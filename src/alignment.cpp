//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

//#define TEST_DETECT 1
//#define TEST_VITERBI 1

#include <fstream>
#include "detect.h"
#include <math.h>
#include <utility>
#include <stdlib.h>
#include <limits>
#include "common.h"
#include "event_handling.h"
#include "../fast5/include/fast5.hpp"
#include "../pod5-file-format/c++/pod5_format/c_api.h"
#include "alignment.h"
#include "error_handling.h"
#include "probability.h"
#include "htsInterface.h"
#include "pod5.h"
#include "fast5.h"
#include "config.h"
#include <mutex>
#include <algorithm>

static const char *help=
"align: DNAscent executable that generates a BrdU- and EdU-aware event alignment.\n"
"To run DNAscent align, do:\n"
"   DNAscent align -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.align\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -o,--output               path to output file that will be generated.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  -m,--maxReads             maximum number of reads to consider,\n"
"  -q,--quality              minimum mapping quality (default is 20),\n"
"  -l,--length               minimum read length in bp (default is 100),\n"
"  --tail-source             which recovered 3' tail signal (past the last full\n"
"                             reference k-mer) to print, for reads that pass the\n"
"                             tail-capture gate (forward-mapped, alignment reaches\n"
"                             the true reference end, no 3' soft clip): none,\n"
"                             source1, source2, or both (default is none).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	bool capReads;
	int minQ, maxReads;
	int minL;
	unsigned int threads;
	TailCaptureMode tailSourceMode;
};


TailCaptureMode parseTailSourceArg( std::string s ){

	if (s == "none") return TailCaptureMode::NONE;
	else if (s == "source1") return TailCaptureMode::SOURCE1;
	else if (s == "source2") return TailCaptureMode::SOURCE2;
	else if (s == "both") return TailCaptureMode::BOTH;
	else throw InvalidOption( "--tail-source " + s );
}

Arguments parseAlignArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent detect." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent detect." << std::endl;
		exit(EXIT_FAILURE);
	}

	Arguments args;

	/*defaults - we'll override these if the option was specified by the user */
	args.threads = 1;
	args.minQ = 20;
	args.minL = 100;
	args.capReads = false;
	args.maxReads = 0;
	args.tailSourceMode = TailCaptureMode::NONE;

	/*parse the command line arguments */

	for ( int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-b" or flag == "--bam" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.bamFilename = strArg;
			i+=2;
		}
		else if ( flag == "-r" or flag == "--reference" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.referenceFilename = strArg;
			i+=2;
		}
		else if ( flag == "-t" or flag == "--threads" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.threads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-q" or flag == "--quality" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.minQ = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-l" or flag == "--length" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.minL = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-i" or flag == "--index" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.indexFilename = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.outputFilename = strArg;
			i+=2;
		}
		else if ( flag == "-m" or flag == "--maxReads" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.capReads = true;
			args.maxReads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "--tail-source" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.tailSourceMode = parseTailSourceArg( strArg );
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	if (args.outputFilename == args.indexFilename or args.outputFilename == args.referenceFilename or args.outputFilename == args.bamFilename) throw OverwriteFailure();

	return args;
}


//counts how many TAIL_SOURCE1/TAIL_SOURCE2/NEXT100EVT lines eventalign() appended to one
//read's output block - used only for the --tail-source terminal summary in align_main.
//This necessarily reflects what the current --tail-source mode actually captured/printed
//for this run, not what signal exists in principle (e.g. under --tail-source source1, a
//read's Source 2 count here will always read zero even if that read has real Source 2
//signal, because it was never collected this run).
void countTailLines(const std::string &out, long &n1, long &n2, long &n3){

	n1 = 0;
	n2 = 0;
	n3 = 0;
	size_t pos = 0;
	while ( (pos = out.find('\n', pos)) != std::string::npos ){

		size_t lineStart = pos + 1;
		if (out.compare(lineStart, 12, "TAIL_SOURCE1") == 0) n1++;
		else if (out.compare(lineStart, 12, "TAIL_SOURCE2") == 0) n2++;
		else if (out.compare(lineStart, 10, "NEXT100EVT") == 0) n3++;
		pos = lineStart;
	}
}


std::string tailCaptureModeToString(TailCaptureMode m){

	switch(m){
		case TailCaptureMode::NONE: return "none";
		case TailCaptureMode::SOURCE1: return "source1";
		case TailCaptureMode::SOURCE2: return "source2";
		case TailCaptureMode::BOTH: return "both";
		default: return "unknown";
	}
}


//min/median/mean/max of a per-read tail-length distribution, for the end-of-run
//--tail-source summary. Empty input (no reads had that source this run) is reported
//explicitly rather than printing meaningless statistics.
void printTailLengthSummary(const std::string &label, std::vector<long> lengths){

	if (lengths.empty()){
		std::cout << "  " << label << ": 0 reads" << std::endl;
		return;
	}
	std::sort(lengths.begin(), lengths.end());
	double sum = 0.;
	for (long v : lengths) sum += v;
	double mean = sum / lengths.size();
	long median = lengths[lengths.size() / 2];
	std::cout << "  " << label << ": " << lengths.size() << " reads, min=" << lengths.front()
	           << " median=" << median << " mean=" << mean << " max=" << lengths.back() << std::endl;
}


//called from both align_main exit points (the early return under -m/--maxReads, and the
//normal end-of-BAM return) so the summary always prints regardless of which one is hit.
void printTailCaptureSummary(TailCaptureMode mode, long readsPassingTailGate,
                              long readsWithSource1, long readsWithSource2, long readsWithNext100,
                              long totalSource1Samples, long totalSource2Samples, long totalNext100Samples,
                              const std::vector<long> &source1Lengths,
                              const std::vector<long> &source2Lengths,
                              const std::vector<long> &next100Lengths){

	std::cout << std::endl << "Tail-source capture summary (--tail-source "
	          << tailCaptureModeToString(mode) << "):" << std::endl;
	std::cout << "  reads passing tail gate (forward-mapped, reaches true ref end, no 3' soft clip): "
	          << readsPassingTailGate << std::endl;
	std::cout << "  reads with Source 1 signal printed this run: " << readsWithSource1
	          << ", total Source 1 samples: " << totalSource1Samples << std::endl;
	std::cout << "  reads with Source 2 signal printed this run: " << readsWithSource2
	          << ", total Source 2 samples: " << totalSource2Samples << std::endl;
	std::cout << "  reads with NEXT100EVT signal printed this run: " << readsWithNext100
	          << ", total NEXT100EVT samples (event-mean granularity, <=100 events/read): " << totalNext100Samples << std::endl;
	printTailLengthSummary("Source 1 tail length per read", source1Lengths);
	printTailLengthSummary("Source 2 tail length per read", source2Lengths);
	printTailLengthSummary("NEXT100EVT tail length per read", next100Lengths);
}


inline double lnVecMax(std::vector<double> v){

	double maxVal = v[0];
	for (size_t i = 1; i < v.size(); i++){
		if (lnGreaterThan( v[i], maxVal )){
			maxVal = v[i];
		}
	}
	return maxVal;
}


inline int lnArgMax(std::vector<double> v){

	double maxVal = v[0];
	int maxarg = 0;

	for (size_t i = 1; i < v.size(); i++){
		if (lnGreaterThan( v[i], maxVal )){
			maxVal = v[i];
			maxarg = i;
		}
	}
	return maxarg;
}


std::pair< double, std::vector< std::string > > builtinViterbi( std::vector <double> &observations,
				std::string &sequence,
				PoreParameters scalings,
				bool flip){
	
	//HMM transition probabilities	
	double externalD2D = eln(Pore_Substrate_Config.HMM_config.externalD2D);
	double externalD2M1 = eln(Pore_Substrate_Config.HMM_config.externalD2M);
	double externalI2M1 = eln(Pore_Substrate_Config.HMM_config.externalI2M);
	double externalM12D = eln(Pore_Substrate_Config.HMM_config.externalM2D);
	double internalM12I = eln(Pore_Substrate_Config.HMM_config.internalM2I);
	double internalI2I = eln(Pore_Substrate_Config.HMM_config.internalI2I);

	//transition probabilities that change on a per-read basis
	double internalM12M1 = eln(1. - (1./scalings.eventsPerBase));
	double externalM12M1 = eln(1.0 - externalM12D - internalM12I - internalM12M1);
	double externalM12M1orD = lnSum( externalM12M1, externalM12D );
	double externalOrInternalM12M1 = lnSum( externalM12M1, internalM12M1 );
	int maxindex;

	unsigned int k = Pore_Substrate_Config.kmer_len;

	size_t n_states = sequence.length() - k + 1;

	//pre-compute 6mer indices
	std::vector<unsigned int> kmerIndices;
	kmerIndices.reserve(n_states);
	for (size_t i = 0; i < n_states; i++){
		std::string kmer = sequence.substr(i, k);
		if (flip) std::reverse(kmer.begin(),kmer.end());
		kmerIndices.push_back(kmer2index(kmer,k));
	}

	std::vector< std::vector< size_t > > backtraceS( 3*n_states, std::vector< size_t >( observations.size() + 1 ) ); /*stores state indices for the Viterbi backtrace */
	std::vector< std::vector< size_t > > backtraceT( 3*n_states, std::vector< size_t >( observations.size() + 1 ) ); /*stores observation indices for the Viterbi backtrace */

	//reserve 0 for start
	ssize_t D_offset = 0;
	ssize_t M_offset = n_states;
	ssize_t I_offset = 2*n_states;

	std::vector< double > I_curr(n_states, NAN), D_curr(n_states, NAN), M_curr(n_states, NAN), I_prev(n_states, NAN), D_prev(n_states, NAN), M_prev(n_states, NAN);
	double start_curr = NAN, start_prev = 0.0;

	double matchProb, insProb;

	/*-----------INITIALISATION----------- */
	//transitions from the start state
	D_prev[0] = lnProd( start_prev, externalM12D );
	backtraceS[0 + D_offset][0] = -1;
	backtraceT[0 + D_offset][0] = 0;

	//account for transitions between deletion states before we emit the first observation
	for ( unsigned int i = 1; i < n_states; i++ ){

		D_prev[i] = D_prev[i-1] + externalD2D;
		backtraceS[i + D_offset][0] = i -1 + D_offset;
		backtraceT[i + D_offset][0] = 0;
	}


	/*-----------RECURSION----------- */
	/*complexity is O(T*N^2) where T is the number of observations and N is the number of states */
	double level_mu, level_sigma;
	for ( unsigned int t = 0; t < observations.size(); t++ ){

		std::fill( I_curr.begin(), I_curr.end(), NAN );
		std::fill( M_curr.begin(), M_curr.end(), NAN );
		std::fill( D_curr.begin(), D_curr.end(), NAN );

		std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmerIndices[0]];

		//uncomment if you scale model
		//level_mu = (scalings.shift + scalings.scale * meanStd.first);
		//level_sigma = meanStd.second;

		//uncomment if you scale events
		level_mu = meanStd.first;
		level_sigma = meanStd.second;

		matchProb = eln( normalPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
		//matchProb = eln( cauchyPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
		insProb = 0.0; //log(1) = 0; set probability equal to 1 and use transition probability as weighting

		//to the base 1 insertion
		I_curr[0] = lnVecMax({I_prev[0] + internalI2I + insProb ,
			                  M_prev[0] + internalM12I + insProb,
							  start_prev + internalM12I + insProb
		                     });
		maxindex = lnArgMax({I_prev[0] + internalI2I + insProb,
					         M_prev[0] + internalM12I + insProb,
							 start_prev + internalM12I + insProb*0.001
						     });
		switch(maxindex){
			case 0:
				backtraceS[0 + I_offset][t+1] = 0 + I_offset;
				backtraceT[0 + I_offset][t+1] = t;
				break;
			case 1:
				backtraceS[0 + I_offset][t+1] = 0 + M_offset;
				backtraceT[0 + I_offset][t+1] = t;
				break;
			case 2:
				backtraceS[0 + I_offset][t+1] = -1;
				backtraceT[0 + I_offset][t+1] = t;
				break;
			default:
				std::cout << "problem" << std::endl;
				exit(EXIT_FAILURE);
		}

		//to the base 1 match
		M_curr[0] = lnVecMax({M_prev[0] + internalM12M1 + matchProb,
							  start_prev + externalOrInternalM12M1 + matchProb
							 });
		maxindex = lnArgMax({M_prev[0] + internalM12M1 + matchProb,
							 start_prev + externalOrInternalM12M1 + matchProb
							});
		switch(maxindex){
			case 0:
				backtraceS[0 + M_offset][t+1] = 0 + M_offset;
				backtraceT[0 + M_offset][t+1] = t;
				break;
			case 1:
				backtraceS[0 + M_offset][t+1] = -1;
				backtraceT[0 + M_offset][t+1] = t;
				break;
			default:
				std::cout << "problem" << std::endl;
				exit(EXIT_FAILURE);
		}

		//to the base 1 deletion
		D_curr[0] = lnProd( NAN, externalM12D );  //start to D
		backtraceS[0 + D_offset][t+1] = -1;
		backtraceT[0 + D_offset][t+1] = t + 1;


		//the rest of the sequence
		for ( unsigned int i = 1; i < n_states; i++ ){

			insProb = 0.0; //log(1) = 0; set probability equal to 1 and use transition probability as weighting

			//get model parameters
			std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmerIndices[i]];

			//uncomment if you scale model
			//level_mu = scalings.shift + scalings.scale * meanStd.first;
			//level_sigma = meanStd.second;

			//uncomment if you scale events
			level_mu = meanStd.first;
			level_sigma = meanStd.second;

			matchProb = eln( normalPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
			//matchProb = eln( cauchyPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
			
			//to the insertion
			I_curr[i] = lnVecMax({I_prev[i] + internalI2I + insProb,
								  M_prev[i] + internalM12I + insProb
								 });
			maxindex = lnArgMax({I_prev[i] + internalI2I + insProb ,
							     M_prev[i] + internalM12I + insProb
								});
			switch(maxindex){
				case 0:
					backtraceS[i + I_offset][t+1] = i + I_offset;
					backtraceT[i + I_offset][t+1] = t;
					break;
				case 1:
					backtraceS[i + I_offset][t+1] = i + M_offset;
					backtraceT[i + I_offset][t+1] = t;
					break;
				default:
					std::cout << "problem" << std::endl;
					exit(EXIT_FAILURE);
			}

			//to the match
			M_curr[i] = lnVecMax({I_prev[i-1] + externalI2M1 + matchProb,
								   M_prev[i-1] + externalM12M1 + matchProb ,
								   M_prev[i] + internalM12M1 + matchProb ,
								   D_prev[i-1] + externalD2M1 + matchProb
								   });
			maxindex = lnArgMax({I_prev[i-1] + externalI2M1 + matchProb,
							   M_prev[i-1] + externalM12M1 + matchProb,
							   M_prev[i] + internalM12M1 + matchProb,
							   D_prev[i-1] + externalD2M1 + matchProb
			});
			switch(maxindex){
				case 0:
					backtraceS[i + M_offset][t+1] = i - 1 + I_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				case 1:
					backtraceS[i + M_offset][t+1] = i - 1 + M_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				case 2:
					backtraceS[i + M_offset][t+1] = i + M_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				case 3:
					backtraceS[i + M_offset][t+1] = i - 1 + D_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				default:
					std::cout << "problem" << std::endl;
					exit(EXIT_FAILURE);
			}
		}

		for ( unsigned int i = 1; i < n_states; i++ ){

			//to the deletion
			D_curr[i] = lnVecMax({M_curr[i-1] + externalM12D,
				                  D_curr[i-1] + externalD2D
			                     });
			maxindex = lnArgMax({ M_curr[i-1] + externalM12D,
                               D_curr[i-1] + externalD2D
			                  });
			switch(maxindex){
				case 0:
					backtraceS[i + D_offset][t+1] = i - 1 + M_offset;
					backtraceT[i + D_offset][t+1] = t + 1;
					break;
				case 1:
					backtraceS[i + D_offset][t+1] = i - 1 + D_offset;
					backtraceT[i + D_offset][t+1] = t + 1;
					break;
				default:
					std::cout << "problem" << std::endl;
					exit(EXIT_FAILURE);
			}
		}

		I_prev = I_curr;
		M_prev = M_curr;
		D_prev = D_curr;
		start_prev = start_curr;
	}

	/*
	for (int i = 0; i < backtraceS.size(); i++){
		for (int j = 0; j < backtraceS[i].size(); j++){
			std::cout << backtraceS[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << "----------------------------------------------------------------------" << std::endl;
	*/

	/*-----------TERMINATION----------- */
	double viterbiScore = NAN;
	viterbiScore = lnVecMax( {D_curr.back() , // + eln( 1.0 ) which is 0 //D to end
							  M_curr.back() + externalM12M1orD,//M to end
							  I_curr.back() + externalI2M1 //I to end
							 });
	//std::cout << "Builtin Viterbi score: " << viterbiScore << std::endl;


	//figure out where to go from the end state
	maxindex = lnArgMax({ D_curr.back() , //+ eln( 1.0 ) which is 0
	                   M_curr.back() + externalM12M1orD ,
					   I_curr.back() + externalI2M1
	                   });

	ssize_t  traceback_new;
	ssize_t  traceback_old;
	ssize_t  traceback_t = observations.size();
	switch(maxindex){
		case 0:
			traceback_old = D_offset + n_states - 1;
			break;
		case 1:
			traceback_old = M_offset + n_states - 1;
			break;
		case 2:
			traceback_old = I_offset + n_states - 1;
			break;
		default:
			std::cout << "problem" << std::endl;
			exit(EXIT_FAILURE);
	}

#if TEST_VITERBI
std::cerr << "Starting traceback..." << std::endl;
std::cerr << "Number of events: " << observations.size() << std::endl;
#endif

	std::vector<std::string> stateIndices;
	stateIndices.reserve(observations.size());
	while (traceback_old != -1){

		traceback_new = backtraceS[ traceback_old ][ traceback_t ];
		traceback_t = backtraceT[ traceback_old ][ traceback_t ];

		if (traceback_old < M_offset){ //Del

			//std::cout << "D " << traceback_old << std::endl;
			stateIndices.push_back(std::to_string(traceback_old) + "_D");

		}
		else if (traceback_old < I_offset){ //M

			//std::cout << "M " << traceback_old - M_offset << std::endl;
			stateIndices.push_back(std::to_string(traceback_old- M_offset) + "_M");
		}
		else { //I

			//std::cout << "I " << traceback_old - I_offset << std::endl;
			stateIndices.push_back(std::to_string(traceback_old - I_offset) + "_I");

		}
		traceback_old = traceback_new;
	}
	std::reverse( stateIndices.begin(), stateIndices.end() );

#if TEST_VITERBI
std::cout << "Traceback terminated." << std::endl;
#endif

	return std::make_pair( viterbiScore,stateIndices);
}


bool hasSoftClipAtReferenceThreePrime(bam1_t *record){

	//BAM CIGAR is always stored in reference/genome coordinate order regardless of the
	//read's mapped strand, so the reference's 3' end corresponds to the LAST cigar
	//operation for a forward-mapped read and the FIRST cigar operation for a reverse-
	//mapped read.
	uint32_t n_cigar = record -> core.n_cigar;
	if (n_cigar == 0) return false;

	const uint32_t *cigar = bam_get_cigar(record);
	int op = bam_is_rev(record) ? bam_cigar_op(cigar[0]) : bam_cigar_op(cigar[n_cigar - 1]);

	return op == BAM_CSOFT_CLIP;
}


bool referenceDefined(std::string &readSnippet){

	//make sure the read snippet is fully defined as A/T/G/C in reference
	unsigned int As = 0, Ts = 0, Cs = 0, Gs = 0;
	for ( std::string::iterator i = readSnippet.begin(); i < readSnippet.end(); i++ ){

		switch( *i ){
			case 'A' :
				As++;
				break;
			case 'T' :
				Ts++;
				break;
			case 'G' :
				Gs++;
				break;
			case 'C' :
				Cs++;
				break;
		}
	}
	if ( readSnippet.length() != (As + Ts + Gs + Cs) ){
		return false;
	}
	else return true;
}


void eventalign( DNAscent::read &r, unsigned int totalWindowLength, TailCaptureMode tailMode){

	int readHead = 0;
	int runningInsertions = 0;
	int runningDeletions = 0;

	unsigned int k = Pore_Substrate_Config.kmer_len;

	r.humanReadable_eventalignOut = ">" + r.readID + " " + r.referenceMappedTo + " " + std::to_string(r.refStart) + " " + std::to_string(r.refEnd) + " " + r.strand + "\n";

	//raw scaled signal observed after the last reference-anchored 9-mer window - the pore
	//model needs a full k-mer of context to centre a call, so nothing beyond that point ever
	//gets a reference coordinate. Two distinct mechanisms feed signal in here (see the two
	//collection sites below), so they're kept in separate buffers rather than one combined
	//pool: Source 1 (terminal-window trailing insertions from the reference-based Viterbi
	//HMM) and Source 2 (events the rough, basecall-level aligner trimmed before eventalign
	//ever saw them, often much larger than Source 1 - see TAIL_SIGNAL_CAPTURE.md). Each
	//source is only collected if tailMode actually asks for it, so requesting just one
	//source skips the other's collection work entirely rather than gathering-then-discarding.
	bool wantSource1 = (tailMode == TailCaptureMode::SOURCE1 or tailMode == TailCaptureMode::BOTH);
	bool wantSource2 = (tailMode == TailCaptureMode::SOURCE2 or tailMode == TailCaptureMode::BOTH);
	std::vector<double> tailSignalSource1;
	std::vector<double> tailSignalSource2;

	//PI-proposed anchor: eventIndeces[lastM_ev] from the LAST window processed below -
	//the raw r.events index of the last event the HMM actually matched to a reference
	//k-mer (not merely the last event the rough aligner touched, which is what Source 2
	//uses and can be later if the terminal window's HMM assigned any trailing events as
	//"I" rather than extending the match). Overwritten unconditionally every iteration,
	//so after the loop below it holds whatever the LAST iteration computed.
	size_t lastMatchRawIdx = 0;
	bool haveLastMatchRawIdx = false;

	unsigned int reference_index = 0;
	while ( reference_index < r.referenceSeqMappedTo.size() - k + 1){

		double insRate = (runningInsertions) / double(readHead+1);
		if (insRate > 0.2){
			r.QCpassed = false;
			return;
		}

		//adjust so we can get the last bit of the read if it doesn't line up with the windows nicely
		unsigned int basesToEnd = r.referenceSeqMappedTo.size() - reference_index;
		unsigned int windowLength = std::min(basesToEnd, totalWindowLength);
	
		//find good breakpoints
		std::string break1, break2;
		if (basesToEnd > 1.5*totalWindowLength){

			std::string breakSnippet = (r.referenceSeqMappedTo).substr(reference_index, 1.5*windowLength);

			bool isDefined = referenceDefined(breakSnippet);
			if (not isDefined){
				reference_index += windowLength;
				continue;
			}

			for (unsigned int i = windowLength; i < 1.5*windowLength - k - 1; i++){

				std::string kmer = breakSnippet.substr(i,k);
				std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmer2index(kmer, k)];

				std::string kmer_back = breakSnippet.substr(i-1,k);
				std::pair<double,double> meanStd_back = Pore_Substrate_Config.pore_model[kmer2index(kmer_back, k)];

				std::string kmer_front = breakSnippet.substr(i+1,k);
				std::pair<double,double> meanStd_front = Pore_Substrate_Config.pore_model[kmer2index(kmer_front, k)];

				double gap1 = std::abs(meanStd.first - meanStd_front.first);
				double gap2 = std::abs(meanStd.first - meanStd_back.first);

				if (gap1 > 0.75 and gap2 > 0.75){
					break1 = breakSnippet.substr(i-1,k);
					break2 = breakSnippet.substr(i+1,k);
					windowLength = i + k;
					break;
				}
			}
		}

		std::string readSnippet = (r.referenceSeqMappedTo).substr(reference_index, windowLength);

		bool isDefined = referenceDefined(readSnippet);
		if (not isDefined){

			reference_index += windowLength;
			continue;
		}

		std::vector< double > eventSnippet_means;
		std::vector< event > eventSnippet;
		//PI-proposed tracking: same length as eventSnippet_means/eventSnippet, holding
		//each entry's index into r.events (the read's full, 5'->3' event vector) - lets
		//us translate a LOCAL position within this window's Viterbi output (e.g. lastM_ev
		//below) back into a raw index in the full per-read event vector, which eventSnippet
		//itself discards.
		std::vector< size_t > eventIndeces;

		//get the events that correspond to the read snippet
		bool firstMatch = true;
		for ( unsigned int j = readHead; j < (r.eventAlignment).size(); j++ ){

			//if an event has been aligned to a position in the window, add it
			if ( (r.refToQuery)[reference_index] <= (r.eventAlignment)[j].second and (r.eventAlignment)[j].second < (r.refToQuery)[reference_index + windowLength - k + 1] ){

				if (firstMatch){
					readHead = j;
					firstMatch = false;
				}

				double event_mean = (r.events)[(r.eventAlignment)[j].first].mean;

				//guard on bad signal
				if (0. < event_mean and event_mean < 250.){
					eventSnippet_means.push_back(event_mean);
					eventSnippet.push_back((r.events)[(r.eventAlignment)[j].first]);
					eventIndeces.push_back((r.eventAlignment)[j].first);
				}
			}

			//stop once we get to the end of the window
			if ( (r.eventAlignment)[j].second >= (r.refToQuery)[reference_index + windowLength - k + 1] ) break;
		}
	
		//flag large insertions
		int querySpan = (r.refToQuery)[reference_index + windowLength - k + 1] - (r.refToQuery)[reference_index];
		assert(querySpan >= 0);
		int referenceSpan = windowLength - k + 1;
		int indelScore = querySpan - referenceSpan;

		//pass on this window if we have a deletion
		if ( eventSnippet_means.size() < 2){

			reference_index += windowLength;
			continue;
		}

		//calculate where we are on the assembly - if we're a reverse complement, we're moving backwards down the reference genome
		int reference_coord;
		if ( r.isReverse ) reference_coord = r.refEnd - reference_index - k/2;
		else reference_coord = r.refStart + reference_index + k/2;
		//std::cout << "Event sizes: " << eventSnippet_means.size() << " " << eventSnippet.size() << std::endl;
		std::pair< double, std::vector<std::string> > builtinAlignment = builtinViterbi( eventSnippet_means, readSnippet, r.scalings, false);

		std::vector< std::string > stateLabels = builtinAlignment.second;
		size_t lastM_ev = 0;
		size_t lastM_ref = 0;

		size_t evIdx = 0;

		//grab the index of the last match so we don't print insertions where we shouldn't
		for (size_t i = 0; i < stateLabels.size(); i++){

			std::string label = stateLabels[i].substr(stateLabels[i].find('_')+1);
			int pos = std::stoi(stateLabels[i].substr(0,stateLabels[i].find('_')));

			if (label == "M"){
				lastM_ev = evIdx;
				lastM_ref = pos;
			}

			if (label != "D") evIdx++; //silent states don't emit an event
		}

		//translate this window's local lastM_ev back into a raw r.events index - see
		//lastMatchRawIdx's declaration above the outer while loop.
		if (not eventIndeces.empty()){
			lastMatchRawIdx = eventIndeces[lastM_ev];
			haveLastMatchRawIdx = true;
		}

		//true only if this window's last match reaches the last valid k-mer in the whole
		//read - equivalently, true only if the outer while loop is about to terminate for
		//good after this window. windowLength reaching basesToEnd is NOT sufficient on its
		//own: if Viterbi D-labels some trailing states in a window that nominally reached
		//the end, lastM_ref falls short, reference_index only advances partway, and the
		//outer loop runs again on the leftover - readHead does not skip past that leftover's
		//events, so a genuinely-final-looking window here can still be followed by another
		//real window that needs those same trailing events back.
		bool reachedFinalKmer = (reference_index + lastM_ref + 1 >= r.referenceSeqMappedTo.size() - k + 1);

		//do a second pass to print the alignment
		evIdx = 0;
		for (size_t i = 0; i < stateLabels.size(); i++){

			std::string label = stateLabels[i].substr(stateLabels[i].find('_')+1);
	        int pos = std::stoi(stateLabels[i].substr(0,stateLabels[i].find('_')));

			//silent states don't emit an event
			if (label == "D"){
				runningDeletions++;
				continue;
			}

			std::string kmerStrand = (r.referenceSeqMappedTo).substr(reference_index + pos, k);

			//calculate the kmer on the strand as well as the reference coordinate for this event
			unsigned int event_coord;
			std::string kmerRef;
			if (r.isReverse){
				event_coord = reference_coord - pos - 1;
				kmerRef = reverseComplement(kmerStrand);
			}
			else{
				event_coord = reference_coord + pos;
				kmerRef = kmerStrand;
			}
			
			//calculate the query (0-based) index from this reference position
			unsigned int event_indexRef = reference_index + pos + k/2;
			unsigned event_indexQuery = r.refToQuery.at(event_indexRef);
			
			if (label == "M"){
				std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmer2index(kmerStrand, k)];

				for (unsigned int idx_raw = 0; idx_raw < eventSnippet[evIdx].raw.size(); idx_raw++){
					//std::cout << evIdx << " " << eventSnippet.size() << " " << idx_raw << " " << eventSnippet[evIdx].raw.size() << std::endl;
					double scaledEvent = (eventSnippet[evIdx].raw[idx_raw] - r.scalings.shift) / r.scalings.scale;

					if (r.refCoordToCalls.count(event_coord) > 0){
						r.humanReadable_eventalignOut += std::to_string(event_coord) 
							      + "\t" + kmerRef 
							      + "\t" + std::to_string(scaledEvent) 
							      + "\t" + kmerStrand 
							      + "\t" + std::to_string(meanStd.first) 
							      + "\t" + std::to_string(r.refCoordToCalls.at(event_coord).first) 
							      + "\t" + std::to_string(r.refCoordToCalls.at(event_coord).second) 					          
							      + "\n";
					}
					else{
						r.humanReadable_eventalignOut += std::to_string(event_coord) 
							      + "\t" + kmerRef 
							      + "\t" + std::to_string(scaledEvent) 
							      + "\t" + kmerStrand 
							      + "\t" + std::to_string(meanStd.first) 
							      + "\n";
						r.addSignal(kmerStrand, event_coord, event_indexQuery, event_indexRef, scaledEvent, indelScore);
					}

				}

			}
			else if (label == "I"){

				//insertions after the last match don't get printed against a reference
				//coordinate because normally the next window picks them up. For the
				//terminal window there is no next window - Viterbi still routed these
				//events somewhere (there's no other state left for it to put them in),
				//it's just that "somewhere" is signal from past the end of the modelled
				//region, so route it to the TAIL pool instead of dropping it.
				bool isTrailingInsertion = (evIdx >= lastM_ev);

				if ( (not isTrailingInsertion) or (reachedFinalKmer and wantSource1) ){

					for (unsigned int idx_raw = 0; idx_raw < eventSnippet[evIdx].raw.size(); idx_raw++){
						double scaledEvent = (eventSnippet[evIdx].raw[idx_raw] - r.scalings.shift) / r.scalings.scale;

						if (isTrailingInsertion) tailSignalSource1.push_back(scaledEvent);
						else r.humanReadable_eventalignOut += std::to_string(event_coord) + "\t" + kmerRef + "\t" + std::to_string(scaledEvent) + "\t" + std::string(k, 'N') + "\t" + "0" + "\n";
					}
					runningInsertions++;
				}
			}

			evIdx ++;
		}

		//go again starting at posOnRef + lastM_ref using events starting at readHead + lastM_ev
		readHead += lastM_ev + 1;
		reference_index += lastM_ref + 1;
	}

	//a second, distinct pool of tail signal: events that never made it into r.eventAlignment
	//at all. The adaptive banded alignment (event_handling.cpp, adaptive_banded_simple_event_align)
	//picks a best-scoring final event index and "trims" (silently discards) anything after it
	//rather than force-matching it to the last k-mer if that scores worse than trimming -
	//exactly the situation we'd expect if a modified base right at the end of the molecule
	//makes the tail signal look unlike the last unmodified k-mer. Those events are otherwise
	//invisible to this function, since everything above only ever walks r.eventAlignment.
	if ( wantSource2 and r.eventAlignment.size() > 0 ){

		size_t lastAlignedEvent = r.eventAlignment.back().first;
		for ( size_t evi = lastAlignedEvent + 1; evi < r.events.size(); evi++ ){

			double event_mean = r.events[evi].mean;
			if ( not (0. < event_mean and event_mean < 250.) ) continue; //same signal guard used when gathering eventSnippet above

			for ( unsigned int idx_raw = 0; idx_raw < r.events[evi].raw.size(); idx_raw++ ){
				double scaledEvent = (r.events[evi].raw[idx_raw] - r.scalings.shift) / r.scalings.scale;
				tailSignalSource2.push_back(scaledEvent);
			}
		}
	}

	//a third, independent pool anchored at lastMatchRawIdx - the raw r.events index of
	//the last event the HMM actually matched to a reference k-mer, in the last window
	//processed above (PI-proposed anchor, eventIndeces[lastM_ev]). This is NOT always the
	//same as Source 2's boundary (r.eventAlignment.back().first): if the terminal
	//window's HMM assigned some trailing events as "I" (Source 1) rather than extending
	//the match, those events sit BETWEEN lastMatchRawIdx and r.eventAlignment.back().first,
	//so this window starts earlier than Source 2's and folds any Source 1 signal into its
	//first entries. Reported at raw EVENT granularity (one value per event.mean) rather
	//than expanded raw ADC samples, capped at the next 100 events.
	std::vector<double> tailSignalNext100Events;
	if ( (wantSource1 or wantSource2) and haveLastMatchRawIdx ){

		size_t stop = std::min( r.events.size(), lastMatchRawIdx + 1 + 100 );
		for ( size_t evi = lastMatchRawIdx + 1; evi < stop; evi++ ){

			double event_mean = r.events[evi].mean;
			if ( not (0. < event_mean and event_mean < 250.) ) continue; //same signal guard used elsewhere

			double scaledEvent = (event_mean - r.scalings.shift) / r.scalings.scale;
			tailSignalNext100Events.push_back(scaledEvent);
		}
	}

	//Each buffer is independently in chronological order. No reference coordinate is
	//meaningful for either source, so lines lead with a non-numeric marker (existing parsers
	//that blindly do int(fields[0]) on every line need to skip these first). wantSource1/
	//wantSource2 already encode exactly which mode(s) ask for which buffer, so printing just
	//mirrors collection above. Source 1 is written before Source 2 when both are requested,
	//matching the old combined buffer's chronological construction (terminal-window
	//insertions first, then rough-alignment-trimmed events). NEXT100EVT deliberately does
	//NOT start with "TAIL" - existing scripts key off a bare "TAIL" prefix to mean "legacy
	//combined tail signal" and would otherwise silently double-count this pool as more of
	//that same signal.
	if (wantSource1){
		for ( size_t i = 0; i < tailSignalSource1.size(); i++ ){
			r.humanReadable_eventalignOut += "TAIL_SOURCE1\t" + std::to_string(i) + "\t" + std::to_string(tailSignalSource1[i]) + "\n";
		}
	}
	if (wantSource2){
		for ( size_t i = 0; i < tailSignalSource2.size(); i++ ){
			r.humanReadable_eventalignOut += "TAIL_SOURCE2\t" + std::to_string(i) + "\t" + std::to_string(tailSignalSource2[i]) + "\n";
		}
	}
	if (wantSource1 or wantSource2){
		for ( size_t i = 0; i < tailSignalNext100Events.size(); i++ ){
			r.humanReadable_eventalignOut += "NEXT100EVT\t" + std::to_string(i) + "\t" + std::to_string(tailSignalNext100Events[i]) + "\n";
		}
	}

	r.QCpassed = true;
}


int align_main( int argc, char** argv ){

	Arguments args = parseAlignArguments( argc, argv );

	//load DNAscent index
	std::map< std::string, IndexEntry > readID2path;
	parseIndex( args.indexFilename, readID2path );

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( args.referenceFilename );

	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	//load the bam
	std::cout << "Opening bam file... ";
	htsFile *bam_fh_cr = sam_open((args.bamFilename).c_str(), "r");
	if (bam_fh_cr == NULL) throw IOerror(args.bamFilename);

	//load the header
	bam_hdr_t *bam_hdr_cr = sam_hdr_read(bam_fh_cr);
	std::cout << "ok." << std::endl;

	//open a log file
	std::cout << "Opening log file... ";
	std::string logFilename = strip_extension(args.outputFilename);
	logFilename += ".align.log";
	std::ofstream logfile(logFilename);
	if (logfile.is_open()) std::cout << "ok." << std::endl;
	else throw IOerror(logFilename);
	std::mutex mtx;

	//initialise progress
	int numOfRecords = 0, prog = 0, failed = 0;
	countRecords( bam_fh_cr, bam_hdr_cr, numOfRecords, args.minQ, args.minL );
	if (args.capReads){
		numOfRecords = std::min(numOfRecords,args.maxReads);
	}
	progressBar pb(numOfRecords,true);
	bam_hdr_destroy(bam_hdr_cr);
	hts_close(bam_fh_cr);

	pod5_init();

	int failedEvents = 0;

	//accumulated only inside the #pragma omp critical block below, alongside the existing
	//outFile/prog/pb updates - see the terminal summary printed at the end of this function.
	long readsPassingTailGate = 0;
	long readsWithSource1 = 0, readsWithSource2 = 0, readsWithNext100 = 0;
	long totalSource1Samples = 0, totalSource2Samples = 0, totalNext100Samples = 0;
	std::vector<long> source1Lengths, source2Lengths, next100Lengths;

	unsigned int maxBufferSize;
	std::vector< bam1_t * > buffer;
	if ( args.threads <= 4 ) maxBufferSize = args.threads;
	else maxBufferSize = 4*(args.threads);

	htsFile *bam_fh = sam_open((args.bamFilename).c_str(), "r");
	if (bam_fh == NULL) throw IOerror(args.bamFilename);
	bam_hdr_t *bam_hdr = sam_hdr_read(bam_fh);
	bam1_t *itr_record = bam_init1();
	int result = sam_read1(bam_fh, bam_hdr, itr_record);

	//int readCount = 0;


	while(result >= 0){
	
		bam1_t *record = bam_dup1(itr_record);

		//add the record to the buffer if it passes the user's criteria, otherwise destroy it cleanly
		int mappingQual = record -> core.qual;
		int refStart,refEnd;
		getRefEnd(record,refStart,refEnd);
		int queryLen = record -> core.l_qseq;

		if ( mappingQual >= args.minQ and refEnd - refStart >= args.minL and queryLen != 0 ){

			buffer.push_back(record);
		}
		else{
			bam_destroy1(record);
		}

		result = sam_read1(bam_fh, bam_hdr, itr_record);
				
		//if we've filled up the buffer with reads, compute them in parallel
		if (buffer.size() >= maxBufferSize or (buffer.size() > 0 and result == -1 ) ){

			#pragma omp parallel for schedule(dynamic) shared(buffer,Pore_Substrate_Config,args,prog,failed) num_threads(args.threads)
			for (unsigned int i = 0; i < buffer.size(); i++){

				DNAscent::read r(buffer[i], bam_hdr, readID2path, reference);

				if (r.refMismatch){

					std::cerr << std::endl << "Error: contig '" << r.referenceMappedTo << "' found in BAM file but not in the reference genome." << std::endl;
					std::cerr << "Please check that the reference genome passed with -r matches the one used for alignment." << std::endl;
					exit(1);
				}
				if (r.missing){

					std::lock_guard<std::mutex> lock(mtx);
					logfile << "ReadID " << r.readID << " missing from index. Skipping." << std::endl;
					prog++;
					continue;
				}

				if (r.missing){

					std::cerr << "ReadID " << r.readID << " missing from index. Skipping." << std::endl;
					prog++;
					continue;
				}

				const char *ext = get_ext(r.filename.c_str());
				
				if (strcmp(ext,"pod5") == 0){
					pod5_getSignal(r);
				}
				else if (strcmp(ext,"fast5") == 0){
					fast5_getSignal(r);
				}

				bool useFitPoreModel = false;
				normaliseEvents(r, useFitPoreModel);

				//catch reads with rough event alignments that fail the QC
				if ( r.eventAlignment.size() == 0 ){

					failed++;
					prog++;
					continue;
				}

				//only trust signal past the end of the modelled region if this read is
				//forward-mapped, its alignment actually reaches the true 3' end of the
				//reference contig, and there's no 3' soft clip - a 3' soft clip here would
				//mean the physical molecule kept going past our reference (e.g. extra
				//ligated sequence), so what follows isn't unambiguously "past position N
				//of our construct."
				bool readPassesTailGate = (r.strand == "fwd")
				                          and (r.refEnd == (int) reference.at(r.referenceMappedTo).size())
				                          and not hasSoftClipAtReferenceThreePrime(r.record);

				//the gate above is a hard read-level requirement regardless of what the
				//user asked for on the command line - a read that fails it never gets any
				//tail signal printed, no matter which --tail-source mode is in effect.
				TailCaptureMode effectiveTailMode = readPassesTailGate ? args.tailSourceMode : TailCaptureMode::NONE;

				eventalign(r, Pore_Substrate_Config.windowLength_align, effectiveTailMode);

				if (not r.QCpassed){
					failed++;
					prog++;
					continue;
				}

				#pragma omp critical
				{
					outFile << r.humanReadable_eventalignOut;
					prog++;
					pb.displayProgress( prog, failed, failedEvents );

					if (readPassesTailGate){

						readsPassingTailGate++;
						long n1 = 0, n2 = 0, n3 = 0;
						countTailLines(r.humanReadable_eventalignOut, n1, n2, n3);
						if (n1 > 0){
							readsWithSource1++;
							totalSource1Samples += n1;
							source1Lengths.push_back(n1);
						}
						if (n2 > 0){
							readsWithSource2++;
							totalSource2Samples += n2;
							source2Lengths.push_back(n2);
						}
						if (n3 > 0){
							readsWithNext100++;
							totalNext100Samples += n3;
							next100Lengths.push_back(n3);
						}
					}
				}
			}
			buffer.clear();
		}
		pb.displayProgress( prog, failed, failedEvents );
		if (args.capReads and prog > args.maxReads){
			bam_destroy1(itr_record);
			bam_hdr_destroy(bam_hdr);
			hts_close(bam_fh);
			printTailCaptureSummary(args.tailSourceMode, readsPassingTailGate, readsWithSource1,
			                        readsWithSource2, readsWithNext100, totalSource1Samples,
			                        totalSource2Samples, totalNext100Samples,
			                        source1Lengths, source2Lengths, next100Lengths);
			return 0;
		}
	}
	bam_destroy1(itr_record);
	bam_hdr_destroy(bam_hdr);
	hts_close(bam_fh);
	std::cout << std::endl;
	pod5_terminate();
	logfile.close();
	printTailCaptureSummary(args.tailSourceMode, readsPassingTailGate, readsWithSource1,
	                        readsWithSource2, readsWithNext100, totalSource1Samples,
	                        totalSource2Samples, totalNext100Samples,
	                        source1Lengths, source2Lengths, next100Lengths);
	return 0;
}
