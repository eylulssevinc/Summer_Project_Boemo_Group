//----------------------------------------------------------
// Copyright 2026 University of Cambridge
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

/*** Gillespie simulation of the following beacon calculus model:
 *
 * fast = 100000;
 * m = 2000;
 * v = 1.4;
 * fr = 0.00001;
 *
 * Fp[i] = {chr![i],fast}.[i <= m] -> {~chr?[i+1],v}.Fp[i+1];
 * Fm[i] = {chr![i],fast}.[i >= 0] -> {~chr?[i-1],v}.Fm[i-1];
 * Ori[i] = {fire,fr}.(Fm[i] || Fp[i]) + {chr?[i],fast};
 * SeedOris[i] = [i <= m] -> {seed, fast}.(Ori[i] || SeedOris[i+1]);
 *
 * SeedOris[0];
 *
 * Each kb position on the chromosome has an origin (Ori) that can
 * either fire at rate fr (spawning two replication forks) or be
 * passively replicated by an incoming fork (chr handshake).  Forks
 * move at rate v kb/min and terminate at chromosome boundaries or
 * when they encounter an already-replicated position.
 ***/

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <numeric>
#include "common.h"
#include "error_handling.h"
#include "htsInterface.h"
#include "iod.h"

#define CHR_LEN 2000 // Assumes reads of longer than 1 Mb are fleetingly rare and can be ignored for the purposes of IOD estimation


static const char *help=
"IOD: DNAscent executable that estimates inter-origin distance.\n"
"To run DNAscent IOD, do:\n"
"   DNAscent IOD -l /path/to/leftForks_DNAscent_forksense.bed -r /path/to/rightForks_DNAscent_forksense.bed -a /path/to/BrdU_DNAscent_forkSense.bed -d /path/to/detectOutput.bam -o /path/to/output.IOD \n"
"Required arguments are:\n"
"  -l,--left                 path to leftForks file from forkSense detect with `bed` extension,\n"
"  -r,--right                path to rightFork file from forkSense detect with `bed` extension,\n"
"     --origin               path to origin file from forkSense detect with `bed` extension,\n"
"     --termination          path to termination file from forkSense detect with `bed` extension,\n"
"  -d,--detect               path to output from detect with `detect` or `bam` extension,\n"
"  -o,--output               path to output file or directory for IOD.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


double wassersteinDistance(std::vector<int> a, std::vector<int> b) {

	std::sort(a.begin(), a.end());
	std::sort(b.begin(), b.end());

	size_t n = a.size();
	size_t m = b.size();

	double totalDist = 0.0;
	size_t i = 0, j = 0;
	double prevCDF_a = 0.0, prevCDF_b = 0.0;
	double prevVal = 0.0;
	bool first = true;

	while (i < n || j < m) {

		double val;
		if (i < n && (j >= m || a[i] <= b[j])) {
			val = a[i];
		} else {
			val = b[j];
		}

		if (!first) {
			totalDist += std::abs(prevCDF_a - prevCDF_b) * (val - prevVal);
		}
		first = false;

		while (i < n && a[i] == val) {
			prevCDF_a += 1.0 / n;
			i++;
		}
		while (j < m && b[j] == val) {
			prevCDF_b += 1.0 / m;
			j++;
		}
		prevVal = val;
	}

	return totalDist;
}


struct Fork {
	int position;
	int direction; // +1 (plus strand) or -1 (minus strand)
};

struct ForkCall {
	int EdU_start;
	int EdU_end;
	int BrdU_start;
	int BrdU_end;
	int direction; // +1 (plus strand) or -1 (minus strand)
};

struct MatchedForkPair {
	ForkCall leftFork;
	ForkCall rightFork;
};


struct SimulationResult {
	std::vector<double> leftForkTimes;
	std::vector<double> rightForkTimes;
	std::vector<int> firedOriginPositions;
	double completionTime;
};


struct Arguments {
	std::string lForkInput;
	std::string rForkInput;
	std::string originInput;
	std::string terminationInput;
	std::string DetectInput;
	std::string output;
    bool specifiedLeft = false;
    bool specifiedRight = false;
    bool specifiedOrigin = false;
    bool specifiedTermination = false;
    bool humanReadable = false;
	bool specifiedOutput = false;
    bool specifiedDetect = false;
	int threads = 1;
};


Arguments parseIODArguments( int argc, char** argv ) {

    if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent IOD." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}
 	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){
 		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}

	Arguments args;

	for ( int i = 1; i < argc; ){

 		std::string flag( argv[ i ] );

 		if ( flag == "-l" or flag == "--left" ){

 			if (i == argc-1) throw TrailingFlag(flag);

            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());

            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);

			args.lForkInput = strArg;
			i+=2;
			args.specifiedLeft = true;
		}
        else if ( flag == "-r" or flag == "--right" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.rForkInput = strArg;
			i+=2;
			args.specifiedRight = true;
        }
        else if ( flag == "--origin" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.originInput = strArg;
			i+=2;
			args.specifiedOrigin = true;
        }
		else if ( flag == "--termination" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.terminationInput = strArg;
			i+=2;
			args.specifiedTermination = true;
        }
        else if ( flag == "-h" or flag == "--help" ){
            std::cout << help << std::endl;
            exit(EXIT_SUCCESS);
		}
        else if ( flag == "-d" or flag == "--detect" ){
 		
 			if (i == argc-1) throw TrailingFlag(flag);		
 		
 			std::string strArg( argv[ i + 1 ] );
 			
 			const char *ext = get_ext(strArg.c_str());
				
			if (strcmp(ext,"bam") == 0){
				args.humanReadable = false;
			}
			else if (strcmp(ext,"detect") == 0){
				args.humanReadable = true;
            }
			else{
				throw InvalidExtension(ext);
			}
            args.specifiedDetect = true;
			
			args.DetectInput = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){
		
			if (i == argc-1) throw TrailingFlag(flag);		
		
 			std::string strArg( argv[ i + 1 ] );
            args.output = strArg;
			i+=2;
			args.specifiedOutput = true;
        }
		else if (flag == "-t" or flag == "--threads") {

			if (i == argc-1) throw TrailingFlag(flag);
			std::string strArg( argv[ i + 1 ] );
			args.threads = std::stoi(strArg);
			i += 2; // Skip the flag and its argument
		}
        else throw InvalidOption( flag );
	}
    if ( !args.specifiedOutput or !args.specifiedOrigin or !args.specifiedTermination or !args.specifiedDetect or (!args.specifiedLeft and !args.specifiedRight) ) {
        std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent IOD." << std::endl;
        exit(EXIT_FAILURE);
    }

    return args;
}


SimulationResult runGillespie(int m, double v, double fr, std::mt19937 &gen) {

	// State: each position is either replicated or unreplicated (has a live Ori)
	std::vector<bool> replicated(m + 1, false);
	std::vector<Fork> forks;
	std::vector<int> firedOrigins;
	int nUnreplicated = m + 1;

	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	double t = 0.0;

    std::vector<double> leftForkTimes(m + 1, -1.0);
    std::vector<double> rightForkTimes(m + 1, -1.0);

	while (nUnreplicated > 0 || !forks.empty()) {

		int nForks = static_cast<int>(forks.size());

		// Propensities
		double fireRate = nUnreplicated * fr;  // any Ori can fire
		double forkRate = nForks * v;           // any fork can advance 1 kb
		double totalRate = fireRate + forkRate;

		if (totalRate <= 0.0) break;

		// Time to next event (exponentially distributed)
		double u1 = uniform(gen);
		if (u1 == 0.0) u1 = std::numeric_limits<double>::min();
		double dt = -std::log(u1) / totalRate;
		t += dt;

		// Select event
		double r = uniform(gen) * totalRate;

		if (r < fireRate) {

			// --- Origin firing event ---
			// Select the k-th unreplicated position
			int k = static_cast<int>(r / fr);
			if (k >= nUnreplicated) k = nUnreplicated - 1;

			int count = 0;
			int firePos = -1;
			for (int i = 0; i <= m; i++) {
				if (!replicated[i]) {
					if (count == k) {
						firePos = i;
						break;
					}
					count++;
				}
			}

			// Ori[firePos] fires: position becomes replicated, spawn Fm || Fp
			replicated[firePos] = true;
			nUnreplicated--;
			leftForkTimes[firePos] = t;
			rightForkTimes[firePos] = t;
			firedOrigins.push_back(firePos);

			forks.push_back({firePos, +1}); // Fp
			forks.push_back({firePos, -1}); // Fm

		} else {

			// --- Fork movement event ---
			int k = static_cast<int>((r - fireRate) / v);
			if (k >= nForks) k = nForks - 1;

			int newPos = forks[k].position + forks[k].direction;

			if (newPos < 0 || newPos > m || replicated[newPos]) {
				// Fork termination: boundary or collision (no Ori to sync with)
				forks[k] = forks.back();
				forks.pop_back();
			} else {
				// chr![newPos] syncs with Ori[newPos]: passive replication
				replicated[newPos] = true;
				nUnreplicated--;
				forks[k].position = newPos;
                if (forks[k].direction == -1) {
                    leftForkTimes[newPos] = t;
                } else {
                    rightForkTimes[newPos] = t;
                }
			}
		}
	}

	return {leftForkTimes, rightForkTimes, firedOrigins, t};
}


std::pair<std::vector<ForkCall>, std::vector<ForkCall>> parseForkCalls(
	const SimulationResult &result,
	double pulseStart,
	double EdUpulseEnd,
	double pulseEnd,
	int readStart,
	int readEnd) {

	std::vector<ForkCall> leftForkCalls;
	std::vector<ForkCall> rightForkCalls;

	int n = std::min(readEnd + 1, static_cast<int>(result.leftForkTimes.size()));

	// Parse left fork tracks (direction = -1)
	// Contiguous segments where leftForkTimes >= 0 each correspond to one left fork
	int i = n - 1;
	while (readStart < i) {

		// Skip positions that weren't replicated by a left fork
		if (result.leftForkTimes[i] < 0) {
			i--;
			continue;
		}

		// Get contiguous segment of left fork replication (one fork track)
		int segStart = i;
		while (i >= readStart && result.leftForkTimes[i] >= 0) {
			i--;
		}
		int segEnd = i + 1;

		// Divide the segment into EdU and BrdU incorporation based on pulse lengths
		int eduStart = -1, eduEnd = -1, brduStart = -1, brduEnd = -1;
		for (int j = segStart; j >= segEnd; j--) {
			double t = result.leftForkTimes[j];
			if (t >= pulseStart && t < EdUpulseEnd) {
				if (eduStart == -1) eduStart = j;
				eduEnd = j;
			}
			if (t >= EdUpulseEnd && t < pulseEnd) {
				if (brduStart == -1) brduStart = j;
				brduEnd = j;
			}
		}

		// This length check ensures that this would actually be called as a fork by DNAscent
		int EdULength = (eduStart != -1 && eduEnd != -1) ? std::abs(eduEnd - eduStart + 1) : 0;
		int BrdULength = (brduStart != -1 && brduEnd != -1) ? std::abs(brduEnd - brduStart + 1) : 0;

		if (eduStart != -1 && brduStart != -1 && EdULength > 2 && BrdULength > 2) {
			ForkCall call;
			call.EdU_start = eduStart;
			call.EdU_end = eduEnd;
			call.BrdU_start = brduStart;
			call.BrdU_end = brduEnd;
			call.direction = -1;
			leftForkCalls.push_back(call);
		}
	}

	// Parse right fork tracks (direction = +1)
	n = std::min(readEnd + 1, static_cast<int>(result.rightForkTimes.size()));
	i = readStart;
	while (i < n) {

		// Skip positions that weren't replicated by a right fork
		if (result.rightForkTimes[i] < 0) {
			i++;
			continue;
		}

		// Get contiguous segment of right fork replication (one fork track)
		int segStart = i;
		while (i < n && result.rightForkTimes[i] >= 0) {
			i++;
		}
		int segEnd = i - 1;

		// Divide the segment into EdU and BrdU incorporation based on pulse lengths
		int eduStart = -1, eduEnd = -1, brduStart = -1, brduEnd = -1;
		for (int j = segStart; j <= segEnd; j++) {
			double t = result.rightForkTimes[j];
			if (t >= pulseStart && t < EdUpulseEnd) {
				if (eduStart == -1) eduStart = j;
				eduEnd = j;
			}
			if (t >= EdUpulseEnd && t < pulseEnd) {
				if (brduStart == -1) brduStart = j;
				brduEnd = j;
			}
		}

		// This length check ensures that this would actually be called as a fork by DNAscent
		int EdULength = (eduStart != -1 && eduEnd != -1) ? std::abs(eduEnd - eduStart + 1) : 0;
		int BrdULength = (brduStart != -1 && brduEnd != -1) ? std::abs(brduEnd - brduStart + 1) : 0;

		if (eduStart != -1 && brduStart != -1 && EdULength > 2 && BrdULength > 2) {
			ForkCall call;
			call.EdU_start = eduStart;
			call.EdU_end = eduEnd;
			call.BrdU_start = brduStart;
			call.BrdU_end = brduEnd;
			call.direction = +1;
			rightForkCalls.push_back(call);
		}
	}

	return {leftForkCalls, rightForkCalls};
}


std::vector<MatchedForkPair> matchForkCalls(
	const std::vector<ForkCall> &leftForkCalls,
	const std::vector<ForkCall> &rightForkCalls) {

	// Merge all calls into a single list
	std::vector<ForkCall> allCalls;
	allCalls.insert(allCalls.end(), leftForkCalls.begin(), leftForkCalls.end());
	allCalls.insert(allCalls.end(), rightForkCalls.begin(), rightForkCalls.end());

	// Sort by leftmost labelled position
	std::sort(allCalls.begin(), allCalls.end(), [](const ForkCall &a, const ForkCall &b) {
		int aLeft = std::min(a.EdU_start, a.BrdU_start);
		int bLeft = std::min(b.EdU_start, b.BrdU_start);
		return aLeft < bLeft;
	});

	// Match adjacent pairs where a left fork is immediately followed by a right fork
	std::vector<MatchedForkPair> matches;
	for (size_t i = 0; i + 1 < allCalls.size(); i++) {
		if (allCalls[i].direction == -1 && allCalls[i + 1].direction == +1) {
			matches.push_back({allCalls[i], allCalls[i + 1]});
		}
	}

	return matches;
}


std::vector<int> runPulseChase(SimulationResult &result, std::mt19937 &gen, std::vector<int> &readLengths) {

	double EdUpulseTime = 5.0; // Time of EdU pulse in minutes
	double BrdUpulseTime = 10.0; // Time of BrdU pulse in minutes
	
	if (result.completionTime <= EdUpulseTime + BrdUpulseTime) return {};
	std::uniform_real_distribution<double> uniform(0.0, result.completionTime - EdUpulseTime - BrdUpulseTime);
	double pulseStart = uniform(gen);
	double EdUpulseEnd = pulseStart + EdUpulseTime;
	double pulseEnd = EdUpulseEnd + BrdUpulseTime;
	

    std::uniform_int_distribution<> readDist(0, readLengths.size()-1);
    int readIndex = readDist(gen);
	int midpoint = CHR_LEN / 2;
	int sampledReadLength = readLengths[readIndex]; // in kb
	int readStart = midpoint - sampledReadLength / 2;
	int readEnd = midpoint + sampledReadLength / 2;
	if (readStart < 0 or readEnd >= CHR_LEN) {
		return {};
	}

	auto forkCalls = parseForkCalls(result, pulseStart, EdUpulseEnd, pulseEnd, readStart, readEnd);
	auto matchedPairs = matchForkCalls(forkCalls.first, forkCalls.second);

	std::vector<int> forkDistanceTravelled;
	for (const auto &pair : matchedPairs) {

		// Approximate the origin as the midpoint between EdU start sites, as DNAscent would do
		int originPos = pair.leftFork.EdU_start + (pair.rightFork.EdU_start - pair.leftFork.EdU_start)/2;
		int leftDistance = originPos - pair.leftFork.BrdU_end;
		int rightDistance = pair.rightFork.BrdU_end - originPos;

		forkDistanceTravelled.push_back(std::max(leftDistance, rightDistance));
	}

	return forkDistanceTravelled;
}


std::vector<double> parseForkBed(std::string fileInput, std::vector<std::string> &ignoreIDs) {
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
	std::vector<double> forkLengths;
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::vector< std::string > columns = split(line);               
            std::string readID = columns[3];

            if ( std::find(ignoreIDs.begin(), ignoreIDs.end(), readID) == ignoreIDs.end() ) {

                int pulse5Prime = std::stoi(columns[1]);
                int pulse3Prime = std::stoi(columns[2]);
                int read5Prime = std::stoi(columns[4]);
                int read3Prime = std::stoi(columns[5]);

				// ignore forks near the end of the read
				if (std::abs(pulse5Prime - read5Prime) < 3000) continue; 
				if (std::abs(read3Prime - pulse3Prime) < 3000) continue; 
				if (std::abs(pulse5Prime - read3Prime) < 3000) continue; 
				if (std::abs(pulse3Prime - read5Prime) < 3000) continue; 

				forkLengths.push_back( (pulse3Prime - pulse5Prime)/1000.); // in kb
			}
		}
	}

    file.close();
    return forkLengths;
}


void getIgnoreIDs(std::string fileInput, std::vector<std::string> &ignoreIDs) {
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::vector< std::string > columns = split(line);               
            std::string readID = columns[3];

			 if ( std::find(ignoreIDs.begin(), ignoreIDs.end(), readID) == ignoreIDs.end() ) {

				 ignoreIDs.push_back(readID);
			}
		}
	}
    file.close();
}


std::vector<int> calcOriginForkDistances(
	const std::string &originFile,
	const std::string &leftForkFile,
	const std::string &rightForkFile) {

	// Parse origin BED: readID -> origin midpoint
	std::map<std::string, int> originMidpoints;
	{
		std::ifstream file(originFile);
		if (!file.is_open()) {
			std::cerr << "Error: Could not open file " << originFile << "\n";
			exit(EXIT_FAILURE);
		}
		std::string line;
		while (std::getline(file, line)) {
			if (!line.empty() && line[0] != '#') {
				std::vector<std::string> columns = split(line);
				std::string readID = columns[3];
				int start = std::stoi(columns[1]);
				int end = std::stoi(columns[2]);
				originMidpoints[readID] = (start + end) / 2;
			}
		}
	}

	// Parse left fork BED: readID -> fork tips (start coordinate, leftward-moving)
	std::map<std::string, std::vector<int>> leftForkTips;
	{
		std::ifstream file(leftForkFile);
		if (!file.is_open()) {
			std::cerr << "Error: Could not open file " << leftForkFile << "\n";
			exit(EXIT_FAILURE);
		}
		std::string line;
		while (std::getline(file, line)) {
			if (!line.empty() && line[0] != '#') {
				std::vector<std::string> columns = split(line);
				std::string readID = columns[3];
				int tip = std::stoi(columns[1]); // start = leading edge for left fork
				leftForkTips[readID].push_back(tip);
			}
		}
	}

	// Parse right fork BED: readID -> fork tips (end coordinate, rightward-moving)
	std::map<std::string, std::vector<int>> rightForkTips;
	{
		std::ifstream file(rightForkFile);
		if (!file.is_open()) {
			std::cerr << "Error: Could not open file " << rightForkFile << "\n";
			exit(EXIT_FAILURE);
		}
		std::string line;
		while (std::getline(file, line)) {
			if (!line.empty() && line[0] != '#') {
				std::vector<std::string> columns = split(line);
				std::string readID = columns[3];
				int tip = std::stoi(columns[2]); // end = leading edge for right fork
				rightForkTips[readID].push_back(tip);
			}
		}
	}

	// For each origin, calculate distance to the nearest fork tip
	std::vector<int> distances;
	for (const auto &origin : originMidpoints) {
		const std::string &readID = origin.first;
		int midpoint = origin.second;

		auto leftIt = leftForkTips.find(readID);
		int minLeftDist = -1;
		if (leftIt != leftForkTips.end()) {
			minLeftDist = std::abs(midpoint - leftIt->second[0]);
			for (size_t j = 1; j < leftIt->second.size(); j++) {
				int d = std::abs(midpoint - leftIt->second[j]);
				if (d < minLeftDist) minLeftDist = d;
			}
		}

		auto rightIt = rightForkTips.find(readID);
		int minRightDist = -1;
		if (rightIt != rightForkTips.end()) {
			minRightDist = std::abs(rightIt->second[0] - midpoint);
			for (size_t j = 1; j < rightIt->second.size(); j++) {
				int d = std::abs(rightIt->second[j] - midpoint);
				if (d < minRightDist) minRightDist = d;
			}
		}

		// Take the tip of the fork that was furthest away from the midpoint
		distances.push_back((std::max(minLeftDist, minRightDist))/1000.); // convert to kb
	}

	return distances;
}


void detectReadlength(std::string filename, std::vector<int> &readLengths) {

    std::ifstream detectFile(filename);
    std::string line;
    int prog = 0;

	while( std::getline( detectFile, line ) ){

        if (line.substr(0,1) == "#" or line.length() == 0) continue; //ignore header and blank lines
        if ( line.substr(0,1) == ">" ){

            std::vector< std::string > columns = split(line);
            assert(columns.size() == 5); //well-formed detect header

            int refStart = std::stoi(columns[2]);
            int refEnd = std::stoi(columns[3]);
            
            readLengths.push_back((refEnd - refStart)/1000); // convert to kb

            prog ++;
            if (prog % 1000 == 0) std::cout << "\rProcessed " << prog << " reads..." << std::flush;	
        }
    }
	std::cout << " Done." << std::endl;
    detectFile.close();
}


void bamReadlength(std::string filename, std::vector<int> &readLengths) {

    htsFile *bam_fh = sam_open(filename.c_str(), "r");
    if (bam_fh == NULL) throw IOerror(filename);

    bam_hdr_t *bam_hdr = sam_hdr_read(bam_fh);
    bam1_t *itr_record = bam_init1();

    int prog = 0;

    while(sam_read1(bam_fh, bam_hdr, itr_record) >= 0){
          
        int refStart,refEnd;
        getRefEnd(itr_record,refStart,refEnd);
        readLengths.push_back((refEnd - refStart)/1000); // convert to kb

        prog ++;
        if (prog % 1000 == 0) std::cout << "\rProcessed " << prog << " reads..." << std::flush;	
    }
	std::cout << " Done." << std::endl;
	bam_destroy1(itr_record);				
	bam_hdr_destroy(bam_hdr);
	hts_close(bam_fh);
}


int iod_main(int argc, char** argv) {

	Arguments args = parseIODArguments(argc, argv);

	// Get readIDs to ignore (those with an origin or termination call) so that we get an accurate measure of fork speed
	std::vector<std::string> ignoreIDs;
	getIgnoreIDs(args.originInput, ignoreIDs);
	getIgnoreIDs(args.terminationInput, ignoreIDs);

	// Calculate median fork speed from bed files, ignoring any forks near ends of reads or reads with an origin or termination call
	std::vector<double> forkLengths = parseForkBed(args.lForkInput, ignoreIDs);
	std::vector<double> rightForkLengths = parseForkBed(args.rForkInput, ignoreIDs);
	forkLengths.insert(forkLengths.end(), rightForkLengths.begin(), rightForkLengths.end());
	double avgForkSpeed = vectorMean(forkLengths) / 15.; // convert to kb/min
	std::cout << "Mean fork speed = " << std::fixed << std::setprecision(3) << avgForkSpeed << " kb/min" << std::endl;

	// Calculate distance between fork tip and DNAscent origin call
	std::vector<int> data_originForkDistances = calcOriginForkDistances(args.originInput, args.lForkInput, args.rForkInput);

    // Parse detect output
    std::vector<int> readLengths;

    if (args.humanReadable) {// is .detect format

        detectReadlength(args.DetectInput, readLengths);
    }
    else{// is modbam
        
        bamReadlength(args.DetectInput, readLengths);       
    }

	int nSims = 100000; 
	unsigned int seed = 42;

	// Objective function: run simulations at a given fr and return Wasserstein distance to data
	// Also computes the median IOD for that firing rate
	auto objective = [&](double fr, double &medianIOD) -> double {
		std::vector<int> sim_originForkDistances;
		std::vector<int> sim_IODs;

		#pragma omp parallel for num_threads(args.threads) schedule(dynamic)
		for (int sim = 0; sim < nSims; sim++) {

			std::mt19937 gen(seed + sim);
			SimulationResult result = runGillespie(CHR_LEN, avgForkSpeed, fr, gen);
			std::vector<int> OriginForkDistances;
			std::vector<int> iods;
			
			for (size_t i = 0; i < 100; i++){
			
				std::vector<int> distances = runPulseChase(result, gen, readLengths);
				OriginForkDistances.insert(OriginForkDistances.end(), distances.begin(), distances.end());

				// Compute IODs from fired origin positions
				std::vector<int> origins = result.firedOriginPositions;
				std::sort(origins.begin(), origins.end());
				for (size_t j = 1; j < origins.size(); j++) {
					iods.push_back(origins[j] - origins[j - 1]);
				}
			}

			#pragma omp critical
			{
				sim_originForkDistances.insert(sim_originForkDistances.end(), OriginForkDistances.begin(), OriginForkDistances.end());
				sim_IODs.insert(sim_IODs.end(), iods.begin(), iods.end());
			}
		}
		medianIOD = vectorMedian(sim_IODs);
		return wassersteinDistance(sim_originForkDistances, data_originForkDistances);
	};

	// Golden section search to find fr that minimises Wasserstein distance
	double frLow = 1e-7;
	double frHigh = 1e-3;
	double tol = 1e-7;
	const double invPhi = (std::sqrt(5.0) - 1.0) / 2.0; // ~0.618

	double x1 = frHigh - invPhi * (frHigh - frLow);
	double x2 = frLow + invPhi * (frHigh - frLow);
	double iod1, iod2;
	double f1 = objective(x1, iod1);
	double f2 = objective(x2, iod2);

	std::cerr << "Starting parameter search for optimal firing rate in [" << frLow << ", " << frHigh << "]..." << std::endl;

	int iteration = 1;
	while ((frHigh - frLow) > tol) {

		if (f1 < f2) {
			frHigh = x2;
			x2 = x1;
			f2 = f1;
			iod2 = iod1;
			x1 = frHigh - invPhi * (frHigh - frLow);
			f1 = objective(x1, iod1);
		} else {
			frLow = x1;
			x1 = x2;
			f1 = f2;
			iod1 = iod2;
			x2 = frLow + invPhi * (frHigh - frLow);
			f2 = objective(x2, iod2);
		}

		double bestFr = (x1 + x2) / 2.0;  // or use the point with better objective value
		double bestIOD = (f1 < f2) ? iod1 : iod2;
		std::cerr << "[Iteration " << std::setw(2) << iteration << "] "
				  << "fr = " << std::scientific << std::setprecision(3) << std::setw(9) << bestFr
				  << "  bracket = [" << std::setw(9) << frLow << ", " << std::setw(9) << frHigh << "]"
				  << "  IOD = " << std::fixed << std::setw(3) << std::setprecision(0) << bestIOD
				  << "  W = " << std::fixed << std::setprecision(3) << std::setw(5) << std::min(f1, f2) << std::endl;
		iteration++;
	}

	double bestFr = (frLow + frHigh) / 2.0;
	double bestW = std::min(f1, f2);
	double bestIOD = (f1 < f2) ? iod1 : iod2;

	// Estimate uncertainty via local curvature of the objective
	// Evaluate at three points around the optimum to get curvature in fr-space,
	// then convert uncertainty to IOD-space using the IOD values at those points
	double h = bestFr * 0.05; // 5% step for finite differences
	double iodLeft, iodCentre, iodRight;
	double fLeft = objective(bestFr - h, iodLeft);
	double fCentre = objective(bestFr, iodCentre);
	double fRight = objective(bestFr + h, iodRight);
	double curvature = (fLeft - 2.0 * fCentre + fRight) / (h * h);


	std::cerr << "Parameter search complete." << std::endl << std::endl;
	std::cerr << "Summary:" << std::endl;
	std::cerr << "--------- " << std::endl;
	std::cerr << "Optimal fr = " << std::scientific << std::setprecision(3) << bestFr << std::endl;
	std::cerr << "Wasserstein distance = " << std::fixed << std::setprecision(3) << bestW << std::endl;
	std::cerr << "Median IOD = " << std::fixed << std::setprecision(0) << bestIOD << " kb" << std::endl;

	if (curvature > 0) {
		// Uncertainty in fr-space
		double frUncertainty = 1.96 / std::sqrt(curvature);

		// Convert fr uncertainty to IOD uncertainty using dIOD/dfr estimated from the three evaluation points
		double dIOD_dfr = (iodRight - iodLeft) / (2.0 * h);
		double iodUncertainty = std::abs(dIOD_dfr) * frUncertainty;

		std::cerr << "Curvature at minimum: " << std::scientific << std::setprecision(3) << curvature << std::endl;
		std::cerr << "95% confidence interval: [" << std::fixed << std::setprecision(0) << bestIOD - iodUncertainty << " kb , " << std::fixed << std::setprecision(0) << bestIOD + iodUncertainty << " kb]" << std::endl;
	} else {
		std::cerr << "Warning: curvature is non-positive (" << std::scientific << std::setprecision(3) << curvature << "), uncertainty could not be estimated." << std::endl;
	}

	return 0;
}
