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
#include <sstream>
#include <ctime>
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include "htsInterface.h"
#include "iod.h"

#define CHR_LEN 1000 // Assumes reads of longer than 1 Mb are fleetingly rare and can be ignored for the purposes of IOD estimation
double MODEL_RES = 10.; // 100 bp resolution

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
"     --tPulse1              duration (in minutes) of the first analogue pulse,\n"
"     --tPulse2              duration (in minutes) of the second analogue pulse,\n"
"  -o,--output               path to output file or directory for IOD.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


double wasserstein(std::vector<double> a, std::vector<double> b) {

	std::sort(a.begin(), a.end());
	std::sort(b.begin(), b.end());

	size_t n = a.size();
	size_t m = b.size();

	// 1-Wasserstein distance: integral of |F_a(x) - F_b(x)| dx
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
			double diff = prevCDF_a - prevCDF_b;
			totalDist += std::abs(diff) * (val - prevVal);
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


double wassersteinPresorted(const std::vector<double> &a, const std::vector<double> &b) {

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
			double diff = prevCDF_a - prevCDF_b;
			totalDist += std::abs(diff) * (val - prevVal);
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
	std::vector<int> terminationPositions;
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
	bool specifiedPulse1 = false;
	bool specifiedPulse2 = false;
	int threads = 1;
	double EdUpulseLength; // in minutes
	double BrdUpulseLength; // in minutes
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
		else if ( flag == "--tPulse1" ){

			if (i == argc-1) throw TrailingFlag(flag);
			std::string strArg( argv[ i + 1 ] );
			args.EdUpulseLength = std::stod(strArg);
			i += 2; // Skip the flag and its argument
            args.specifiedPulse1 = true;
        }
        else if ( flag == "--tPulse2" ){
			if (i == argc-1) throw TrailingFlag(flag);
			std::string strArg( argv[ i + 1 ] );
			args.BrdUpulseLength = std::stod(strArg);
			i += 2; // Skip the flag and its argument
            args.specifiedPulse2 = true;
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
    if ( !args.specifiedOutput or !args.specifiedPulse1 or !args.specifiedPulse2 or !args.specifiedOrigin or !args.specifiedTermination or !args.specifiedDetect or (!args.specifiedLeft and !args.specifiedRight) ) {
        std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent IOD." << std::endl;
        exit(EXIT_FAILURE);
    }

    return args;
}


SimulationResult runGillespie(int m, double v, double fr, std::mt19937 &gen) {

	//map to model resolution
	m = m * MODEL_RES;
	v = v * MODEL_RES;
	fr = fr / MODEL_RES;

	// State: each position is either replicated or unreplicated (has a live Ori)
	std::vector<bool> replicated(m + 1, false);
	std::vector<Fork> forks;
	std::vector<int> firedOrigins;
	std::vector<int> terminationPositions;
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
				// Fork termination: boundary or collision
				forks[k] = forks.back();
				forks.pop_back();
				if (newPos >= 0 && newPos <= m && replicated[newPos]) {
					// Record termination position
					terminationPositions.push_back(newPos);
				}
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

	return {leftForkTimes, rightForkTimes, firedOrigins, terminationPositions, t};
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
	int BrdUReqLength = 3*MODEL_RES; // minimum length of BrdU track to call a fork (ensures that this would likely be called as a fork by DNAscent's read end QC)
	int EdUReqLength = 3*MODEL_RES; // minimum length of EdU track to call a fork (ensures that this would likely be called as a fork by DNAscent's read end

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

		// check that BrdU track extends at least 2 kb from read start, otherwise this would likely be filtered out by DNAscent's read end QC
		if (eduStart != -1 && brduStart != -1 && EdULength >= EdUReqLength && BrdULength >= BrdUReqLength && std::abs(eduStart - readEnd) >= 2*MODEL_RES) { 
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

		if (eduStart != -1 && brduStart != -1 && EdULength >= EdUReqLength && BrdULength >= BrdUReqLength && std::abs(eduStart - readStart) >= 2*MODEL_RES) {
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


std::pair<std::vector<double>, std::vector<double>> runPulseChase(SimulationResult &result, std::mt19937 &gen, std::vector<int> &readLengths, double EdUpulseTime, double BrdUpulseTime) {

	if (result.completionTime <= EdUpulseTime + BrdUpulseTime) return {};
	std::uniform_real_distribution<double> uniform(0.0, result.completionTime - EdUpulseTime - BrdUpulseTime);
	double pulseStart = uniform(gen);
	double EdUpulseEnd = pulseStart + EdUpulseTime;
	double pulseEnd = EdUpulseEnd + BrdUpulseTime;
	

    std::uniform_int_distribution<> readDist(0, readLengths.size()-1);
    int readIndex = readDist(gen);
	int midpoint = CHR_LEN / 2 * MODEL_RES;
	int sampledReadLength = readLengths[readIndex] * MODEL_RES; // in kb
	int readStart = midpoint - sampledReadLength / 2;
	int readEnd = midpoint + sampledReadLength / 2;
	if (readStart < 0 or readEnd >= CHR_LEN * MODEL_RES) {
		return {};
	}

	auto forkCalls = parseForkCalls(result, pulseStart, EdUpulseEnd, pulseEnd, readStart, readEnd);

	std::vector<ForkCall> &leftForks = forkCalls.first;
	std::vector<ForkCall> &rightForks = forkCalls.second;
	int totalForks = leftForks.size() + rightForks.size();

	if (totalForks == 0) return {};

	std::vector<double> singleForkDists; // behind-distances from single-fork reads
	std::vector<double> forkPairDists;   // EdU start-to-start distances from multi-fork reads

	if (totalForks == 1) {

		// Single-fork read: behind-distance from EdU start to read end
		const ForkCall &fork = leftForks.empty() ? rightForks[0] : leftForks[0];
		double dist;
		if (fork.direction == +1) {
			dist = static_cast<double>( (fork.EdU_start - readStart) / MODEL_RES );
		} else {
			dist = static_cast<double>( (readEnd - fork.EdU_start) / MODEL_RES );
		}
		singleForkDists.push_back(dist);
	} else {
		// Multi-fork read: match right-moving forks with left-moving forks from the same origin
		// Sort right forks by EdU_start position (ascending)
		std::sort(rightForks.begin(), rightForks.end(), [](const ForkCall &a, const ForkCall &b) { return a.EdU_start < b.EdU_start; });
		// Sort left forks by EdU_start position (ascending)
		std::sort(leftForks.begin(), leftForks.end(), [](const ForkCall &a, const ForkCall &b) { return a.EdU_start < b.EdU_start; });

		// For each right fork, find the nearest left fork whose EdU_start is to the left (diverging pair from same origin)
		std::vector<bool> leftUsed(leftForks.size(), false);
		for (const auto &rf : rightForks) {
			int bestIdx = -1;
			int bestGap = std::numeric_limits<int>::max();
			for (size_t li = 0; li < leftForks.size(); li++) {
				if (leftUsed[li]) continue;
				int gap = rf.EdU_start - leftForks[li].EdU_start;
				if (gap >= 0 && gap < bestGap) {
					bestGap = gap;
					bestIdx = static_cast<int>(li);
				}
			}
			if (bestIdx >= 0) {
				leftUsed[bestIdx] = true;
				double d = static_cast<double>(bestGap) / MODEL_RES; // convert to kb
				forkPairDists.push_back(d);
			}
		}
	}

	return {singleForkDists, forkPairDists};
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

				assert(read5Prime <= pulse3Prime && pulse3Prime <= read3Prime);
				assert(read5Prime <= pulse5Prime && pulse5Prime <= read3Prime);

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


std::pair<std::vector<double>, std::vector<double>> calcForkBehindDistances(
	const std::string &leftForkFile,
	const std::string &rightForkFile) {

	struct ForkInfo {
		int eduStart;   // 5' end of the pulse (start of EdU track, origin side)
		int direction;  // +1 (rightward) or -1 (leftward)
		int readStart;  // 5' read boundary
		int readEnd;    // 3' read boundary
	};

	// Group forks by readID
	std::map<std::string, std::vector<ForkInfo>> readForks;

	// Parse left fork BED: eduStart is pulse3Prime (rightmost, the origin-side/trailing edge for a leftward-moving fork)
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
				int pulse5Prime = std::stoi(columns[1]);
				int pulse3Prime = std::stoi(columns[2]);
				int read5Prime = std::stoi(columns[4]);
				int read3Prime = std::stoi(columns[5]);
				if (std::abs(pulse3Prime - read3Prime) < 2000) continue;
				readForks[readID].push_back({pulse3Prime, -1, read5Prime, read3Prime});
			}
		}
	}

	// Parse right fork BED: eduStart is pulse5Prime (leftmost, the origin-side/trailing edge for a rightward-moving fork)
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
				int pulse5Prime = std::stoi(columns[1]);
				int pulse3Prime = std::stoi(columns[2]);
				int read5Prime = std::stoi(columns[4]);
				int read3Prime = std::stoi(columns[5]);
				if (std::abs(pulse5Prime - read5Prime) < 2000) continue;
				readForks[readID].push_back({pulse5Prime, +1, read5Prime, read3Prime});
			}
		}
	}

	// Compute distances: single-fork reads get behind-distance, multi-fork reads get matched fork-pair EdU start distances
	std::vector<double> singleForkDists; // behind-distances from single-fork reads
	std::vector<double> forkPairDists;   // EdU start-to-start distances from multi-fork reads
	for (const auto &entry : readForks) {
		const std::vector<ForkInfo> &forks = entry.second;
		if (forks.empty()) continue;

		if (forks.size() == 1) {
			
			// Single-fork read: distance from EdU start to read end behind the fork
			const ForkInfo &fork = forks[0];
			if (fork.direction == +1) {
				singleForkDists.push_back(static_cast<double>(fork.eduStart - fork.readStart) / 1000.0);
			} else {
				singleForkDists.push_back(static_cast<double>(fork.readEnd - fork.eduStart) / 1000.0);
			}

		} else {
			// Multi-fork read: match right-moving forks with left-moving forks from the same origin
			std::vector<const ForkInfo*> rightForks, leftForks;
			for (const auto &f : forks) {
				if (f.direction == +1) rightForks.push_back(&f);
				else leftForks.push_back(&f);
			}
			// Sort by eduStart position
			std::sort(rightForks.begin(), rightForks.end(), [](const ForkInfo *a, const ForkInfo *b) { return a->eduStart < b->eduStart; });
			std::sort(leftForks.begin(), leftForks.end(), [](const ForkInfo *a, const ForkInfo *b) { return a->eduStart < b->eduStart; });

			// For each right fork, find the nearest left fork whose EdU start is to the left (diverging pair)
			std::vector<bool> leftUsed(leftForks.size(), false);
			for (const auto *rf : rightForks) {
				int bestIdx = -1;
				int bestGap = std::numeric_limits<int>::max();
				for (size_t li = 0; li < leftForks.size(); li++) {
					if (leftUsed[li]) continue;
					int gap = rf->eduStart - leftForks[li]->eduStart;
					if (gap >= 0 && gap < bestGap) {
						bestGap = gap;
						bestIdx = static_cast<int>(li);
					}
				}
				if (bestIdx >= 0) {
					leftUsed[bestIdx] = true;
					forkPairDists.push_back(static_cast<double>(bestGap) / 1000.0); // convert to kb
				}
			}
		}
	}

	return {singleForkDists, forkPairDists};
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
	double avgForkSpeed = vectorMean(forkLengths) / (args.EdUpulseLength + args.BrdUpulseLength); // convert to kb/min
	std::cout << "Mean fork speed = " << std::fixed << std::setprecision(3) << avgForkSpeed << " kb/min" << std::endl;

	// Calculate distances: single-fork behind-distances and multi-fork pair distances
	auto dataDists = calcForkBehindDistances(args.lForkInput, args.rForkInput);
	std::vector<double> data_singleForkDists = dataDists.first;
	std::vector<double> data_forkPairDists = dataDists.second;
	std::cout << "Single-fork reads: " << data_singleForkDists.size() << " distances" << std::endl;
	std::cout << "Multi-fork reads: " << data_forkPairDists.size() << " fork-pair distances" << std::endl;

	// Precompute means for normalising the two Wasserstein components
	double meanSingleFork = vectorMean(data_singleForkDists);
	double meanForkPair = data_forkPairDists.size() > 0 ? vectorMean(data_forkPairDists) : 1.0;
	double dataFracSingle = static_cast<double>(data_singleForkDists.size()) / (data_singleForkDists.size() + data_forkPairDists.size());

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
	std::pair<std::vector<double>, std::vector<double>> *saveSimDist = nullptr;
	std::vector<std::tuple<double, double, double>> landscapePoints; // (fr, IOD, W)
	const size_t MAX_LANDSCAPE_SAMPLES = 50000;

	struct LandscapeEntry {
		std::vector<double> sortedSingleFork;
		std::vector<double> sortedForkPair;
		double medianIOD;
		double fracSingle;
	};
	std::map<double, LandscapeEntry> landscapeSimDists; // fr -> stored sim distributions + median IOD

	// Objective function: run simulations at a given fr and return combined Wasserstein distance to data
	// Uses the sum of 1D Wasserstein distances for single-fork and fork-pair distributions
	auto objective = [&](double fr, double &medianIOD) -> double {
		std::vector<double> sim_singleForkDists;
		std::vector<double> sim_forkPairDists;
		std::vector<int> sim_IODs;

		#pragma omp parallel for num_threads(args.threads) schedule(dynamic)
		for (int sim = 0; sim < nSims; sim++) {

			std::mt19937 gen(seed + sim);
			SimulationResult result = runGillespie(CHR_LEN, avgForkSpeed, fr, gen);
			std::vector<double> singleForkDists;
			std::vector<double> forkPairDists;
			std::vector<int> iods;
			
			for (size_t i = 0; i < 100; i++){
			
				auto distances = runPulseChase(result, gen, readLengths, args.EdUpulseLength, args.BrdUpulseLength);
				singleForkDists.insert(singleForkDists.end(), distances.first.begin(), distances.first.end());
				forkPairDists.insert(forkPairDists.end(), distances.second.begin(), distances.second.end());

				// Compute IODs from fired origin positions
				std::vector<int> origins = result.firedOriginPositions;
				std::sort(origins.begin(), origins.end());
				for (size_t j = 1; j < origins.size(); j++) {
					iods.push_back( (origins[j] - origins[j - 1]) / MODEL_RES ); // convert to kb
				}
			}

			#pragma omp critical
			{
				sim_singleForkDists.insert(sim_singleForkDists.end(), singleForkDists.begin(), singleForkDists.end());
				sim_forkPairDists.insert(sim_forkPairDists.end(), forkPairDists.begin(), forkPairDists.end());
				sim_IODs.insert(sim_IODs.end(), iods.begin(), iods.end());
			}
		}
		medianIOD = vectorMedian(sim_IODs);

		// Store subsampled, sorted sim distributions for bootstrap CI
		{
			auto subsample = [](const std::vector<double> &src, size_t maxN) -> std::vector<double> {
				std::vector<double> out;
				if (src.size() > maxN) {
					out.resize(maxN);
					std::mt19937 subGen(42);
					std::uniform_int_distribution<size_t> idx(0, src.size() - 1);
					for (size_t s = 0; s < maxN; s++) out[s] = src[idx(subGen)];
				} else {
					out = src;
				}
				std::sort(out.begin(), out.end());
				return out;
			};

			LandscapeEntry entry;
			entry.sortedSingleFork = subsample(sim_singleForkDists, MAX_LANDSCAPE_SAMPLES);
			entry.sortedForkPair = subsample(sim_forkPairDists, MAX_LANDSCAPE_SAMPLES);
			entry.medianIOD = medianIOD;
			double nSimTotal = sim_singleForkDists.size() + sim_forkPairDists.size();
			entry.fracSingle = (nSimTotal > 0) ? sim_singleForkDists.size() / nSimTotal : 0.0;
			landscapeSimDists[fr] = std::move(entry);
		}

		if (saveSimDist) *saveSimDist = {sim_singleForkDists, sim_forkPairDists};

		// Combined objective: normalised Wasserstein distances + single-fork fraction penalty
		double w1 = (!sim_singleForkDists.empty() && !data_singleForkDists.empty()) ? wasserstein(sim_singleForkDists, data_singleForkDists) : 0.0;
		double w2 = (!sim_forkPairDists.empty() && !data_forkPairDists.empty()) ? wasserstein(sim_forkPairDists, data_forkPairDists) : 0.0;
		double simTotal = sim_singleForkDists.size() + sim_forkPairDists.size();
		double simFracSingle = (simTotal > 0) ? sim_singleForkDists.size() / simTotal : 0.0;
		double fracPenalty = std::abs(simFracSingle - dataFracSingle) / dataFracSingle;
		double w = w1 / meanSingleFork + w2 / meanForkPair + fracPenalty;
		landscapePoints.push_back(std::make_tuple(fr, medianIOD, w));
		return w;
	};

	// Golden section search to find fr that minimises Wasserstein distance
	// Search in log10(fr) space so that probe points are evenly distributed
	// across orders of magnitude
	double frLow = 1e-7;
	double frHigh = 1e-3;
	double logFrLow = std::log10(frLow);
	double logFrHigh = std::log10(frHigh);
	double tol = 1e-7;
	double logTol = std::log10((frLow + tol) / frLow); // tolerance in log-space

	const double invPhi = (std::sqrt(5.0) - 1.0) / 2.0; // ~0.618

	double lx1 = logFrHigh - invPhi * (logFrHigh - logFrLow);
	double lx2 = logFrLow + invPhi * (logFrHigh - logFrLow);
	double iod1, iod2;
	double f1 = objective(std::pow(10.0, lx1), iod1);
	double f2 = objective(std::pow(10.0, lx2), iod2);

	std::cerr << "Starting parameter search for optimal firing rate in [" << frLow << ", " << frHigh << "]..." << std::endl;

	int iteration = 1;
	while ((logFrHigh - logFrLow) > logTol) {

		if (f1 < f2) {
			logFrHigh = lx2;
			lx2 = lx1;
			f2 = f1;
			iod2 = iod1;
			lx1 = logFrHigh - invPhi * (logFrHigh - logFrLow);
			f1 = objective(std::pow(10.0, lx1), iod1);
		} else {
			logFrLow = lx1;
			lx1 = lx2;
			f1 = f2;
			iod1 = iod2;
			lx2 = logFrLow + invPhi * (logFrHigh - logFrLow);
			f2 = objective(std::pow(10.0, lx2), iod2);
		}

		double bestFr = std::pow(10.0, (lx1 + lx2) / 2.0);
		double bestIOD = (f1 < f2) ? iod1 : iod2;
		std::cerr << "[Iteration " << std::setw(2) << iteration << "] "
				  << "fr = " << std::scientific << std::setprecision(3) << std::setw(9) << bestFr
				  << "  bracket = [" << std::setw(9) << std::pow(10.0, logFrLow) << ", " << std::setw(9) << std::pow(10.0, logFrHigh) << "]"
				  << "  IOD = " << std::fixed << std::setw(3) << std::setprecision(0) << bestIOD
				  << "  W = " << std::fixed << std::setprecision(3) << std::setw(5) << std::min(f1, f2) << std::endl;
		iteration++;
	}

	double bestFr = std::pow(10.0, (logFrLow + logFrHigh) / 2.0);
	double bestW = std::min(f1, f2);
	double bestIOD = (f1 < f2) ? iod1 : iod2;

	// Evaluate at the optimum to get the simulation distribution for output
	std::pair<std::vector<double>, std::vector<double>> bestSimDist;
	saveSimDist = &bestSimDist;
	double iodAtBest;
	double wAtBest = objective(bestFr, iodAtBest);
	saveSimDist = nullptr;
	bestIOD = iodAtBest;
	bestW = wAtBest;

	// The golden section search clusters most evaluation points near the
	// optimum, leaving large gaps elsewhere.  The bootstrap needs coverage
	// across the full firing-rate range so that resampled optima can be
	// located accurately rather than snapping to sparse points.
	{
		std::cerr << "Building landscape for confidence interval estimation..." << std::endl;
		int savedNSims = nSims;
		nSims = 20000;

		double logBestFr = std::log10(bestFr);
		double logFrMin  = std::log10(1e-5);
		double logFrMax  = std::log10(1e-2);

		// Broad grid across full range
		const int nBroad = 10;
		for (int i = 0; i < nBroad; i++) {
			double logFr = logFrMin + i * (logFrMax - logFrMin) / (nBroad - 1);
			bool exists = false;
			for (const auto &lp : landscapeSimDists) {
				if (std::abs(std::log10(lp.first) - logFr) < 0.05) { exists = true; break; }
			}
			if (!exists) {
				double iod;
				objective(std::pow(10.0, logFr), iod);
			}
		}

		// Fine grid: +/- 1 decade around optimum
		const int nFine = 20;
		double fineLogLow  = std::max(logFrMin,  logBestFr - 1.0);
		double fineLogHigh = std::min(logFrMax, logBestFr + 1.0);
		for (int i = 0; i < nFine; i++) {
			double logFr = fineLogLow + i * (fineLogHigh - fineLogLow) / (nFine - 1);
			bool exists = false;
			for (const auto &lp : landscapeSimDists) {
				if (std::abs(std::log10(lp.first) - logFr) < 0.02) { exists = true; break; }
			}
			if (!exists) {
				double iod;
				objective(std::pow(10.0, logFr), iod);
			}
		}

		nSims = savedNSims;
		std::cerr << "Landscape has " << landscapeSimDists.size() << " points for bootstrap CI." << std::endl;
	}

	// Bootstrap CI: resample the observed data and find the best landscape point for each resample
	const int nBootstrap = 10000;
	std::vector<double> bootstrapIODs;
	bootstrapIODs.reserve(nBootstrap);

	std::mt19937 bootGen(42);

	// Helper to compute smoothed bootstrap parameters for a distribution
	struct SmoothParams {
		double mean;
		double smoothBW;
		double shrinkFactor;
	};
	auto computeSmoothParams = [](const std::vector<double> &data) -> SmoothParams {
		SmoothParams p;
		p.mean = 0.0;
		for (const auto &v : data) p.mean += v;
		p.mean /= data.size();
		double sd = 0.0;
		for (const auto &v : data) sd += (v - p.mean) * (v - p.mean);
		sd = std::sqrt(sd / (data.size() - 1));
		p.smoothBW = sd * std::pow(static_cast<double>(data.size()), -0.2) * 0.5;
		p.shrinkFactor = 1.0 / std::sqrt(1.0 + (p.smoothBW * p.smoothBW) / (sd * sd));
		return p;
	};

	SmoothParams sp1 = computeSmoothParams(data_singleForkDists);
	SmoothParams sp2 = data_forkPairDists.size() > 1 ? computeSmoothParams(data_forkPairDists) : SmoothParams{0, 0, 1};

	for (int b = 0; b < nBootstrap; b++) {

		// Resample the data fraction (as if resampling reads, some of which are single-fork vs multi-fork)
		size_t nTotal = data_singleForkDists.size() + data_forkPairDists.size();
		std::binomial_distribution<size_t> fracDist(nTotal, dataFracSingle);
		double bootFracSingle = static_cast<double>(fracDist(bootGen)) / nTotal;

		// Resample single-fork data with replacement and smooth
		std::vector<double> bootSingle(data_singleForkDists.size());
		{
			std::uniform_int_distribution<size_t> dist(0, data_singleForkDists.size() - 1);
			std::normal_distribution<double> noise(0.0, sp1.smoothBW);
			for (size_t i = 0; i < data_singleForkDists.size(); i++) {
				double val = data_singleForkDists[dist(bootGen)] + noise(bootGen);
				bootSingle[i] = std::max(0.0, sp1.mean + (val - sp1.mean) * sp1.shrinkFactor);
			}
			std::sort(bootSingle.begin(), bootSingle.end());
		}

		// Resample fork-pair data with replacement and smooth
		std::vector<double> bootPair;
		if (!data_forkPairDists.empty()) {
			bootPair.resize(data_forkPairDists.size());
			std::uniform_int_distribution<size_t> dist(0, data_forkPairDists.size() - 1);
			std::normal_distribution<double> noise(0.0, sp2.smoothBW);
			for (size_t i = 0; i < data_forkPairDists.size(); i++) {
				double val = data_forkPairDists[dist(bootGen)] + noise(bootGen);
				bootPair[i] = std::max(0.0, sp2.mean + (val - sp2.mean) * sp2.shrinkFactor);
			}
			std::sort(bootPair.begin(), bootPair.end());
		}

		// Compute combined Wasserstein distance at each landscape point for this bootstrap sample
		std::vector<std::pair<double, double>> bootLandscape; // (log10(fr), W)
		std::vector<std::pair<double, double>> iodLandscape;  // (log10(fr), IOD)
		for (const auto& lp : landscapeSimDists) {
			double logFr = std::log10(lp.first);
			double w1 = wassersteinPresorted(lp.second.sortedSingleFork, bootSingle);
			double w2 = (!bootPair.empty() && !lp.second.sortedForkPair.empty()) ? wassersteinPresorted(lp.second.sortedForkPair, bootPair) : 0.0;
			double fracPenalty = (bootFracSingle > 0) ? std::abs(lp.second.fracSingle - bootFracSingle) / bootFracSingle : 0.0;
			bootLandscape.push_back({logFr, w1 / meanSingleFork + w2 / meanForkPair + fracPenalty});
			iodLandscape.push_back({logFr, lp.second.medianIOD});
		}

		// Find the index of the minimum W
		size_t minIdx = 0;
		for (size_t i = 1; i < bootLandscape.size(); i++) {
			if (bootLandscape[i].second < bootLandscape[minIdx].second) minIdx = i;
		}

		// Fit a quadratic to the minimum and its neighbours to interpolate
		double bootIOD;
		if (minIdx > 0 && minIdx < bootLandscape.size() - 1) {
			double x0 = bootLandscape[minIdx - 1].first, w0 = bootLandscape[minIdx - 1].second;
			double x1 = bootLandscape[minIdx].first,     w1 = bootLandscape[minIdx].second;
			double x2 = bootLandscape[minIdx + 1].first, w2 = bootLandscape[minIdx + 1].second;

			// Vertex of parabola through (x0,w0), (x1,w1), (x2,w2)
			double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
			double a = (x2 * (w1 - w0) + x1 * (w0 - w2) + x0 * (w2 - w1)) / denom;

			if (a > 0) {
				double b = (x2 * x2 * (w0 - w1) + x1 * x1 * (w2 - w0) + x0 * x0 * (w1 - w2)) / denom;
				double xStar = -b / (2.0 * a);

				// Clamp to the bracket to avoid extrapolation
				xStar = std::max(x0, std::min(x2, xStar));

				// Linearly interpolate IOD at xStar
				double iod0 = iodLandscape[minIdx - 1].second;
				double iod1_val = iodLandscape[minIdx].second;
				double iod2 = iodLandscape[minIdx + 1].second;
				if (xStar <= x1) {
					double t = (x1 - x0 > 0) ? (xStar - x0) / (x1 - x0) : 0.0;
					bootIOD = iod0 + t * (iod1_val - iod0);
				} else {
					double t = (x2 - x1 > 0) ? (xStar - x1) / (x2 - x1) : 0.0;
					bootIOD = iod1_val + t * (iod2 - iod1_val);
				}
			} else {
				bootIOD = iodLandscape[minIdx].second;
			}
		} else {
			bootIOD = iodLandscape[minIdx].second;
		}

		bootstrapIODs.push_back(bootIOD);
	}

	std::sort(bootstrapIODs.begin(), bootstrapIODs.end());
	double ciLow = bootstrapIODs[static_cast<int>(std::floor(0.025 * nBootstrap))];
	double ciHigh = bootstrapIODs[static_cast<int>(std::ceil(0.975 * nBootstrap)) - 1];

	std::cerr << "Parameter search complete." << std::endl << std::endl;
	std::cerr << "Summary:" << std::endl;
	std::cerr << "--------- " << std::endl;
	std::cerr << "Optimal fr = " << std::scientific << std::setprecision(3) << bestFr << std::endl;
	std::cerr << "Wasserstein distance = " << std::fixed << std::setprecision(3) << bestW << std::endl;
	std::cerr << "Median IOD = " << std::fixed << std::setprecision(0) << bestIOD << " kb" << std::endl;
	std::cerr << "95% confidence interval: [" << std::fixed << std::setprecision(0) << ciLow << " kb , " << std::fixed << std::setprecision(0) << ciHigh << " kb]" << std::endl;

	// Write output file
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d/%m/%Y %H:%M:%S");
	auto startTime = oss.str();
	std::ofstream outFile(args.output);
	outFile << "#DetectFile " << args.DetectInput << "\n";
	outFile << "#ForkFiles " << args.lForkInput << " " << args.rForkInput << "\n";
	outFile << "#OriginFile " << args.originInput << "\n";
	outFile << "#TerminationFile " << args.terminationInput << "\n";
	outFile << "#SystemStartTime " + startTime + "\n";
	outFile << "#Software " << std::string(getExePath()) << "\n";
	outFile << "#Version " << std::string(VERSION) << "\n";
	outFile << "#Commit " << std::string(getGitCommit()) << "\n";
	outFile << "#MeanForkSpeed " << std::fixed << std::setprecision(3) << avgForkSpeed << "\n";
	outFile << "#OptimalFiringRate " << std::scientific << std::setprecision(6) << bestFr << "\n";
	outFile << "#WassersteinDistance " << std::fixed << std::setprecision(6) << bestW << "\n";
	outFile << "#MedianIOD " << std::fixed << std::setprecision(1) << bestIOD << "\n";
	outFile << "#95ConfidenceInterval " << std::fixed << std::setprecision(1) << ciLow << " " << ciHigh << "\n";
	outFile << ">DataSingleForkDistances:\n";
	for (const auto& val : data_singleForkDists) {
		outFile << val << "\n";
	}
	outFile << ">SimSingleForkDistances:\n";
	for (const auto& val : bestSimDist.first) {
		outFile << val << "\n";
	}
	outFile << ">DataForkPairDistances:\n";
	for (const auto& val : data_forkPairDists) {
		outFile << val << "\n";
	}
	outFile << ">SimForkPairDistances:\n";
	for (const auto& val : bestSimDist.second) {
		outFile << val << "\n";
	}

	// Write landscape from points collected during the search (sorted by fr)
	std::sort(landscapePoints.begin(), landscapePoints.end());
	outFile << ">Landscape:\n";
	for (const auto& pt : landscapePoints) {
		outFile << std::scientific << std::setprecision(6) << std::get<0>(pt) << "\t"
				<< std::fixed << std::setprecision(1) << std::get<1>(pt) << "\t"
				<< std::fixed << std::setprecision(4) << std::get<2>(pt) << "\n";
	}
	outFile.close();

	return 0;
}
