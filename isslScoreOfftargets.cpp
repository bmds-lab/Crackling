/*

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o isslScoreOfftargets isslScoreOfftargets.cpp -O3 -std=c++11 -fopenmp -mpopcnt -Iparallel_hashmap

*/

#include "cfdPenalties.h"

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>
#include <chrono>
#include <bitset>
#include <iostream>
#include <climits>
#include <stdio.h>
#include <cstring>
#include <omp.h>
#include <phmap.h>
#include <map>

using namespace std;

size_t seqLength, seqCount, sliceWidth, sliceCount, offtargetsCount, scoresCount;

vector<uint8_t> nucleotideIndex(256);
vector<char> signatureIndex(4);

size_t getFileSize(const char *path)
{
    struct stat64 statBuf;
    stat64(path, &statBuf);
    return statBuf.st_size;
}

uint64_t sequenceToSignature(const char *ptr)
{
    uint64_t signature = 0;
    for (size_t j = 0; j < seqLength; j++) {
        signature |= (uint64_t)(nucleotideIndex[*ptr]) << (j * 2);
        ptr++;
    }
    return signature;
}

string signatureToSequence(uint64_t signature)
{
    string sequence = string(seqLength, ' ');
    for (size_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s [issltable] [query file] [max distance] [score-threshold] [score-method]\n", argv[0]);
        exit(1);
    }
    
    nucleotideIndex['A'] = 0;
    nucleotideIndex['C'] = 1;
    nucleotideIndex['G'] = 2;
    nucleotideIndex['T'] = 3;
    signatureIndex[0] = 'A';
    signatureIndex[1] = 'C';
    signatureIndex[2] = 'G';
    signatureIndex[3] = 'T';
	
    int maxDist = atoi(argv[3]);
    double threshold = atof(argv[4]);
	string scoreMethod = argv[6];

    FILE *fp = fopen(argv[1], "rb");
    vector<size_t> slicelistHeader(6);
    if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), fp) == 0) {
		fprintf(stderr, "Error reading index: header invalid\n");
		return 1;
	}
	
	offtargetsCount = slicelistHeader[0];
    seqLength       = slicelistHeader[1];
    seqCount        = slicelistHeader[2];
    sliceWidth      = slicelistHeader[3];
    sliceCount      = slicelistHeader[4];
    scoresCount     = slicelistHeader[5];
    
    size_t sliceLimit = 1 << sliceWidth;
    
	// read in the precalculated scores	
	//map<uint64_t, double> precalculatedScores;
	phmap::flat_hash_map<uint64_t, double> precalculatedScores;

	for (int i = 0; i < scoresCount; i++) {
		uint64_t mask = 0;
		double score = 0.0;
		fread(&mask, sizeof(uint64_t), 1, fp);
		fread(&score, sizeof(double), 1, fp);
		
		precalculatedScores.insert(pair<uint64_t, double>(mask, score));
	}
	
	// Load in all of the off-target sites
    vector<uint64_t> offtargets(offtargetsCount);
    if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, fp) == 0) {
		fprintf(stderr, "Error reading index: loading off-target sequences failed\n");
		return 1;
	}
		
	// Create enough 1-bit "seen" flags for the off-targets
	// We only want to score a candidate guide against an off-target once.
	// The least-significant bit represents the first off-target
	// 0 0 0 1   0 1 0 0   would indicate that the 3rd and 5th off-target have been seen.
	// The CHAR_BIT macro tells us how many bits are in a byte (C++ >= 8 bits per byte)
	uint64_t numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;

	
	vector<size_t> allSlicelistSizes(sliceCount * sliceLimit);
    vector<uint64_t> allSignatures(seqCount * sliceCount);
    
    if (fread(allSlicelistSizes.data(), sizeof(size_t), allSlicelistSizes.size(), fp) == 0) {
		fprintf(stderr, "Error reading index: reading slice list sizes failed\n");
		return 1;
	}
	
    if (fread(allSignatures.data(), sizeof(uint64_t), allSignatures.size(), fp) == 0) {
		fprintf(stderr, "Error reading index: reading slice contents failed\n");
		return 1;
	}
    fclose(fp);
    vector<vector<uint64_t *>> sliceLists(sliceCount, vector<uint64_t *>(sliceLimit));
    
    {
        uint64_t *offset = allSignatures.data();
        for (size_t i = 0; i < sliceCount; i++) {
            for (size_t j = 0; j < sliceLimit; j++) {
                size_t idx = i * sliceLimit + j;
                sliceLists[i][j] = offset;
                offset += allSlicelistSizes[idx];
            }
        }
    }
    
    size_t seqLineLength = seqLength + 1;
    size_t fileSize = getFileSize(argv[2]);
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    fp = fopen(argv[2], "rb");
    vector<char> queryDataSet(fileSize);
    vector<uint64_t> querySignatures(queryCount);
    vector<double> querySignatureMitScores(queryCount);
    vector<double> querySignatureCfdScores(queryCount);

    if (fread(queryDataSet.data(), fileSize, 1, fp) < 1) {
        fprintf(stderr, "Failed to read in query file.\n");
        exit(1);
    }
    fclose(fp);

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i < queryCount; i++) {
            char *ptr = &queryDataSet[i * seqLineLength];
            uint64_t signature = sequenceToSignature(ptr);
            querySignatures[i] = signature;
        }
    }

    #pragma omp parallel
    {
        unordered_map<uint64_t, unordered_set<uint64_t>> searchResults;
		vector<uint64_t> offtargetToggles(numOfftargetToggles);
    
		uint64_t * offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;

        #pragma omp for
        for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {

			auto searchSignature = querySignatures[searchIdx];

			double totScoreMit = 0.0;
			double totScoreCfd = 0.0;
            int numOffTargetSitesScored = 0;

            double maximum_sum = (10000.0 - threshold*100) / threshold;
			bool checkNextSlice = true;
			
            for (size_t i = 0; i < sliceCount; i++) {
                uint64_t sliceMask = sliceLimit - 1;
                int sliceShift = sliceWidth * i;
                sliceMask = sliceMask << sliceShift;
                auto &sliceList = sliceLists[i];
                
                uint64_t searchSlice = (searchSignature & sliceMask) >> sliceShift;
                
                size_t idx = i * sliceLimit + searchSlice;
                
                size_t signaturesInSlice = allSlicelistSizes[idx];
                uint64_t *sliceOffset = sliceList[searchSlice];
				
                for (size_t j = 0; j < signaturesInSlice; j++) {
                    auto signatureWithOccurrencesAndId = sliceOffset[j];
                    auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFull;

					uint64_t xoredSignatures = searchSignature ^ offtargets[signatureId];
					uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAull;
					uint64_t oddBits = xoredSignatures & 0x5555555555555555ull;
					uint64_t mismatches = (evenBits >> 1) | oddBits;
					int dist = __builtin_popcountll(mismatches);

					// Begin calculating MIT score
					if (!scoreMethod.compare("mit")) {
						if (dist > 0 && dist <= maxDist) {
							uint64_t seenOfftargetAlready = 0;
							uint64_t * ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
							if (i > 0) {
								seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;
							}
							
							if (!seenOfftargetAlready) {
								uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));
								totScoreMit += precalculatedScores[mismatches] * (double)occurrences;
								
								if (totScoreMit > maximum_sum) {
									checkNextSlice = false;
									break;
								}
								
								*ptrOfftargetFlag |= (1ULL << (signatureId % 64));
								numOffTargetSitesScored += occurrences;
							}
						}
						
					} 
					
					// Begin calculating CFD score
					else if (!scoreMethod.compare("cfd")) {
						if (dist == 1) {
							totScoreCfd = 1.0;
						}
						else if (dist <= maxDist) {
							uint64_t seenOfftargetAlready = 0;
							uint64_t * ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
							if (i > 0) {
								seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;
							}
							
							if (!seenOfftargetAlready) {
								uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));
								double cfdScore = cfdPamPenalties[0b1010]; // PAM: NGG
								
								for (size_t pos = 0; pos < 20; pos++) {
									size_t mask = pos << 4;
									
									// Create the mask to look up the position-identity score
									// In Python... c2b is char to bit
									// 	mask = pos << 4
									// 	mask |= c2b[sgRNA[pos]] << 2
									// 	mask |= c2b[revcom(offTaret[pos])]
									
									// Find identity at `pos` for search signature
									// example: find identity in pos=2
									// 	Recall ISSL is inverted, hence:
									//              3'-  T  G  C  C  G  A -5'
									//	start		    11 10 01 01 10 00 	
									//	3UL << pos*2    00 00 00 11 00 00 
									//  and			    00 00 00 01 00 00
									//  shift		    00 00 00 00 01 00
									uint64_t searchSigIdentityPos = searchSignature;
									searchSigIdentityPos &= (3UL << (pos * 2));
									searchSigIdentityPos = searchSigIdentityPos >> (pos * 2); 
									searchSigIdentityPos = searchSigIdentityPos << 2;

									// Find identity at `pos` for offtarget
									// example: find identity in pos=2
									// 	Recall ISSL is inverted, hence:
									//              3'-  T  G  C  C  G  A -5'
									//	start		    11 10 01 01 10 00 	
									//	3UL<<pos*2      00 00 00 11 00 00 
									//  and			    00 00 00 01 00 00
									//  shift		    00 00 00 00 00 01
									//  rev comp 3UL    00 00 00 00 00 10 (done below)
									uint64_t offtargetIdentityPos = offtargets[signatureId];
									offtargetIdentityPos &= (3UL << (pos * 2));
									offtargetIdentityPos = offtargetIdentityPos >> (pos * 2); 

									// Complete the mask
									// reverse complement (^3UL) `offtargetIdentityPos` here
									mask = (mask | searchSigIdentityPos | (offtargetIdentityPos ^ 3UL));

									if (searchSigIdentityPos >> 2 != offtargetIdentityPos) {
										cfdScore *= cfdPosPenalties[mask];
									}
									
								}
								totScoreCfd += cfdScore;
								
								if (totScoreCfd > maximum_sum) {
									checkNextSlice = false;
									break;
								}
							
								*ptrOfftargetFlag |= (1ULL << (signatureId % 64));
								numOffTargetSitesScored += occurrences;
							}
						}
					}
                }
				
				if (!checkNextSlice)
					break;
            }
			
			querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
			querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);

			memset(offtargetToggles.data(), 0, sizeof(uint64_t)*offtargetToggles.size());
        }
		
    }
	
	for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
		auto querySequence = signatureToSequence(querySignatures[searchIdx]);
		printf("%s\t", querySequence.c_str());
		if (!scoreMethod.compare("cfd"))
			printf("%f\n", querySignatureCfdScores[searchIdx]);
		if (!scoreMethod.compare("mit"))
			printf("%f\n", querySignatureMitScores[searchIdx]);
	}

    return 0;
}