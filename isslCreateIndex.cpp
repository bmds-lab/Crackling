/*

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o isslCreateIndex isslCreateIndex.cpp -O3 -std=c++11 -fopenmp -mpopcnt

*/


#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <map>

using namespace std;

size_t seqLength;
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

// combinations
vector<uint64_t> computeMasksTwoBit(int seqLength, int mismatches) {
	vector<uint64_t> masks;
	
	// there are more positions than mismatches
	if (mismatches < seqLength) {
		
		if (mismatches > 0) {
			// some mismatches across a long sequence
			for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches -1)) {
				// orig: masks.push_back((1 << (seqLength - 1)) + mask);
				masks.push_back((1LLU << (seqLength - 1)*2) + mask);
			}
			for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches)) {
				masks.push_back(mask);
			}
		} else {
			// no mismatches at all
			masks.push_back(0LLU);
		}
		
	// mismatches >= seqLength
	// every position is going to be a mismatch
	// eg: mismatches=4 and seqLength=4, we want: 10 10 10 10
	} else {
		// orig: masks.push_back((1LLU << seqLength) - 1);
		uint64_t tempMask = 0;
		for (int i = 0; i < seqLength; i++) {
			tempMask |= (1LLU << i*2);
		}
		masks.push_back(tempMask);
	}
	return masks;
}

double single_score(int* mismatch_array, int length) {
    int i;
    double T1=1.0, T2, T3, d=0.0, score;
    double M[] = {0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583};

    /* 1st term */
    for(i=0; i<length; ++i)
            T1 = T1*(1.0-M[mismatch_array[i]]);

    /* 2nd term */
    if(length==1)
            d = 19.0;
    else {
            for(i=0; i<length-1; ++i)
                    d += mismatch_array[i+1]-mismatch_array[i];
            d = d/(length-1);
    }
    T2 = 1.0 / ((19.0-d)/19.0 * 4.0 + 1);

    /* 3rd term */
    T3 = 1.0 / (length*length);

    /* Total score */
    score = T1*T2*T3*100;
    return score;
}

double sscore(uint64_t xoredSignatures)
{
    int mismatch_array[20], m = 0;
    for (size_t j = 0; j < seqLength; j++) {
        if ((xoredSignatures >> (j * 2)) & 0x3) {
            mismatch_array[m++] = j;
        }
    }
    if (m == 0) return 0.0;
    return single_score(mismatch_array, m);
}

int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr, "Usage: %s [offtargetSites.txt] [sequence length] [slice width (bits)] [sissltable]\n", argv[0]);
        exit(1);
    }
    size_t fileSize = getFileSize(argv[1]);
    
    FILE *fp = fopen(argv[1], "rb");
    seqLength = atoi(argv[2]);
    if (seqLength > 32) {
        fprintf(stderr, "Sequence length is greater than 32, which is the maximum supported currently\n");
        exit(1);
    }
    size_t seqLineLength = seqLength + 1; // '\n'
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "fileSize: %zu\n", fileSize);
        fprintf(stderr, "Error: file does is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t seqCount = fileSize / seqLineLength;
    fprintf(stderr, "Number of sequences: %zu\n", seqCount);
    
    
    nucleotideIndex['A'] = 0;
    nucleotideIndex['C'] = 1;
    nucleotideIndex['G'] = 2;
    nucleotideIndex['T'] = 3;
    signatureIndex[0] = 'A';
    signatureIndex[1] = 'C';
    signatureIndex[2] = 'G';
    signatureIndex[3] = 'T';
    
    size_t globalCount = 0;
    
    vector<uint64_t> seqSignatures;
    vector<uint32_t> seqSignaturesOccurrences;
    
	size_t myCount = 0;
    {
        vector<char> entireDataSet(fileSize);

        if (fread(entireDataSet.data(), fileSize, 1, fp) < 1) {
            fprintf(stderr, "Failed to read in file.\n");
            exit(1);
        }
        fclose(fp);

		size_t i = 0;
		size_t offtargetId = 0;
		while (i < seqCount) {
			char *ptr = &entireDataSet[i * seqLineLength];
			
			uint64_t signature = sequenceToSignature(ptr);

			// check how many times the off-target appears
			// (assumed the list is sorted)
			uint32_t occurrences = 1;
			while (memcmp(ptr, ptr + (seqLineLength * occurrences), seqLength) == 0) {
				occurrences++;
			}

			seqSignatures.push_back(signature);
			seqSignaturesOccurrences.push_back(occurrences);
			
			myCount++;
			if (myCount % 10000 == 0) {
				fprintf(stderr, "%zu/%zu : %zu\n", myCount, seqCount, i);
			}
			
			i += occurrences;
		}
    
    }
	printf("Finished counting occurrences, now constructing index...\n");
    size_t sliceWidth = atoi(argv[3]);
    size_t sliceLimit = 1 << sliceWidth;
    size_t sliceCount = (seqLength * 2) / sliceWidth;
    size_t offtargetsCount = myCount;
	
    vector<vector<vector<uint64_t>>> sliceLists(sliceCount, vector<vector<uint64_t>>(sliceLimit));
    
    #pragma omp parallel for
    for (size_t i = 0; i < sliceCount; i++) {
        uint64_t sliceMask = sliceLimit - 1;
        int sliceShift = sliceWidth * i;
        sliceMask = sliceMask << sliceShift;
        auto &sliceList = sliceLists[i];
        
		uint32_t signatureId = 0;
        for (uint64_t signature : seqSignatures) {
			uint32_t occurrences = seqSignaturesOccurrences[signatureId];
            uint8_t sliceVal = (signature & sliceMask) >> sliceShift;
			
			uint64_t seqSigIdVal = (((uint64_t)occurrences) << 32) | (uint64_t)signatureId;
			sliceList[sliceVal].push_back(seqSigIdVal);
			signatureId++;
        }
    }
    
	printf("Finished constructing index, now precalculating scores...\n");
	
	// Precalculate all the scores
	map<uint64_t, double> precalculatedScores;
	
	int maxDist = seqLength * 2 / sliceWidth - 1;
	size_t scoresCount = 0;
	
	for (int i = 1; i <= maxDist; i++) {
		vector<uint64_t> tempMasks;
		tempMasks = computeMasksTwoBit(20, i);
		for (auto mask : tempMasks) {
			double score = sscore(mask);
			precalculatedScores.insert(pair<uint64_t, double>(mask, score));
			scoresCount++;
		}
	}
	
	printf("Finished calculating scores, now preparing to write to disk...\n");
	
    fp = fopen(argv[4], "wb");
    vector<size_t> slicelistHeader;
    slicelistHeader.push_back(offtargetsCount);
    slicelistHeader.push_back(seqLength);
    slicelistHeader.push_back(seqCount);
    slicelistHeader.push_back(sliceWidth);
    slicelistHeader.push_back(sliceCount);
    slicelistHeader.push_back(scoresCount);
	
	
	// write the header
    fwrite(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), fp);

	// write the precalculated scores
	for (auto const& x : precalculatedScores) {
		fwrite(&x.first, sizeof(uint64_t), 1, fp);
		fwrite(&x.second, sizeof(double), 1, fp);
	}

	// write the offtargets
	fwrite(seqSignatures.data(), sizeof(uint64_t), seqSignatures.size(), fp);
	
    for (size_t i = 0; i < sliceCount; i++) { // 5
        for (size_t j = 0; j < sliceLimit; j++) { // 256
            size_t sz = sliceLists[i][j].size(); 
            fwrite(&sz, sizeof(size_t), 1, fp);
        }
    }

    for (size_t i = 0; i < sliceCount; i++) { // 5
        for (size_t j = 0; j < sliceLimit; j++) { // 256
            fwrite(sliceLists[i][j].data(), sizeof(uint64_t), sliceLists[i][j].size(), fp); // vector
        }
    }
	
	printf("Writing to disk...\n");
    
    fclose(fp);
    printf("Done.\n");
    return 0;
}
