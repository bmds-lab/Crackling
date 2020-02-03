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

#include <chrono>    // for high_resolution_clock

using namespace std;

size_t seqLength, seqCount, sliceWidth, sliceCount;

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

/**************************************/
/**                    single_score()                    **/
/**************************************/

/*
 Computes the score of a single off-target site
 Input: array containing all the mismatches between this site and the target, length of the array
 Output: score
 */

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

int distance(uint64_t xoredSignatures)
{
    uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAull;
    uint64_t oddBits = xoredSignatures & 0x5555555555555555ull;
    return __builtin_popcountll((evenBits >> 1) | oddBits);
}

int main(int argc, char **argv)
{
    
    long startTimer, endTimer;
    struct timeval timecheck;
    
    gettimeofday(&timecheck, NULL);
    startTimer = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    
    
    
    if (argc < 3) {
        fprintf(stderr, "Usage: %s [sissltable] [query file] [max distance]\n", argv[0]);
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

    //auto start = std::chrono::high_resolution_clock::now();

    FILE *fp = fopen(argv[1], "rb");
    vector<size_t> slicelistHeader(4);
    fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), fp);
    seqLength = slicelistHeader[0];
    seqCount = slicelistHeader[1];
    sliceWidth = slicelistHeader[2];
    sliceCount = slicelistHeader[3];
    
    size_t sliceLimit = 1 << sliceWidth;
    vector<size_t> allSlicelistSizes(sliceCount * sliceLimit);
    vector<uint64_t> allSignatures(seqCount * sliceCount);
    
    fread(allSlicelistSizes.data(), sizeof(size_t), allSlicelistSizes.size(), fp);
    fread(allSignatures.data(), sizeof(uint64_t), allSignatures.size(), fp);
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
    
		//auto startB = std::chrono::high_resolution_clock::now();
        #pragma omp for
        for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            auto searchSignature = querySignatures[searchIdx];
            auto &searchSigs = searchResults[searchSignature];
			
			double totscore = 0.0;
            int numOffTargetSitesScored = 0;
            double threshold = 75.0;
            double maximum_sum = (10000.0 - threshold*100) / threshold;
			bool checkNextSlice = true;
			
            for (size_t i = 0; i < sliceCount; i++) { //  && checkNextSlice
                uint64_t sliceMask = sliceLimit - 1;
                int sliceShift = sliceWidth * i;
                sliceMask = sliceMask << sliceShift;
                auto &sliceList = sliceLists[i];
                
                uint64_t searchSlice = (searchSignature & sliceMask) >> sliceShift;
                
                size_t idx = i * sliceLimit + searchSlice;
                
                size_t signaturesInSlice = allSlicelistSizes[idx];
                uint64_t *sliceOffset = sliceList[searchSlice];
                
                for (size_t j = 0; j < signaturesInSlice; j++) {
                    auto signature = sliceOffset[j];
                    int dist = distance(searchSignature ^ signature);
                    if (dist <= maxDist && searchSigs.insert(signature).second) {
                        //searchSigs.insert(signature);
						
						totscore += sscore(searchSignature ^ signature);
						
						if (totscore > maximum_sum) {
							checkNextSlice = false;
							break;
						}
                    }
                }
				
				if (!checkNextSlice)
					break;
            }
			
			auto querySequence = signatureToSequence(searchSignature);
			#pragma omp critical
			{
				printf("%s\t", querySequence.c_str());
				printf("%f\n", 10000.0 / (100.0 + totscore));
			}
        }
		
    }
    
    
    gettimeofday(&timecheck, NULL);
    endTimer = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
    
    printf("%ld milliseconds elapsed\n", (endTimer - startTimer));
    
    
    
    return 0;
}
