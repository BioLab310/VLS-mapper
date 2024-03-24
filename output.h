//
// Created by spoony on 23-3-11.
//

#ifndef MAPPER_OUTPUT_H
#define MAPPER_OUTPUT_H

#include <string>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <deque>
#include <queue>
#include <string.h>

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "assert.h"

using namespace std;

typedef struct {
    samFile *output_sam_file;
    sam_hdr_t *sam_header;
    vector <bam1_t> samOutputVec;
} OutputQueue;


void initialize_output_queue(const char *output_file_path, unordered_map<string, string>& RefInfo, OutputQueue* output_queue, int argc,char *argv[]);
void output_sam_header(const char *output_file_path, unordered_map<string, string> RefInfo, OutputQueue* output_queue, int argc,char *argv[]);
void output_Queue(OutputQueue * a);
void generate_BAM_t(string readName,uint16_t flag, int refId, long refstart,uint8_t MapQ,\
vector<uint32_t>& cigar, char a, int b, int c, string readSeq, string readQual, string options,bam1_t* r);
#endif //MAPPER_OUTPUT_H
