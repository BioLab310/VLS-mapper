#ifndef MAPPER_BASIC_H
#define MAPPER_BASIC_H

#endif //MAPPER_BASIC_H

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
#include <stack>
#include <float.h>


#include "output.h"
#include "edlib/include/edlib.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#define isNum(c) (isdigit(c)?c-48:(c=='e'?10:(c=='.'?11:(c=='-'?12:(c=='+'?13:-1)))))

using namespace std;

struct ShortMapResult {
    string ReadName = "";
    string MapFlag = "";
    string RefName = "";
    string MapPosition = "";
    string MapQuality = "";
    string CIGAR = "";
    string MapSame = "";
    string MatePosition = "";
    string SeqMapSize = "";
    string CompleteReadSeq = "";
    string AsciiSeqQuality = "";
    string OptionalArea = "";
};

struct SplitReadFragment {
    int StartPosition = -1;
    int SplitId = -1;
    int SplitReadLength = 0;
    bool IsLastSplit = false;
    vector<ShortMapResult> SplitReadMapResult;
};

struct AnchorPoint {
    string RefName="";
    bool isReverse=false;
    long ReadPosSubRefPos = 0;
    long StartPosInRead = 0;
    long StartPosInRef = 0;
    int segment=0;
};

struct Chain {
    vector<int> chainSequence;
    int length=0;
    long beginPosition=0;
    long endPosition=0;
    long totalDiff=0;
    long shift=0;
    int segment=0;
};

struct Clique {
    vector<int> maximalClique;
    double score=0.0;
};

struct AlignPosition {
    long ReadPos;
    long RefPos;
};

struct BaseLevelPos {
    string RefName="";
    bool dir=true;
    vector<AlignPosition> position;
    double score=0.0;
    long RefStart=0;
    long AnchorNum=0;
    long AnchorNumRev=0;
};

class ChainCmpByRefPos {
public:
    bool operator()(AnchorPoint &a, AnchorPoint &b) {
        return a.StartPosInRef < b.StartPosInRef;
    }
};

class PosCmpByReadPos {
public:
    bool operator()(AlignPosition &a, AlignPosition &b) {
        return a.RefPos < b.RefPos;
    }
};

struct BlPosNum {
    bool operator()(BaseLevelPos& a, BaseLevelPos& b){
        return a.position.size()>b.position.size();
    }
};

void SplitReadSet(double sampleProportion, string readPath);

double StrToNum(string s);
double SamtoolsToGetErrorRate();

double binomial(int N, int k, double p);
double binomialAccumulator(int N, int k, double p);
double CalReadSegLen(double errorRate, double allowRate, double mapRate);

vector<ShortMapResult> SeparateSamInfo();

unordered_map<string, vector<SplitReadFragment>> SplitRead(int SplitLength, string readPath);

string CharTo16bit(string seq);
char * CharTo16Bit(string seq);
string CharToComple(string seq);

string &replace_all(string &src, const string &old_value, const string &new_value);

string getReadName(string ReadFragmentName);

int getReadFragmentOrder(string ReadFragmentName);

unordered_map<string, string> ReadOrRefInfo(string file);

BaseLevelPos chain(vector<SplitReadFragment> OneRead, int SplitLength, int e);
bool checkConnected(AnchorPoint a, AnchorPoint b, int e);
void dfsToFindChain(vector<vector<int>>& AdjList, vector<int>& temp, vector<Chain>& chainResult, vector<int>& used, int index);
vector<Clique> maximalClique(int n, int m, vector<vector<int>> MP);
void dfsToMaximalClique(int d, int an, int sn, int nn, vector<vector<int>>& some, vector<vector<int>>& none, vector<vector<int>>& all, vector<vector<int>> MP, int ans, int n, int m, vector<Clique>& res);
Clique maxClique(vector<Clique> cliques, vector<Chain> chainResult, vector<AnchorPoint> anchors, bool dir);
Clique calculateCliqueScore(Clique& clique, vector<Chain> chainResult, vector<AnchorPoint> anchors, bool dir);
vector<vector<long>> sortClique(Clique& clique, vector<Chain> chainResult, bool dir);
bool checkChainCanAddClique(int chain, vector<vector<long>> cover, vector<Chain> chainResult, Clique& temp, bool dir);
BaseLevelPos getCliqueInfo(Clique maxScoreClique, vector<Chain> chainResult, vector<AnchorPoint> anchors, bool dir, string RefName);
vector<AlignPosition> LIS(vector<AlignPosition>& arr);
vector<AlignPosition> LDS(vector<AlignPosition>& arr);

uint32_t getIandM(string s);
vector<uint32_t> handleCIGAR(string BLresult, uint32_t ReadLen, int DIFF);

vector<BaseLevelPos> chainLDS(vector<SplitReadFragment> OneRead, int SplitLength, int d, long long e);
int calPre(vector<AnchorPoint>& v, int idx, unsigned d);
int calNext(vector<AnchorPoint>& v, int idx, unsigned d);
void getPreLevel(vector<unsigned>& preRes, vector<unsigned>& preLevel, vector<bool>& usedPre, vector<unsigned>& prePath);
void getNextLevel(vector<unsigned>& nextRes, vector<unsigned>& nextLevel, vector<bool>& usedPre, vector<unsigned>& nextPath);

EdlibAlignResult edlibBLCommon(char* read, char* ref);
EdlibAlignResult edlibBLHeadOrTail(char* read, char* ref);
string reverseCIGAR(char* cigar);

