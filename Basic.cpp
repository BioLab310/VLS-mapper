#include "Basic.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Waddress-of-temporary"
#pragma clang diagnostic ignored "-Wc++11-narrowing"

void SplitReadSet(double sampleProportion, string readPath) {
    double SplitReadStartTime = clock();
    ifstream finPre(readPath);
    string EveryLine;
    unsigned int countAll = 0;
    while (getline(finPre, EveryLine)) {
        if (EveryLine[0] == '>') {
            ++countAll;
        }
    }
    finPre.close();

    ifstream fin(readPath);
    ofstream fout("readSample.fasta");
    unsigned int position = 0;
    unsigned int countSample = 0;
    unsigned int sampleNumber = (int) (countAll * sampleProportion);
    while (getline(fin, EveryLine)) {
        if (EveryLine[0] == '>') {
            ++position;
            if (position % ((int) (1 / sampleProportion)) == 0) {
                ++countSample;
                fout << EveryLine << endl;
                while (getline(fin, EveryLine)) {
                    if (EveryLine[0] != '>') {
                        fout << EveryLine << endl;
                    } else {
                        ++position;
                        break;
                    }
                }
                if (countSample >= sampleNumber) {
                    break;
                }
            }
        }
    }
    fin.close();
    fout.close();
    double SplitReadEndTime = clock();
    cout << "The run time of splitting readset is: "
         << (double) (SplitReadEndTime - SplitReadStartTime) / CLOCKS_PER_SEC
         << "s" << std::endl;
}

double StrToNum(string s) {
    //字符串转数字，包括整数、小数和科学记数法
    int i, j, k, negative = 0;
    double n = 0;
    string s1, s2;

    if (s.empty()) return 0;
    if (s[0] == '-') negative = 1; //设置负数标记
    if (s[0] == '+' || s[0] == '-') s = s.substr(1, s.size());
    //---------------
    for (i = 0; i < s.size(); i++) //排除不需要的字符
        if (isNum(s[i]) == -1) return pow(-1.1, 1.1);
    if (s[0] == 'e' || s[0] == '.' || s[s.size() - 1] == 'e' || s[s.size() - 1] == '.')
        return pow(-1.1, 1.1); //排除 e或. 出现在首尾
    i = -1;
    j = 0;
    while ((i = s.find('.', ++i)) != s.npos) j++;
    if (j > 1) return pow(-1.1, 1.1); //排除多个小数点
    i = -1;
    j = 0;
    while ((i = s.find('e', ++i)) != s.npos) j++;
    if (j > 1) return pow(-1.1, 1.1); //排除多个字母e
    if (s.find('e') == s.npos) //没有e时排除加减
        if (s.find('+') != s.npos || s.find('-') != s.npos) return pow(-1.1, 1.1);
    //---------------
    if ((i = s.find('e')) != s.npos) {
        s1 = s.substr(0, i); //尾数部分
        s2 = s.substr(i + 1, s.size()); //阶码
        if (s2[0] == '+') s2 = s2.substr(1, s2.size()); //阶码为正数，去掉+
        if (s2.find('.') != s2.npos) return pow(-1.1, 1.1); //阶码不准出现小数
        n = StrToNum(s1) * pow(10.0, StrToNum(s2)); //尾数和阶码分别递归调用
        return negative ? -n : n;
    }
    i = 0;
    k = 1;
    if ((i = s.find('.')) != s.npos) {
        for (j = i + 1; j < s.length(); j++, k++)
            n += isNum(s[j]) / pow(10.0, (double) k);
        n += StrToNum(s.substr(0, i));  //整数部分递归调用
    } else
        for (j = 0; j < s.size(); j++)
            n = n * 10 + isNum(s[j]);

    return negative ? -n : n; //负数返回-n
}

double SamtoolsToGetErrorRate() {
    const char *SamToBam = "./samtools view -bS readSample.sam > readSample.bam";
    const char *stats = "./samtools stats --threads 10 readSample.bam > readSample.txt";
    const char *OutStatsResult = "grep \"^SN\" readSample.txt";
    system(SamToBam);
    system(stats);
    freopen("samInfo.txt", "w", stdout);
    system(OutStatsResult);
    freopen("/dev/tty", "w", stdout);
    ifstream fin("samInfo.txt");
    string EveryLine;
    string errorPre = "error rate:";
    string originErrorRate = "";
    while (getline(fin, EveryLine)) {
        if (EveryLine.find(errorPre) != string::npos) {
            int i = EveryLine.find("error rate:") + errorPre.size();
            while (isspace(EveryLine[i])) {
                ++i;
            }
            while (!isspace(EveryLine[i])) {
                originErrorRate += EveryLine[i];
                ++i;
            }
            cout << "Sample error rate: " << originErrorRate << endl;
            break;
        }
    }
    fin.close();
    double errorRate = StrToNum(originErrorRate);
    return errorRate;
}

double binomial(int N, int k, double p) {
    if (N == 0 && k == 0) {
        return 1;
    }
    if (N < 0 || k < 0) {
        return 0;
    }
    if (k > N) {
        return 0;
    }

    std::deque<std::vector<double>> deqVecDou;
    deqVecDou.push_back(std::vector<double>{1, 0});

    double q = 1;
    std::vector<double> lieOne(N + 1);
    for (size_t i = 0; i != N + 1; ++i) {
        lieOne[i] = q;
        q *= 1 - p;
    }

    double pp = 1 - p;
    for (size_t i = 1; i != N + 1; ++i) {
        if (i < k) {
            deqVecDou.push_back(std::vector<double>(i + 2, 0));
            deqVecDou[1][0] = lieOne[i];
            for (size_t j = 1; j != i + 1; ++j) {
                deqVecDou[1][j] = pp * deqVecDou[0][j] + p * deqVecDou[0][j - 1];
            }
            deqVecDou.pop_front();
        } else {
            deqVecDou.push_back(std::vector<double>(k + 2, 0));
            deqVecDou[1][0] = lieOne[i];
            for (size_t j = 1; j != k + 1; ++j) {
                deqVecDou[1][j] = pp * deqVecDou[0][j] + p * deqVecDou[0][j - 1];
            }
            deqVecDou.pop_front();
        }
    }
    double result = deqVecDou[0][k];
    deqVecDou.clear();
    lieOne.clear();
    return result;
}

double binomialAccumulator(int N, int k, double p) {
    double res = 0.0;
    for (int i = 0; i <= k; ++i) {
        res += binomial(N, i, p);
    }
    return res;
}

double CalReadSegLen(double errorRate, double allowRate, double mapRate) {
    int ReadSegLen = 0;
    int left = 0, right = 1000;
    while (left <= right) {
        int mid = left + ((right - left) >> 1);
        double rate = binomialAccumulator(mid, (int) (allowRate * mid), errorRate);
        if (rate <= mapRate) {
            right = mid - 1;
        } else {
            left = mid + 1;
            ReadSegLen = mid;
        }
    }
    double SplitReadEndTime = clock();
    cout << "Read split length: " << ReadSegLen << endl;
    return ReadSegLen;
}

unordered_map<string, vector<SplitReadFragment>> SplitRead(int SplitLength, string readPath) {
    double SplitReadStartTime = clock();
    ifstream fin(readPath);
    ofstream fout("readSegment.fasta");
    string EveryLine;
    string ReadName;
    unordered_map<string, vector<SplitReadFragment>> SplitReadResult;

    getline(fin, EveryLine);
    while (!EveryLine.empty()) {
        if (EveryLine[0] == '>') {
            ReadName = "";
            for (int i = 1; i < EveryLine.size() && EveryLine[i] != ' '; ++i) {
                ReadName += EveryLine[i];
            }
            getline(fin, EveryLine);
        }
        else {
            string read = "";
            while (!EveryLine.empty() && EveryLine[0] != '>') {
                read += EveryLine;
                getline(fin, EveryLine);
            }
            int CountFragmentNumber = 0;
            bool still = true;
            for (int i = 0; i < read.size();) {
                SplitReadFragment OneFragment;
                string SplitName = ReadName + '.' + to_string(CountFragmentNumber);
                OneFragment.SplitId = CountFragmentNumber;
                OneFragment.StartPosition = 0 + (int)(SplitLength * 0.5) * CountFragmentNumber;
                string FragmentSeq;
                for (int j = 0; j < SplitLength && i < read.size(); ++j, ++i) {
                    FragmentSeq += read[i];
                }

                if (i == read.size()) {
                    OneFragment.IsLastSplit = true;
                    still=false;
                    OneFragment.SplitReadLength = EveryLine.size() - (int)(CountFragmentNumber * SplitLength * 0.5);
                } else {
                    OneFragment.SplitReadLength = SplitLength;
                }
                SplitReadResult[ReadName].emplace_back(OneFragment);
                ++CountFragmentNumber;

                fout << '>' << SplitName << endl;
                fout << FragmentSeq << endl;
                FragmentSeq = "";
                if(still) {
                    i-=(int)(SplitLength * (1-0.5));
                }
            }


        }
    }
    fin.close();
    fout.close();
    double SplitReadEndTime = clock();
    cout << "The run time of splitting read is: "
         << (double) (SplitReadEndTime - SplitReadStartTime) / CLOCKS_PER_SEC
         << "s" << std::endl;
    return SplitReadResult;
}

vector<ShortMapResult> SeparateSamInfo() {
    double SeparateSamStartTime = clock();
    ifstream f("readSegment.sam");
    string EveryLine;
    vector<ShortMapResult> SeparateResult;

    while (getline(f, EveryLine)) {
        if (EveryLine[0] == '@') {
            continue;
        }

        ShortMapResult OneResult;

        int i = 0;

        //Query name of th read or the read pair
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.ReadName += EveryLine[i];
            ++i;
        }
        ++i;

        //Bitwise flag
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.MapFlag += EveryLine[i];
            i++;
        }
        ++i;

        //Reference sequence name
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.RefName += EveryLine[i];
            i++;
        }
        ++i;

        //1-Based leftmost position of clipped alignment
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.MapPosition += EveryLine[i];
            i++;
        }
        ++i;

        //Mapping quality
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.MapQuality += EveryLine[i];
            i++;
        }
        ++i;

        //Compact Idiosyncratic Gapped Alignment Report
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.CIGAR += EveryLine[i];
            i++;
        }
        ++i;

        //'=' if same as ref
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.MapSame += EveryLine[i];
            i++;
        }
        ++i;

        //1-Based leftmost mate position
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.MatePosition += EveryLine[i];
            i++;
        }
        ++i;

        //Inferred insert size
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.SeqMapSize += EveryLine[i];
            i++;
        }
        ++i;

        //Complete read sequence
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.CompleteReadSeq += EveryLine[i];
            i++;
        }
        ++i;

        //ASCII format sequence quality
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.AsciiSeqQuality += EveryLine[i];
            i++;
        }
        ++i;

        //Optional area
        while (i < EveryLine.size() && !isspace(EveryLine[i])) {
            OneResult.OptionalArea += EveryLine[i];
            i++;
        }
        ++i;

        if (OneResult.MapFlag != "4") {
            SeparateResult.push_back(OneResult);
        }

        EveryLine = "";

        //break;
    }
    f.close();
    double SeparateSamEndTime = clock();
    cout << "The run time of splitting sam is: "
         << (double) (SeparateSamEndTime - SeparateSamStartTime) / CLOCKS_PER_SEC << "s" << std::endl;
    return SeparateResult;
}

string &replace_all(string &src, const string &old_value, const string &new_value) {
    // 每次重新定位起始位置，防止上轮替换后的字符串形成新的old_value
    for (string::size_type pos(0); pos != string::npos; pos += new_value.length()) {
        if ((pos = src.find(old_value, pos)) != string::npos) {
            src.replace(pos, old_value.length(), new_value);
        } else break;
    }
    return src;
}

string CharTo16bit(string seq) {
    string res = seq;
    replace_all(res, "a", "\\x01");
    replace_all(res, "c", "\\x02");
    replace_all(res, "g", "\\x04");
    replace_all(res, "t", "\\x08");
    replace_all(res, "A", "\\x01");
    replace_all(res, "C", "\\x02");
    replace_all(res, "G", "\\x04");
    replace_all(res, "T", "\\x08");
    return res;
}

char * CharTo16Bit(string seq) {
    char *res=(char*) malloc((seq.size()+1) * sizeof(char));
    for(uint32_t i=0; i<seq.size(); ++i) {
        switch(seq[i]) {
            case 'A':
            case 'a':
                res[i]='\x01';
                break;
            case 'C':
            case 'c':
                res[i]='\x02';
                break;
            case 'G':
            case 'g':
                res[i]='\x04';
                break;
            case 'T':
            case 't':
                res[i]='\x08';
                break;
        }
    }
    res[seq.size()]='\0';
    return res;
}

string CharToComple(string seq) {
    string res = seq;
    replace_all(res, "a", "T");
    replace_all(res, "c", "G");
    replace_all(res, "g", "C");
    replace_all(res, "t", "A");
    replace_all(res, "A", "T");
    replace_all(res, "C", "G");
    replace_all(res, "G", "C");
    replace_all(res, "T", "A");
    return res;
}

string getReadName(string ReadFragmentName) {
    string ReadName = "";
    for (int i = ReadFragmentName.size() - 1; i >= 0; --i) {
        if (ReadFragmentName[i] != '.') {
            continue;
        }
        ReadName = ReadFragmentName.substr(0, i);
        break;
    }
    return ReadName;
}

int getReadFragmentOrder(string ReadFragmentName) {
    int order = 0;
    for (int i = ReadFragmentName.size() - 1; i >= 0; --i) {
        if (ReadFragmentName[i] != '.') {
            continue;
        }
        order = stoi(ReadFragmentName.substr(i + 1));
        break;
    }
    return order;
}

unordered_map<string, string> ReadOrRefInfo(string file) {
    ifstream fin(file);

    string EveryLine = "";
    string ReadName = "";
    string ReadSeq = "";
    unordered_map<string, string> res;
    getline(fin, EveryLine);
    while (!EveryLine.empty()) {
        if (EveryLine[0] == '>') {
            ReadName = "";
            for (int i = 1; i < EveryLine.size() && EveryLine[i] != ' '; ++i) {
                ReadName += EveryLine[i];
            }
            getline(fin, EveryLine);
        } else {
            while (!EveryLine.empty() && EveryLine[0] != '>') {
                ReadSeq += EveryLine;
                getline(fin, EveryLine);
            }
            res[ReadName] = ReadSeq;
            ReadName = "";
            ReadSeq = "";
        }
    }
    fin.close();

    return res;
}

vector<BaseLevelPos> chainLDS(vector<SplitReadFragment> OneRead, int SplitLength, int d, long long e) {
    unordered_map<string, vector<AnchorPoint>> UseRefNameToDivide;
    for (int i = 0; i < OneRead.size(); ++i) {
        for (int j = 0; j < OneRead[i].SplitReadMapResult.size(); ++j) {
            AnchorPoint anchor;
            if (stoi(OneRead[i].SplitReadMapResult[j].MapFlag) & 16) {
                anchor.isReverse = true;
            }
            anchor.segment = OneRead[i].SplitReadLength;
            anchor.StartPosInRead = (int)(i * SplitLength * 0.5);
            anchor.StartPosInRef = stoi(OneRead[i].SplitReadMapResult[j].MapPosition);
            anchor.ReadPosSubRefPos = anchor.StartPosInRead - anchor.StartPosInRef;
            anchor.RefName = OneRead[i].SplitReadMapResult[j].RefName;
            UseRefNameToDivide[anchor.RefName].emplace_back(anchor);
        }
    }

    priority_queue<BaseLevelPos, vector<BaseLevelPos>, BlPosNum> Best5;

    vector<BaseLevelPos> OneReadResAll;
    for (auto OneRef: UseRefNameToDivide) {
        cout<<"-------OneRef-------"<<endl;
        vector<AnchorPoint> forwardAnchor;
        vector<AnchorPoint> reverseAnchor;
        for (auto Anchor: OneRef.second) {
            if (Anchor.isReverse) {
                reverseAnchor.emplace_back(Anchor);
            } else {
                forwardAnchor.emplace_back(Anchor);
            }
        }
        int n = forwardAnchor.size();
        int m = reverseAnchor.size();
        cout<<n<<"++++++++++++++++++++++++++++++++"<<m<<endl;
        sort(forwardAnchor.begin(), forwardAnchor.end(), ChainCmpByRefPos());
        sort(reverseAnchor.begin(), reverseAnchor.end(), ChainCmpByRefPos());

//---------------------------------Adj---------------------------------------
        vector<vector<int>> AdjList(n, vector<int>());
        vector<vector<int>> AdjListReverse(m, vector<int>());

        vector<int> indegree(n, 0);
        vector<int> outdegree(n, 0);
        vector<int> alldegree(n, 0);

        vector<int> indegreeReverse(m, 0);
        vector<int> outdegreeReverse(m, 0);
        vector<int> alldegreeReverse(m, 0);

        // AdjList
        for (int i = 0; i < n; ++i) {
            AdjList[i].emplace_back(i);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (checkConnected(forwardAnchor[i], forwardAnchor[j], e)) {
                    ++indegree[j];
                    ++outdegree[i];
                    AdjList[i].emplace_back(j);
                }
            }
        }
        int maxDegree=0;
        vector<int> maxDegreeIndex;
        for (int i = 0; i < n; ++i) {
            alldegree[i]=indegree[i]+outdegree[i];
            if(alldegree[i]>maxDegree) {
                maxDegree=alldegree[i];
                maxDegreeIndex.clear();
                maxDegreeIndex.emplace_back(i);
            } else if(alldegree[i]==maxDegree) {
                maxDegreeIndex.emplace_back(i);
            }
        }

        // AdjListReverse
        for (int i = m - 1; i >= 0; --i) {
            AdjListReverse[i].emplace_back(i);
        }
        for (int i = m - 1; i >= 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                if (checkConnected(reverseAnchor[i], reverseAnchor[j], e)) {
                    ++indegreeReverse[j];
                    ++outdegreeReverse[i];
                    AdjListReverse[i].emplace_back(j);
                }
            }
        }
        int maxDegreeReverse=0;
        vector<int> maxDegreeIndexReverse;
        for (int i = 0; i < m; ++i) {
            alldegreeReverse[i]=indegreeReverse[i]+outdegreeReverse[i];
            if(alldegreeReverse[i]>maxDegreeReverse) {
                maxDegreeReverse=alldegreeReverse[i];
                maxDegreeIndexReverse.clear();
                maxDegreeIndexReverse.emplace_back(i);
            } else if(alldegreeReverse[i]==maxDegreeReverse) {
                maxDegreeIndexReverse.emplace_back(i);
            }
        }
//---------------------------------Adj---------------------------------------


//        ----------forward----------
        vector<unsigned> preRes(n, n);
        vector<bool> usedPre(n, false);
        for(unsigned  i=0; i<n; ++i) {
            unsigned preIdx=calPre(forwardAnchor, i, d);
            usedPre[preIdx]=true;
            preRes[i]=preIdx;
        }
        vector<unsigned> preLevel(n, 1);
        vector<unsigned> prePath(n, n);
        for(unsigned i=0; i<n; ++i) {
            prePath[i]=i;
        }
        getPreLevel(preRes, preLevel, usedPre, prePath);

        vector<unsigned> nextRes(n, n);
        vector<bool> usedNext(n, false);
        for(unsigned  i=0; i<n; ++i) {
            unsigned nextIdx=calNext(forwardAnchor, i, d);
            usedNext[nextIdx]=true;
            nextRes[i]=nextIdx;
        }
        vector<unsigned> nextLevel(n, 1);
        vector<unsigned> nextPath(n, n);
        for(unsigned i=0; i<n; ++i) {
            nextPath[i]=i;
        }
        getNextLevel(nextRes, nextLevel, usedNext, nextPath);
        vector<unsigned> preAddNext(n, 0);

//----------------------level------------------------------
        unsigned maxLevel=0;
        vector<unsigned> maxIdx;
        for(unsigned i=0; i<n; ++i) {
            preAddNext[i]=preLevel[i] + nextLevel[i];
            if(preAddNext[i]>maxLevel) {
                maxLevel=preAddNext[i];
                maxIdx.clear();
                maxIdx.emplace_back(i);
            } else if(preAddNext[i]==maxLevel) {
                maxIdx.emplace_back(i);
            }
        }
        cout<<"maxLevel:"<<maxLevel<<endl;
        cout<<"maxIdx.size:"<<maxIdx.size()<<endl;
        vector<vector<int>> forwardRes(maxIdx.size());
        for(unsigned i=0; i<maxIdx.size(); ++i) {
            unsigned temp=maxIdx[i];
            forwardRes[i].emplace_back(temp);
            while(prePath[temp]!=temp) {
                forwardRes[i].emplace_back(prePath[temp]);
                temp=prePath[temp];
            }
            temp=maxIdx[i];
            while(nextPath[temp]!=temp) {
                forwardRes[i].emplace_back(nextPath[temp]);
                temp=nextPath[temp];
            }
            std::sort(forwardRes[i].begin(), forwardRes[i].end());
        }
        cout<<"-----forward-----"<<endl;
        vector<int> forFinRes;
        int forFinSize=0;
        for(unsigned i=0; i<maxIdx.size(); ++i) {
            if(forwardRes[i].size()>forFinSize) {
                forFinSize=forwardRes[i].size();
                forFinRes=forwardRes[i];
            }
        }
//----------------------level------------------------------

//----------------------Adj------------------------------
//        cout<<"maxDegree:"<<maxDegree<<endl;
//        cout<<"maxDegreeIndex.size:"<<maxDegreeIndex.size()<<endl;
//        vector<vector<int>> forwardRes(maxDegreeIndex.size());
//        for(unsigned i=0; i<maxDegreeIndex.size(); ++i) {
//            unsigned temp=maxDegreeIndex[i];
//            forwardRes[i].emplace_back(temp);
//            while(prePath[temp]!=temp) {
//                forwardRes[i].emplace_back(prePath[temp]);
//                temp=prePath[temp];
//            }
//            temp=maxDegreeIndex[i];
//            while(nextPath[temp]!=temp) {
//                forwardRes[i].emplace_back(nextPath[temp]);
//                temp=nextPath[temp];
//            }
//            std::sort(forwardRes[i].begin(), forwardRes[i].end());
//        }
//        cout<<"-----forward-----"<<endl;
//        vector<int> forFinRes;
//        int forFinSize=0;
//        for(unsigned i=0; i<maxDegreeIndex.size(); ++i) {
//            if(forwardRes[i].size()>forFinSize) {
//                forFinSize=forwardRes[i].size();
//                forFinRes=forwardRes[i];
//            }
//        }
//----------------------Adj------------------------------

        cout<<">> ";
        for(int i=0; i<forFinSize; ++i) {
            cout<<forFinRes[i]<<"   ";
        }
        cout<<endl;

//        ----------reverse----------
        vector<unsigned> preResRev(m, m);
        vector<bool> usedPreRev(m, false);
        for(unsigned i=0; i<m; ++i) {
            unsigned preIdx=calPre(reverseAnchor, i, d);
            usedPreRev[preIdx]=true;
            preResRev[i]=preIdx;
        }
        vector<unsigned> preLevelRev(m, 1);
        vector<unsigned> prePathRev(m, m);
        for(unsigned i=0; i<m; ++i) {
            prePathRev[i]=i;
        }
        getPreLevel(preResRev, preLevelRev, usedPreRev, prePathRev);
        vector<unsigned> nextResRev(m, m);
        vector<bool> usedNextRev(m, false);
        for(unsigned  i=0; i<m; ++i) {
            unsigned nextIdx=calNext(reverseAnchor, i, d);
            usedNextRev[nextIdx]=true;
            nextResRev[i]=nextIdx;
        }
        vector<unsigned> nextLevelRev(m, 1);
        vector<unsigned> nextPathRev(m, m);
        for(unsigned i=0; i<m; ++i) {
            nextPathRev[i]=i;
        }
        getNextLevel(nextResRev, nextLevelRev, usedNextRev, nextPathRev);
        vector<unsigned> preAddNextRev(m, 0);

//----------------------level------------------------------
        unsigned maxLevelRev=0;
        vector<unsigned> maxIdxRev;
        for(unsigned i=0; i<m; ++i) {
            preAddNextRev[i]=preLevelRev[i]+nextLevelRev[i];
            if(preAddNextRev[i]>maxLevelRev) {
                maxLevelRev=preAddNextRev[i];
                maxIdxRev.clear();
                maxIdxRev.emplace_back(i);
            } else if(preAddNextRev[i]==maxLevelRev) {
                maxIdxRev.emplace_back(i);
            }
        }
        cout<<"maxLevelRev:"<<maxLevelRev<<endl;
        cout<<"maxIdxRev.size:"<<maxIdxRev.size()<<endl;
        vector<vector<int>> reverseRes(maxIdxRev.size());
        for(unsigned i=0; i<maxIdxRev.size(); ++i) {
            unsigned temp=maxIdxRev[i];
            reverseRes[i].emplace_back(temp);
            while(prePathRev[temp]!=temp) {
                reverseRes[i].emplace_back(prePathRev[temp]);
                temp=prePathRev[temp];
            }
            temp=maxIdxRev[i];
            while(nextPathRev[temp]!=temp) {
                reverseRes[i].emplace_back(nextPathRev[temp]);
                temp=nextPathRev[temp];
            }
            std::sort(reverseRes[i].begin(), reverseRes[i].end());
        }
        cout<<"-----reverse-----"<<endl;
        vector<int> revFinRes;
        int revFinSize=0;
        for(unsigned i=0; i<maxIdxRev.size(); ++i) {
            if(reverseRes[i].size()>revFinSize) {
                revFinSize=reverseRes[i].size();
                revFinRes=reverseRes[i];
            }
        }
//----------------------level------------------------------

//----------------------Adj------------------------------
//        cout<<"maxDegreeRev:"<<maxDegreeReverse<<endl;
//        cout<<"maxDegreeIndexRev.size:"<<maxDegreeIndexReverse.size()<<endl;
//        vector<vector<int>> reverseRes(maxDegreeIndexReverse.size());
//        for(unsigned i=0; i<maxDegreeIndexReverse.size(); ++i) {
//            unsigned temp=maxDegreeIndexReverse[i];
//            reverseRes[i].emplace_back(temp);
//            while(prePathRev[temp]!=temp) {
//                reverseRes[i].emplace_back(prePathRev[temp]);
//                temp=prePathRev[temp];
//            }
//            temp=maxDegreeIndexReverse[i];
//            while(nextPathRev[temp]!=temp) {
//                reverseRes[i].emplace_back(nextPathRev[temp]);
//                temp=nextPathRev[temp];
//            }
//            std::sort(reverseRes[i].begin(), reverseRes[i].end());
//        }
//        cout<<"-----reverse-----"<<endl;
//        vector<int> revFinRes;
//        int revFinSize=0;
//        for(unsigned i=0; i<maxDegreeIndexReverse.size(); ++i) {
//            if(reverseRes[i].size()>revFinSize) {
//                revFinSize=reverseRes[i].size();
//                revFinRes=reverseRes[i];
//            }
//        }
//----------------------Adj------------------------------

        cout<<">> ";
        for(int i=0; i<revFinSize; ++i) {
            cout<<revFinRes[i]<<"   ";
        }
        cout<<endl;

//        revFinRes
//        forFinRes
        BaseLevelPos OneRefRes;
        OneRefRes.AnchorNum=n;
        OneRefRes.AnchorNumRev=m;
//        cout<<"11"<<endl;
        if(forFinRes.size()>=revFinRes.size()) {
            OneRefRes.RefName=OneRef.first;
            OneRefRes.dir=true;
//            cout<<"22"<<endl;
            for(int ii=0; ii<forFinRes.size(); ++ii) {
                AlignPosition temp;
                temp.ReadPos=forwardAnchor[forFinRes[ii]].StartPosInRead;
                temp.RefPos=forwardAnchor[forFinRes[ii]].StartPosInRef;
                OneRefRes.position.emplace_back(temp);
            }


            if(OneRefRes.position.size()<3 && forwardAnchor.size()>10) {
                ofstream spOut("spOut.txt", ios::app);
                spOut<<"----------------------------------"<<endl;
                spOut<<"forward"<<endl;
                spOut<<OneRefRes.position.size()<<"     and     "<<forwardAnchor.size()<<endl;
                for(int i=0; i<forwardAnchor.size(); ++i) {
                    spOut<<forwardAnchor[i].StartPosInRead<<"    "<<forwardAnchor[i].StartPosInRef<<endl;
                }
                spOut<<endl;
            }


//            cout<<"33"<<endl;
        } else {
            OneRefRes.RefName=OneRef.first;
            OneRefRes.dir=false;
//            cout<<"44"<<endl;
            for(int ii=0; ii<revFinRes.size(); ++ii) {
                AlignPosition temp;
                temp.ReadPos=reverseAnchor[revFinRes[ii]].StartPosInRead;
                temp.RefPos=reverseAnchor[revFinRes[ii]].StartPosInRef;
                OneRefRes.position.emplace_back(temp);
            }


            if(OneRefRes.position.size()<3 && reverseAnchor.size()>10) {
                ofstream spOut("spOut.txt", ios::app);
                spOut<<"----------------------------------"<<endl;
                spOut<<"reverse"<<endl;
                spOut<<OneRefRes.position.size()<<"     and     "<<reverseAnchor.size()<<endl;
                for(int i=0; i<reverseAnchor.size(); ++i) {
                    spOut<<reverseAnchor[i].StartPosInRead<<"    "<<reverseAnchor[i].StartPosInRef<<endl;
                }
                spOut<<endl;
            }


//            cout<<"55"<<endl;
        }
//        cout<<"66"<<endl;

        if(Best5.size()<5) {
            Best5.push(OneRefRes);
        } else if(OneRefRes.position.size()>Best5.top().position.size()) {
            Best5.pop();
            Best5.push(OneRefRes);
        }

//        if(OneReadRes.position.size()<OneRefRes.position.size()) {
//            OneReadResSec=OneReadRes;
//            OneReadRes=OneRefRes;
//        } else if(OneReadResSec.position.size()<OneRefRes.position.size()) {
//            OneReadResSec=OneRefRes;
//        }
//        cout<<"77"<<endl;
    }
    OneReadResAll.resize(Best5.size());
    for(int i=Best5.size()-1; i>=0; --i) {
        OneReadResAll[i]=Best5.top();
        Best5.pop();
    }
//
//    vector<int> store;
//    ofstream zz0("outAnchInfo.txt", ios::app);
//    vector<BaseLevelPos> AlignPosArr;
//    AlignPosArr=OneReadResAll;
//
//    vector<double> chainBest5(AlignPosArr.size());
//    chainBest5[0]=100;
//    double ave=100;
//    for(int i=1; i<chainBest5.size(); ++i) {
//        chainBest5[i]=(double)AlignPosArr[i].position.size() / (double)AlignPosArr[0].position.size() * 100;
//        ave+=chainBest5[i];
//    }
//    ave/=5;
//    double zz=0.0;
//    for(int i=0; i<chainBest5.size(); ++i) {
//        zz+=pow(chainBest5[i]-ave, 2);
//    }
//    for(int i=chainBest5.size(); i<5; ++i) {
//        zz+=pow(ave, 2);
//    }
//
//    zz=sqrt(zz)/4;
//
//    int choseFLAG=0;
//    int maxAnchNum=0;
//    if(zz==0) {
//        for(int i=0; i<chainBest5.size(); ++i) {
//            zz0<<AlignPosArr[i].position.size()<<" dir"<<AlignPosArr[i].dir<<
//               "("<<AlignPosArr[i].AnchorNum<<"+"<<AlignPosArr[i].AnchorNumRev<<")"<<"  "<<endl;
//            if(true) {
//                vector<AnchorPoint> forwardAnchor;
//                vector<AnchorPoint> reverseAnchor;
//                for (auto Anchor: UseRefNameToDivide[AlignPosArr[i].RefName]) {
//                    if (Anchor.isReverse) {
//                        reverseAnchor.emplace_back(Anchor);
//                    } else {
//                        forwardAnchor.emplace_back(Anchor);
//                    }
//                }
//                int n = forwardAnchor.size();
//                int m = reverseAnchor.size();
//
//                sort(forwardAnchor.begin(), forwardAnchor.end(), ChainCmpByRefPos());
//                sort(reverseAnchor.begin(), reverseAnchor.end(), ChainCmpByRefPos());
//                zz0<<"forward"<<endl;
//                for(int i=0; i<n; ++i) {
//                    zz0<<forwardAnchor[i].StartPosInRead<<"     "<<forwardAnchor[i].StartPosInRef<<endl;
//                }
//                zz0<<"reverse"<<endl;
//                for(int i=0; i<m; ++i) {
//                    zz0<<reverseAnchor[i].StartPosInRead<<"     "<<reverseAnchor[i].StartPosInRef<<endl;
//                }
//            }
//        }
//        zz0<<endl;
//        zz0<<endl;
//        zz0<<endl;
//    }
//    store.emplace_back(zz);


//    OneReadResAll.emplace_back(OneReadRes);
//    OneReadResAll.emplace_back(OneReadResSec);
    return OneReadResAll;
}

int calPre(vector<AnchorPoint>& v, int idx, unsigned d) {
    //find the idx of the last fragment
    if(idx==0){
        return idx;
    }
    int pos_tmp=0;
    int i;
    unsigned r=5000;    // 500 1000 1500
    for(i=idx-1; i>=0; --i) {
        if(v[i].StartPosInRead!=v[idx].StartPosInRead &&
                ((!v[i].isReverse && v[i].StartPosInRead<v[idx].StartPosInRead) ||
                 (v[i].isReverse && v[i].StartPosInRead>v[idx].StartPosInRead))) {
            pos_tmp=v[i].StartPosInRef;
            if(abs(abs(v[i].StartPosInRead-v[i].StartPosInRef)-abs(v[idx].StartPosInRead-v[idx].StartPosInRef))>r){
                i=-1;
            }
            break;
        }
    }
    if(pos_tmp>d) {
        pos_tmp-=d;
    } else {
        pos_tmp=0;
    }

    int rr=i>=0?i:v.size();
    for(int j=i; j>=0 && v[j].StartPosInRef>=pos_tmp; --j) {
        unsigned tmp=abs(abs(v[j].StartPosInRead-v[idx].StartPosInRead)-abs(v[j].StartPosInRef-v[idx].StartPosInRef));
        if(tmp<r &&
                ((!v[j].isReverse && v[j].StartPosInRead<v[idx].StartPosInRead) ||
                 (v[j].isReverse && v[j].StartPosInRead>v[idx].StartPosInRead))) {
            r=tmp;
            rr=j;
        }
    }
    return rr;
}

int calNext(vector<AnchorPoint>& v, int idx, unsigned d) {
    int n=v.size();
    if(idx==n-1){
        return idx;
    }
    unsigned pos_tmp=0;
    unsigned r=5000;
    uint32_t i;
    for(i=idx+1; i<n; ++i) {
        if(v[i].StartPosInRead!=v[idx].StartPosInRead &&
                ((!v[i].isReverse && v[i].StartPosInRead>v[idx].StartPosInRead) ||
                 (v[i].isReverse && v[i].StartPosInRead<v[idx].StartPosInRead))) {
            pos_tmp=v[i].StartPosInRef;
            if(abs(abs(v[i].StartPosInRead-v[i].StartPosInRef)-abs(v[idx].StartPosInRead-v[idx].StartPosInRef))>r){
                i=n;
            }
            break;
        }
    }

    pos_tmp+=d;
//    if(pos_tmp>d) {
//        pos_tmp-=d;
//    }
//    else {
//        pos_tmp=0;
//    }


    int rr=i;
    for(uint32_t j=i; j<n && v[j].StartPosInRef<=pos_tmp; ++j) {
        unsigned tmp=abs(abs(v[j].StartPosInRead-v[idx].StartPosInRead)-abs(v[j].StartPosInRef-v[idx].StartPosInRef));
        if(tmp<r &&
                ((!v[j].isReverse && v[j].StartPosInRead>v[idx].StartPosInRead) ||
                 (v[j].isReverse && v[j].StartPosInRead<v[idx].StartPosInRead))) {
            r=tmp;
            rr=j;
        }
    }

    return rr;
}

void getPreLevel(vector<unsigned>& preRes, vector<unsigned>& preLevel, vector<bool>& usedPre, vector<unsigned>& prePath) {
    int n = preLevel.size();
    for(int i = n - 1; i > 0; --i) {
        if(preLevel[i] + 1 > preLevel[preRes[i]] && preRes[i] != n) {
            preLevel[preRes[i]] = preLevel[i] + 1;
            prePath[preRes[i]] = i;
        }

//        if(usedPre[i]) {
//            continue;
//        }
//        int idx=i;
//        while(preRes[idx]!=idx) {
//            if(preLevel[preRes[idx]]<preLevel[idx]+1) {
//                preLevel[preRes[idx]]=preLevel[idx]+1;
//                idx=preRes[idx];
//                prePath[preRes[idx]]=idx;
//            } else {
//                break;
//            }
//        }
    }
}

void getNextLevel(vector<unsigned>& nextRes, vector<unsigned>& nextLevel, vector<bool>& usedNext, vector<unsigned>& nextPath) {
    int n=nextLevel.size();
    for(int i=0; i<n-1; ++i) {
        if(nextLevel[i]+1>nextLevel[nextRes[i]] && nextRes[i]!=n) {
            nextLevel[nextRes[i]]=nextLevel[i]+1;
            nextPath[nextRes[i]]=i;
        }

//        if(usedNext[i]) {
//            continue;
//        }
//        int idx=i;
//        while(nextRes[idx]!=idx) {
//            if(nextLevel[nextRes[idx]]<nextLevel[idx]+1) {
//                nextLevel[nextRes[idx]]=nextLevel[idx]+1;
//                idx=nextRes[idx];
//                nextPath[nextRes[idx]]=idx;
//            } else {
//                break;
//            }
//        }
    }
}

BaseLevelPos chain(vector<SplitReadFragment> OneRead, int SplitLength, int e) {
    unordered_map<string, vector<AnchorPoint>> UseRefNameToDivide;
//    vector<vector<AnchorPoint>> ChainResultSet;
//    ofstream clearMaxClique("testMaxClique.txt", std::ios::out | std::ios::trunc);
//    clearMaxClique.close();
//    ofstream clearFinal("testFinal.txt", std::ios::out | std::ios::trunc);
//    clearFinal.close();
//    ofstream fout("testMaxClique.txt", ios::app);
//    ofstream foutfinal("testFinal.txt", ios::app);
//    fout << "-----------------oneread-----------------" << endl;
//    foutfinal << "-----------------oneread-----------------" << endl;
    for (int i = 0; i < OneRead.size(); ++i) {
        for (int j = 0; j < OneRead[i].SplitReadMapResult.size(); ++j) {
            AnchorPoint anchor;
//            if (((stoi(OneRead[i].SplitReadMapResult[j].MapFlag) >> 3) & 1) ||
//                ((stoi(OneRead[i].SplitReadMapResult[j].MapFlag) >> 4) & 1)) {
//                anchor.isReverse = true;
//            }
            if (stoi(OneRead[i].SplitReadMapResult[j].MapFlag) & 16) {
                anchor.isReverse = true;
            }
            anchor.segment = SplitLength;
            anchor.StartPosInRead = (int)(i * SplitLength * 0.5);
            anchor.StartPosInRef = stoi(OneRead[i].SplitReadMapResult[j].MapPosition);
            anchor.ReadPosSubRefPos = anchor.StartPosInRead - anchor.StartPosInRef;
            anchor.RefName = OneRead[i].SplitReadMapResult[j].RefName;
            UseRefNameToDivide[anchor.RefName].emplace_back(anchor);
        }
    }
    BaseLevelPos OneReadRes;
    for (auto OneRef: UseRefNameToDivide) {
//        fout << OneRef.first << endl;
        vector<AnchorPoint> forwardAnchor;
        vector<AnchorPoint> reverseAnchor;
        for (auto Anchor: OneRef.second) {
            if (Anchor.isReverse) {
                reverseAnchor.emplace_back(Anchor);
            } else {
                forwardAnchor.emplace_back(Anchor);
            }
        }
        int n = forwardAnchor.size();
        int m = reverseAnchor.size();
        sort(forwardAnchor.begin(), forwardAnchor.end(), ChainCmpByRefPos());
        sort(reverseAnchor.begin(), reverseAnchor.end(), ChainCmpByRefPos());
        vector<vector<int>> AdjList(n, vector<int>());
        vector<vector<int>> AdjListReverse(m, vector<int>());
        vector<int> indegree(n, 0);
        vector<int> indegreeReverse(m, 0);
        // AdjList
        for (int i = 0; i < n; ++i) {
            AdjList[i].emplace_back(i);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (checkConnected(forwardAnchor[i], forwardAnchor[j], e)) {
                    ++indegree[j];
                    AdjList[i].emplace_back(j);
                }
            }
        }
        // AdjListReverse
        for (int i = m - 1; i >= 0; --i) {
            AdjListReverse[i].emplace_back(i);
        }
        for (int i = m - 1; i >= 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                if (checkConnected(reverseAnchor[i], reverseAnchor[j], e)) {
                    ++indegreeReverse[j];
                    AdjListReverse[i].emplace_back(j);
                }
            }
        }

        vector<Chain> chainResult;
        vector<Chain> chainResultRe;
        vector<int> indegree0;
        vector<int> indegree0Reverse;
        vector<int> originUsed(n, 0);
        vector<int> originUsedReverse(m, 0);
        for (int i = 0; i < n; ++i) {
            if (indegree[i] == 0) {
                indegree0.emplace_back(i);
                originUsed[i] = 1;
            }
        }

        for (int i = m - 1; i >= 0; --i) {
            if (indegreeReverse[i] == 0) {
                indegree0Reverse.emplace_back(i);
                originUsedReverse[i] = 1;
            }
        }

        // ------------------------forward---------------------------
        // find chain
        for (int i = 0; i < indegree0.size(); ++i) {
            vector<int> used = originUsed;
            vector<int> temp;
            temp.emplace_back(indegree0[i]);
            dfsToFindChain(AdjList, temp, chainResult, used, indegree0[i]);
        }
//        fout << "forward" << endl;;

        // chain Info
        for (int i = 0; i < chainResult.size(); ++i) {
            chainResult[i].segment = SplitLength;
            chainResult[i].length = chainResult[i].chainSequence.size();
            chainResult[i].beginPosition = forwardAnchor[chainResult[i].chainSequence[0]].StartPosInRef;
            chainResult[i].endPosition = forwardAnchor[chainResult[i].chainSequence[chainResult[i].length -
                                                                                    1]].StartPosInRef;
            chainResult[i].shift = forwardAnchor[chainResult[i].chainSequence[0]].ReadPosSubRefPos;
            for (int j = 1; j < chainResult[i].length; ++j) {
                chainResult[i].totalDiff += abs(forwardAnchor[chainResult[i].chainSequence[j]].ReadPosSubRefPos -
                                                forwardAnchor[chainResult[i].chainSequence[j - 1]].ReadPosSubRefPos);
            }
        }
        // vector<Chain> chainResult;
        int chainNumber = chainResult.size();
        int edgeNumber = 0;
        vector<vector<int>> MP(chainNumber + 1, vector<int>(chainNumber + 1));
        for (int i = 0; i < chainNumber; ++i) {
            for (int j = i + 1; j < chainNumber; ++j) {
                if (abs(chainResult[i].shift - chainResult[j].shift) < 10 * e) {
                    MP[i + 1][j + 1] = MP[j + 1][i + 1] = 1;
                    ++edgeNumber;
                }
            }
        }


        vector<Clique> maximalCliqueResult = maximalClique(chainNumber, edgeNumber, MP);
//        for (int i = 0; i < maximalCliqueResult.size(); ++i) {
//            fout << "one clique" << maximalCliqueResult[i].maximalClique.size() << endl;
//            for (int j = 0; j < maximalCliqueResult[i].maximalClique.size(); ++j) {
//                for (int k = 0; k < chainResult[j].chainSequence.size(); ++k) {
//                    fout << chainResult[j].chainSequence[k] << "    ";
//                }
//                fout << endl;
//            }
//        }

        Clique ForwardMaxScoreClique = maxClique(maximalCliqueResult, chainResult, forwardAnchor, true);
//        if (ForwardMaxScoreClique.score != 0) {
//            fout << "--ForwardMaxScoreClique--" << "    " << ForwardMaxScoreClique.maximalClique.size() << "    "
//                 << ForwardMaxScoreClique.score << endl;
//        }
//
//        for (int i = 0; i < ForwardMaxScoreClique.maximalClique.size(); ++i) {
//            for (int j = 0; j < chainResult[ForwardMaxScoreClique.maximalClique[i]].chainSequence.size(); ++j) {
//                fout << chainResult[ForwardMaxScoreClique.maximalClique[i]].chainSequence[j] << "    ";
//            }
//            fout << endl;
//        }

        // ------------------------reverse---------------------------
        for (int i = 0; i < indegree0Reverse.size(); ++i) {
            vector<int> used = originUsedReverse;
            vector<int> temp;
            temp.emplace_back(indegree0Reverse[i]);
            dfsToFindChain(AdjListReverse, temp, chainResultRe, used, indegree0Reverse[i]);
        }
//        fout << "reverse" << endl;

        // chain Info
        for (int i = 0; i < chainResultRe.size(); ++i) {
            chainResultRe[i].segment = SplitLength;
            chainResultRe[i].length = chainResultRe[i].chainSequence.size();
            chainResultRe[i].beginPosition = reverseAnchor[chainResultRe[i].chainSequence[0]].StartPosInRef;
            chainResultRe[i].endPosition = reverseAnchor[chainResultRe[i].chainSequence[chainResultRe[i].length -
                                                                                        1]].StartPosInRef;
            chainResultRe[i].shift = reverseAnchor[chainResultRe[i].chainSequence[0]].ReadPosSubRefPos;
            for (int j = 1; j < chainResultRe[i].length; ++j) {
                chainResultRe[i].totalDiff += abs(reverseAnchor[chainResultRe[i].chainSequence[j]].ReadPosSubRefPos -
                                                  reverseAnchor[chainResultRe[i].chainSequence[j -
                                                                                               1]].ReadPosSubRefPos);
            }
        }

        int chainNumberRe = chainResultRe.size();
        int edgeNumberRe = 0;
        vector<vector<int>> MPRe(chainNumberRe + 1, vector<int>(chainNumberRe + 1));
        for (int i = 0; i < chainNumberRe; ++i) {
            for (int j = i + 1; j < chainNumberRe; ++j) {
                if (abs(chainResultRe[i].shift - chainResultRe[j].shift) < 10 * e) {
                    MPRe[i + 1][j + 1] = MPRe[j + 1][i + 1] = 1;
                    ++edgeNumberRe;
                }
            }
        }

        vector<Clique> maximalCliqueResultRe = maximalClique(chainNumberRe, edgeNumberRe, MPRe);
//        for (int i = 0; i < maximalCliqueResultRe.size(); ++i) {
//            fout << "one clique" << maximalCliqueResultRe[i].maximalClique.size() << endl;
//            for (int j = 0; j < maximalCliqueResultRe[i].maximalClique.size(); ++j) {
//                for (int k = 0; k < chainResultRe[j].chainSequence.size(); ++k) {
//                    fout << chainResultRe[j].chainSequence[k] << "    ";
//                }
//                fout << endl;
//            }
//        }

        Clique ReverseMaxScoreClique = maxClique(maximalCliqueResultRe, chainResultRe, reverseAnchor, false);
//        if (ReverseMaxScoreClique.score != 0) {
//            fout << "--ReverseMaxScoreClique--" << "    " << ReverseMaxScoreClique.maximalClique.size() << "    "
//                 << ReverseMaxScoreClique.score << endl;
//        }
//
//        for (int i = 0; i < ReverseMaxScoreClique.maximalClique.size(); ++i) {
//            for (int j = 0; j < chainResultRe[ReverseMaxScoreClique.maximalClique[i]].chainSequence.size(); ++j) {
//                fout << chainResultRe[ReverseMaxScoreClique.maximalClique[i]].chainSequence[j] << "    ";
//            }
//            fout << endl;
//        }

        BaseLevelPos OneRefRes;
        if (ForwardMaxScoreClique.score >= ReverseMaxScoreClique.score) {
            OneRefRes = getCliqueInfo(ForwardMaxScoreClique, chainResult, forwardAnchor, true, OneRef.first);
        } else {
            OneRefRes = getCliqueInfo(ReverseMaxScoreClique, chainResultRe, reverseAnchor, false, OneRef.first);
        }
        if (OneReadRes.score < OneRefRes.score) {
            OneReadRes = OneRefRes;
        }
    }
//    fout.close();
//
//    foutfinal << "RefName: " << OneReadRes.RefName << endl;
//    foutfinal << "Score: " << OneReadRes.score << endl;
//    foutfinal << "dir: " << OneReadRes.dir << endl;
//    foutfinal << "position number: " << OneReadRes.position.size() << endl;
//    for (int i = 0; i < OneReadRes.position.size(); ++i) {
//        foutfinal << OneReadRes.position[i].RefPos << " " << OneReadRes.position[i].ReadPos << "        ";
//    }
//    foutfinal << endl;
//    foutfinal.close();

    return OneReadRes;
}

bool checkConnected(AnchorPoint a, AnchorPoint b, int e) {
    if (a.StartPosInRead < b.StartPosInRead && abs(a.ReadPosSubRefPos - b.ReadPosSubRefPos) < e) {
        return true;
    } else {
        return false;
    }
}

void dfsToFindChain(vector<vector<int>> &AdjList, vector<int> &temp, vector<Chain> &chainResult, vector<int> &used,
                    int index) {
    if (AdjList[index].size() == 1) {
        Chain chainTemp;
        chainTemp.chainSequence = temp;
        chainResult.emplace_back(chainTemp);
        return;
    }
    for (int i = 1; i < AdjList[index].size(); ++i) {
        if (used[AdjList[index][i]] == 1) {
            continue;
        } else {
            used[AdjList[index][i]] = 1;
            temp.emplace_back(AdjList[index][i]);
            dfsToFindChain(AdjList, temp, chainResult, used, AdjList[index][i]);
            temp.pop_back();
        }
    }
}

vector<Clique> maximalClique(int n, int m, vector<vector<int>> MP) {
    // MP:connectivity
    vector<vector<int>> some(n + 100, vector<int>(n + 100));
    vector<vector<int>> none(n + 100, vector<int>(n + 100));
    vector<vector<int>> all(n + 100, vector<int>(n + 100));
    vector<Clique> res;
    int ans = 0;
    for (int i = 0; i < n + 100; ++i) {
        some[1][i] = i + 1;
    }
    dfsToMaximalClique(1, 0, n, 0, some, none, all, MP, ans, n, m, res);
    return res;
}

void dfsToMaximalClique(int d, int an, int sn, int nn, vector<vector<int>> &some, vector<vector<int>> &none,
                        vector<vector<int>> &all, vector<vector<int>> MP, int ans, int n, int m, vector<Clique> &res) {
    if (!sn && !nn) {
        ++ans;
        Clique temp;
        for (int i = 0; i < 100; ++i) {
            if (all[d][i] != 0) {
                temp.maximalClique.emplace_back(all[d][i] - 1);
            }
        }
        if (temp.maximalClique.size() != 0) {
            res.emplace_back(temp);
        }
    }
    int u = some[d][0];  //选取Pivot结点
    for (int i = 0; i < sn; ++i) {
        int v = some[d][i];
        if (MP[u][v]) continue;
        //如果是邻居结点，就直接跳过下面的程序，进行下一轮的循环。显然能让程序运行下去的，只有两种，一种是v就是u结点本身，另一种是v不是u的邻居结点。
        for (int j = 0; j < an; ++j) {
            all[d + 1][j] = all[d][j];
        }
        all[d + 1][an] = v;
        int tsn = 0, tnn = 0;
        for (int j = 0; j < sn; ++j) if (MP[v][some[d][j]]) some[d + 1][tsn++] = some[d][j];
        for (int j = 0; j < nn; ++j) if (MP[v][none[d][j]]) none[d + 1][tnn++] = none[d][j];
        dfsToMaximalClique(d + 1, an + 1, tsn, tnn, some, none, all, MP, ans, n, m, res);
        some[d][i] = 0, none[d][nn++] = v;
    }
}

Clique maxClique(vector<Clique> cliques, vector<Chain> chainResult, vector<AnchorPoint> anchors, bool dir) {
    int n = cliques.size();
    Clique res;
    for (int i = 0; i < n; ++i) {
        Clique curr = calculateCliqueScore(cliques[i], chainResult, anchors, dir);
        if (curr.score > res.score) {
            res = curr;
        }
    }
    return res;
}

Clique calculateCliqueScore(Clique &clique, vector<Chain> chainResult, vector<AnchorPoint> anchors, bool dir) {
    vector<vector<long>> cover = sortClique(clique, chainResult, dir);
    Clique res;
    res.score = 0.0;
    int n = clique.maximalClique.size();
    for (int i = 0; i < n; ++i) {
        checkChainCanAddClique(i, cover, chainResult, res, dir);
    }
    for (int i = 0; i < res.maximalClique.size(); ++i) {
        res.maximalClique[i] = clique.maximalClique[res.maximalClique[i]];
    }
    return res;
}

vector<vector<long>> sortClique(Clique &clique, vector<Chain> chainResult, bool dir) {
    int n = clique.maximalClique.size();
    vector<vector<long>> cover(n, vector<long>(3));
    for (int i = 0; i < n; ++i) {
        if (dir) {
            cover[i][0] = chainResult[clique.maximalClique[i]].beginPosition;
            cover[i][1] =
                    chainResult[clique.maximalClique[i]].endPosition + chainResult[clique.maximalClique[i]].segment;
        } else {
            cover[i][0] = chainResult[clique.maximalClique[i]].endPosition;
            cover[i][1] =
                    chainResult[clique.maximalClique[i]].beginPosition + chainResult[clique.maximalClique[i]].segment;
        }
        cover[i][2] = abs(chainResult[clique.maximalClique[i]].endPosition -
                          chainResult[clique.maximalClique[i]].beginPosition) +
                      chainResult[clique.maximalClique[i]].segment;
    }
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (cover[j][2] < cover[j + 1][2]) {
                swap(clique.maximalClique[j], clique.maximalClique[j + 1]);
                swap(cover[j], cover[j + 1]);
            }
        }
    }
    return cover;
}

bool checkChainCanAddClique(int chain, vector<vector<long>> cover, vector<Chain> chainResult, Clique &temp, bool dir) {
    temp.maximalClique.emplace_back(chain);
    int n = temp.maximalClique.size();

    vector<vector<long>> totalCover;
    totalCover.emplace_back(cover[temp.maximalClique[0]]);
    for (int i = 1; i < n; ++i) {
        bool flag = true;
        for (int j = 0; j < totalCover.size(); ++j) {
            if (cover[temp.maximalClique[i]][0] >= totalCover[j][1] ||
                cover[temp.maximalClique[i]][1] <= totalCover[j][0]) {
                continue;
            } else if (cover[temp.maximalClique[i]][0] >= totalCover[j][0] &&
                       cover[temp.maximalClique[i]][1] <= totalCover[j][1]) {
                flag = false;
                break;
            } else if (cover[temp.maximalClique[i]][0] >= totalCover[j][0] &&
                       cover[temp.maximalClique[i]][1] >= totalCover[j][1]) {
                totalCover[j][1] = cover[temp.maximalClique[i]][1];
                flag = false;
            } else if (cover[i][0] <= totalCover[j][0] && cover[i][1] <= totalCover[j][1]) {
                totalCover[j][0] = cover[temp.maximalClique[i]][0];
                flag = false;
            }
        }
        if (flag) {
            totalCover.emplace_back(cover[temp.maximalClique[i]]);
        }
    }

    // if already calculate cover not overlap
    long coverLen = 0;
    for (int i = 0; i < totalCover.size(); ++i) {
        coverLen += abs(totalCover[i][1] - totalCover[i][0]);
    }
    double newScore = 1.0 * coverLen / (temp.maximalClique.size() + 50);
    if (newScore > temp.score) {
        temp.score = newScore;
        return true;
    } else {
        temp.maximalClique.pop_back();
        return false;
    }

}

BaseLevelPos getCliqueInfo(Clique maxScoreClique, vector<Chain> chainResult, vector<AnchorPoint> anchors, bool dir, string RefName) {
    BaseLevelPos res;
    res.dir = dir;
    res.score = maxScoreClique.score;
    res.RefName = RefName;
    unordered_set<int> pointSet;
    for (int i = 0; i < maxScoreClique.maximalClique.size(); ++i) {
        for (int j = 0; j < chainResult[maxScoreClique.maximalClique[i]].chainSequence.size(); ++j) {
            pointSet.emplace(chainResult[maxScoreClique.maximalClique[i]].chainSequence[j]);
        }
    }
    for (int anchorIndex: pointSet) {
        AlignPosition temp;
        temp.ReadPos = anchors[anchorIndex].StartPosInRead;
        temp.RefPos = anchors[anchorIndex].StartPosInRef;
        res.position.emplace_back(temp);
    }
    sort(res.position.begin(), res.position.end(), PosCmpByReadPos());
    if(dir) {
        res.position=LIS(res.position);
    } else {
        res.position=LDS(res.position);
    }
    return res;
}

vector<AlignPosition> LIS(vector<AlignPosition>& arr) {
    // write code here
    // 第一步：利用贪心+二分求最长递增子序列长度
    vector<AlignPosition> res;//组里面存放递增子序列
    vector<int> maxLen;
    if (arr.size() < 1) return res;
    res.emplace_back(arr[0]);  // 注：emplace_back(val)作用同push_back，效率更高
    maxLen.emplace_back(1);
    for (int i = 1; i < arr.size(); ++i) {
        if (arr[i].ReadPos > res.back().ReadPos) {
            res.emplace_back(arr[i]);
            maxLen.emplace_back(res.size());
        } else {
            // lower_bound(begin, end, val)包含在<algorithm>中
            // 它的作用是返回有序数组begin..end中第一个大于等于val的元素的迭代器
//            int pos = lower_bound(res.begin(), res.end(), ResPosCmp) - res.begin();
            int pos=res.size();
            for(int ii=0; ii<res.size(); ++ii) {
                if(res[ii].ReadPos>=arr[i].ReadPos) {
                    pos=ii;
                    break;
                }
            }
            if(pos==res.size()) {
                res.emplace_back(arr[i]);
            } else {
                res[pos] = arr[i];
            }
            maxLen.emplace_back(pos + 1);
        }
    }
    // 第二步：填充最长递增子序列
    for (int i = arr.size() - 1, j = res.size(); j > 0; --i) {
        if (maxLen[i] == j) {
            res[--j] = arr[i];
        }
    }
    return res;
}

vector<AlignPosition> LDS(vector<AlignPosition>& arr) {
    // write code here
    // 第一步：利用贪心+二分求最长递增子序列长度
    vector<AlignPosition> res;//组里面存放递增子序列
    vector<int> maxLen;
    if (arr.size() < 1) return res;
    res.emplace_back(arr[0]);  // 注：emplace_back(val)作用同push_back，效率更高
    maxLen.emplace_back(1);
    for (int i = 1; i < arr.size(); ++i) {
        if (arr[i].ReadPos < res.back().ReadPos) {
            res.emplace_back(arr[i]);
            maxLen.emplace_back(res.size());
        } else {
            // lower_bound(begin, end, val)包含在<algorithm>中
            // 它的作用是返回有序数组begin..end中第一个大于等于val的元素的迭代器
//            int pos = lower_bound(res.begin(), res.end(), ResPosCmp) - res.begin();
            int pos=res.size();
            for(int ii=0; ii<res.size(); ++ii) {
                if(res[ii].ReadPos<=arr[i].ReadPos) {
                    pos=ii;
                    break;
                }
            }
            if(pos==res.size()) {
                res.emplace_back(arr[i]);
            } else {
                res[pos] = arr[i];
            }
            maxLen.emplace_back(pos + 1);
        }
    }
    // 第二步：填充最长递增子序列
    for (int i = arr.size() - 1, j = res.size(); j > 0; --i) {
        if (maxLen[i] == j) {
            res[--j] = arr[i];
        }
    }
    return res;
}

uint32_t getIandM(string s) {
    long long CCount=0;
    char op=' ';
    uint32_t num=0;
    for(uint32_t i=0; i<s.size(); ++i) {
        if((int)s[i]>=48 && (int)s[i]<=57) {
            num=num*10+(int)s[i]-48;
        } else {
            op=s[i];
            if(op=='M' || op=='I') {
                CCount+=num;
            }
            num=0;
            op=' ';
        }
    }
    return CCount;
}

vector<uint32_t> handleCIGAR(string BLresult, uint32_t ReadLen, int DIFF) {
    long CCount=0;
    vector<uint32_t> res;
    uint32_t num=0;
    char op=' ';
    uint32_t prenum=0;
    char preop='_';
    for(uint32_t i=0; i<BLresult.size(); ++i) {
        if((int)BLresult[i]>=48 && (int)BLresult[i]<=57) {
            num=num*10+(int)BLresult[i]-48;
        } else {
            op=BLresult[i];
            if(op=='M' || op=='I' || op=='S') {
                CCount+=num;
            }
            if(op==preop) {
                num+=prenum;
                res.pop_back();
            }
            prenum=num;
            preop=op;
            num<<=4;
            switch(op) {
                case 'M':
                    num |= 0;
                    break;
                case 'I':
                    num |= 1;
                    break;
                case 'D':
                    num |= 2;
                    break;
                case 'S':
                    num |= 4;
                    break;
            }
            res.emplace_back(num);
            num=0;
            op=' ';
        }
    }
    if(CCount!=ReadLen) {
        ++DIFF;
    }

    return res;
}

EdlibAlignResult edlibBLCommon(char* read, char* ref) {
    EdlibAlignResult result = edlibAlign(read, strlen(read), ref, strlen(ref), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
//        printf("%d\n", result.editDistance);
//        printf("%d\n", result.alignmentLength);
//        printf("%d\n", result.endLocations[0]);
        return result;
    } else {
        cout<<"error?"<<endl;
        int errorTemp;
        cin>>errorTemp;
    }
//    char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
//    free(cigar);
//    edlibFreeAlignResult(result);
}

EdlibAlignResult edlibBLHeadOrTail(char* read, char* ref) {
    EdlibAlignResult result = edlibAlign(read, strlen(read), ref, strlen(ref), edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
//        printf("%d\n", result.editDistance);
//        printf("%d\n", result.alignmentLength);
//        printf("%d\n", result.endLocations[0]);
        return result;
    } else {
        cout<<"error?"<<endl;
        int errorTemp;
        cin>>errorTemp;
    }
}

string reverseCIGAR(char* cigar) {
    string cigarTemp=cigar;
    stack<string> stk;
    string res="";
    uint32_t num=0;
    char op=' ';

    for(uint32_t i=0; i<cigarTemp.size(); ++i) {
        if((int)cigarTemp[i]>=48 && (int)cigarTemp[i]<=57) {
            num=num*10+(int)cigar[i]-48;
        } else {
            op=cigarTemp[i];
            string temp=to_string(num)+op;
            stk.push(temp);
            num=0;
            op=' ';
        }
    }
    while(!stk.empty()) {
        res+=stk.top();
        stk.pop();
    }
    return res;
}

#pragma clang diagnostic pop