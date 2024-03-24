#include <iostream>
#include <string>
#include <stdlib.h>

#include "Basic.h"

using namespace std;

int main(int argc,char *argv[]) {

    string readPath="pb_human.fasta";
    string refPath="05_hg38.fa";
    string outPath="test.sam";
    for(int32_t i=1;i<argc;i=i+2) {
        if (argv[i][0] == '-' && argv[i][1] == 'i') {
            readPath = argv[i + 1];
        }
        if (argv[i][0] == '-' && argv[i][1] == 'I') {
            refPath = argv[i + 1];
        }
        if (argv[i][0] == '-' && argv[i][1] == 'o') {
            outPath = argv[i + 1];
        }
    }

    string refPre="";
    for(int i=0; i<refPath.size(); ++i) {
        if(refPath[i]!='.') {
            refPre+=refPath[i];
        } else {
            break;
        }
    }


    SplitReadSet(0.1, readPath);
    cout << "Sample read success" << endl;

//    //minimap2
//    string Minimap2="./minimap2 -a ";
//    Minimap2+=refPath;
//    Minimap2+=" readSample.fasta > readSample.sam";
//    system(Minimap2.c_str());
//    cout << "Minimap2 success" << endl;

    //samtools
    double errorRate = SamtoolsToGetErrorRate();
    cout << "Samtools to get errorRate success" << endl;

    //Calculate read segment length
    double SplitLength = CalReadSegLen(errorRate, 0.12, 0.3);
    cout << "Calculate read segment length success" << endl;

    //Split read
    SplitLength=50;
    cout<<"SplitLen:"<<SplitLength<<endl;
    unordered_map<string, vector<SplitReadFragment>> SplitReadInfo = SplitRead(SplitLength, readPath);
    cout << "Split read success" << endl;

    //GEM-mapper
//    string GemIndex = "./gem-indexer -i ";
//    GemIndex+=refPath;
//    GemIndex+=" -o ";
//    GemIndex+=refPre;

    string GemMapper = "./gem-mapper -I ";
    GemMapper+=refPre;
    GemMapper+=".gem";
    GemMapper+=" -i readSegment.fasta -o readSegment.sam";
//    system(GemIndex);
    system(GemMapper.c_str());
    cout << "Gem-mapper success" << endl;

//    //BWA
////    const char *BWAIndex = "./bwa -index 05_hg38.fa";
////    const char *BWAMapper = "./bwa mem -x pacbio 05_hg38.fa readSegment.fasta > readSegment.sam";
//    const char *BWAMapper1 = "./bwa aln 05_hg38.fa readSegment.fasta > readSegment.sai";
//    const char *BWAMapper2 = "./bwa samse 05_hg38.fa readSegment.sai readSegment.fasta > readSegment.sam";
////    bwa aln ref.fa reads.fq > reads.sai; bwa samse ref.fa reads.sai reads.fq > aln-se.sam
//
////    system(BWAIndex);
//    system(BWAMapper1);
//    system(BWAMapper2);
//    cout << "BWA success" << endl;


    //Separate SAM
    vector<ShortMapResult> SplitReadMapInfo = SeparateSamInfo();
    cout << "Separate sam success" << endl;

    //Combination
    for (const auto &SSI: SplitReadMapInfo) {
        SplitReadInfo[getReadName(SSI.ReadName)][getReadFragmentOrder(SSI.ReadName)].SplitReadMapResult.push_back(SSI);
    }
    cout << "Combination success" << endl;

    ofstream foutAnchor("anchor.txt");
    ofstream fout2("Read_AnchorCount.txt");
    unordered_map<int, int> countEveryRead;
    unordered_map<string, int> ReadAndAnchorCount;
    for(const auto item : SplitReadInfo) {
        int count=0;
        for(int i=0; i<item.second.size(); ++i) {
            count+=item.second[i].SplitReadMapResult.size();
        }

        if(count<20) {
            foutAnchor<<"OneReadless20-------------------"<<count<<endl;
            unordered_map<string, int> anchorInRef;
            for(int i=0; i<item.second.size(); ++i) {
                for(int j=0; j<item.second[i].SplitReadMapResult.size(); ++j) {
                    ++anchorInRef[item.second[i].SplitReadMapResult[j].RefName];
                }
            }
            for(const auto one : anchorInRef) {
                foutAnchor<<one.first<<"    "<<one.second<<endl;
            }
        }

        ++countEveryRead[count];
        ReadAndAnchorCount[item.first]=count;
        fout2<<item.first<<"    "<<count<<endl;
    }
    foutAnchor.close();

    int less10=0;
    int _10_20=0;
    int less20=0;
    int _20_40=0;
    int _40_60=0;
    int _60_80=0;
    int _80_100=0;
    int _100_120=0;
    int more120=0;

    int totalCount=0;
    ofstream fout("countEveryRead.txt");

    for(const auto item : countEveryRead) {
        totalCount+=item.second;
        if(item.first<20) {
            less20+=item.second;
            if(item.first<10) {
                less10+=item.second;
            } else {
                _10_20+=item.second;
            }
        } else if(item.first>=20 && item.first<40) {
            _20_40+=item.second;
        } else if(item.first>=40 && item.first<60) {
            _40_60+=item.second;
        } else if(item.first>=60 && item.first<80) {
            _60_80+=item.second;
        } else if(item.first>=80 && item.first<100) {
            _80_100+=item.second;
        } else if(item.first>=100 && item.first<120) {
            _100_120+=item.second;
        } else {
            more120+=item.second;
        }

        fout<<item.first<<"     "<<item.second<<endl;

    }
    fout<<endl;
    fout<<"totalCount   "<<totalCount<<endl;
    fout<<endl;
    fout<<"less10   "<<less10<<endl;
    fout<<"_10_20   "<<_10_20<<endl;
    fout<<endl;
    fout<<"less20   "<<less20<<endl;
    fout<<"_20_40   "<<_20_40<<endl;
    fout<<"_40_60   "<<_40_60<<endl;
    fout<<"_60_80   "<<_60_80<<endl;
    fout<<"_80_100   "<<_80_100<<endl;
    fout<<"_100_120   "<<_100_120<<endl;
    fout<<"more120   "<<more120<<endl;

    fout.close();
    cout<<"over"<<endl;
//    int stopNow;
//    cin>>stopNow;


    cout<<"collect info"<<endl;
    unordered_map<string, string> ReadInfo = ReadOrRefInfo("pb_human.fasta");
    unordered_map<string, string> RefInfo = ReadOrRefInfo("05_hg38.fa");
    cout<<"collect info over"<<endl;

    unordered_map<string, int> RefIndex;
    uint32_t refi=0;
    for(auto OneRef:RefInfo) {
        RefIndex[OneRef.first]=refi;
        ++refi;
    }

    OutputQueue outputQueue;
    initialize_output_queue(outPath.c_str(), RefInfo, &outputQueue, argc, argv);

    double ChainAndBaseLevelStartTime = clock();
    cout << "Chain and baselevel begin" << endl;

    bam1_t r;
    r.data=0;
    r.m_data=0;

    int diff1=0;
    int diff2=0;
    int diff=0;
    int DIFF=0;
    int diffSingle=0;

    vector<double> store;
    double maxZ=DBL_MIN;
    double minZ=DBL_MAX;
    ofstream zz0("zz0.txt");
    ofstream AnchOut("AnchOut.txt");
    ofstream smallFA("smallRead.fasta");

    for (auto OneRead: SplitReadInfo) {
        cout<<"----------OneRead---------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
        vector<BaseLevelPos> AlignPosArr;
        AlignPosArr=chainLDS(OneRead.second, SplitLength, 1500, 500);

        AnchOut<<AlignPosArr.size()<<endl;
        for(int i=0; i<AlignPosArr.size(); ++i) {
            AnchOut<<AlignPosArr[i].position.size()<<"  ";
        }
        AnchOut<<endl;

        vector<double> chainBest5(AlignPosArr.size());
        chainBest5[0]=100;
        double ave=100;
        for(int i=1; i<chainBest5.size(); ++i) {
            chainBest5[i]=(double)AlignPosArr[i].position.size() / (double)AlignPosArr[0].position.size() * 100;
            ave+=chainBest5[i];
        }
        ave/=5;
        double zz=0.0;
        for(int i=0; i<chainBest5.size(); ++i) {
            zz+=pow(chainBest5[i]-ave, 2);
        }
        for(int i=chainBest5.size(); i<5; ++i) {
            zz+=pow(ave, 2);
        }

        zz=sqrt(zz)/4;

        int choseFLAG=0;
        int maxAnchNum=0;
        if(zz==0) {
            for(int i=0; i<chainBest5.size(); ++i) {
                zz0<<AlignPosArr[i].position.size()<<" dir"<<AlignPosArr[i].dir<<
                "("<<AlignPosArr[i].AnchorNum<<"+"<<AlignPosArr[i].AnchorNumRev<<")"<<"  ";
                if((AlignPosArr[i].dir && AlignPosArr[i].AnchorNum>maxAnchNum) ||
                (!AlignPosArr[i].dir && AlignPosArr[i].AnchorNum>maxAnchNum)) {

                }
            }
            zz0<<endl;
            for(int i=0; i<chainBest5.size(); ++i) {
                zz0<<chainBest5[i]<<"   ";
            }
            zz0<<endl;
            zz0<<endl;
            smallFA<<">"<<OneRead.first<<endl;
            smallFA<<ReadInfo[OneRead.first]<<endl;
        }
        store.emplace_back(zz);

//        if(zz<61) {
//            AnchNumOut<<ReadAndAnchorCount[OneRead.first]<<endl;
//        }

        if(zz>maxZ) {
            maxZ=zz;
        }
        if(zz<minZ) {
            minZ=zz;
        }


//        for(int i=0; i<AlignPosArr.size(); ++i) {
            BaseLevelPos AlignPos = AlignPosArr[choseFLAG];
//        BaseLevelPos AlignPos = chain(OneRead.second, SplitLength, 500);
            int PosNumber=AlignPos.position.size();
//            if(zz!=0) {
//                PosNumber=0;
//            }
            cout<<"PosNumber:"<<PosNumber<<endl;
            long CCread=0;
            string totalCIGAR="";
            // forward
            if(AlignPos.dir && PosNumber>0) {
                // start
                // cout<<"1"<<endl;
                // cout<<"seq"<<endl;
                string ReadStart=ReadInfo[OneRead.first].substr(0, AlignPos.position[0].ReadPos);
                string RefStart=RefInfo[AlignPos.RefName].substr(max((long)0, AlignPos.position[0].RefPos-AlignPos.position[0].ReadPos), min(AlignPos.position[0].RefPos, (long)(AlignPos.position[0].ReadPos*1.5)));
                // cout<<ReadStart.size()<<" "<<RefStart.size()<<endl;
                CCread+=ReadStart.size();
                reverse(ReadStart.begin(), ReadStart.end());
                reverse(RefStart.begin(), RefStart.end());

                char* readS=(char*)ReadStart.c_str();
                char* refS=(char*)RefStart.c_str();
                // cout<<"edlib"<<endl;
                // cout<<strlen(readS)<<"  "<<strlen(refS)<<endl;
                EdlibAlignResult resHead= edlibBLHeadOrTail(readS, refS);

                // result.editDistance
                // result.alignmentLength
                // result.endLocations[0]
                // result.startLocations[0]
                // cout<<"cigar"<<endl;
                char* cigarHead = edlibAlignmentToCigar(resHead.alignment, resHead.alignmentLength, EDLIB_CIGAR_STANDARD);
                // cout<<"other"<<endl;
                string addMHead="";
                if(ReadStart.size()!= getIandM((string)cigarHead)) {
                    ++diffSingle;
                    if(ReadStart.size()>getIandM((string)cigarHead)) {
                        addMHead+= to_string(ReadStart.size()-getIandM((string)cigarHead));
                        addMHead+='M';
                    } else {
                        cout<<"ReadSeq.size()<getIandM((string)cigarHead)"<<endl;
                        int stop;
                        cin>>stop;
                    }
                }

                totalCIGAR+=addMHead;
                string cigarHeadRev=reverseCIGAR(cigarHead);
                totalCIGAR+=cigarHeadRev;

                AlignPos.RefStart=max((long)1, AlignPos.position[0].RefPos-resHead.endLocations[0]);

                free(cigarHead);
                edlibFreeAlignResult(resHead);
//            free(readS);
//            free(refS);

                // common
                // cout<<"2"<<endl;
                for(int i=0; i<PosNumber-1; ++i) {
                    // cout<<"seq"<<endl;
                    string ReadSeq = ReadInfo[OneRead.first].substr(min(AlignPos.position[i].ReadPos, AlignPos.position[i + 1].ReadPos), max(AlignPos.position[i].ReadPos, AlignPos.position[i + 1].ReadPos)-min(AlignPos.position[i].ReadPos, AlignPos.position[i + 1].ReadPos));
                    string RefSeq = RefInfo[AlignPos.RefName].substr(min(AlignPos.position[i].RefPos, AlignPos.position[i + 1].RefPos), max(AlignPos.position[i].RefPos, AlignPos.position[i + 1].RefPos)-min(AlignPos.position[i].RefPos, AlignPos.position[i + 1].RefPos));
                    CCread+=ReadSeq.size();

                    char* read=(char*)ReadSeq.c_str();
                    char* ref=(char*)RefSeq.c_str();
                    // cout<<"edlib"<<endl;
                    // cout<<strlen(read)<<"  "<<strlen(ref)<<endl;
                    EdlibAlignResult resCommon=edlibBLCommon(read, ref);

                    // result.editDistance
                    // result.alignmentLength
                    // result.endLocations[0]
                    // result.startLocations[0]
                    // cout<<"cigar"<<endl;
                    char* cigarCommon = edlibAlignmentToCigar(resCommon.alignment, resCommon.alignmentLength, EDLIB_CIGAR_STANDARD);
                    // cout<<"other"<<endl;
                    string addMCommon="";
                    if(ReadSeq.size()!= getIandM((string)cigarCommon)) {
                        ++diffSingle;
                        if(ReadSeq.size()>getIandM((string)cigarCommon)) {
                            addMCommon+= to_string(ReadSeq.size()-getIandM((string)cigarCommon));
                            addMCommon+='M';
                        } else {
                            cout<<"ReadSeq.size()<getIandM((string)cigarCommon)"<<endl;
                            int stop;
                            cin>>stop;
                        }
                    }
                    totalCIGAR+=(string)cigarCommon;
                    totalCIGAR+=addMCommon;
                    free(cigarCommon);
                    edlibFreeAlignResult(resCommon);
//                free(read);
//                free(ref);
                }

                // end
                // cout<<"3"<<endl;
                if(AlignPos.position[PosNumber-1].ReadPos<ReadInfo[OneRead.first].size()) {
                    // cout<<"seq"<<endl;
                    string ReadEnd=ReadInfo[OneRead.first].substr
                            (AlignPos.position[PosNumber-1].ReadPos,
                             ReadInfo[OneRead.first].size()-AlignPos.position[PosNumber-1].ReadPos);
                    string RefEnd=RefInfo[AlignPos.RefName].substr
                            (AlignPos.position[PosNumber-1].RefPos,
                             min((long)(RefInfo[AlignPos.RefName].size()-AlignPos.position[PosNumber-1].RefPos),
                                 (long)((ReadInfo[OneRead.first].size()-AlignPos.position[PosNumber-1].ReadPos)*1.5)));
                    CCread+=ReadEnd.size();

                    char* readE=(char*)ReadEnd.c_str();
                    char* refE=(char*)RefEnd.c_str();
                    // cout<<"edlib"<<endl;
                    // cout<<strlen(readE)<<"  "<<strlen(refE)<<endl;
                    EdlibAlignResult resTail=edlibBLHeadOrTail(readE, refE);

                    // result.editDistance
                    // result.alignmentLength
                    // result.endLocations[0]
                    // result.startLocations[0]
                    // cout<<"cigar"<<endl;
                    char* cigarTail = edlibAlignmentToCigar(resTail.alignment, resTail.alignmentLength, EDLIB_CIGAR_STANDARD);
                    // cout<<"other"<<endl;
                    string addMTail="";
                    if(ReadEnd.size()!= getIandM((string)cigarTail)) {
                        ++diffSingle;
                        if(ReadEnd.size()>getIandM((string)cigarTail)) {
                            addMTail+= to_string(ReadEnd.size()-getIandM((string)cigarTail));
                            addMTail+='M';
                        } else {
                            cout<<"ReadEnd.size()<getIandM((string)cigarTail)"<<endl;
                            int stop;
                            cin>>stop;
                        }
                    }

                    totalCIGAR+=(string)cigarTail;
                    totalCIGAR+=addMTail;

                    free(cigarTail);
                    edlibFreeAlignResult(resTail);
//                free(readE);
//                free(refE);
                }
            }
                // reverse
            else if(PosNumber>0) {
                // start
                // cout<<"4"<<endl;
                if(AlignPos.position[0].ReadPos<ReadInfo[OneRead.first].size()) {
                    // cout<<"seq"<<endl;
                    string ReadEnd=ReadInfo[OneRead.first].substr
                            (AlignPos.position[0].ReadPos,
                             ReadInfo[OneRead.first].size()-AlignPos.position[0].ReadPos);
                    string RefEnd=RefInfo[AlignPos.RefName].substr
                            (max((long)0, (long)(AlignPos.position[0].RefPos-(ReadInfo[OneRead.first].size()-AlignPos.position[0].ReadPos)*1.5+SplitLength)),
                             min((long)(AlignPos.position[0].RefPos+SplitLength), (long)((ReadInfo[OneRead.first].size()-AlignPos.position[0].ReadPos)*1.5)));
//                min(RefInfo[AlignPos.RefName].size()-AlignPos.position[0].RefPos, ReadInfo[OneRead.first].size()-AlignPos.position[0].ReadPos)
                    CCread+=ReadEnd.size();
                    ReadEnd=CharToComple(ReadEnd);
//                reverse(ReadEnd.begin(), ReadEnd.end());
//                reverse(ReadEnd.begin(), ReadEnd.end());
                    reverse(RefEnd.begin(), RefEnd.end());
                    char* readE=(char*)ReadEnd.c_str();
                    char* refE=(char*)RefEnd.c_str();
                    // cout<<"edlib"<<endl;
                    // cout<<strlen(readE)<<"  "<<strlen(refE)<<endl;
                    EdlibAlignResult resHead= edlibBLHeadOrTail(readE, refE);

                    // result.editDistance
                    // result.alignmentLength
                    // result.endLocations[0]
                    // result.startLocations[0]
                    // cout<<"cigar"<<endl;
                    char* cigarHead = edlibAlignmentToCigar(resHead.alignment, resHead.alignmentLength, EDLIB_CIGAR_STANDARD);
                    // cout<<"other"<<endl;
                    string addMHead="";
                    if(ReadEnd.size()!= getIandM((string)cigarHead)) {
                        ++diffSingle;
                        if(ReadEnd.size()>getIandM((string)cigarHead)) {
                            addMHead+= to_string(ReadEnd.size()-getIandM((string)cigarHead));
                            addMHead+='M';
                        } else {
                            // cout<<"ReadSeq.size()<getIandM((string)cigarHead)"<<endl;
                            int stop;
                            cin>>stop;
                        }
                    }

                    totalCIGAR+=addMHead;

                    string cigarHeadRev=reverseCIGAR(cigarHead);
                    totalCIGAR+=cigarHeadRev;

                    AlignPos.RefStart=max((long)1, (long)(AlignPos.position[0].RefPos-resHead.endLocations[0]));

                    free(cigarHead);
                    edlibFreeAlignResult(resHead);
//                free(readE);
//                free(refE);
                }
                // common
                // cout<<"5"<<endl;
                for(int i=0; i<PosNumber-1; ++i) {
                    // cout<<"seq"<<endl;
                    string ReadSeq = ReadInfo[OneRead.first].substr(min(AlignPos.position[i].ReadPos, AlignPos.position[i + 1].ReadPos), max(AlignPos.position[i].ReadPos, AlignPos.position[i + 1].ReadPos)-min(AlignPos.position[i].ReadPos, AlignPos.position[i + 1].ReadPos));
                    string RefSeq = RefInfo[AlignPos.RefName].substr(min(AlignPos.position[i].RefPos, AlignPos.position[i + 1].RefPos)+SplitLength, max(AlignPos.position[i].RefPos, AlignPos.position[i + 1].RefPos)-min(AlignPos.position[i].RefPos, AlignPos.position[i + 1].RefPos));
                    ReadSeq=CharToComple(ReadSeq);
                    reverse(ReadSeq.begin(), ReadSeq.end());
                    CCread+=ReadSeq.size();

                    char* read=(char*)ReadSeq.c_str();
                    char* ref=(char*)RefSeq.c_str();
                    // cout<<"edlib"<<endl;
                    // cout<<strlen(read)<<"  "<<strlen(ref)<<endl;
                    EdlibAlignResult resCommon=edlibBLCommon(read, ref);

                    // result.editDistance
                    // result.alignmentLength
                    // result.endLocations[0]
                    // result.startLocations[0]
                    // cout<<"cigar"<<endl;
                    char* cigarCommon = edlibAlignmentToCigar(resCommon.alignment, resCommon.alignmentLength, EDLIB_CIGAR_STANDARD);
                    // cout<<"other"<<endl;
                    string addMCommon="";
                    if(ReadSeq.size()!= getIandM((string)cigarCommon)) {
                        ++diffSingle;
                        if(ReadSeq.size()>getIandM((string)cigarCommon)) {
                            addMCommon+= to_string(ReadSeq.size()-getIandM((string)cigarCommon));
                            addMCommon+='M';
                        } else {
                            cout<<"ReadSeq.size()<getIandM((string)cigarCommon)"<<endl;
                            int stop;
                            cin>>stop;
                        }
                    }
                    totalCIGAR+=(string)cigarCommon;
                    totalCIGAR+=addMCommon;

                    free(cigarCommon);
                    edlibFreeAlignResult(resCommon);
//                free(read);
//                free(ref);
                }

                // end
                // cout<<"6"<<endl;
                // cout<<"seq"<<endl;
                string ReadStart=ReadInfo[OneRead.first].substr(0, AlignPos.position[PosNumber-1].ReadPos);
                string RefStart=RefInfo[AlignPos.RefName].substr
                        (min((long)(AlignPos.position[PosNumber-1].RefPos+SplitLength), (long)RefInfo[AlignPos.RefName].size()),
                         min((long)(RefInfo[AlignPos.RefName].size()-AlignPos.position[PosNumber-1].RefPos-SplitLength), (long)(AlignPos.position[PosNumber-1].ReadPos*1.5)));
                CCread+=ReadStart.size();
                ReadStart=CharToComple(ReadStart);
                reverse(ReadStart.begin(), ReadStart.end());

                char* readS=(char*)ReadStart.c_str();
                char* refS=(char*)RefStart.c_str();
                // cout<<"edlib"<<endl;
                // cout<<strlen(readS)<<"  "<<strlen(refS)<<endl;
                EdlibAlignResult resTail=edlibBLHeadOrTail(readS, refS);

                // result.editDistance
                // result.alignmentLength
                // result.endLocations[0]
                // result.startLocations[0]
                // cout<<"cigar"<<endl;
                char* cigarTail = edlibAlignmentToCigar(resTail.alignment, resTail.alignmentLength, EDLIB_CIGAR_STANDARD);
                // cout<<"other"<<endl;
                string addMTail="";
                if(ReadStart.size()!= getIandM((string)cigarTail)) {
                    ++diffSingle;
                    if(ReadStart.size()>getIandM((string)cigarTail)) {
                        addMTail+= to_string(ReadStart.size()-getIandM((string)cigarTail));
                        addMTail+='M';
                    } else {
                        cout<<"ReadEnd.size()<getIandM((string)cigarTail)"<<endl;
                        int stop;
                        cin>>stop;
                    }
                }

                totalCIGAR+=(string)cigarTail;
                totalCIGAR+=addMTail;

                free(cigarTail);
                edlibFreeAlignResult(resTail);
//            free(readS);
//            free(refS);
            }

            uint16_t flag=4;
            if(PosNumber<1) {
                flag=4;
            } else if(!AlignPos.dir) {
                flag=16;
            } else {
                flag=1;
            }
            uint8_t MapQ=255;

            if(CCread!=ReadInfo[OneRead.first].size() && CCread!=0) {
                ++diff;
                if(flag==1) {
                    ++diff1;
                } else {
                    ++diff2;
                }
            }
            vector<uint32_t> CIGAR= handleCIGAR(totalCIGAR, ReadInfo[OneRead.first].size(), DIFF);

            generate_BAM_t(OneRead.first, flag, RefIndex[AlignPos.RefName], AlignPos.RefStart, MapQ,\
CIGAR, '*', 0, 0, ReadInfo[OneRead.first], "", "",&r);

            outputQueue.samOutputVec.emplace_back(r);

            if(outputQueue.samOutputVec.size()>=500) {
                output_Queue(&outputQueue);
            }
//        }

    }


    sort(store.begin(), store.end());
    ofstream storeOut("storeOut.txt");
    for(int i=0; i<store.size(); ++i) {
        storeOut<<store[i]<<endl;
    }

    // cout<<"store start"<<endl;
    // cout<<maxZ<<"   "<<minZ<<endl;
    vector<int> storeRes(11, 0);
    for(int i=0; i<store.size(); ++i) {
        ++storeRes[(int)(store[i]/(maxZ-minZ)*10)];
    }
    for(int i=0; i<storeRes.size(); ++i) {
        cout<<storeRes[i]<<endl;
    }
    // cout<<"store end"<<endl;

//    BaseLevel.close();
    output_Queue(&outputQueue);
    // cout<<diff<<endl;
    // cout<<diff1<<endl;
    // cout<<diff2<<endl;
    // cout<<DIFF<<endl;
    // cout<<diffSingle<<endl;

    double ChainAndBaseLevelEndTime = clock();
    cout << "The run time of chain and baselevel align is: "
         << (double) (ChainAndBaseLevelEndTime - ChainAndBaseLevelStartTime) / CLOCKS_PER_SEC << "s" << std::endl;
    cout << "Chain and baselevel align success" << endl;

    return 0;
}
