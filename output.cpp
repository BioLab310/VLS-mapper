//
// Created by spoony on 23-3-11.
//
#include "output.h"

//void initialize_output_queue(const char *output_file_path, \
//const SequenceBatch *reference_sequence_batch, \
//int num_mapping_threads, size_t max_queue_size, OutputQueue *output_queue) {
void initialize_output_queue(const char *output_file_path, unordered_map<string, string>& RefInfo, OutputQueue *output_queue, int argc, char* argv[]) {
    output_queue->output_sam_file = sam_open_format(output_file_path, "w", NULL);
    output_queue->sam_header = sam_hdr_init();
    output_sam_header(output_file_path, RefInfo, output_queue, argc, argv);
    fprintf(stderr, "Initialize output queue successfully!\n");
}

void output_sam_header(const char *output_file_path, unordered_map<string, string> RefInfo, OutputQueue *output_queue, int argc, char* argv[]) {
    //sam_hdr_t sam_header;
    output_queue->sam_header->n_targets = RefInfo.size();
    output_queue->sam_header->target_len = (uint32_t*) malloc(output_queue->sam_header->n_targets * sizeof(uint32_t));
    output_queue->sam_header->target_name = (char**) malloc(output_queue->sam_header->n_targets * sizeof(char*));
    kstring_t header_kstr = {0, 0, NULL};
    // TODO(Haowen): add PG later
    ksprintf(&header_kstr, "@PG\tID:mapper\tPN:mapper\tVN:1.0\tCL:");
    for (int i = 0; i < argc; ++i) {
      ksprintf(&header_kstr, " %s", argv[i]);
    }
    ksprintf(&header_kstr, "\n");
    uint32_t ri=0;
    for(auto OneRef:RefInfo) {
        output_queue->sam_header->target_len[ri] = RefInfo[OneRef.first].size();
        output_queue->sam_header->target_name[ri] =(char*)malloc(sizeof(char)*(OneRef.first.size()+1));
        output_queue->sam_header->target_name[ri][OneRef.first.size()]='\0';
        for(uint32_t ii=0;ii<OneRef.first.size();ii++)
        {
            output_queue->sam_header->target_name[ri][ii]=OneRef.first[ii];
        }

//        = (char *)OneRef.first.c_str();
        ksprintf(&header_kstr, "@SQ\tSN:%s\tLN:%d\n", output_queue->sam_header->target_name[ri], output_queue->sam_header->target_len[ri]);
        ++ri;
    }

    output_queue->sam_header->l_text = ks_len(&header_kstr);
    output_queue->sam_header->text = ks_str(&header_kstr);
    output_queue->sam_header->sdict = NULL;
    output_queue->sam_header->hrecs = NULL;
    output_queue->sam_header->ref_count = 1;
    int htslib_err = sam_hdr_write(output_queue->output_sam_file, output_queue->sam_header);
    assert(htslib_err == 0);
}

//void init_BAM_t(bam1_t* r)
//{
//    r->
//}

/*! @typedef
 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  id
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  mempolicy  memory handling policy, see bam_set_mempolicy()
 */

/*! @typedef
 @abstract Structure for core alignment information.
 @field  pos     0-based leftmost coordinate
 @field  tid     chromosome ID, defined by sam_hdr_t
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_extranul length of extra NULs between qname & cigar (for alignment)
 @field  flag    bitwise flag
 @field  l_qname length of the query name
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by sam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template
 @field  isize   observed template length ("insert size")
 */

    void generate_BAM_t(string readName, uint16_t flag, int refId, long refStart, uint8_t MapQ,\
vector<uint32_t>& cigar, char a, int b, int c, string readSeq, string readQual, string options, bam1_t* r) {
    r->core.pos=refStart;
    r->core.tid=refId;
    r->core.qual=MapQ;
    r->core.flag=flag;

    r->core.l_extranul=0;
    if((readName.size()+1)%4!=0) {
        r->core.l_extranul=4-((readName.size()+1)%4);
    }
    r->core.l_qname=readName.size()+1+r->core.l_extranul;

    r->core.n_cigar=cigar.size();
    r->core.l_qseq=readSeq.size();
    r->core.mtid=-1;
    r->core.mpos=-1;
    r->core.isize=0;

//    r->l_data=r->core.l_qname+(r->core.n_cigar<<2)+((r->core.l_qseq+1)>>1)+1;
    r->l_data=r->core.l_qname+(r->core.n_cigar<<2)+((r->core.l_qseq+1)>>1)+r->core.l_qseq;
    // r->l_data=r->core.l_qname+(r->core.n_cigar<<2)+((r->core.l_qseq+1)>>1);
    r->data=(uint8_t*)malloc(r->l_data*sizeof(uint8_t));
    r->m_data=r->l_data;

    // readName
    char * p_tmp=0;
    p_tmp=bam_get_qname(r);
    for(uint32_t ii=0;ii<readName.size();ii++)
    {
        p_tmp[ii]=readName[ii];
    }

    bam_get_qname(r)[readName.size()]='\0';

    //CIGAR
    memcpy(bam_get_cigar(r), &(cigar[0]), r->core.n_cigar*sizeof(uint32_t));

    //readSeq
    uint8_t *seq=bam_get_seq(r);
    for(size_t i = 0; i < readSeq.size(); ++i) {
        bam_set_seqi(seq, i, seq_nt16_table[(uint8_t)readSeq[i]]);
    }

    //readQual
    uint8_t * p_tmp1=bam_get_qual(r);
//    p_tmp1[0]='*';
    for(uint32_t ii=0;ii<readSeq.size();ii++)
    {
//        p_tmp1[ii]='\0';
        p_tmp1[ii]=readSeq[ii];
    }

    for(int i=0; i<readSeq.size(); ++i) {
        p_tmp1[i]-=33;
    }
//    NM:i:2322	ms:i:18279	AS:i:18198	nn:i:0	tp:A:P	cm:i:342	s1:i:3187	s2:i:0	de:f:0.1405	rl:i:0
//    bam_aux_update_str(sam_alignment, "MD", ks_len(MD_tag) + 1, ks_str(MD_tag));
//    bam_aux_update_int(r, "NM", 2322);
//    bam_aux_update_int(r, "ms", 100);
//    bam_aux_update_int(r, "AS", 13198);
//    bam_aux_update_int(r, "nn", 100);
//    bam_aux_update_int(r, "cm", 100);
//    bam_aux_update_int(r, "s1", 100);
//    bam_aux_update_int(r, "s2", 100);
//    bam_aux_update_int(r, "rl", 100);
}

void output_Queue(OutputQueue *a)
{
    if(a->samOutputVec.size()>0) {
        for (size_t si = 0; si < a->samOutputVec.size(); ++si) {
            int htslib_err = sam_write1(a->output_sam_file, a->sam_header, &(a->samOutputVec[si]));
            assert(htslib_err >= 0);
        }
        a->samOutputVec.clear();
    }
}