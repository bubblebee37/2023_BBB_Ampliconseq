#!/usr/bin/python
import sys

filename_VR_freq = sys.argv[1]
filename_out = filename_VR_freq.replace('.nVR_freq', '') + '.pVR_freq'

trans_tbl = dict()

## http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
AAs   = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
Base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
Base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
Base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
for i in range(0,len(AAs)):
    trans_tbl['%s%s%s'%(Base1[i],Base2[i],Base3[i])] = AAs[i]


def translate(tmp_nseq):
    rv = []
    for i in range(0,len(tmp_nseq)-2,3):
        tmp_codon = tmp_nseq[i:i+3]
        if( not trans_tbl.has_key(tmp_codon) ):
            rv.append('*')
        else:
            rv.append( trans_tbl[tmp_codon] )
    return ''.join(rv)

pseq_freq = dict()
pseq2nseq = dict()
f_VR = open(filename_VR_freq, 'r')
f_VR.readline()
for line in f_VR:
    tokens = line.strip().split("\t")
    count_matched = int(tokens[5])
    count_uniq = int(tokens[6])
    tmp_nseq = tokens[0]
    tmp_pseq = translate(tmp_nseq)

    if count_uniq == 0:
        continue
    
    if tmp_pseq not in pseq_freq:
        pseq_freq[tmp_pseq] = 0
        pseq2nseq[tmp_pseq] = dict()

    pseq_freq[tmp_pseq] += count_uniq
    pseq2nseq[tmp_pseq][tmp_nseq] = count_uniq
f_VR.close()

pseq_list = sorted(pseq_freq.keys(), key=pseq_freq.get, reverse=True)
pseq_total = sum(pseq_freq.values())
pseq_count_list = sorted([x for x in pseq_freq.values()], reverse=True)
pseq_over3 = sum([x for x in pseq_freq.values() if x >= 3])

f_out = open(filename_out, 'w')
f_out.write('#pseq\tRank\tCount\tpct_total\tpct_over3\tnseq_count\tnseq_rep\tnseq_rep_count\n')
for tmp_pseq in pseq_list:
    tmp_nseq_count = len(pseq2nseq[tmp_pseq])
    tmp_nseq_rep = sorted(pseq2nseq[tmp_pseq].keys(), key=pseq2nseq[tmp_pseq].get)[-1]

    if pseq_freq[tmp_pseq] < 3:
        continue

    tmp_rank = pseq_count_list.index(pseq_freq[tmp_pseq]) + 1
    f_out.write("%s\t%d\t%d\t%.2f\t%.2f\t%d\t%s\t%d\n" %
                       (tmp_pseq, tmp_rank, pseq_freq[tmp_pseq], 
                        pseq_freq[tmp_pseq]*100.0/pseq_total, 
                        pseq_freq[tmp_pseq]*100.0/pseq_over3, 
                        tmp_nseq_count, tmp_nseq_rep, 
                        pseq2nseq[tmp_pseq][tmp_nseq_rep]) )
f_out.close()
