#!/usr/bin/env python3
import pysam
import sys

filename_bam = sys.argv[1]
filename_out = filename_bam.replace('.bam', '') + '.nVR_freq'

umi_len = 6
VR_len = 36

read2VR = dict()
VR_freq = dict()

save = pysam.set_verbosity(0)
bamfile = pysam.AlignmentFile(filename_bam, "rb")
pysam.set_verbosity(save)

for s in bamfile:
    q_name = s.query_name
    t_name = s.reference_name
    q_seq = s.query_sequence
    c_tuples = s.cigartuples

    if c_tuples == None:
        continue
    
    if q_name not in read2VR:
        read2VR[q_name] = dict()
        read2VR[q_name]['R1'] = {'umi': '', 'VR_seq': '', 'is_multi':False}
        read2VR[q_name]['R2'] = {'umi': '', 'VR_seq': '', 'is_multi':False}
    
    if s.is_reverse == True and s.is_read1 == True and t_name.endswith('_flank3'):
        umi_3 = ''
        (tmp_umi_code, tmp_umi_len) = c_tuples[-1]
        if tmp_umi_code == 4 and tmp_umi_len == umi_len:
            umi_3 = q_seq[umi_len*-1:]
        
        if umi_3 == '':
            continue
       
        tmp_VR_seq = ''
        (tmp_VR_code, tmp_VR_len) = c_tuples[0]
        if tmp_VR_code == 4:
            tmp_VR_start = tmp_VR_len - VR_len
            tmp_VR_seq = q_seq[tmp_VR_start:tmp_VR_len]
        
        if read2VR[q_name]['R1']['umi'] != '':
            read2VR[q_name]['R1']['is_multi'] = True
            continue
        
        read2VR[q_name]['R1']['umi'] = umi_3
        read2VR[q_name]['R1']['VR_seq'] = tmp_VR_seq

        if tmp_VR_seq not in VR_freq:
            VR_freq[tmp_VR_seq] = {'total': 0, 'match': 0, 'R1_missed': 0, 'R2_missed':0, 'mismatch':0, 'umi':[]}
        VR_freq[tmp_VR_seq]['total'] += 1

    elif s.is_reverse == False and s.is_read2 == True and t_name.endswith('_flank5'):
        umi_5 = ''
        (tmp_umi_code, tmp_umi_len) = c_tuples[0]
        if tmp_umi_code == 4 and tmp_umi_len == umi_len:
            umi_5 = q_seq[:umi_len]
        
        # if UMI is not found, skip the read.
        if umi_5 == '':
            continue

        tmp_VR_start = umi_len
        for (tmp_code, tmp_len) in c_tuples:
            if tmp_code == 0:   # match
                tmp_VR_start += tmp_len
            elif tmp_code == 1: # insert
                tmp_VR_start += tmp_len
            elif tmp_code == 2: # delete
                pass

        tmp_flank_5 = q_seq[umi_len:tmp_VR_start]
        tmp_VR_seq = q_seq[tmp_VR_start:tmp_VR_start+VR_len]
        flank_3 = q_seq[tmp_VR_start+VR_len:]

        if read2VR[q_name]['R2']['umi'] != '':
            read2VR[q_name]['R2']['is_multi'] = True
            continue

        read2VR[q_name]['R2']['umi'] = umi_5
        read2VR[q_name]['R2']['VR_seq'] = tmp_VR_seq

        if tmp_VR_seq not in VR_freq:
            VR_freq[tmp_VR_seq] = {'total': 0, 'match': 0, 'R1_missed': 0, 'R2_missed':0, 'mismatch':0, 'umi':[]}
        VR_freq[tmp_VR_seq]['total'] += 1


VR_counts = {'match': 0, 'mismatch': 0, 'R1_missed': 0, 'R2_missed': 0, 'empty': 0}

for q_name, tmp in read2VR.items():
    VR_R1 = tmp['R1']['VR_seq']
    VR_R2 = tmp['R2']['VR_seq']
    VR_umi = '%s.%s' % (tmp['R1']['umi'], tmp['R2']['umi'])

    if VR_R1 != '':
        if VR_R1 == VR_R2:
            VR_freq[VR_R1]['match'] += 1
            VR_freq[VR_R1]['umi'].append(VR_umi)
            VR_counts['match'] += 1

        elif VR_R2 == '':
            VR_freq[VR_R1]['R1_missed'] += 1
            VR_counts['R1_missed'] += 1

        else:
            VR_freq[VR_R1]['mismatch'] += 1
            VR_freq[VR_R2]['mismatch'] += 1
            VR_counts['mismatch'] += 1

    else:   # VR_R1 == ''
        if VR_R2 != '':
            VR_freq[VR_R2]['R2_missed'] += 1
            VR_counts['R2_missed'] += 1
        else:
            VR_counts['empty'] += 1

sys.stderr.write('Write %s\n' % filename_out)
f_out = open(filename_out, 'w')
f_out.write('#VRseq\tTotal\tR1_missed\tR2_missed\tMismatch\tMatched\tUniqMatched\n')
for tmp_VR, tmp in VR_freq.items():
    tmp_umi_count = len(list(set(tmp['umi'])))
    if tmp_VR == '':
        continue

    f_out.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\n' % 
                (tmp_VR, tmp['total'], tmp['R1_missed'], tmp['R2_missed'], 
                 tmp['mismatch'], tmp['match'], tmp_umi_count))
f_out.close()
