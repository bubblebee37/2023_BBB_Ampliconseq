#Normalized ratio와 rank가 포함된 전체 excel 파일 만들기
#파일명을 아래대로 미리 지정해둘 것

import pandas as pd
import sys
c1 = 'chip1'
c2 = 'chip2'
l_1 = 'input'
l_2 = 'input-new'
t = 'well'
m = 'invivo'
amp_1 = 1
amp_2 = 2
amp_3 = 3
folder = '/home/kyungha/Downloads/BBB_phage_display/Data/210217_Sequencing/Fastq/Trimmed'
C1_1 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, c1, amp_1), sep='\t')
C1_2 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, c1, amp_2), sep='\t')
C1_3 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, c1, amp_3), sep='\t')
C2_1 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, c2, amp_1), sep='\t')
C2_2 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, c2, amp_2), sep='\t')
C2_3 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, c2, amp_3), sep='\t')
T_1 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, t, amp_1), sep='\t')
T_2 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, t, amp_2), sep='\t')
T_3 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, t, amp_3), sep='\t')
L1_1 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, l_1, amp_1), sep='\t')
L1_2 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, l_1, amp_2), sep='\t')
L1_3 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, l_1, amp_3), sep='\t')
L2_1 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, l_2, amp_1), sep='\t')
L2_2 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, l_2, amp_2), sep='\t')
L2_3 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, l_2, amp_3), sep='\t')
M_1 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, m, amp_1), sep='\t')
M_2 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, m, amp_2), sep='\t')
M_3 = pd.read_csv('{}/{}-{}-.PhageM13.pVR_freq'.format(folder, m, amp_3), sep='\t')

C1_1_r = C1_1.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
C1_2_r = C1_2.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
C1_3_r = C1_3.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
C2_1_r = C2_1.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
C2_2_r = C2_2.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
C2_3_r = C2_3.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
T_1_r = T_1.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
T_2_r = T_2.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
T_3_r = T_3.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
L1_1_r = L1_1.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
L1_2_r = L1_2.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
L1_3_r = L1_3.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
L2_1_r = L2_1.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
L2_2_r = L2_2.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
L2_3_r = L2_3.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
M_1_r = M_1.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
M_2_r = M_2.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)
M_3_r = M_3.drop(['Rank', 'Count', 'nseq_count', 'pct_over3', 'nseq_rep', 'nseq_rep_count'], axis = 1)

m_1 = pd.merge(L1_1_r, L1_2_r, on = '#pseq', how = 'outer')
m_2 = pd.merge(m_1, L1_3_r, on = '#pseq', how = 'outer')
m_3 = pd.merge(m_2, L2_1_r, on = '#pseq', how = 'outer')
m_4 = pd.merge(m_3, L2_2_r, on = '#pseq', how = 'outer')
m_5 = pd.merge(m_4, L2_3_r, on = '#pseq', how = 'outer')
m_6 = pd.merge(m_5, T_1_r, on = '#pseq', how = 'outer')
m_7 = pd.merge(m_6, T_2_r, on = '#pseq', how = 'outer')
m_8 = pd.merge(m_7, T_3_r, on = '#pseq', how = 'outer')
m_9 = pd.merge(m_8, C1_1_r, on = '#pseq', how = 'outer')
m_10 = pd.merge(m_9, C1_2_r, on = '#pseq', how = 'outer')
m_11 = pd.merge(m_10, C1_3_r, on = '#pseq', how = 'outer')
m_12 = pd.merge(m_11, C2_1_r, on = '#pseq', how = 'outer')
m_13 = pd.merge(m_12, C2_2_r, on = '#pseq', how = 'outer')
m_14 = pd.merge(m_13, C2_3_r, on = '#pseq', how = 'outer')
m_15 = pd.merge(m_14, M_1_r, on = '#pseq', how = 'outer')
m_16 = pd.merge(m_15, M_2_r, on = '#pseq', how = 'outer')
m_17 = pd.merge(m_16, M_3_r, on = '#pseq', how = 'outer')

m_17.columns = ["Protein_sequence", "Input_1", "Input_2", "Input_3", "Input_new_1", "Input_new_2", "Input_new_3", "Transwell_1", "Transwell_2", "Transwell_3", "Chip_1um_1", "Chip_1um_2", "Chip_1um_3", "Chip_2um_1", "Chip_2um_2", "Chip_2um_3", "In_vivo_1", "In_vivo_2", "In_vivo_3"]

#전체 protein sequence를 포함한 dataframe 지정 --> m_pseq_r_total
m_17 = m_17.fillna(0)
m_pseq_r_total = m_14.sort_values(by=["In_vivo_3"], axis=0, ascending=False)
