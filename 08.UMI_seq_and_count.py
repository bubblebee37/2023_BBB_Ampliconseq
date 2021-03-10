#UMI_list 담긴 파일 만들고, pie chart 
import pandas as pd
import sys
#아래 적은 format 대로 파일명을 미리 저장해둘 것 / 파일 경로 입력 : ~ 채우기
c1 = 'chip1'
c2 = 'chip2'
l_1 = 'input'
#새 lot number의 phage로 만든 library sample
l_2 = 'input-new'
t = 'well'
m = 'invivo'
amp_1 = 1
amp_2 = 2
amp_3 = 3
folder = '/home/kyungha/Downloads/BBB_phage_display/Data/~'

C1_1 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, c1, amp_1), sep='\t')
C1_2 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, c1, amp_2), sep='\t')
C1_3 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, c1, amp_3), sep='\t')
C2_1 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, c2, amp_1), sep='\t')
C2_2 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, c2, amp_2), sep='\t')
C2_3 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, c2, amp_3), sep='\t')
T_1 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, t, amp_1), sep='\t')
T_2 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, t, amp_2), sep='\t')
T_3 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, t, amp_3), sep='\t')
L1_1 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, l_1, amp_1), sep='\t')
L1_2 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, l_1, amp_2), sep='\t')
L1_3 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, l_1, amp_3), sep='\t')
L2_1 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, l_2, amp_1), sep='\t')
L2_2 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, l_2, amp_2), sep='\t')
L2_3 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, l_2, amp_3), sep='\t')
M_1 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, m, amp_1), sep='\t')
M_2 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, m, amp_2), sep='\t')
M_3 = pd.read_csv('{}/{}-{}-.PhageM13.nVR_freq'.format(folder, m, amp_3), sep='\t')

#UMI_frequency를 정리할 파일의 list를 뽑기 (~에 지정한 csv 파일 넣기)
a = ~['#VRseq'].tolist()
b = ~['UniqMatched'].tolist()
c = {key:value for key,value in zip(a, b)}

def translate_rna(s):
    codon = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N", 
             "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T", 
             "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S", 
             "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I", 
             "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H", 
             "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P", 
             "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R", 
             "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L", 
             "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D", 
             "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A", 
             "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G", 
             "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V", 
             "UAA":"_", "UAC":"Y", "UAG":"_", "UAU":"T", 
             "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S", 
             "UGA":"_", "UGC":"C", "UGG":"W", "UGU":"C", 
             "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}

    l = [codon.get(s[n:n+3], 'X') for n in range(0, len(s), 3)]
    return "".join(l)

#같은 protein으로 나왔지만 다른 염기 서열을 갖는 것 합쳐서 frequency 다시 계산
prot_f = {}
for key, val in c.items():
    re_key = key.replace('T','U')
    #print(re_key)
    #print(translate_rna(re_key),val)
    prot=translate_rna(re_key)
    if prot in prot_f:
        prot_f[prot] += val
    else:
        prot_f[prot] = val
        
#~에 원하는 sample명과 amplification number 넣기 (예 : input_new 등)
fw = open('~_{}.txt'.format(~), 'w')
fw.write('#Protein_sequence\tUMI\n')
dic_overlap = {k: v for k, v in sorted(prot_f.items(), key=lambda item: item[1])}
for k in list(dic_overlap.keys())[::-1]:
    v = dic_overlap[k]
    print(k, v)
    fw.write(k+'\t'+str(v)+'\n')
    
#List형에서 frequency가 높은 순에서 낮은 순으로 정렬
def f2(x):
    return x[1]
order=sorted(c.items(), key=(lambda x:x[1]), reverse=True)


#Dictionary형으로 다시 바꿔줌
order_1=dict(order)
for key, val in order_1.items():
    print(key)
    print(val)

#Frequency가 0인 VR sequence를 전부 제외함 --> pie chart에서 제외하기 위해서!
order_d={}
for key, val in order_1.items():
    if val != 0:
        order_d[key] = val        
#UMI List 뽑아서 frequency 갯수 세서 차지하는 비율 확인
seq_num_for_UMI_freq = order_d.values()
list_UMI_freq = list(seq_num_for_UMI_freq)

count={}
for i in list_UMI_freq:
    try: count[i] += 1
    except: count[i]=1
print(count)
umi_in_seq_val = order_d.values()
umi_count = list(umi_in_seq_val)
print(umi_count)

#Pie chart의 mean value 구하기
import numpy as np
average = np.array(umi_count)
np.mean(average)
key = count.keys()
count_k_label = list(key)
value = count.values()
count_v_ratio = list(value)
print(count_k_label)
print(count_v_ratio)

#Pie chart로 나타내기
import matplotlib.pyplot as plt
#colors에 추가로 색 코드 더 넣어서 label 갯수에 맞추기
colors = ['#1f77b4', '#d62728', '#8c564b', '#e377c2','#9467bd', '#7f7f7f']
wedgeprops={'width': 0.7, 'edgecolor': 'w', 'linewidth': 5}
plt.pie(c1_ratio, labels=c1_label, autopct='%.1f%%', startangle=260, counterclock=False, colors=colors, wedgeprops=wedgeprops, pctdistance=0.65, labeldistance=1.05)
plt.title('Ratio of UMI frequency in each sequence', fontsize=40)
plt.rcParams["figure.figsize"] = (20, 20)
plt.text(1,1,'Mean: {:.2f}'.format(average.mean()), fontsize=30)
plt.rcParams['font.size']=15
plt.show()
