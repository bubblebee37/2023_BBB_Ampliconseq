#log로 나타내기 위해 value 0 값을 다 0.0001로 변경
m_pseq_scatter = m_pseq_r_total.replace(0, 0.0001)

i1_s = m_pseq_scatter["Input_1"].tolist()
i2_s = m_pseq_scatter["Input_2"].tolist()
i3_s = m_pseq_scatter["Input_3"].tolist()
i_n_1_s = m_pseq_scatter["Input_new_1"].tolist()
i_n_2_s = m_pseq_scatter["Input_neW_2"].tolist()
i_n_3_s = m_pseq_scatter["Input_new_3"].tolist()
t1_s = m_pseq_scatter["Transwell_1"].tolist()
t2_s = m_pseq_scatter["Transwell_2"].tolist()
t3_s = m_pseq_scatter["Transwell_3"].tolist()
c1_1_s = m_pseq_scatter["Chip_1um_1"].tolist()
c1_2_s = m_pseq_scatter["Chip_1um_2"].tolist()
c1_3_s = m_pseq_scatter["Chip_1um_3"].tolist()
c2_1_s = m_pseq_scatter["Chip_2um_1"].tolist()
c2_2_s = m_pseq_scatter["Chip_2um_2"].tolist()
c2_3_s = m_pseq_scatter["Chip_2um_3"].tolist()
m1_s = m_pseq_scatter["In_vivo_1"].tolist()
m2_s = m_pseq_scatter["In_vivo_2"].tolist()
m3_s = m_pseq_scatter["In_vivo_3"].tolist()

import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from scipy.stats import gaussian_kde

#Data를 log10 값으로 변경함 : ~에 list화된 비교할 sample data를 넣음 
x = list(np.log10(~))
y = list(np.log10(~))
#Figure size 지정
plt.rcParams['figure.figsize'] = [8,6]

# plt.scatter 내부 리스트 바꾸기, format 함수 바꾸기
#amp_1 = '1st'
#amp_2 = '2nd'
#amp_3 = '3rd'
plt.scatter(x, y, s = 20, alpha=0.45, color='red')
#plt.show()

#(4, 4)는 나눈 칸 구획으로 변경 가능 --> 5칸씩 분할하고 싶으면 (5, 5)로 변경
hist, xbins, ybins, im  = plt.hist2d(x, y, (4, 4), alpha=.4, cmap=plt.cm.jet)
# 각 셀의 텍스트 표기 : 위치 지정
for i in range(len(ybins)-1):
    for j in range(len(xbins)-1):
        plt.text(xbins[j]+0.55,ybins[i]+0.55, str(int(hist.T[i,j])), fontsize=13, color="w", ha="center", va="center", fontweight="bold")

# x, y의 축의 위치 및 범위 설정 (간격 지정)
plt.xticks(np.arange(-4, ~, 1), fontsize = 12)
plt.yticks(np.arange(-4, ~, 1), fontsize = 12)
plt.title('Enriched peptide sequence comparison between samples', fontsize = 15)
plt.xlabel('Log10(normalized ratio)', fontsize = 13)
plt.ylabel('Log10(normalized ratio)', fontsize = 13)
plt.colorbar()
plt.show()
