#Trimming 후 bowtie2를 사용해서 sequence alignment 정보를 담은 SAM 파일로 변경함
#PhageM13 amplicon fasta 파일로 bowtie2를 사용해 indexing을 하고 나서 수행할 것 (.1.bt2 ~ .4.bt2, .rev.1.bt2, .rev.2.bt2 파일이 생성되야 함)

for FQ1 in $(ls *_trim_1P)
do
  FQ2=${FQ1/_trim_1P/_trim_2P}
  OUT=${FQ1/_trim_1P}".PhageM13.sam"
  echo $OUT
  bowtie2 --local -x PhageM13 -1 $FQ1 -2 $FQ2 -S $OUT
done
