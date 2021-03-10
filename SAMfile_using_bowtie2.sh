#Trimming 후 bowtie2를 사용해서 sequence alignment 정보를 담은 SAM 파일로 변경함
for FQ1 in $(ls *_trim_1P)
do
  FQ2=${FQ1/_trim_1P/_trim_2P}
  OUT=${FQ1/_trim_1P}".PhageM13.sam"
  echo $OUT
  bowtie2 --local -x PhageM13 -1 $FQ1 -2 $FQ2 -S $OUT
done
