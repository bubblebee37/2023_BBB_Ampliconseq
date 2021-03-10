NUM_THREADS=4
DIR_CONDA=$HOME/anaconda3
#사용하는 adapter에 따라 다른 fasta 파일 넣기 (Nextera, Truseq 등)
FA_ADAPTER=$(find $DIR_CONDA | grep -m1 NexteraPE-PE.fa)
echo "Adapter: "$FA_ADAPTER
cp $FA_ADAPTER .

#아래에 파일경로 지정 (~ 부분)
for FQ1 in $(ls /home/kyungha/Downloads/~/*R1.fastq.gz)
do
        FQ2=${FQ1/_R1/_R2}
        OUT=${FQ1/_R1.fastq.gz}"_trim"
        echo $FQ1 $FQ2 $OUT
#ILLUMINACLIP: 뒤에 해당되는 adapter 파일명 넣기, 그외 pair-end (PE) 또는 single-end (SE)에 따라 다른 option 넣기
        trimmomatic PE -validatePairs -threads $NUM_THREADS -summary $OUT".summary" $FQ1 $FQ2 -baseout $OUT \
   ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
done
