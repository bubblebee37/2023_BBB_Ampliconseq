#SAM 파일 내에서 각 조건 별 replicate를 a, b를 붙여서 파일명 변경할 것 --> replicate를 통합한 BAM 
for SAM1 in $(ls *a.PhageM13.sam)
do
  SAM2=${SAM1/a.PhageM13/b.PhageM13}
  MERGED=${SAM1/a.PhageM13.sam/}".PhageM13.bam"
  echo $MERGED
  samtools merge $MERGED $SAM1 $SAM2
done
