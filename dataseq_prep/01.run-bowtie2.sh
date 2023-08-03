#!/bin/bash
for FQ1 in $(ls *_trim_1P)
do
  FQ2=${FQ1/_trim_1P/_trim_2P}
  OUT=${FQ1/_trim_1P}".PhageM13.sam"
  echo $OUT
  bowtie2 --local -x PhageM13_amplicon -1 $FQ1 -2 $FQ2 -S $OUT
done
