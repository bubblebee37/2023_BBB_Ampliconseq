#!/bin/bash
for SAM1 in $(ls *a.PhageM13.sam)
do
  SAM2=${SAM1/a.PhageM13/b.PhageM13}
  MERGED=${SAM1/a.PhageM13.sam/}".PhageM13.bam"
  echo $MERGED
  samtools merge $MERGED $SAM1 $SAM2
done
