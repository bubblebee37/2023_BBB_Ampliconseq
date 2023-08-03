#!/bin/bash
for BAM in $(ls *.bam)
do
  echo $BAM
  ./phage_bam-to-nVR_freq.py $BAM
done
