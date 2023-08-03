# BBBpeptidePredictor
Version: 2023-08-03 Author: Kyungha Kim (bubblebee37@unist.ac.kr)

Introduction
BBBpeptidePredictor predicts specific BBB-penetrating peptide sequences. The enriched 12-mer length of the peptide is from a combinatorial library of random 12-mer M13 phage display peptide library. The peptide sequence was converted to the one-hot encoded matrix which is used as input of the penetration-predicting ML model based on the binary classification. The BBB-penetration probability and the class of 12-mer peptide could be predicted. The convolutional neural network (CNN) model was fitted using Keras (ver. 2.4.0).

Step 0. Preparation of the input 36-bp DNA sequence dataset from amplicon sequencing data and converting to the one-hot encoded matrix
First, Trim the paired raw FASTQ reads. (/scripts/trim_seq.sh) Use conda to install the trimmomatic. (conda -c bioconda install trimmomatic).


Rename the files as "_R[12].raw.fastq" (i.e. SNU201903_Gliblastoma+e_GBL-67-tumor_R1.raw.fastq).
Map the reads with BWA MEM, and make 'sort-index-rmdup-indexed' BAM file using samtools (/scripts/run-bwa_p.sh).

Step 1. 
Make a config directory to run Strelka2 (/scripts/conf-strelka2-somatic.sh).
Run a Strelka2 (/scripts/run-strelka2.sh).
Rename VCF file and compress the file (/scripts/extract-strelka2-vcf.sh).
Step 2. Variant Filtering
Filter InDel from VCF with 3-8 length range
/vcf/filter-vcf-InDel.py
Filter by Read Depth.
/vcf/filter-vcf-DP.py for Strelka2 somatic variants
/vcf/filter-vcf-AD.py for Strelka2 germline variants
Filter by variant scores (optional)
/vcf/filter-vcf-strelka2SomaticEVC.py
Filter by homozygosity (optional)
/vcf/filter-vcf-Homo.py
Step 3. crRNA Design
Make var_flank_seq from vcf and reference fasta.
/crRNA/make-var_flank_seq.py
Make crRNA fasta from var_flank_seq.
/crRNA/var_flank_seq-to-crRNA <.var_flank_seq>
Make alternative short fragments for all variants.
/crRNA/make-var_flank_seq.py
/crRNA/var_flank_seq-to-alt_fasta.py
Run exonerate to check the off-targets
exonerate
exonerate <alt seq from Step 3.3 above>
Select crRNAs without off target candidates
/crRNA/select-crRNA.py
Step 4. Report
Copy /report/make-region_tview_script.py to the directory, and run it with designed crRNA FASTA file. It will generate bash commands to run samtools tview.
/report/make-region_tview_script.py > run.sh
Check the location of BAM files and the referene FASTA file.
Run the script produced by the above step. It will generate tview files for each crRNA target site.
bash run.sh
Make the summary report by running /report/make-somatic-tview-report.txt.
