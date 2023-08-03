# BBBpeptidePredictor
Version: 2023-08-03

Author: Kyungha Kim (bubblebee@unist.ac.kr)

## Introduction
BBBpeptidePredictor is a tool designed to predict specific blood-brain barrier(BBB)-penetrating peptide sequences. These peptides are amplified 12-mer sequences obtained from a combinatorial library of random 12-mer M13 phage display peptides. The peptide sequences are converted into one-hot encoded matrices, which serve as input for the penetration-predicting machine learning model based on binary classification. This model allows the prediction of BBB-penetration probabilities and the corresponding class of the 12-mer peptide.

The prediction model is built using a Convolutional Neural Network (CNN) and implemented using Keras (ver. 2.4.0).

## Step 0. Pre-process the sequencing data
1. Use conda to install the trimmomatic, bowtie2, samtools.
````
    $ conda install -c bioconda trimmomatic
    $ conda install -c bioconda bowtie2
    $ conda install -c bioconda samtools
````

2. Make index file to search the variable regions. (fasta file format)
To index the variable region in the phage amplicon library, use the flanking region sequence around the variable region (for flanking sequence confirmation, please refer to the protocol or manual of the specific phage library product used).

* i.e.
```fasta
>PhageM13_flank5
ACCGATACAA TTAAAGGCTC CTTTTGGAGC CTTTTTTTTG GAGATTTTCA ACGTGAAAAA ATTATTATTC GCAATTCCTT TAGTGGTACC TTTCTATTCT CACTCT
>PhageM13_flank3
GGTGGAGG TTCGGCCGAA ACTGTTGAAA GTTGTTTAGC AAAATCCCAT ACAGAAAATC ATTACTAACG TCTGGAAAGA CGACAA
```

3. Rename the fastq files as "a_[12].fastq" and "b_[12].fastq" for duplicates.
* i.e. Sample-a_1.fastq, Sample-a_2.fastq for duplicate 1 and Sample-b_1.fastq, Sample-b_2.fastq for duplicate 2

4. Run the bash files in '/dataseq_prep/' sequentially by 02.merge-sam.sh  Through this process, sequence trimming is followed by alignment and file merging using bowtie2 and samtools. One bam file from duplicates is created. (i.e. Sample.PhageM13.bam)
   
5. Download the python codes ('/dataseq_prep/phage_bam-to-nVR_freq.py' & '/dataseq_prep/translate-nVR_freq.py')

6. Next, execute the codes using '/dataseq_prep/03.make-nVR_freq.sh' and '/dataseq_prep/04.make-pVR_freq.sh'. Add the flanking region indexing file path to the bash files.
* Output of 03.make-nVR_freq.sh : sample_name.nVR_freq (Count of identified DNA sequence)
* Output of 04.make-pVR_freq.sh : sample_name.pVR_freq (After translating DNA sequences into peptide sequences, counts and frequencies are displayed.)

## Step 1. Make the input dataset composed of filtered peptide sequences
This step is to remove the outlier sequences and make the input dataset composed of 12-mer target sequences.
Add the '.pVR_freq' file path to the jupyter notebook code, '/dataseq_prep/~.ipynb'. Outlier sequences in all samples are filtered out. Finally, you can make the input dataset to predict the penetration of peptide sequences.

## Step 2. Predict the peptide penetration through the BBB using the CNN model
This step is to load the input dataset to a finalized 10-fold cross-validated CNN model and predict the BBB penetration probability of peptide sequences.
By Add the '.csv' input file path to '/ML_pred/~.ipynb'. The peptide sequences can be converted to the one-hot encoded matrix format and the predicted class and score of each peptide would be derived.
* Class 1 means penetrating peptides.
* Class 0 means Not-penetrating peptides.
