#먼저, fastq.gz 확장명인 raw file을 따로 다른 폴더에 옮겨둔다.

for GZIP in $(ls *.fastq.gz)
do
        FASTQ=${GZIP/.fastq.gz}".fastq"
        echo $FASTQ
        gzip -d $GZIP
done
