for file in *.fastq.gz
do
# get file basename
basename=${file:0:10}
mkdir $basename
# run FastQC
fastqc -o $basename $file
echo "======FINISH!======"
done