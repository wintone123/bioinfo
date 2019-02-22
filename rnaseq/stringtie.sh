for file in SRR*.bam
do
basename=${file:0:6}
# check if file exists
if [[ -f $basename.gtf ]]
    then
        echo "$basename.gtf"
        echo "======FILE EXIST!======"
        continue
    else
        # 1st round alignment 
        (set -x; stringtie -p 2 -G chr.gff3 -l $basename -o $basename.gtf $file)
        echo "======FINISH!======"
    fi
done

name_string=""
break=" "
# make input file name
for file in SRR*.gtf
do
    name_string=$file$break$name_string
done
# merge gtf file
(set -x; stringtie --merge -p 2 -G chr.gff3 -o merge.gtf $name_string)
echo "======FINISH!======"

for file in SRR*.bam
do
basename=${file:0:6}
# check if file exists
if [[ -f ballgown/$basename/$basename.gtf ]]
    then
        echo "ballgown/$basename/$basename.gtf"
        echo "======FILE EXIST!======"
        continue
    else
        # 2nd round alignment
        (set -x; stringtie -e -B -p 2 -G merge.gtf -o ballgown/$basename/$basename.gtf $file)
        echo "======FINISH======"
    fi
done

echo "======DONE!======"