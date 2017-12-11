#!/bin/bash

module load bedtools/2.17.0
endSeqDir='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools'
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src:$PATH
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/utils:$PATH

if [[ $# -eq 0 ]] ; then
    echo "Please provide Samples_File and counts_Matrix files as first and second argument, respectively"
    exit 0
fi

case "$1" in
    1) echo 'you gave 1' ;;
    *) echo 'OK' ;;
esac


SamplesFile=$1
CountMatrix=$2

NoS=`grep fastq $SamplesFile|wc -l`

while read line
do
#SAMPLEfastq=`head -$x $SamplesFile|sed 's/\//\t/g' |awk '{print $NF}'`

SAMPLE_NAME=${line%.fastq}

cut -f1-6 $CountMatrix|bedtools intersect -s -wa -wb -a - -b "$SAMPLE_NAME"_Atrimmed_sorted_stranded_read_count.bed|cut -f7-|sort -k1,1 -k2,2n > "$SAMPLE_NAME"_TruePAS_Dinamic.bed;
python $endSeqDir/endSeq_tools.py makeBedgraph -bd "$SAMPLE_NAME"_TruePAS_Dinamic.bed;

rm "$SAMPLE_NAME"_TruePAS_Dinamic.bed;

done < $SamplesFile

module unload bedtools/2.17.0
