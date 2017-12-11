#!/bin/bash

module load java/1.7.0_25
module load samtools/0.0.19
module load R/3.0.1
module load bedtools/2.17.0
module load python/2.7.5


endSeqDir='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools'
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src:$PATH
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/utils:$PATH

utilsdir='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/utils'
WORKINGDIR=$@
#This directory needs to contain "file_contains_sample_fastqs.txt" as well as a file for condition assignments of corresponding samples "file_contains_sample_conditions.txt".

#make a special bed format. 
#Note: this merge_bed_file uses a compressed version of bed file (such as first 3 columns are merged.)

while read line
do
Fastq=$line
basename=${Fastq%.fastq}

awk '{print $1"_"$2"_"$3"_"$6"\t"$5 }' "$basename"_Atrimmed_sorted_stranded_CPM_TruePeaks_3UTR.bed|uniq > \
"$basename"_Atrimmed_sorted_stranded_CPM_TruePeaks_3UTR.bed.tmp

echo "$basename"_Atrimmed_sorted_stranded_CPM_TruePeaks_3UTR.bed.tmp >> $WORKINGDIR/samples.tmp

done < $WORKINGDIR/file_contains_sample_fastqs.txt

x=`wc -l $WORKINGDIR/samples.tmp|awk '{print $1}'`
echo "Number of samples:" $x

c=$(for a in `seq $x`;do let b=$a+6; printf ",""%s"$b ; done)
o=$(for a in `seq $x`;do printf ",""sum" ; done)

let g=$x+7

$utilsdir/merge_bed_files.pl $WORKINGDIR/samples.tmp | \
tail -n+2 |\
awk -F _ '{print $1"\t"$2"\t"$3"\t""ID""\t""1""\t"$4}'|\
sort -k1,1 -k2,2g |\
clusterBed -s -d 40 |\
awk '{OFS="\t"; $(NF)=$1"_"$(NF); print $0}' |\
groupBy -g $g -c 1,2,3,4,5,6$c -o first,min,max,concat,sum,first$o |\
awk '{OFS="\t"; $5=$1; print $0}' |\
cut -f2- > $WORKINGDIR/mergedPeaks.bed

echo "Done"

$endSeqDir/src/annotate.sh $WORKINGDIR/mergedPeaks.bed human 1000

#this t if the column in which 3UTR annotation isoform id is stored. for ex: NM_153254_utr3_8_0_chr1_1120523_f
let t=$x+10
#parse the id and only get NM number: NM_153254
cut -f$t $WORKINGDIR/mergedPeaks.bed_3UTR.bed |cut -d "_" -f1-2 > $WORKINGDIR/nms.tmp

#Look up which NM number associates with which Gene symbol and store this in a temp file
awk 'FNR==NR{a[$1]=$0;next}{print $0"\t", a[$1]}' \
$endSeqDir/data/Annotation/human/hg19genome_reference.rsem.geneNames_transcriptsNames_Len_GCcontent_for_isoforms.txt $WORKINGDIR/nms.tmp |\
cut -f1,5 > $WORKINGDIR/genes.tmp

#paste this with the original annotated bed file
let f=$x+6
cut -f1-$f $WORKINGDIR/mergedPeaks.bed_3UTR.bed|paste -d "\t" - $WORKINGDIR/genes.tmp > $WORKINGDIR/mergedPeaks.bed_3UTR_withGeneNames.bed

rm $WORKINGDIR/*.tmp






#Make a file for statistical testing




