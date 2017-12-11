#!/bin/bash

##$1 = an input bed file to annotate
##$2 = organism (optional, default is human)
##$3 = downstream region length of known 3'UTRs to annotates APAs.
down=$3
dataDIR='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/data'
#If reads are coming from opposite strand of annotated gene;
#We should have an option for other annotation regions
#SAMPLE_NAME=`basename ${1%_Atrimmed_sorted_stranded_CPM_TruePeaks_GB.bed}`
SAMPLE_NAME=`basename ${1%.bed}`

awk 'down='"$down"'{if($6 == "+")print $1"\t"$2"\t"$3+down"\t"$4"\t"$5"\t"$6; else print $1"\t"$2-down"\t"$3"\t"$4"\t"$5"\t"$6}' $dataDIR/Annotation/$2/3UTR.bed|\
bedtools intersect -s -wa -wb -a $1 -b - > ${1%_GB.bed}_3UTR.bed
echo "track name='"$SAMPLE_NAME"_3UTR_TruePeaks' description='User Supplied Track'" > ${1%_GB.bed}_3UTR_GB.bed
cut -f1-9 ${1%_GB.bed}_3UTR.bed >> ${1%_GB.bed}_3UTR_GB.bed

awk '{if($6 == "+")print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $dataDIR/Annotation/$2/5UTR.bed|\
bedtools intersect -s -wa -wb -a $1 -b - > ${1%_GB.bed}_5UTR.bed
echo "track name='"$SAMPLE_NAME"_5UTR_TruePeaks' description='User Supplied Track'" > ${1%_GB.bed}_5UTR_GB.bed
cut -f1-9 ${1%_GB.bed}_5UTR.bed >> ${1%_GB.bed}_5UTR_GB.bed

awk '{if($6 == "+")print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $dataDIR/Annotation/$2/coding_exons.bed|\
bedtools intersect -s -wa -wb -a $1 -b - > ${1%_GB.bed}_coding_exons.bed
echo "track name='"$SAMPLE_NAME"_Exons_TruePeaks' description='User Supplied Track'" > ${1%_GB.bed}_coding_exons_GB.bed
cut -f1-9 ${1%_GB.bed}_coding_exons.bed >> ${1%_GB.bed}_coding_exons_GB.bed

awk '{if($6 == "+")print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $dataDIR/Annotation/$2/Introns_exclusive.bed|\
bedtools intersect -s -wa -wb -a $1 -b - > ${1%_GB.bed}_Introns.bed
echo "track name='"$SAMPLE_NAME"_Introns_TruePeaks' description='User Supplied Track'" > ${1%_GB.bed}_Introns_GB.bed
cut -f1-9 ${1%_GB.bed}_Introns.bed >> ${1%_GB.bed}_Introns_GB.bed


