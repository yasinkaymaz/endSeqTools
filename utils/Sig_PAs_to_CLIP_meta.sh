#!/bin/bash

module load bedtools/2.22.0
module load samtools/0.0.19
module load python/2.7.9
export PATH=/home/yk42w/project/share/bin_sync/bbmap:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/deepTools2.0/bin:$PATH
export PYTHONPATH=/project/umw_jeffrey_bailey/share/bin_sync/deepTools2.0/lib/python2.7/site-packages/
chromsizesfile='/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/chromSizes.txt'

SigPAFile=$1
#SigPAFile=Significantly_DE_pA_counts_between_shC_shTDP_padj0.01.txt
#Modified_PACountTable=$2
Modified_PACountTable=mergedPeaks_annotated_Count_Modified.bed

#CLIPrep1=340_01.basedon_340_01.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.encode.bed
#CLIPrep2=340_02.basedon_340_02.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.encode.bed

CLIPrep1=K562_TDP-CLIP_hg19_ENCFF785PDJ.bed
CLIPrep2=K562_TDP-CLIP_hg19_ENCFF489VNW.bed


for clip in $CLIPrep1 $CLIPrep2;
do

if [ -f ${clip%.bed}_Rv_peaks_Filtered.bigWig ]
then
	echo "bigwig files are present"
else

awk '{OFS="\t"; if($6 == "+")  print $1,$2,$3,$8}' $clip |\
sort -k1,1 -k2,2g |\
bedtools merge  -c 4 -o sum -i - |awk '($4 > 1)' > ${clip%.bed}_Fw_peaks_Filtered.bedgraph

/project/umw_jeffrey_bailey/share/bin_sync/ucsc/bedGraphToBigWig \
${clip%.bed}_Fw_peaks_Filtered.bedgraph \
$chromsizesfile \
${clip%.bed}_Fw_peaks_Filtered.bigWig

awk '{OFS="\t"; if($6 == "-")  print $1,$2,$3,$8}' $clip |\
sort -k1,1 -k2,2g |\
bedtools merge  -c 4 -o sum -i - |awk '($4 > 1)' > ${clip%.bed}_Rv_peaks_Filtered.bedgraph

/project/umw_jeffrey_bailey/share/bin_sync/ucsc/bedGraphToBigWig \
${clip%.bed}_Rv_peaks_Filtered.bedgraph \
$chromsizesfile \
${clip%.bed}_Rv_peaks_Filtered.bigWig

fi

done


for DE in UP DOWN;
do

for strand in Fw Rv;
do

if [ "$DE" == "UP" ]
then
	if [ "$strand" == "Fw" ]
	then

	echo $strand $DE
	tail -n+2 $SigPAFile |\
	awk '($3 > 0 )' |\
	cut -f1|\
	sed 's/"//g'|\
	awk 'FNR==NR{a[$1]=$0;next}{print $0"\t", a[$1]}' $Modified_PACountTable -|\
	awk '($8 == "+")' |\
	awk '{OFS="\t"; print $3,$4,$5+1,$1,"1",$8}' > "${SigPAFile}"_"$strand"_$DE.bed
	else

	echo $strand $DE
	tail -n+2 $SigPAFile |\
	awk '($3 > 0 )' |\
	cut -f1|\
	sed 's/"//g'|\
	awk 'FNR==NR{a[$1]=$0;next}{print $0"\t", a[$1]}' $Modified_PACountTable -|\
	awk '($8 == "-")' |\
	awk '{OFS="\t"; print $3,$4,$5+1,$1,"1",$8}' > "${SigPAFile}"_"$strand"_$DE.bed

	fi

else

	if [ "$strand" == "Fw" ]
	then

	echo $strand $DE
	tail -n+2 $SigPAFile |\
	awk '($3 < 0 )' |\
	cut -f1|\
	sed 's/"//g'|\
	awk 'FNR==NR{a[$1]=$0;next}{print $0"\t", a[$1]}' $Modified_PACountTable -|\
	awk '($8 == "+")' |\
	awk '{OFS="\t"; print $3,$4,$5+1,$1,"1",$8}' > "${SigPAFile}"_"$strand"_$DE.bed
	else

	echo $strand $DE
	tail -n+2 $SigPAFile |\
	awk '($3 < 0 )' |\
	cut -f1|\
	sed 's/"//g'|\
	awk 'FNR==NR{a[$1]=$0;next}{print $0"\t", a[$1]}' $Modified_PACountTable -|\
	awk '($8 == "-")' |\
	awk '{OFS="\t"; print $3,$4,$5+1,$1,"1",$8}' > "${SigPAFile}"_"$strand"_$DE.bed

	fi

fi


computeMatrix reference-point \
-S ${CLIPrep1%.bed}_"$strand"_peaks_Filtered.bigWig \
   ${CLIPrep2%.bed}_"$strand"_peaks_Filtered.bigWig \
-R "${SigPAFile}"_"$strand"_$DE.bed \
--referencePoint center \
--beforeRegionStartLength 5000 \
--afterRegionStartLength 5000 \
--binSize 100 \
--outFileName "${SigPAFile}"_CLIP_"$strand"_peaks.mtrx

#profile
plotProfile \
--matrixFile "${SigPAFile}"_CLIP_"$strand"_peaks.mtrx \
--numPlotsPerRow 1 \
--plotWidth 20 \
--plotHeight 10 \
--averageType "mean" \
--regionsLabel "TPDP43-CLIP" \
--refPointLabel "'$strand' PA sites" \
--yAxisLabel "mean Peak intensity" \
--plotTitle "CLIPseq-TDP43 peaks around PA " \
--verbose \
--outFileName "${SigPAFile}"_"$DE"_"${clip%.bed}"_"$strand"_peaks.pdf
#--yMax 2000 \
done
done



