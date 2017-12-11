#!/bin/bash

## $1=cluster distance cs
## $2=inputBed
## $3 NB predictions file (*.nbpred)
## outputClusterBed (this should be same name as the input name)
SAMPLE_NAME=`basename ${2%_Atrimmed_sorted_stranded_read_count.bed}`

array=( '<' '>' )
for x in "${array[@]}"
do

if [ "$x" = ">" ]
then
#outPeaks=${2:0:-4}.TruePeaks
outPeaks=${2%.bed}.TruePeaks
#outPeaksBed=${2:0:-15}_CPM_TruePeaks_GB.bed
outPeaksBed=${2%_read_count.bed}_CPM_TruePeaks_GB.bed
#outPeaksBg=${2:0:-15}_CPM_TruePeaks_GB.bedgraph
outPeaksBg=${2%_read_count.bed}_CPM_TruePeaks_GB.bedgraph
echo "Calling True peaks!"
else
#outPeaks=${2:0:-4}.FalsePeaks
outPeaks=${2%.bed}.FalsePeaks
#outPeaksBed=${2:0:-15}_CPM_FalsePeaks_GB.bed
outPeaksBed=${2%_read_count.bed}_CPM_FalsePeaks_GB.bed
#outPeaksBg=${2:0:-15}_CPM_FalsePeaks_GB.bedgraph
outPeaksBg=${2%_read_count.bed}_CPM_FalsePeaks_GB.bedgraph
echo "Calling False peaks!"
fi

awk '{x1=$1; x2=$2; x3=$3; x4=$4; gsub("\"","",x1); gsub("\"","",x2);  gsub("\"","",x3); gsub("\"","",x4); print x1"\t"x2"\t"x3"\t"x4}' $3|\
tail -n+2|\
paste -d "\t" $2 - |\
cut -f1-6,10|\
clusterBed -s -d $1 -i stdin|\
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"_"$8"\t"$7}'|\
groupBy -g 7 -c 1,2,2,3,3,4,5,5,5,5,5,5,5,6,8 -o collapse,collapse,min,collapse,max,collapse,collapse,max,min,sum,median,mean,stdev,collapse,mean|expandCols -c 2,3,5,7,8,15 |\
awk '{if($16 '$x' 0.05) print $0}'|\
awk '{if($8 == $9)printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.02f\t%.02f\t%s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15)}'|\
awk '{print $2"\t"$4"\t"$6"\t"$7"\t"$8"\t"$15"\t"$3"\t"$5"\t"$1"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$6-$4}'|\
sort -k1,1 -k2,2g |\
groupBy -g 9 -c 1,2,3,4,5,6,7,8,10,11,12,13,14,15,16 -o first,first,first,first,first,first,first,last,first,first,first,first,first,first,first |\
awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$5"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' > $outPeaks

total_c=`awk '{x=x+$5}END{print x}' $outPeaks`

echo "track name='"$SAMPLE_NAME"' description='Bed format'" > $outPeaksBed
cut -f1-8 $outPeaks|\
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5*1000000/'$total_c'"\t"$6"\t"$7"\t"$8"\t""255,0,0"}' >> $outPeaksBed

echo "browser position chr8:128748315-128753680" > $outPeaksBg
echo "browser hide all" >> $outPeaksBg
echo "browser pack refGene encodeRegions" >> $outPeaksBg
echo "browser full altGraph" >> $outPeaksBg
echo "browser full altGraph" >> $outPeaksBg
echo "track type=bedGraph name=$SAMPLE_NAME description='BedGraph format' visibility=full color=200,100,0 altColor=0,100,200 priority=20" >> $outPeaksBg
awk '{x1=$1; x2=$2; x3=$3; x4=$4; gsub("\"","",x1); gsub("\"","",x2);  gsub("\"","",x3); gsub("\"","",x4); print x1"\t"x2"\t"x3"\t"x4}' $3|\
tail -n+2|\
paste -d "\t" $2 -|\
cut -f1-6,10|\
clusterBed -s -d $1 -i stdin|\
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"_"$8"\t"$7}'|\
groupBy -g 7 -c 1,2,2,3,3,4,5,5,5,5,5,5,5,6,8 -o collapse,collapse,min,collapse,max,collapse,collapse,max,min,sum,median,mean,stdev,collapse,mean|\
awk '{if($16 '$x' 0.05) print $0}'|\
expandCols -c 2,3,5,7,8,15 |\
awk '{if($15 == "+")print $2"\t"$3"\t"$5"\t"$8*1000000/'$total_c'; else print $2"\t"$3"\t"$5"\t"$15$8*1000000/'$total_c'}' >> $outPeaksBg

#run peak stats
#mkdir peak_stats #give correct directory
#histogram of true/false peak distrubutions

done
