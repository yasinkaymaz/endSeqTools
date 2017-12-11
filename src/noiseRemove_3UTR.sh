module load bedtools/2.17.0
module unload R/3.1.0
module load R/3.0.1
dataDIR='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/data'
srcDIR='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src'

##$1= input bed reads_count.bed
##$2= reads_count.nbpred
##$3=cluster distance cs
clusterdistance=$3

##$4 = organism (optional, default is human)
species=$4

SAMPLE_NAME=`basename ${1%_Atrimmed_sorted_stranded_read_count.bed}`

#Make a bedgraph for visualization of GenomeBrowser
rm ${1%.bed}.bedgraph.gz
TrackNAME=$SAMPLE_NAME'_noisy'
echo "browser position chr8:128748315-128753680" > ${1%.bed}.bedgraph
echo "browser hide all" >> ${1%.bed}.bedgraph
echo "browser pack refGene encodeRegions" >> ${1%.bed}.bedgraph
echo "browser full altGraph" >> ${1%.bed}.bedgraph
echo "browser full altGraph" >> ${1%.bed}.bedgraph
echo "track type=bedGraph name="$TrackNAME" description='"$TrackNAME"_TP' visibility=full color=200,100,0 altColor=0,100,200 priority=20" >> ${1%.bed}.bedgraph
awk '{print $1"\t"$2"\t"$3"\t"$6$5}' $1 >> ${1%.bed}.bedgraph
gzip ${1%.bed}.bedgraph

### Combine read count file and NB predictions file. Assign each site to a gene name. Calculate mean count per gene.
#tail -n+2|sed 's/\"//g'|\

#awk '{x2=$2; gsub("\"","",x2); print x2}' $2 |\
cut -f2 $2 |\
paste -d "\t" $1 - |\
bedtools intersect -s -wa -wb -a - -b $dataDIR/Annotation/$species/genespace7 |\
sort -k11,11 |\
groupBy -g 11 -c 1,2,3,4,5,6,7,8,9,10,11,5 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,mean |\
expandCols -c 2,3,4,5,6,7,8,9,10,11,12 |\
cut -f 2-8,12,13 > ${1%.bed}_genespace_lambda.bed
### Use mean count per gene as the lambda for Poisson test. Calculate pvalues of each site. Filter background noise using Poisson distribution. Filter out sites that are possibly Internal Priming Events.
Rscript $srcDIR/ppois.R ${1%.bed}_genespace_lambda.bed ${1%.bed}_genespace_lambda_ppois.bed

#Make a bedgraph for visualization of GenomeBrowser
rm ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph.gz
TrackNAME=$SAMPLE_NAME'_Noise_IPE_removed'
echo "browser position chr8:128748315-128753680" > ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph
echo "browser hide all" >> ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph
echo "browser pack refGene encodeRegions" >> ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph
echo "browser full altGraph" >> ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph
echo "browser full altGraph" >> ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph
echo "track type=bedGraph name="$TrackNAME" description="$TrackNAME" visibility=full color=200,100,0 altColor=0,100,200 priority=20" >> ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph
awk '{ print $1"\t"$2"\t"$3"\t"$6$5}' ${1%.bed}_genespace_lambda_ppois.bed >> ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph
gzip ${1%.bed}_genespace_lambda_ppois_filtered.bedgraph

#Group PA sites in a close proximity and sum their counts.
clusterBed -s -d $clusterdistance -i ${1%.bed}_genespace_lambda_ppois.bed |\
awk '{$10=$1"_"$10"_"$6; OFS="\t"; print $0}'|\
groupBy -g 10 -c 1,2,3,5,6,8 -o first,min,max,sum,first,first |\
awk '{OFS="\t"; print $2,$3,$4,$1,$5,$6,$7}' > ${1%.bed}_genespace_lambda_ppois_groupedSites.bed

rm ${1%.bed}_genespace_lambda_ppois_groupedSites_GB.bed.gz
TrackNAME=$SAMPLE_NAME'_PAsites'
echo "track name='"$TrackNAME"' description='Bed format'" > ${1%.bed}_genespace_lambda_ppois_groupedSites_GB.bed
cut -f1-6 ${1%.bed}_genespace_lambda_ppois_groupedSites.bed >> ${1%.bed}_genespace_lambda_ppois_groupedSites_GB.bed
gzip ${1%.bed}_genespace_lambda_ppois_groupedSites_GB.bed

#Calculate the percent usage of each PA sites within each gene.
sort -k7,7 ${1%.bed}_genespace_lambda_ppois_groupedSites.bed|\
groupBy -g 7 -c 1,2,3,4,5,6,7,5 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,sum|\
expandCols -c 2,3,4,5,6,7,8|\
cut -f2-|\
awk '{print $0"\t"(100*$5/$8)}'|\
awk '{OFS="\t"; print $1,$2,$3,$4,$9,$6}' > ${1%.bed}_genespace_lambda_ppois_groupedSites_PCT.bed

# Create stats log file:
SitesFound=`wc -l ${1%.bed}_genespace_lambda_ppois_groupedSites.bed|awk '{print $1}'`
Genes=`cut -f7 ${1%.bed}_genespace_lambda_ppois_groupedSites.bed|sort|uniq|wc -l`
total_c=`awk '{x=x+$5}END{print x}' ${1%.bed}_genespace_lambda_ppois_groupedSites.bed`
echo "Total number of PA sites found is $SitesFound" > ${1%.bed}_PA_sites.stats
echo "Total number of Genes to which these PA sites assigned to is $Genes" >> ${1%.bed}_PA_sites.stats
echo "Library Size (Total read counts of PA sites) is $total_c" >> ${1%.bed}_PA_sites.stats

#This is how to run this script
#/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src/noiseRemove_3UTR.sh N0_Atrimmed_sorted_stranded_read_count.bed N0_Atrimmed_sorted_stranded_read_count.nbpred 20 human

