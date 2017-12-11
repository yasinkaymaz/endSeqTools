#!/bin/bash
#BSUB -n 2
#BSUB -R rusage[mem=16000]
#BSUB -R "span[hosts=1]"
#BSUB -q long
#BSUB -W 8:00
#BSUB -o "/project/umw_jeffrey_bailey/OTHERS/std_out/%J.out"
#BSUB -e "/project/umw_jeffrey_bailey/OTHERS/std_err/%J.err"


module load java/1.7.0_25
module load samtools/0.0.19
module load tophat/2.0.9
module load bowtie2/2-2.1.0
module load R/3.0.1
module load bedtools/2.17.0
module load python/2.7.5
module load fastqc/0.10.1

endSeqDir='/home/yk42w/codes/endSeq_Tools'
export PATH=/home/yk42w/codes/endSeq_Tools/src:$PATH
export PATH=/home/yk42w/codes/endSeq_Tools/utils:$PATH
#TO DO:
#----------
# 1) PROVIDE AN APPROPRIATE FILE "/Full-directory/file_contains_sample_fastqs.txt". This file should contain each fastq file for the samples with a unique name and full directories at each line.
#	For example; file_contains_sample_fastqs.txt file should look like this;
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N0.fastq
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N4.fastq
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N6.fastq
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N7.fastq


# 2) Modify the script below; give the correct "/Full-directory/" for "file_contains_sample_fastqs.txt" file as well as *.gtf file for ex; gencode.v19.annotation.gtf.
# 3) SUBMIT THIS SCRIPT BY DOING THIS: bsub < MakeSampleTable.sh

#python $endSeqDir/endSeq_tools.py makeSamplesTable -sf /Full-directory/file_contains_sample_fastqs.txt --gtf /Full-directory/gencode.v19.annotation.gtf
DIR='/project/umw_jeffrey_bailey/OTHERS/LINGTAO/new_bed_nbpred3'

python $endSeqDir/endSeq_tools.py makeSamplesTable -sf $DIR/file_contains_sample_fastqs_YK.txt --gtf /project/umw_jeffrey_bailey/share/Homo_sapiens/Gencode/v19/gencode.v19.annotation.gtf

python $endSeqDir/endSeq_tools.py annotatePeaks -bd mergedPeaks_CPM.bed

#Rscript $endSeqDir/utils/PASseq_chisquire_beta_v4.R $DIR/mergedPeaks_annotated_CPM-3.bed $DIR/switch_Test.txt 6 8 $DIR

#In summary, this function will
#1) Collect all non-zero genomic locations from all samples into a pool,
#2) Test each location if significantly high expressed enough to be distinguished from background noise using a Poisson model (p-val cutoff 0.05),
#3) Eliminate possible Internal priming events based on Naive Bayes classifier probability function (cut off for classification 0.001),
#4) Exclude PA sites of which cumulative read counts do not exceed to 100 reads (this is sum of all samples at that particular site so the threshold is not that high sample wise)
#5) Normalize each site expression by the total number of reads comprising real PA sites (after filtrations) in each sample separately (Counts Per Million),
#6) Calculate percent PA site usage within gene for each PA site.
#7) Output sample matrices in three different format after annotating with gene names;
#        mergedPeaks_annotated_Count.bed
#        mergedPeaks_annotated_CPM.bed
#        mergedPeaks_annotated_Percent.bed
