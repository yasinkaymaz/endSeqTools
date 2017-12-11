#!/bin/bash
#BSUB -n 2
#BSUB -R rusage[mem=3000]
#BSUB -R "span[hosts=1]"
#BSUB -q short
#BSUB -W 04:00
#BSUB -o "/project/umw_jeffrey_bailey/yk42w/std_out/%J.out"
#BSUB -e "/project/umw_jeffrey_bailey/yk42w/std_err/%J.err"
#BSUB -J arraytest[1-14]

#TO DO:
#----------
# 1) PROVIDE AN APPROPRIATE FILE "/Full-directory/file_contains_sample_fastqs.txt". This file should contain each fastq file for the samples with a unique name and full directories at each line. 
#	For example; file_contains_sample_fastqs.txt file should look like this;
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N0.fastq
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N4.fastq
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N6.fastq
#	/project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N7.fastq

# 2) Modify the script below; give the correct "/Full-directory/" for "file_contains_sample_fastqs.txt" file.
 
# 3) PLEASE MODIFY arraytest[1-X] BASED ON THE NUMBER OF SAMPLES in "file_contains_sample_fastqs.txt" file. X is the total number of samples. This will allow to submit one job to the GHPCC and run all of the samples in PARALEL independently. The output files will be stored in the same directory of fastq files.

# 4) MODIFY MEMORY (-R) AND JOB TIME (-W) ACCORDING TO YOUR DATA SIZE. If you have more than 10M reads per sample, it would be wise to increase time (-W) to 240:00 (10days) and memory (-R) to 16000 with at least 4 cores (-n). Especially, NBclassifier function consumes a lot of time!

# 5) SUBMIT THIS SCRIPT BY DOING THIS: bsub < PASseq_pipe.sh

#load modules
module load bedtools/2.17.0

endSeqDir='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools'
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src:$PATH
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/utils:$PATH

#Configure
index_id=`echo $LSB_JOBINDEX`
#Fastq=`head -$index_id /project/umw_jeffrey_bailey/OTHERS/LINGTAO/new_bed_nbpred3/ReRunNBpred.txt | tail -1`
Fastq=`head -$index_id /project/umw_jeffrey_bailey/OTHERS/LINGTAO/new_bed_nbpred3/file_contains_sample_fastqs_3.txt | tail -1`

basename=${Fastq%.fastq}
File_DIR=`dirname $Fastq`
echo $Fastq
echo $basename

cd /project/umw_jeffrey_bailey/OTHERS/LINGTAO/new_bed_nbpred3/
mkdir tmp.files/
## Split the bed files into chroms and process in paralel.
for i in `cut -f1 "$basename"_Atrimmed_sorted_stranded_read_count.bed|sort |uniq`;
do 
	cat "$basename"_"$i".bed >>  "$basename"_Atrimmed_sorted_stranded_read_countS.bed;
	tail -n+2 "$basename"_"$i".nbpred >> "$basename"_Atrimmed_sorted_stranded_read_count.nbpred;
	mv "$basename"_"$i".* tmp.files/;
done

mv "$basename"_Atrimmed_sorted_stranded_read_count.bed tmp.files/
mv "$basename"_Atrimmed_sorted_stranded_read_countS.bed "$basename"_Atrimmed_sorted_stranded_read_count.bed
#python $endSeqDir/endSeq_tools.py NBclassifier --bed "$basename"_Atrimmed_sorted_stranded_read_count.bed
#echo "NBC done!"

#python $endSeqDir/endSeq_tools.py clusterPeaks --bed "$basename"_Atrimmed_sorted_stranded_read_count.bed --NBprobs "$basename"_Atrimmed_sorted_stranded_read_count.nbpred
#echo "Clustering done!"
#python $endSeqDir/endSeq_tools.py annotatePeaks --bed "$basename"_Atrimmed_sorted_stranded_CPM_TruePeaks_GB.bed
#echo "Peak annotation done!"

echo "Finished!"

module unload python/2.7.5
module unload java/1.7.0_25
module unload samtools/0.0.19
module unload tophat/2.0.9
module unload bowtie2/2-2.1.0
module unload R/3.0.1
module unload bedtools/2.17.0
module unload fastqc/0.10.1

