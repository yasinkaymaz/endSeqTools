#!/bin/bash
#BSUB -n 2
#BSUB -R rusage[mem=16000]
#BSUB -R "span[hosts=1]"
#BSUB -q long
#BSUB -W 240:00
#BSUB -o "/project/umw_jeffrey_bailey/OTHERS/std_out/%J.out"
#BSUB -e "/project/umw_jeffrey_bailey/OTHERS/std_err/%J.err"
#BSUB -J arraytest[1-20]

#TO DO:
#----------
# 1) PROVIDE AN APPROPRIATE FILE "/Full-directory/file_contains_sample_fastqs.txt". This file should contain each fastq file for the samples with a unique name and full directories at each line.
#       For example; file_contains_sample_fastqs.txt file should look like this;
#       /project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N0.fastq
#       /project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N4.fastq
#       /project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N6.fastq
#       /project/umw_jeffrey_bailey/OTHERS/Ami/rerun/N7.fastq

# 2) Modify the script below; give the correct "/Full-directory/" for "file_contains_sample_fastqs.txt" file.

# 3) PLEASE MODIFY arraytest[1-X] BASED ON THE NUMBER OF SAMPLES in "file_contains_sample_fastqs.txt" file. X is the total number of samples. This will allow to submit one job to the GHPCC and run all of the samples in PARALEL independently. The output files will be stored in the same directory of fastq files.

# 4) MODIFY MEMORY (-R) AND JOB TIME (-W) ACCORDING TO YOUR DATA SIZE. If you have more than 10M reads per sample, it would be wise to increase time (-W) to 240:00 (10days) and memory (-R) to 16000 with at least 4 cores (-n). Especially, NBclassifier function consumes a lot of time!

# 5) SUBMIT THIS SCRIPT BY DOING THIS: bsub < PASseq_pipe_mouse.sh

#load modules
module load java/1.7.0_25
module load samtools/0.0.19
module load tophat/2.0.9
module load bowtie2/2-2.1.0
module load R/3.0.1
module load bedtools/2.17.0
module load python/2.7.5

endSeqDir='/home/yk42w/codes/endSeqTools'
export PATH=/home/yk42w/codes/endSeqTools/src:$PATH
export PATH=/home/yk42w/codes/endSeqTools/utils:$PATH
PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/utils:$PATH

#Configure
index_id=`echo $LSB_JOBINDEX`
Fastq=`head -$index_id /project/umw_jeffrey_bailey/OTHERS/Ami/rerun/file_contains_sample_fastqs.txt | tail -1`
basename=${Fastq%.fastq}
echo $Fastq
echo $basename
#PASseq pipeline
python $endSeqDir/endSeq_tools.py trimPolyAstrecth --fastq $Fastq
echo "Reads trimmed!"
python $endSeqDir/endSeq_tools.py genomeAlign --fastq "$basename"_Atrimmed.fastq --mouse
echo "Alignment done!"
python $endSeqDir/endSeq_tools.py samProcess --sam "$basename"_Atrimmed.sam
echo "Sam processed!"
python $endSeqDir/endSeq_tools.py genomeCov --bam "$basename"_Atrimmed_sorted.bam --re 3 --mouse
echo "Coverage calculated!"
python $endSeqDir/endSeq_tools.py NBclassifier --bed "$basename"_Atrimmed_sorted_stranded_read_count.bed --mouse
echo "NBC done!"
python $endSeqDir/endSeq_tools.py clusterPeaks --bed "$basename"_Atrimmed_sorted_stranded_read_count.bed --NBprobs "$basename"_Atrimmed_sorted_stranded_read_count.nbpred --mouse
echo "Clustering done!"
python $endSeqDir/endSeq_tools.py annotatePeaks --bed "$basename"_Atrimmed_sorted_stranded_CPM_TruePeaks_GB.bed --mouse
echo "Peak annotation done!"

echo "Finished!"

module unload python/2.7.5
#module unload python3/3.3.2
module unload java/1.7.0_25
module unload samtools/0.0.19
module unload tophat/2.0.9
module unload bowtie2/2-2.1.0
module unload R/3.0.1
module unload bedtools/2.17.0
