#!/bin/bash

#load modules
module load java/1.7.0_25
module load samtools/0.0.19
module load tophat/2.0.9
module load bowtie2/2-2.1.0
module load R/3.0.1
module load bedtools/2.17.0
module load python/2.7.5
module load fastqc/0.10.1

endSeqDir='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools'
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src:$PATH
export PATH=/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/utils:$PATH

python $endSeqDir/endSeq_tools.py NBclassifier --bed $1
echo "NBC done!"

