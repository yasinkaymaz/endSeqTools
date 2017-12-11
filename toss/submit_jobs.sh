#!/bin/bash

while read line
do

Fastq=$line

bsub ./PASseq_pipe.sh $line

done < ./file_contains_sample_fastqs.txt

