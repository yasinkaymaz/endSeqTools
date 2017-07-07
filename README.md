# endSeqTools
A pipeline for end sequencing data sets. 

USAGE: 
endSeqTools.py [-h] subcommand [suboptions] 
-----
DESCRIPTION: 
A collection of tools to analyze end sequencing data such as PASseq, 3Pseq etc.
-----------
SUBCOMMANDS:
-----------
makeSamplesTable:   Combines all bed files and returns a table with peaks.

NBclassifier:   Detecting internally priming events with NB classifier.

trimPolyTstrecth:   This function trims polyT stretch that is present in the sequence reads pulled down with oligo-dA.

trimPolyAstrecth:   This function trims polyA stretch that is present in the sequence reads pulled down with oligo-dT.

switchTest:   Performs PAS switch test using called peaks.

makeBedgraph:   Creates a bedgraph file for a given bed file.

clusterPeaks:   Cluster peak locations.

genomeCov:   Calculates coverage levels of peaks in a strand specific manner.

genomeAlign:   Aligns reads to genomic sequence with Bowtie2 or BWA.

samProcess:   Processing of sam alignment file: Sam2Bam, sortBam, indexBam.

annotatePeaks:   Annotates peaks with genomic locations such as 3'UTR and returns only annotated intervals.

HELP:
----
-h/-help   short / long  subcommand descriptions 

For specific options: endSeqTools.py [subcommand] --help 


