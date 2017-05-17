#!/usr/bin/python

import os.path
import sys
import os, time, subprocess
import textwrap
import argparse
import glob

argparser = argparse.ArgumentParser(description='Process some integers.')

subcommands={}
def main(args):
    if len(args)  == 0 or args[0] in ["h", "help", "-h", "--h", "--help","-help"] :
        verbosity= 'shortDesc'
        if args[0] in ["help" , "--help", "-help"]:
            verbosity = 'longDesc'
        program_name=os.path.basename(__file__)
        print "\n","USAGE:",program_name, "[-h] subcommand [suboptions]","\n", "-----"
        print "DESCRIPTION: A collection of tools to analyze end sequencing data such as PASseq, 3Pseq etc."
        print "-----------"
        print "SUBCOMMANDS:"
        print "-----------"
        for k in subcommands.keys():
            text=subcommands[k][verbosity]
            text= textwrap.dedent(text)
            if text:
                text =  "%s:   %s " %(k, text )
                print textwrap.fill(text, 75, initial_indent='', subsequent_indent='                    ')
        print "\n","HELP:"
        print "----"
        print "-h/-help   short / long  subcommand descriptions","\n"
        print "For specific options:", program_name, "[subcommand] --help", "\n"
    elif args[0] == 'pydoc':
        os.system( "pydoc " + os.path.abspath(__file__) )
    elif args[0] in subcommands.keys():
        globals()[args[0]](args[1:])
    else:
        print "unknown subcommand (" + args[0] + ") use -h for list of subcommands!"
        sys.exit(-1)
    sys.exit(0)
    
    return

def exportTools():
    """
    This is for exporting external tools to the environment.
    """
    #Bowtie
    #BWA
    #Samtools
    #BedTools
    #Rscript
    
    return

def run_command(params,shell=True):
    cmd=' '.join(params)
    proc = subprocess.Popen(cmd,stdin=None,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=shell)
    (odata,edata) = proc.communicate()
    
    if len(edata)>0:
        edata=0
    return [edata,odata]

shortDescText="Trimming polyA tail within sequence reads."
longDescText="""This function trims polyA stretch that is present in the sequence reads pulled down with oligo-dT."""
subcommands['trimPolyAstrecth'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def trimPolyAstrecth(args):
    """
    
    #Input: Fastq file, required
    #Arguments: NA
    """
    #TODO:    1. Export cutadapt to the environment
    argparser.add_argument("-f", "--fastq", required=True, help='Raw read sequences for polyA trimming in fastq format. Please provide with full directory.')
    args=argparser.parse_args(args=args)
    fastq = glob.glob(args.fastq)[0]
    outDir=os.path.dirname(fastq)
    print "input fastq file is: ", glob.glob(args.fastq)[0]
    print "Output directory is: ", outDir
    params = ['cutadapt', '-b', 'AAAAAAAAAAA', '-e', '0.1', '-O', '5', '-o', '%s' %fastq[:-6]+'_Atrimmed.fastq', '%s' %fastq]
    run_command(params,shell=True)
    print "Trimming polyA tails..."
    print "Done!"
    print "Trimmed output fastq file is:", fastq[:-6]+'_Atrimmed.fastq'
    print "Quality checking for trimmed reads..."
    params=['fastqc', '%s' %fastq[:-6]+'_Atrimmed.fastq']
    run_command(params,shell=True)
    print "Done! Please check fastqc folder for quality report."
    return


shortDescText="Trimming polyT tail within sequence reads."
longDescText="""This function trims polyT stretch that is present in the sequence reads pulled down with oligo-dA."""
subcommands['trimPolyTstrecth'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def trimPolyTstrecth(args):
    """
    
    #Input: Fastq file, required
    #Arguments: NA
    """
    #TODO:    1. Export cutadapt to the environment
    argparser.add_argument("-f", "--fastq", required=True, help='Raw read sequences for polyT trimming in fastq format. Please provide with full directory.')
    args=argparser.parse_args(args=args)
    fastq = glob.glob(args.fastq)[0]
    outDir=os.path.dirname(fastq)
    print "input fastq file is: ", glob.glob(args.fastq)[0]
    print "Output directory is: ", outDir
    params = ['cutadapt', '-g', 'NNNNTTTTTTTTTT','-n','3','-g', 'TTTTTTTTTT','-n','10','-g', 'AAAAAAAAAAA','-n','10','-a', 'TTTTTTTTTT','-n','10','-a', 'AAAAAAAAAAA','-n','10', '-b','TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG', '-e', '0.1', '-O', '4', '-m', '10', '-o', '%s' %fastq[:-6]+'_Ttrimmed.fastq', '%s' %fastq]
    run_command(params,shell=True)
    print "Trimming polyT tails..."
    print "Done!"
    print "Trimmed output fastq file is:", fastq[:-6]+'_Ttrimmed.fastq'
    print "Quality checking for trimmed reads..."
    params=['fastqc', '%s' %fastq[:-6]+'_Ttrimmed.fastq']
    run_command(params,shell=True)
    print "Done! Please check fastqc folder for quality report."
    return


shortDescText="Aligns reads to genome."
longDescText=""" Aligns reads to genomic sequence with Bowtie2 or BWA. """
subcommands['genomeAlign'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def genomeAlign(args):
    """
    This function maps reads in Fastq format to a target genome using Bowtie2
    #Input: Fastq file
    #Arguments: -
    #User option: bowtie2, bwa 
    #Output(s): samfile, unmapped reads file
    """
    #TODO: 1. needs bowtie2 exported to environment --> export PATH=/home/baileyj1/share/bin_sync/bowtie2-2.1.0/
    #    2. Provide option for other aligner such as BWA? STAR?
    #    3. change the output directory of excluded un-mapped reads to save space.
    #    4. Include option for other genomes such as mouse. DONE
    #    5. Provide statistics and plots
    #    6. Fix the output directory of unmapped.fastz.gz
    
    argparser.add_argument("-f", "--fastq", required=True, help='Raw read sequences in fastq format for mapping to genome. Please provide with full directory.')
    argparser.add_argument("--mouse", action='store_true', required=False, help='Choose target organism genome as mouse to map read sequences (mm10). Please add this argument if your organism is mouse (Default is human).')
    argparser.add_argument("-m","--multimappers",  action='store', type=int, required=False, dest="MM", help='Choose the limit of reported multimappers, integer (Default is 1).')
    argparser.add_argument("-p","--multithread",  action='store', type=int, required=False, dest="MT", help='Choose the number of processors, integer (Default is 2).')

    args=argparser.parse_args(args=args)
    
    fastq = glob.glob(args.fastq)[0]
    print "Input read file is: ", fastq
    outDir=os.path.dirname(fastq)
    print outDir
    basename=fastq[:-6]
    samFile = basename+'.sam'
    if args.mouse:
        genome='%s/data/Sequence/mouse/Bowtie2Index/genome' %install_dir
    else:
        genome='%s/data/Sequence/human/Bowtie2Index/genome' %install_dir

    if args.MM:
        mm=args.MM
    else:
        mm=1

    if args.MT:
        mt=args.MT
    else:
        mt=2

    print '%(dir)s/unmapped.fastq.gz' %{"dir":outDir}
    params = ['bowtie2','-p','%s' %mt, '-x', '%s' %genome,'-k','%s' %mm, '--un-gz', '%(dir)s/unmapped.fastq.gz' %{"dir":outDir}, '-U', '%(fq1)s' %{"fq1":fastq}, '-S', '%s' %samFile]
    
    run_command(params,shell=True)
    print "Output genome alignment file is: ", samFile
    
    return

shortDescText="Post-alignment processing of sam files."
longDescText=""" Processing of sam alignment file: Sam2Bam, sortBam, indexBam. """
subcommands['samProcess'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }


def samProcess(args):
    """
    This function converts sam to bam (filters the low mapping quality reads), sorts bam, and indexes bam.
    #Input: samfile
    #Arguments: -
    #User options: Thr_MQ, default is 20
    #Output(s): filtered-sorted.bam/bai
    """
    #TODO: 1.Provide statistics and plots
    argparser.add_argument("-s", "--sam", required=True, help='A sequence alignment file in sam format is required. Please provide with full directory.')
    argparser.add_argument("-q", "--mq_threshold", action='store', type=int, required=False, dest="Thr_MQ", help="Mapping quality threshold for filtering non-uniquely mapping reads, integer. (Default is 20)")
    args=argparser.parse_args(args=args)

    samFile=glob.glob(args.sam)[0]
    outDir=os.path.dirname(samFile)
    basename=samFile.strip().split("/")[-1][:-4]
    bamFile = outDir+"/"+basename+'.bam'
    sortedBamFile = outDir+"/"+basename+'_sorted.bam'
   
    if args.Thr_MQ:
        Thr_MQ=args.Thr_MQ
    else:
        Thr_MQ=20
    print "Output directory is: ", outDir
    print "Processing sam file and converting to bam file..."
    print "Mapping quality threshold is: ", Thr_MQ
    params = ['samtools', 'view', '-q', '%s' %Thr_MQ, '-bS', '%s' %samFile, '>', '%s' %bamFile]
    run_command(params,shell=True)
    print "Sorting sam file..."
    params = ['samtools', 'sort', '%s' %bamFile, '%s' %bamFile[:-4]+'_sorted']
    run_command(params,shell=True)
    params = ['samtools', 'flagstat', '%s' %bamFile, '>', '%s' %bamFile[:-4]+'_alignment_stats.txt']
    run_command(params,shell=True)
    print "Indexing sorted bam file..."
    params = ['samtools', 'index', '%s' %sortedBamFile]
    run_command(params,shell=True)
    print "Done!"
    return


def excludeChr():
    """
    #Input: 
    #Arguments: 
    #User options: 
    #Output(s): 
    """
    return

def genomCov2bed(GenCov, strand, Bed):
    """
    #Input: Output file of genomeCoverage from BedTools
    #Arguments: 
    #User options: 
    #Output(s): Bedformat
    """
    fileGenCov = open('%s' %GenCov, "r")
    fileBed = open('%s' %Bed, "w")
    peakid = 0
    for linex in fileGenCov:
        peakid = peakid +1
        end = int(linex.strip().split("\t")[1]) + 1
        if strand == 'minus':
            strd = '-'
        else:
            strd = '+'
        #This step filters single coverage locations.
        if int(float(linex.strip().split("\t")[2])) > 1:
            fileBed.write(linex.strip().split("\t")[0]+"\t"+linex.strip().split("\t")[1]+"\t"+str(end)+"\t"+'%s' %strand+str(peakid)+"\t"+linex.strip().split("\t")[2]+"\t"+'%s' %strd+"\n")
        else:
            pass
    peakid = 0
    fileGenCov.close()
    fileBed.close()

    
    return

shortDescText="Calculates coverage levels of peaks in a strand specific manner."
longDescText=""" Calculates coverage levels of peaks in a strand specific manner. """
subcommands['genomeCov'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def genomeCov(args):
    """
    #Input: 
    #Arguments: 
    #User options: mouse or human (default human). User optionally can provide bam file with bam argument.
    #Output(s):
    """
    #TODO: 1. 
    #      2. Provide statistics and plots
    #      3. 
    #      4. it does not work if you dont provide bam file. Fix!
    #	   5. Give an option to switch strand for RNaseH protocol or 3Pseq etc.
    #Calculate Coverage for + strand
    #Note: Because RNaseH protocol produces reads from opposite strand of original RNA, coverage should be calculated at the 5'-end of the reads.

    argparser.add_argument("-b", "--bam", required=True, help='A sequence alignment file in bam format is required. Please provide with full directory.')
    argparser.add_argument("--mouse", action='store_true', required=False, help='Choose target organism genome as mouse to map read sequences (mm10). Please add this argument if your organism is mouse (Default is human).')
    argparser.add_argument("-re", "--read_end", action='store', type=int, required=True, dest="re", help="Read end site to calculate coverage. For example; because RNaseH protocol produces reads from opposite strand of original RNA, coverage should be calculated at the 5'-end of the reads. (Default is 3 for 3-prime end).")
    args=argparser.parse_args(args=args)

    sortedBamFile = glob.glob(args.bam)[0]
    outDir=os.path.dirname(sortedBamFile)
    basename=sortedBamFile.strip().split("/")[-1][:-4]
    plusStrGenCov = outDir+"/"+basename+'_pos_read_count.gencov'
    minStrGenCov = outDir+"/"+basename+'_min_read_count.gencov'
    plusStrBed = outDir+"/"+basename+'_pos_read_count.bed'
    minStrBed = outDir+"/"+basename+'_min_read_count.bed'
    mergedBed = outDir+"/"+basename+'_stranded_read_count.bed'
    mergedBedtmp = outDir+"/"+basename+'_stranded_read_count.bed.tmp'

    if args.mouse:
        genome='%s/data/Sequence/mouse/Bowtie2Index/genome.fa' %install_dir
    else:
        genome='%s/data/Sequence/human/Bowtie2Index/genome.fa' %install_dir
    re=args.re

    print "Scanning for plus strand..."
    params = ['genomeCoverageBed', '-ibam', '%s' %sortedBamFile, '-g',\
                                '%s' %genome, '-dz', '-%s' %re, '-strand', '+', '>',\
                                '%s' %plusStrGenCov]
    run_command(params,shell=True)
    #Convert output to bed
    strand = 'plus'
    print "Bed file is being created..."
    genomCov2bed(plusStrGenCov, strand, plusStrBed)
    
    #Calculate Coverage for - strand
    print "Scanning for minus strand..."
    params = ['genomeCoverageBed', '-ibam', '%s' %sortedBamFile, '-g',\
                                '%s' %genome, '-dz', '-%s' %re, '-strand', '-', '>',\
                                '%s' %minStrGenCov]
    run_command(params,shell=True)
    #Convert output to bed
    strand = 'minus'
    print "Bed file is being created..."
    genomCov2bed(minStrGenCov, strand, minStrBed)
    
    #merge POS strand NEG strand files
    params = ['cat', '%s' %plusStrBed,\
               '%s' %minStrBed, '>', \
               '%s' %mergedBed]
    run_command(params,shell=True)

    params = ['grep', '-v', 'chrM', '%s' %mergedBed, '>', '%s' %mergedBedtmp]
    run_command(params,shell=True)

    os.remove(mergedBed)
    
    params = ['mv', '%s' %mergedBedtmp, '%s' %mergedBed]
    run_command(params,shell=True)

    print "A merged bed file storing coverages is reported: ", mergedBed
    os.remove(plusStrGenCov)
    os.remove(minStrGenCov)
    os.remove(plusStrBed)
    os.remove(minStrBed)
    
    return

shortDescText="Detecting internally priming events with NB classifier."
longDescText=""" Detecting internally priming events with NB classifier. """
subcommands['NBclassifier'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def NBclassifier(args):
    """
    #Run Naive Bayes classification in order to compute probability of a given PAS peak being an internal priming event. 
    This step runs an external R script in ./src/ called "cleanUpdTSeq.R".
    This step takes too long. Time consuming!
    #Input: 
    #Arguments: 
    #User options: 
    #Output(s):
    """
    #TODO: 1. Provide statistics and plots
    #      2. Needs Rscript/R exported to environment
    #      3. R script misses first interval in the bed file. Fix it.

    argparser.add_argument("-bd", "--bed", required=True, help='A file stores candidate locations in bed format is required. Please provide with full directory.')
    argparser.add_argument("--mouse", action='store_true', required=False, help='Choose target organism genome as mouse to map read sequences (mm10). Please add this argument if your organism is mouse (Default is human).')
    args=argparser.parse_args(args=args)
    
    mergedBed = glob.glob(args.bed)[0]
    outDir=os.path.dirname(mergedBed)
    basename=mergedBed.strip().split("/")[-1][:-4]

    print 'The input bed file for NBc is: %s' %mergedBed
#    NB_output = outDir+"/"+basename+'_stranded_read_count.nbpred'
    NB_output = outDir+"/"+basename+'.nbpred'
    print 'The output file 0f NBc will be: %s' %NB_output

    if args.mouse:
        organism='mouse'
    else:
        organism='human'
    
    print "Calculating the internal priming event probabilities... Please be patient! This step is time consuming!"
    params = ['Rscript', '--vanilla', '%s' %cleanUpdTSeq_R, '%s' %mergedBed, '%s' %NB_output, organism]
    run_command(params,shell=True)
    
    return


shortDescText="Cluster peak locations."
longDescText=""" Cluster peak locations. """
subcommands['clusterPeaks'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def clusterPeaks(args):
    """
    Cluster reads & calculate peak stats & combine with NB probs.
    This step runs an external bash script in ./src/ called "cluster.sh".
    #Input: 
    #Arguments: 
    #User options: 
    #Output(s):
    """
    #TODO: 1.Provide statistics and plots
    #      2.Needs bedtools exported to environment
    #      3.Provide options to user for choosing annotation. Create annotation for mouse
    #      4.You can annotate peaks with annotate() function
    argparser.add_argument("-bd", "--bed", required=True, help='A file stores candidate locations in bed format is required. Please provide with full directory.')
    argparser.add_argument("-cd", "--clust_dist", action='store', type=int, required=False, dest="cd", help="Clustering distance for adjacent intervals to determine peak locations, integer (Default is 10nt).")
    argparser.add_argument("-np", "--NBprobs", required=True, help='A file stores candidate locations with their probabilities being true event is required. Please provide with full directory.')
    argparser.add_argument("--mouse", action='store_true', required=False, help='Choose target organism genome as mouse to map read sequences (mm10). Please add this argument if your organism is mouse (Default is human).')
    args=argparser.parse_args(args=args)

    mergedBed = glob.glob(args.bed)[0]
    outDir=os.path.dirname(mergedBed)
    basename=mergedBed.strip().split("/")[-1][:-4]
    NB_output = outDir+"/"+basename+'.nbpred'

    if args.cd:
        cd=args.cd
    else:
        cd=20
    print cluster_sh
    print mergedBed
    print NB_output
    params = ['%s' %cluster_sh, '%s'%cd, '%s' %mergedBed, '%s' %NB_output]
    run_command(params,shell=True)
    
    return

def calculateLibsize(bedfile2librarysize):
	libSize=0
	for line in bedfile2librarysize:
		libSize = libSize + int(line.strip().split("\t")[4])
	return float(libSize)/1000000


shortDescText="Creates a bedgraph file for a given bed file."
longDescText=""" Creates a bedgraph file for a given bed file. """
subcommands['makeBedgraph'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def makeBedgraph(args):
    """
    #Input: 
    #Arguments: 
    #User options: 
    #Output(s): bed
    """
    #TODO: 
    argparser.add_argument("-bd", "--bed", required=True, help='A file stores candidate locations in bed format is required. Expression in raw read count. Please provide with full directory.')
    argparser.add_argument("--mouse", action='store_true', required=False, help='Choose target organism genome as mouse to map read sequences (mm10). Please add this argument if your organism is mouse (Default is human).')
    args=argparser.parse_args(args=args)
    file_bed = glob.glob(args.bed)[0]
    outDir=os.path.dirname(file_bed)
    basename=file_bed.strip().split("/")[-1][:-4]
    BedGraph_output = outDir+"/"+basename+'.bedgraph'
    bed = open('%s' %file_bed)
    libSize = calculateLibsize(bed)
    print "The library size of this file is :",libSize

    bed = open('%s' %file_bed)
    file_correctbedgraph = open('%s' %BedGraph_output, "w")

    file_correctbedgraph.write("browser position chr8:128748315-128753680"+"\n"+"browser hide all"+"\n"+"browser pack refGene encodeRegions"+"\n"+"browser full altGraph"+"\n"+"track type=bedGraph name="'%s' %basename+" description='BedGraph format' visibility=full color=200,100,0 altColor=0,100,200 priority=20"+"\n"+"#chrom"+"\t"+"start"+"\t"+"end"+"\t"+'%s'%basename+"\n")
    for line in bed:
	normCount = 0
    	if line.startswith("#"):
        	pass
        elif line.strip().split("\t")[5] == "-":
		normCount = float(int(line.strip().split("\t")[4]))/float(libSize)
        	file_correctbedgraph.write(line.strip().split("\t")[0]+"\t"+line.strip().split("\t")[1]+"\t"+line.strip().split("\t")[2]+"\t"+"-"+str(normCount)+"\n")
##        	file_correctbedgraph.write(line.strip().split("\t")[0]+"\t"+line.strip().split("\t")[1]+"\t"+line.strip().split("\t")[2]+"\t"+"-"+str(line.strip().split("\t")[4])+"\n")
        else:
                normCount = float(int(line.strip().split("\t")[4]))/float(libSize)
                file_correctbedgraph.write(line.strip().split("\t")[0]+"\t"+line.strip().split("\t")[1]+"\t"+line.strip().split("\t")[2]+"\t"+str(normCount)+"\n")
##                file_correctbedgraph.write(line.strip().split("\t")[0]+"\t"+line.strip().split("\t")[1]+"\t"+line.strip().split("\t")[2]+"\t"+str(line.strip().split("\t")[4])+"\n")
    print "Bedgraph file has been created!"
    bed.close()
    file_correctbedgraph.close()


    return

shortDescText="Combines all bed files of called peaks."
longDescText=""" Combines all bed files and returns a table with peaks. """
subcommands['makeSamplesTable'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def makeSamplesTable(args):
    """
    This step runs an external bash script in ./src/ called "make_samples_matrix.sh".
    #Input: 
    #Arguments: 
    #User options: 
    #Output(s):
    """
    #TODO: 1. Provide statistics and plots
    argparser.add_argument("-sf", "--Sfile", required=True, help='A file stores samples in fastq with full directories is required. Please provide with full directory.')
    argparser.add_argument("-g", "--gtf", required=True, help='This file is the gene annotation file in GTF format from which gene boundaries will be extracted. Please provide with a full directory.')
    args=argparser.parse_args(args=args)

    SamplesFile = glob.glob(args.Sfile)[0]
    GTFFile = glob.glob(args.gtf)[0]

    params = ['%s' %makeSamplesTable_sh, '%s' %SamplesFile, '%s' %GTFFile]
    run_command(params,shell=True)
    
    return


shortDescText="Performs PAS switch test using called peaks."
longDescText=""" Performs PAS switch test using called peaks. """
subcommands['switchTest'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def switchTest(args):
    """
    This step runs an external bash script in ./utils/ called "PASseq_chisquire_beta_v3.R".
    #Input: 
    #Arguments: 
    #User options: 
    #Output(s):
    """
    #TODO: 1. Provide statistics and plots
    argparser.add_argument("-sf", "--Sfile", required=True, help='A file stores samples in fastq with full directories is required. Please provide with full directory.')
    argparser.add_argument("-g", "--gtf", required=True, help='This file is the gene annotation file in GTF format from which gene boundaries will be extracted. Please provide with a full directory.')
    args=argparser.parse_args(args=args)

    SamplesFile = glob.glob(args.Sfile)[0]
    GTFFile = glob.glob(args.gtf)[0]

    params = ['%s' %SwitchTest_R_script, '%s' %SamplesFile, '%s' %GTFFile]
    run_command(params,shell=True)
    
    return






shortDescText="Annotates peaks with genomic locations such as 3'UTR."
longDescText=""" Annotates peaks with genomic locations such as 3'UTR and returns only annotated intervals."""
subcommands['annotatePeaks'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }


def annotatePeaks(args):
    """
    This step runs an external bash script in ./src/ called "annotate.sh".
    #Input: 
    #Arguments: 
    #User options: 
    #Output(s):
    """
    #TODO: 1. Provide statistics and plots
    argparser.add_argument("--region", default='3UTR', choices=['5UTR', 'Intronic', 'CodingExon','3UTR', '3UTR+1Kb_down'], required=False, help='Please choose a genomic region to annotate the peak locations with (Default is 3UTR for PASseq).')
    argparser.add_argument("-bd", "--bed", required=True, help='A file stores True Peaks locations in bed format to annotate is required. Please provide with full directory.')
    argparser.add_argument("-dd", "--down_dist", action='store', type=int, required=False, dest="downDist", help="Additional distance range downstream of 3'UTR regions to annotate peaks, integer (Default is 1000, 1Kb).")
    argparser.add_argument("--mouse", action='store_true', required=False, help='Choose target organism genome as mouse to map read sequences (mm10). Please add this argument if your organism is mouse (Default is human).')
    args=argparser.parse_args(args=args)

    truePeaksBed = glob.glob(args.bed)[0]
    outDir=os.path.dirname(truePeaksBed)
    
    basename=truePeaksBed.strip().split("/")[-1][:-4]
    
    if args.downDist:
        downDist=args.downDist
    else:
        downDist=1000
    print downDist
    if args.mouse:
        organism='mouse'
    else:
        organism='human'
    
    params = ['%s' %annotate_sh, '%s' %truePeaksBed, '%s' %organism, '%s' %downDist]
    run_command(params,shell=True) 
    
    return

#################################
### Expression quantification ###
#################################

class inputfiles:
    'Common base class for input files'
    file=0
    def __init__(self,name):
        self.name=name
        inputfiles.name += 1
    
class outputfiles:
    'Common base class for output files'
    file=0
    def __init__(self,name):
        self.name=name
        outputfiles.name += 1
        
####
install_dir='/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools'

cleanUpdTSeq_R = '%s/src/cleanUpdTSeq.R' %install_dir
cluster_sh = '%s/src/cluster.sh' %install_dir
annotate_sh = '%s/src/annotate.sh' %install_dir
makeSamplesTable_sh = '%s/src/make_samples_matrix.sh' %install_dir
SwitchTest_R_script = '%s/utils/PASseq_chisquire_beta_v3.R' %install_dir
####



if __name__ == "__main__":
    if len (sys.argv)==1:
        sys.argv.append("--help")  #if no command then it is a cry for help    
    main(sys.argv[1:])
