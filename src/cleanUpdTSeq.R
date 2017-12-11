# Rscript for cleanUpdTseq Naive Bayes classifier
library(cleanUpdTSeq)
args = commandArgs (T)
InputBedFile = args[1]
NBPredictions = args[2]
organism = args[3]

testSet <- read.table(InputBedFile, header=FALSE)
peaks <- BED2GRangesSeq(testSet, withSeq=FALSE)

if(organism == "mouse"){
library(BSgenome.Mmusculus.UCSC.mm10)
testSet.NaiveBayes <- buildFeatureVector(peaks, BSgenomeName=Mmusculus,
                                          upstream=40, downstream=30,
                                          wordSize=6, alphabet=c("ACGT"),
                                          sampleType="unknown",
                                          replaceNAdistance=30,
                                          method="NaiveBayes",
                                          ZeroBasedIndex=1, fetchSeq=TRUE)

}else{
library("BSgenome.Hsapiens.UCSC.hg19")
testSet.NaiveBayes <- buildFeatureVector(peaks, BSgenomeName=Hsapiens,
                                          upstream=40, downstream=30,
                                          wordSize=6, alphabet=c("ACGT"),
                                          sampleType="unknown",
                                          replaceNAdistance=30,
                                          method="NaiveBayes",
                                          ZeroBasedIndex=1, fetchSeq=TRUE)

}

data(data.NaiveBayes)

predictTestSet(data.NaiveBayes$Negative, 
	data.NaiveBayes$Positive, 
	testSet.NaiveBayes=testSet.NaiveBayes, 
	outputFile=NBPredictions, 
	assignmentCutoff=0.5)


