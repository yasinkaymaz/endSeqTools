#Chi-Square test for differential PolyA usage
options("scipen"=100, "digits"=5)

library(cleanUpdTSeq)
args = commandArgs (T)
InputBedFile = args[1]
Outfile = args[2]
NoC1=as.numeric(args[3])
NoC2=as.numeric(args[4])
#temp;

NoS=NoC1+NoC2
#NoS=15
#NoC1=6
#NoC2=9

#setwd("/home/yasinkaymaz/Dropbox/PASseq")
#data <- read.delim("mergedPeaks_annotated_CPM_subset.bed",header=FALSE)
#data <- read.delim("mergedPeaks_annotated_CPM.bed",header=FALSE)
data <- read.delim(InputBedFile,header=FALSE)
head(data)
data <- data[rowSums(data[,8:(7+NoS)]) > 0,  ]
#Filter low count PA sites!

genes <- data[,NoS+8]
gene.table <- as.data.frame(table(genes))
genelist <- noquote(droplevels(gene.table[which(gene.table$Freq >1),]))
rownames(genelist) <- genelist$genes
tested.Genes <- NULL
for (i in genelist$genes){ 
  nPA=genelist[i,]$Freq
#  print(paste("Gene",i,"has",nPA,"PolyA sites",sep = " "))
  #print(data[which(data[,NoS+8] == i  ),8:(7+NoC1)])
  #print(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)])
  GeneTotalC1 <- sum(data[which(data[,NoS+8] == i  ),8:(7+NoC1)])
  GeneTotalC2 <- sum(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)])
#  print(paste("Total normalized count of gene",i,"in condition 1 is",GeneTotalC1))
#  print(paste("Total normalized count of gene",i,"in condition 2 is",GeneTotalC2))
  #print(rowSums(data[which(data[,NoS+8] == i  ),8:(7+NoC1)]))
  #print(rowSums(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)]))
  PAsC1 <- rowSums(data[which(data[,NoS+8] == i  ),8:(7+NoC1)])
  PAsC2 <- rowSums(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)])
  
    for (j in 1:nPA){
#      print(paste("Rowsum of ",j,"th PA site of gene",i,"in condition 1 is", PAsC1[j]))
#      print(paste("Rowsum of ",j,"th PA site of gene",i,"in condition 2 is", PAsC2[j]))
      
      x <- round(matrix(c(PAsC1[j], PAsC2[j], GeneTotalC1-PAsC1[j], GeneTotalC2-PAsC2[j]), byrow = TRUE, 2, 2))
      
      #print(x)
      cqtest <- chisq.test(x)
#      print(cqtest$p.val)
      #pval.colum <- noquote(paste("V",(11+NoS),sep=""))
      #print(pval.colum)
      tested.PAsites <- as.data.frame(data[which(data[,NoS+8] == i  ),c(1,2,3,4,5,6,(8+NoS),(9+NoS),(10+NoS))][j,])
      tested.PAsites$MeanCond1 <- PAsC1[j]/NoC1
      tested.PAsites$MeanCond2 <- PAsC2[j]/NoC2
      tested.PAsites$Chi.pvalue <- cqtest$p.val
      tested.PAsites$Fis.pvalue <- fisher.test(x)$p.value
      
      tested.Genes <- rbind(tested.Genes, tested.PAsites)
      
    }
    
}
tested.Genes$Chi.padj <- p.adjust(tested.Genes$Chi.pvalue, n=length(tested.Genes[,1]), method="BH")
tested.Genes$Fis.padj <- p.adjust(tested.Genes$Fis.pvalue, n=length(tested.Genes[,1]), method="BH")


#write.table(tested.Genes[order(tested.Genes$Chi.padj),], file="switchtest.txt")
write.table(tested.Genes[order(tested.Genes$Chi.padj),], file=Outfile)





#for (i in 1:length(data[,1])){#for each of the line
#  if (data[i,10] >1) {#if the number of peaks are more than 1
#  tab <- matrix(as.integer(c(data[i,11:14])), ncol=2)#Store these columns in a matrix: Cond1PeakMean	Cond2PeakMean	Cond1RestofPeaks	Cond2RestofPeaks
#  data[i, 15] <- fisher.test(tab)$p.value #Do a fisher's test on these values and store the pvalue on col 15
#  }
#  else {#otherwise
#  data[i, 15] <- t.test(tab)$p.value# do a student's T-test and store the pvalue on col 15
#  #data[i, 15] <- c('NA')
#  }
#}

#data[,16] <- p.adjust(data$V15, n=length(data$V15), method="BH")#Then adjust the calculated pvalues on col16 for multiple hypothesis testing

#for (i in 1:length(data[,1])){#then for each line of the file
#  if(data[i,16] < 0.05){#if the calculated pvalue is less then alpha
#    data[i, 17] <- c('Significant')#assign a column wih "Significant"
#    data[i, 18] <- log((as.numeric(data[i,11])/as.numeric(data[i,12])), 2)#and a column with the log-base2 Fold ratio value to col 18
#  } else {#otherwise
#    data[i, 17] <- c('Not Sig')#assign 'Not Significant'
#    data[i, 18] <- c('NA')#and NA
#  }
#}

#write the output with additional columns (+4 columns)
#write.table(data, file="golub_All_Union_closerPeaksClustered_final.bed.annotated_NoRedundants.out2_comparePeaks2_significance", sep="\t")


#Alternative
#library("data.table")
#data.t <- data.table(data)
#data <- read.table("mergedPeaks_annotated.bed",header=FALSE)

