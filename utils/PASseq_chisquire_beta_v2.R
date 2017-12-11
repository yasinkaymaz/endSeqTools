#Chi-Square test for differential PolyA usage
options("scipen"=100, "digits"=5)

library(cleanUpdTSeq)
args = commandArgs (T)
InputBedFile = args[1]
Outfile = args[2]
NoC1=as.numeric(args[3])
NoC2=as.numeric(args[4])
NoS=NoC1+NoC2

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
#for (i in c("ACO2")){ 
  nPA=genelist[i,]$Freq
  GeneTotalC1 <- sum(data[which(data[,NoS+8] == i  ),8:(7+NoC1)])
  GeneTotalC2 <- sum(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)])
  PAsC1 <- rowSums(data[which(data[,NoS+8] == i  ),8:(7+NoC1)])
  PAsC2 <- rowSums(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)])
  
  for (j in 1:nPA){
    x <- round(matrix(c(PAsC1[j]/NoC1, PAsC2[j]/NoC2, (GeneTotalC1-PAsC1[j])/NoC1, (GeneTotalC2-PAsC2[j])/NoC2), byrow = TRUE, 2, 2))
    print(x)
    cqtest <- chisq.test(x)
    tested.PAsites <- as.data.frame(data[which(data[,NoS+8] == i  ),c(1,2,3,4,5,6,(8+NoS),(9+NoS),(10+NoS))][j,])
    tested.PAsites$Num_Total_PA_sites_of_Gene <- nPA
    tested.PAsites$PA_MeanCond1 <- PAsC1[j]/NoC1
    tested.PAsites$PA_MeanCond2 <- PAsC2[j]/NoC2
    tested.PAsites$Rest_MeanCond1 <- (GeneTotalC1-PAsC1[j])/NoC1
    tested.PAsites$Rest_MeanCond2 <- (GeneTotalC2-PAsC2[j])/NoC2
    tested.PAsites$Chi2.statistics  <- cqtest$statistic
    tested.PAsites$Chi2.pvalue <- cqtest$p.val
    tested.PAsites$Fis.pvalue <- fisher.test(x)$p.value
    
    tested.Genes <- rbind(tested.Genes, tested.PAsites)
    
  }
  
}


tested.Genes$Chi.padj <- p.adjust(tested.Genes$Chi2.pvalue, n=length(tested.Genes[,1]), method="BH")
tested.Genes$Fis.padj <- p.adjust(tested.Genes$Fis.pvalue, n=length(tested.Genes[,1]), method="BH")

#write.table(tested.Genes[order(tested.Genes$Chi2.padj),], file="switchtest.txt")
#write.table(tested.Genes[which(tested.Genes$Chi2.padj <0.05),], file="Outfile.switch.test.txt",sep="\t")
write.table(tested.Genes[order(tested.Genes$Chi.padj),], file=Outfile,sep="\t")




