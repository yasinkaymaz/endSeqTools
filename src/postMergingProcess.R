
library(stats)
args = commandArgs (T)
InputMeatMapData = args[1]
Outfile = args[2]
NoS= as.numeric(args[3])

SumCol=NoS+1
SumCol=NoS+2

data <- read.table(InputMeatMapData)

data$VSumCol <- rowSums(data[,2:NoS])
data$VMeanCol <- rowMeans(data[,2:NoS])


write.table(data[which(data$VSumCol > 100 ),1:NoS], file=Outfile, quote=FALSE, sep="\t",row.names=FALSE,col.names=FALSE)



