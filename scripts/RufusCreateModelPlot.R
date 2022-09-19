#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	          stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
	            print(args[1])
}

library(ggplot2)
data <- read.table(args[1],header=TRUE,  skip=5)


myCon = file(description = args[1], open="r", blocking = TRUE)
min = as.numeric(readLines(myCon, n = 1)) # Read one line from the connection.
cutoff = as.numeric(readLines(myCon, n = 1) )
genomesize = as.numeric(readLines(myCon, n = 1))
diploid = as.numeric(readLines(myCon, n = 1))
close(myCon)
#cutoff <- as.numeric(args[2])
#diploid <-as.numeric(args[3])
haploid <-diploid/2
trans=0.8
print(paste0(trimws(args[1]),".pdf"))
pdf(paste0(trimws(args[1]),".pdf"), width=6, height =2)
ggplot(data=data)+
  geom_line(aes(x=K, y=RawCount), size=1, color = "black", alpha=trans)+
  geom_line(aes(x=K, y=ModelSum), size=1, color="red", alpha=trans)+
  geom_line(aes(x=K, y=ErrorModel), color="yellow", alpha=trans)+
  geom_line(aes(x=K, y=X1x), color="green", alpha=trans)+
  geom_line(aes(x=K, y=X2x), color="#00008B", alpha=trans)+
  geom_line(aes(x=K, y=X3x), color="#0000CD", alpha=trans)+
  geom_line(aes(x=K, y=X4x), color="#0000FF", alpha=trans)+
  geom_line(aes(x=K, y=X5x), color="#4169E1", alpha=trans)+
  scale_y_log10(limits = c(1, max(data$RawCount)))+
  scale_x_continuous(limits=c(2,max(data$K)), breaks=seq(0,max(data$K), max(data$K)/10))+ 
  geom_vline(xintercept=cutoff, color="red", alpha=0.5)+
  geom_vline(xintercept=haploid, color="green", alpha=0.5)+
  geom_vline(xintercept=diploid, color="blue", alpha=0.5)+
  xlab("Kmer depth")+
  ylab("Frequency")
dev.off()
  
