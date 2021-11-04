#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	          stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
	            print(args[1])
}

library(ggplot2)
data <- read.table(args[1],header=FALSE,  skip=6)
trans=0.8
print(paste0(trimws(args[1]),".pdf"))
pdf(paste0(trimws(args[1]),".pdf"), width=6, height =2)
ggplot(data=data)+
  geom_line(aes(x=V1, y=V4), size=1, color = "black", alpha=trans)+
  geom_line(aes(x=V1, y=V5), size=1, color="red", alpha=trans)+
  geom_line(aes(x=V1, y=V6), color="green", alpha=trans)+
  geom_line(aes(x=V1, y=V7), color="#00008B", alpha=trans)+
  geom_line(aes(x=V1, y=V8), color="#0000CD", alpha=trans)+
  geom_line(aes(x=V1, y=V9), color="#0000FF", alpha=trans)+
  geom_line(aes(x=V1, y=V10), color="#4169E1", alpha=trans)+
  xlim(0,500)+
  xlab("Kmer depth")+
  ylab("Frequency")
dev.off()
  
