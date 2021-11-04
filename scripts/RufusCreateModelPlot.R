#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	          stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
	            print(args[1])
}

library(ggplot2)
data <- read.table(args[1],header=TRUE,  skip=5)
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
  scale_y_log10(limits = c(10000, max(data$RawCount)))+
  xlim(2,max(data$K))+
  xlab("Kmer depth")+
  ylab("Frequency")
dev.off()
  
