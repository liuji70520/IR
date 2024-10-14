args<-commandArgs(TRUE)
a <- args[1]
b <- args[2]

setwd(b)

data<-read.table(paste0(a,".txt"),sep = "\t")

tick <-rownames(data)

data<-t(data)

colnames(data)<-NULL

c <- ncol(data)

library(ggseqlogo)
library(ggplot2)

csl <- make_col_scheme(chars = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"), cols = c("#D62839","#D62839","#D62839","#D62839","#D62839","#F7B32B","#F7B32B","#F7B32B","#F7B32B","#F7B32B","#109648","#109648","#109648","#109648","#109648","#255C99","#255C99","#255C99","#255C99","#255C99"))

png(paste0(a,".png"),width=100*c,height=270)
ggseqlogo(data, method="prob",col_scheme=csl) + scale_x_discrete(limits=tick) + ylab("Allele Frequency") + xlab("Potential Mutation Site")
dev.off()
