x=read.table("LRT_AlloVsSym.txt",sep="\t",head=T)
pdf("PvalueHist_Geo.pdf")
hist(x$PValue,nclass=50,ylim=c(0,3000),main="AlloVsSym",xlab="PValue",ylab="Number of transcripts")
dev.off()

x=read.table("LRT_FVsM.txt",sep="\t",head=T)
pdf("PvalueHist_Sex.pdf")
hist(x$PValue,nclass=50,ylim=c(0,3000),main="FVsM",xlab="PValue",ylab="Number of transcripts")
dev.off()
