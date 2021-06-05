genesDE_sex=read.table("topTags_FVsM.txt",sep="\t",head=T)$genes
genesDE_sexA=read.table("topTags_ALLOPATRY_FVsM.txt",sep="\t",head=T)$genes
genesDE_sexS=read.table("topTags_SYMPATRY_FVsM.txt",sep="\t",head=T)$genes
genesDE_geo=read.table("topTags_AlloVsSym.txt",sep="\t",head=T)$genes
genesDE_geoF=read.table("topTags_FEMALES_SymVsAllo.txt",sep="\t",head=T)$genes
genesDE_geoM=read.table("topTags_MALES_SymVsAllo.txt",sep="\t",head=T)$genes
genesDE=unique(c(genesDE_sex,genesDE_sexA,genesDE_sexS,genesDE_geo,genesDE_geoF,genesDE_geoM))


library(plotly)
x=read.csv("CPM_normalized_counts.csv")
#retain only DE genes
good=match(genesDE,x$genes)
good <- good[!is.na(good)]

x=x[good,]
#remove poor quality sample
x2=log2(x[,c(2:9,11:18)])

d=dist(t(x2))
colrs=c("red","red","orange","red","orange","red","orange","red","blue","cyan","blue","cyan","cyan","blue","blue","cyan")

mds=cmdscale(d,eig=FALSE,k=2)
pdf("MDS_AllDEGenes.pdf")
plot(mds,col=colrs,pch=19,cex=2,xlab="MDS Axis 1",ylab="MDS Axis 2")
text(mds[,1],mds[,2],labels=c(29357:29364,29366:29373))
dev.off()

mds=cmdscale(d,eig=TRUE,k=3)
df=as.data.frame(mds$points)
plot_ly(df,x=~V1,y=~V2,z=~V3,color=colrs)
