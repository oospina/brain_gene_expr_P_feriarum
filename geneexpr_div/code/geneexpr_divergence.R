library(ape)
library(phytools)
library(castor)

#get candidate list
cands=scan("CandidateTranscriptsFinal.csv",what=character())

#get key containing indivdual IDs and Group membership
key=read.csv("IndGroups.csv",head=F)

#get DE list
#genesDE_sex=read.table("topTags_FVsM.txt",sep="\t",head=T)$genes
#genesDE_sexA=read.table("topTags_ALLOPATRY_FVsM.txt",sep="\t",head=T)$genes
#genesDE_sexS=read.table("topTags_SYMPATRY_FVsM.txt",sep="\t",head=T)$genes
#genesDE_geo=read.table("topTags_AlloVsSym.txt",sep="\t",head=T)$genes
#genesDE_geoF=read.table("topTags_FEMALES_SymVsAllo.txt",sep="\t",head=T)$genes
#genesDE_geoM=read.table("topTags_MALES_SymVsAllo.txt",sep="\t",head=T)$genes
#genesDE=unique(c(genesDE_sex,genesDE_sexA,genesDE_sexS,genesDE_geo,genesDE_geoF,genesDE_geoM))

#get 7 Candidates that are differentially expressed
#candDEExpr=read.csv("CPM_normalized_counts_7CandDE.csv",head=F)
#candDE=candDEExpr[,1]


#GET NORMALIZED COUNTS
expr=read.csv("CPM_normalized_counts.csv",head=T)
geneNames=expr[,1]
nGenes=length(geneNames)

#remove column corresponding to gene name
#remove poor quality sample

expr=expr[,c(2:9,11:18)]

#standardize each gene's expression level to mean=0 and sd=1
expr2=expr-apply(expr,1,mean)
expr3=expr2/apply(expr2,1,sd)

expr=expr3

#Identify non-candidate genes
set1="NonCandidate"
isSet1=!is.element(geneNames,cands)
indexesSet1=(1:nGenes)[isSet1]
nSet1=length(indexesSet1)

#Identify candidate genes
set2="CandidateGenes"
isSet2=is.element(geneNames,cands)
indexesSet2=(1:nGenes)[isSet2]
nSet2=length(indexesSet2)


#Obtain the sample group information
indNames=colnames(expr)
groupNames=c("AM","AF","SM","SF")
samples=read.csv("samples.csv")

#combine inds within each of four groups by averaging
idsAM=samples$library[samples$situation=="Allo" & samples$sex=="M"]
idsAF=samples$library[samples$situation=="Allo" & samples$sex=="F"]
idsSM=samples$library[samples$situation=="Sym" & samples$sex=="M"]
idsSF=samples$library[samples$situation=="Sym" & samples$sex=="F"]

isAM=is.element(indNames,idsAM)
isAF=is.element(indNames,idsAF)
isSM=is.element(indNames,idsSM)
isSF=is.element(indNames,idsSF)

exprAM=apply(expr[,isAM],1,mean)
exprAF=apply(expr[,isAF],1,mean)
exprSM=apply(expr[,isSM],1,mean)
exprSF=apply(expr[,isSF],1,mean)

exprAvg=cbind(AM=exprAM,AF=exprAF,SM=exprSM,SF=exprSF)


#construct trees from groups
trees1=nj(dist(t(exprAvg[isSet1,])))
trees2=nj(dist(t(exprAvg[isSet2,])))
brls1=trees1$edge.length[order(trees1$edge[,2])]
brls2=trees2$edge.length[order(trees2$edge[,2])]


#plot expression trees
pdf("ExpressionTrees_Groups.pdf")
plot(trees1,type="unrooted",main=set1)
plot(trees2,type="unrooted",main=set2)
dev.off()

#compute ratios
ratios=matrix(0,1,8)
colnames(ratios)=c("AFvAM1","SMvSF1","SMvAM1","AFvSF1","AFvAM2","SMvSF1","SMvAM2","AFvSF1")

ratios[1]=log2(brls1[2]/brls1[1]) #AFvAM1
ratios[2]=log2(brls1[3]/brls1[4]) #SMvSF1
ratios[3]=log2(brls1[3]/brls1[1]) #SMvAM1
ratios[4]=log2(brls1[2]/brls1[4]) #AFvSF1
ratios[5]=log2(brls2[2]/brls2[1]) #AFvAM2
ratios[6]=log2(brls2[3]/brls2[4]) #SMvSF1
ratios[7]=log2(brls2[3]/brls2[1]) #SMvAM2
ratios[8]=log2(brls2[2]/brls2[4]) #AFvSF1

#generate bootstrap distributions
ratiosB=matrix(0,1000,8)
colnames(ratiosB)=colnames(ratios)

for(i in 1:1000){
	bootChosen1=sample(indexesSet1,nSet1,replace=T)
	bootChosen2=sample(indexesSet2,nSet2,replace=T)
	t1=nj(dist(t(exprAvg[bootChosen1,])))
	t2=nj(dist(t(exprAvg[bootChosen2,])))
	b1=t1$edge.length[order(t1$edge[,2])]
	b2=t2$edge.length[order(t2$edge[,2])]
	ratiosB[i,1]=log2(b1[2]/b1[1]) #AFvAM1
	ratiosB[i,2]=log2(b1[3]/b1[4]) #SMvSF1
	ratiosB[i,3]=log2(b1[3]/b1[1]) #SMvAM1
	ratiosB[i,4]=log2(b1[2]/b1[4]) #AFvSF1
	ratiosB[i,5]=log2(b2[2]/b2[1]) #AFvAM2
	ratiosB[i,6]=log2(b2[3]/b2[4]) #SMvSF1
	ratiosB[i,7]=log2(b2[3]/b2[1]) #SMvAM2
	ratiosB[i,8]=log2(b2[2]/b2[4]) #AFvSF1	
}

#determine 95% bootstrap confidence intervals
boot95=matrix(0,2,8)
colnames(boot95)=colnames(ratiosB)
rownames(boot95)=c("2.5%","9.75%")
for(i in 1:8){
	boot95[,i]=as.numeric(quantile(ratiosB[,i],probs=c(0.025,0.975)))
}
pdf("RatiosPlot.pdf")
plot(1:8,ratios,ylim=range(boot95))
text(1:8,ratios,labels=colnames(boot95))
for(i in 1:8){
	lines(c(i,i),boot95[,i])
}
dev.off()

#SIGNIFICANCE TESTS
#construct tree for each gene separately
trees=list()
brls=matrix(0,nGenes,5)	#first four coorespond to groups in order of groupNames, last brl is central branch
colnames(brls)=c(groupNames,"central")
for(i in 1:nGenes){
	trees[[i]]=nj(dist(exprAvg[i,]))
	brls[i,]=(trees[[i]]$edge.length[order(trees[[i]]$edge[,2])])
}

diffs=matrix(0,nGenes,4)
colnames(diffs)=c("AFvAM","SMvSF","SMvAM","AFvSF")
diffs[,1]=brls[,2]-brls[,1] #AFvAM
diffs[,2]=brls[,3]-brls[,4] #SMvSF
diffs[,3]=brls[,3]-brls[,1] #SMvAM
diffs[,4]=brls[,2]-brls[,4] #AFvSF

#conduct randomization tests to see if each of the 8 comparisions is significantly differnt than zero
#test stat is mean branch length difference
testStats=rep(0,8) #corresponds to 8 comparisons of ratios
testStats[1]=mean(diffs[indexesSet1,1])
testStats[2]=mean(diffs[indexesSet1,2])
testStats[3]=mean(diffs[indexesSet1,3])
testStats[4]=mean(diffs[indexesSet1,4])
testStats[5]=mean(diffs[indexesSet2,1])
testStats[6]=mean(diffs[indexesSet2,2])
testStats[7]=mean(diffs[indexesSet2,3])
testStats[8]=mean(diffs[indexesSet2,4])

nullDists=matrix(0,1000,8)
colnames(nullDists)=colnames(ratios)
for(i in 1:1000){
	tempIndexes=sample(c(indexesSet1,indexesSet2))
	tempIndexes1=tempIndexes[1:nSet1]
	tempIndexes2=tempIndexes[(nSet1+1):(nSet1+nSet2)]
	nullDists[i,1]=mean(diffs[tempIndexes1,1])
	nullDists[i,2]=mean(diffs[tempIndexes1,2])
	nullDists[i,3]=mean(diffs[tempIndexes1,3])
	nullDists[i,4]=mean(diffs[tempIndexes1,4])
	nullDists[i,5]=mean(diffs[tempIndexes2,1])
	nullDists[i,6]=mean(diffs[tempIndexes2,2])
	nullDists[i,7]=mean(diffs[tempIndexes2,3])
	nullDists[i,8]=mean(diffs[tempIndexes2,4])	
}
pvalues=rep(0,8)
for(i in 1:8){
	pvalues[i]=mean(nullDists[,i]>=testStats[i])
}
testResultIsZero=cbind(TestStat=testStats,PValue=pvalues)
rownames(testResultIsZero)=colnames(ratios)

#conduct randomization teset to see if branch length differences between groups are more extreme in sympatry than allopatry
#test stat is difference of mean branch length differences (cand - noncand)
testStats2=rep(0,4) #corresponds to 8 comparisons of ratios
testStats2[1]=mean(diffs[indexesSet2,1])-mean(diffs[indexesSet1,1])
testStats2[2]=mean(diffs[indexesSet2,2])-mean(diffs[indexesSet1,2])
testStats2[3]=mean(diffs[indexesSet2,3])-mean(diffs[indexesSet1,3])
testStats2[4]=mean(diffs[indexesSet2,4])-mean(diffs[indexesSet1,4])

nullDists2=matrix(0,1000,4)
for(i in 1:1000){
	tempIndexes=sample(c(indexesSet1,indexesSet2))
	tempIndexes1=tempIndexes[1:nSet1]
	tempIndexes2=tempIndexes[(nSet1+1):(nSet1+nSet2)]
	nullDists2[i,1]=mean(diffs[tempIndexes2,1])-mean(diffs[tempIndexes1,1])
	nullDists2[i,2]=mean(diffs[tempIndexes2,2])-mean(diffs[tempIndexes1,2])
	nullDists2[i,3]=mean(diffs[tempIndexes2,3])-mean(diffs[tempIndexes1,3])
	nullDists2[i,4]=mean(diffs[tempIndexes2,4])-mean(diffs[tempIndexes1,4])
}
pvalues2=rep(0,4)
for(i in 1:4){
	pvalues2[i]=mean(nullDists2[,i]>=testStats2[i])
}
testResultIsEnhanced=cbind(TestStat=testStats2,PValue=pvalues2)
rownames(testResultIsEnhanced)=c("AFvAM","SMvSF","SMvAM","AFvSF")

