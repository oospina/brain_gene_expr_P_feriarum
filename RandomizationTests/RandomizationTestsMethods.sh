
########
print("Are there more genes differentially expressed by Geography than by Sex?",quote=F)
geoOnly=rep(c(1,0),241)
sexOnly=rep(c(0,1),43)
GS_Both=rep(c(1,1),14)
neither=rep(c(0,0),47956)

all=t(matrix(c(geoOnly,sexOnly,GS_Both,neither),nrow=2))
nRow=nrow(all)

nGeo=sum(all[,1]==1 & all[,2]==0)+sum(all[,1] & all[,2]==1)
nSex=sum(all[,1]==0 & all[,2]==1)+sum(all[,1] & all[,2]==1)

testStat=nGeo/nSex
print(paste("Test stat (nGeo/nSex)=",testStat,sep=""),quote=F)
nullDist=rep(0,10000)

for(i in 1:10000){
	#randomize values across cells within a row (mix up geo and sex)
	switch=sample(c(T,F),nRow,replace=T)
	rAll=all
	rAll[switch,1]=all[switch,2]
	rAll[switch,2]=all[switch,1]

	rnGeo=sum(rAll[,1]==1 & rAll[,2]==0)+sum(rAll[,1] & rAll[,2]==1)
	rnSex=sum(rAll[,1]==0 & rAll[,2]==1)+sum(rAll[,1] & rAll[,2]==1)

	nullDist[i]=rnGeo/rnSex


}
pvalue=mean(nullDist>=testStat)
if(pvalue==0){
	print("p-value<0.0001")
}else{
	print(paste("p-value=",pvalue,sep=""),quote=F)
}

########
print("Are there more gene differentially expressed by Geography in females than in males?",quote=F)
#Values from Figure 4
geoF=rep(c(1,0),170+4+3)
geoM=rep(c(0,1),52+3+4)
geoB=rep(c(1,1),19)
neither=rep(c(0,0),48254-(170+4+3+52+3+4+19))

all=t(matrix(c(geoF,geoM,geoB,neither),nrow=2))
nRow=nrow(all)

nGeoF=sum(all[,1]==1 & all[,2]==0)+sum(all[,1] & all[,2]==1)
nGeoM=sum(all[,1]==0 & all[,2]==1)+sum(all[,1] & all[,2]==1)

testStat=nGeoF/nGeoM
print(paste("Test stat (nGeoF/nGeoM)=",testStat,sep=""),quote=F)
nullDist=rep(0,10000)

for(i in 1:10000){
	#randomize values across cells within a row (mix up geo and sex)
	switch=sample(c(T,F),nRow,replace=T)
	rAll=all
	rAll[switch,1]=all[switch,2]
	rAll[switch,2]=all[switch,1]

	rnGeoF=sum(rAll[,1]==1 & rAll[,2]==0)+sum(rAll[,1] & rAll[,2]==1)
	rnGeoM=sum(rAll[,1]==0 & rAll[,2]==1)+sum(rAll[,1] & rAll[,2]==1)

	nullDist[i]=rnGeoF/rnGeoM
}
pvalue=mean(nullDist>=testStat)
if(pvalue==0){
	print("p-value<0.0001")
}else{
	print(paste("p-value=",pvalue,sep=""),quote=F)
}


########
print("Are there more genes differentially expressed by Sex in sympatry than in allopatry?",quote=F)
#Values from Figure 4
sexS=rep(c(1,0),4+25+4)
sexA=rep(c(0,1),3+17+3)
sexB=rep(c(1,1),1)
neither=rep(c(0,0),48254-(4+25+4+3+17+3))

all=t(matrix(c(sexS,sexA,sexB,neither),nrow=2))
nRow=nrow(all)

nSexS=sum(all[,1]==1 & all[,2]==0)+sum(all[,1] & all[,2]==1)
nSexA=sum(all[,1]==0 & all[,2]==1)+sum(all[,1] & all[,2]==1)

testStat=nSexS/nSexA
print(paste("Test stat (nSexS/nSexA)=",testStat,sep=""),quote=F)
nullDist=rep(0,10000)

for(i in 1:10000){
	#randomize values across cells within a row (mix up geo and sex)
	switch=sample(c(T,F),nRow,replace=T)
	rAll=all
	rAll[switch,1]=all[switch,2]
	rAll[switch,2]=all[switch,1]

	rnSexS=sum(rAll[,1]==1 & rAll[,2]==0)+sum(rAll[,1] & rAll[,2]==1)
	rnSexA=sum(rAll[,1]==0 & rAll[,2]==1)+sum(rAll[,1] & rAll[,2]==1)

	nullDist[i]=rnSexS/rnSexA
}
pvalue=mean(nullDist>=testStat)
if(pvalue==0){
	print("p-value<0.0001")
}else{
	print(paste("p-value=",pvalue,sep=""),quote=F)
}

#MODULE TESTS
library(qvalue)
x=read.table("ModuleSummary.txt",sep="\t",head=T)
mods=x$module
nMods=length(mods)
nPaths=x$nPaths
nSynap=x$nSynap

#correct for multiple tests
pvalues=c(x$pvalGeo,x$pvalSex)
qvalues=qvalue(pvalues)$qvalues

#create expanded table for geo
isSignificant=qvalues[1:nMods]<0.05
Module=NULL
isSig=NULL
isSyn=NULL
for(i in 1:nMods){
	Module=c(Module,rep(mods[i],nPaths[i]))
	isSig=c(isSig,rep(isSignificant[i],nPaths[i])+0)
	isSyn=c(isSyn,rep(1,nSynap[i]),rep(0,nPaths[i]-nSynap[i])+0)
}
y=cbind(Module,isSig,isSyn)

write.table(y,file="ModuleSummaryExpanded_Geo.txt",quote=F,row.names=F,sep="\t")



x=read.table("ModuleSummaryExpanded_Geo.txt",sep="\t",head=T)

mods=unique(x$Module)
nMods=length(mods)
sigs=rep(F,nMods)
for(i in 1:nMods){sigs[i]="1" %in% x$isSig[x$Module==mods[i]]}

testStats=rep(-1,nMods)
pvalues2=rep(-1,nMods)

#############################################
print("Are synaptic pathways significantally concentrated in module X?",quote=F)
#module 1...
for(i in 1:nMods){
	print("------------")
	print(paste("Module ",mods[i],sep=""),quote=F)
	#compute test statistic, number of synaptic pathway in module
	testStat=sum(x$isSyn[x$Module==mods[i]])
	testStats[i]=testStat
	print(paste("     TestStat = ",testStat,sep=""),quote=F)
	#generate null distribution
	nullDist=10000
	for(j in 1:10000){
		#recompute test stat after randomizing synaptic pathway designation
		nullDist[j]=sum(sample(x$isSyn)[x$Module==mods[i]])
	}
	#compute p-value
	pval=mean(nullDist>=testStat)
	pvalues2[i]=pval
	print(paste("     p-value = ",pval,sep=""),quote=F)
	print("------------",quote=F)
}
print(" ",quote=F)
print(" ",quote=F)

library(qvalue)
qvalues2=qvalue(pvalues2)$qvalues

results=cbind(mods,pvalues,qvalues,testStats,pvalues2,qvalues2)
write.table(results,file="SynapConcModulesTestResultsA.txt",sep="\t",quote=F)

###############################
print("Does synap come up more in the significant modules than you would expect by chance?",quote=F)
#compute number of synaptic pathways for each module
print("     Number of synaptic pathways per module:",quote=F)
nSyn=rep(0,nMods)
for(i in 1:nMods){
	nSyn[i]=sum(x$isSyn[x$Module==mods[i]])
	print(paste("          ",mods[i],"=",nSyn[i],sep=""),quote=F)
}
#compute test statistic, avgNSynInSigModules - avgNSigInNonSigModules
testStat=mean(nSyn[sigs])-mean(nSyn[!sigs])
print(paste("     TestStat = ",testStat,sep=""),quote=F)

#generate the null distribution
nullDist=10000
for(j in 1:10000){
	#recompute test stat after randomizing synaptic pathway designation
	randIsSyn=sample(x$isSyn)
	#recompute nSyn for each module
	nSyn=rep(0,nMods)
	for(i in 1:nMods){
		nSyn[i]=sum(randIsSyn[x$Module==mods[i]])
	}
	nullDist[j]=mean(nSyn[sigs])-mean(nSyn[!sigs])
}
#compute p-value
pval=mean(nullDist>=testStat)
print(paste("     p-value = ",pval,sep=""),quote=F)
print(" ",quote=F)
print(" ",quote=F)


###############################
print("Are there more significant modules than expected by chance?",quote=F)
#GEO
#compute test statistics , number of significant modules
testStat=nSigByGeo
print(paste("     TestStatGeo = ",testStat,sep=""),quote=F)
#Generate null distribution
nullDist=rbinom(10000,22,0.05)
#compute pvalue
pval=mean(nullDist>=testStat)
print(paste("     p-value = ",pval,sep=""),quote=F)
print(" ",quote=F)
#SEX
#compute test statistics , number of significant modules
testStat=nSigBySex
print(paste("     TestStatSex = ",testStat,sep=""),quote=F)
#Generate null distribution
nullDist=rbinom(10000,22,0.05)
#compute pvalue
pval=mean(nullDist>=testStat)
print(paste("     p-value = ",pval,sep=""),quote=F)
print(" ",quote=F)
print(" ",quote=F)

