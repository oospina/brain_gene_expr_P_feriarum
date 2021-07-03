x=read.table("ModuleSummaryExpanded_Geo.txt",sep="\t",head=T)

mods=unique(x$Module)
nMods=length(mods)
sigs=rep(F,nMods)
for(i in 1:nMods){sigs[i]="1" %in% x$isSig[x$Module==mods[i]]}

uMods=unique(mods)

print("     Number of synaptic pathways per module:",quote=F)
nSyn=rep(0,nMods)
sigMark=c("","*")
for(i in 1:nMods){
	nSyn[i]=sum(x$isSyn[x$Module==mods[i]])
	print(paste("          ",mods[i],"=",nSyn[i],sigMark[sigs[i]+1],sep=""),quote=F)
}
print(" ",quote=F)
print(" ",quote=F)


print("QUESTION STATED IN PAPER: Are there more synaptic genes in the significant modules than in the non-significant modules?")


print("INTENDED QUESTION: Are the synaptic pathways concentrated in the significant modules?",quote=F)
#compute number of synaptic pathways for each module
testStat=sum(x$isSig & x$isSyn)
print(paste("     TestStat = sum(synap pathways in signif modules) = ",testStat,sep=""),quote=F)

#generate the null distribution
nullDist=10000
for(j in 1:10000){
	randIsSyn=sample(x$isSyn)
	nullDist[j]=sum(x$isSig & randIsSyn)
}
#compute p-value
pval=mean(nullDist>=testStat)
print(paste("Total number of synaptic pathways expected in significant modules = ",mean(nullDist)))
print(paste("     p-value = ",pval,sep=""),quote=F)
print(" ",quote=F)
print(" ",quote=F)



###############################
print("QUESTION ANSWERED BY REVIEWERS METHOD: Is the number of synaptic modules per module higher for the significant modules compared to the nonsignificant modules?")
#compute number of synaptic pathways for each module
print("     Number of synaptic pathways per module:",quote=F)
nSyn=rep(0,nMods)
for(i in 1:nMods){
	nSyn[i]=sum(x$isSyn[x$Module==mods[i]])
}
#compute test statistic, avgNSynInSigModules - avgNSigInNonSigModules
testStat=mean(nSyn[sigs])-mean(nSyn[!sigs])
print(paste("     TestStat = mean(synap pathways in signif modules) - mean(synap pathways in nonsignif modules) = ",testStat,sep=""),quote=F)

#generate the null distribution
nullDist=10000
for(j in 1:10000){
	#recompute test stat after selecting a randomly chosen four modules
	chosen=rep(FALSE,nMods)
	chosen[sample(1:nMods,4)]=TRUE
	nullDist[j]=mean(nSyn[chosen])-mean(nSyn[!chosen])
}
#compute p-value
pval=mean(nullDist>=testStat)
print(paste("     p-value = ",pval,sep=""),quote=F)
print(" ",quote=F)
print(" ",quote=F)


