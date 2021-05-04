## FOR INSTALLATION OF WGCNA AND DEPENDENCIES:
# install.packages("fastcluster", "dynamicTreeCut", "robust")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#	install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("GO.db", "preprocessCore", "impute"))

# orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg")
# orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6))
# packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="")
# BiocManager::install(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))

# install.packages("~/Downloads/WGCNA_1.67.tgz")

###################

library("edgeR")
library("WGCNA")
library("miceadds")
library(gplots) # For network heatmap colors
library("stringr") # For replacing of DEgenes column in Cytoscape network

## Both of the following options seem to be important for WGCNA
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

setwd("~/Desktop/wgcna_RNASeq_TwoFourTreatsDEs/")

## Vector of EXCLUDED SAMPLES
exclude <- c("I29365")

## Read and prepare data as if it were to be analyzed by edgeR
fercounts <- read.csv("sumRawCounts.csv")
fercounts <- na.omit(fercounts)
rownames(fercounts) <- fercounts[, 1]
fercounts <- fercounts[, -1]
fercounts <- fercounts[, -grep(exclude, colnames(fercounts))]
sampleMeta <- read.csv("samples.csv", header = T)
sampleMeta <- sampleMeta[-grep(exclude, sampleMeta[ , 1]), ]

## Create DGEList object with factors and groups
rownames(sampleMeta) <- sampleMeta[, 1]
sample.ids <- apply(sampleMeta[, 1:3], 1, paste, collapse = "_")
indivs <- c()
for (i in 1:length(colnames(fercounts))) {
	loop <- grep(colnames(fercounts[i]), sample.ids, value = T)
	indivs[[length(indivs) + 1]] <- loop
}

situation <- gsub("^[[:alnum:]]*_", "", indivs)
situation <- gsub("_[[:alnum:]]*$", "", situation)
sex <- gsub("\\w*_", "", indivs)
Group <- factor(paste(situation, sex, sep = "."))
dgList <- DGEList(counts = fercounts, group = Group)

## Add "Group" (situation x sex) column to sampleMeta
sampleMeta <- cbind(sampleMeta, Group)

## perform removal of low (near zero) count genes and TMM normalization
keep <- rowMeans(cpm(dgList, log=F, prior.count=2)) >= 1 
dgList_NormTMM <- DGEList(dgList$counts[keep, ], group = Group, lib.size = colSums(dgList$counts[keep, ]))
dgList_NormTMM <- calcNormFactors(dgList_NormTMM, method = "TMM") ## Trimmed Mean M-values
dgList_NormTMM <- cpm(dgList_NormTMM, log=T, prior.count=2)

## Keep top % genes (eg. 25% for quantile 0.75) with highest variance
## Specifify quantile:
quant <- c(0.90)
dgList_NormTMMVarFil <- cbind(dgList_NormTMM, apply(dgList_NormTMM, 1, var))
colnames(dgList_NormTMMVarFil)[length(colnames(dgList_NormTMMVarFil))] <- "variance"
dgList_NormTMM <- dgList_NormTMMVarFil[dgList_NormTMMVarFil[, 17] >= quantile(dgList_NormTMMVarFil[, 17], quant), ]
dgList_NormTMM <- dgList_NormTMM[, 1:16]
write.table(dgList_NormTMM, paste("topGenesMostVariable_", quant*100, "percent.txt", sep=""), quote=F)  ## Change name accordingly

## Transpose data for WGCNA
fercountsTr <- as.data.frame(t(dgList_NormTMM))

## Check data for missing entries and zerp-Var genes
gsg = goodSamplesGenes(fercountsTr, verbose = 3)
gsg$allOK  ## If "TRUE" we can proceed with analysis. Otherwise, filter data.

## Outlier sample detection using clustering... This is not gene clustering
numTraits <- sampleMeta[, -c(1:6, 8:10, 12:13)]
sampleTree <- hclust(dist(fercountsTr), method = "average")
traitColors <- numbers2colors(numTraits, signed = FALSE)
pdf("sampleOutlierDetect_hclust.pdf")
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) ## When no "trait" data is available
plotDendroAndColors(sampleTree, traitColors, groupLabels=names(numTraits), main = "Sample clustering to detect outliers")
dev.off()

## Choose a set of soft-thresholding powers and call the network topology analysis function. Then plot results (Scale-free topology fit index as a function of the soft-thresholding power)
## Select a scale-free topology fit R^2 thershold
RsquaredCut = c(0.9)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(fercountsTr, powerVector = powers, verbose = 5, blockSize=20000, RsquaredCut=RsquaredCut) ## blockSize controls how many genes in each block to analyze. Ideally, one single block ( i.e. blockSize > length(colnames(ExprMatrix)) )
pdf("scaleFreeTopology_vs_SFT.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red")
abline(h=RsquaredCut, col="red")
dev.off()
pdf("MeanConnectivity_vs_SFT.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

## Estimate and save adjancency, TOM and Dissimilarity matrices
adjacency <- adjacency(fercountsTr, power=sft$powerEstimate)
write.table(adjacency, paste("adjacencyMatrix_sft_", sft$powerEstimate, ".txt", sep=""))
TOMmatrix = TOMsimilarity(adjacency)
write.table(TOMmatrix, "TOMSimilarity.txt")
dissTOM = 1-TOMmatrix
write.table(dissTOM, "disSimilarityTOM.txt")

## Cluster genes!
## Select a cutHeight value to cluster branches
cutHeight = c(0.99)
geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(paste("geneTree_cutHeight", cutHeight, ".pdf", sep=""))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
abline(h=cutHeight,col="red")
dev.off()
# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = c(40)
deepSplit = c(3)
# Module identification using dynamic tree cut: 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = deepSplit, pamRespectsDendro = FALSE,  minClusterSize = minModuleSize, cutHeight = cutHeight)
write.table(as.data.frame(table(dynamicMods)), paste("geneClusterSizes_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".txt", sep=""), row.names=F, quote=F)

# Convert numeric lables into colors 
dynamicColors = labels2colors(dynamicMods)
write.table(as.data.frame(table(dynamicColors)), paste("geneClusterSizes_Colors_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".txt", sep=""), row.names=F, quote=F)
# Plot the dendrogram and colors underneath
pdf(paste("geneTree_wColors_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".pdf", sep=""))
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

## Calculate eigengenes 
MEList = moduleEigengenes(fercountsTr, colors = dynamicColors) 
MEs = MEList$eigengenes 
## Calculate dissimilarity of module eigengenes 
MEDiss = 1-cor(MEs)
## Cluster module eigengenes 
METree = hclust(as.dist(MEDiss), method = "average")
pdf(paste("eigenGeneTree_wColors_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".pdf", sep=""), height=6, width=12)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

## Call an automatic merging function 
## Select a cluster merge threshold:
MEDissThres = c(0.3)
merge = mergeCloseModules(fercountsTr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
# The merged module colors 
mergedColors = merge$colors
write.table(as.data.frame(table(mergedColors)), paste("geneMergedClusterSizes_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".txt", sep=""), row.names=F, quote=F)
# Eigengenes of the new merged modules: 
mergedMEs = merge$newMEs
pdf(file = paste("moduleDendrogram_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".pdf", sep=""), wi = 9, he = 6) 
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) 
dev.off()

## Visualize dissimlarity between modules in a heatmap
plotTOM = dissTOM^7
diag(plotTOM) = NA
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
tiff(paste("module_DendrogramHeatmap_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".tif", sep=""), width = 1500, height = 1500, res=300)
TOMplot(plotTOM, geneTree, mergedColors, main = "Network heatmap plot, all genes", col=myheatcol)
dev.off()

geneColors <- (cbind(rownames(adjacency), dynamicColors))
write.table(as.data.frame(geneColors), paste("geneColors_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".txt", sep=""), row.names=F, quote=F)

## MODULE-TRAIT CORRELATIONS - WITH UN-MERGED MODULES
nGenes <- ncol(fercountsTr)
nSamples <- nrow(fercountsTr)
MEs0 = MEList$eigengenes
MEs0 = orderMEs(MEs0)
traits <-  model.matrix(~situation + sex, data=sampleMeta)
traits <- as.data.frame(traits)[2:3]
moduleTraitCor = cor(MEs0, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor) 

pdf(file = paste("module-Trat_Corr_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".pdf", sep=""), wi = 6, he = 10) 
par(mar = c(2, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = c("Situation", "Sex"), yLabels = names(MEs0), ySymbols = names(MEs0), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text =0.5, zlim = c(-1,1), main = paste("Module-trait relationships"),xLabelsAngle=0)
dev.off()

## MODULE-TRAIT CORRELATIONS - WITH MERGED MODULES
nGenes <- ncol(fercountsTr)
nSamples <- nrow(fercountsTr)
MEsMerged = mergedMEs
MEsMerged = orderMEs(MEsMerged)
traits <-  model.matrix(~situation + sex, data=sampleMeta)
traits <- as.data.frame(traits)[2:3]
moduleTraitCor = cor(MEsMerged, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor) 

pdf(file = paste("moduleMerged-Trat_Corr_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".pdf", sep=""), wi = 6, he = 10) 
par(mar = c(2, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = c("Situation", "Sex"), yLabels = names(MEsMerged), ySymbols = names(MEsMerged), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text =0.5, zlim = c(-1,1), main = paste("Module-trait relationships"),xLabelsAngle=0)
dev.off()

## Choose hub genes/module
hubGenes <- chooseTopHubInEachModule(fercountsTr, mergedColors, power=sft$powerEstimate, type="unsigned")
write.table(hubGenes, paste("hubGenes_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".txt", sep=""), quote=F, col.names=F, sep="\t")

## SAVE SESSION AT THIS POINT

## Export modules to Cytoscape
## Select modules to extract and threshold (correlation):
for( ModCol in names(table(mergedColors))) { ## OPENING LOOP THROUGH MODULE COLORS

selModules <- ModCol ## FOR LOOPING THROUGH COLORS

# selModules <- as.vector(levels(as.factor(dynamicColors))) ## ALL MODULES
# selModules <- c("darkred", "red", "purple", "black", "lightyellow", "blue", "turquoise")
adjThreshold <- c("0.00")
## Path to DE gene IDs:
deGenesfpath <- "~/Desktop/AlanCounts_DEGenes_ALL_Two-FourCompars.txt"
annotsfpath <- "~/Desktop/ref_Transcriptome_annotationFiles/sprotAnnots_Trinotate_UniqueGenes_SYMBOLS.txt"
socMatfpath <- "~/Desktop/Social_MatePref_genes_inNetwork_UNIPROTSYMBOLS.txt"
sPGIEGfpath <- "~/Desktop/SPG_IEG_genes_inNetwork_UNIPROTSYMBOLS.txt"
annots <- read.csv(annotsfpath)
deGenes <- read.csv(deGenesfpath)
deGenes <- as.vector(deGenes[[1]])
socMatGenes <- read.csv(socMatfpath)
sPGIEGGenes <- read.csv(sPGIEGfpath)

genes <- names(fercountsTr)
inModule <- is.finite(match(mergedColors, selModules))
modGenes <- genes[inModule]
modTOM = TOMmatrix[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)
deGenesNet <- unlist(dimnames(modTOM)[1]) %in% deGenes ## Finds TOM genes in the list of DE genes
deGenesNet <- str_replace(deGenesNet, "FALSE", "NoDE")
deGenesNet <- str_replace(deGenesNet, "TRUE", "DE")

annotsNet <- unlist(dimnames(modTOM)[1]) %in% annots[[1]] ## Finds TOM genes in the list of DE genes
annotAlt <- c()
for(i in 1:length(annotsNet)) {
	idTemp <- as.vector(dimnames(modTOM)[[1]][i])
	if(annotsNet[i] == "TRUE") {
		symbolTemp <- annots[grep(idTemp, as.vector(annots[[1]])), 2]
		annotAlt[i] <- as.vector(symbolTemp)
	}
	else {
		annotAlt[i] <- idTemp
	}
}

socMatNet <- annotAlt %in% socMatGenes[[1]] ## Finds Social/Mate Pref genes in the list of annotated genes in-network
socMatNet <- str_replace(socMatNet, "TRUE", "SocMat")
socMatNet <- str_replace(socMatNet, "FALSE", "NoSocMat")
sPGIEGNet <- annotAlt %in% sPGIEGGenes[[1]] ## Finds Social/Mate Pref genes in the list of annotated genes in-network
sPGIEGNet <- str_replace(sPGIEGNet, "TRUE", "SPGIEG")
sPGIEGNet <- str_replace(sPGIEGNet, "FALSE", "NoSPGIEG")

selModules <- paste(selModules, collapse="-")
classGenesNet <- cbind(deGenesNet, socMatNet, sPGIEGNet, mergedColors[inModule])
classGenesNet <- apply(classGenesNet, 1, paste, collapse="_")
TOMCytos = exportNetworkToCytoscape(
modTOM, edgeFile = paste("module_nets/edgesTOM_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, "_", "corrThres_", adjThreshold, selModules, ".txt", sep=""),
# modTOM, edgeFile = paste("module_nets/edgesTOM_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, "_", "corrThres_", adjThreshold, "ALLMODULES", ".txt", sep=""),
nodeFile = paste("module_nets/nodesTOM_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, "_", "corrThres_", adjThreshold, selModules, ".txt", sep=""), 
# nodeFile = paste("module_nets/nodesTOM_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, "_", "corrThres_", adjThreshold, "ALLMODULES", ".txt", sep=""),
weighted = TRUE, threshold = adjThreshold, nodeNames = modGenes, altNodeNames = annotAlt, 
nodeAttr=classGenesNet)
# nodeAttr = mergedColors[inModule])
} ## CLOSING BRACKET FOR LOOPING THROUGH MODULE COLORS


## g:Profiler Enrichment Analysis
# install.packages("gprofiler2")
library("gprofiler2")
library("plyr")

NetGenesUniProtNamesfPath <- "~/Desktop/networkgenes_GeneNamesFromUniprot.txt"
NetGenesUniProtNames <- read.table(NetGenesUniProtNamesfPath, header=F)
customBG <- as.vector(NetGenesUniProtNames[[2]])

for( ModCol in names(table(mergedColors))) { ## OPENING LOOP THROUGH MODULE COLORS

ModSelected <- ModCol ## FOR LOOPING THROUGH COLORS
adjThreshold <- c("0.00")

# ModSelected <- c("yellow")
organism=c('xtropicalis', 'hsapiens')
for(org in organism){
	ModGeneList <- read.table(paste("module_nets/nodesTOM_sft6_ModSize40_cutHeight0.99_deepSplit3_corrThres_", adjThreshold, ModSelected, ".txt", sep=""), skip=1)
	ModGeneList <- as.vector(ModGeneList[,2])
	ModGeneListMask <- NetGenesUniProtNames[[1]] %in% ModGeneList
	ModGeneListNames <- as.vector(NetGenesUniProtNames[[2]])[ModGeneListMask]
	gostres <- gost(ModGeneListNames, organism=org,  exclude_iea=T, correction_method="fdr", sources=c("GO", "KEGG", "TF"), custom_bg=customBG, user_threshold = 0.1)
	dfTest <- c(empty(gostres))
	if(is.na(dfTest)){
		gostres_pValSort <-  as.data.frame(gostres$result[order(gostres$result$p_value), ][, 1:13])	
		write.table(gostres_pValSort, file=paste("module_nets/enrichmentTable_", ModSelected, "_", org, ".csv", sep=""), quote=F, sep=",")
	}else{print("No table file generated")}
}
} ## CLOSING BRACKET FOR LOOPING THROUGH MODULE COLORS

