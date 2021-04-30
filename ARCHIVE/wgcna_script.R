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

setwd("~/Desktop/wgcna_RNASeq/")

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
write.table(as.data.frame(table(dynamicColors)), paste("geneClusterSizes_ColorsModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, ".txt", sep=""), row.names=F, quote=F)
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

## Export modules to Cytoscape
## Select modules to extract and threshold (correlation):
selModules <- c("brown")
adjThreshold <- c("0.3")
## Path to DE gene IDs:
deGenesfpath <- "~/Desktop/AlanCounts_DEGenes_ALL.txt"
annotsfpath <- "~/Desktop/Trinotate_Genes_FirstAnnot_wSYMBOLS.txt"
annots <- read.csv(annotsfpath)
deGenes <- read.csv(deGenesfpath)
deGenes <- as.vector(deGenes[[1]])

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
		annotAlt[i] <- symbolTemp
	}
	else {
		annotAlt[i] <- idTemp
	}
}

TOMCytos = exportNetworkToCytoscape(
modTOM, edgeFile = paste("edgesTOM_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, "_", "corrThres_", adjThreshold, selModules, ".txt", sep=""),
nodeFile = paste("nodesTOM_sft", sft$powerEstimate, "_ModSize", minModuleSize, "_cutHeight", cutHeight, "_deepSplit", deepSplit, "_", "corrThres_", adjThreshold, selModules, ".txt", sep=""), 
weighted = TRUE, threshold = adjThreshold, nodeNames = modGenes, altNodeNames = annotAlt, 
nodeAttr=deGenesNet)
# nodeAttr = mergedColors[inModule])
