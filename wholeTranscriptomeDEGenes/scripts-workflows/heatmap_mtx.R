##
# This function takes a DGEList and the result table of a contrast test, and 
# outputs a matrix for ComplexHeatmaps. The genes to be plotted are extracted
# from the contrast table.
#
# @param dgelist, the DGEList object.
# @param degenes, a vector with gene names matching names of genes in the DGEList.
# @return, a matrix for ComplexHeatmaps.
#
#
heatmap_mtx <- function(dgelist=NULL, contrasts=NULL){
  degenes <- contrasts$table$genes
  de_counts <- dgelist$counts[is.element(unlist(dgelist$genes), degenes), ]
  de_genes <- dgelist$genes[is.element(unlist(dgelist$genes), degenes), ]
  logcpm <- cpm(de_counts, prior.count=2, log=TRUE)
  z_logcpm <- t(scale(t(logcpm), scale=T, center=T))
  z_logcpm <- merge(z_logcpm, contrasts$table, by=0)
  z_logcpm <- z_logcpm[order(z_logcpm$logFC), ]
  z_logcpm <- z_logcpm[, -((length(colnames(z_logcpm))-4):length(colnames(z_logcpm)))]
  rownames(z_logcpm) <- z_logcpm$genes
  z_logcpm <- z_logcpm[,-c(1, length(colnames(z_logcpm)))]
  z_logcpm <- as.matrix(z_logcpm)

  return(z_logcpm)
}