# 
# write.csv(new.data0,file = "my_counts.csv",row.names = TRUE,quote = FALSE)
# write.csv(new.sample.info,file = "input_data/sample_info.csv",row.names = FALSE,quote = FALSE)
# write.csv(bac_traits,file = "input_data/traits.csv",row.names = T,quote = FALSE)
setwd("data/20220623/NP/Yu2021NaturePlants-master/input_data/")
library(DESeq2)
library(WGCNA)
enableWGCNAThreads(nThreads = 4)
data0<-read.table("my_counts.csv",
                  header=T,
                  row.names=1,
                  sep=",")
head(data0)
sample_metadata <- read.csv(file = "sample_info.csv",
                            header = TRUE)
head(sample_metadata)

bac_traits<-read.csv("input_data/traits.csv",row.names = 1)
head(bac_traits)


dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-20],
                                         colData = sample_metadata,
                                         design = ~ Zone)
mcols(dataExpr_deseq)$basepairs = data0$geneLength
fpkm_matrix = fpkm(dataExpr_deseq)
datExpr = t(log2(fpkm_matrix+1))

dim(datExpr)

nSamples = nrow(datExpr)  # 统计样品数目
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))  #按列（基因）取方差
no.missingdatExpr=as.vector(apply(is.na(as.matrix(datExpr)),2, sum) )#按列（基因）统计缺失数目
KeepGenes= variancedatExpr>0 & no.missingdatExpr<0.1*nSamples  # 保留方差不等于0，且缺失低于10%的基因
table(KeepGenes)# 过滤统计

datExpr=datExpr[, KeepGenes]

setwd("input_data/")
sampleTree = hclust(dist(datExpr), method = "average")
pdf(file = "1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.5);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

power=sft$powerEstimate
power

## long time big mem
## multiple cpu
#TOM = TOMsimilarityFromExpr(datExpr, power = power)
## 这一步用linux系统的电脑设置多核心还挺快的 enableWGCNAThreads(nThreads = 4)，
## 我设置12核心很快就得到了结果
## 但是我在自己的win系统上设置多核好像没有作用
load("../TOM.rdata")
dissTOM = 1-TOM 
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, 
     xlab="", 
     sub="", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, 
     hang = 0.04);
dev.off()

dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = 2, 
                            pamRespectsDendro = FALSE,
                            minClusterSize = 30)

dynamicColors = labels2colors(dynamicMods)

pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, 
                    dynamicColors, 
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEDissThres=0.25

merge = mergeCloseModules(datExpr, 
                          dynamicColors, 
                          cutHeight = MEDissThres, 
                          verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, 
                    cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05)  
dev.off()
write.table(merge$oldMEs,file="oldMEs.txt");
write.table(merge$newMEs,file="newMEs.txt");

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

table(rownames(MEs) == rownames(bac_traits))
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("module-traits-bacteria-order.pdf", width = 10, height = 10)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bac_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
