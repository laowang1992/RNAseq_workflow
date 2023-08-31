#!/bin/env Rscript

###########################################################################
########################### A WGCNA EXAMPLE ###############################
###########################################################################

library(WGCNA)

## -----准备输入文件-------------------------------------------------------------------
f_exp <- './data/genes.TMM.average.txt'
f_de_id <- './data/DE_gene_id.txt'
f_trait <- './data/Traits.txt'

## 创建输出文件夹
# 创建图片文件夹
#dir.create("Plots")
# 创建表格文件夹
#dir.create("Tables")

## -----导入表达矩阵------------------------------------------------------------------

# 不要将文件中的字符串转换为因子
options(stringsAsFactors = FALSE)

# 读取表达数据
exp0 <- read.delim(f_exp, row.names=1)
de_id <- read.delim(f_de_id, header=F)

# 过滤掉在所有样品中表达之和小于10 的基因, RNA-Seq 建议
# exp0 <- exp0[rowSums(exp0) > 10, ]

# 筛选出差异表达基因
exp0 <- exp0[de_id$V1, ]

# 对数，RNA-Seq 建议
exp0 <- log2(exp0 + 1)

# 转置
datExpr0 <- as.data.frame(t(exp0))

datExpr0[1:4,1:4]


## 样本聚类，检验样本

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = "Plots/sampleClustering.pdf", width = 6, height = 4);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()

# 保留左右样本
datExpr <- datExpr0


## ----寻找最佳 beta 值-------------------------------------------------------

# 开启多线程模式
enableWGCNAThreads(nThreads = 20)
# MACOS 下需要改为 
# disableWGCNAThreads()

# 通过对 power 的多次迭代，确定最佳 power
powers <- 1:20
sft <- pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 5,
                        networkType = "unsigned"
                        )
# 保存选择系数的表格
write.table(as.data.frame(sft$ fitIndices),"Tables/softPowerSelection.xls", sep="\t", quote=F)
# 计算出来的最佳 β 存放在
sft$powerEstimate

# 画图
pdf(file = "Plots/softPowerSelction.pdf", width = 9, height = 5);
par(mfrow = c(1,2))
cex1 <- 0.9
#  R2 ~ soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.80,col="green")
# Mean connectivity ~ soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


## ----构建共表达网络/聚类分析-------------------------------------------------------

# 选择软阈值，可以多尝试，这里先选用sft$powerEstimate
#softPower = sft$powerEstimate
softPower = 11

# 是否要进行自动模块聚类识别
auto <- FALSE
if (auto) {
  net <- blockwiseModules(
    # 0.输入数据
    datExpr, 
    
    # 1. 计算相关系数
    corType = "pearson", # 相关系数算法，pearson|bicor
    
    # 2. 计算邻接矩阵
    power = softPower, # 前面得到的 soft power
    networkType = "unsigned", # unsigned | signed | signed hybrid
    
    # 3. 计算 TOM 矩阵
    TOMType = "unsigned", # none | unsigned | signed
    saveTOMs = TRUE, # 是否保存
    saveTOMFileBase = "TOM", # .... 重命名
    
    # 4. 聚类并划分模块
    deepSplit = 2, # 0|1|2|3|4, 值越大得到的模块就越多越小
    minModuleSize = 30,
    
    # 5. 合并相似模块
    ## 5.1 计算模块特征向量（module eigengenes， MEs），即第一个主成分（1st PC）
    ## 5.2 计算 MEs 与 datTrait 之间的相关性
    ## 5.3 对距离小于 mergeCutHeight （1-cor）的模块进行合并
    mergeCutHeight = 0.25, 
  
    # 其他参数
    maxBlockSize = 20000, 
    numericLabels = TRUE, # 以数字命名模块
    nThreads = 10 # 0 则使用所有可用线程
    )
  # 查看每个模块包含基因数目
  table(net$colors)
  ## -----模块可视化----------------------------------------------
  
  # 模块名称修改为颜色
  pdf(file = paste0("Plots/","softPower_",softPower,"_autoDetection.pdf"), width = 9, height = 5);
  moduleColors <- labels2colors(net$colors)
  # 同时绘制聚类图和模块颜色
  plotDendroAndColors(
    net$dendrograms[[1]], 
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, 
    addGuide = TRUE)
  dev.off()
   
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  save(MEs, moduleLabels, moduleColors, geneTree, 
    file = "allSample-02-networkConstruction-auto.sft_6.RData")
}

# 手动构建模块
networkType = "unsigned"
adjacency = adjacency(datExpr, power = softPower);
cat(date(),"\t我们要计算 TOMsimilarity 了，可能要很久\n")
TOM = TOMsimilarity(adjacency);
cat(date(),"\t终于算完了，用了多长时间呢？后面应该很快了\n")

dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf(file = "Plots/geneClusteringTreeBasedOnDissTOM.pdf", width = 9, height = 5);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
  labels = FALSE, hang = 0.04)
dev.off()


# 多个deepSplit
mColorh=NULL;
minModuleSize = 30;
for (ds in 0:3){
	tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
						minClusterSize = (minModuleSize-3*ds), cutHeight = 0.99,
						deepSplit = ds, distM = dissTOM);
	mColorh=cbind(mColorh,labels2colors(tree$labels));
}
# pdf("Module_choices.pdf", height=10,width=25);
pdf(file = "Plots/deepSplitChoices.pdf", width = 15, height = 10);
plotDendroAndColors(geneTree, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);
dev.off()
# 查看每个deepSplit对应的模块数目
str(apply(mColorh,2,table))

# 查看之后，决定用 deepSplit=2

# -------------------
finalDs = 2;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = finalDs, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

# 
#------------------------
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
pdf(file = "Plots/geneDendrogramAndModuleColors.pdf", width = 9, height = 5);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# 补充代码，关闭绘图设备，再下载图片到本地可查看pdf结果
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
# sizeGrWindow(7, 6)
pdf(file = "Plots/clusteringOfModuleEigegenes.pdf", width = 14, height = 6);
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# 限定最小差异为0.10
MEDissThres = 0.10	
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")	# 画一条线 方便查看
# 限定最小差异为0.10
MEDissThres = 0.20
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "green")	# 画一条线 方便查看
# 限定最小差异为0.10
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "blue")	# 画一条线 方便查看

# 补充代码，关闭绘图设备，再下载图片到本地可查看pdf结果
dev.off()

##########
#	以下是新增加的部分，直接查看模块切分结果，再决定
########
merge0.05 = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.05, verbose = 3)
merge0.10 = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.10, verbose = 3)
merge0.15 = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.15, verbose = 3)
merge0.20 = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.20, verbose = 3)
merge0.25 = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
# The merged module colors
mergedColors0.05 = merge0.05$colors;
mergedColors0.10 = merge0.10$colors;
mergedColors0.15 = merge0.15$colors;
mergedColors0.20 = merge0.20$colors;
mergedColors0.25 = merge0.25$colors;
# Eigengenes of the new merged modules:
# mergedMEs = merge$newMEs;

# sizeGrWindow(12, 9)
pdf(file = "Plots/geneDendrogramAndModuleColorsAndDyamicColorsChoose.pdf", wi = 12, he = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors0.05,mergedColors0.10,mergedColors0.15,mergedColors0.20,mergedColors0.25),
					# dynamicColors 为原来的颜色
					# mergedColors 为合并后的颜色
					c("Dynamic Tree Cut",  "MEDissThres 0.05","MEDissThres 0.10","MEDissThres 0.15", "MEDissThres 0.20", "MEDissThres 0.25"),
					dendroLabels = FALSE, hang = 0.03,
					addGuide = TRUE, guideHang = 0.05)
dev.off()
# 输出模块数目比较？
moduleCounts<-data.frame(cutHeight=c(0.05,0.10,0.15,0.20,0.25),moduleCounts=c(length(table(mergedColors0.05)),length(table(mergedColors0.10)),length(table(mergedColors0.15)),length(table(mergedColors0.20)),length(table(mergedColors0.25))))

write.table(moduleCounts,
	file="Tables/cutHeightChoose.xls",
	row.names = FALSE,
	col.names = FALSE,
	quote = FALSE
	)

# 确定的参数
MEDissThres = 0.25

# Call an automatic merging function	# 进行模块合并
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# sizeGrWindow(12, 9)
pdf(file = "Plots/geneDendrogramAndModuleColorsAndDyamicColors.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    # dynamicColors 为原来的颜色
					# mergedColors 为合并后的颜色
					c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

######
#	创建参数文件
######
system(paste0("touch Tables/UsedParamers","_softPower_",softPower,"_networkType_",networkType,"_deepSplit_",finalDs,"_minModuleSize_",minModuleSize,"_MEDissThres_",MEDissThres))


#######
#	输出模块对应的GID
#########
# softPower = 11
gids <- colnames(datExpr)
intModules <- gsub("ME","",names(mergedMEs))
for(module in intModules){
	modGenes <- gids[mergedColors == module]
	fileName<-paste("Tables/ModuleGID_",module,".txt",sep="")
	write.table(as.data.frame(modGenes),
		file=fileName,
		row.names = FALSE,
		col.names = FALSE,
		quote = FALSE
		)
}

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "allSample_ave-networkConstruction-stepByStep.RData")

##########
#	有些人可能想要module和sample的相关性图片
##########
# 计算样本相关系数
# datCor<-cor(t(datExpr))
datCor<-cor(t(datExpr),use = 'pairwise.complete.obs')
# 
nGenes = ncol(datExpr)
nSambles = nrow(datExpr)
# 
MEs0 = moduleEigengenes(datExpr,mergedColors)$eigengenes
MEs = orderMEs(MEs0)
moduleSampleCor = cor(MEs,datCor,use='p')
moduleSamplePvalue = corPvalueStudent(moduleSampleCor,nSambles)
# 绘图
# sizeGrWindow(20,16)
pdf(file = "Plots/moduleSampleCorrelationGreenWhiteRed.pdf", wi = 30, he = 15)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleSampleCor, 2), "\n(",
                           signif(moduleSamplePvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleSampleCor)
par(mar = c(6, 15, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleSampleCor,
               xLabels = colnames(datCor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-Sample relationships"))
dev.off()

# 改变色彩
pdf(file = "Plots/moduleSampleCorrelationBlueWhiteRed.pdf", wi = 30, he = 15)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleSampleCor, 2), "\n(",
                           signif(moduleSamplePvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleSampleCor)
par(mar = c(6, 15, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleSampleCor,
               xLabels = colnames(datCor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-Sample relationships"))
dev.off()

##############
# 如果有性状数据
############
# 导入traits信息
traits <- read.table(f_trait,sep="\t",header=T,row.names=1)
# 这里需注意有个大坑，traits的样本顺序需要和MEs中的样本顺序一致
# 因此，需要重新排序
traits <- traits[rownames(MEs), ]

moduleTraitCor = cor(MEs,traits,use='p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nSambles)
# 绘图
# sizeGrWindow(20,16)
pdf(file = "Plots/moduleTraitCorrelationGreenWhiteRed.pdf", wi = 8, he = 12)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 15, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-Traits relationships"))
dev.off()
#####
pdf(file = "Plots/moduleTraitCorrelationBlueWhiteRed.pdf", wi = 8, he = 12)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 15, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-Traits relationships"))
dev.off()



# 输出模块与样本的相关系数
write.table(moduleSampleCor,"Tables/moduleSampleCor.xls", sep="\t", quote=F)
write.table(moduleSamplePvalue,"Tables/moduleSampleCorPvalue.xls", sep="\t", quote=F)



#########
#	输出数据，用于可视化
##############
#######
#	做一些另外的东西 比如模块表达量
#########

# 查看模块之间的相关系数
datME<-moduleEigengenes(datExpr,mergedColors)$eigengenes
# 
MEcor = signif(cor(datME,use="p"),2);
write.table(MEcor,file="Tables/ModuleCorelation.xls",sep="\t")

# 查看聚类图 keeps track of the sign 
dissimME = (1-MEcor)/2
hclustdatME = hclust(as.dist(dissimME),method="average")
# 绘制
sizeGrWindow(9,9)
pdf(file = "Plots/ClusteringOfModuleEigengenes1.pdf", wi = 8, he = 5)
par(mfrow=c(1,1))
plot(hclustdatME,main="Clustering tree based of the module eigengenes")
dev.off()

# Pairwise scatter plots of the samples along the module eigengenes
sizeGrWindow(8,9)
pdf(file = "Plots/moduleEigengenesCorrelation.pdf", wi = 12, he = 9)
# plotMEpairs(datME,y=y) # y对应一个性状 -- 这个y是怎么得来的？
# plotMEpairs(datME,y=moduleSampleCor) # 失败 应该是对应 moduleSampleCor 的某一列
# 直接只绘制模块之间的
plotMEpairs(datME)
dev.off() 

## 保存导出模块的信息
save(mergedColors, datExpr, TOM, file = "export_network.RData")
## 导出模块信息
allModule <- names(table(mergedColors))
for (modules in allModule) {
  probes = colnames(datExpr)
  inModule = is.finite(match(mergedColors, modules));
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  # 输出所有模块的，暂时注释掉
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("Tables/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("Tables/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 nodeAttr = mergedColors[inModule])
}
