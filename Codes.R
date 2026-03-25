getwd()
setwd("~/GSE171306")


library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(SingleR)
library(clustree)
library(cowplot)
library(tidyverse)
library(SCpubr)
library(GSVA)
library(GSEABase)
library(hdf5r)         #读入h5文件需要引用奥
library(plyr)
library(harmony)


library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(clustree)
library(cowplot)
library(tidyverse)
library(SCpubr)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(ggpubr)



packageVersion("Seurat")
packageVersion("Matrix")




ids <- "KIRC_GSE171306_expression.h5"

data1 = Read10X_h5("KIRC_GSE171306_expression.h5")
afcol=gsub(".h5","",ids)
#创建seurat对象
# min.cell：每个feature至少在多少个细胞中表达
# min.features：每个细胞中至少有多少个feature被检测到
#nFeature_RNA是每个细胞中检测到的基因数量
#nCount_RNA是细胞内检测到的分子总数
#nFeature_RNA过低，表示该细胞可能已死/将死或是空液滴。高nCount_RNA和/或nFeature_RNA表明“细胞”实际上可以是两个或多个细胞。
#结合线粒体基因（percent.mt）与核糖体基因（percent.rb）除去异常值，即可除去大多数双峰/死细胞/空液滴，因此它们过滤是常见的预处理步骤
af<-CreateSeuratObject(counts =data1, 
                       project = afcol, 
                       min.cells = 3,
                       min.features = 200)
af$Type=Idents(af)


af[["percent.mt"]] <- PercentageFeatureSet(af, pattern = "^MT-")
af[["percent.rb"]] <- PercentageFeatureSet(af, pattern = "^RP")
mask1 <- af$nCount_RNA >= 1000
mask2 <- af$nFeature_RNA >= 200 & af$nFeature_RNA <= 10000
mask3 <- af$percent.mt <= 20
mask4<-af$percent.rb<= 20
mask <- mask1 & mask2 & mask3 & mask4
af <- af[, mask]




# QC指标使用小提琴图可视化,ncol为图片排列列数
pdf(file = "01.vlnplot.pdf",width = 12,height = 5)
VlnPlot(af, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)+scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))
dev.off()

# 指标之间的相关性
plot1 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.mt")+ RotatedAxis()
plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
#组图
pdf(file = "01.corqc.pdf",width =12,height = 5)
plot1+plot2+plot3+plot_layout(ncol = 3)      #plot_layout，patchwork函数，指定一行有几个图片
dev.off()

##标准化,使用LogNormalize方法
af <- NormalizeData(af, normalization.method = "LogNormalize", scale.factor = 10000)

## 鉴定高变基因
# 高变基因：在一些细胞中表达高，另一些细胞中表达低的基因
# 变异指标： mean-variance relationship
# 返回2000个高变基因，用于下游如PCA降维分析。
af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)

# 提取前10的高变基因
top10 <- head(VariableFeatures(af), 10)
top10


# 展示高变基因
plot1 <- VariableFeaturePlot(af)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "01.topgene.pdf",width =7,height = 6)
plot2                   #plot_layout，patchwork函数，指定一行有几个图片
dev.off()

## 归一化
# 归一化处理：每一个基因在所有细胞中的均值变为0，方差标为1，对于降维来说是必需步骤
# 归一化后的值保存在：af[["RNA"]]@scale.data
af <- ScaleData(af)

# 可以选择全部基因归一化
all.genes <- rownames(af)
af <- ScaleData(af, features = all.genes)

##########0.3.part3 降维(绘制原始分布)##########################
# PCA降维，用前面1500个高变基因，可以使用features改变用于降维的基因集
af <- Seurat::RunPCA(af, features = VariableFeatures(object = af))
af <- Seurat::RunTSNE(af,dims = 1:20)
pdf(file = "02.rawtsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right") #top为图列位置最上方，除此之外还有right、left、bottom(意思同英文)
dev.off()
pdf(file = "02.rawpca.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "pca",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

colaa=distinctColorPalette(100)
pdf(file = "02.raw.tsne.split.pdf",width =8,height =7)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.1
) #选择展示的主成分，这边是PC2与PC1
dev.off()


af$orig.ident <- substr(rownames(af@meta.data),1,10)
#####################################（选做 harmony 去批次与降维）################################################################################################
af <- RunHarmony(af, group.by.vars = "orig.ident")
###########################################################################################################
#####################################################################
#######################################0.3.part3 矫正后结果可视化
#################################################################
# PCA降维，用前面1500个高变基因，可以使用features改变用于降维的基因集
pdf(file = "03.harmony.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "harmony",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
af <- Seurat::RunTSNE(af,dims = 1:20,reduction ='harmony')
pdf(file = "03.tsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$orig.ident))
pdf(file = "03.tsne.split.pdf",width =12,height = 7)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "orig.ident",pt.size =0.1
) #选择展示的主成分，这边是PC2与PC1
dev.off()


#########################################################################################################################################
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$Type))
# 前两个PC特征基因可视化
VizDimLoadings(af, dims = 1:2, reduction = "pca")
#热图可视化前15个PC
pdf(file = "04.pc_heatmap.pdf",width =7.5,height = 9)
DimHeatmap(af, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()
##确定使用PC个数
# each PC essentially representing a ‘metafeature’
af <- JackStraw(af, num.replicate = 100)
af <- ScoreJackStraw(af, dims = 1:20)
pdf(file = "04.jackstrawplot.pdf",width =7.5,height = 5.5)
JackStrawPlot(af, dims = 1:20)
dev.off()
pdf(file = "04.ElbowPlot.pdf",width =5,height = 4)
ElbowPlot(af,ndims = 30,reduction = "harmony")
dev.off()
pdf(file = "04.ElbowPlot_pca.pdf",width =5,height = 4)
ElbowPlot(af,ndims = 30,reduction = "pca")
dev.off()


#选择PC
afPC=7
##对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图
#af <- FindNeighbors(af, dims = 1:15)
#设置不同的分辨率，观察分群效果，dim为PCA选择的主成分数
af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
  af=FindClusters(af, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(af@meta.data[,grep("RNA_snn_res",colnames(af@meta.data))],2,table)

p2_tree=clustree(af@meta.data, prefix = "RNA_snn_res.")
pdf(file = "04.clustertree.pdf",width =12,height =10)
p2_tree
dev.off()



# 聚类并最优化
# resolution参数：值越大，细胞分群数越多，根据前面进行选择
# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
# Optimal resolution often increases for larger datasets. 


#选择分辨率进行降维
af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
af <- FindClusters(af, resolution = 1.2)

# 查看聚类数ID
head(Idents(af), 5)

# 查看每个类别多少个细胞
head(af@meta.data)
table(af@meta.data$seurat_clusters)
# 鉴定各个细胞集群的标志基因only.pos：只保留上调差异表达的基因
af.markers <- FindAllMarkers(af, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "05.cluster_markers.csv")

## 将细胞在低维空间可视化UMAP/tSNE
af <- RunUMAP(af, dims = 1:afPC, reduction = "harmony")
af <- RunTSNE(af, dims = 1:afPC, reduction = "harmony")

# 可视化UMAP/tSNE
pdf(file = "05-cluster.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "05-cluster.TSEN.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


###################################0.5.part5 singleR细胞注释#############################################
refdata=celldex::HumanPrimaryCellAtlasData()    # 网络为公用可能加载不好，用手机开热点吧
#saveRDS(refdata, "MouseRNAseqData_refdata.rds")
refdata
head(colnames(refdata))
head(rownames(refdata))

# 查看共有多少种细胞类型
unique(refdata@colData@listData[["label.main"]])
# 使用的数据为标化后的数据
testdata <- GetAssayData(af, slot="data")
dim(testdata)
testdata[1:30,1:4]
clusters <- af@meta.data$seurat_clusters
refdata <- singleR_use[refdata]$ref_min
table(clusters)
table(refdata[,1])
cellpred <- SingleR(test = testdata,  
                    ref = refdata, 
                    labels = refdata$label.main,
                    method = "cluster", 
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

str(cellpred,max.level = 3)
metadata <- cellpred@metadata
head(metadata)

celltype = data.frame(ClusterID = rownames(cellpred), 
                      celltype = cellpred$labels, 
                      stringsAsFactors = F)
celltype
write.csv(celltype, "07.singleR.celltype_anno_SingleR.csv")
# 打分热图上面的注释结果需要校正
pdf(file = "07-singleR.pdf",width =7.5,height = 5.5)
p = plotScoreHeatmap(cellpred, clusters = rownames(cellpred), order.by = "cluster")
p
dev.off()
#########sigleR注释后结果可视化
newLabels=cellpred$labels
names(newLabels)=levels(af)
af=RenameIdents(af, newLabels)
# 可视化UMAP/tSNE
pdf(file = "07-scRNA.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3,pt.size = 0.1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "07-scRNA.TSEN.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3,pt.size = 0.1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


af$Celltype=Idents(af)





###############################################
#########################################################
################人工注释###########################################################################################################
#实用网址：
#https://www.thermofisher.cn/cn/zh/home/life-science/cell-analysis/cell-analysis-learning-center/immunology-at-work.html
#http://xteam.xbio.top/CellMarker/
#https://www.jianshu.com/p/15dddefc7038
#https://toppgene.cchmc.org/enrichment.jsp
#genes <- list("Cancer stem cells" = c("PROM1","CD34","CD90"),
#              "Monocyte"=c("CD14")
#              "M1 macrophages" = c("CD16", "FCGR3B","FCGR1A"),
#              "M2 macrophages" = c("MSR1","CD163","MRC1","CSF1R"),
#              "CD8+ T cells" = c("GZMK", "CD8A","CD8B"),
#              "CD4+ memory cells" = c("IL7R", "CD27","CCR7"),
#              "B cells" = c("CD79A", "CD37"),
#              "Regulatory T cells" = c("LAG3", "ITGA2","FOXP3","HELIOS","NRP1"),
#              "NK cells" = c("CD160","NKG7", "GNLY", "CD247", "CCL3", "GZMB"),
#              "Fibroblasts" = c("FGF7", "MME"),
#              "Endothelial cells" = c("PECAM1", "VWF"),
#              "Neurons" = c("ENO2")
#              "Epithelial_cells"=c("cD24","CDH1","CLDN4")
#)
af <- FindClusters(af, resolution = 1.2)
genes <- list("Cancer stem cells" = c("PROM1","CD34","CD90"),
              "Monocyte"=c("CD14","FCN1","APOBEC3A","THBS1"),
              "M1 macrophages" = c("CD16", "FCGR3B","FCGR1A"),
              "M2 macrophages" = c("MSR1","CD163","MRC1","CSF1R"),
              "Macrophage" = c("CD163","CD68","CD14","FCGR3A"),
              "Exhausted CD8+ T cell" =c("ABCG1", "ACP5", "ACSL"),
              "CD8+ T cells" = c("GZMK", "CD8A","CD8B"),
              "CD4+ T ells" = c("CD4", "AQP3", "GPR183"),
              "B cells" = c("CD79A", "CD37"),
              "Regulatory T cells" = c("LAG3", "ITGA2","FOXP3","HELIOS","NRP1"),
              "NK cells" = c("CD160","NKG7", "GNLY", "CD247", "CCL3", "GZMB"),
              "Fibroblasts" = c("COL1A2", "MYL9","COL1A1","MMP2"),
              "Endothelial cells" = c("PECAM1", "VWF","CHGB"),
              "Neurons" = c("ENO2"),
              "Epithelial_cells"=c("cD24","CDH1","CLDN4"),
              "Mast cell" = c("KIT","CPA3","TPSAB1"),
              "T profile" = c("MKI67"),
              "Malignant" = c("CDH1","CD24","MYC"),
              "Epi" =c("EPCAM","KRT18","KLK3", "ACPP", "MSMB")
              
)
pdf(file = "08.ann_cluster_marker.pdf",width =33,height = 15)
do_DotPlot(sample = af,features = genes,dot.scale = 12,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()
#人工注释
table(af@active.ident) # 看一下现在的cluster数量
celltype   # 看一下singleR软件提供的注释
ann.ids <- c("CD8+ T cell",  #cluster0
             "Epithelial (Maligant)",  #cluster1
             "CD4+ T cell",      #以下按顺序操作
             "CD8+ T cell",
             "Endothelial",
             "NK cell",
             "Macrophage",
             "NK cell",
             "Fibroblast",
             "Epithelial (Maligant)",
             "Endothelial",
             "Epithelial (Maligant)",
             "Monocyte",
             "Macrophage",
             "Epithelial (Maligant)",
             "Macrophage",
             "Mast cell",
             "B cell",
             "NK cell",
             "Fibroblast",
             "Epithelial (Maligant)",
             "Endothelial",
             "Macrophage",
             "Endothelial"
)
length(ann.ids)    # 看看长度对不对
afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens
af$Celltype=Idents(af)


#########人工注释后结果可视化
# 可视化UMAP/tSNE
pdf(file = "08-ann.scRNA.UMAP.pdf",width =6.5,height = 4.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 0.1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "08-ann.scRNA.TSEN.pdf",width =6.5,height = 4.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 0.1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()




#########用SCP美化作图
library(SCP)

##经典
pdf(file = "08-ann.scRNA.UMAP(美化 2).pdf",width =7,height = 5)
CellDimPlot(af,group.by="Celltype",reduction="UMAP")
dev.off()


##小箭头
pdf(file = "08-ann.scRNA.UMAP(美化)(clusters).pdf",width =7,height = 5)
CellDimPlot(af,group.by="Celltype",reduction="UMAP",theme_use="theme_blank")
dev.off()


###clusters
pdf(file = "08-ann.scRNA.UMAP(美化)(clusters).pdf",width =7,height = 5)
CellDimPlot(af,group.by="seurat_clusters",reduction="UMAP",theme_use="theme_blank")
dev.off()


pdf(file = "08-ann.scRNA.UMAP(美化)(clusters 2).pdf",width =7,height = 5)
CellDimPlot(af,group.by="seurat_clusters",reduction="UMAP")
dev.off()



##########用scRNAtoolVis美化##########

library(scRNAtoolVis)

colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
            '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
            '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
            '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
            '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
            '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
            '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
            '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
            '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
            '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
            '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
            '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')

colors2 <- c('#d4de9c','#94c58f','#86c7b4','#9cd2ed','#a992c0',
             '#ea9994','#f2c396','#bb82b1')

##cluster
pdf(file = "08-ann.scRNA.UMAP(美化)(clusters 3).pdf",width =7,height = 5)
clusterCornerAxes(object = af, 
                  reduction = 'umap',
                  pSize=0.1,
                  clusterCol = "seurat_clusters",
                  noSplit = T,
                  cellLabel = T,
                  cellLabelSize = 3.5,
                  show.legend = T) +
  scale_color_manual(values=colors)
dev.off()


##细胞群
pdf(file = "08-ann.scRNA.UMAP(美化 3).pdf",width =7,height = 5)
clusterCornerAxes(object = af, 
                  reduction = 'umap',
                  pSize=0.1,
                  clusterCol = "Celltype",
                  noSplit = T,
                  cellLabel = T,
                  cellLabelSize = 3.5,
                  show.legend = T) +
  scale_color_manual(values=colors)
dev.off()




######## Feature plot
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)

gene<-c("CD8A","CD8B","NKG7","GNLY","CD4","GPR183",
        "CD37","CD79A","CPA3","TPSAB1","MYL9","VWF",
        "PECAM1","KRT18","CD24","FCN1","CD163","CD68")


pdf(file = "08-ann.scRNA.UMAP(feature plot).pdf",width =20,height = 20)
FeaturePlot(af,features = gene,cols = c("#4b6aa8", "white", "#9d3b62")) & NoAxes()
dev.off()



FeaturePlot(af,features = "TNFAIP3",cols = c("#4b6aa8", "white", "#9d3b62")) & NoAxes()




############## Circle umap ################
library(plot1cell)


cors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
          '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
          '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
          '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
          '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
          '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
          '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
          '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
          '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
          '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
          '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
          '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')



cors_group <- c("#3674a2","#8d689d") #定义颜色



# Prepare data for ploting 准备圈图数据
circ_data <- prepare_circlize_data(af, scale = 0.8)
# plot and save figures
pdf(file ='scRNA_celltype_circle_umap.pdf', width = 4.5, height = 4.5)
plot_circlize(circ_data,do.label = T, pt.size = 0.05, col.use = cors[1:10], bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
add_track(circ_data, group = "orig.ident", colors = cors_group, track_num = 2)
dev.off()





# 鉴定各个细胞集群的标志基因only.pos：只保留上调差异表达的基因
af.markers <- FindAllMarkers(af, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "08.cell_markers.csv")








#######  AUCELL
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
library(monocle)
library(tidyverse)
library(CellChat)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(foreach)
library(KernSmooth)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(pheatmap)
library(ggheatmap)
library(reshape2)
library(SummarizedExperiment)
library(ggpubr)
library(ggrepel)
library(SeuratWrappers)
library(SeuratData)
library(copykat)
library(limma)
library(GSEABase)
library(GSVA)
library(readxl)

setwd("~/triapotosis")

#读入
pbmc <- af
DefaultAssay(pbmc) <- "RNA"
#标准化
pbmc <- NormalizeData(pbmc)
dim(pbmc)

#获取表达矩阵
expr = as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "data"))
#获取平均表达量
#expr2 = as.matrix(AverageExpression(pbmc, assays = "RNA", layer = "data", group.by = "Is_Double")[[1]])
dim(expr)

#筛选
expr <- expr[rowSums(expr)>0,]
dim(expr)

#读入基因集，将你想要打分的基因集，放入genesets.xlsx
#xlsx中第二列需为对基因集的描述，如果没有，设为na
geneSets2 = read_excel("triaptosis.xlsx",col_names = F)

#导出为gmt文件
write.table(geneSets2,file = "my_genesets.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

#读入gmt文件
geneSets=getGmt("my_genesets.gmt", geneIdType=SymbolIdentifier())

#AUCell，排序,计算评分
cells_rankings <- AUCell_buildRankings(expr)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores = 10)
gsvaResult <- getAUC(cells_AUC)
dev.off()

#导出
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="ssgseaOut.txt", sep="\t", quote=F, col.names=F)

#添加到meta.data中
pbmc <- AddMetaData(pbmc, metadata = t(gsvaResult))

#画图
p <- FeaturePlot(pbmc, features = rownames(gsvaResult),reduction = "umap", cols = c("#C3DBD9", "#990000"))
ggsave("Mono_Macro triaptosis_score_umap2.pdf", p, width = 5.5, height = 5)



p <- FeaturePlot(scRNA_m3, features = rownames(gsvaResult),reduction = "umap", cols = c("#C3DBD9", "#990000")) +
  NoAxes()+ labs(x="UMAP1",y="UMAP2") +
  theme_dr() + theme(panel.grid = element_blank())


p2 <- SCP::FeatureDimPlot(pbmc, features = rownames(gsvaResult), reduction = "umap", theme_use="theme_blank")
ggsave("aucell_scp.pdf", p2, width = 5, height = 5)



library(scRNAtoolVis)

scRNAtoolVis::featureCornerAxes(object = pbmc, reduction = "umap",
                                features = "Triaptosis")






####### Cytotrace
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(monocle)
library(ggsci)


#### CytoTRACE ####
## 目的：判断分化潜能最强的簇为干细胞，作为monocle3的起点
library(CytoTRACE)

scRNAsub <- scRNA_epi
##提取表型文件##提取no##提取表型文件##提取no表型文件
table(scRNAsub)
phe <- scRNAsub$Celltype
phe = as.character(phe)
names(phe) <- rownames(scRNAsub@meta.data)

##提取表达矩阵
mat <- as.matrix(scRNAsub@assays$RNA$counts)
mat[1:4,1:4]


##运行CytoTRACE
results <- CytoTRACE(mat = mat)
plotCytoGenes(results, numOfGenes = 10,outputDir = "cytoTRACE")

plotCytoTRACE(results,
              phenotype = phe,
              gene = "MKI67",
              outputDir = "cytoTRACE")
#结果显示Epi3为干性最强的簇，作为monocle3的起点

##回贴到seurat,个性化绘图
scRNAsub$CytoTRACE <- results$CytoTRACE
scRNA_epi$CytoTRACE <- results$CytoTRACE



p <- SCP::FeatureDimPlot(scRNAsub,features = "CytoTRACE", reduction = "umap",
                         ,theme_use="theme_blank")
ggsave("Epi_CytoTRACE_SCP.pdf",p,width = 5.5,height = 5)






######### Metabolism
library(Seurat)
library(stringr)
library(dplyr)
library(future)
library(future.apply)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(harmony)
library(clustree)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(org.Hs.eg.db)
library(tidyverse)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(AUCell)
library(ggsci)


cors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
          '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
          '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
          '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
          '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
          '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
          '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
          '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
          '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
          '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
          '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
          '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')

cors_group <- c("#3C5488FF","#BB0021FF") #定义颜色



####代谢评估 #### 
library(scMetabolism)
library(ggplot2)
library(rsvd)

#######V4运行
#scRNA<-sc.metabolism.Seurat(obj = scRNA_epi, method = "AUCell", 
#                            imputation = F, metabolism.type = "KEGG")


#######V5的函数
sc.metabolism.SeuratV5 <- function (obj,method="VISION",imputation=F,ncores=2,
                                    metabolism.type="KEGG")
{
  countexp<-GetAssayData(obj,layer='counts')
  countexp<-data.frame(as.matrix(countexp))
  signatures_KEGG_metab<-system.file("data","KEGG_metabolism_nc.gmt",
                                     package="scMetabolism")
  signatures_REACTOME_metab<-system.file("data","REACTOME_metabolism.gmt",
                                         package="scMetabolism")
  if(metabolism.type=="KEGG"){
    gmtFile<-signatures_KEGG_metab
    cat("Your choice is: KEGG\n")
  }
  if(metabolism.type=="REACTOME"){
    gmtFile<-signatures_REACTOME_metab
    cat("Your choice is: REACTOME\n")
  }
  if(imputation==F){
    countexp2<-countexp
  }
  if(imputation==T){
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- row.names(countexp)
  }
  cat("Start quantify the metabolism activity...\n")
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
    vis <- Vision(scaled_counts,signatures=gmtFile)
    options(mc.cores=ncores)
    vis<-analyze(vis)
    signature_exp<-data.frame(t(vis@SigScores))
  }
  if(method=="AUCell"){
    library(AUCell)
    library(GSEABase)
    cells_rankings<-AUCell_buildRankings(as.matrix(countexp2),
                                         nCores=ncores,plotStats=F)
    geneSets<-getGmt(gmtFile)
    cells_AUC<-AUCell_calcAUC(geneSets,cells_rankings)
    signature_exp<-data.frame(getAUC(cells_AUC))
  }
  if(method=="ssGSEA"){
    library(GSVA)
    library(GSEABase)
    geneSets<-getGmt(gmtFile)
    gsva_es<-gsva(as.matrix(countexp2),geneSets,method=c("ssgsea"),
                  kcdf=c("Poisson"),parallel.sz=ncores)
    signature_exp<-data.frame(gsva_es)
  }
  if(method=="ssGSEA"){
    library(GSVA)
    library(GSEABase)
    geneSets<-getGmt(gmtFile)
    gsva_es<-gsva(as.matrix(countexp2),geneSets,method=c("gsva"),
                  kcdf=c("Poisson"),parallel.sz=ncores)
    signature_exp<-data.frame(gsva_es)
  }
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  obj@assays$METABOLISM$score <- signature_exp
  obj
}



####运行V5的sc metabolism
res <-sc.metabolism.SeuratV5(obj = scRNA_epi,
                             method = "AUCell",#VISION、AUCell、ssgsea和gsva
                             imputation=F, ncores=2,
                             metabolism.type = "KEGG") #KEGG和REACTOME



# 提取代谢评估通路
input.pathway <- rownames(res@assays[["METABOLISM"]][["score"]])


# 糖代谢：1-15；脂代谢：18-32；氨基酸代谢：35-54；核苷酸代谢：33-34
input.pathway <- input.pathway[1:54]



pdf(file='scRNA_celltype_scMet_三大代谢2_dotplot.pdf',width = 10,height = 15)
DotPlot.metabolism(obj = res,
                   pathway = input.pathway, phenotype = "Celltype", norm = "y")+
  labs(x=NULL,y=NULL) + 
  theme(axis.text = element_text(size = 10, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size = 1))+ #轴标签
  scale_color_gradientn(colours = c('#3B4992FF',"#8491B4FF",'#F39B7FFF','#BB0021FF')) #颜色
dev.off()



#将代谢评分添加到meta.data中
METABOLISM <- res@assays$METABOLISM$score
colnames(METABOLISM) <- rownames(res@meta.data)
identical(colnames(METABOLISM), rownames(res@meta.data))
res@meta.data <- cbind(res@meta.data,t(METABOLISM))
saveRDS(METABOLISM, file = "METABOLISM.rds")
res@assays$METABOLISM$score <- METABOLISM



library(tidydr)
library(ggrepel)
p <- FeaturePlot(res,features = "Glycerophospholipid metabolism",
                 cols = c("#C3DBD9", "#990000")) + 
  theme_dr() + theme(panel.grid=element_blank())
ggsave("scRNA_Glycerophospholipid metabolism_umap.pdf",width = 5.5,height = 5)

VlnPlot(res, pt.size = 0, cols = cors,
        features = "Arginine biosynthesis")+
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white")+
  labs(title = "Arginine biosynthesis", x="",y="")+
  theme(legend.position = 'None')
ggsave("scRNA_ Arginine biosynthesis2.pdf",width = 6,height = 4)


library(SCP)
p <- SCP::FeatureDimPlot(res,"Glycerophospholipid metabolism",reduction = "umap",theme_use="theme_blank")
ggsave("Glycerophospholipid metabolism_scp.pdf",p,width = 5.5,height = 5)


cors <- c('#d4de9c','#94c58f','#86c7b4','#9cd2ed','#a992c0',
          '#ea9994','#f2c396','#bb82b1')







## 代谢差异分析 ####
res[["METABOLISM"]] <- CreateAssayObject(counts = METABOLISM)
DefaultAssay(res) <- 'METABOLISM'
Idents(res)='Epi 6'
DEG_scMet = FindMarkers(object = res,
                        only.pos = F, #only.pos改为T则只输出高表达gene
                        ident.1 = 'Epi 6',ident.2 = c('Epi 1','Epi 2'),
                        test.use = "LR",
                        min.cells.group = 1,
                        min.cells.feature = 1,
                        min.pct = 0,
                        logfc.threshold = 0)
DEG_scMet$pathway <- rownames(DEG_scMet)

## Cohen's d函数
cohens_d <- function(x, y, conf.level = 0.95) {
  pooled_std <- sqrt(((length(x)-1) * var(x) + (length(y)-1) * var(y)) / (length(x) + length(y) - 2))
  d <- (mean(x) - mean(y)) / pooled_std
  
  df <- length(x) + length(y) - 2
  nc <- d * sqrt((length(x) * length(y)) / (length(x) + length(y)))
  ci_lower <- qt((1 - conf.level) / 2, df, ncp = nc) / sqrt((length(x) * length(y)) / (length(x) + length(y)))
  ci_upper <- qt(1 - (1 - conf.level) / 2, df, ncp = nc) / sqrt((length(x) * length(y)) / (length(x) + length(y)))
  return(list(d = d, ci_lower = ci_lower, ci_upper = ci_upper))
}
cells_MEL <- rownames(scRNA_Epi@meta.data)[scRNA_Epi$tumor_type == "Tumor Epi"]
cells_HC <- rownames(scRNA_Epi@meta.data)[scRNA_Epi$tumor_type == "Normal Epi"]
METABOLISM_MEL <- METABOLISM[,cells_MEL]
METABOLISM_HC <- METABOLISM[,cells_HC]
# Cohen's d 的值约为 0.2 表示小效应,约为 0.5 表示中等效应,约为 0.8 表示大效应
for (id in rownames(DEG_scMet)) {
  A <- as.numeric(METABOLISM_MEL[id, ])
  B <- as.numeric(METABOLISM_HC[id, ])
  c_d <- cohens_d(B, A)
  DEG_scMet[id, 'cohens_d'] <- c_d
  DEG_scMet[id, 'ci_lower'] <- c_d$ci_lower; DEG_scMet[id, 'ci_upper'] <- c_d$ci_upper
}

DEG_scMet_sig = subset(DEG_scMet, p_val_adj < 0.05 & abs(cohens_d) > 0.25)
DEG_scMet_sig$pathway <- factor(DEG_scMet_sig$pathway, 
                                levels = DEG_scMet_sig[order(DEG_scMet_sig$cohens_d), "pathway"])
write.csv(DEG_scMet_sig, 'Epi_METABOLISM_DEG.csv', row.names = F)

# 森林图
ggplot(DEG_scMet_sig, aes(x = cohens_d, y = reorder(pathway, cohens_d))) +
  geom_point(size = 1.5, color = "black") +  # 缩小圆形点，颜色设置为蓝色
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "#3777ac") +  # 水平误差线
  geom_vline(xintercept = 0, linetype = "dashed", color = "#1f2d6f") +  # 无效应线
  theme_bw() +  # 最小化主题
  labs(title = NULL, x = "Cohen's d (Effect Size)", y = NULL) +
  xlim(-1, 0.6) +  # 限制 x 轴范围
  theme(
    panel.grid = element_blank(),  # 隐藏网格线
    axis.text = element_text(size = 10),  # 设置坐标轴文本大小和颜色
    axis.title = element_text(size = 10)
  )
ggsave("Epi_tumor_normal_scMet_forest.pdf", width = 6.5, height = 5.5)







######### monocle3

#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
library(monocle)
library(tidyverse)
library(CellChat)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(foreach)
library(KernSmooth)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(pheatmap)
library(ggheatmap)
library(reshape2)
library(SummarizedExperiment)
library(ggpubr)
library(ggrepel)
library(SeuratWrappers)
library(SeuratData)
library(copykat)
library(limma)
library(GSEABase)
library(GSVA)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringi)	
library(GOplot)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(monocle3)
#monocle3 version 1.3.7


sco <- scRNA_epi



#查看meta.data
colnames(sco@meta.data)
#画图
p1 = DimPlot(sco, group.by = "Celltype", label = T)
ggsave("1.UMAP.pdf", plot = p1, width = 6.5, height = 5)

#获取表达矩阵
data <- GetAssayData(sco, assay = 'RNA', layer  = 'counts')
#细胞注释信息
cell_metadata <- sco@meta.data
#基因名
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
#创建一个新的 cell_data_set 对象
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
#预处理cds
cds <- preprocess_cds(cds, method = "PCA")
#降维
cds <- reduce_dimension(cds, reduction_method = "UMAP",preprocess_method = 'PCA')
#画图
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters", show_trajectory_graph = FALSE) + ggtitle('cds.umap')

#将seurat对象的UMAP导入
int.embed <- Embeddings(sco, reduction = "umap")
#排序
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP),]
#导入
cds@int_colData$reducedDims$UMAP <- int.embed
#画图
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters", show_trajectory_graph = FALSE) + ggtitle('seurat.umap')
p = p1|p2
ggsave("2.Reduction_Compare.pdf", plot = p, width = 10, height = 5)

#聚类分区，不同分区的细胞会进行单独的轨迹分析
cds <- cluster_cells(cds)
#画图
p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + ggtitle("partition")
ggsave("3.Cluster_Partition.pdf", plot = p1, width = 6, height = 5)

#构建细胞轨迹
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))
#画图
p = plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave("4.Trajectory.pdf", plot = p, width = 6, height = 5)

#发育轨迹（拟时序）排列细胞。
cds <- order_cells(cds)
#画图
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE) +
  
  ggsave("5.Trajectory_Pseudotime.pdf", plot = p, width = 5.5, height = 5)
#保存结果
saveRDS(cds, file = "6.cds.rds")

#寻找拟时轨迹差异基因
#根据低维嵌入和主图测试基因的差异表达
#Monocle3 引入了一种寻找此类基因的新方法，它借鉴了空间相关性分析中的一项强大技术--莫兰 I 检验。
#莫兰 I 是一种多向、多维空间自相关的测量方法。
#该统计量可以告诉您，轨迹上邻近位置的细胞对被测基因的表达水平是否相似（或不同）。
#虽然皮尔逊相关性和莫兰 I 的范围都在-1 到 1 之间，但对莫兰 I 的解释略有不同：
#+1表示邻近细胞的表达完全相似；0表示没有相关性；-1表示邻近细胞没有相关性。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
#导出
write.csv(Track_genes, "7.Track_genes_all.csv", row.names = F)
#保存
saveRDS(Track_genes, "8.Track_genes_all.rds")

#拟时基因热图
genes <- row.names(subset(Track_genes, morans_I > 0.25))
#选取合适的clusters
num_clusters = 3
#画图
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
mycol = rev(RColorBrewer::brewer.pal(11, "Spectral"))
p = ComplexHeatmap::Heatmap(
  pt.matrix, name = "z-score", show_row_names = T, show_column_names = F,
  col = circlize::colorRamp2(seq(from=-2,to=2,length=11), mycol),
  row_names_gp = gpar(fontsize = 6), row_title_rot= 0, km = num_clusters, 
  cluster_rows = TRUE, cluster_row_slices = FALSE, cluster_columns = FALSE,use_raster=F
)
pdf("9.PseudotimeGenes_heatmap.pdf", width = 8, height = 8)
plot(p)
dev.off()

#基因展示,选取morans_I排名前12的基因，可自定义
genes = rownames(top_n(Track_genes,n=12,morans_I))
#画图
p <- plot_genes_in_pseudotime(cds[genes,], color_cells_by="celltype_rename", min_expr=0.5, ncol = 3)
ggsave("10.Genes_Jitterplot.pdf", plot = p, width = 8, height = 6)
p <- plot_cells(cds, genes=genes, show_trajectory_graph=FALSE,
                label_cell_groups=FALSE,  label_leaves=FALSE)
p$facet$params$ncol <- 3
ggsave("11.Genes_Featureplot.pdf", plot = p, width = 12, height = 9)

#寻找共表达基因模块
genelist <- row.names(subset(Track_genes, morans_I > 0.1))
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 6)
table(gene_module$module)
write.csv(gene_module, "12.PseudotimeGenes_Module.csv", row.names = F)

#热图
cell_group <- tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$celltype_rename)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
ggsave("13.PseudotimeGenes_Module.pdf", plot = p, width = 8, height = 6.5)

#散点图
ME.list <- lapply(1:length(unique(gene_module$module)), function(i){subset(gene_module,module==i)$id})
names(ME.list) <- paste0("module", 1:length(unique(gene_module$module)))
sco <- AddModuleScore(sco, features = ME.list, name = "module")
p <- FeaturePlot(sco, features = paste0("module", 1:length(unique(gene_module$module))), ncol = 2)
p <- p + plot_layout()&scale_color_viridis_c(option = "C")
ggsave("14.PseudotimeGenes_ModuleScore1.pdf", plot = p, width = 10, height = 8)








###### Cell chat
setwd("~/triapotosis")

library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(stringr)
library(CellChat)



####筛选目标组
scRNA_as <- subset(af,Celltype == "Astrocytes")

##根据氮代谢高低分组
scRNA_as$tse_group <- ifelse(scRNA_as$`Selenium Metabolism and Selenoproteins`>median(scRNA_as$`Selenium Metabolism and Selenoproteins`), 
                             'High_Astrocytes', 'Low_Astrocytes')


####准备工作####
#读取文件
scRNA <- af


##合并上皮细胞与其他细胞
table(scRNA$Celltype)
library(dplyr)
scRNA$cellid <- rownames(scRNA@meta.data) #原来的seurat metadata中添加一列celllid
scRNA_Epi$cellid <- rownames(scRNA_Epi@meta.data) #亚群的metadata中添加一列celllid

scRNA_Epi$subcelltype <- scRNA_Epi$tri_group

metadata <- scRNA_Epi@meta.data[,c('cellid','subcelltype')]
metadata <- left_join(x = scRNA@meta.data, y = metadata, by = "cellid")

metadata$merge_celltype <- ifelse(is.na(metadata$subcelltype), 
                                  as.character(metadata$Celltype), as.character(metadata$subcelltype))
table(metadata$merge_celltype)
#形成新的一列，分别是原来seurat的细胞注释和亚群的分类注释
scRNAsub <- scRNA[,metadata$cellid]
scRNAsub$merge_celltype <- metadata$merge_celltype#将合并的注释添加到原来的seurat meta即可
Idents(scRNAsub) <- "merge_celltype"
table(Idents(scRNAsub))

#将细胞按组织类型拆分

tumordata<-scRNAsub



####二、tumor组织####
####细胞通讯计算####
library(CellChat)
meta =tumordata@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(tumordata@assays$RNA$data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "merge_celltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net <- subsetCommunication(cellchat)
write.csv(df.net,"Epi_Communication.csv")

####细胞通讯数量强度####
#aggregateNet
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
#查看细胞通讯的数量/权重矩阵
mat <- cellchat@net$weight
mat <- cellchat@net$count
#绘制细胞通讯图（Ligand-receptor）
pdf(file=paste0('Epi_to_other_netVisual_circle_count.pdf'),width = 5.5,height = 5.5)
netVisual_circle(mat, 
                 sources.use = c(6,7),
                 targets.use = c(1:5,8:11), 
                 arrow.size = 0.2,
                 remove.isolate = F, #去掉孤立点
                 weight.scale = T, edge.weight.max = max(mat))
dev.off()


####信号通路通讯分析####
cellchat@netP$pathways   #信号通路查看
pathways.show <- c('VEGF')   #以'MIF'信号通路展示为例
levels(cellchat@idents)   #查看细胞亚群及factor顺序
vertex.receiver = c(6,7)  #左侧列展示感兴趣的亚群
#层级图（Hierarchy plot）
pdf(file=paste0('Epi_VEGF_hierarchy.pdf'),width = 10,height = 7)
netVisual_aggregate(cellchat,                #左侧列展示感兴趣的亚群
                    layout = c('hierarchy'), #"circle", "hierarchy", "chord"
                    signaling = pathways.show,
                    vertex.receiver = vertex.receiver)
dev.off()

#计算配受体对在目标信号通路中的贡献条形图
pdf(file=paste0('Epi_VEGF_contribution.pdf'),width = 4,height = 3)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
#提取细胞对
pairLR.CXCL <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] #以贡献度top1的配受体对为例
pairLR.CXCL
LR.show
#Hierarchy plot
pdf(file=paste0('Epi_VEGF_hierarchy_2.pdf'),width = 10,height = 7)
netVisual_individual(cellchat,
                     layout = c('hierarchy'),
                     signaling = pathways.show, #目标信号通路
                     pairLR.use = LR.show, #目标配受体对
                     vertex.receiver = vertex.receiver) #感兴趣的细胞亚群
dev.off()


####配体-受体通讯情况####
levels(cellchat@idents) #查看有哪些细胞类型
cellchat@netP$pathways   #信号通路查看
#sources.use为发出信号的细胞，targets.use为接受信号的细胞
#c()里面的数字与levels(cellchat@idents)中细胞的位次对应
##指定信号通路绘制气泡图
pdf(file=paste0('Epi_to_others_netVisual_bubble.pdf'),width = 7,height = 6)
netVisual_bubble(cellchat,
                 sources.use = c(6,7),
                 targets.use = c(1:5,8:11), 
                 #thresh = 0.01,
                 signaling = cellchat@netP$pathways[1:10], #指定信号通路
                 remove.isolate = T)+
  scale_color_gradientn(colors=c("#3C5488FF", "#f0e2a3", "#BB0021FF")) +
  theme_bw() + 
  theme(axis.text=element_text(size=8, colour = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()

#指定配受体对绘制气泡图
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("ANNEXIN","CCL","CXCL","VEGF")) #确定在目标信号通路中有重要作用的配受体对
pairLR.use
netVisual_bubble(cellchat,
                 sources.use = c(6,7),
                 targets.use = c(1:5,8:11), 
                 #thresh = 0.01,
                 pairLR.use = pairLR.use,
                 remove.isolate = T)+
  scale_color_gradientn(colors=c("#3C5488FF", "#f0e2a3", "#BB0021FF")) +
  theme_bw() + 
  theme(axis.text=element_text(size=8, colour = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) + coord_flip()

#参与目标信号通路的基因在各细胞亚群的表达分布展示
pdf(file=paste0('GC_celltype_MIF_genes_vlnplot.pdf'),width = 8,height = 4)
plotGeneExpression(cellchat, signaling = 'MIF', type = 'violin',color.use = cors) #小提琴图
dev.off()


#信号通路发出接收强度
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf(file=paste0('EPi_cellchat_VEGF_signaling.pdf'),width = 10,height = 3.5)
netAnalysis_signalingRole_network(cellchat, signaling = 'VEGF', width = 16, height = 3, font.size = 10)
dev.off()

h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf(file=paste0('Epi_cellchat_signaling_heatmap.pdf'),width = 11,height = 8)
h1 + h2
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
pdf(file=paste0('Epi_cellchat_signaling.pdf'),width = 5,height = 4.5)
gg1
dev.off()




########## Enrichment
library(tidyverse)
library(clusterProfiler)
library(corrplot)
library(ggforce)
library(RColorBrewer)




gene_df<-bitr(
  marker_epi6_sig$gene,         #######富集所用基因
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb="org.Hs.eg.db")%>%
  mutate(ENTREZID=as.integer(ENTREZID))


colnames(gene_df)[1] <- 'gene'
gene_df <- merge(gene_df,corr_data,by = 'gene')
gene_pos<-dplyr::filter(gene_df,cor>0.5)
gene_neg<-dplyr::filter(gene_df,cor<0)%>% dplyr::filter(abs(cor)>0.5)

gene_pos <- gene_df


kegg_pos <-enrichKEGG(gene_pos$ENTREZID,organism='hsa',keyType='kegg')
hsa <- read.delim('https://rest.kegg.jp/list/pathway/hsa',sep="\t",
                  col.names=c('ID','Name'))%>%
  mutate(Name=str_replace(Name,' - Homo sapiens \\(human\\)',''))
gene_info<-read.delim('https://rest.kegg.jp/list/hsa',sep='\t')%>%
  dplyr::select(1,4)%>%
  rename_with(~c('hsa','info'))%>%
  mutate(symbol=str_extract(info,"(^[\\w-]+)[;,]",group=1))
pathway<-read.delim('https://rest.kegg.jp/link/pathway/hsa',sep='\t',
                    col.names=c('hsa','ID'))%>%
  mutate(ID=str_remove(ID,'path:'))%>%
  inner_join(hsa)%>%
  inner_join(gene_info)%>%
  dplyr::select(c('Name','symbol'))
colnames(kegg_category)[4] <- 'Description'

kegg_pos<-enricher(gene_pos$gene,TERM2GENE=pathway)%>%
  as.data.frame()%>%
  dplyr::select(-ID)%>%
  dplyr::inner_join(kegg_category,.)%>%
  rowwise()%>%
  mutate(GeneRatio=eval(parse(text=GeneRatio)),
         type=1,Change='Positive')





go_pos<-enrichGO(gene_pos$gene,OrgDb="org.Hs.eg.db",
                 keyType='SYMBOL',ont='ALL')


data<-as.data.frame(go_pos)%>%
  dplyr::select('Description','ONTOLOGY','GeneRatio','qvalue')%>%
  mutate(len=str_length(Description))%>%
  dplyr::filter(len<120)%>%
  rowwise()%>%
  mutate(GeneRatio=round(eval(parse(text=GeneRatio)),3)*100)%>%
  arrange(qvalue)%>%
  group_by(ONTOLOGY)%>%
  mutate(ID=1:n())%>%
  top_n(5,wt=-ID)%>%{
    kegg<-as.data.frame(kegg_pos)%>%
      dplyr::select('Description','GeneRatio','qvalue')%>%
      rowwise()%>%
      mutate(GeneRatio=round(eval(parse(text=GeneRatio)),3)*100,
             ONTOLOGY='KEGG')%>%
      arrange(qvalue)%>%
      group_by(ONTOLOGY)%>%
      mutate(ID=1:n())%>%
      top_n(5,wt=-ID)
    rbind(.,kegg)
  }%>%
  mutate(ONTOLOGY=factor(ONTOLOGY,levels=c("BP","CC","MF","KEGG")),
         Description=str_wrap(Description,width=60),
         Description=factor(Description,levels=rev(Description)))



color <- c("#ea9994","#F2c396","#9cd2ed","#94c58f")
p <- ggplot(data,aes(-log10(qvalue),ID))+
  geom_col(aes(y=Description,fill=ONTOLOGY),alpha=0.5,show.legend=F)+
  geom_line(aes(x=GeneRatio,y=ID,group=1),col='black',size = 1,orientation="y",show.legend = F)+
  geom_point(aes(x=GeneRatio,y=ID,fill = ONTOLOGY),size=4,color='black',shape = 21,show.legend=F)+
  facet_wrap("ONTOLOGY",scales="free", nrow = 2)+
  scale_fill_manual(values=color)+
  labs(x="GeneRatio(%) and -log10(FDR)",y=NULL)+
  scale_x_continuous(limits=c(0,max(-log10(data$qvalue))+1))+
  theme(strip.text = element_text(size = 12,face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12,color="black"),
        axis.title.x = element_text(size = 14,color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black",
                                    size = 0.8,linetype = 1)) +
  scale_y_discrete(labels=function(y) str_wrap(y, width=30))

ggsave("Epi6_GOKEGG.pdf",width = 8,height = 5)
