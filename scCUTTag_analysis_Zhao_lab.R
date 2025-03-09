library(ArchR)
ArchR::installExtraPackages()
addArchRGenome("hg38") # hg19, mm9, mm10
addArchRThreads(threads = 5)
set.seed(5)


#######################################################
###   step1:  create arrow file   #####################

### input file
inputFiles <- getInputFiles(paths = "./Data/") ## postfix: *.fragments.tsv.gz
ArrowFiles <- ArchR::createArrowFiles(inputFiles = inputFiles,
                                      sampleNames = names(inputFiles), 
                                      filterTSS = 0,
                                      filterFrags = 1000,
                                      minFrags = 400,
                                      maxFrags = 20000,
                                      addTileMat = TRUE, 
                                      addGeneScoreMat = TRUE,
                                      TileMatParams = list(tileSize=5000), force=T
)
### delete doublet
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1 
)

########################################################
###   step2: create ArchRProject   #####################

proj_1 <- ArchRProject(ArrowFiles = ArrowFiles,outputDirectory = "H3K27ac",copyArrows = TRUE)
getAvailableMatrices(proj_1) 
proj_1 <- filterDoublets(proj_1)

##########################################################################
###   step3: LSI dimensionality reduction analysis   #####################

proj_1 <- addIterativeLSI(ArchRProj = proj_1,useMatrix = "TileMatrix",name = "IterativeLSI",iterations = 1,force=TRUE,varFeatures = 25000, dimsToUse = 1:30)
proj_1 <- addClusters(input = proj_1,reducedDims = "IterativeLSI",method = "Seurat",name = "Clusters",maxClusters = 40,resolution = 1,force = TRUE)
proj_1 <- addUMAP(ArchRProj = proj_1, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, force = T)
# proj_1 <- addTSNE(ArchRProj = proj_1, reducedDims = "IterativeLSI",name = "TSNE", perplexity = 30, force = T)
plotEmbedding(ArchRProj = proj_1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
plotEmbedding(ArchRProj = proj_1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
# plotEmbedding(ArchRProj = proj_1, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
# plotEmbedding(ArchRProj = proj_1, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
saveRDS(proj_1,"./rds/proj_LSI_reduction_raw_Int.rds") ## iteration1 var25000
proj_1 <- readRDS("./rds/proj_LSI_reduction_raw_Int.rds")

### QC analysis
library(ggplot2)
library(ggthemes)
clusters_order <- names(inputFiles)
##plotFragmentSizes to visualized fragment size
plotFragmentSizes(ArchRProj = proj_1)
ggsave("./plot/fragment_size_distribution.pdf", width = 5, height = 4)
##TSS enrichment
plotGroups(ArchRProj = proj_1,groupBy = "Sample",colorBy = "cellColData",name = "TSSEnrichment",plotAs = "ridges")
ggsave("./plot/TSS_enrichment_distribution.pdf", width = 5, height = 4)
## log10(UMI)
df <- getCellColData(proj_1, select = c("Sample", "log10(nFrags)", "TSSEnrichment"))
df$Sample <- factor(df$Sample,levels=names(inputFiles))
df = as.data.frame(df)
ggplot(data=df,aes(x=Sample,y=df[,2])) + 
  geom_violin(aes(fill=Sample)) + 
  geom_boxplot(width=0.1,outlier.shape = NA) +
  ylim(2.5,4.5) +
  theme_few() + 
  NoLegend() + ylab("Unique reads\nper cell (log10)") + xlab("")
ggsave(plot=p4, filename = './plot/UMI_boxplot.pdf',width=4,height = 3)

## Remove clusters with fewer than 100 cells
idxPass <- which(proj_1$Clusters %in% c("C13", "C14", "C15", "C16", "C17", "C18", "C19", "C24", "C6")) ## iteration1 var25000
cellsPass <- proj_1$cellNames[idxPass]
proj_2 = proj_1[cellsPass, ]
proj_2 <- addIterativeLSI(ArchRProj = proj_2, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 3, force=TRUE,varFeatures = 25000, dimsToUse = 1:50)
proj_2 <- addClusters(input = proj_2, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 1, force = TRUE)
proj_2 <- addUMAP(ArchRProj = proj_2, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 40, minDist = 0.5, force = T)
plotEmbedding(ArchRProj = proj_2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
plotEmbedding(ArchRProj = proj_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
saveRDS(proj_2,"./rds/proj_LSI_reduction_filt_low_cluster_Int.rds")

###############################################################
###   step4: Harmony remove batch effect  #####################

proj_3 <- addHarmony(ArchRProj = proj_2,reducedDims = "IterativeLSI",name = "Harmony",groupBy = "Sample")
proj_3 <- addClusters(input = proj_3,reducedDims = "Harmony",method = "Seurat",name = "Harmony_Clusters",maxClusters = 40,resolution = 1,force = T)
proj_3 <- addUMAP(ArchRProj = proj_3, reducedDims = "Harmony",  name = "Harmony_UMAP", nNeighbors = 40, minDist = 0.5, force = T)
proj_3 <- addTSNE(ArchRProj = proj_3, reducedDims = "Harmony", name = "Harmony_TSNE", perplexity = 30, force = T)
plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "Sample", embedding = "Harmony_UMAP")
plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "Harmony_Clusters", embedding = "Harmony_UMAP")
saveRDS(proj_3,"./rds/proj_LSI_reduction_after_harmony_Int.rds")
proj_3 = readRDS("./rds/proj_LSI_reduction_after_harmony_Int.rds")

##################################################################
###   step5: cells annotation with scRNA-seq   ###################

proj_4 <- addImputeWeights(proj_3)
## 鉴定不同cluster的标记基因
markersGS <- getMarkerFeatures(
  ArchRProj = proj_4, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Harmony_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05")
markerDf = data.frame()
markerDf_top20 = data.frame()
for (cluster in names(markerList)){
  clusterDF = as.data.frame(markerList[[cluster]])
  if (nrow(clusterDF) > 0){
    clusterDF$cluster = cluster
    clusterDF = clusterDF[order(-clusterDF$Log2FC), ]
    clusterDF_top20 = clusterDF[1:20, ]
    markerDf = rbind(markerDf, clusterDF)
    markerDf_top20 = rbind(markerDf_top20, clusterDF_top20)
  }
}
write.table(markerDf, "./markers/GeneScoreMatrix_total_markers_info.txt", sep = "\t", row.names = F, quote = F)
write.table(markerDf_top20, "./markers/GeneScoreMatrix_top20_markers_info.txt", sep = "\t", row.names = F, quote = F)


library(Seurat) ## V4 version
rna <- readRDS("../GSE138794/rds/GSE138794_five_samples_final_with_harmony.rds")
Idents(rna) <- "celltype"
proj_4 <- addGeneIntegrationMatrix(
  ArchRProj = proj_4,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_ArchR",
  reducedDims = "Harmony", ## 
  seRNA = rna,
  # groupList = groupList, ## 配对样本时设置
  addToArrow = F,
  groupRNA = "celltype",
  # useImputation = F,
  nameCell = "predictedCell_ArchR",
  nameGroup = "predictedGroup_ArchR",
  nameScore = "predictedScore_ArchR"
)
getAvailableMatrices(proj_4)
plotEmbedding(ArchRProj = proj_4, colorBy = "cellColData", name = "predictedGroup_ArchR", embedding = "Harmony_UMAP")
plotEmbedding(ArchRProj = proj_4, colorBy = "cellColData", name = "Harmony_Clusters", embedding = "Harmony_UMAP")

markerlist = c("CD2", "CD3D", "CD3E", "CD3G", ## T cells
               "CD14", "AIF1", "FCER1G", "TYROBP", ## macrophages
               "TF", "MOG", "MAG", "PLP1", ## oligodendrocytes ("FCGR3B")
               "AQP4", "GFAP", "MLC1", "ALDH1L1",  ## astrocytes ("GFAP")
               "CLDN5", "VWF", "CD34", "ABCG2",'RAMP2',  ## endothelials
               "PTPRZ1", "SOX2", "OLIG2", ## tumor cell
               "CHI3L1", "DLL3", "PDGFRA", "SMOC1", "CD44") ## GSC

pp <- plotEmbedding(
  ArchRProj = proj_4, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerlist, 
  embedding = "Harmony_UMAP",
  imputeWeights = getImputeWeights(proj_4)
)

ggsave(plot=p6$TYROBP, filename = './plot/macrophages_GeneScoreMatrix_TYROBP_UMAP_reduction_with.pdf',width=6,height = 5)
ggsave(plot=p6$PTPRZ1, filename = './plot/tumor_cells_GeneScoreMatrix_OLIG2_UMAP_reduction_with.pdf',width=6,height = 5)
ggsave(plot=p6$AQP4, filename = './plot/AST_GeneScoreMatrix_AQP4_UMAP_reduction_with.pdf',width=6,height = 5)
ggsave(plot=p6$MAG, filename = './plot/oligodendrocyte_GeneScoreMatrix_MAG_UMAP_reduction_with.pdf',width=6,height = 5)
# ggsave(plot=p6$DLL3, filename = './plot/GSC_GeneScoreMatrix_DLL3_UMAP_reduction_with.pdf',width=6,height = 5)
# ggsave(plot=p6$PDGFRA, filename = './plot/GSC_GeneScoreMatrix_PDGFRA_UMAP_reduction_with.pdf',width=6,height = 5)
# ggsave(plot=p6$SMOC1, filename = './plot/GSC_GeneScoreMatrix_SMOC1_UMAP_reduction_with.pdf',width=6,height = 5)

new_type = getCellColData(proj_4, select = c("Harmony_Clusters", "predictedGroup_ArchR"))
new_type$cellType = new_type$predictedGroup_ArchR
new_type[new_type$Harmony_Clusters %in% c("C6", "C7"), "cellType"] = "endothelials"
new_type[new_type$Harmony_Clusters%in% c("C1", "C2", "C3", "C4", "C5"), "cellType"] = "macrophages"
new_type[new_type$Harmony_Clusters %in% c("C8", "C9", "C10"), "cellType"] = "oligodendrocytes"
new_type[new_type$Harmony_Clusters %in% c("C17"), "cellType"] = "astrocytes"
new_type[new_type$Harmony_Clusters %in% c("C11", "C12", "C13", "C14", "C15", "C16", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25"), "cellType"] = "tumor cells"


proj_4$Celltype = new_type$cellType
plotEmbedding(ArchRProj = proj_4, colorBy = "cellColData", name = "Celltype", embedding = "Harmony_UMAP")
plotEmbedding(ArchRProj = proj_4, colorBy = "cellColData", name = "predictedGroup_ArchR", embedding = "Harmony_UMAP")
ggsave(filename = "./plot/Harmony_UMAP_reduction_with_celltype_label.pdf", width = 5, height = 6)
saveRDS(proj_4,"./rds/proj_LSI_reduction_after_harmony_GeneScores_Annotations_Int.rds")
proj_4 <- readRDS("./rds/proj_LSI_reduction_after_harmony_GeneScores_Annotations_Int.rds")


#############################################################
###   step6: re-cluster for tumor cells ######################

GSCsample <- BiocGenerics::which(proj_4$Celltype == "tumor cells")
cellsSample <- proj_4$cellNames[GSCsample]
proj_GSC <- proj_4[cellsSample, ]
proj_GSC <- addImputeWeights(proj_GSC)
getAvailableMatrices(proj_GSC)
## LSI降维分析
proj_GSC <- addClusters(
  input = proj_GSC,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters",
  # maxClusters = 25,
  resolution = 0.4,force = T ## 0.4
)
proj_GSC <- addUMAP(
  ArchRProj = proj_GSC, 
  reducedDims = "Harmony", 
  name = "UMAP_GSC", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)
plotEmbedding(ArchRProj = proj_GSC, colorBy = "cellColData", name = "Harmony_Clusters", embedding = "UMAP_GSC")
plotEmbedding(ArchRProj = proj_GSC, colorBy = "cellColData", name = "Sample", embedding = "UMAP_GSC")
saveRDS(proj_GSC,"./rds/proj_GSC_cell_Int.rds")

#########################
# meta module analysis ##
library(msigdbr)
library(scrabble)
CellInfo=proj_GSC@cellColData
CellInfo$Harmony_Clusters = factor(CellInfo$Harmony_Clusters)
Counts <- getMatrixFromProject(proj_GSC, useMatrix = "GeneScoreMatrix")
counts = assay(Counts)
rownames(counts) = rowData(Counts)[,5]
counts = as.data.frame(counts)
#Read the Genelist file #from Neftel et al
Subtype.Markers=read.csv("../Final/SuvaGBMsubtype genelist.csv")
Subtype.Markers = Subtype.Markers[Subtype.Markers$gene %in% rownames(counts), ]
#get scores
subtypelist=Subtype.Markers%>% split(x = .$gene, f = .$module)
colCenter = function(m, by = 'mean') {
  m = as.matrix(m)
  if (by == 'mean')  by = T
  else if (by == 'median') by = matrixStats::colMedians(m)
  else stopifnot(is.numeric(by) & length(by) == ncol(m))
  scale(m, center = by, scale = F)
}
sc=scrabble::score(counts, groups=subtypelist, binmat = NULL, bins = NULL, controls = NULL, 
                   bin.control = F, center = T, nbin = 30, n = 100, replace = T)
write.csv(sc,"HumanGliomacellstates_GSC_cells.csv")
sc=read.csv("HumanGliomacellstates_GSC_cells.csv",row.names = 1)
#get coordinates
h=scrabble::hierarchy (sc, quadrants = c("AC","OPC","MES","NPC"), log.scale = T)
# make plots
##per cluster
cols = c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767", "#90D5E4", "#89C75F")
for( i in 1:11){
  CellInfo[CellInfo$Harmony_Clusters == paste0("C", i), "Clustercolor"] = cols[i]
}
Clustercolor <- CellInfo$Clustercolor
names(Clustercolor ) <- row.names(CellInfo)

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

#Make the big plot--> all Clusters
groups=CellInfo[,c("Clustercolor","Harmony_Clusters")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
matrix$Harmony_Clusters[is.na(matrix$Harmony_Clusters)] <- "Other"
row.names(matrix)=matrix$Row.names
x=matrix$Clustercolor
y=matrix$Harmony_Clusters
col=x[!duplicated(x)]
names(col)=y[!duplicated(y)]
matrix=matrix[,-1]
title="All GSC Clusters"
p0 <- ggplot(matrix, aes(x = X,
                         y =Y,color=factor(Harmony_Clusters)))+geom_point()+geom_point(data = subset(matrix, Harmony_Clusters !="Other"))+ 
  scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Harmony_Clusters")+theme(legend.position = "none")+ 
  theme(legend.text=element_text(size=15),legend.title =element_text(size=15,face="bold") )+ 
  theme(axis.title.x = element_text(hjust = 0, vjust=-2, colour="black",size=10,face="bold"))+ 
  theme(axis.title.y = element_text(hjust = 0, vjust=4, colour="black",size=10,face="bold"))+  
  scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
  theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  geom_hline(yintercept=0, color = "black", size=0.5)+
  geom_vline(xintercept=0, color = "black", size=0.5)+
  annotate("rect", xmin = xmin, xmax = xmin+1.6, ymin = ymax-0.6, ymax = ymax, fill= "black")  + 
  annotate("text",x = xmin+0.8, y = ymax-0.3,label = "Mes-Like",color="white",fontface="bold",size=4)+ 
  annotate("rect", xmin = xmax-1.6, xmax = xmax, ymin = ymax-0.6, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-0.8, y = ymax-0.3,label = "NPC-Like",color="white",fontface="bold",size=4)+
  annotate("rect", xmin = xmin, xmax = xmin+1.6, ymin = ymin+0.6, ymax = ymin, fill= "black")  + 
  annotate("text",x = xmin+0.8, y = ymin+0.3,label = "AC-Like",color="white",fontface="bold",size=4)+ 
  annotate("rect", xmin = xmax-1.6, xmax =xmax, ymin = ymin+0.6, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-0.8, y = ymin+0.3,label = "OPC-Like",color="white",fontface="bold",size=4)  
P0=ggplotGrob(p0)
##make the small plots--> one per Cluster
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:11) {
  ClusterMD=CellInfo[CellInfo$Harmony_Clusters==paste0("C", i),]
  groups=ClusterMD[,c("Clustercolor","Harmony_Clusters")]
  title=paste0("C", i)
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
  matrix$Harmony_Clusters=as.character(matrix$Harmony_Clusters)
  matrix$Harmony_Clusters[is.na(matrix$Harmony_Clusters)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Clustercolor
  y=matrix$Harmony_Clusters
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  P=ggplot(matrix, aes(x = X,y =Y, color=factor(Harmony_Clusters)))+
    geom_point()+geom_point(data = subset(matrix, Harmony_Clusters !="Other"), alpha=0.5)+ 
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("text",x = xmin+0.8, y = ymax-0.3,label = "Mes-Like",color="black",fontface="bold",size=4)+ 
    annotate("text",x = xmax-0.8, y = ymax-0.3,label = "NPC-Like",color="black",fontface="bold",size=4)+
    annotate("text",x = xmin+0.8, y = ymin+0.3,label = "AC-Like",color="black",fontface="bold",size=4)+ 
    annotate("text",x = xmax-0.8, y = ymin+0.3,label = "OPC-Like",color="black",fontface="bold",size=4)  
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)

pdf(file ="./plot/HumanGliomacellstates_subtype_butterfly_by_GSC_Cluster_transferent.pdf", height = 12, width =18)
grid.arrange(grobs=Final, widths = c(1,1,1,1),layout_matrix = rbind(c(1, 2, 3, 4),c(5,6,7,8),c(9,10,11,12))) 
dev.off()

##########################
# markers GO enrichment ##

run_enrichment = function(UP,pro='test', wd = "./"){
  # library(KEGG.db)
  library(clusterProfiler)
  library(DOSE)
  library(stringr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(enrichplot)
  gene.up <- bitr(UP, fromType = "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)
  gene_up=unique(gene.up$ENTREZID)
  ###############   UP genes
  bp.up = enrichGO(gene=gene_up, OrgDb = org.Hs.eg.db,ont = "BP",readable=TRUE,
                   pvalueCutoff = 0.1,pAdjustMethod = "BH", qvalueCutoff = 0.1)
  cc.up = enrichGO(gene=gene_up, OrgDb = org.Hs.eg.db,ont = "CC",readable=TRUE,
                   pvalueCutoff = 0.1,pAdjustMethod = "BH", qvalueCutoff = 0.1)
  mf.up = enrichGO(gene=gene_up, OrgDb = org.Hs.eg.db,ont = "MF",readable=TRUE,
                   pvalueCutoff = 0.1,pAdjustMethod = "BH", qvalueCutoff = 0.1)
  k.up <- enrichKEGG(gene=gene_up, organism = "hsa", pvalueCutoff = 0.1,
                     pAdjustMethod = "BH", qvalueCutoff = 0.1, minGSSize = 1,use_internal_data =F)
  k.up=DOSE::setReadable(k.up, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.table(bp.up@result,paste0(wd, "/", pro,'_genes_BP_pathway.txt'), row.names = F, quote = F, sep = "\t")
  write.table(cc.up@result,paste0(wd, "/", pro,'_genes_CC_pathway.txt'), row.names = F, quote = F, sep = "\t")
  write.table(mf.up@result,paste0(wd, "/", pro,'_genes_MF_pathway.txt'), row.names = F, quote = F, sep = "\t")
  write.table(k.up@result,paste0(wd, "/", pro,'_genes_KEGG_pathway.txt'), row.names = F, quote = F, sep = "\t")
}

## marker heatmap
markersGS_GSC <- getMarkerFeatures(
  ArchRProj = proj_GSC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Cluster_merge",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS_GSC, cutOff = "FDR <= 0.01")
for (cluster in paste0("C", 1:11)){
  cluster_markers = markerList[[cluster]][, "name"]
  run_enrichment(cluster_markers, pro=cluster, wd="./GSC_cluster/enrichment")
}
pathway_for_plot = read.table("./GSC_cluster/enrichment/top10_pathway_for plot.txt", header = T, sep = "\t")
pathway_for_plot$ES = sapply(1:nrow(pathway_for_plot), function(x) 
{gr = as.numeric(strsplit(pathway_for_plot[x, "GeneRatio"], "/")[[1]]); gr1 = gr[1]; gr2 = gr[2]; 
br = as.numeric(strsplit(pathway_for_plot[x, "BgRatio"], "/")[[1]]); br1 = br[1]; br2 = br[2];
(gr1 * br2)/(gr2 * br1)})
ggplot(pathway_for_plot, aes(x=ES, y=Description))+
  geom_point(shape=21, aes(size=Count, fill=-log10(qvalue)))+
  scale_fill_gradient(low="#329DBE", high="#EE6F52") +
  facet_grid(cluster ~ ., scales = "free_y") + 
  # scale_y_discrete(limits=rev(shCtrl_Runx1_specific_top10$Description)) +
  labs(y="",title="enrichment pathways")+
  theme_bw() +
  theme(axis.text = element_text(colour = "#000000"),
        axis.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave("./GSC_cluster/enrichment/merged_top_BP_pathway.pdf", width = 6, height = 6)

###############################
## sample & cells distribution

metaInfo = data.frame()
for (cluster in paste0("C", 1:11)){
  temp = as.data.frame(table(GSCmetaData[GSCmetaData$Cluster_merge==cluster, "Sample"]))
  temp$cluster = cluster
  metaInfo = rbind(metaInfo, temp)
}
library(plyr)
metaInfo.ce = ddply(metaInfo, "cluster", transform, percet = Freq/sum(Freq) * 100)
ggplot(metaInfo.ce, aes(x=cluster, y=percet, fill=Var1)) +
  geom_bar(stat="identity") +
  labs(y="Percentage (%)", x = "") +
  scale_x_discrete(limits = paste0("C", 1:11))+
  scale_fill_brewer(palette = "Set3") +
  # scale_fill_manual(values = c("#7FC97F", "#FDC086", "#D58A91"))+
  theme_bw() +
  theme(axis.text = element_text(colour = "#000000"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave("plot/GSC_cluster_sample_contribution.pdf", width = 5, height = 5)

cellNum = as.data.frame(table(CellInfo$Cluster_merge))
cellNum$Var1 = factor(cellNum$Var1, levels = paste0("C", 1:11))
library(ggbreak)
ggplot(cellNum, aes(x=Freq, y=Var1, fill=Var1)) +
  geom_bar(stat="identity") +
  labs(y="", x = "No. of cells") +
  scale_y_discrete(limits = paste0("C", 1:11))+
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(colour = "#000000"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        legend.position = "none")
ggsave("plot/GSC_cluster_cell_number.pdf", width = 5, height = 3)

##########################################
# extract GSC cluster for SE analysis ####

GSCmetaData = proj_GSC@cellColData
for(sample in unique(GSCmetaData$Sample)){
  sampleCell = GSCmetaData[GSCmetaData$Sample ==sample, ]
  sampleDF = data.frame(cluster = sampleCell$Cluster_merge, barcodes = rownames(sampleCell))
  write.table(sampleDF, paste0("./GSC_cluster/", sample, "_GSC_cluster_barcode.csv"), col.names = F, row.names = F, quote = F, sep = ",")
}


#########################################################
################### track plot ##########################
fragments.ls <- lapply(names(inputFiles),function(x){
    fragments <- inputFiles[[x]]
    fragments <- rtracklayer::import(fragments,format = 'bed')
    fragments$label = paste0(x, "#", fragments$name)
    fragments
})
fragments <- Reduce(function(x,y) c(x,y), fragments.ls)
fragments <- sort(fragments)
saveRDS(fragments, "./rds/fragments.rds")
fragments = readRDS("./rds/fragments.rds")
chrom.sizes <- read.table(url('http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'),sep="\t",stringsAsFactors = FALSE)
chrom.sizes <- chrom.sizes[1:24,]
  
exportBW <- function(object,cluster,fragments){
  # if(class(object) == "Seurat"){
  #   cells <- rownames(object@meta.data[object@active.ident == cluster,])
  # }
  cells <- rownames(object[object$Cluster_merge==cluster,])
  new_read <- GRanges(seqnames = chrom.sizes[,1], 
                      ranges =IRanges(start = as.numeric(chrom.sizes[,2]),
                                      width=1),
                      name = rep("in_silico_extra_read",dim(chrom.sizes)[1]),
                      score = rep(0,dim(chrom.sizes)[1])
  )
  fragments.x <- fragments$label %in% cells
  fragments.x <- fragments[fragments.x]
  fragments.x <- c(fragments.x,new_read)
  
  coverage.x <- coverage(fragments.x)
  coverage.x <- coverage.x/(length(fragments.x)/1e6)
  rtracklayer::export.bw(object = coverage.x,paste0("GSC_cluster/GSC_cluster_",cluster,".bw"))
}

lapply(unique(GSCmetaData$Cluster_merge),function(x){
  exportBW(GSCmetaData,x,fragments)
})
  

library(Gviz)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(Signac)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bw.files <- paste0("GSC_cluster/GSC_cluster_C",1:11, ".bw")
bw.ls <- lapply(bw.files, function(x) {x <- rtracklayer::import(x);x})

marker.regions <- c(
  StringToGRanges("chr3:126094754-126107427",sep = c(":","-")),             # SLC41A3 10 ## c1
  StringToGRanges("chr19:35507259-35519487",sep = c(":","-")),             # 8  ## c1
  StringToGRanges("chr19:34587729-34602601",sep = c(":","-")),             # 8  ## c1
  StringToGRanges("chr7:55562748-55581074",sep = c(":","-")),             # VOPP1 50  ## c2
  StringToGRanges("chr15:43249712-43259506",sep = c(":","-")),              # TGM5 , 7  ## c3
  StringToGRanges("chr20:63694743-63703953",sep = c(":","-")),              # 10  ## c3
  StringToGRanges("chr19:6765022-6776617",sep=c(":","-")),                 # SH2D3A, 3.5 ## c4
  StringToGRanges("chr6:31554641-31568310",sep=c(":","-")),                    # 5  ## c4
  StringToGRanges("chr11:128691489-128699840",sep=c(":","-")),             # 4.5  ## c4
  StringToGRanges("chr5:1877522-1887735",sep = c(":","-")),               # IRX4, 2.5  ## c5
  StringToGRanges("chr16:49277995-49285916",sep = c(":","-")),            # 8.5  ## c5
  StringToGRanges("chr8:24950673-24960587",sep = c(":","-")),             # 4.5  ## c5
  StringToGRanges("chr7:54405632-54415596",sep = c(":","-")),               # , 5  ## c6
  StringToGRanges("chr7:55257822-55266880",sep = c(":","-")),               # 6    ## c6
  StringToGRanges("chr7:54765150-54771624",sep = c(":","-")),               # 6    ## c6
  StringToGRanges("chr7:26530890-26540594",sep=c(":","-")),                # , 2.5  ## c7
  StringToGRanges("chr13:21580749-21592374",sep=c(":","-")),               # 3.5    ## c7
  StringToGRanges("chr12:119003101-119014717",sep=c(":","-")),             # 5      ## c7
  StringToGRanges("chr1:153625279-153630654",sep = c(":","-")),             # S100A1, 26   ## c8
  StringToGRanges("chr19:35283715-35301285",sep = c(":","-")),              # 4.5   ## c8
  StringToGRanges("chr17:81121030-81137552",sep = c(":","-")),              # 31    ## c8
  StringToGRanges("chr1:151836995-151842965",sep = c(":","-")),             # 34    ## c8
  StringToGRanges("chr3:134673829-134687470",sep = c(":","-")),             # 5     ## c8
  StringToGRanges("chr17:46727759-46738547",sep = c(":","-")),            # NSF, 13  ## c9 
  StringToGRanges("chr1:160163856-160175832",sep = c(":","-")),           # 13   ## c9  
  StringToGRanges("chr7:50331064-50344640",sep = c(":","-")),             # 25   ## c9 
  StringToGRanges("chr8:127139607-127171142",sep = c(":","-")),            # 30
  StringToGRanges("chr8:127005501-127021484",sep = c(":","-")),            # 12
  StringToGRanges("chr8:124446632-124457740",sep = c(":","-")),            # 51
  StringToGRanges("chr1:203338889-203353383",sep = c(":","-")),            # 71
  StringToGRanges("chr1:203298114-203311157",sep = c(":","-")),            # 222
  StringToGRanges("chr1:204185296-204200267",sep = c(":","-"))           # 116
)             

to.plot.ls <- lapply(bw.ls,function(x){
  x <- subsetByOverlaps(x,sort(marker.regions))
  x
})

# Colors defined by CTcolors
CTcolors = c("#D51F26", "#90D5E4", "#89C75F", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767")

p <- lapply(seq(marker.regions),function(y){
  ylimits = c(0,ceiling(max(subsetByOverlaps(do.call('c',bw.ls),marker.regions[y])$score) /10) * 10)
  bw.tracks <- lapply(seq(bw.ls),function(x){
    track <- DataTrack(range = bw.ls[[x]],chromosome = as.character(seqnames(marker.regions[y])),
                       from = start(marker.regions[y]), to = end(marker.regions[y]),
                       type="polygon",showTitle=T,col.axis="black",
                       background.title="transparent",col.baseline="black",
                       col.mountain="transparent",fill.mountain=c(CTcolors[x],CTcolors[x]),ylim=ylimits[y],yTicksAt = ylimits[y])
    track
  })
  
  myAxisTrack <- GenomeAxisTrack(col="black")
  grtrack <- GeneRegionTrack(txdb,showTitle=T,col.axis="black",background.title="transparent",col.baseline="black",
                             chromosome = as.character(seqnames(marker.regions[y])),from = start(marker.regions[y]), to = end(marker.regions[y]),
                             stacking = 'squish',col='black',fill='red')
  return(c(myAxisTrack,bw.tracks,grtrack))
})

pdf(file = 'plot/GSC_pseudobulk_cluster_H3K27.pdf',width = 2,height=15)
sapply(seq(marker.regions),function(x){
  plotTracks(p[[x]],chromosome=as.character(seqnames(marker.regions[x])),from = start(marker.regions[x]), to = end(marker.regions[x]), ylim=range(0,ylimits[x]),
             scale=5000,min.width=1,min.distance=1,mergeGroups= TRUE,lwd=0.5,col.line='black',window=2000,sizes = c(1,rep(1,length(bw.ls)),3))
  
})
dev.off()


############################################################
################### bw plot ################################

fetchFragmentsByRegion <- function(fragments, region) { 
  if(length(region) > 1){
    stop("Only one region can be used")
  }
  if(class(region)[1] == "GRanges"){
    seqlevelsStyle(region) <- "UCSC"
    chr    = as.character(seqnames(region))
    from  = as.numeric(start(region))
    to    = as.numeric(end(region))
  }
  if(class(region)[1] == "character"){
    region <- strsplit(region,":|-")[[1]]
    
    chr    = as.character(region[1])
    from  = as.numeric(region[2])
    to    = as.numeric(region[3])
  }
  # Extract fragments overalpping the region
  fragments.interesting <- fragments[seqnames(fragments) == chr]
  fragments.interesting <- fragments.interesting[start(fragments.interesting) > from]
  fragments.interesting <- fragments.interesting[end(fragments.interesting) < to]
  return(fragments.interesting)
}

filterCells <- function(object,fragments,fraction){
  # Make fragments per cell table
  fragments.summary <- sort(table(c(fragments$name,rownames(object))),decreasing = TRUE) - 1
  
  # Split cells by cluster - return fraction of cells per cluster
  fragments.summary.ls <-  lapply(levels(object$Clusters),function(x){
    barcodes <- rownames(object[object$Clusters==x, ])
    fragments.summary.x <- fragments.summary[names(fragments.summary) %in% barcodes]
    fragments.summary.x <- sort(fragments.summary.x,decreasing = TRUE)
    head(x = fragments.summary.x, n = length(barcodes) * fraction)
  })
  
  # unlist and merge cella names
  cells <- unlist(lapply(fragments.summary.ls,names))
  return(cells)
}

fragmentsToMatrix <- function(fragments,window){
  bins <- tile(range(fragments),width=window)[[1]]
  fragments.ls <- split(fragments,as.factor(fragments$name))
  seqlevels(bins,pruning.mode = "coarse") <- names(coverage(fragments.ls[[1]]))
  fragments.ls.coverage <- lapply(fragments.ls,function(x){
    x <- coverage(x)
    x.bins <- binnedAverage(bins = bins,numvar = x, "coverage")
    as.numeric(x.bins$coverage)
  })
  fragments.matrix <- t(do.call("rbind",fragments.ls.coverage))
  rownames(fragments.matrix) <- start(bins)
  return(fragments.matrix)
}

plotHeatmapBW <- function(object,
                          fragmentsGrange,
                          interestingRegion,
                          window = 5000,
                          annotFilter,
                          fraction = 0.1,
                          clusters_order = NULL,
                          cells = NULL,...
) {
  library(pheatmap)
  if(is.null(clusters_order)){clusters_order  <- levels(object)}
  if(is.null(cells)){cells.to.plot            <- sample(rownames(object),fraction * nrow(object))}
  
  fragments.interesting <- fragmentsGrange  
  
  
  fragments.interesting.matrix <- fragmentsToMatrix(fragments    = fragments.interesting,
                                                    window       = window)
  
  
  fragments.interesting.matrix <- fragments.interesting.matrix[,colnames(fragments.interesting.matrix) %in% cells]
  
  
  dummy.mat <- matrix(data = 0, nrow = dim(fragments.interesting.matrix)[1],ncol = length(cells.to.plot[!cells.to.plot %in% colnames(fragments.interesting.matrix)]))
  colnames(dummy.mat) <- cells.to.plot[!cells.to.plot %in% colnames(fragments.interesting.matrix)]
  rownames(dummy.mat) <- rownames(fragments.interesting.matrix)
  
  fragments.interesting.matrix <- cbind(fragments.interesting.matrix,dummy.mat)
  
  fragments.interesting.matrix.binary <- apply(fragments.interesting.matrix,2,as.logical)
  fragments.interesting.matrix.binary <- Matrix::Matrix(apply(fragments.interesting.matrix.binary,2,as.numeric),sparse = TRUE)
  
  rownames(fragments.interesting.matrix) <- as.numeric(rownames(fragments.interesting.matrix))
  rownames(fragments.interesting.matrix.binary) <- as.numeric(rownames(fragments.interesting.matrix))
  
  
  ####################################################
  round_coordinates = 250000
  
  
  
  clusters_annotations <- sample(factor(as.matrix(object[colnames(fragments.interesting.matrix.binary),])[,"Clusters"],levels = clusters_order))
  
  annotation_row       <- data.frame(clusters_annotations[order(clusters_annotations)])
  colnames(annotation_row) <- "Cluster"
  
  xlabels  <- as.numeric(rownames(fragments.interesting.matrix))
  xlabels[xlabels %% round_coordinates != 0] <- ""
  
  p1 <- pheatmap(t(fragments.interesting.matrix.binary)[rownames(annotation_row),],
                 cluster_rows=FALSE,
                 cluster_cols=FALSE,
                 show_rownames=FALSE,
                 show_colnames=TRUE,
                 labels_col = xlabels,
                 annotation_row = annotation_row,
                 annotation_legend= TRUE,
                 annotation_names_row = FALSE,
                 color = colorRampPalette(c("white","black"))(2),
                 gaps_row= cumsum(table(annotation_row$Cluster)[clusters_order]),
                 border_color = NA)
  
  return(p1)
}

fragments <- readRDS("./rds/fragments.rds")
fragments$name <- fragments$label

interesting_region <- Signac::StringToGRanges("chr17:48510001-48609832",sep=c(":","-")) ## HOXB3
window_bin = 250

library(AnnotationFilter)
flt <- AnnotationFilterList(TxBiotypeFilter("protein_coding"),
                            GRangesFilter(interesting_region),
                            logicOp = c("&"))
fragments.interesting <- fetchFragmentsByRegion(fragments = fragments,
                                                region    = interesting_region)


proj_GSC <- readRDS("./rds/proj_GSC_cell_Int.rds")
sampleInfo <- as.data.frame(getCellColData(proj_GSC, select = c("Celltype", "Clusters")))
sampleInfo$Clusters <- factor(sampleInfo$Clusters, levels = paste0("C", 1:11))

fragments.interesting.GSC <- fragments.interesting[fragments.interesting$label %in% rownames(sampleInfo)]
fraction              <- length(unique(fragments.interesting.GSC$label)) / nrow(sampleInfo)
cells.to.plot         <- filterCells(object = sampleInfo,
                                     fragments = fragments.interesting.GSC,
                                     fraction=fraction)
p <- plotHeatmapBW(object = sampleInfo,
                   fragmentsGrange = fragments.interesting.GSC,
                   cells = cells.to.plot,
                   interestingRegion = interesting_region,
                   window = window_bin,
                   annotFilter = flt,
                   fraction = fraction,
                   clusters_order=levels(sampleInfo$Clusters))

pdf("plot/GSC_HOXB3_pseudobulk_cluster_H3K27_heatmap.pdf",height = 5, width = 5)
p
dev.off()

