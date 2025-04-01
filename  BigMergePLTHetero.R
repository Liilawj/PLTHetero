session_info <- sessionInfo()
print(session_info)
capture.output(session_info, file = "session_info.txt")##uploaded in GitHub.
library(Matrix)
library(vctrs)
library(matrixStats) 
library(Rcpp)
library(dplyr) 
library(dbplyr) 
library(BiocGenerics) 
library(scales)
library(multtest)
library(Seurat)
library(patchwork)
library(R.utils)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(igraph)
library(openxlsx)
library(monocle)
library(rstudioapi)
current_doc <- getSourceEditorContext()
content <- current_doc$contents
setDocumentContents(paste(new_content, collapse = "\n"), current_doc$id)
obj <- readRDS("~/E-MTAB-10026/obj.rds")
CovintegratedAnnorerun2SR <- readRDS("~/Covint.filtered.clusteredAllmoreres.anno.rds")
BigMerge <- merge(CovintegratedAnnorerun2SR, y = obj, add.cell.ids = c("GSE149689","EMTAB10026"), project = "BigMerge")
BigMPLT = BigMerge[,BigMerge$predicted.celltypewithnew %in% c("Platelet")]
BigMPLT@meta.data
dim(BigMPLT)
set.seed(1)
BigMPLT <- NormalizeData(BigMPLT, normalization.method = "LogNormalize", scale.factor = 10000)
BigMPLT <- FindVariableFeatures(BigMPLT, selection.method = "vst", nfeatures = 2000)
BigMPLT <- ScaleData(BigMPLT,features = rownames(BigMPLT))
BigMPLT <- RunPCA(BigMPLT, npcs = 50, verbose = FALSE)
DimHeatmap(BigMPLT,dims=1:10,cells=420,balanced=TRUE)
DimHeatmap(BigMPLT,dims=11:20,cells=420,balanced=TRUE)
VizDimLoadings(BigMPLT,dims = 1:11, reduction="pca")
DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",dims = 10:11,
        order="M")
DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",
        dims = 1:2,shuffle = T)+
  NoLegend()|
  DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",
          dims = 3:4,shuffle = T)+
  NoLegend()|
  DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",
          dims = 5:6,shuffle = T)+
  NoLegend()|
  DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",
          dims = 7:8,shuffle = T)+
  NoLegend()|
  DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",
          dims = 8:9,shuffle = T)+
  NoLegend()|
  DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",
          dims = 9:10,shuffle = T)+
  NoLegend()|
  DimPlot(BigMPLT,reduction="pca",group.by = "MSH.idents",
          dims = 11:12,shuffle = T)
print(PLTCov@reductions$pca)
set.seed(1)
BigMPLT <- FindNeighbors(
  BigMPLT, dims = 1:20)
BigMPLT <- FindClusters(
  BigMPLT, resolution = 0.5)
BigMPLT <- RunUMAP(
  BigMPLT, dims = 1:20,return.model = T)
BigMPLT <- RunTSNE(
  BigMPLT, dims = 1:20)
DimPlot(BigMPLT, reduction = "tsne",group.by = "orig.ident",
        label = T,pt.size=0.5)
DimPlot(BigMPLT, reduction = "tsne",group.by = "sample_id",
        label = F,pt.size=0.5)+
  NoLegend()
DimPlot(BigMPLT, reduction = "tsne",group.by = "MSH.idents",
        label = F,pt.size=0.5)
BigMPLTintegrated <- BigMPLT
library(harmony)
BigMPLTintegrated=BigMPLTintegrated %>% RunHarmony("studyID", plot_convergence = TRUE,lambda=0.5)
DimHeatmap(BigMPLTintegrated,dims=1:10,cells=500,balanced=TRUE)
DimHeatmap(BigMPLTintegrated,dims=11:20,cells=500,balanced=TRUE)
DimHeatmap(BigMPLTintegrated,dims=21:30,cells=500,balanced=TRUE)
DimHeatmap(BigMPLTintegrated,dims=31:40,cells=500,balanced=TRUE)
set.seed(1)
BigMPLTintegrated <- BigMPLTintegrated %>% 
  RunUMAP(reduction = "harmony", dims = 1:20,return.model = T) %>% 
  RunTSNE(reduction = "harmony", dims = 1:20)%>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
library(DoubletFinder)
sweep.res.list_kidney <- 
  paramSweep_v3(BigMPLTintegrated,PC=1:10)
sweep.stats_kidney <- 
  summarizeSweep(sweep.res.list_kidney,GT=F)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
pk_bcmvn <- bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)] %>% as.character() %>% as.numeric()
annotations <- BigMPLTintegrated$predicted.celltypePLTres0.7
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*length(BigMPLTintegrated$predicted.celltypePLTres0.7))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
BigMPLTintegrated <- doubletFinder_v3(
  BigMPLTintegrated,PCs=1:10,
  pN=0.25,pK=pk_bcmvn,
  nExp=nExp_poi,reuse.pANN = F)
BigMPLTintegrated@meta.data[,"DF_hi.lo"] <- 
  BigMPLTintegrated@meta.data$DF.classification
DimPlot(BigMPLTintegrated,reduction="umap",group.by="DF.classifications_0.25_0.01_1118")
table(BigMPLTintegrated$DF.classifications_0.25_0.01_1118)
Doublet <- BigMPLTintegrated@meta.data
BigMPLTintegratedsingle <- BigMPLTintegrated %>% select(rownames(Doublet %>% filter(DF.classifications_0.25_0.01_1118 != "Doublet")))
DimPlot(BigMPLTintegrated, reduction = "umap", group.by = "ID",
        pt.size=0.0001,label = F)+
  NoLegend()+
  theme(axis.line = element_blank(),axis.ticks = element_blank(),
        axis.text = element_blank())
DimPlot(BigMPLTintegrated, reduction = "tsne", group.by = "orig.ident",
        pt.size=0.0001,label = F)+
  NoLegend()+
  theme(axis.line = element_blank(),axis.ticks = element_blank(),
        axis.text = element_blank())
DimPlot(BigMPLTintegrated, reduction = "umap",
        group.by = "RNA_snn_res.0.5",split.by = "MSH.idents",
        label = TRUE,pt.size=0.5,repel = TRUE)+
  theme(axis.line = element_blank(),axis.ticks = element_blank(),
        axis.text = element_blank())
DimPlot(BigMPLTintegrated, reduction = "umap",
        group.by = "MSH.idents",
        label = TRUE,pt.size=0.5,repel = TRUE)+
  theme(axis.line = element_blank(),axis.ticks = element_blank(),
        axis.text = element_blank())
DimPlot(BigMPLTintegrated, reduction = "tsne", 
        group.by = "RNA_snn_res.0.5",split.by = "MSH.idents",
        pt.size=0.5,label = TRUE,repel = TRUE)+
  theme(axis.line = element_blank(),axis.ticks = element_blank(),
        axis.text = element_blank())
set.seed(1)
PLTCov <- FindClusters(
  PLTCov, resolution = 0.7)
PLTCov <- RunUMAP(
  PLTCov, dims = 1:20,return.model = T)
PLTCov <- FindClusters(
  PLTCov, resolution = 1) 
PLTCov <- RunUMAP(
  PLTCov, dims = 1:20,return.model = T)
PLTCov <- FindClusters(
  PLTCov, resolution = 1.5)
PLTCov <- RunUMAP(
  PLTCov, dims = 1:20,return.model = T)
PLTCov <- FindClusters(
  PLTCov, resolution = 0.3)
PLTCov <- RunUMAP(
  PLTCov, dims = 1:20,return.model = T)
PLTCov <- FindClusters(
  PLTCov, resolution = 0.17)
PLTCov <- RunUMAP(
  PLTCov, dims = 1:20,return.model = T)
PLTCov <- FindClusters(
  PLTCov, resolution = 0.05)
PLTCov <- RunUMAP(
  PLTCov, dims = 1:20,return.model = T)
library(clustree)
p_tree <- clustree(PLTCov@meta.data, prefix = "RNA_snn_res.")
DimPlot(PLTCov, reduction = "tsne",group.by = "RNA_snn_res.0.5",
        label = T,pt.size=0.5)
DimPlot(PLTCov, reduction = "tsne",group.by = "RNA_snn_res.0.5",
        split.by = "MSH.idents",
        label = T,pt.size=0.5)
DimPlot(PLTCov, reduction = "umap",group.by = "RNA_snn_res.0.5",
        label = T,pt.size=0.5)
DimPlot(PLTCov, reduction = "umap",group.by = "RNA_snn_res.0.5",
        split.by = "MSH.idents",
        label = T,pt.size=0.5)
PLT.cluster@meta.data$Worst_Clinical_Status <- 
  factor(PLT.cluster$Worst_Clinical_Status,
         levels =c("Death","Severe","Critical ","Moderate",
                   "Mild","Asymptomatic","Non-covid","LPS",
                   "Healthy","nan"))
DimPlot(PLT.cluster, reduction = "ref.umap",cols = cols,
        group.by = "predicted.celltype", label = TRUE,
        split.by = "Worst_Clinical_Status",
        label.size = 3, repel = TRUE) + 
  NoLegend() + 
  ggtitle("Query transferred labels")
DimPlot(PLTCov, reduction = "tsne",
        group.by = "RNA_snn_res.0.7",label = T,pt.size=0.3)
DimPlot(PLTCov, reduction = "tsne",
        group.by = "RNA_snn_res.1",label = T,pt.size=0.3)
DimPlot(PLTCov, reduction = "tsne",
        group.by = "RNA_snn_res.1.5",label = T,pt.size=0.3)
DimPlot(PLTCov, reduction = "tsne",
        group.by = "RNA_snn_res.1.5",
        split.by = "cond.idents",label = T,pt.size=1)
DimPlot(PLTCov, reduction = "tsne",
        group.by = "RNA_snn_res.1",
        split.by = "cond.idents",label = T,pt.size=1)

Idents(BigMPLTintegrated) <- BigMPLTintegrated$PLT.idents
markers.f3 <- FindAllMarkers(BigMPLTintegrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
DT::datatable(markers.f3)
library(dplyr)
PLTtop10 <- markers.f3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
PLTtop100 <- markers.f3 %>% group_by(cluster)%>% top_n(n = 100, wt = avg_log2FC) 
PLTtop20 <- markers.f3 %>% group_by(cluster)%>% top_n(n = 20, wt = avg_log2FC) 
BiocManager::install("ComplexHeatmap",force=T)
BiocManager::install("Mfuzz",force=T)
devtools::install_github("junjunlab/jjAnno")
library(ComplexHeatmap)
library(Mfuzz)
library(jjAnno)
library(yulab.utils)
library(monocle)
BiocManager::install("monocle")
install_zip("~/Downloads/ClusterGVis-main.zip")
install_zip("~/Downloads/ClusterGVis-0.1.1.zip")
library(ClusterGVis)
library(org.Hs.eg.db)
library(ggplot2)
devtools::install_github("junjunlab/scRNAtoolVis")
devtools::install_github("sajuukLyu/ggunchull")
library(scRNAtoolVis)
markers.f3 <- read.table("~/PLT07unsupmarkers.txt",sep='\t')
PLTtop100 <- markers.f3 %>% group_by(cluster)%>% top_n(n = 100, wt = avg_log2FC)
PLTtop10 <- markers.f3 %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)
Idents(BigMPLTintegrated) <- BigMPLTintegrated$predicted.celltypePLTres0.7
input.data3 <- prepareDataFromscRNA(
  BigMPLTintegrated,
  group.by = "predicted.celltypePLTres0.7",
  diffData = PLTtop100,showAverage = T)
input.data3modi2 <- prepareDataFromscRNA(
  BigMPLTintegrated,
  diffData = PLTtop100,showAverage = T)
visCluster(input.data3modi2,plot.type = "both")
visCluster(input.data3,plot.type = "both")
visCluster(input.data3modifinal,plot.type = "both",
           line.side = "left",column_names_rot=45,
           markGenes = PLTtop10$gene,
           cluster.order=c("3","7","2","4","9","6","1","5","8"))
write.csv(input.data3modi2$long.res,"long.res.csv")
write.csv(input.data3modi2$wide.res,"wide.res.csv")
input.data3modifinal <- input.data3modi2
input.data3modifinal$wide.res <- read.csv("wide.res.csv")
str(input.data3)
dim(input.data3modi2$long.res)
dim(input.data3modi2$wide.res)
levels(as.factor(input.data3$long.res$cluster))
levels(as.factor(input.data3$long.res$cluster_name))
levels(as.factor(input.data3$long.res$cell_type)) 
input.data3modi <- input.data3
input.data3modi$long.res$cluster_name <- input.data3$long.res$cell_type
input.data3modi$long.res$cluster <- input.data3$long.res$cell_type
input.data3modi$wide.res$cluster <- input.data3$long.res$cell_type[1:674]
enrichGO <- enrichCluster(input.data3modifinal, 
                          OrgDb = org.Hs.eg.db,
                          type = "BP", organism = "hsa",
                          pvalueCutoff = 0.05,topn = 10, seed = 1)
levels(as.factor(input.data3modifinal$wide.res$cluster))
levels(as.factor(input.data3modifinal$long.res$cluster))
levels(as.factor(input.data3modifinal$long.res$cluster_name))
enrichKEGG <- enrichCluster(input.data3modifinal, 
                            OrgDb = org.Hs.eg.db,
                          type = "KEGG", organism = "hsa",
                          pvalueCutoff = 0.05,topn = 10, seed = 1)
enrichGO$ratio <- as.numeric(enrichGO$ratio)
enrichKEGG <- enrichGO
str(enrichGO)
write.csv(enrichGO,"enrichGO.csv")
length(levels(as.factor(enrichGO$group)))
levels(as.factor(enrichGO$group))
palette = c("1Grays","2Light Grays","3Blues2","4Blues3",
            "5Purples2","6Purples3","7Reds2","8Reds3","9Greens2","10Greens3")
palette = c("Blues2","Reds2","Light Grays","Blues3",
            "Greens2","Purples3","Grays","Purples2","Reds3","Greens3")
lapply(seq_along(unique(enrichGO$group)), function(x){
  tmp <- enrichGO |> dplyr::filter(group == unique(enrichGO$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  tmp$Description <- factor(tmp$Description,levels = tmp$Description)
  p <-
    ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue),y = Description,fill = -log10(pvalue)),
             width = 0.75) +
    geom_line(aes(x = log10(ratio),y = as.numeric(Description)),color = "grey50") +
    geom_point(aes(x = log10(ratio),y = Description),size = 3,color = "orange") +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~.,name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  tmp.kg <- enrichKEGG |> dplyr::filter(group == unique(enrichKEGG$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  tmp.kg$Description <- factor(tmp.kg$Description,levels = tmp.kg$Description)
  pk <-
    ggplot(tmp.kg) +
    geom_segment(aes(x = 0,xend = -log10(pvalue),y = Description,yend = Description),
                 lty = "dashed",linewidth = 0.75) +
    geom_point(aes(x = -log10(pvalue),y = Description,color = -log10(pvalue)),size = 5) +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    colorspace::scale_color_binned_sequential(palette = palette[x]) +
    ylab("") + xlab("-log10(pvalue)")
  cb <- cowplot::plot_grid(plotlist = list(p,pk))
  return(cb)
}) -> gglist
length(gglist)
gglist[[4]]
enrichGO$group
names(gglist) <- paste("C",c(1,2,3,4,5,6,7,9,8),sep = "")
names(gglist) <- paste("C",c(3,7,2,4,9,6,1,5,8),sep = "")
pdf('~/heatmap.pdf',height = 15,width = 25,onefile = F)
visCluster(object = input.data3modifinal,plot.type = "both",
           line.side = "left",column_names_rot = 45,
           markGenes = PLTtop10$gene,
           cluster.order=c("3","7","2","4","9","6","1","5","8"),
           ggplot.panel.arg = c(5,0.5,32,"grey90",NA),
           markGenes=F,
           show_row_dend = F,
           annoTerm.data = enrichGO)
dev.off()
visCluster(input.data3modifinal,plot.type = "both",
           line.side = "left",column_names_rot=45,
           markGenes = PLTtop10$gene,
           cluster.order=c("3","7","2","4","9","6","1","5","8"))
pdf('bigHmap.pdf',height = 20,width = 22,onefile = F)
visCluster(object = input.data3modifinal,
           plot.type = "both",
           line.side = "left",
           column_names_rot = 45,
           markGenes = PLTtop10$gene,
           cluster.order=c("3","7","2","4","9","6","1","5","8"),
           ggplot.panel.arg = c(5,2,32,"grey90",NA),
           gglist = gglist,
           show_row_dend = F,
           annoTerm.data = enrichGO)
dev.off()
str(gglist)
saveRDS(gglist,"~/LargeHmap/gglist.rds")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
markers <- read.table("~/BigMPLTres07markers.txt",sep='\t')
all.top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
top100 <- all.top100[!duplicated(all.top100$gene),]
select_genes_all = split(top100$gene,top100$cluster)
markers <- read.table("~/BigMPLTres07markers.txt",sep='\t')
markers2 <- markers[!duplicated(markers$gene),]
select_genes_all = split(markers2$gene,markers2$cluster)
clu9.sig.genes <- select_genes_all$"9"
clu9.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu9.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu9.Ssig_entrezId)
clu9.go_Ssiggenes <- enrichGO(gene = clu9.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu9.go_Ssiggenes_results <- as.data.frame(clu9.go_Ssiggenes@result)
colnames(clu9.go_Ssiggenes_results)
clu9.go_Ssiggenes_results <- clu9.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu9.go_Ssiggenes_results$geneID <- 
  str_replace(clu9.go_Ssiggenes_results$geneID,"/",",")
names(clu9.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu9.go_Ssiggenes_results)
heatplot(clu9.go_Ssiggenes)
clu9.go_Ssiggenes@result$Description
clu9.go_Ssiggenes@result$pvalue
cnetplot(clu9.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu9.go_Ssiggenes_lowP <- 
  subset(clu9.go_Ssiggenes,clu9.go_Ssiggenes@result$pvalue < 0.05)
dim(clu9.go_Ssiggenes_lowP)
clu9.go_Ssiggenes_lowP
str(input.data3modifinal)
saveRDS(input.data3modifinal,"~/LargeHmap/input.data3modifinal.rds")
clu9.sig.genes <- c("EIF1", "MT1X", "HBG2", "GPX4", "IFITM3", "FTL", "NAP1L1", "IFI27", "OAZ1", "AQP10", "TIMP1", "CLEC1B", "GP9", "FTH1", "CMTM5", "HLA-A", "HIST1H2AC", "H3F3A", "GSTO1", "CTSD", "SOD2", "SVIP", "SLC39A3", "TLN1", "RAB32", "EIF2AK1", "STOM", "NAA38", "MT-CYB", "FAH", "HLA-C", "YWHAH", "MT-CO2", "TGFB1", "EMP3", "ETFA", "MT-ATP6", "GAS2L1", "MT-CO1", "MAP1A", "PARK7", "MT-ND4", "LY6G6F", "IFI27L2", "IFI6", "RGS18", "C9orf16", "C12orf76", "PF4", "SRP14", "CCS", "ITGA2B", "STXBP2", "MT-CO3", "ARPC2", "MGLL", "GTF3C6", "CAP1", "FAXDC2", "NDUFS5", "TRAPPC5", "DAP", "ASAH1", "UQCR11", "PNP", "MT-ND2", "RAB4A", "SSR4", "SMIM5", "TMSB4X", "RHOC", "ISCA1", "PPBP", "AP2S1", "ARHGAP6", "C12orf75", "CAPN2", "TMBIM6", "MTURN", "SEC14L1", "DERA", "MT-ND3", "GNAS", "EGFL7", "CD63", "PCMT1", "GLA", "TERF2IP", "ISCU", "VPS28", "MT-ND1", "PLEKHO1", "MYL9", "ACTG1", "ARPC1B", "ACTB", "SH3BGRL3", "MYL6", "SERF2")
clu9.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu9.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu9.Ssig_entrezId)
clu9.go_Ssiggenes <- enrichGO(gene = clu9.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu9.go_Ssiggenes_lowP <- 
  subset(clu9.go_Ssiggenes,clu9.go_Ssiggenes@result$pvalue < 0.05)
dim(clu9.go_Ssiggenes_lowP)
clu9.go_Ssiggenes_lowP
write.csv(clu9.go_Ssiggenes_lowP,"~/LargeHmap/clu9.go_Ssiggenes_lowP.csv")
enrichGO
str(enrichGO)
write.csv(enrichGO,"~/LargeHmap/enrichGO.csv")
enrichGO <- read.csv("~/LargeHmap/enrichGO.csv")
function (object = NULL, diffData = NULL, showAverage = TRUE, 
    cells = NULL, group.by = "ident", assays = "RNA", slot = "data", 
    scale.data = TRUE, cluster.order = NULL, keep.uniqGene = TRUE, 
    sep = "_") 
{
    markerGene <- unique(diffData$gene)
    if (showAverage == TRUE) {
        mean_gene_exp <- Seurat::AverageExpression(object, features = markerGene, 
            group.by = group.by, assays = assays, slot = slot) %>% 
            data.frame() %>% as.matrix()
        name1 <- gsub(pattern = paste0(assays, ".", sep = ""), 
            replacement = "", colnames(mean_gene_exp))
        colnames(mean_gene_exp) <- gsub(pattern = "\\.", replacement = " ", 
            name1)
        colnames(mean_gene_exp) <- levels(Seurat::Idents(object))
        if (scale.data == TRUE) {
            mean_gene_exp <- t(scale(t(mean_gene_exp)))
        }
        if (!is.null(cluster.order)) {
            mean_gene_exp <- mean_gene_exp[, cluster.order]
        }
        geneMode = "average"
    }
    else {
        cell.order <- data.frame(cell.id = names(Seurat::Idents(object)), 
            cell.ident = Seurat::Idents(object))
        if (is.null(cluster.order)) {
            cell.order$cell.ident <- factor(cell.order$cell.ident, 
                levels = levels(Seurat::Idents(object)))
        }
        else {
            cell.order$cell.ident <- factor(cell.order$cell.ident, 
                levels = cluster.order)
        }
        cell.order <- cell.order[order(cell.order$cell.ident), 
            ]
        getassy <- Seurat::GetAssayData(object = object, slot = slot)[features = markerGene, 
            cells = NULL, drop = FALSE] %>% as.matrix()
        id.order <- match(cell.order$cell.id, colnames(getassy))
        getassy <- getassy[, id.order]
        colnames(getassy) <- paste(colnames(getassy), cell.order$cell.ident, 
            sep = "|")
        mean_gene_exp <- getassy
        if (scale.data == TRUE) {
            mean_gene_exp <- t(scale(t(mean_gene_exp)))
        }
        geneMode = "all"
    }
    merMat <- data.frame(mean_gene_exp, check.names = FALSE) %>% 
        tibble::rownames_to_column(., var = "gene")
    cinfo.gene <- diffData[, c("cluster", "gene")]
    cn <- unique(cinfo.gene$cluster)
    wide.res <- purrr::map_df(seq_along(cn), function(x) {
        tmp <- cinfo.gene[which(cinfo.gene$cluster == cn[x]), 
            ]
        tmp2 <- merMat[which(merMat$gene %in% tmp$gene), ] %>% 
            dplyr::mutate(cluster = as.character(x))
        return(tmp2)
    })
    if (keep.uniqGene == TRUE) {
        wide.res <- wide.res %>% dplyr::distinct(., gene, .keep_all = TRUE)
        geneType = paste("unique", sep, sep = "|")
    }
    else {
        wide.res <- wide.res %>% dplyr::mutate(., gene = make.unique(gene, 
            sep = sep))
        geneType = paste("nounique", sep, sep = "|")
    }
    df <- reshape2::melt(wide.res, id.vars = c("cluster", "gene"), 
        variable.name = "cell_type", value.name = "norm_value")
    df$cluster_name <- paste("cluster ", df$cluster, sep = "")
    if (showAverage == FALSE) {
        df$cell_type <- sapply(strsplit(as.character(df$cell_type), 
            split = "\\|"), "[", 2)
    }
    cl.info <- data.frame(table(wide.res$cluster)) %>% dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>% 
        dplyr::arrange(Var1)
    id <- unique(df$cluster_name)
    df <- purrr::map_df(seq_along(id), function(x) {
        tmp <- df %>% dplyr::filter(cluster_name == id[x])
        tmp %>% dplyr::mutate(cluster_name = paste(cluster_name, 
            " (", cl.info$Freq[x], ")", sep = ""))
    })
    df$cluster_name <- factor(df$cluster_name, levels = paste("cluster ", 
        cl.info$Var1, " (", cl.info$Freq, ")", sep = ""))
    return(list(wide.res = wide.res, long.res = df, type = "scRNAdata", 
        geneMode = geneMode, geneType = geneType))
}
source("getGOTerm.R")
GO_DATA <- get_GO_data("org.Hs.eg.db","ALL","SYMBOL")
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
Idents(BigMPLTintegrated) <- BigMPLTintegrated$predicted.celltypePLTres0.7
markers <- FindAllMarkers(BigMPLTintegrated, 
                          test.use = "wilcox",
                          only.pos = TRUE, 
                          min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"~/BigMPLTres07markers.txt",sep='\t')
markers <- read.table("~/BigMPLTres07markers.txt",sep='\t')
DT::datatable(markers)
library(dplyr)
all.top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
all.top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
all.top10
write.table(all.top10,"PLTusingConventionalPC1-11/unSupmarkers/PLTCovMSH.cluster_markersTop10res1.5.MAST.txt",sep='\t')
DoHeatmap(BigMPLTintegrated, features = all.top10$gene) + NoLegend()
all.top10 <- markers %>% group_by(cluster) %>% top_n(10,avg_log2FC)
top10 <- all.top10[!duplicated(all.top10$gene),]
top100 <- all.top100[!duplicated(all.top100$gene),]
top100$cluster <- factor(top100$cluster,levels=c("9","8","7","6","5",
                                               "4","3","2","1","0"))
BigMPLTintegrated$predicted.celltypePLTres0.7 <- factor(
  BigMPLTintegrated$predicted.celltypePLTres0.7,
  levels=c("9","8","7","6","5","4","3","2","1","0"))
select_genes_all = split(top10$gene,top10$cluster)
select_genes_all = split(top100$gene,top100$cluster)
str(select_genes_all)
saveRDS(select_genes_all,"select_genes_all.rds")
library(ggplot2)
DefaultAssay(BigMPLTintegrated) <- "RNA"
DotPlot(object = BigMPLTintegrated, 
        features=select_genes_all,
        assay="RNA",dot.scale = 10,
        cols=c("lightgrey", "darkred"))+
  theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5,size=5),axis.text.y=element_text(vjust=0.5,hjust=0.5,size=7))
DoHeatmap(BigMPLTintegrated, 
        features=top10$gene,
        assay="RNA")+
  theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5,size=5),
        axis.text.y=element_text(vjust=0.5,hjust=0.5,size=7))
top10
library(ggsci)
library(scales)
mycol = c(pal_npg()(10))
# mycol = c(pal_npg()(8),"
#          "
#          "
DimPlot(BigMPLTintegrated,group.by = "predicted.celltypePLTres0.7",
        label = T,cols = mycol,pt.size = 0.5) +
  theme(legend.position = "none",
        plot.title = element_blank()) +
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  theme(axis.line = element_blank(),axis.title = element_blank())
DimPlot(BigMPLTintegrated,group.by = "PLT.idents",
        label = T,cols = mycol,pt.size = 0.5) +
  theme(legend.position = "none",
        plot.title = element_blank()) +
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  theme(axis.line = element_blank(),axis.title = element_blank())
select_genes_all <- readRDS("select_genes_all.rds")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu2.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu2.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu2.Ssig_entrezId)
clu2.go_Ssiggenes <- enrichGO(gene = clu2.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu2.go_Ssiggenes_results <- as.data.frame(clu2.go_Ssiggenes@result)
colnames(clu2.go_Ssiggenes_results)
clu2.go_Ssiggenes_results <- clu2.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu2.go_Ssiggenes_results$geneID <- 
  str_replace(clu2.go_Ssiggenes_results$geneID,"/",",")
names(clu2.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu2.go_Ssiggenes_results)
heatplot(clu2.go_Ssiggenes)
clu2.go_Ssiggenes@result$Description
clu2.go_Ssiggenes@result$pvalue
cnetplot(clu2.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu2.go_Ssiggenes_lowP <- 
  subset(clu2.go_Ssiggenes,clu2.go_Ssiggenes@result$pvalue < 0.05)
dim(clu2.go_Ssiggenes_lowP)
str(clu2.go_Ssiggenes_lowP)
write.table(clu2.go_Ssiggenes_lowP,"~/LargeHmap/clu2.go_Ssiggenes_lowP_top100.txt")
clu2.go_Ssiggenes_lowP
clu2.kk <- enrichKEGG(clu2.Ssig_entrezId, 
                       organism = "hsa",
                      keyType = "kegg", pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)
dim(clu2.kk)
barplot(clu2.kk)
enrich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$pvalue,decreasing=FALSE)[1:40],]
  data4plot$BgRatio <-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][2])
    })
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p <- p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(pvalue))) +
    scale_size(range=c(2,13))
  
  pr <- pbubble + 
  scale_colour_gradient(low="",
    labs(color=expression(-log[10](pvalue)),size="observed.gene.count",
         x="Richfactor",y="term.description",title="Enrichment Process")
  
  pr <- pr + theme_bw()
  pr
}
p1 <- enrich2plot(clu2.kk@result)+
  scale_color_gradient(low="blue",high="red")+
  theme_bw()
p2 <- enrich2plot(clu2.go_Ssiggenes@result) +
  scale_color_gradient(low="blue",high="red")+
  theme_bw()
clu2.go_SsiggenesOrdered <- clu2.go_Ssiggenes[order(clu2.go_Ssiggenes$GeneRatio),]
clu2.go_SsiggenesOrdered$ONTOLOGY <- factor(clu2.go_SsiggenesOrdered$ONTOLOGY,levels=clu2.go_SsiggenesOrdered$ONTOLOGY)
enrich2plot(clu2.go_SsiggenesOrdered)
library(patchwork)
p2 + p1 + plot_layout(ncol=2)
p1 <- enrich2plot(clu45.kk@result)+
  scale_color_gradient(low="blue",high="red")+
  theme_bw()+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size=15),axis.text.y = element_text(size=27),
        axis.title = element_text(size=27))+
  theme(legend.title = element_text(size=27),legend.text = element_text(size=20))+
  scale_y_discrete(position="right")
p2 <- enrich2plot(clu45.go_Ssiggenes@result) +
  scale_color_gradient(low="blue",high="red")+
  theme_bw()+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size=15),axis.text.y = element_text(size=27),
        axis.title = element_text(size=27))+
  theme(legend.position = "none")
p2
p1
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu0.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu0.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu0.Ssig_entrezId)
clu0.go_Ssiggenes <- enrichGO(gene = clu0.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu0.go_Ssiggenes_results <- as.data.frame(clu0.go_Ssiggenes@result)
colnames(clu0.go_Ssiggenes_results)
clu0.go_Ssiggenes_results <- clu0.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu0.go_Ssiggenes_results$geneID <- 
  str_replace(clu0.go_Ssiggenes_results$geneID,"/",",")
names(clu0.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu0.go_Ssiggenes_results)
heatplot(clu0.go_Ssiggenes)
clu0.go_Ssiggenes@result$Description
clu0.go_Ssiggenes@result$pvalue
cnetplot(clu0.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu0.go_Ssiggenes_lowP <- 
  subset(clu0.go_Ssiggenes,clu0.go_Ssiggenes@result$pvalue < 0.05)
dim(clu0.go_Ssiggenes_lowP)
write.table(clu0.go_Ssiggenes_lowP,"~/LargeHmap/clu0.go_Ssiggenes_lowP_top100.txt")
clu0.go_Ssiggenes_lowP
enrich2plot(clu0.go_Ssiggenes@result) +
  scale_color_gradient(low="blue",high="red")+
  theme_bw()
enrich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$pvalue,decreasing=FALSE)[1:40],]
  data4plot$BgRatio <-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][2])
    })
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p <- p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(pvalue))) +
    scale_size(range=c(2,13))
  
  pr <- pbubble + scale_colour_gradient(low="",
    labs(color=expression(-log[10](pvalue)),size="observed.gene.count",
         x="Richfactor",y="term.description",title="Enrichment Process")
  
  pr <- pr + theme_bw()
  pr
}
p1 <- enrich2plot(clu2.kk@result)+
  scale_color_gradient(low="blue",high="red")+
  theme_bw()
p2 <- enrich2plot(clu2.go_Ssiggenes@result) +
  scale_color_gradient(low="blue",high="red")+
  theme_bw()
clu2.go_SsiggenesOrdered <- clu2.go_Ssiggenes[order(clu2.go_Ssiggenes$GeneRatio),]
clu2.go_SsiggenesOrdered$ONTOLOGY <- factor(clu2.go_SsiggenesOrdered$ONTOLOGY,levels=clu2.go_SsiggenesOrdered$ONTOLOGY)
enrich2plot(clu2.go_SsiggenesOrdered)
library(patchwork)
p2 + p1 + plot_layout(ncol=2)
p1 <- enrich2plot(clu45.kk@result)+
  scale_color_gradient(low="blue",high="red")+
  theme_bw()+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size=15),axis.text.y = element_text(size=27),
        axis.title = element_text(size=27))+
  theme(legend.title = element_text(size=27),legend.text = element_text(size=20))+
  scale_y_discrete(position="right")
p2 <- enrich2plot(clu45.go_Ssiggenes@result) +
  scale_color_gradient(low="blue",high="red")+
  theme_bw()+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size=15),axis.text.y = element_text(size=27),
        axis.title = element_text(size=27))+
  theme(legend.position = "none")
p2
p1
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu1.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu1.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu1.Ssig_entrezId)
clu1.go_Ssiggenes <- enrichGO(gene = clu1.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu1.go_Ssiggenes_results <- as.data.frame(clu1.go_Ssiggenes@result)
colnames(clu1.go_Ssiggenes_results)
clu1.go_Ssiggenes_results <- clu1.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu1.go_Ssiggenes_results$geneID <- 
  str_replace(clu1.go_Ssiggenes_results$geneID,"/",",")
names(clu1.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu1.go_Ssiggenes_results)
heatplot(clu1.go_Ssiggenes)
clu1.go_Ssiggenes@result$Description
clu1.go_Ssiggenes@result$pvalue
cnetplot(clu1.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu1.go_Ssiggenes_lowP <- 
  subset(clu1.go_Ssiggenes,clu1.go_Ssiggenes@result$pvalue < 0.05)
dim(clu1.go_Ssiggenes_lowP)
write.table(clu1.go_Ssiggenes_lowP,"~/LargeHmap/clu1.go_Ssiggenes_lowP_top100.txt")
clu1.go_Ssiggenes_lowP
enrich2plot(clu1.go_Ssiggenes@result) +
  scale_color_gradient(low="blue",high="red")+
  theme_bw()
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu3.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu3.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu3.Ssig_entrezId)
clu3.go_Ssiggenes <- enrichGO(gene = clu3.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu3.go_Ssiggenes_results <- as.data.frame(clu3.go_Ssiggenes@result)
colnames(clu3.go_Ssiggenes_results)
clu3.go_Ssiggenes_results <- clu3.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu3.go_Ssiggenes_results$geneID <- 
  str_replace(clu3.go_Ssiggenes_results$geneID,"/",",")
names(clu3.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu3.go_Ssiggenes_results)
heatplot(clu3.go_Ssiggenes)
clu3.go_Ssiggenes@result$Description
clu3.go_Ssiggenes@result$pvalue
cnetplot(clu3.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu3.go_Ssiggenes_lowP <- 
  subset(clu3.go_Ssiggenes,clu3.go_Ssiggenes@result$pvalue < 0.05)
dim(clu3.go_Ssiggenes_lowP)
write.table(clu3.go_Ssiggenes_lowP,"~/LargeHmap/clu3.go_Ssiggenes_lowP_top100.txt")
clu1.go_Ssiggenes_lowP
enrich2plot(clu3.go_Ssiggenes@result) +
  scale_color_gradient(low="blue",high="red")+
  theme_bw()
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu4.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu4.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu4.Ssig_entrezId)
clu4.go_Ssiggenes <- enrichGO(gene = clu4.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu4.go_Ssiggenes_results <- as.data.frame(clu4.go_Ssiggenes@result)
colnames(clu4.go_Ssiggenes_results)
clu4.go_Ssiggenes_results <- clu4.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu4.go_Ssiggenes_results$geneID <- 
  str_replace(clu4.go_Ssiggenes_results$geneID,"/",",")
names(clu4.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu4.go_Ssiggenes_results)
heatplot(clu4.go_Ssiggenes)
clu4.go_Ssiggenes@result$Description
clu4.go_Ssiggenes@result$pvalue
cnetplot(clu4.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu4.go_Ssiggenes_lowP <- 
  subset(clu4.go_Ssiggenes,clu4.go_Ssiggenes@result$pvalue < 0.05)
dim(clu4.go_Ssiggenes_lowP)
write.table(clu4.go_Ssiggenes_lowP,"~/LargeHmap/clu4.go_Ssiggenes_lowP_top100.txt")
clu4.go_Ssiggenes_lowP
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu8.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu8.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu8.Ssig_entrezId)
clu8.go_Ssiggenes <- enrichGO(gene = clu8.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu8.go_Ssiggenes_results <- as.data.frame(clu8.go_Ssiggenes@result)
colnames(clu8.go_Ssiggenes_results)
clu8.go_Ssiggenes_results <- clu8.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu8.go_Ssiggenes_results$geneID <- 
  str_replace(clu8.go_Ssiggenes_results$geneID,"/",",")
names(clu8.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu8.go_Ssiggenes_results)
heatplot(clu8.go_Ssiggenes)
clu8.go_Ssiggenes@result$Description
clu8.go_Ssiggenes@result$pvalue
cnetplot(clu8.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu8.go_Ssiggenes_lowP <- 
  subset(clu8.go_Ssiggenes,clu8.go_Ssiggenes@result$pvalue < 0.05)
dim(clu8.go_Ssiggenes_lowP)
write.table(clu8.go_Ssiggenes_lowP,"~/LargeHmap/clu8.go_Ssiggenes_lowP_top100.txt")
clu8.go_Ssiggenes_lowP
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu7.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu7.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu7.Ssig_entrezId)
clu7.go_Ssiggenes <- enrichGO(gene = clu7.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu7.go_Ssiggenes_results <- as.data.frame(clu7.go_Ssiggenes@result)
colnames(clu7.go_Ssiggenes_results)
clu7.go_Ssiggenes_results <- clu7.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu7.go_Ssiggenes_results$geneID <- 
  str_replace(clu7.go_Ssiggenes_results$geneID,"/",",")
names(clu7.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu7.go_Ssiggenes_results)
heatplot(clu7.go_Ssiggenes)
clu7.go_Ssiggenes@result$Description
clu7.go_Ssiggenes@result$pvalue
cnetplot(clu7.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu7.go_Ssiggenes_lowP <- 
  subset(clu7.go_Ssiggenes,clu7.go_Ssiggenes@result$pvalue < 0.05)
dim(clu7.go_Ssiggenes_lowP)
write.table(clu7.go_Ssiggenes_lowP,"~/LargeHmap/clu7.go_Ssiggenes_lowP_top100.txt")
clu7.go_Ssiggenes_lowP
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu5.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu5.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu5.Ssig_entrezId)
clu5.go_Ssiggenes <- enrichGO(gene = clu5.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu5.go_Ssiggenes_results <- as.data.frame(clu5.go_Ssiggenes@result)
colnames(clu5.go_Ssiggenes_results)
clu5.go_Ssiggenes_results <- clu5.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu5.go_Ssiggenes_results$geneID <- 
  str_replace(clu5.go_Ssiggenes_results$geneID,"/",",")
names(clu5.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu5.go_Ssiggenes_results)
heatplot(clu5.go_Ssiggenes)
clu5.go_Ssiggenes@result$Description
clu5.go_Ssiggenes@result$pvalue
cnetplot(clu5.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu5.go_Ssiggenes_lowP <- 
  subset(clu5.go_Ssiggenes,clu5.go_Ssiggenes@result$pvalue < 0.05)
dim(clu5.go_Ssiggenes_lowP)
write.table(clu5.go_Ssiggenes_lowP,"~/LargeHmap/clu5.go_Ssiggenes_lowP_top100.txt")
clu5.go_Ssiggenes_lowP
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu6.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu6.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu6.Ssig_entrezId)
clu6.go_Ssiggenes <- enrichGO(gene = clu6.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu6.go_Ssiggenes_results <- as.data.frame(clu6.go_Ssiggenes@result)
colnames(clu6.go_Ssiggenes_results)
clu6.go_Ssiggenes_results <- clu6.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu6.go_Ssiggenes_results$geneID <- 
  str_replace(clu6.go_Ssiggenes_results$geneID,"/",",")
names(clu6.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu6.go_Ssiggenes_results)
heatplot(clu6.go_Ssiggenes)
clu6.go_Ssiggenes@result$Description
clu6.go_Ssiggenes@result$pvalue
cnetplot(clu6.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu6.go_Ssiggenes_lowP <- 
  subset(clu6.go_Ssiggenes,clu6.go_Ssiggenes@result$pvalue < 0.05)
dim(clu6.go_Ssiggenes_lowP)
write.table(clu6.go_Ssiggenes_lowP,"~/LargeHmap/clu6.go_Ssiggenes_lowP_top100.txt")
clu6.go_Ssiggenes_lowP
library(dplyr)
markers <- read.table("~/BigMPLTres07markers.txt",sep='\t')
all.top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
top100 <- all.top100[!duplicated(all.top100$gene),]
select_genes_all = split(top100$gene,top100$cluster)
markers <- read.table("~/BigMPLTres07markers.txt",sep='\t')
markers2 <- markers[!duplicated(markers$gene),]
select_genes_all = split(markers2$gene,markers2$cluster)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
clu0.sig.genes <- select_genes_all$"0"
clu1.sig.genes <- select_genes_all$"1"
clu2.sig.genes <- select_genes_all$"2"
clu3.sig.genes <- select_genes_all$"3"
clu4.sig.genes <- select_genes_all$"4"
clu5.sig.genes <- select_genes_all$"5"
clu6.sig.genes <- select_genes_all$"6"
clu7.sig.genes <- select_genes_all$"7"
clu8.sig.genes <- select_genes_all$"8"
clu9.sig.genes <- select_genes_all$"9"
clu9.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                             keys=clu9.sig.genes,
                             keytype = "SYMBOL",
                             column="ENTREZID")
length(clu9.Ssig_entrezId)
clu9.go_Ssiggenes <- enrichGO(gene = clu9.Ssig_entrezId, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",ont="ALL",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5,readable = TRUE)
clu9.go_Ssiggenes_results <- as.data.frame(clu9.go_Ssiggenes@result)
colnames(clu9.go_Ssiggenes_results)
clu9.go_Ssiggenes_results <- clu9.go_Ssiggenes_results[,c(1,2,3,9,7)]
clu9.go_Ssiggenes_results$geneID <- 
  str_replace(clu9.go_Ssiggenes_results$geneID,"/",",")
names(clu9.go_Ssiggenes_results) <- c("Category","ID","term","Genes","adj_pval")
dim(clu9.go_Ssiggenes_results)
heatplot(clu9.go_Ssiggenes)
clu9.go_Ssiggenes@result$Description
clu9.go_Ssiggenes@result$pvalue
cnetplot(clu9.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
clu9.go_Ssiggenes_lowP <- 
  subset(clu9.go_Ssiggenes,clu9.go_Ssiggenes@result$pvalue < 0.05)
dim(clu9.go_Ssiggenes_lowP)
write.table(clu9.go_Ssiggenes_lowP,"~/LargeHmap/clu9.go_Ssiggenes_lowP_top100.txt")
clu9.go_Ssiggenes_lowP
lowPpathways <- rbind(clu0.go_Ssiggenes_lowP,clu1.go_Ssiggenes_lowP,
                      clu2.go_Ssiggenes_lowP,clu3.go_Ssiggenes_lowP,
                      clu4.go_Ssiggenes_lowP,clu5.go_Ssiggenes_lowP,
                      clu6.go_Ssiggenes_lowP,clu7.go_Ssiggenes_lowP,
                      clu8.go_Ssiggenes_lowP,clu9.go_Ssiggenes_lowP)
lowPpathways
dim(clu0.go_Ssiggenes_lowP)
dim(clu1.go_Ssiggenes_lowP)
dim(clu2.go_Ssiggenes_lowP)
dim(clu3.go_Ssiggenes_lowP)
dim(clu4.go_Ssiggenes_lowP)
dim(clu5.go_Ssiggenes_lowP)
dim(clu6.go_Ssiggenes_lowP)
dim(clu7.go_Ssiggenes_lowP)
dim(clu8.go_Ssiggenes_lowP)
dim(clu9.go_Ssiggenes_lowP)
write.table(lowPpathways,"~/LargeHmap/BigMPLTres07GO_top100.txt")
generatio <- cbind(clu0.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu1.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu2.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu3.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu4.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu5.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu6.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu7.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu8.go_Ssiggenes_lowP$GeneRatio[1:50],
                   clu9.go_Ssiggenes_lowP$GeneRatio[1:50])
dim(generatio)
str(generatio)
write.table(generatio,"~/LargeHmap/generatio.txt")
AllintegratedAnnorerunSR <- readRDS("3rd-PartA/stepone/intHarmony/Allint.filtered.clusteredAllres.anno.rds")
Idents(AllintegratedAnnorerunSR) <- AllintegratedAnnorerunSR$RNA_snn_res.5
markers.ori <- FindAllMarkers(AllintegratedAnnorerunSR, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers.ori,"3rd-PartA/steptwoAnno/ModifylablesbySR/AllintegratedAnnorerunSR_markersres5.txt",sep='\t')
DT::datatable(markers.ori)
library(dplyr)
all.top10 <- markers.ori %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
all.top10
write.table(all.top10,"3rd-PartA/steptwoAnno/ModifylablesbySR/AllintegratedAnnorerunSR_markersTop10byres5.txt",sep='\t')
library(dplyr)
all.top50 <- markers.ori %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) 
all.top50
write.table(all.top50,"3rd-PartA/steptwoAnno/ModifylablesbySR/AllintegratedAnnorerunSR_markersTop30byres5.txt",sep='\t')
DoHeatmap(AllintegratedAnnorerunSR, features = all.top10$gene) + NoLegend()
select_genes_all = split(all.top10$gene,all.top10$cluster)
DotPlot(object = AllintegratedAnnorerunSR,features=select_genes_all,assay="RNA",
        cols = c("lightgrey", "darkred"),group.by = "SRsortlabelsmain",
        dot.scale=10,col.min=1)+
  theme(axis.text.x.bottom =element_text(angle=45,vjust=0.5,hjust=0.5,size=12),
        axis.text.y=element_text(vjust=0.5,hjust=0.5,size=17),
  )
HPCA@colData$label.ont
HPCA <- HumanPrimaryCellAtlasData()
hpca
getwd()
cran.packages.need <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                   "gghalves", "Seurat", "SeuratObject", "methods",
                   "devtools", "BiocManager","data.table","doParallel",
                   "doRNG")
if(!requireNamespace(cran.packages.need,quietly = T)){
  install.packages(cran.packages.need,ask=F,update=F)
}
bioconductor.pac.need <- c("GSEABase", "AUCell", "SummarizedExperiment",
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree",
                           "Nebulosa")
if(!requireNamespace(bioconductor.pac.need)){
  BiocManager::install(bioconductor.pac.need,ask=F,update=F)
}
if(!requireNamespace("UCell",quietly = T)){
  devtools::install_github("carmonalab/UCell")
}
if(!requireNamespace("irGSEA",quietly = T)){
   devtools::install_github("chuiqin/irGSEA")
  }
BiocManager::install("Nebulosa")
library(irGSEA)
library(UCell)
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
levels(as.factor(Idents(BigMPLTintegrated)))
Idents(BigMPLTintegrated) <- BigMPLTintegrated$PLT.idents
BigMPLTintegrated.Enscore <- irGSEA.score(
  BigMPLTintegrated,assay="RNA",slot="data",seeds=1,ncores=1,min.cells = 3,
  min.feature = 0,custom = F, geneset = NULL,msigdb = T, species = "Homo sapiens",
  category = "H", subcategory = NULL, geneid = "symbol",
  method = c("AUCell","UCell","singscore","ssgsea"),kcdf = "Gaussian")
msigdbr::msigdbr_collections
Seurat::Assays(BigMPLTintegrated.Enscore)
str(BigMPLTintegrated.Enscore)
saveRDS(BigMPLTintegrated.Enscore,"~/BigMPLTintegrated.Enscore.rds")
result.dge <- irGSEA.integrate(BigMPLTintegrated.Enscore,
                               group.by = "PLT.idents",
                               method =c("AUCell","UCell","singscore","ssgsea"))
result.dge
class(result.dge)
saveRDS(result.dge,"~/BigMPLTintegrated-Enrich-result.dge.rds")
write.csv(result.dge$RRA,"~/BigMPLTintegrated-Enrich-RRA.csv")
write.csv(result.dge$UCell,"~/BigMPLTintegrated-Enrich-UCell.csv")
write.csv(result.dge$AUCell,"~/BigMPLTintegrated-Enrich-AUCell.csv")
write.csv(result.dge$singscore,"~/BigMPLTintegrated-Enrich-singscore.csv")
write.csv(result.dge$ssgsea,"~/BigMPLTintegrated-Enrich-ssgsea.csv")
result.dge.RRA.sorted <- read.csv("~/BigMPLTintegrated-Enrich-RRA.csv")
library(openxlsx)
result.dge.UCell.sorted <- read.csv("~/BigMPLTintegrated-Enrich-UCell.csv")
result.dge.AUCell.sorted <- read.csv("~/BigMPLTintegrated-Enrich-AUCell.csv")
result.dge.singscore.sorted <- read.csv("~/BigMPLTintegrated-Enrich-singscore.csv")
result.dge.ssgsea.sorted <- read.csv("~/BigMPLTintegrated-Enrich-ssgsea.csv")
result.dge.sorted <- result.dge
result.dge.sorted$RRA <- result.dge.RRA.sorted
result.dge.sorted$UCell <- result.dge.UCell.sorted
result.dge.sorted$AUCell <- result.dge.AUCell.sorted
result.dge.sorted$singscore <- result.dge.singscore.sorted
result.dge.sorted$ssgsea <- result.dge.ssgsea.sorted
# cols = c(pal_npg()(8),"
#          "
#          "
hmap1 <- irGSEA.heatmap(result.dge.sorted,method="RRA",
                       heatmap.width = 15,rowname.fointsize = 10,
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("darkred","darkblue"))
                       # # significance.color = c("white","
hmap1
result.dge.sorted$RRA$pvalue
hmap2 <- irGSEA.heatmap(result.dge.sorted,method="UCell",
                       heatmap.width = 15,rowname.fointsize = 10,
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("darkred","darkblue"),
                       # significance.color = c("white","
hmap2
hmap3 <- irGSEA.heatmap(result.dge.sorted,method="AUCell",
                       heatmap.width = 15,rowname.fointsize = 10,
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("darkred","darkblue"),
                       # significance.color = c("white","
hmap3
hmap4 <- irGSEA.heatmap(result.dge.sorted,method="singscore",
                       heatmap.width = 15,rowname.fointsize = 10,
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("darkred","darkblue"),
                       # significance.color = c("white","
hmap4
hmap5 <- irGSEA.heatmap(result.dge.sorted,method="ssgsea",
                       heatmap.width = 15,rowname.fointsize = 10,
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("darkred","darkblue"),
                       # significance.color = c("white","
hmap5
bubbleplot1 <- irGSEA.bubble(result.dge.sorted,,method="RRA",
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("blue","red"),
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot1
bubbleplot2 <- irGSEA.bubble(result.dge.sorted,,method="UCell",
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                        direction.color=c("blue","red"),
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot2
bubbleplot3 <- irGSEA.bubble(result.dge.sorted,,method="AUCell",
                       cluster.color =c("
                                        "
                                        "
                                        "
                                        "
                        direction.color=c("blue","red"),
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot3
bubbleplot4 <- irGSEA.bubble(result.dge.sorted,,method="singscore",
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                        direction.color=c("blue","red"),
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot4
bubbleplot5 <- irGSEA.bubble(result.dge.sorted,,method="ssgsea",
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("blue","red"),
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot5
scatterplot <- irGSEA.density.scatterplot(BigMPLTintegrated.Enscore,
                                          method="AUCell",
                                          show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE",
                                          reduction="tsne")
scatterplot
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
BigMPLTintegrateduse = BigMPLTintegrated[,!(BigMPLTintegrated$predicted.celltypePLTres0.7 %in% c("9"))]
table(BigMPLTintegrateduse$predicted.celltypePLTres0.7)
BigMPLTintegrateduse$predicted.celltypePLTres0.7 <- factor(BigMPLTintegrateduse$predicted.celltypePLTres0.7,levels = c("8","7","6","5","4","3","2","1","0"))
table(BigMPLTintegrateduse$predicted.celltypePLTres0.7)
BigMerge <- readRDS("~/BigMergeProjectedAll.rds")
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
BigMergereU <- BigMerge
MannualBasiccellNaming=data.frame(
  Allcell=BigMerge$predicted.celltypewithnew,
  Basiccell=BigMerge$predicted.celltype)
write.table(MannualBasiccellNaming,"~/MannualBasiccellNaming.txt")
library(openxlsx)
MannualBasiccellNaming <- read.xlsx("~/MannualBasiccellNaming.xlsx")
BigMergereU$predicted.celltypeBasicwithnewreU <- 
  MannualBasiccellNaming$Basiccell
length(BigMergereU@meta.data$predicted.celltypeBasicwithnewreU)
str(BigMPLTintegrated@meta.data$predicted.celltypePLTres0.7)
library(tidyverse)
PLT.idents <- BigMPLTintegrated$predicted.celltypePLTres0.7
PLT.idents <- str_replace_all(PLT.idents,"0","cluA")
PLT.idents <- str_replace_all(PLT.idents,"1","cluB")
PLT.idents <- str_replace_all(PLT.idents,"2","cluC")
PLT.idents <- str_replace_all(PLT.idents,"3","cluD")
PLT.idents <- str_replace_all(PLT.idents,"4","cluE")
PLT.idents <- str_replace_all(PLT.idents,"5","cluF")
PLT.idents <- str_replace_all(PLT.idents,"6","cluG")
PLT.idents <- str_replace_all(PLT.idents,"7","cluH")
PLT.idents <- str_replace_all(PLT.idents,"8","cluI")
PLT.idents <- str_replace_all(PLT.idents,"9","cluJ")
BigMPLTintegrated@meta.data$PLT.idents <- PLT.idents
levels(as.factor(BigMPLTintegrated@meta.data$PLT.idents))
table(BigMPLTintegrated@meta.data$PLT.idents)
table(BigMPLTintegrated@meta.data$predicted.celltypePLTres0.7)
BigMergereU@meta.data$CCCreU <- 
  BigMergereU@meta.data$predicted.celltypeBasicwithnewreU
n = 1
for (i in (1:695699)){
  if (BigMergereU@meta.data$CCCreU[i] == "Platelet"){
    BigMergereU@meta.data$CCCreU[i] <-
      BigMPLTintegrated@meta.data$PLT.idents[n]
    n = n+1
  }
  i = i+1
}
table(BigMergereU@meta.data$CCCreU)
table(BigMPLTintegrated@meta.data$predicted.celltypePLTres0.7)
table(BigMPLTintegrated@meta.data$PLT.idents)
levels(as.factor(BigMergereU@meta.data$CCCreU))
saveRDS(BigMergereU,"~/BigMergereU.rds")
saveRDS(BigMPLTintegrated,"~/BigMPLTintegratedProjected.rds")
BigMPLTintegrated$PLT.idents
library(tidyverse)
library(cowplot)
library(patchwork)
library(devtools)
library(igraph)
library(WGCNA)
   
library(scWGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())
set.seed(1)
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
PLTCovintegrated = 
  BigMPLTintegrated[,(BigMPLTintegrated$Covidsplit 
                      %in% c("use"))]
table(BigMPLTintegrated$Covidsplit)
table(PLTCovintegrated$predicted.celltypePLTres0.7)
PLTCovintegrateduse = 
  PLTCovintegrated[,!(PLTCovintegrated$PLT.idents 
                      %in% c("cluJ"))]
table(PLTCovintegrateduse$predicted.celltypePLTres0.7)
PLTCovintegrateduse$PLT.idents <- factor(PLTCovintegrateduse$PLT.idents,levels = c("cluA","cluB","cluC","cluD","cluE","cluF","cluG","cluH","cluI"))
table(PLTCovintegrateduse$PLT.idents)
PLTCovintegrateduse$predicted.celltypePLTres0.7 <- factor(PLTCovintegrateduse$predicted.celltypePLTres0.7,levels = c("0","1","2","3","4","5","6","7","8"))
table(PLTCovintegrateduse$predicted.celltypePLTres0.7)
PLTCovintegrateduse$PLT.idents
PLTCovintegrateduse$predicted.celltypePLTres0.7
dim(PLTCovintegrateduse)
dim(PLTCovintegrated)
DefaultAssay(PLTCovintegrateduse)
DefaultDimReduc(PLTCovintegrateduse)
Idents(PLTCovintegrateduse) <- PLTCovintegrateduse$PLT.idents
PLTCovintegrateduse.forWGCNA <- SetupForWGCNA(
  PLTCovintegrateduse,gene_select="fraction",
  fraction=0.05,
  wgcna_name="Clu0.7Frac0.05",
  group.by="PLT.idents")
length(PLTCovintegrateduse.forWGCNA@misc$Clu0.7Frac0.05$wgcna_genes)
PLTCovintegrateduse.forWGCNA@misc$active_wgcna
table(PLTCovintegrateduse.forWGCNA$predicted.celltypePLTres0.7)
table(PLTCovintegrateduse.forWGCNA$MSH.idents)
PLTCovintegrateduse.forWGCNA3 <- MetacellsByGroups(
  PLTCovintegrateduse.forWGCNA,
  group.by=c("PLT.idents","MSH.idents"),
  k = 20,
  reduction="ref.pca",
  max_shared=7,
  ident.group='PLT.idents'
)
PLTCovintegrateduse.forWGCNA3 <- 
  NormalizeMetacells(PLTCovintegrateduse.forWGCNA3)
table(PLTCovintegrateduse.forWGCNA3$metacell_grouping)
metacell_PLTCovintegrateduse3 <- GetMetacellObject(PLTCovintegrateduse.forWGCNA3)
metacell_PLTCovintegrateduse3@meta.data
metacell_PLTCovintegrateduse3
PLTCovintegrateduse.forWGCNA3
levels(as.factor(PLTCovintegrateduse.forWGCNA3$metacell_grouping))
table(PLTCovintegrateduse.forWGCNA3$PLT.idents,PLTCovintegrateduse.forWGCNA3$MSH.idents)
saveRDS(PLTCovintegrateduse.forWGCNA3,"~/PLTCovintegrateduse.forWGCNAv1.rds")
table(PLTCovintegrateduse.forWGCNA$predicted.celltypePLTres0.7)
DefaultAssay(PLTCovintegrateduse.forWGCNA)
DefaultDimReduc(PLTCovintegrateduse.forWGCNA)
Idents(PLTCovintegrateduse.forWGCNA) <- PLTCovintegrateduse.forWGCNA$PLT.idents
length(PLTCovintegrateduse.forWGCNA$PLT.idents)
length(PLTCovintegrateduse.forWGCNA$ID)
PLTCovintegrateduse.forWGCNAv2refpca <- MetacellsByGroups(
  PLTCovintegrateduse.forWGCNA,
  group.by=c("PLT.idents","ID"),
  k = 20,
  reduction="ref.pca",
  max_shared=7,
  ident.group='PLT.idents'
)
PLTCovintegrateduse.forWGCNAv2refpca <- 
  NormalizeMetacells(PLTCovintegrateduse.forWGCNAv2refpca)
table(PLTCovintegrateduse.forWGCNAv2$metacell_grouping)
metacell_PLTCovintegratedusev2refpca <- GetMetacellObject(PLTCovintegrateduse.forWGCNAv2refpca)
metacell_PLTCovintegratedusev2@meta.data
metacell_PLTCovintegratedusev2refpca
PLTCovintegrateduse.forWGCNAv2refpca
levels(as.factor(PLTCovintegrateduse.forWGCNAv2refpca$metacell_grouping))
saveRDS(PLTCovintegrateduse.forWGCNAv2refpca,"~/PLTCovintegrateduse.forWGCNAfrac0005v2refpca.rds")
table(PLTCovintegrateduse.forWGCNA$predicted.celltypePLTres0.7)
table(PLTCovintegrateduse.forWGCNA$MSH.idents)
PLTCovintegrateduse.forWGCNA$MSH.idents <-
  factor(PLTCovintegrateduse.forWGCNA$MSH.idents,levels = c("H","A","M","S","nan"))
PLTCovintegrateduse.forWGCNA3 <- MetacellsByGroups(
  PLTCovintegrateduse.forWGCNA,
  group.by=c("PLT.idents","Covidsplit"),
  k = 20,
  max_shared=7,
  ident.group='PLT.idents'
)
PLTCovintegrateduse.forWGCNA <- 
  NormalizeMetacells(PLTCovintegrateduse.forWGCNA)
metacell_PLTCovintegrateduse <- GetMetacellObject(PLTCovintegrateduse.forWGCNA)
metacell_PLTCovintegrateduse
PLTCovintegrateduse.forWGCNA
levels(as.factor(PLTCovintegrateduse.forWGCNA3$metacell_grouping))
k
table(PLTCovintegrateduse.forWGCNA$PLT.idents,PLTCovintegrateduse.forWGCNA$ID)
table(metacell_PLTCovintegrateduse$PLT.idents,metacell_PLTCovintegrateduse$ID)
saveRDS(PLTCovintegrateduse.forWGCNA,"~/PLTCovintegrateduse.forWGCNA.rds")
DefaultAssay(PLTCovintegrateduse.forWGCNAv2refpca)
DefaultDimReduc(PLTCovintegrateduse.forWGCNAv2refpca)
Idents(PLTCovintegrateduse.forWGCNAv2refpca) <- PLTCovintegrateduse.forWGCNAv2refpca$PLT.idents
PLTCovintegrateduse.forWGCNAv2refpca <- SetDatExpr(
  PLTCovintegrateduse.forWGCNAv2refpca, 
  group_name=c("cluA","cluB","cluC","cluD",
               "cluE","cluF"
               ),
  group.by = "PLT.idents",
  assay='RNA',use_metacells = T,
)
Idents(PLTCovintegrateduse.forWGCNAv2refpca)
allowWGCNAThreads(nThreads = 8)
PLTCovintegrateduse.forWGCNAv2refpca <- TestSoftPowers(
  PLTCovintegrateduse.forWGCNAv2refpca,
  powers=c(seq(1,10,by=1),seq(12,30,by=2)),
  networkType = "signed",
)
plot_list <- PlotSoftPowers(PLTCovintegrateduse.forWGCNAv2refpca,
                            point_size = 5,
                            text_size=3)
wrap_plots(plot_list,ncol=2)
power_table <- GetPowerTable(PLTCovintegrateduse.forWGCNAv2refpca)
head(power_table)
assignInNamespace("quantile.default",quantile.default,ns="stats",
                  pos="package:stats")
assign("quantile.default",quantile.default,
                  pos="package:stats")
lockBinding("quantile.default",as.environment("package:stats"))
stats:::quantile.default
package.version("stats")
PLTCovintegrateduse.forWGCNAv2refpca <- ConstructNetwork(
  PLTCovintegrateduse.forWGCNAv2refpca,
  soft_power=8,
  SetDatExpr=F,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  tom_outdir = "TOM",
  tom_name = "v2",
  overwrite_tom = T,
  replaceMissingAdjacencies=TRUE,na.rm=TRUE
)
PlotDendrogram(PLTCovintegrateduse.forWGCNAv2refpca,main="hdWGCNA Dendrogram")
saveRDS(PLTCovintegrateduse.forWGCNAv2refpca,"~/PLTCovintegrateduse.forWGCNAv2.rds")
PLTCovintegrateduse.forWGCNAv2refpca <- readRDS("~/PLTCovintegrateduse.forWGCNAv2.rds")
PLTCovintegrateduse.forWGCNA3 <- SetDatExpr(
  PLTCovintegrateduse.forWGCNA3, 
  group_name=c("cluA","cluB","cluC","cluD",
               "cluE","cluF","cluG","cluH","cluI"),
  group.by = "PLT.idents",
  assay='RNA',use_metacells = T,
  slot='data'
)
Idents(PLTCovintegrateduse.forWGCNA3)
PLTCovintegrateduse.forWGCNA3 <- TestSoftPowers(
  PLTCovintegrateduse.forWGCNA3,
  powers=c(seq(1,10,by=1),seq(12,30,by=2)),
  networkType = "signed",
  use_metacells = TRUE,
  setDatExpr =F,
  group_name=c("cluA","cluB","cluC","cluD",
               "cluE","cluF","cluG","cluH","cluI"),
  group.by = "PLT.idents",
)
plot_list <- PlotSoftPowers(PLTCovintegrateduse.forWGCNA3,
                            point_size = 5,
                            text_size=3)
wrap_plots(plot_list,ncol=2)
power_table <- GetPowerTable(PLTCovintegrateduse.forWGCNA3)
head(power_table)
PLTCovintegrateduse.forWGCNA3 <- ConstructNetwork(
  PLTCovintegrateduse.forWGCNA3,
  soft_power=8,
  SetDatExpr=F,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  tom_outdir = "TOM",
  tom_name = "v1",
  overwrite_tom = T,
  replaceMissingAdjacencies=T,na.rm=T
)
PlotDendrogram(PLTCovintegrateduse.forWGCNA3,main="hdWGCNA Dendrogram")
GetTOM
saveRDS(PLTCovintegrateduse.forWGCNA3,"PLTusingConventionalPC1-11/WGCNA/PLTCovMSH.cluster.forWGCNA.rds") 
PLTCovintegrateduse.forWGCNAv2refpca <-
  ScaleData(PLTCovintegrateduse.forWGCNAv2refpca,
            features=VariableFeatures(PLTCovintegrateduse.forWGCNAv2refpca))
GetModules(PLTCovintegrateduse.forWGCNAv2refpca)
DefaultAssay(PLTCovintegrateduse.forWGCNAv2refpca)
DefaultDimReduc(PLTCovintegrateduse.forWGCNAv2refpca)
Idents(PLTCovintegrateduse.forWGCNAv2refpca) <- PLTCovintegrateduse.forWGCNAv2refpca$PLT.idents
PLTCovintegrateduse.forWGCNAv2refpca <- ModuleEigengenesnew(
  PLTCovintegrateduse.forWGCNAv2refpca,
  scale.model.use = "linear",
  assay=NULL,
  pc_dim=1,
  group.by.vars = "PLT.idents"
)
hMEs <- GetMEs(PLTCovintegrateduse.forWGCNAv2refpca,wgcna_name = "Clu0.7Frac0.05")
MEs <- GetMEs(PLTCovintegrateduse.forWGCNAv2refpca,wgcna_name = "Clu0.7Frac0.05",harmonized = F)
PLTCovintegrateduse.forWGCNAv2refpca@misc
colnames(PLTCovintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05$MEs)
colnames(PLTCovintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05$hMEs)
PLTCovintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05$hMEs[1:3,]
ComputeModuleEigengenenew <- 
  function (seurat_obj, cur_mod, modules, group.by.vars = NULL, 
  verbose = TRUE, vars.to.regress = NULL, scale.model.use = "linear", 
  pc_dim = 1, assay = NULL, wgcna_name = NULL,...) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name %>% 
    as.character()
  if (CheckSeurat5()) {
    X <- SeuratObject::LayerData(seurat_obj, layer = "counts", 
      assay = assay)[cur_genes, ]
    X_dat <- SeuratObject::LayerData(seurat_obj, layer = "data", 
      assay = assay)[cur_genes, ]
  }
  else {
    X <- Seurat::GetAssayData(seurat_obj, slot = "counts", 
      assay = assay)[cur_genes, ]
    X_dat <- Seurat::GetAssayData(seurat_obj, slot = "data", 
      assay = assay)[cur_genes, ]
  }
  cur_seurat <- CreateSeuratObject(X, assay = assay, meta.data = seurat_obj@meta.data)
  if (CheckSeurat5()) {
    cur_seurat <- SetAssayData(cur_seurat, layer = "data", 
      new.data = X_dat, assay = assay)
  }
  else {
    cur_seurat <- SetAssayData(cur_seurat, slot = "data", 
      new.data = X_dat, assay = assay)
  }
  if (is.null(vars.to.regress)) {
    cur_seurat <- ScaleData(cur_seurat, features = rownames(cur_seurat), 
      model.use = scale.model.use)
  }
  else if (all(vars.to.regress %in% colnames(seurat_obj@meta.data))) {
    cur_seurat <- ScaleData(cur_seurat, features = rownames(cur_seurat), 
      model.use = scale.model.use, vars.to.regress = vars.to.regress)
  }
  else {
    stop(paste0("Some variables specified in vars.to.regress are not found in seurat_obj@meta.data"))
  }
  if (CheckSeurat5()) {
    cur_expr <- SeuratObject::GetAssayData(cur_seurat, layer = "data")
  }
  else {
    cur_expr <- Seurat::GetAssayData(cur_seurat, slot = "data")
  }
  expr <- Matrix::t(cur_expr)
  averExpr <- Matrix::rowSums(expr)/ncol(expr)
  cur_pca <- RunPCAnew(cur_seurat, features = cur_genes, 
    reduction.key = paste0("pca", cur_mod), verbose = verbose, 
    ...)@reductions$pca
  pc <- cur_pca@cell.embeddings[, pc_dim]
  pc_loadings <- cur_pca@feature.loadings[, pc_dim]
  pca_cor <- cor(averExpr, pc,use = 'pairwise.complete.obs')
  if (!is.null(group.by.vars)) {
    seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(embeddings = cur_pca@cell.embeddings, 
      assay = assay)
    cur_harmony <- harmony::RunHarmony(seurat_obj, group.by.vars = group.by.vars, 
      reduction.use = "ME", verbose = verbose, assay.use = assay, 
      ...)@reductions$harmony
    ha <- cur_harmony@cell.embeddings[, pc_dim]
    ha_loadings <- cur_pca@feature.loadings[, pc_dim]
    if (pca_cor < 0) {
      cur_harmony@cell.embeddings[, pc_dim] <- -ha
      ha_loadings <- -ha_loadings
    }
    seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(embeddings = cur_harmony@cell.embeddings, 
      assay = assay)
    seurat_obj <- SetMELoadings(seurat_obj, loadings = ha_loadings, 
      harmonized = TRUE, wgcna_name = wgcna_name)
  }
  if (pca_cor < 0) {
    cur_pca@cell.embeddings[, pc_dim] <- -pc
    pc_loadings <- -pc_loadings
  }
  seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(embeddings = cur_pca@cell.embeddings, 
    assay = assay)
  seurat_obj <- SetMELoadings(seurat_obj, loadings = pc_loadings, 
    harmonized = FALSE, wgcna_name = wgcna_name)
  seurat_obj
}
ModuleEigengenesnew <- 
function (seurat_obj, group.by.vars = NULL, modules = NULL, vars.to.regress = NULL, scale.model.use = "linear", verbose = TRUE, assay = NULL, pc_dim = 1, wgcna_name = NULL,...) {
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }
    CheckWGCNAName(seurat_obj, wgcna_name)
    harmonized <- !is.null(group.by.vars)
    if (harmonized & !any(grepl("ScaleData", seurat_obj@commands))) {
        stop("Need to run ScaleData before running ModuleEigengenes with group.by.vars option.")
    }
    if (is.null(assay)) {
        assay <- DefaultAssay(seurat_obj)
    }
    me_list <- list()
    harmonized_me_list <- list()
    seurat_obj <- SetMELoadings(seurat_obj, loadings = c(""), 
        harmonized = FALSE, wgcna_name = wgcna_name)
    if (harmonized) {
        seurat_obj <- SetMELoadings(seurat_obj, loadings = c(""), 
            harmonized = TRUE, wgcna_name = wgcna_name)
    }
    if (is.null(modules)) {
        modules <- GetModules(seurat_obj, wgcna_name)
        projected <- FALSE
    }
    else {
        projected <- TRUE
    }
    mods <- levels(modules$module)
    mods_loop <- mods
    for (cur_mod in mods_loop) {
        print(cur_mod)
        seurat_obj <- ComputeModuleEigengenenew(seurat_obj = seurat_obj, 
            cur_mod = cur_mod, modules = modules, group.by.vars = group.by.vars, 
            vars.to.regress = vars.to.regress, scale.model.use = scale.model.use, 
            verbose = verbose, pc_dim = pc_dim, assay = assay, 
            wgcna_name = wgcna_name)
        cur_me <- seurat_obj@reductions$ME@cell.embeddings[, 
            pc_dim]
        me_list[[cur_mod]] <- cur_me
        if (harmonized) {
            cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[, 
                pc_dim]
            harmonized_me_list[[cur_mod]] <- cur_harmonized_me
        }
    }
    me_df <- do.call(cbind, me_list)
    if (!projected) {
        me_df <- WGCNA::orderMEs(me_df)
    }
    seurat_obj <- SetMEs(seurat_obj, me_df, harmonized = FALSE, 
        wgcna_name)
    if (!is.null(group.by.vars)) {
        hme_df <- do.call(cbind, harmonized_me_list)
        if (!projected) {
            hme_df <- WGCNA::orderMEs(hme_df)
        }
        seurat_obj <- SetMEs(seurat_obj, hme_df, harmonized = TRUE, 
            wgcna_name)
    }
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
    modules$module <- factor(as.character(modules$module), levels = mods)
    seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)
    seurat_obj@reductions$ME <- NULL
    seurat_obj@reductions$ME_harmony <- NULL
    seurat_obj
}
library(WGCNA)
table(PLTCovintegrateduse.forWGCNAv2refpca$MSH.idents)
PLTCovintegrateduse.forWGCNAv2refpca$ConditionsH_M_S <- 
  as.factor(PLTCovintegrateduse.forWGCNAv2refpca$MSH.idents)
PLTCovintegrateduse.forWGCNAv2refpca$ConditionsH_M_S<- 
  factor(PLTCovintegrateduse.forWGCNAv2refpca$ConditionsH_M_S,levels=c("H","A","M","S"))
PLTCovintegrateduse.forWGCNAv2refpca$SEX <- 
  PLTCovintegrateduse.forWGCNAv2refpca$Sex
PLTCovintegrateduse.forWGCNAv2refpca$SEX[1:2842] <- 
  PLTCovintegrateduse.forWGCNAv2refpca$sex.idents[1:2842]
levels(as.factor(PLTCovintegrateduse.forWGCNAv2refpca$SEX))
library(tidyverse)
SEX.idents <- PLTCovintegrateduse.forWGCNAv2refpca$SEX 
SEX.idents <- str_replace_all(SEX.idents,"Female","F")
SEX.idents <- str_replace_all(SEX.idents,"Male","M")
PLTCovintegrateduse.forWGCNAv2refpca@meta.data$SEX <- SEX.idents
levels(as.factor(PLTCovintegrateduse.forWGCNAv2refpca@meta.data$SEX))
PLTCovintegrateduse.forWGCNAv2refpca$SexF_M <- 
  as.factor(PLTCovintegrateduse.forWGCNAv2refpca$SEX)
PLTCovintegrateduse.forWGCNAv2refpca$SexF_M <- 
  factor(PLTCovintegrateduse.forWGCNAv2refpca$SexF_M,levels = c("F","M"))
PLTCovintegrateduse.forWGCNAv2refpca$Age_increase <- 
  as.numeric(PLTCovintegrateduse.forWGCNAv2refpca$age.idents)
phenolist <- c("ConditionsH_M_S","SexF_M")
str(PLTCovintegrateduse.forWGCNAv2refpca@meta.data[,phenolist])
PLTCovintegrateduse.forWGCNAv2refpca <- ModuleTraitCorrelation(
  PLTCovintegrateduse.forWGCNAv2refpca,
  traits = phenolist,
  features ="MEs",
  cor_method = "pearson",
  group.by = "PLT.idents"
)
PLTCovMSH.cluster.forWGCNA2 <- ModuleTraitCorrelation(
  PLTCovMSH.cluster.forWGCNA,
  traits = phenolist,
  features ="hMEs",
  cor_method = "pearson",
  group.by = "RNA_snn_res.0.05"
)
moduleCor <- GetModuleTraitCorrelation(PLTCovintegrateduse.forWGCNAv2refpca)
names(moduleCor)
moduleCor2 <- GetModuleTraitCorrelation(PLTCovMSH.cluster.forWGCNA2)
names(moduleCor2)
PlotModuleTraitCorrelation(
  PLTCovintegrateduse.forWGCNAv2refpca,
  label="fdr",
  label_symbol = "stars",
  text_size = 2,
  text_digits = 2,
  text_color = "white",
  # #high_color = "
  # #mid_color = "
  # #low_color = "
  plot_max = 0.2,
  combine=T 
)
PlotModuleTraitCorrelation(
  PLTCovMSH.cluster.forWGCNA2,
  label="fdr",
  label_symbol = "stars",
  text_size = 2,
  text_digits = 2,
  text_color = "white",
  # #high_color = "
  # #mid_color = "
  # #low_color = "
  plot_max = 0.2,
  combine=T 
)
phenolist2 <- c("Conditions","SexF_M","Age_increase")
PLTCovMSH.cluster.forWGCNA.3 <- ModuleTraitCorrelation(
  PLTCovMSH.cluster.forWGCNA,
  traits = phenolist2,
  features ="hMEs",
  cor_method = "pearson",
  group.by = "RNA_snn_res.0.5"
)
moduleCor <- GetModuleTraitCorrelation(PLTCovMSH.cluster.forWGCNA.3)
names(moduleCor)
PlotModuleTraitCorrelation(
  PLTCovMSH.cluster.forWGCNA.3,
  label="fdr",
  label_symbol = "stars",
  text_size = 3,
  text_digits = 2,
  text_color = "black",
  # #low_color = "
  # #mid_color = "
  # #high_color = "
  plot_max = 0.2,
  combine=T 
)
phenolist2 <- c("Conditions","SexF_M","Age_increase")
PLTCovMSH.cluster.forWGCNA.4 <- ModuleTraitCorrelation(
  PLTCovMSH.cluster.forWGCNA,
  traits = phenolist2,
  features ="hMEs",
  cor_method = "pearson"
)
moduleCor <- GetModuleTraitCorrelation(PLTCovMSH.cluster.forWGCNA.4)
names(moduleCor)
PlotModuleTraitCorrelation(
  PLTCovMSH.cluster.noA.forWGCNA.4,
  label="fdr",
  label_symbol = "stars",
  text_size = 3,
  text_digits = 2,
  text_color = "black",
  # #low_color = "
  # #mid_color = "
  # #high_color = "
  plot_max = 0.2,
  combine=T 
)+
  theme(axis.text.x = element_text(size=10),
        axis.text.y.right = element_text(size=10),
        axis.text.y.left = element_text(size=10))
phenolist2 <- c("Conditions","SexF_M","Age_increase")
PLTCovMSH.cluster.forWGCNA <- ModuleTraitCorrelation(
  PLTCovMSH.cluster.forWGCNA,
  traits = phenolist2,
  features ="MEs",
  cor_method = "pearson",
  group.by = "RNA_snn_res.0.05"
)
PlotModuleTraitCorrelation(
  PLTCovMSH.cluster.forWGCNA,
  label="fdr",
  label_symbol = "stars",
  text_size = 3,
  text_digits = 2,
  text_color = "black",
  #low_color = "
  #mid_color = "
  #high_color = "
  plot_max = 0.2,
  combine=T 
)+
  theme(axis.text.x = element_text(size=10),
        axis.text.y.right = element_text(size=10),
        axis.text.y.left = element_text(size=10))
PLTCovintegrateduse.forWGCNAv2refpca <- ModuleConnectivity(
  PLTCovintegrateduse.forWGCNAv2refpca,
  group_name = "PLT.idents",
  corFnc = "bicor",
  corOptions = "use='p'",
  harmonized = TRUE,
  assay = "RNA",
  slot = "data",
)
PLTCovintegrateduse.forWGCNAv2refpca <- ResetModuleNames(
  PLTCovintegrateduse.forWGCNAv2refpca,
  new_name="PLT-M"
)
modules <- GetModules(PLTCovintegrateduse.forWGCNAv2refpca)
print(levels(modules$module))
head(modules[,1:6])
dim(modules)
library(hdWGCNA)
library(conflicted)
conflict_prefer("select","dplyr")
p <- PlotKMEs(PLTCovintegrateduse.forWGCNAv2refpca,
              ncol=3,
              n_hubs=20,
              text_size=2,
              plot_widths = c(3,2)
              )
p
hubgenetop100 <- GetHubGenes(PLTCovintegrateduse.forWGCNAv2refpca,n_hubs = 100)
head(hubgenetop100)
hubgenetop100
write.table(hubgenetop100,"~/hubgenetop100.txt")
saveRDS(PLTCovintegrateduse.forWGCNAv2refpca,"~/PLTCovintegrateduse.forWGCNAv2TOMrefpcaMESpca.rds")
PLTCovintegrateduse.forWGCNAv2refpca <-
  readRDS("~/PLTCovintegrateduse.forWGCNAv2TOMrefpcaMESpca.rds")
PLTCovintegrateduse.forWGCNAv2refpca <- ModuleExprScore(
  PLTCovintegrateduse.forWGCNAv2refpca,
  n_genes=100,
  method = "Seurat"
)
PLTCovintegrateduse.forWGCNAv2refpca <- ModuleExprScore(
  PLTCovintegrateduse.forWGCNAv2refpca,
  n_genes = 50,
  method = "UCell"
)
PLTCovMSH.cluster.forWGCNA
plot_list1 <- ModuleFeaturePlot(
  PLTCovintegrateduse.forWGCNAv2refpca,
  reduction="ref.umap",
  features = "hMEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 1,
  label_legend = F,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list1,ncol=3)
plot_list2 <- ModuleFeaturePlot(
  PLTCovMSH.cluster.forWGCNA,
  features = "scores",
  ucell=TRUE,
  point_size = 0.3,
  reduction="tsne"
)
wrap_plots(plot_list2,ncol=3)
plot_list3 <- ModuleFeaturePlot(
  PLTCovMSH.cluster.forWGCNA,
  reduction="tsne",
  features = "MEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 0.5,
  label_legend = F,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list3,ncol=3)
plot_list1 <- ModuleFeaturePlot(
  PLTCovintegrateduse.forWGCNAv2refpca,
  reduction="ref.umap",
  features = "hMEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 1,
  label_legend = F,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list1,ncol=3)
plot_list2 <- ModuleFeaturePlot(
  PLTCovintegrateduse.forWGCNAv2refpca,
  features = "scores",
  ucell=TRUE,
  point_size = 0.3,
  reduction="ref.umap"
)
wrap_plots(plot_list2,ncol=3)
plot_list3 <- ModuleFeaturePlot(
  PLTCovintegrateduse.forWGCNAv2refpca,
  reduction="ref.umap",
  features = "MEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 0.5,
  label_legend = F,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list3,ncol=3)
ModuleCorrelogram(PLTCovintegrateduse.forWGCNAv2refpca,
                  exclude_grey = T,
                  features ="MEs")
PLTCovintegrateduse.forWGCNAv2refpca.b <- PLTCovintegrateduse.forWGCNAv2refpca
Idents(PLTCovintegrateduse.forWGCNAv2refpca.b) <- PLTCovintegrateduse.forWGCNAv2refpca.b$PLT.idents
levels(as.factor(Idents(PLTCovintegrateduse.forWGCNAv2refpca.b)))
markers_forhere <- Seurat::FindAllMarkers(
  PLTCovintegrateduse.forWGCNAv2refpca.b,
  only.pos=T,
  logfc.threshold = 1
) 
library(GeneOverlap)
overlap_MwithMarker <- OverlapModulesDEGs(
  PLTCovintegrateduse.forWGCNAv2refpca.b,
  deg_df = markers_forhere,
  fc_cutoff=1
)
overlap_MwithMarker
OverlapDotPlot(
  overlap_MwithMarker,
  plot_var="odds_ratio")+
  ggtitle("Modules mapping to markers of PLT-subclusters")
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
library(dplyr)
theme_set(theme_cowplot())
library(enrichR)
library(GeneOverlap)
set.seed(1)
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021',"GO_Molecular_Function_2021")
PLTCovMSH.cluster.forWGCNA <- RunEnrichr(
  PLTCovMSH.cluster.forWGCNA,
  dbs=dbs,
  max_genes=100
)
str(PLTCovMSH.cluster.forWGCNA@misc$fourthClu1.5Frac0.05withAres05$enrichr_table)
enrich_df <- PLTCovMSH.cluster.forWGCNA@misc$fourthClu1.5Frac0.05withAres05$enrichr_table
enrich_df
write.table(enrich_df,"PLTusingConventionalPC1-11/WGCNA/fourth/enrich_df.txt")
EnrichrBarPlot(
  PLTCovMSH.cluster.forWGCNA,
  outdir="PLTusingConventionalPC1-11/WGCNA/fourth/enrich_plots",
  n_terms=10,
  plot_size=c(5,7),
  logscale=T
)
EnrichrDotPlot(
  PLTCovMSH.cluster.forWGCNA,
  mods="all",
  database = "GO_Biological_Process_2021",
  n_terms=5,
  break_ties = F
)
EnrichrDotPlot(
  PLTCovMSH.cluster.forWGCNA,
  mods="PLT-M3",
  database = "GO_Biological_Process_2021",
  n_terms=10
)
EnrichrDotPlot(
  PLTCovMSH.cluster.forWGCNA,
  mods="PLT-M6",
  database = "GO_Biological_Process_2021",
  n_terms=10
)
EnrichrDotPlot <- function (seurat_obj, database, mods = "all", n_terms = 3, break_ties = TRUE, logscale = TRUE, wgcna_name = NULL, ...) 
{
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }
    modules <- GetModules(seurat_obj, wgcna_name)
    if (mods == "all") {
        mods <- levels(modules$module)
        mods <- mods[mods != "grey"]
    }
    enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)
    mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct
    enrichr_df$color <- mod_colors[match(enrichr_df$module, mod_colors$module), 
        "color"]
    wrapText <- function(x, len) {
        sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), 
            USE.NAMES = FALSE)
    }
    plot_df <- enrichr_df %>% subset(db == database & module %in% 
        mods) %>% group_by(module) %>% top_n(n_terms, wt = Combined.Score)
    if (break_ties) {
        plot_df <- do.call(rbind, lapply(plot_df %>% group_by(module) %>% 
            group_split, function(x) {
            x[sample(n_terms), ]
        }))
    }
    plot_df$Term <- wrapText(plot_df$Term, 45)
    plot_df$module <- factor(as.character(plot_df$module), levels = levels(modules$module))
    plot_df <- arrange(plot_df, module)
    plot_df$Term <- factor(as.character(plot_df$Term), levels = unique(as.character(plot_df$Term)))
    if (logscale) {
        plot_df$Combined.Score <- log(plot_df$Combined.Score)
        lab <- "Enrichment\nlog(combined score)"
        x <- 0.2
    }
    else {
        lab <- "Enrichment\n(combined score)"
        x <- 5
    }
    p <- plot_df %>% ggplot(aes(x = module, y = Term)) + geom_point(aes(size = Combined.Score), 
        color = plot_df$color) + RotatedAxis() + ylab("") + xlab("") + 
        labs(size = lab) + scale_y_discrete(limits = rev) + ggtitle(database) + 
        theme(plot.title = element_text(hjust = 0.5), axis.line.x = element_blank(), 
            axis.line.y = element_blank(), panel.border = element_rect(colour = "black", 
                fill = NA, size = 1))
    p
}
wgcna_name <- PLTCovMSH.cluster.forWGCNA@misc$active_wgcna
mods <- levels(modules$module)
mods <- mods[mods != "grey"]
enrichr_df <- GetEnrichrTable(PLTCovMSH.cluster.forWGCNA, wgcna_name)
mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct
enrichr_df$color <- mod_colors[match(enrichr_df$module, mod_colors$module),"color"]
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
plot_df <- enrichr_df %>% subset(db == database & module %in% mods) %>% group_by(module) %>% top_n(n_terms, wt = Combined.Score)
if (break_ties) {
    plot_df <- do.call(rbind, lapply(plot_df %>% group_by(module) %>% 
        group_split, function(x) {
        x[sample(n_terms), ]
    }))
}
plot_df$Term <- wrapText(plot_df$Term, 45)
plot_df$module <- factor(as.character(plot_df$module), levels = levels(modules$module))
plot_df <- arrange(plot_df, module)
plot_df$Term <- factor(as.character(plot_df$Term), levels = unique(as.character(plot_df$Term)))
if (logscale) {
    plot_df$Combined.Score <- log(plot_df$Combined.Score)
    lab <- "Enrichment\nlog(combined score)"
    x <- 0.2
}
else {
    lab <- "Enrichment\n(combined score)"
    x <- 5
}
p <- plot_df %>% ggplot(aes(x = module, y = Term)) + 
  geom_point(aes(size = Combined.Score), 
        color = plot_df$color) + RotatedAxis() + ylab("") + xlab("") + 
        labs(size = lab) + scale_y_discrete(limits = rev) + ggtitle(database) + 
        theme(plot.title = element_text(hjust = 0.5), axis.line.x = element_blank(), 
            axis.line.y = element_blank(), panel.border = element_rect(colour = "black", 
                fill = NA, size = 1))
p
ModuleNetworkPlot(PLTCovMSH.cluster.forWGCNA)
HubGeneNetworkPlot(
  PLTCovMSH.cluster.forWGCNA,
  mods="all",
  n_hubs=10,
  n_other=50,
  edge_prop=0.75
)
set.seed(1)
HubGeneNetworkPlot(
  PLTCovMSH.cluster.forWGCNA,
  mods="all",
  n_hubs=10,
  n_other=20,
  edge_prop=0.75,
  vertex.label.cex = 1
)
g <- HubGeneNetworkPlot(
  PLTCovMSH.cluster.forWGCNA,
  return_graph = T
)
g
PLTCovMSH.cluster.forWGCNA <- RunModuleUMAP(
  PLTCovMSH.cluster.forWGCNA,
  n_hubs=20,
  n_neighbors = 15,
  min_dist = 0.1
)
umap_table <- GetModuleUMAP(PLTCovMSH.cluster.forWGCNA)
ggplot(umap_table,aes(x=UMAP1,y=UMAP2))+
  geom_point(
    color=umap_table$color,
    size=umap_table$kME*2
  )+
  umap_theme()
ModuleUMAPPlot(
  PLTCovMSH.cluster.forWGCNA,
  edge.alpha = 0.25,
  sample_edges = T,
  edge_prop = 0.1,
  label_hubs = 2,
  keep_grey_edges = F
)
PLTCovMSH.cluster.forWGCNA <- RunModuleUMAP(
  PLTCovMSH.cluster.forWGCNA,
  n_hubs = 10,
  n_neighbors=15,
  min_dist=0.1,
  supervised=TRUE,
  target_weight=0.5
)
umap_df <- GetModuleUMAP(PLTCovMSH.cluster.forWGCNA)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,
   size=umap_df$kME*2
  ) +
  umap_theme()
ModuleUMAPPlot(
  PLTCovMSH.cluster.forWGCNA,
  edge.alpha = 0.25,
  sample_edges = T,
  edge_prop = 0.1,
  label_hubs = 2,
  keep_grey_edges = F
)
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
BigMPLTintegrateduse = 
  BigMPLTintegrated[,!(BigMPLTintegrated$PLT.idents 
                      %in% c("cluJ"))]
table(BigMPLTintegrateduse$PLT.idents)
BigMPLTintegrateduse$PLT.idents
dim(BigMPLTintegrateduse)
dim(BigMPLTintegrated)
DefaultAssay(BigMPLTintegrateduse)
DefaultDimReduc(BigMPLTintegrateduse)
Idents(BigMPLTintegrateduse) <- BigMPLTintegrateduse$PLT.idents
BigMPLTintegrateduse.forWGCNA <- SetupForWGCNA(
  BigMPLTintegrateduse,gene_select="fraction",
  fraction=0.05,
  wgcna_name="Clu0.7Frac0.05v2",
  group.by="PLT.idents")
length(BigMPLTintegrateduse.forWGCNA@misc$Clu0.7Frac0.05v2$wgcna_genes)
BigMPLTintegrateduse.forWGCNA@misc$active_wgcna
DefaultAssay(BigMPLTintegrateduse.forWGCNA)
DefaultDimReduc(BigMPLTintegrateduse.forWGCNA)
Idents(BigMPLTintegrateduse.forWGCNA) <- BigMPLTintegrateduse.forWGCNA$PLT.idents
length(BigMPLTintegrateduse.forWGCNA$PLT.idents)
length(BigMPLTintegrateduse.forWGCNA$ID)
BigMPLTintegrateduse.forWGCNAv2refpca <- MetacellsByGroups(
  BigMPLTintegrateduse.forWGCNA,
  group.by=c("PLT.idents","ID"),
  k = 20,
  reduction="ref.pca",
  max_shared=7,
  ident.group='PLT.idents'
)
BigMPLTintegrateduse.forWGCNAv2refpca <- 
  NormalizeMetacells(BigMPLTintegrateduse.forWGCNAv2refpca)
table(BigMPLTintegrateduse.forWGCNAv2refpca$metacell_grouping)
metacell_BigMPLTintegratedusev2refpca <- GetMetacellObject(BigMPLTintegrateduse.forWGCNAv2refpca)
metacell_BigMPLTintegratedusev2refpca@meta.data
metacell_BigMPLTintegratedusev2refpca
BigMPLTintegrateduse.forWGCNAv2refpca
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca$metacell_grouping))
DefaultAssay(BigMPLTintegrateduse.forWGCNAv2refpca)
DefaultDimReduc(BigMPLTintegrateduse.forWGCNAv2refpca)
Idents(BigMPLTintegrateduse.forWGCNAv2refpca) <-
  BigMPLTintegrateduse.forWGCNAv2refpca$PLT.idents
BigMPLTintegrateduse.forWGCNAv2refpca <- SetDatExpr(
  BigMPLTintegrateduse.forWGCNAv2refpca, 
  group_name=c("cluA","cluB","cluC","cluD",
               "cluE","cluF"
               ),
  group.by = "PLT.idents",
  assay='RNA',use_metacells = T
)
dim(BigMPLTintegrateduse.forWGCNAv2refpca)
Idents(BigMPLTintegrateduse.forWGCNAv2refpca)
allowWGCNAThreads(nThreads = 8)
BigMPLTintegrateduse.forWGCNAv2refpca <- TestSoftPowers(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  powers=c(seq(1,10,by=1),seq(12,30,by=2)),
  networkType = "signed",
)
plot_list <- PlotSoftPowers(BigMPLTintegrateduse.forWGCNAv2refpca,
                            point_size = 5,
                            text_size=3)
wrap_plots(plot_list,ncol=2)
power_table <- GetPowerTable(BigMPLTintegrateduse.forWGCNAv2refpca)
head(power_table)
BigMPLTintegrateduse.forWGCNAv2refpca <- ConstructNetwork(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  soft_power=7,
  SetDatExpr=F,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  tom_outdir = "TOM",
  tom_name = "v2",
  overwrite_tom = T,
  replaceMissingAdjacencies=TRUE,na.rm=TRUE
)
PlotDendrogram(BigMPLTintegrateduse.forWGCNAv2refpca,main="hdWGCNA Dendrogram")
saveRDS(BigMPLTintegrateduse.forWGCNAv2refpca,"~/BigMPLTintegrateduse.forWGCNAv2refpcauntilTOM.rds")
BigMPLTintegrateduse.forWGCNAv2refpca <-
  ScaleData(BigMPLTintegrateduse.forWGCNAv2refpca,
            features=VariableFeatures(BigMPLTintegrateduse.forWGCNAv2refpca))
table(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$module)
conflicts_prefer(WGCNA::cor)
BigMPLTintegrateduse.forWGCNAv2refpca <- ModuleEigengenesnew(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  scale.model.use = "linear",
  assay=NULL,
  pc_dim=1,
  group.by.vars = "PLT.idents"
)
sum(BigMPLTintegrateduse.forWGCNAv2refpca@assays$RNA@meta.features$vst.variance==0)
table(is.na(BigMPLTintegrateduse.forWGCNAv2refpca@assays$RNA@counts@x))
table(is.na(BigMPLTintegrateduse.forWGCNAv2refpca@assays$RNA@data@x))
BigMPLTintegrateduse.forWGCNAv2refpca@assays$RNA@data@x[is.na(BigMPLTintegrateduse.forWGCNAv2refpca@assays$RNA@data@x)] <- 0
BigMPLTintegrateduse.forWGCNAv2refpca <- ModuleEigengenes(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  scale.model.use = "linear",
  assay=NULL,
  pc_dim=1,
  group.by.vars = "PLT.idents"
)
hMEs <- GetMEs(BigMPLTintegrateduse.forWGCNAv2refpca,wgcna_name = "Clu0.7Frac0.05v2")
MEs <- GetMEs(BigMPLTintegrateduse.forWGCNAv2refpca,
              wgcna_name = "Clu0.7Frac0.05v2",harmonized = F)
BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules
colnames(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$MEs)
colnames(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$hMEs)
BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05$hMEs[1:3,]
ComputeModuleEigengenenew <- 
  function (seurat_obj, cur_mod, modules, group.by.vars = NULL, 
  verbose = TRUE, vars.to.regress = NULL, scale.model.use = "linear", 
  pc_dim = 1, assay = NULL, wgcna_name = NULL,...) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name %>% 
    as.character()
  if (CheckSeurat5()) {
    X <- SeuratObject::LayerData(seurat_obj, layer = "counts", 
      assay = assay)[cur_genes, ]
    X_dat <- SeuratObject::LayerData(seurat_obj, layer = "data", 
      assay = assay)[cur_genes, ]
  }
  else {
    X <- Seurat::GetAssayData(seurat_obj, slot = "counts", 
      assay = assay)[cur_genes, ]
    X_dat <- Seurat::GetAssayData(seurat_obj, slot = "data", 
      assay = assay)[cur_genes, ]
  }
  cur_seurat <- CreateSeuratObject(X, assay = assay, meta.data = seurat_obj@meta.data)
  if (CheckSeurat5()) {
    cur_seurat <- SetAssayData(cur_seurat, layer = "data", 
      new.data = X_dat, assay = assay)
  }
  else {
    cur_seurat <- SetAssayData(cur_seurat, slot = "data", 
      new.data = X_dat, assay = assay)
  }
  if (is.null(vars.to.regress)) {
    cur_seurat <- ScaleData(cur_seurat, features = rownames(cur_seurat), 
      model.use = scale.model.use)
  }
  else if (all(vars.to.regress %in% colnames(seurat_obj@meta.data))) {
    cur_seurat <- ScaleData(cur_seurat, features = rownames(cur_seurat), 
      model.use = scale.model.use, vars.to.regress = vars.to.regress)
  }
  else {
    stop(paste0("Some variables specified in vars.to.regress are not found in seurat_obj@meta.data"))
  }
  if (CheckSeurat5()) {
    cur_expr <- SeuratObject::GetAssayData(cur_seurat, layer = "data")
  }
  else {
    cur_expr <- Seurat::GetAssayData(cur_seurat, slot = "data")
  }
  expr <- Matrix::t(cur_expr)
  averExpr <- Matrix::rowSums(expr)/ncol(expr)
  cur_pca <- Seurat::RunPCA(cur_seurat, features = cur_genes, 
    reduction.key = paste0("pca", cur_mod), verbose = verbose, 
    ...)@reductions$pca
  pc <- cur_pca@cell.embeddings[, pc_dim]
  pc_loadings <- cur_pca@feature.loadings[, pc_dim]
  pca_cor <- cor(averExpr, pc,use = 'pairwise.complete.obs')
  if (!is.null(group.by.vars)) {
    seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(embeddings = cur_pca@cell.embeddings, 
      assay = assay)
    cur_harmony <- harmony::RunHarmony(seurat_obj, group.by.vars = group.by.vars, 
      reduction.use = "ME", verbose = verbose, assay.use = assay, 
      ...)@reductions$harmony
    ha <- cur_harmony@cell.embeddings[, pc_dim]
    ha_loadings <- cur_pca@feature.loadings[, pc_dim]
    if (pca_cor < 0) {
      cur_harmony@cell.embeddings[, pc_dim] <- -ha
      ha_loadings <- -ha_loadings
    }
    seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(embeddings = cur_harmony@cell.embeddings, 
      assay = assay)
    seurat_obj <- SetMELoadings(seurat_obj, loadings = ha_loadings, 
      harmonized = TRUE, wgcna_name = wgcna_name)
  }
  if (pca_cor < 0) {
    cur_pca@cell.embeddings[, pc_dim] <- -pc
    pc_loadings <- -pc_loadings
  }
  seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(embeddings = cur_pca@cell.embeddings, 
    assay = assay)
  seurat_obj <- SetMELoadings(seurat_obj, loadings = pc_loadings, 
    harmonized = FALSE, wgcna_name = wgcna_name)
  seurat_obj
}
ModuleEigengenesnew <- 
function (seurat_obj, group.by.vars = NULL, modules = NULL, vars.to.regress = NULL, scale.model.use = "linear", verbose = TRUE, assay = NULL, pc_dim = 1, wgcna_name = NULL,...) {
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }
    CheckWGCNAName(seurat_obj, wgcna_name)
    harmonized <- !is.null(group.by.vars)
    if (harmonized & !any(grepl("ScaleData", seurat_obj@commands))) {
        stop("Need to run ScaleData before running ModuleEigengenes with group.by.vars option.")
    }
    if (is.null(assay)) {
        assay <- DefaultAssay(seurat_obj)
    }
    me_list <- list()
    harmonized_me_list <- list()
    seurat_obj <- SetMELoadings(seurat_obj, loadings = c(""), 
        harmonized = FALSE, wgcna_name = wgcna_name)
    if (harmonized) {
        seurat_obj <- SetMELoadings(seurat_obj, loadings = c(""), 
            harmonized = TRUE, wgcna_name = wgcna_name)
    }
    if (is.null(modules)) {
        modules <- GetModules(seurat_obj, wgcna_name)
        projected <- FALSE
    }
    else {
        projected <- TRUE
    }
    mods <- levels(modules$module)
    mods_loop <- mods
    for (cur_mod in mods_loop) {
        print(cur_mod)
        seurat_obj <- ComputeModuleEigengenenew(seurat_obj = seurat_obj, 
            cur_mod = cur_mod, modules = modules, group.by.vars = group.by.vars, 
            vars.to.regress = vars.to.regress, scale.model.use = scale.model.use, 
            verbose = verbose, pc_dim = pc_dim, assay = assay, 
            wgcna_name = wgcna_name)
        cur_me <- seurat_obj@reductions$ME@cell.embeddings[, 
            pc_dim]
        me_list[[cur_mod]] <- cur_me
        if (harmonized) {
            cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[, 
                pc_dim]
            harmonized_me_list[[cur_mod]] <- cur_harmonized_me
        }
    }
    me_df <- do.call(cbind, me_list)
    if (!projected) {
        me_df <- WGCNA::orderMEs(me_df)
    }
    seurat_obj <- SetMEs(seurat_obj, me_df, harmonized = FALSE, 
        wgcna_name)
    if (!is.null(group.by.vars)) {
        hme_df <- do.call(cbind, harmonized_me_list)
        if (!projected) {
            hme_df <- WGCNA::orderMEs(hme_df)
        }
        seurat_obj <- SetMEs(seurat_obj, hme_df, harmonized = TRUE, 
            wgcna_name)
    }
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
    modules$module <- factor(as.character(modules$module), levels = mods)
    seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)
    seurat_obj@reductions$ME <- NULL
    seurat_obj@reductions$ME_harmony <- NULL
    seurat_obj
}
BigMPLTintegrateduse.forWGCNAv2refpca <- 
  readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.rds")
BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2
library(WGCNA)
table(BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents)
BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_M_S <- 
  as.factor(BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents)
BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_M_S<- 
  factor(BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_M_S,levels=c("H","A","M","S"))
BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_m_S <- 
  as.factor(BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents)
BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_m_S<- 
  factor(BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_m_S,levels=c("H","M","S"))
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$ConditionsH_m_S))
BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_AM_S <- 
  as.factor(BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents)
library(tidyverse)
Cond.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_AM_S
Cond.idents <- str_replace_all(Cond.idents,"A","am")
Cond.idents <- str_replace_all(Cond.idents,"M","am")
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$ConditionsH_AM_S <- Cond.idents
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$ConditionsH_AM_S))
BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_AM_S <- 
  factor(BigMPLTintegrateduse.forWGCNAv2refpca$ConditionsH_AM_S,levels = c("H","am","S"))
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$ConditionsH_AM_S))
BigMPLTintegrateduse.forWGCNAv2refpca$SEX <- 
  BigMPLTintegrateduse.forWGCNAv2refpca$Sex
BigMPLTintegrateduse.forWGCNAv2refpca$SEX[1:2842] <- 
  BigMPLTintegrateduse.forWGCNAv2refpca$sex.idents[1:2842]
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca$SEX))
library(tidyverse)
SEX.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$SEX 
SEX.idents <- str_replace_all(SEX.idents,"Female","F")
SEX.idents <- str_replace_all(SEX.idents,"Male","M")
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$SEX <- SEX.idents
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$SEX))
BigMPLTintegrateduse.forWGCNAv2refpca$SexF_M <- 
  as.factor(BigMPLTintegrateduse.forWGCNAv2refpca$SEX)
BigMPLTintegrateduse.forWGCNAv2refpca$SexF_M <- 
  factor(BigMPLTintegrateduse.forWGCNAv2refpca$SexF_M,levels = c("F","M"))
table(BigMPLTintegrateduse.forWGCNAv2refpca$age.idents)
table(BigMPLTintegrateduse.forWGCNAv2refpca$Age_interval)
BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase <- 
  as.character(BigMPLTintegrateduse.forWGCNAv2refpca$Age_interval)
BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase[1:2842] <- 
  BigMPLTintegrateduse.forWGCNAv2refpca$age.idents[1:2842]
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase))
library(tidyverse)
AGE.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase 
AGE.idents <- str_replace_all(AGE.idents,"38","(30, 39]")
AGE.idents <- str_replace_all(AGE.idents,"46","(40, 49]")
AGE.idents <- str_replace_all(AGE.idents,"54","(50, 59]")
AGE.idents <- str_replace_all(AGE.idents,"61","(60, 69]")
AGE.idents <- str_replace_all(AGE.idents,"62","(60, 69]")
AGE.idents <- str_replace_all(AGE.idents,"63","(60, 69]")
AGE.idents <- str_replace_all(AGE.idents,"67","(60, 69]")
AGE.idents <- str_replace_all(AGE.idents,"73","(70, 79]")
AGE.idents <- str_replace_all(AGE.idents,"82","(80, 89]")
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$Age_increase <- AGE.idents
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$Age_increase))
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$Age_increase <-
  as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$Age_increase)
phenolist <- c("ConditionsH_M_S","SexF_M","Age_increase")
phenolist <- c("ConditionsH_m_S","SexF_M","Age_increase") 
str(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data[,phenolist])
BigMPLTintegrateduse.forWGCNAv2refpca1 <- ModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  traits = phenolist,
  features ="MEs",
  cor_method = "pearson",
  group.by = "PLT.idents"
)
moduleCor1 <- GetModuleTraitCorrelation(BigMPLTintegrateduse.forWGCNAv2refpca1)
names(moduleCor1)
PlotModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca1,
  label="fdr",
  label_symbol = "stars",
  text_size = 3,
  text_digits = 2,
  text_color = "black",
  #low_color = "
  #mid_color = "
  #high_color = "
  plot_max = 0.15,
  combine=T 
)+
  theme(axis.text.x = element_text(size=10),
        axis.text.y.right = element_text(size=10),
        axis.text.y.left = element_text(size=10))
BigMPLTintegrateduse.forWGCNAv2refpca2 <- ModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  traits = phenolist,
  features ="MEs",
  cor_method = "pearson",
  group.by = "Covidsplit"
)
moduleCor2 <- GetModuleTraitCorrelation(BigMPLTintegrateduse.forWGCNAv2refpca2)
names(moduleCor2)
PlotModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca2,
  label="fdr",
  label_symbol = "stars",
  text_size = 3,
  text_digits = 2,
  text_color = "black",
  #low_color = "
  #mid_color = "
  #high_color = "
  plot_max = 0.15,
  combine=T 
)+
  theme(axis.text.x = element_text(size=10),
        axis.text.y.right = element_text(size=10),
        axis.text.y.left = element_text(size=10))
table(BigMPLTintegrateduse.forWGCNAv2refpca$Covidsplit)
str(BigMPLTintegrateduse.forWGCNAv2refpca$Covidsplit)
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$Covidsplit <-
  as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$Covidsplit)
phenolist1 <- c("ConditionsH_M_S","SexF_M","Age_increase","Covidsplit")
str(BigMPLTintegrateduse.forWGCNAv2refpca@meta.data[,phenolist1])
BigMPLTintegrateduse.forWGCNAv2refpca3 <- ModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  traits = phenolist1,
  features ="MEs",
  cor_method = "pearson",
  group.by = "PLT.idents"
)
moduleCor3 <- GetModuleTraitCorrelation(BigMPLTintegrateduse.forWGCNAv2refpca3)
names(moduleCor3)
PlotModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca3,
  label="fdr",
  label_symbol = "stars",
  text_size = 3,
  text_digits = 2,
  text_color = "black",
  #low_color = "
  #mid_color = "
  #high_color = "
  plot_max = 0.15,
  combine=T 
)+
  theme(axis.text.x = element_text(size=10),
        axis.text.y.right = element_text(size=10),
        axis.text.y.left = element_text(size=10))
phenolist1 <- c("ConditionsH_M_S","SexF_M","Age_increase","Covidsplit")
BigMPLTintegrateduse.forWGCNAv2refpca3 <- ModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  traits = phenolist1,
  features ="MEs",
  cor_method = "pearson",
  group.by = "PLT.idents"
)
moduleCor3 <- GetModuleTraitCorrelation(BigMPLTintegrateduse.forWGCNAv2refpca3)
names(moduleCor3)
PlotModuleTraitCorrelation(
  BigMPLTintegrateduse.forWGCNAv2refpca3,
  label="fdr",
  label_symbol = "stars",
  text_size = 3,
  text_digits = 2,
  text_color = "black",
  #low_color = "
  #mid_color = "
  #high_color = "
  plot_max = 0.15,
  combine=T 
)+
  theme(axis.text.x = element_text(size=10),
        axis.text.y.right = element_text(size=10),
        axis.text.y.left = element_text(size=10))
BigMPLTintegrateduse.forWGCNAv2refpca <- ModuleConnectivity(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  group_name = "PLT.idents",
  corFnc = "bicor",
  corOptions = "use='p'",
  harmonized = TRUE,
  assay = "RNA",
  slot = "data",
)
sort(colnames(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$datExpr))[1:100]
str(colnames(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$datExpr))
library(tidyverse)
tempGenenames <- colnames(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$datExpr)
tempGenenames <- str_replace_all(tempGenenames,"AB-*","")
sort(tempGenenames)[1:100]
colnames(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$datExpr) <- 
  tempGenenames
sort(colnames(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$datExpr))[1:100]
str(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$gene_name)
library(tidyverse)
temp <- BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$gene_name
temp <- str_replace_all(temp,"AB-*","")
sort(temp)[1:100]
BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$gene_name_old <-
  BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$gene_name
BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$gene_name <- temp
sort(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$gene_name)[1:100]
BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$wgcna_modules$gene_name[1:100]
BigMPLTintegrateduse.forWGCNAv2refpca <- ResetModuleNames(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  new_name="PLT-M"
)
modules <- GetModules(BigMPLTintegrateduse.forWGCNAv2refpca)
print(levels(modules$module))
head(modules[,1:6])
dim(modules)
library(hdWGCNA)
library(conflicted)
conflict_prefer("select","dplyr")
p <- PlotKMEs(BigMPLTintegrateduse.forWGCNAv2refpca,
              ncol=3,
              n_hubs=20,
              text_size=2,
              plot_widths = c(3,2)
              )
p
hubgenetop100 <- GetHubGenes(BigMPLTintegrateduse.forWGCNAv2refpca,n_hubs = 100)
head(hubgenetop100)
hubgenetop100
write.table(hubgenetop100,"~/hubgenetop100BigM-anotherday2.txt")
saveRDS(BigMPLTintegrateduse.forWGCNAv2refpca,"~/PLTCovintegrateduse.forWGCNAv2TOMrefpcaMESpca.anotherday2.rds")
PLTCovintegrateduse.forWGCNAv2refpca <- ModuleExprScore(
  PLTCovintegrateduse.forWGCNAv2refpca,
  n_genes=100,
  method = "Seurat"
)
PLTCovintegrateduse.forWGCNAv2refpca <- ModuleExprScore(
  PLTCovintegrateduse.forWGCNAv2refpca,
  n_genes = 50,
  method = "UCell"
)
PLTCovMSH.cluster.forWGCNA
plot_list1 <- ModuleFeaturePlot(
  PLTCovintegrateduse.forWGCNAv2refpca,
  reduction="ref.umap",
  features = "hMEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 1,
  label_legend = F,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list1,ncol=3)
plot_list2 <- ModuleFeaturePlot(
  PLTCovMSH.cluster.forWGCNA,
  features = "scores",
  ucell=TRUE,
  point_size = 0.3,
  reduction="tsne"
)
wrap_plots(plot_list2,ncol=3)
plot_list3 <- ModuleFeaturePlot(
  PLTCovMSH.cluster.forWGCNA,
  reduction="tsne",
  features = "MEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 0.5,
  label_legend = F,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list3,ncol=3)
plot_list3 <- ModuleFeaturePlot(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  reduction="ref.umap",
  features = "MEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 0.5,
  label_legend = F,
  raster_dpi = 420,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list3,ncol=3)
BigMPLTintegrateduse.forWGCNAv2refpca <- ModuleExprScore(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  n_genes=100,
  method = "Seurat"
)
BigMPLTintegrateduse.forWGCNAv2refpca <- ModuleExprScore(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  n_genes = 50,
  method = "UCell"
)
plot_list2 <- ModuleFeaturePlot(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  features = "scores",
  ucell=TRUE,
  point_size = 0.3,
  reduction="ref.umap",
  restrict_range = F
)
wrap_plots(plot_list2,ncol=3)
BigMPLTintegrateduse.forWGCNAv2refpca <-
  readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
plot_list3 <- ModuleFeaturePlot(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  reduction="ref.umap",
  features = "MEs",
  order_points = T,
  restrict_range = T,
  point_size = 0.5,
  alpha = 0.5,
  label_legend = F,
  raster_dpi = 420,
  raster_scale = 1,
  plot_ratio = 1,
  title=T
)
wrap_plots(plot_list3,ncol=3)
ModuleCorrelogram(BigMPLTintegrateduse.forWGCNAv2refpca,
                  exclude_grey = T,
                  features ="MEs")
library(conflicted)
conflict_prefer("select","dplyr")
p <- Seurat::DotPlot(BigMPLTintegrateduse.forWGCNAv2refpca,
             features=modules$gene_name[1:10],
             group.by="PLT.idents")
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red',mid='grey95',low='blue')
p
BigMPLTintegrateduse.forWGCNAv2refpca.b <- BigMPLTintegrateduse.forWGCNAv2refpca
Idents(BigMPLTintegrateduse.forWGCNAv2refpca.b) <- BigMPLTintegrateduse.forWGCNAv2refpca.b$PLT.idents
levels(as.factor(Idents(BigMPLTintegrateduse.forWGCNAv2refpca.b)))
markers_forhere <- Seurat::FindAllMarkers(
  BigMPLTintegrateduse.forWGCNAv2refpca.b,
  only.pos=T,
  logfc.threshold = 1
) 
library(GeneOverlap)
overlap_MwithMarker <- OverlapModulesDEGs(
  BigMPLTintegrateduse.forWGCNAv2refpca.b,
  deg_df = markers_forhere,
  fc_cutoff=1
)
overlap_MwithMarker
OverlapDotPlot(
  overlap_MwithMarker,
  plot_var="odds_ratio")+
  ggtitle("Modules mapping to markers of PLT-subclusters")
saveRDS(BigMPLTintegrateduse.forWGCNAv2refpca,"~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
library(dplyr)
theme_set(theme_cowplot())
library(enrichR)
library(GeneOverlap)
set.seed(1)
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021',"GO_Molecular_Function_2021")
conflicts_prefer(BiocGenerics::intersect)
BigMPLTintegrateduse.forWGCNAv2refpca <- RunEnrichr(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  dbs=dbs,
  max_genes=100
)
str(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$enrichr_table)
enrich_df <- BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$enrichr_table
dim(enrich_df)
levels(as.factor(BigMPLTintegrateduse.forWGCNAv2refpca@misc$Clu0.7Frac0.05v2$enrichr_table$module))
write.table(enrich_df,"~/BigMFrac005v2enrich_df.txt")
EnrichrBarPlot(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  outdir="~/WGCNAenrich_plots",
  n_terms=10,
  plot_size=c(5,7),
  logscale=T
)
EnrichrDotPlot(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  mods="all",
  database = "GO_Biological_Process_2021",
  n_terms=5,
  break_ties = F
)
EnrichrDotPlot(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  mods="all",
  database = 'GO_Cellular_Component_2021',
  n_terms=5,
  break_ties = F
)
EnrichrDotPlot(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  mods="all",
  database = "GO_Molecular_Function_2021",
  n_terms=5,
  break_ties = F
)
RunEnrichrnew <- function (seurat_obj, dbs = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", 
    "GO_Molecular_Function_2021"), max_genes = 100, wait = TRUE, 
    wait_time = 5, wgcna_name = NULL, ...) 
{
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }
    CheckWGCNAName(seurat_obj, wgcna_name)
    if (!is.numeric(wait_time)) {
        stop(paste0("wait_time must be a numeric."))
    }
    if (wait_time > 60 | wait_time < 1) {
        stop(paste0("Invalid value selected for wait_time, must be greater than 0 and less than 60."))
    }
    modules <- GetModules(seurat_obj, wgcna_name)
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
    combined_output <- data.frame()
    n <- length(mods)+1
    for (i in 1:n) {
        cur_mod <- mods[i]
        if (max_genes != Inf) {
            cur_info <- subset(modules, module == cur_mod)
            cur_info <- cur_info[, c("gene_name", paste0("kME_", 
                cur_mod))]
            cur_genes <- top_n(cur_info, max_genes) %>% .$gene_name %>% 
                as.character
        }
        else {
            cur_genes <- subset(modules, module == cur_mod) %>% 
                .$gene_name %>% as.character
        }
        enriched <- enrichR::enrichr(cur_genes, dbs)
        if (wait) {
            Sys.sleep(wait_time)
        }
        for (db in names(enriched)) {
            cur_df <- enriched[[db]]
            if (nrow(cur_df) > 1) {
                cur_df$db <- db
                cur_df$module <- cur_mod
                combined_output <- rbind(combined_output, cur_df)
            }
        }
    }
    seurat_obj <- SetEnrichrTable(seurat_obj, combined_output, 
        wgcna_name)
    seurat_obj
}
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
group1 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(MSH.idents == 'S') %>% rownames
group2 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(MSH.idents != 'S') %>% rownames
group1 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(PLT.idents == 'cluC') %>% rownames
group2 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(PLT.idents != 'cluC') %>% rownames
group1 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(Covidsplit == 'use') %>% rownames
group2 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(Covidsplit != 'use') %>% rownames
head(group1)
DMEs <- FindDMEs(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
package_info("hdWGCNA")
theme_set(theme_cowplot())
PlotDMEsLollipop(
  BigMPLTintegrateduse.forWGCNAv2refpcatry, 
  DMEs, 
  wgcna_name='Clu0.7Frac0.05v2',
  pvalue = "p_val_adj"
)
ggforestplot
usethis::create_github_token()
devtools::install_local("~/Downloads/ggforestplot-master.zip")
library(ggforestplot)
wgcna_name
PlotDMEsVolcano(
  BigMPLTintegrateduse.forWGCNAv2refpcatry,
  DMEs,
)
BigMPLTintegrateduse.forWGCNAv2refpcatry <- BigMPLTintegrateduse.forWGCNAv2refpca
BigMPLTintegrateduse.forWGCNAv2refpcatry@misc$Clu0.7Frac0.05v2$wgcna_modules <- BigMPLTintegrateduse.forWGCNAv2refpcatry@misc$Clu0.7Frac0.05v2$wgcna_modules[,1:10]
dim(GetModules(BigMPLTintegrateduse.forWGCNAv2refpca, "Clu0.7Frac0.05v2"))
dim(GetModules(BigMPLTintegrateduse.forWGCNAv2refpcatry, "Clu0.7Frac0.05v2"))
str(mo2)
mo2
dim(mo2)
sum(mo2$gene_name=="")
table(is.na(mo2$gene_name))
table(is.na(mo2$module))
table(is.na(mo2$color))
table(is.na(mo2$kME_grey))
table(is.na(mo2$'kME_PLT-M1'))
table(is.na(mo2$'kME_PLT-M2'))
table(is.na(mo2$'kME_PLT-M3'))
table(is.na(mo2$'kME_PLT-M4'))
table(is.na(mo2$'kME_PLT-M5'))
table(is.na(mo2$'kME_PLT-M6'))
table(is.na(mo2))
BigMPLTintegrateduse.forWGCNAv2refpca@assays$RNA@data@x[is.na(BigMPLTintegrateduse.forWGCNAv2refpca@assays$RNA@data@x)] <- 0
levels(as.factor(PLTCovMSH.cluster.forWGCNA$cond.idents))
group1 <- PLTCovMSH.cluster.forWGCNA@meta.data %>% subset(CCClable1 == 'immPLT' & cond.idents == "S") %>% rownames
group2 <- PLTCovMSH.cluster.forWGCNA@meta.data %>% subset(CCClable1 == 'immPLT' & cond.idents != "S") %>% rownames
head(group1)
str(group1)
DMEs <- FindDMEs(
  PLTCovMSH.cluster.forWGCNA,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
package_info("hdWGCNA")
theme_set(theme_cowplot())
PlotDMEsLollipop(
  PLTCovMSH.cluster.forWGCNA, 
  DMEs, 
  wgcna_name='fourthClu1.5Frac0.05withAres05',
  pvalue = "p_val_adj"
)
wgcna_name
PlotDMEsVolcano(
  PLTCovMSH.cluster.forWGCNA,
  DMEs,
)
PlotDMEsLollipop <- function(
  seurat_obj,
  DMEs,
  wgcna_name,
  group.by=NULL,
  comparison= NULL,
  pvalue,
  avg_log2FC = 'avg_log2FC'
){
  if (!require("ggforestplot")) {
    print('Missing package: ggforestplot')
    print('Installing package: ggforestplot')
    devtools::install_github("NightingaleHealth/ggforestplot")
  }
  if(!(pvalue %in% colnames(DMEs))){
  stop('Selected pvalue is not found in DMEs dataframe column names.')
  }
  if(missing(wgcna_name) || !(wgcna_name %in% names(seurat_obj@misc))){
  stop('Please provide wgcna_name or the selected wgcna_name is not found in seurat_obj@misc.')
  }
  modules <- GetModules(seurat_obj, wgcna_name) %>% 
    subset(module != 'grey') %>% 
    mutate(module=droplevels(module))
  if (!missing(group.by) & !missing(comparison)) {
    comparisons <- comparison
    if(!(all(comparisons %in% DMEs[[group.by]]))){
    stop('Not all selected comparisons are not found in DMEs[[group.by]] or the comparison column, DMEs[[group.by]], is not correctly supplied.')
    }
    plot_list <- list()
    for(cur_comp in comparisons){
        print(cur_comp)
        cur_DMEs <- subset(DMEs, DMEs[[group.by]] == cur_comp)
        cur_title <- cur_comp
        p <- PlotLollipop(modules, cur_DMEs, pvalue, avg_log2FC = 'avg_log2FC')
        p <- p + ggtitle(cur_title) + NoLegend() +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=cur_DMEs)
        plot_list[[cur_comp]] <- p
      }
  } else if (missing(group.by) && !missing(comparison)) {
    stop('The group.by column is not provided in the DMEs data, and comparison cannot be found.')
  } else if (!missing(group.by) && missing(comparison)) {
    if (!(group.by %in% names(DMEs))) {
    stop('The group.by column is not found in the DMEs data.')
    }
    comparisons <- unique(DMEs[[group.by]])
    plot_list <- list()
    for(cur_comp in comparisons){
        print(cur_comp)
        cur_DMEs <- subset(DMEs, DMEs[[group.by]] == cur_comp)
        cur_title <- cur_comp
        p <- PlotLollipop(modules, cur_DMEs, pvalue, avg_log2FC = 'avg_log2FC')
        p <- p + ggtitle(cur_title) + NoLegend() +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=cur_DMEs)
        plot_list[[cur_comp]] <- p
      }
  } else{
    print('Please be aware comparison group/groups are not provided, which may casue an ERROR. PlotDMEsLollipop function will automatically assume all values are within the same group.')
    plot_list <- list()
    cur_DMEs <- DMEs
    p <- PlotLollipop(modules, cur_DMEs, pvalue, avg_log2FC = 'avg_log2FC')
    p <- p + NoLegend() +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=cur_DMEs)
    plot_list <- p
  }
  return(plot_list)
}
PlotLollipop <- function(
  modules,
  cur_DMEs,
  pvalue,
  avg_log2FC = 'avg_log2FC'
){
    cur_DMEs$shape <- ifelse(cur_DMEs[[pvalue]] < 0.05, 21, 4)
    cur_DMEs <- cur_DMEs %>% arrange(avg_log2FC, descending=TRUE)
    cur_DMEs$module <- factor(as.character(cur_DMEs$module), levels=as.character(cur_DMEs$module))
    n_genes <- table(modules$module)
    cur_DMEs$n_genes <- as.numeric(n_genes[as.character(cur_DMEs$module)])
    mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct
    cp <- mod_colors$color; names(cp) <- mod_colors$module
    p <- cur_DMEs %>%
    ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
    geom_vline(xintercept=0, color='black') +
    geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
    geom_point() +
    geom_point(shape=cur_DMEs$shape, color='black', fill=NA) +
    scale_color_manual(values=cp, guide='none') +
    ylab('') +
    xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=0.5, face='plain', size=10)
    )
    return(p)
}
FindDMEs <- function (seurat_obj, barcodes1, barcodes2, harmonized = TRUE, 
    wgcna_name = NULL, add_missing = FALSE, test.use = "wilcox", 
    only.pos = FALSE, logfc.threshold = 0, min.pct = 0, verbose = FALSE, 
    pseudocount.use = 0, ...) 
{
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }
    if (!(all(barcodes1 %in% colnames(seurat_obj)))) {
        stop("Some barcodes in barcodes1 not found in colnames(seurat_obj).")
    }
    if (!(all(barcodes2 %in% colnames(seurat_obj)))) {
        stop("Some barcodes in barcodes2 not found in colnames(seurat_obj).")
    }
    if (length(intersect(barcodes1, barcodes2)) > 0) {
        stop("Some barcodes overlap in barcodes1 and barcodes2")
    }
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
    MEs <- MEs[, colnames(MEs) != "grey"]
    print(dim(MEs))
    print(colnames(MEs))
    MEs[MEs < 0] <- 0
    MEs <- t(MEs)
    ME_assay <- Seurat::CreateAssayObject(MEs)
    DMEs <- FindMarkers(ME_assay, cells.1 = barcodes1, cells.2 = barcodes2, 
        slot = "counts", test.use = test.use, only.pos = only.pos, 
        logfc.threshold = logfc.threshold, min.pct = min.pct, 
        verbose = verbose, pseudocount.use = pseudocount.use, 
        ...)
    DMEs$module <- rownames(DMEs)
    if (add_missing) {
        missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
        for (cur_mod in missing_mods) {
            DMEs[cur_mod, ] <- NA
            DMEs[cur_mod, "module"] <- cur_mod
        }
    }
    DMEs
}
colnames(PLTCovMSH.cluster.forWGCNA)
MEs <- GetMEs(PLTCovMSH.cluster.forWGCNA)
MEs <- MEs[, colnames(MEs) != "grey"]
print(dim(MEs))
print(colnames(MEs))
MEs[MEs < 0] <- 0
MEs <- t(MEs)
ME_assay <- Seurat::CreateAssayObject(MEs)
DMEs <- FindMarkers(ME_assay, cells.1 = group1, cells.2 = group2, 
slot = "counts")
DMEs$module <- rownames(DMEs)
if (add_missing) {
  missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
  for (cur_mod in missing_mods) {
    DMEs[cur_mod, ] <- NA
    DMEs[cur_mod, "module"] <- cur_mod
  }
}
DMEs
clusters <- c("cluA","cluB","cluC","cluD","cluE",
              "cluF","cluG","cluH","cluI","cluJ")
DMEs2 <- data.frame()
for(cur_cluster in clusters){
  group1 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(MSH.idents == "S") %>% rownames
  group2 <- BigMPLTintegrateduse.forWGCNAv2refpca@meta.data %>% subset(MSH.idents != "S") %>% rownames
  cur_DMEs <- FindDMEs(
    BigMPLTintegrateduse.forWGCNAv2refpca,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    pseudocount.use=0.01
  )
  cur_DMEs$cluster <- cur_cluster
  DMEs2 <- rbind(DMEs2, cur_DMEs)
}
modules <- GetModules(PLTCovMSH.cluster.forWGCNA)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
plot_df <- DMEs2
plot_df$module <- factor(as.character(plot_df$module), levels=mods)
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')
p <- plot_df %>% 
  ggplot(aes(y=cluster, x=module, fill=avg_log2FC)) +
  geom_tile() 
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 
p <- p + 
  scale_fill_gradient2(low='purple', mid='black', high='yellow') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()
p
levels(as.factor(PLTCovMSH.cluster.forWGCNA$RNA_snn_res.1.5))
clusters <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
DMEs2 <- data.frame()
for(cur_cluster in clusters){
  group1 <- PLTCovMSH.cluster.forWGCNA@meta.data %>% subset(RNA_snn_res.1.5 == cur_cluster & cond.idents == "S") %>% rownames
  group2 <- PLTCovMSH.cluster.forWGCNA@meta.data %>% subset(RNA_snn_res.1.5 == cur_cluster & cond.idents != "S") %>% rownames
  cur_DMEs <- FindDMEs(
    PLTCovMSH.cluster.forWGCNA,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    pseudocount.use=0.01
  )
  cur_DMEs$cluster <- cur_cluster
  DMEs2 <- rbind(DMEs2, cur_DMEs)
}
modules <- GetModules(PLTCovMSH.cluster.forWGCNA)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
plot_df <- DMEs2
plot_df$module <- factor(as.character(plot_df$module), levels=mods)
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')
p <- plot_df %>% 
  ggplot(aes(y=cluster, x=module, fill=avg_log2FC)) +
  geom_tile() 
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 
p <- p + 
  scale_fill_gradient2(low='purple', mid='black', high='yellow') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()
p
DMEs_all <- FindAllDMEs(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  group.by = 'PLT.idents'
)
head(DMEs_all)
p <- PlotDMEsVolcano(
  BigMPLTintegrateduse.forWGCNAv2refpcatry,
  DMEs_all,
  plot_labels=T,
  show_cutoff=T
)
p + facet_wrap(~group, ncol=3)
plot_df <- DMEs_all
plot_df$module <- factor(as.character(plot_df$module), levels=mods)
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')
p <- plot_df %>%
  ggplot(aes(y=group, x=module, fill=avg_log2FC)) +
  geom_tile()
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 
p <- p + 
  scale_fill_gradient2(low='purple', mid='black', high='yellow') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()
p
DMEs_all <- FindAllDMEs(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  group.by = 'Covidsplit'
)
head(DMEs_all)
p <- PlotDMEsVolcano(
  BigMPLTintegrateduse.forWGCNAv2refpcatry,
  DMEs_all,
  plot_labels=T,
  show_cutoff=T
)
p + facet_wrap(~group, ncol=3)
plot_df <- DMEs_all
plot_df$module <- factor(as.character(plot_df$module), levels=mods)
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')
p <- plot_df %>%
  ggplot(aes(y=group, x=module, fill=avg_log2FC)) +
  geom_tile()
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 
p <- p + 
  scale_fill_gradient2(low='purple', mid='black', high='yellow') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()
p
DMEs_all <- FindAllDMEs(
  BigMPLTintegrateduse.forWGCNAv2refpca,
  group.by = 'MSH.idents'
)
head(DMEs_all)
p <- PlotDMEsVolcano(
  BigMPLTintegrateduse.forWGCNAv2refpcatry,
  DMEs_all,
  plot_labels=T,
  show_cutoff=T
)
p + facet_wrap(~group, ncol=3)
plot_df <- DMEs_all
plot_df$module <- factor(as.character(plot_df$module), levels=mods)
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')
p <- plot_df %>%
  ggplot(aes(y=group, x=module, fill=avg_log2FC)) +
  geom_tile()
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 
p <- p + 
  scale_fill_gradient2(low='purple', mid='black', high='yellow') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()
p

if(!require(DDRTree))install.packages('DDRTree')
if(!require(pheatmap))install.packages('pheatmap')
if(!require(Seurat))install.packages('Seurat')
if(!require(dplyr))install.packages('dplyr')
library(monocle)
packageVersion('monocle')
devtools::load_all("~/R/x86_64-pc-linux-gnu-library/4.1/monocle")
package.version("monocle")

.libPaths()
library(monocle)
suppressMessages({
  suppressWarnings({library(Seurat)
library(monocle)})
})
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
BigMerge <- readRDS("~/BigMergeProjectedAll.rds")
BigMPLT = BigMerge[,BigMerge$predicted.celltypewithnew 
                   %in% c("Platelet")]
BigMPLT$PLT.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$PLT.idents
BigMPLT$MSH.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents
BigMPLT$SEX <- BigMPLTintegrateduse.forWGCNAv2refpca$SEX
BigMPLT$Age_increase <- BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase
BigMPLT$ID <- BigMPLTintegrateduse.forWGCNAv2refpca$ID
BigMPLT$Covidsplit <- BigMPLTintegrateduse.forWGCNAv2refpca$Covidsplit
BigMPLT$studyID <- BigMPLTintegrateduse.forWGCNAv2refpca$studyID
BigMPLT$predicted.celltypePLTCCClable1 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTCCClable1
BigMPLT$predicted.celltypePLTres0.7 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTres0.7
BigMPLT@misc <-
  BigMPLTintegrateduse.forWGCNAv2refpca@misc
PLTforCT <- BigMPLT
data <- as(as.matrix(PLTforCT@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = PLTforCT@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=10, relative_expr = TRUE)
premycds <- mycds
mycds <- premycds
saveRDS(premycds,"~/premycds0.rds")
      mycds <- detectGenes(mycds, min_expr = 0.1)
      mycds_expressed_genes <- row.names(subset(fData(mycds),
                                                num_cells_expressed >= 5))
      mycds_filteredLow <- mycds[mycds_expressed_genes,]
      disp_table <- dispersionTable(mycds_filteredLow)
      disp.genes <- subset(disp_table, mean_expression >= 0.1 & 
                             dispersion_empirical >= 1 * dispersion_fit)
      length(disp.genes$gene_id)
mycds <- setOrderingFilter(mycds, disp.genes$gene_id)
  plot_ordering_genes(mycds)
  plot_pc_variance_explained(mycds, return_all = F)
mycds <- reduceDimension(mycds, max_components = 2, num_dim = 9, 
                         reduction_method = 'DDRTree')
mycds <- reduceDimension(mycds, max_components = 3, num_dim = 9, 
                         reduction_method = 'DDRTree')
suppressWarnings(mycds <- orderCells(mycds))
saveRDS(mycds,"~/noNormmycdsDDRTree2.rds")
saveRDS(mycds,"~/noNormmycdsDDRTree3.rds")
saveRDS(mycds,"~/noNormmycdsDDRTree4.rds")
library(patchwork)
library(Biobase)
library(multtest)
pData(mycds)
pData
e <- plot_cell_trajectory(mycds, color_by = "PLT.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "MSH.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.3)
e$layers[[1]]$aes_params$alpha <- 0.5
e
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~cond.idents, nrow = 1)
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~RNA_snn_res.1.5, nrow = 5)
mycds <- readRDS("~/noNormmycdsDDRTree3.rds")
library(RColorBrewer)
display.brewer.pal(11,"Set3")
mycol = brewer.pal(11,"Set3")
mycol
newcol = brewer.pal(12,"Set3")
newcol
mycol=c("
        "
        "
cols = c(pal_npg()(8),"
         "
         "
plot_cell_trajectory(mycds, color_by = c("State"),
                     show_branch_points = T,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~Covidsplit, nrow = 5)+
  NoLegend()
plot_cell_trajectory(mycds, color_by = c("State"),
                     show_branch_points = T,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~disease.idents, nrow = 5)+
  NoLegend()
data_df <- t(reducedDimS(mycds))%>%
  as.data.frame()%>%
  select_('Component 1' = 1,'Component 2' = 2)%>%
  rownames_to_column("Cells")%>%
  mutate(pData(mycds)$State,pData(mycds)$Pseudotime,
         pData(mycds)$ID,pData(mycds)$PLT.idents,
         pData(mycds)$disease.idents)
head(data_df)
colnames(data_df) <- c("Cells","Component_1",
                       "Component_2","State",
                       "Pseudotime","ID",
                       "PLT.idents","disease.idents")
reduced_dim_coords <- reducedDimK(mycds)
ica_space_df <- Matrix::t(reduced_dim_coords)%>%
  as.data.frame()%>%
  select_(prin_graph_dim_1=1,prin_graph_dim_2=2)%>%
  mutate(sample_name=rownames(.),sample_state=rownames(.))
dp_mst <- minSpanningTree(mycds)
edge_df <- dp_mst%>%
  igraph::as_data_frame()%>%
  select_(source="from",target="to")%>%
  left_join(ica_space_df%>%select_(source="sample_name",
                                   source_prin_graph_dim_1="prin_graph_dim_1",
                                   source_prin_graph_dim_2="prin_graph_dim_2"),
            by="source")%>%
  left_join(ica_space_df%>%select_(target="sample_name",
                                   target_prin_graph_dim_1="prin_graph_dim_1",
                                   target_prin_graph_dim_2="prin_graph_dim_2"),
            by="target")
str(edge_df)
write.csv(edge_df,"~/edge_df_diseasenew.csv")
ggplot()+
  geom_point_rast(data=data_df,aes(x=Component_1,
                                   y=Component_2,
                                   color=Pseudotime))+
  scale_color_viridis()+
  theme_bw()+
  theme_dr(arrow=grid::arrow(length=unit(0,"inches")))+
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.8,"line"),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    axis.title = element_text(size=15)
  )
plot_cell_trajectory(mycds, color_by = c("State"),
                     show_branch_points = F,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~PLT.idents, nrow = 5)+
  NoLegend()
library(ggpubr)
library(ggsignif)
library(patchwork)
library(tidydr)
library(tidyr)
library(ggpubr)
library(ggforce)
library(ggrastr)
library(viridis)
library(tidyverse)
data_df <- t(reducedDimS(mycds))%>%
  as.data.frame()%>%
  select_('Component 1' = 1,'Component 2' = 2)%>%
  rownames_to_column("Cells")%>%
  mutate(pData(mycds)$State,pData(mycds)$Pseudotime,
         pData(mycds)$ID,pData(mycds)$PLT.idents,
         pData(mycds)$MSH.idents)
head(data_df)
colnames(data_df) <- c("Cells","Component_1",
                       "Component_2","State",
                       "Pseudotime","ID",
                       "PLT.idents","MSH.idents")
reduced_dim_coords <- reducedDimK(mycds)
ica_space_df <- Matrix::t(reduced_dim_coords)%>%
  as.data.frame()%>%
  select_(prin_graph_dim_1=1,prin_graph_dim_2=2)%>%
  mutate(sample_name=rownames(.),sample_state=rownames(.))
dp_mst <- minSpanningTree(mycds)
edge_df <- dp_mst%>%
  igraph::as_data_frame()%>%
  select_(source="from",target="to")%>%
  left_join(ica_space_df%>%select_(source="sample_name",
                                   source_prin_graph_dim_1="prin_graph_dim_1",
                                   source_prin_graph_dim_2="prin_graph_dim_2"),
            by="source")%>%
  left_join(ica_space_df%>%select_(target="sample_name",
                                   target_prin_graph_dim_1="prin_graph_dim_1",
                                   target_prin_graph_dim_2="prin_graph_dim_2"),
            by="target")
str(edge_df)
write.csv(edge_df,"~/edge_df.csv")
getwd()
ggplot()+
  geom_point_rast(data=data_df,aes(x=Component_1,
                                   y=Component_2,
                                   color=Pseudotime))+
  scale_color_viridis()+
  theme_bw()+
  theme_dr(arrow=grid::arrow(length=unit(0,"inches")))+
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.8,"line"),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    axis.title = element_text(size=15)
  )
plot_cell_trajectory(mycds, color_by = c("State"),
                          show_backbone = F,show_state_number = F,
                          show_tree = T,show_branch_points = F) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~MSH.idents, nrow = 5)+
  NoLegend()
plot_cell_trajectory(mycds, color_by = c("State"),
                     show_branch_points = T,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~SEX, nrow = 5)+
  NoLegend()
suppressMessages({
  suppressWarnings({library(Seurat)
library(monocle)})
})
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
BigMerge <- readRDS("~/BigMergeProjectedAll.rds")
BigMPLT = BigMerge[,BigMerge$predicted.celltypewithnew 
                   %in% c("Platelet")]
BigMPLT = BigMPLT[,BigMPLT$Covidsplit
                   %in% c("use")]
BigMPLT$PLT.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$PLT.idents
BigMPLT$MSH.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents
BigMPLT$SEX <- BigMPLTintegrateduse.forWGCNAv2refpca$SEX
BigMPLT$Age_increase <- BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase
BigMPLT$ID <- BigMPLTintegrateduse.forWGCNAv2refpca$ID
BigMPLT$Covidsplit <- BigMPLTintegrateduse.forWGCNAv2refpca$Covidsplit
BigMPLT$studyID <- BigMPLTintegrateduse.forWGCNAv2refpca$studyID
BigMPLT$predicted.celltypePLTCCClable1 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTCCClable1
BigMPLT$predicted.celltypePLTres0.7 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTres0.7
BigMPLT@misc <-
  BigMPLTintegrateduse.forWGCNAv2refpca@misc
PLTforCT <- BigMPLT
data <- as(as.matrix(PLTforCT@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = PLTforCT@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=10, relative_expr = TRUE)
premycds <- mycds
mycds <- premycds
saveRDS(premycds,"~/premycds.rds")
      mycds <- detectGenes(mycds, min_expr = 0.1)
      mycds_expressed_genes <- row.names(subset(fData(mycds),
                                                num_cells_expressed >= 5))
      mycds_filteredLow <- mycds[mycds_expressed_genes,]
      disp_table <- dispersionTable(mycds_filteredLow)
      disp.genes <- subset(disp_table, mean_expression >= 0.1 & 
                             dispersion_empirical >= 1 * dispersion_fit)
      length(disp.genes$gene_id)
mycds <- setOrderingFilter(mycds, disp.genes$gene_id)
  plot_ordering_genes(mycds)
  plot_pc_variance_explained(mycds, return_all = F)
mycds <- reduceDimension(mycds, max_components = 2, num_dim = 9, 
                         reduction_method = 'DDRTree')
mycds <- reduceDimension(mycds, max_components = 3, num_dim = 9, 
                         reduction_method = 'DDRTree')
getwd()
saveRDS(mycds,"BigMergeCovid/mycds2.1.rds")
mycds <- readRDS("BigMergeCovid/mycds2.1.rds")
suppressWarnings(mycds <- orderCells(mycds))
saveRDS(mycds,"~/noNormmycdsDDRTree2.rds")
saveRDS(mycds,"~/noNormmycdsDDRTree3onlycovid.rds")
saveRDS(mycds,"~/noNormmycdsDDRTree4.rds")
library(patchwork)
library(Biobase)
library(multtest)
pData(mycds)
pData
e <- plot_cell_trajectory(mycds, color_by = "PLT.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "MSH.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.3)
e$layers[[1]]$aes_params$alpha <- 0.5
e
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~cond.idents, nrow = 1)
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~RNA_snn_res.1.5, nrow = 5)
library(RColorBrewer)
display.brewer.pal(11,"Set3")
mycol = brewer.pal(11,"Set3")
mycol
newcol = brewer.pal(12,"Set3")
newcol
mycol=c("
        "
        "
cols = c(pal_npg()(8),"
         "
         "
plot_cell_trajectory(mycds, color_by = c("RNA_snn_res.0.7"),
                     show_branch_points = T,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~RNA_snn_res.0.7, nrow = 5)+
  NoLegend()
library(ggpubr)
library(ggsignif)
library(patchwork)
library(tidydr)
library(tidyr)
library(ggpubr)
library(ggforce)
library(ggrastr)
library(viridis)
library(tidyverse)
data_df <- readRDS("/home/jin/ws/covid2/data_dfusing4thmycds.rds")
basic <- ggplot()+
  geom_point_rast(data=data_df,aes(x=Component_1,
                                   y=Component_2,
                                   color=Pseudotime))+
  theme_bw()+
  theme_dr(arrow=grid::arrow(length=unit(0,"inches")))+
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.8,"line"),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    axis.title = element_text(size=15)
  )
basic
Cellratio <- prop.table(table(data_df$State,data_df$PLTclusters),margin=2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("State","PLT.idents","Freq")
ggplot()+
  geom_point_rast(data=data_df,aes(x=Component_1,
                                   y=Component_2,
                                   color=Pseudotime))+
  theme_bw()+
  theme_dr(arrow=grid::arrow(length=unit(0,"inches")))+
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.8,"line"),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    axis.title = element_text(size=15)
  )+
  geom_arc_bar(data=subset(Cellratio,State=='1'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))+
  geom_arc_bar(data=subset(Cellratio,State=='2'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))+
  geom_arc_bar(data=subset(Cellratio,State=='3'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))+
  geom_arc_bar(data=subset(Cellratio,State=='4'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))+
  geom_arc_bar(data=subset(Cellratio,State=='5'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))+
  geom_arc_bar(data=subset(Cellratio,State=='6'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))+
  geom_arc_bar(data=subset(Cellratio,State=='7'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))
write.csv(Cellratio,"Cellratio.csv")
basic+
  geom_arc_bar(data=subset(Cellratio,State=='9'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))
table(mycds$State)
basic+
  geom_arc_bar(data=subset(Cellratio,State=='9'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))+
  scale_colour_manual(values=cols)
basic+
  geom_arc_bar(data=subset(Cellratio,State=='9'),stat="pie",
               aes(x0=2,y0=3.5,r0=0,r=0.8,amount=Freq,fill=PLT.idents))
plot_cell_trajectory(mycds, color_by = c("State"),
                          show_backbone = F,show_state_number = F,
                          show_tree = T,show_branch_points = F) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~MSH.idents, nrow = 5)+
  NoLegend()
plot_cell_trajectory(mycds, color_by = c("State"),
                     show_branch_points = T,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~SEX, nrow = 5)+
  NoLegend()
library(ggthemes)
library(ggsci)
cols = c(pal_npg()(8),"
         "
         "
mycol=c("
        "
        "
ggplot(Cellratio,aes(State,Freq,fill=PLT.idents))+
  geom_bar(stat="identity",position="fill",width = 0.9,just=0.5)+
  scale_fill_manual(values=cols)+
  theme_bw()+
  theme_few()+
  theme_classic()+
  theme(axis.ticks.length=unit(0.01,'cm'))+
  theme(legend.position = "left",legend.text = element_text(size=90),
        legend.title = element_blank(),legend.spacing  = unit(0.2,"cm"), 
        legend.key.size=unit(4,"cm"))+
  theme(legend.box.spacing = unit(0.7,"cm"))+
  theme(axis.text=element_text(size=85,color='black'),
        axis.title.y = element_text(size=85))+
  theme(axis.title.x = element_blank(),axis.line.x = element_blank())
plot_cell_trajectory(mycds, color_by = c("State"),
                     show_branch_points = T,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  NoLegend()
Cellratio$PLT.idents <- factor(Cellratio$PLT.idents,levels = c("8","9","5","2","3","4","1","6","7"))
ggplot(Cellratio,aes(PLT.idents,Freq,fill=State))+
  geom_bar(stat="identity",position="fill",width = 0.9,just=0.5)+
  scale_fill_brewer("blues")+
  theme_bw()+
  theme_few()+
  theme_classic()+
  theme(axis.ticks.length=unit(0.01,'cm'))+
  theme(legend.position = "left",legend.text = element_text(size=90),
        legend.title = element_blank(),legend.spacing  = unit(0.2,"cm"), 
        legend.key.size=unit(4,"cm"))+
  theme(legend.box.spacing = unit(0.7,"cm"))+
  theme(axis.text=element_text(size=85,color='black'),
        axis.title.y = element_text(size=85))+
  theme(axis.title.x = element_blank(),axis.line.x = element_blank()))
mycds <- readRDS("~/noNormmycdsDDRTree3.rds")
data_df <- t(reducedDimS(mycds))%>%
  as.data.frame()%>%
  select_('Component 1' = 1,'Component 2' = 2)%>%
  rownames_to_column("Cells")%>%
  mutate(pData(mycds)$State,pData(mycds)$Pseudotime,
         pData(mycds)$ID,pData(mycds)$cond.idents,
         pData(mycds)$disease.idents)
head(data_df)
colnames(data_df) <- c("Cells","Component_1",
                       "Component_2","State",
                       "Pseudotime","ID",
                       "cond.idents","disease.idents")
reduced_dim_coords <- reducedDimK(mycds)
ica_space_df <- Matrix::t(reduced_dim_coords)%>%
  as.data.frame()%>%
  select_(prin_graph_dim_1=1,prin_graph_dim_2=2)%>%
  mutate(sample_name=rownames(.),sample_state=rownames(.))
dp_mst <- minSpanningTree(mycds)
edge_df <- dp_mst%>%
  igraph::as_data_frame()%>%
  select_(source="from",target="to")%>%
  left_join(ica_space_df%>%select_(source="sample_name",
                                   source_prin_graph_dim_1="prin_graph_dim_1",
                                   source_prin_graph_dim_2="prin_graph_dim_2"),
            by="source")%>%
  left_join(ica_space_df%>%select_(target="sample_name",
                                   target_prin_graph_dim_1="prin_graph_dim_1",
                                   target_prin_graph_dim_2="prin_graph_dim_2"),
            by="target")
str(edge_df)
write.csv(edge_df,"~/edge_df.csv")
plot_cell_trajectory(mycds, color_by = c("State"),
                     show_branch_points = T,show_backbone = T) +
  scale_colour_manual(values=mycol)+
  facet_wrap(~disease.idents, nrow = 5)+
  NoLegend()
Cellratio <- prop.table(table(data_df$State,data_df$disease.idents),margin=2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("State","disease.idents","Freq")
library(ggthemes)
library(ggsci)
ggplot(Cellratio,aes(disease.idents,Freq,fill=State))+
  geom_bar(stat="identity",position="fill",width = 0.9,just=0.5)+
  scale_fill_manual(values=mycol)+
  theme_bw()+
  theme_few()+
  theme_classic()+
  theme(axis.ticks.length=unit(0.01,'cm'))+
   theme(legend.position = "none",legend.text = element_text(size=90),
        legend.title = element_blank(),legend.spacing  = unit(0.2,"cm"),
        legend.key.size=unit(4,"cm"))+
  theme(legend.box.spacing = unit(0.7,"cm"))+
  theme(axis.text=element_text(size=85,color='black'),
        axis.title.y = element_text(size=85))+
  theme(axis.title.x = element_blank(),axis.line.x = element_blank())
suppressMessages({
  suppressWarnings({library(Seurat)
library(monocle)})
})
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
BigMerge <- readRDS("~/BigMergeProjectedAll.rds")
BigMPLT = BigMerge[,BigMerge$predicted.celltypewithnew 
                   %in% c("Platelet")]
BigMPLT$PLT.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$PLT.idents
BigMPLT$MSH.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents
BigMPLT$SEX <- BigMPLTintegrateduse.forWGCNAv2refpca$SEX
BigMPLT$Age_increase <- BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase
BigMPLT$ID <- BigMPLTintegrateduse.forWGCNAv2refpca$ID
BigMPLT$Covidsplit <- BigMPLTintegrateduse.forWGCNAv2refpca$Covidsplit
BigMPLT$studyID <- BigMPLTintegrateduse.forWGCNAv2refpca$studyID
BigMPLT$predicted.celltypePLTCCClable1 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTCCClable1
BigMPLT$predicted.celltypePLTres0.7 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTres0.7
BigMPLT@misc <-
  BigMPLTintegrateduse.forWGCNAv2refpca@misc
PLTforCT <- BigMPLT
data <- as(as.matrix(PLTforCT@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = PLTforCT@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=10, relative_expr = TRUE)
premycds <- mycds
mycds <- premycds
mycds <- detectGenes(mycds, min_expr = 0.1)
mycds_expressed_genes <- row.names(subset(fData(mycds),
                                                num_cells_expressed >= 5))
length(mycds_expressed_genes)
mycds_filteredLow <- mycds[mycds_expressed_genes,]
disp_table <- dispersionTable(mycds_filteredLow)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
length(disp.genes$gene_id)
mycds_filteredLowDisper <- mycds_filteredLow[disp.genes$gene_id,]
set.seed(1)
diff_test_res <- differentialGeneTest(mycds_filteredLowDisper,fullModelFormulaStr = "~PLT.idents")
write.csv(diff_test_res,"~/CT_diff_test_res_for_PLTident.csv")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
length(ordering_genes)
mycds <- setOrderingFilter(mycds, ordering_genes)
  plot_ordering_genes(mycds)
  plot_pc_variance_explained(mycds, return_all = F)
mycds <- reduceDimension(mycds, max_components = 2, num_dim = 9, 
                         reduction_method = 'DDRTree')
suppressWarnings(mycds <- orderCells(mycds))
saveRDS(mycds,"~/noNormmycdsDDTreeusingPLTDE.rds")
library(patchwork)
library(Biobase)
library(multtest)
pData(mycds)
pData
e <- plot_cell_trajectory(mycds, color_by = "PLT.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "MSH.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.3)
e$layers[[1]]$aes_params$alpha <- 0.5
e
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~cond.idents, nrow = 1)
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~RNA_snn_res.1.5, nrow = 5)
suppressMessages({
  suppressWarnings({library(Seurat)
library(monocle)})
})
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
BigMerge <- readRDS("~/BigMergeProjectedAll.rds")
BigMPLT = BigMerge[,BigMerge$predicted.celltypewithnew 
                   %in% c("Platelet")]
BigMPLT$PLT.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$PLT.idents
BigMPLT$MSH.idents <- BigMPLTintegrateduse.forWGCNAv2refpca$MSH.idents
BigMPLT$SEX <- BigMPLTintegrateduse.forWGCNAv2refpca$SEX
BigMPLT$Age_increase <- BigMPLTintegrateduse.forWGCNAv2refpca$Age_increase
BigMPLT$ID <- BigMPLTintegrateduse.forWGCNAv2refpca$ID
BigMPLT$Covidsplit <- BigMPLTintegrateduse.forWGCNAv2refpca$Covidsplit
BigMPLT$studyID <- BigMPLTintegrateduse.forWGCNAv2refpca$studyID
BigMPLT$predicted.celltypePLTCCClable1 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTCCClable1
BigMPLT$predicted.celltypePLTres0.7 <-
  BigMPLTintegrateduse.forWGCNAv2refpca$predicted.celltypePLTres0.7
BigMPLT@misc <-
  BigMPLTintegrateduse.forWGCNAv2refpca@misc
PLTforCT <- BigMPLT
data <- as(as.matrix(PLTforCT@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = PLTforCT@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=10, relative_expr = TRUE)
premycds <- mycds
mycds <- premycds
mycds <- detectGenes(mycds, min_expr = 0.1)
mycds_expressed_genes <- row.names(subset(fData(mycds),num_cells_expressed >= 5))
length(mycds_expressed_genes)
mycds_filteredLow <- mycds[mycds_expressed_genes,]
    pData(mycds_filteredLow)
    mycds_filteredLow@phenoData@data
    table(is.na(mycds_filteredLow@phenoData@data))
    mycds_filteredLow@phenoData@data[is.na(mycds_filteredLow@phenoData@data)] <- 0
    table(is.na(mycds_filteredLow@phenoData@data))
disp_table <- dispersionTable(mycds_filteredLow)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
length(disp.genes$gene_id)
mycds_filteredLowDisper <- mycds_filteredLow[disp.genes$gene_id,]
set.seed(1)
diff_test_res <- differentialGeneTest(mycds_filteredLowDisper,fullModelFormulaStr = "~MSH.idents")
mycds_filteredLowDisper@phenoData@data$MSH.idents[is.na(mycds_filteredLowDisper@phenoData@data$MSH.idents)] <- "nan"
write.csv(diff_test_res,"~/CT_diff_test_res_for_MSHident.csv")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
length(ordering_genes)
mycds <- setOrderingFilter(mycds, ordering_genes)
  plot_ordering_genes(mycds)
  plot_pc_variance_explained(mycds, return_all = F)
mycds <- reduceDimension(mycds, max_components = 2, num_dim = 9, 
                         reduction_method = 'DDRTree')
suppressWarnings(mycds <- orderCells(mycds))
saveRDS(mycds,"~/noNormmycdsDDTreeusingMSHDE.rds")
mycds <- readRDS("~/noNormmycdsDDTreeusingPLTDE.rds")
library(patchwork)
library(Biobase)
library(multtest)
pData(mycds)
pData
e <- plot_cell_trajectory(mycds, color_by = "PLT.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "MSH.idents",cell_size = 0.3)|
  plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.3)
e$layers[[1]]$aes_params$alpha <- 0.5
e
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~cond.idents, nrow = 1)
plot_cell_trajectory(mycds, color_by = c("cond.idents")) +
  facet_wrap(~RNA_snn_res.1.5, nrow = 5)
BEAMplotdata=data.frame(BEAM_res1=BEAM_geme1,BEAM_res2=BEAM_geme2,BEAM_res3=BEAM_geme3,BEAM_res4=BEAM_geme4,BEAM_res5=BEAM_geme5)
BEAMplotdata=cbind(BEAM_res1=BEAM_geme1,BEAM_res2=BEAM_geme2,BEAM_res3=BEAM_geme3,BEAM_res4=BEAM_geme4,BEAM_res5=BEAM_geme5)
BEAM_geme1 <- write.table(BEAM_geme1,"BEAM_geme1.txt")
BEAM_geme2 <- write.table(BEAM_geme2,"BEAM_geme2.txt")
BEAM_geme3 <- write.table(BEAM_geme3,"BEAM_geme3.txt")
BEAM_geme4 <- write.table(BEAM_geme4,"BEAM_geme4.txt")
BEAM_geme5 <- write.table(BEAM_geme5,"BEAM_geme5.txt")
plot_cell_trajectory(mycds_filteredLowDisp,color_by = 'Pseudotime')
BEAM_res1 <- BEAM(mycds_filteredLowDisp, branch_point = 1, cores = 10)
BEAM_res1 <- BEAM(mycds, branch_point = 1, cores = 10)
BEAM_res2 <- BEAM(mycds_filteredLowDisp, branch_point = 2, cores = 10)
mycds_filteredLowDisp <- readRDS("~/noNormmycdsDDRTree_filteredLowDisprerunOrdercell.rds")
mycds_filteredLowDisp@minSpanningTree
mycds@minSpanningTree
mycds3@minSpanningTree
str(BEAM_res)
BEAM_res2 <- readRDS("~/BEAM_res2.rds")
BEAM_res5 <- readRDS("~/BEAM_res5.rds")
BEAM_res1 <- BEAM_res1[order(BEAM_res1$qval),]
BEAM_res1 <- BEAM_res1[,c("gene_short_name", "pval", "qval")]
BEAM_geme1 <- row.names(subset(BEAM_res1,qval < 0.01))
str(BEAM_geme1)
length(BEAM_geme1)
BEAM_res2 <- BEAM_res2[order(BEAM_res2$qval),]
BEAM_res2 <- BEAM_res2[,c("gene_short_name", "pval", "qval")]
BEAM_geme2 <- row.names(subset(BEAM_res2,qval < 0.01))
str(BEAM_geme2)
length(BEAM_geme2)
BEAM_res5 <- BEAM_res5[order(BEAM_res5$qval),]
BEAM_res5 <- BEAM_res5[,c("gene_short_name", "pval", "qval")]
BEAM_geme5 <- row.names(subset(BEAM_res5,qval < 0.01))
str(BEAM_geme5)
length(BEAM_geme5)
plot_genes_branched_heatmap(mycds[BEAM_geme1[1:50],],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme2[1:50],],branch_point = 2,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme3[1:50],],branch_point = 3,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme4[1:50],],branch_point = 4,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme5[1:50],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme1[1:100],],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme2[1:100],],branch_point = 2,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme3[1:100],],branch_point = 3,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme4[1:100],],branch_point = 4,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme5[1:100],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds_filteredLowDisp[BEAM_geme5[1:50],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme1[1:100],],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme2[1:100],],branch_point = 2,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme3[1:100],],branch_point = 3,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme4[1:100],],branch_point = 4,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme5[1:200],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme5,],branch_point = 5,num_clusters = 8,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds.unsup[BEAM_geme[51:100],],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
mycds <- readRDS("~/noNormmycdsDDRTree3rerunOrdercell.rds")
mycds_filteredLowDisp <- readRDS("~/noNormmycdsDDRTree_filteredLowDisprerunOrdercell.rds")
BEAM_res1 <- BEAM(mycds_filteredLowDisp, branch_point = 1, cores = 10)
saveRDS(BEAM_res1,"~/BEAM_res1.rds")
BEAM_res2 <- BEAM(mycds_filteredLowDisp, branch_point = 1, cores = 10)
saveRDS(BEAM_res2,"~/BEAM_res2.rds")
BEAM_res3 <- BEAM(mycds_filteredLowDisp, branch_point = 1, cores = 10)
saveRDS(BEAM_res3,"~/BEAM_res3.rds")
BEAM_res4 <- BEAM(mycds_filteredLowDisp, branch_point = 1, cores = 10)
saveRDS(BEAM_res4,"~/BEAM_res4.rds")
BEAM_res5 <- BEAM(mycds_filteredLowDisp, branch_point = 1, cores = 10)
saveRDS(BEAM_res5,"~/BEAM_res5.rds")
BEAM_res3 <- readRDS("~/BEAM_res3.rds")
BEAM_res5 <- readRDS("~/BEAM_res5.rds")
BEAM_res1 <- readRDS("~/BEAM_res1.rds")
BEAM_res2 <- readRDS("~/BEAM_res2.rds")
BEAM_res4 <- readRDS("~/BEAM_res4.rds")
BEAM_res3 <- BEAM_res3[order(BEAM_res3$qval),]
BEAM_res3 <- BEAM_res3[,c("gene_short_name", "pval", "qval")]
BEAM_geme3 <- row.names(subset(BEAM_res3,qval < 0.01))
str(BEAM_geme3)
length(BEAM_geme3)
BEAM_res5 <- BEAM_res5[order(BEAM_res5$qval),]
BEAM_res5 <- BEAM_res5[,c("gene_short_name", "pval", "qval")]
BEAM_geme5 <- row.names(subset(BEAM_res5,qval < 0.01))
str(BEAM_geme5)
length(BEAM_geme5)
BEAM_res1 <- BEAM_res1[order(BEAM_res1$qval),]
BEAM_res1 <- BEAM_res1[,c("gene_short_name", "pval", "qval")]
BEAM_geme1 <- row.names(subset(BEAM_res1,qval < 0.01))
str(BEAM_geme1)
length(BEAM_geme1)
BEAM_res2 <- BEAM_res2[order(BEAM_res2$qval),]
BEAM_res2 <- BEAM_res2[,c("gene_short_name", "pval", "qval")]
BEAM_geme2 <- row.names(subset(BEAM_res2,qval < 0.01))
str(BEAM_geme2)
length(BEAM_geme2)
BEAM_res4 <- BEAM_res4[order(BEAM_res4$qval),]
BEAM_res4 <- BEAM_res4[,c("gene_short_name", "pval", "qval")]
BEAM_geme4 <- row.names(subset(BEAM_res4,qval < 0.01))
str(BEAM_geme4)
length(BEAM_geme4)
plot_genes_branched_heatmap(mycds[BEAM_geme2[1:50],],branch_point = 2,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme3[1:50],],branch_point = 3,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme4[1:50],],branch_point = 4,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme5[1:50],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme1[1:100],],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme2[1:100],],branch_point = 2,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme3[1:100],],branch_point = 3,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme4[1:100],],branch_point = 4,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme5[1:100],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds_filteredLowDisp[BEAM_geme5[1:50],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme1[1:100],],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme2[1:100],],branch_point = 2,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme3[1:100],],branch_point = 3,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme4[1:100],],branch_point = 4,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
plot_genes_branched_heatmap(mycds[BEAM_geme5[1:200],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
BEAMHmapBP5All <- plot_genes_branched_heatmap(mycds[BEAM_geme5,],branch_point = 5,num_clusters = 8,cores = 1,use_gene_short_name = T,show_rownames = T,return_heatmap=T)
BEAMHmapBP3All <- plot_genes_branched_heatmap(mycds[BEAM_geme3,],branch_point = 3,num_clusters = 8,cores = 1,use_gene_short_name = T,show_rownames = T,return_heatmap=T)
BEAMHmapBP5t100 <- plot_genes_branched_heatmap(mycds[BEAM_geme5[1:100],],branch_point = 5,num_clusters = 5,cores = 1,use_gene_short_name = T,show_rownames = T,return_heatmap=T)
BEAMHmapBP3t100 <- plot_genes_branched_heatmap(mycds[BEAM_geme3[1:100],],branch_point = 3,num_clusters = 5,cores = 1,use_gene_short_name = T,show_rownames = T,return_heatmap=T)
BEAMHmapBP5t50 <- plot_genes_branched_heatmap(mycds[BEAM_geme5[1:50],],branch_point = 5,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T,return_heatmap=T)
BEAMHmapBP3t50 <- plot_genes_branched_heatmap(mycds[BEAM_geme3[1:50],],branch_point = 3,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T,return_heatmap=T)
saveRDS(BEAMHmapBP3t50,"~/BEAMHmapBP3t50.rds")
saveRDS(BEAMHmapBP5t50,"~/BEAMHmapBP5t50.rds")
saveRDS(BEAMHmapBP3t100,"~/BEAMHmapBP3t100.rds")
saveRDS(BEAMHmapBP5t100,"~/BEAMHmapBP5t100.rds")
saveRDS(BEAMHmapBP3All,"~/BEAMHmapBP3All.rds")
saveRDS(BEAMHmapBP5All,"~/BEAMHmapBP5All.rds")
BEAMHmapBP3t100 <- readRDS("~/BEAMHmapBP3t100.rds")
BEAMHmapBP3All <- readRDS("~/BEAMHmapBP3All.rds")
hclusters <- cutree(BEAMHmapBP3All$ph$tree_row,k=8)
range(hclusters)
hclusters <- data.frame(hclusters)
hclusters[,1] <- as.character(hclusters[,1])
colnames(hclusters) <- "Genes_clusters"
table(hclusters)
Hieclu1 <- subset(hclusters,Genes_clusters==1)
dim(Hieclu1)
Hieclu1genes <- rownames(Hieclu1)
Hieclu1genes
Hieclu3 <- subset(hclusters,Genes_clusters==3)
dim(Hieclu3)
Hieclu3genes <- rownames(Hieclu3)
Hieclu3genes
Hieclu7 <- subset(hclusters,Genes_clusters==7)
dim(Hieclu7)
Hieclu7genes <- rownames(Hieclu7)
Hieclu7genes
hclusters <- cutree(BEAMHmapBP3t50$ph$tree_row,k=4)
range(hclusters)
hclusters <- data.frame(hclusters)
hclusters[,1] <- as.character(hclusters[,1])
colnames(hclusters) <- "Genes_clusters"
table(hclusters)
Hieclu1 <- subset(hclusters,Genes_clusters==1)
dim(Hieclu1)
Hieclu1genes <- rownames(Hieclu1)
Hieclu1genes
Hieclu3 <- subset(hclusters,Genes_clusters==3)
dim(Hieclu3)
Hieclu3genes <- rownames(Hieclu3)
Hieclu3genes
Hieclu2 <- subset(hclusters,Genes_clusters==2)
dim(Hieclu2)
Hieclu2genes <- rownames(Hieclu2)
Hieclu2genes
hclusters <- cutree(BEAMHmapBP3t100$ph$tree_row,k=5)
range(hclusters)
hclusters <- data.frame(hclusters)
hclusters[,1] <- as.character(hclusters[,1])
colnames(hclusters) <- "Genes_clusters"
table(hclusters)
Hieclu1 <- subset(hclusters,Genes_clusters==1)
dim(Hieclu1)
Hieclu1genes <- rownames(Hieclu1)
Hieclu1genes
Hieclu2 <- subset(hclusters,Genes_clusters==2)
dim(Hieclu2)
Hieclu2genes <- rownames(Hieclu2)
Hieclu2genes
Hieclu5 <- subset(hclusters,Genes_clusters==5)
dim(Hieclu5)
Hieclu5genes <- rownames(Hieclu5)
Hieclu5genes
Hieclu4 <- subset(hclusters,Genes_clusters==4)
dim(Hieclu4)
Hieclu4genes <- rownames(Hieclu4)
Hieclu4genes
Hieclu3 <- subset(hclusters,Genes_clusters==3)
dim(Hieclu3)
Hieclu3genes <- rownames(Hieclu3)
Hieclu3genes
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("GOplot")
install.packages("stringr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
Hieclu1genes <- as.character(Hieclu1genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu1genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu1.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
Hieclu1.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/top100/Hieclu1.go_siggenes.rds")
Hieclu1.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/All/Hieclu1.go_siggenes.rds")
Hieclu1.go_Ssiggenes_results <- as.data.frame(Hieclu1.go_Ssiggenes@result)
write.csv(Hieclu1.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP3/top100/Hieclu1.go_siggenes_results.csv")
colnames(Hieclu1.go_Ssiggenes_results)
Hieclu1.go_Ssiggenes_results_used <- Hieclu1.go_Ssiggenes_results[,c(1,2,3,9,7)]
Hieclu1.go_Ssiggenes_results_used$geneID <- 
  str_replace(Hieclu1.go_Ssiggenes_results_used$geneID,"/",",")
names(Hieclu1.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(Hieclu1.go_Ssiggenes)
str(Hieclu1.go_Ssiggenes)
Hieclu1.go_Ssiggenes_lowP <- 
  subset(Hieclu1.go_Ssiggenes,Hieclu1.go_Ssiggenes@result$p.adjust < 0.05)
str(Hieclu1.go_Ssiggenes_lowP)
dim(Hieclu1.go_Ssiggenes_lowP)
heatplot(Hieclu1.go_Ssiggenes_lowP)
  Hieclu1.go_Ssiggenes_lowP@result$Description
  Hieclu1.go_Ssiggenes_lowP@result$pvalue
cnetplot(Hieclu1.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
BiocManager::install("clusterProfiler",update = T)
hclu3.kk <- enrichKEGG(Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu3.kk)
barplot(hclu3.kk,showCategory = 24)
Hieclu2genes <- as.character(Hieclu2genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu2genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu2.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu2.go_Ssiggenes,"~/CT/BEAM/enrich/BP3/top100/Hieclu2.go_siggenes.rds")
Hieclu2.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/top100/Hieclu2.go_siggenes.rds")
Hieclu2.go_Ssiggenes_results <- as.data.frame(Hieclu2.go_Ssiggenes@result)
write.csv(Hieclu2.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP3/top100/Hieclu2.go_siggenes_results.csv")
a <- cnetplot(Hieclu2.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
Hieclu5genes <- as.character(Hieclu5genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu5genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu5.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu5.go_Ssiggenes,"~/CT/BEAM/enrich/BP3/top100/Hieclu5.go_siggenes.rds")
Hieclu5.go_Ssiggenes_results <- as.data.frame(Hieclu5.go_Ssiggenes@result)
write.csv(Hieclu5.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP3/top100/Hieclu5.go_siggenes_results.csv")
cnetplot(Hieclu5.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
Hieclu4genes <- as.character(Hieclu4genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu4genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu4.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu4.go_Ssiggenes,"~/CT/BEAM/enrich/BP3/top100/Hieclu4.go_siggenes.rds")
Hieclu4.go_Ssiggenes_results <- as.data.frame(Hieclu4.go_Ssiggenes@result)
write.csv(Hieclu4.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP3/top100/Hieclu4.go_siggenes_results.csv")
cnetplot(Hieclu4.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
Hieclu3genes <- as.character(Hieclu3genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu3genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu3.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu3.go_Ssiggenes,"~/CT/BEAM/enrich/BP3/top100/Hieclu3.go_siggenes.rds")
Hieclu3.go_Ssiggenes_results <- as.data.frame(Hieclu3.go_Ssiggenes@result)
write.csv(Hieclu3.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP3/top100/Hieclu3.go_siggenes_results.csv")
cnetplot(Hieclu3.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
Hieclu7genes <- as.character(Hieclu7genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu7genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu7.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu7.go_Ssiggenes,"~/CT/BEAM/enrich/BP3/All/Hieclu7.go_siggenes.rds")
Hieclu7.go_Ssiggenes_results <- as.data.frame(Hieclu7.go_Ssiggenes@result)
write.csv(Hieclu7.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP3/All/Hieclu7.go_siggenes_results.csv")
cnetplot(Hieclu7.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
Hieclu1.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/All/Hieclu1.go_siggenes.rds")
cnetplot(Hieclu1.go_Ssiggenes,showCategory=20,categorySize="pvalue",circular=F,colorEdge=TRUE,layout="dh")
Hieclu3.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/All/Hieclu3.go_siggenes.rds")
cnetplot(Hieclu3.go_Ssiggenes,showCategory=20,categorySize="pvalue",circular=F,colorEdge=TRUE,layout="dh")
Hieclu1.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/top100/Hieclu1.go_siggenes.rds")
cnetplot(Hieclu1.go_Ssiggenes,showCategory=20,categorySize="pvalue",circular=F,colorEdge=TRUE,layout="dh")
Hieclu2.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/top100/Hieclu2.go_siggenes.rds")
cnetplot(Hieclu2.go_Ssiggenes,showCategory=20,categorySize="pvalue",circular=F,colorEdge=TRUE,layout="dh")
Hieclu3.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/top100/Hieclu3.go_siggenes.rds")
cnetplot(Hieclu3.go_Ssiggenes,showCategory=20,categorySize="pvalue",circular=F,colorEdge=TRUE,layout="dh")
Hieclu4.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/top100/Hieclu4.go_siggenes.rds")
cnetplot(Hieclu4.go_Ssiggenes,showCategory=20,categorySize="pvalue",circular=F,colorEdge=TRUE,layout="dh")
Hieclu5.go_Ssiggenes <-
  readRDS("~/CT/BEAM/enrich/BP3/top100/Hieclu5.go_siggenes.rds")
cnetplot(Hieclu5.go_Ssiggenes,showCategory=20,categorySize="pvalue",circular=F,colorEdge=TRUE,layout="dh")

BEAMHmapBP5t50 <- readRDS("~/BEAMHmapBP5t50.rds")
BEAMHmapBP5t100 <- readRDS("~/BEAMHmapBP5t100.rds")
hclusters <- cutree(BEAMHmapBP5t50$ph$tree_row,k=4)
range(hclusters)
hclusters <- data.frame(hclusters)
hclusters[,1] <- as.character(hclusters[,1])
colnames(hclusters) <- "Genes_clusters"
table(hclusters)
Hieclu1 <- subset(hclusters,Genes_clusters==1)
dim(Hieclu1)
Hieclu1genes <- rownames(Hieclu1)
Hieclu1genes
Hieclu2 <- subset(hclusters,Genes_clusters==2)
dim(Hieclu2)
Hieclu2genes <- rownames(Hieclu2)
Hieclu2genes
Hieclu3 <- subset(hclusters,Genes_clusters==3)
dim(Hieclu3)
Hieclu3genes <- rownames(Hieclu3)
Hieclu3genes
Hieclu4 <- subset(hclusters,Genes_clusters==4)
dim(Hieclu4)
Hieclu4genes <- rownames(Hieclu4)
Hieclu4genes
hclusters <- cutree(BEAMHmapBP5t100$ph$tree_row,k=4)
range(hclusters)
hclusters <- data.frame(hclusters)
hclusters[,1] <- as.character(hclusters[,1])
colnames(hclusters) <- "Genes_clusters"
table(hclusters)
Hieclu1 <- subset(hclusters,Genes_clusters==1)
dim(Hieclu1)
Hieclu1genes <- rownames(Hieclu1)
Hieclu1genes
Hieclu2 <- subset(hclusters,Genes_clusters==2)
dim(Hieclu2)
Hieclu2genes <- rownames(Hieclu2)
Hieclu2genes
Hieclu3 <- subset(hclusters,Genes_clusters==3)
dim(Hieclu3)
Hieclu3genes <- rownames(Hieclu3)
Hieclu3genes
Hieclu4 <- subset(hclusters,Genes_clusters==4)
dim(Hieclu4)
Hieclu4genes <- rownames(Hieclu4)
Hieclu4genes
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("GOplot")
install.packages("stringr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
Hieclu1genes <- as.character(Hieclu1genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu1genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu1.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu1.go_Ssiggenes,"~/CT/BEAM/enrich/BP5/top100/Hieclu1.go_siggenes.rds")
Hieclu1.go_Ssiggenes_results <- as.data.frame(Hieclu1.go_Ssiggenes@result)
write.csv(Hieclu1.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP5/top100/Hieclu1.go_siggenes_results.csv")
colnames(Hieclu1.go_Ssiggenes_results)
Hieclu1.go_Ssiggenes_results_used <- Hieclu1.go_Ssiggenes_results[,c(1,2,3,9,7)]
Hieclu1.go_Ssiggenes_results_used$geneID <- 
  str_replace(Hieclu1.go_Ssiggenes_results_used$geneID,"/",",")
names(Hieclu1.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(Hieclu1.go_Ssiggenes)
str(Hieclu1.go_Ssiggenes)
Hieclu1.go_Ssiggenes_lowP <- 
  subset(Hieclu1.go_Ssiggenes,Hieclu1.go_Ssiggenes@result$p.adjust < 0.05)
str(Hieclu1.go_Ssiggenes_lowP)
dim(Hieclu1.go_Ssiggenes_lowP)
heatplot(Hieclu1.go_Ssiggenes_lowP)
  Hieclu1.go_Ssiggenes_lowP@result$Description
  Hieclu1.go_Ssiggenes_lowP@result$pvalue
cnetplot(Hieclu1.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
cnetplot(Hieclu1.go_Ssiggenes,showCategory=10,
         categorySize="pvalue",
         colorEdge=TRUE,
         circular=F)
hclu3.kk <- enrichKEGG(hclu3.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu3.kk)
barplot(hclu3.kk,showCategory = 24)
Hieclu2genes <- as.character(Hieclu2genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu2genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu2.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu2.go_Ssiggenes,"~/CT/BEAM/enrich/BP5/top100/Hieclu2.go_siggenes.rds")
Hieclu2.go_Ssiggenes_results <- as.data.frame(Hieclu2.go_Ssiggenes@result)
write.csv(Hieclu2.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP5/top100/Hieclu2.go_siggenes_results.csv")
cnetplot(Hieclu2.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
Hieclu3genes <- as.character(Hieclu3genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu3genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu3.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu3.go_Ssiggenes,"~/CT/BEAM/enrich/BP5/top100/Hieclu3.go_siggenes.rds")
Hieclu3.go_Ssiggenes_results <- as.data.frame(Hieclu3.go_Ssiggenes@result)
write.csv(Hieclu3.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP5/top100/Hieclu3.go_siggenes_results.csv")
cnetplot(Hieclu3.go_Ssiggenes,showCategory=10,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
Hieclu4genes <- as.character(Hieclu4genes)
Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=Hieclu4genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
Hieclu4.go_Ssiggenes <- enrichGO(gene = Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(Hieclu4.go_Ssiggenes,"~/CT/BEAM/enrich/BP5/top100/Hieclu4.go_siggenes.rds")
Hieclu4.go_Ssiggenes_results <- as.data.frame(Hieclu4.go_Ssiggenes@result)
write.csv(Hieclu4.go_Ssiggenes_results,"~/CT/BEAM/enrich/BP5/top100/Hieclu4.go_siggenes_results.csv")
cnetplot(Hieclu4.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
mycds <- readRDS("~/noNormmycdsDDRTree3rerunOrdercell.rds")
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.rds")
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
table(mycds$State)
dim(PLTforCT)
dim(mycds)
mycds
dim(t(mycds@assayData$exprs))
dim(BigMPLTintegrated)
length(mycds$State)
BigMPLTintegrated$State <- mycds$State
BigMPLTintegrated$Pseudotime <- mycds$Pseudotime
library(presto)
statesDEgenes <- wilcoxauc(BigMPLTintegrated,"State")
head(statesDEgenes)
dplyr::count(statesDEgenes,group)
state9genes <- statesDEgenes%>%dplyr::filter(group=="9")%>%arrange(desc(auc))%>%dplyr::select(feature,auc)
state9genes
ranks <- deframe(state9genes)
head(ranks)
ranks
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
m_df <- msigdbr(species = "Homo sapiens",category = "H")
fgsea_set <- m_df %>% split(x=.$gene_symbol,f=.$gs_name)
summary(fgsea_set)
fgseaRes <- fgsea(fgsea_set,stats = ranks,nperm=100000)
fgseaRes_results <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaRes_results %>%
  dplyr::select(-leadingEdge,-ES,-nMoreExtreme) %>% 
  arrange(padj) %>% head()
ggplot(fgseaRes_results %>% filter(padj<0.05)%>%head(n=20),aes(reorder(pathway,NES),NES))
plotEnrichment(fgsea_set[["HALLMARK_NOTCH_SIGNALING"]],ranks)+
  labs(title="HALLMARK_NOTCH_SIGNALING")
fgseaRes_results
length(ranks)
DimPlot(BigMPLTintegrated,group.by = "Pseudotime",label=F)+
  NoLegend()
library(clusterProfiler)
devtools::install_github("nicolash2/gggsea")
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library('GSEABase')
library(fgsea)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
Hieclu1genes <- sort(Hieclu1genes,decreasing = T)
Hieclu1genes
hallmark <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
gsea.re1 <- GSEA(ranks,TERM2GENE=hallmark,pvalueCutoff=0.05,pAdjustMethod = "BH")
getwd()
write.table(gsea.re1,"state9wilxaucDEgenes.GSEA.txt")
g1 <- as.data.frame(gsea.re1)
g1 <- g1[order(g1$NES,decreasing = F),]
g1
num <- g1[,c(1,11)]
num <- num %>% separate_rows(core_enrichment,sep="/")%>%group_by(ID)%>%count()
num<- num[match(g1$ID,num$ID),]
sum(num$ID==g1$ID)
g1$Count <- num$n
data <- g1%>%mutate(GeneRatio=Count/setSize)
data$sign <- ifelse(data$NES>0,"activated","suppressed")
data$Description  <- gsub("HALLMARK_","",data$Description)
library(RColorBrewer)
library(wesanderson)
ggplot(data)+
  geom_point(aes(x=GeneRatio,y=Description,colour=p.adjust,size=Count))+
  facet_grid(~sign,scales="free")+
  scale_colour_gradientn(colors=wes_palette("Zissou1",80,type="continuous"))
col_gsea1 <- pal_simpsons()(16)
num1=2
gseaplot(gsea.re1,
         geneSetID = rownames(g1)[4],
          title="",
          color = col_gsea1[4],
          base_size = 14,
          rel_heights=c(1,0.2,0.4),
          subplots = 1:3,
          pvalue_table = F,
          ES_geom="line")
rownames(g1)[4]
gseaplot(gsea.re1,
         geneSetID ="HALLMARK_NOTCH_SIGNALING",
          title="",
          color = col_gsea1[4],
          base_size = 14,
          rel_heights=c(1,0.2,0.4),
          subplots = 1:3,
          pvalue_table = F,
          ES_geom="line")
cran.packages.need <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                   "gghalves", "Seurat", "SeuratObject", "methods",
                   "devtools", "BiocManager","data.table","doParallel",
                   "doRNG")
if(!requireNamespace(cran.packages.need,quietly = T)){
  install.packages(cran.packages.need,ask=F,update=F)
}
bioconductor.pac.need <- c("GSEABase", "AUCell", "SummarizedExperiment",
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree",
                           "Nebulosa")
if(!requireNamespace(bioconductor.pac.need)){
  BiocManager::install(bioconductor.pac.need,ask=F,update=F)
}
if(!requireNamespace("UCell",quietly = T)){
  devtools::install_github("carmonalab/UCell")
}
if(!requireNamespace("irGSEA",quietly = T)){
   devtools::install_github("chuiqin/irGSEA")
  }
BiocManager::install("Nebulosa")
library(irGSEA)
library(UCell)
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
mycds <- readRDS("~/noNormmycdsDDRTree3rerunOrdercell.rds")
BigMPLTintegrated$State <- mycds$State
levels(as.factor(Idents(BigMPLTintegrated)))
Idents(BigMPLTintegrated) <- BigMPLTintegrated$State
BigMPLTintegrated.Enscore <- irGSEA.score(
  BigMPLTintegrated,assay="RNA",slot="data",seeds=1,ncores=1,min.cells = 3,
  min.feature = 0,custom = F, geneset = NULL,msigdb = T, species = "Homo sapiens",
  category = "H", subcategory = NULL, geneid = "symbol",
  method = c("AUCell","UCell","singscore","ssgsea"),kcdf = "Gaussian")
msigdbr::msigdbr_collections
Seurat::Assays(BigMPLTintegrated.Enscore)
str(BigMPLTintegrated.Enscore)
saveRDS(BigMPLTintegrated.Enscore,"~/BigMPLTintegrated.Enscore.byCTstate.rds")
result.dge <- irGSEA.integrate(BigMPLTintegrated.Enscore,
                               group.by = "State",
                               method =c("AUCell","UCell","singscore","ssgsea"))
result.dge
class(result.dge)
saveRDS(result.dge,"~/BigMPLTintegrated-Enrich-result.dge.byCTstate.rds")
write.csv(result.dge$RRA,"~/BigMPLTintegrated-Enrich-RRA.byCTstate.csv")
write.csv(result.dge$UCell,"~/BigMPLTintegrated-Enrich-UCell.byCTstate.csv")
write.csv(result.dge$AUCell,"~/BigMPLTintegrated-Enrich-AUCell.byCTstate.csv")
write.csv(result.dge$singscore,"~/BigMPLTintegrated-Enrich-singscore.byCTstate.csv")
write.csv(result.dge$ssgsea,"~/BigMPLTintegrated-Enrich-ssgsea.byCTstate.csv")
result.dge.RRA.sorted <- read.csv("~/BigMPLTintegrated-Enrich-RRA.byCTstate.csv")
result.dge.UCell.sorted <- read.csv("~/BigMPLTintegrated-Enrich-UCell.byCTstate.csv")
result.dge.AUCell.sorted <- read.csv("~/BigMPLTintegrated-Enrich-AUCell.byCTstate.csv")
result.dge.singscore.sorted <- read.csv("~/BigMPLTintegrated-Enrich-singscore.byCTstate.csv")
result.dge.ssgsea.sorted <- read.csv("~/BigMPLTintegrated-Enrich-ssgsea.byCTstate.csv")
result.dge.sorted <- result.dge
result.dge.sorted$RRA <- result.dge.RRA.sorted
result.dge.sorted$UCell <- result.dge.UCell.sorted
result.dge.sorted$AUCell <- result.dge.AUCell.sorted
result.dge.sorted$singscore <- result.dge.singscore.sorted
result.dge.sorted$ssgsea <- result.dge.ssgsea.sorted
# cols = c(pal_npg()(8),"
#          "
#          "
hmap1 <- irGSEA.heatmap(result.dge.sorted,method="RRA",
                       heatmap.width = 15,rowname.fointsize = 10,
                       direction.color=c("darkred","darkblue"))
                       # significance.color = c("white","
hmap1
result.dge.sorted$RRA$pvalue
hmap2 <- irGSEA.heatmap(result.dge.sorted,method="UCell",
                       heatmap.width = 15,rowname.fointsize = 10,
                       direction.color=c("darkred","darkblue"))
                       # significance.color = c("white","
hmap2
hmap3 <- irGSEA.heatmap(result.dge.sorted,method="AUCell",
                       heatmap.width = 15,rowname.fointsize = 10,
                       direction.color=c("darkred","darkblue"))
                       # significance.color = c("white","
hmap3
hmap4 <- irGSEA.heatmap(result.dge.sorted,method="singscore",
                       heatmap.width = 15,rowname.fointsize = 10,
                       direction.color=c("darkred","darkblue"))
                       # significance.color = c("white","
hmap4
hmap5 <- irGSEA.heatmap(result.dge.sorted,method="ssgsea",
                       heatmap.width = 15,rowname.fointsize = 10,
                       direction.color=c("darkred","darkblue"))
                       # significance.color = c("white","
hmap5
bubbleplot1 <- irGSEA.bubble(result.dge.sorted,,method="RRA",
                       direction.color=c("blue","red"))
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot1
bubbleplot2 <- irGSEA.bubble(result.dge.sorted,,method="UCell",
                        direction.color=c("blue","red"))
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot2
bubbleplot3 <- irGSEA.bubble(result.dge.sorted,,method="AUCell",
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                        direction.color=c("blue","red"))+
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot3
bubbleplot4 <- irGSEA.bubble(result.dge.sorted,,method="singscore",
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                        direction.color=c("blue","red"))+
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot4
bubbleplot5 <- irGSEA.bubble(result.dge.sorted,,method="ssgsea",
                       # cluster.color =c("
                       #                  "
                       #                  "
                       #                  "
                       #                  "
                       direction.color=c("blue","red"))+
                       # significance.color = c("white","
  geom_point(size=10)
bubbleplot5
scatterplot <- irGSEA.density.scatterplot(BigMPLTintegrated.Enscore,
                                          method="AUCell",
                                          show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE",
                                          reduction="tsne")
scatterplot

mycds <- readRDS("~/noNormmycdsDDRTree3.rds")
mycds <- readRDS("~/noNormmycdsDDRTree3rerunOrdercell.rds")
mycds$Pseudotime
mycds$State
levels(mycds$State)
mycds$Cluster
levels(mycds$Cluster)
mycds@assayData$exprs
mycds.S <- subset(mycds,mycds$MSH.idents == "S")
plot_cell_trajectory(mycds.S, color_by = 'as.factor(State)')
mycds_expressed_genes_S <-  
  row.names(subset(fData(mycds.S),
                   num_cells_expressed >= 5))
length(mycds_expressed_genes_S)
mycds_filteredLow_S <- mycds.S[mycds_expressed_genes_S,]
diff_test_res_S <- 
  differentialGeneTest(mycds_filteredLow_S,fullModelFormulaStr = "~sm.ns(Pseudotime)")
saveRDS(diff_test_res_S,"~/CT/DE/SvsM/diff_test_res_S.rds")
library(dplyr)
diff_test_res_S[,c("gene_short_name", "pval", "qval")]%>%head()
sig.gene.S <-row.names(subset(diff_test_res_S, qval < 0.05))
length(sig.gene.S)
plot_genes_in_pseudotime(mycds_filteredLow_S[sig.gene.S[1:10]],
                         color_by = 'Cluster',ncol = 5)
g <- plot_pseudotime_heatmap(mycds_filteredLow_S[sig.gene.S,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
g
diff_test_res_S <- readRDS("~/CT/DE/SvsM/diff_test_res_S.rds")
sig.gene.S <-row.names(subset(diff_test_res_S, qval < 0.05))
class(mycds_filteredLow_S[sig.gene.S,]@assayData$exprs)
S.siggene.matrix <- mycds_filteredLow_S[sig.gene.S,]@assayData$exprs
distS <- dist(S.siggene.matrix, method="euclidean")
str(distS)
saveRDS(distS,"~/CT/DE/SvsM/distS.rds")
clust <- hclust(distS,method="ward.D2")
clust$order
plot(clust, labels=S.siggene.matrix[,2],cex=0.3)
a <- plot_pseudotime_heatmap(mycds_filteredLow_S[sig.gene.S,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap=TRUE) 
k <- (cutree(a$tree_row,k=3)-2) %% 3 +1
levels(k)
length(k)
str(k)
hclu1 <- subset(k,k==1)
length(hclu1)
hclu1.genes <- names(hclu1)
hclu1.genes
hclu2 <- subset(k,k==2)
length(hclu2)
hclu2.genes <- names(hclu2)
hclu2.genes
hclu3 <- subset(k,k==3)
length(hclu3)
str(hclu3)
hclu3.genes <- names(hclu3)
hclu3.genes
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("GOplot")
install.packages("stringr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
hclu3.genes <- as.character(hclu3.genes)
hclu3.genes
hclu3.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu3.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu3.go_Ssiggenes <- enrichGO(gene = hclu3.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
hclu3.go_Ssiggenes_results <- as.data.frame(hclu3.go_Ssiggenes@result)
colnames(hclu3.go_Ssiggenes_results)
hclu3.go_Ssiggenes_results_used <- hclu3.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu3.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu3.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu3.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu3.go_Ssiggenes)
str(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes_lowP <- 
  subset(hclu3.go_Ssiggenes,hclu3.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu3.go_Ssiggenes_lowP)
dim(hclu3.go_Ssiggenes_lowP)
heatplot(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes@result$Description
hclu3.go_Ssiggenes@result$pvalue
cnetplot(hclu3.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu3.kk <- enrichKEGG(hclu3.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu3.kk)
barplot(hclu3.kk,showCategory = 24)
hclu2.genes <- as.character(hclu2.genes)
hclu2.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu2.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu2.kk <- enrichKEGG(hclu2.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu2.kk)
barplot(hclu2.kk,showCategory = 18)
hclu1.genes <- as.character(hclu1.genes)
hclu1.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu1.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu1.kk <- enrichKEGG(hclu1.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu1.kk)
barplot(hclu1.kk,showCategory = 4)
mycds <- readRDS("~/noNormmycdsDDRTree3.rds")
mycds$Pseudotime
mycds$State
levels(mycds$State)
mycds$Cluster
levels(mycds$Cluster)
mycds@assayData$exprs
mycds.M <- subset(mycds,mycds$MSH.idents == "M")
plot_cell_trajectory(mycds.M, color_by = 'as.factor(State)')
mycds_expressed_genes_M <-  row.names(subset(fData(mycds.M),
                                             num_cells_expressed >= 5))
length(mycds_expressed_genes_M)
mycds_filteredLow_M <- mycds.M[mycds_expressed_genes_M,]
diff_test_res_M <- differentialGeneTest(mycds_filteredLow_M,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
str(diff_test_res_M)
saveRDS(diff_test_res_M,"~/CT/DE/SvsM/diff_test_res_M.rds")
library(dplyr)
diff_test_res_M[,c("gene_short_name", "pval", "qval")]%>%head()
sig.gene.M <-row.names(subset(diff_test_res_M, qval < 0.05))
length(sig.gene.M)
plot_genes_in_pseudotime(mycds_filteredLow_M[sig.gene.M[1:10]],
                         color_by = 'Cluster',ncol = 5)
g <- plot_pseudotime_heatmap(mycds_filteredLow_M[sig.gene.M,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
g
class(mycds_filteredLow_S[sig.gene.S,]@assayData$exprs)
M.siggene.matrix <- mycds_filteredLow_M[sig.gene.M,]@assayData$exprs
dist <- dist(M.siggene.matrix, method="euclidean")
saveRDS(dist,"~/CT/DE/SvsM/distM.rds")
clust <- hclust(dist,method="ward.D2")
clust$order
plot(clust, labels=M.siggene.matrix[,2],cex=0.3)
a <- plot_pseudotime_heatmap(mycds_filteredLow_M[sig.gene.M,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap=TRUE)
k <- (cutree(a$tree_row,k=3)-2) %% 3 +1
levels(k)
length(k)
str(k)
hclu1 <- subset(k,k==1)
length(hclu1)
hclu1.genes <- names(hclu1)
hclu1.genes
hclu2 <- subset(k,k==2)
length(hclu2)
hclu2.genes <- names(hclu2)
hclu2.genes
hclu3 <- subset(k,k==3)
length(hclu3)
str(hclu3)
hclu3.genes <- names(hclu3)
hclu3.genes
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("GOplot")
install.packages("stringr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
hclu3.genes <- as.character(hclu3.genes)
hclu3.genes
hclu3.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu3.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu3.go_Ssiggenes <- enrichGO(gene = hclu3.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
hclu3.go_Ssiggenes_results <- as.data.frame(hclu3.go_Ssiggenes@result)
colnames(hclu3.go_Ssiggenes_results)
hclu3.go_Ssiggenes_results_used <- hclu3.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu3.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu3.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu3.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu3.go_Ssiggenes)
str(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes_lowP <- 
  subset(hclu3.go_Ssiggenes,hclu3.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu3.go_Ssiggenes_lowP)
dim(hclu3.go_Ssiggenes_lowP)
heatplot(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes@result$Description
hclu3.go_Ssiggenes@result$pvalue
cnetplot(hclu3.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu3.kk <- enrichKEGG(hclu3.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu3.kk)
barplot(hclu3.kk,showCategory = 24)
hclu2.genes <- as.character(hclu2.genes)
hclu2.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu2.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu2.kk <- enrichKEGG(hclu2.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu2.kk)
barplot(hclu2.kk,showCategory = 18)
hclu1.genes <- as.character(hclu1.genes)
hclu1.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu1.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu1.kk <- enrichKEGG(hclu1.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu1.kk)
barplot(hclu1.kk,showCategory = 4)

mycds <- readRDS("~/noNormmycdsDDRTree3.rds")
mycds@assayData$exprs
install.packages("remotes")
devtools::install_github("jacobheng/cellwrangler")
table(is.na(mycds$Covidsplit))
mycds$Covidsplit[is.na(mycds$Covidsplit)] <- "notuse"
table(is.na(mycds$Covidsplit))
library("cellwrangler")
mycds.Cov <- subset_cds("use",mycds,"Covidsplit")
str(mycds.Cov)
dim(mycds.Cov)
dim(mycds)
table(mycds.Cov$Covidsplit)
mycds.Cov <- detectGenes(mycds.Cov, min_expr = 0.1)
mycds_expressed_genes_Cov <-  row.names(subset(fData(mycds.Cov),
                                               num_cells_expressed >= 5))
length(mycds_expressed_genes_Cov)
mycds_filteredLow_Cov <- mycds.Cov[mycds_expressed_genes_Cov,] 
mycds_expressed_genes 
disp_table <- dispersionTable(mycds_filteredLow_Cov)
disp.genes <- subset(disp_table, 
                     mean_expression >= 0.1 
                     &dispersion_empirical >= 1 * dispersion_fit)
length(disp.genes$gene_id)
mycds_filteredLowDisp_Cov <- mycds_filteredLow_Cov[disp.genes$gene_id,] 
diff_test_res_Cov <- differentialGeneTest(mycds_filteredLowDisp_Cov,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(diff_test_res_Cov,"~/diff_test_res_Cov.csv")
library(dplyr)
diff_test_res_Cov[,c("gene_short_name", "pval", "qval")]%>%head()
sig.gene.Cov <-row.names(subset(diff_test_res_Cov, qval < 0.05))
length(sig.gene.Cov)
write.table(sig.gene.Cov,"~/sig.gene.Cov.txt")
plot_genes_in_pseudotime(mycds_filteredLowDisp_Cov[sig.gene.Cov[1:10]],
                         color_by = 'State',ncol = 5)
plot_pseudotime_heatmap(mycds_filteredLowDisp_Cov[sig.gene.Cov,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap=TRUE)
saveRDS(mycds_filteredLowDisp_Cov,"~/mycds_filteredLowDisp_Cov.rds")
class(mycds_filteredLowDisp_Cov[sig.gene.Cov,]@assayData$exprs)
Cov.siggene.matrix <- mycds_filteredLowDisp_Cov[sig.gene.Cov,]@assayData$exprs
dist <- dist(Cov.siggene.matrix, method="euclidean")
clust <- hclust(dist,method="ward.D2")
clust$order
plot(clust, labels=Cov.siggene.matrix[,2],cex=0.3)
a <- plot_pseudotime_heatmap(mycds_filteredLowDisp_Cov[sig.gene.Cov,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap=TRUE) 
k <- (cutree(a$tree_row,k=3)-2) %% 3 +1
cutree(a$tree_row,k=3)  
levels(k)
length(k)
str(k)
hclu1 <- subset(k,k==1)
length(hclu1)
hclu1.genes <- names(hclu1)
hclu1.genes
hclu2 <- subset(k,k==2)
length(hclu2)
hclu2.genes <- names(hclu2)
hclu2.genes
hclu3 <- subset(k,k==3)
length(hclu3)
str(hclu3)
hclu3.genes <- names(hclu3)
hclu3.genes
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("GOplot")
install.packages("stringr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
hclu3.genes <- as.character(hclu3.genes)
hclu3.genes
hclu3.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu3.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu3.go_Ssiggenes <- enrichGO(gene = hclu3.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(hclu3.go_Ssiggenes,"~/hclu3.go_Covsiggenes.rds")
hclu3.go_Ssiggenes_results <- as.data.frame(hclu3.go_Ssiggenes@result)
write.csv(hclu3.go_Ssiggenes_results,"~/hclu3.go_Covsiggenes_results.csv")
colnames(hclu3.go_Ssiggenes_results)
hclu3.go_Ssiggenes_results_used <- hclu3.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu3.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu3.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu3.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu3.go_Ssiggenes)
str(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes_lowP <- 
  subset(hclu3.go_Ssiggenes,hclu3.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu3.go_Ssiggenes_lowP)
dim(hclu3.go_Ssiggenes_lowP)
heatplot(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes@result$Description
hclu3.go_Ssiggenes@result$pvalue
cnetplot(hclu3.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu3.kk <- enrichKEGG(hclu3.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu3.kk)
barplot(hclu3.kk,showCategory = 24)
hclu2.genes <- as.character(hclu2.genes)
hclu2.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu2.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu2.go_Ssiggenes <- enrichGO(gene = hclu2.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(hclu2.go_Ssiggenes,"~/hclu2.go_Covsiggenes.rds")
hclu2.go_Ssiggenes_results <- as.data.frame(hclu2.go_Ssiggenes@result)
write.csv(hclu2.go_Ssiggenes_results,"~/Covclu2.go_Covsiggenes_results.csv")
colnames(hclu2.go_Ssiggenes_results)
hclu2.go_Ssiggenes_results_used <- hclu2.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu2.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu2.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu2.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu2.go_Ssiggenes)
str(hclu2.go_Ssiggenes)
hclu2.go_Ssiggenes_lowP <- 
  subset(hclu2.go_Ssiggenes,hclu2.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu2.go_Ssiggenes_lowP)
dim(hclu2.go_Ssiggenes_lowP)
heatplot(hclu2.go_Ssiggenes)
hclu2.go_Ssiggenes@result$Description
hclu2.go_Ssiggenes@result$pvalue
cnetplot(hclu2.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu2.kk <- enrichKEGG(hclu2.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu2.kk)
barplot(hclu2.kk,showCategory = 18)
hclu1.genes <- as.character(hclu1.genes)
hclu1.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu1.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu1.go_Ssiggenes <- enrichGO(gene = hclu1.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(hclu1.go_Ssiggenes,"~/hclu1.go_Covsiggenes.rds")
hclu1.go_Ssiggenes_results <- as.data.frame(hclu1.go_Ssiggenes@result)
write.csv(hclu1.go_Ssiggenes_results,"~/Covclu1.go_Covsiggenes_results.csv")
colnames(hclu1.go_Ssiggenes_results)
hclu1.go_Ssiggenes_results_used <- hclu1.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu1.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu1.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu1.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu1.go_Ssiggenes)
str(hclu1.go_Ssiggenes)
hclu1.go_Ssiggenes_lowP <- 
  subset(hclu1.go_Ssiggenes,hclu1.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu1.go_Ssiggenes_lowP)
dim(hclu1.go_Ssiggenes_lowP)
heatplot(hclu1.go_Ssiggenes)
hclu1.go_Ssiggenes@result$Description
hclu1.go_Ssiggenes@result$pvalue
cnetplot(hclu1.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu1.kk <- enrichKEGG(hclu1.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu1.kk)
barplot(hclu1.kk,showCategory = 4)
mycds <- readRDS("~/noNormmycdsDDRTree3.rds")
mycds@assayData$exprs
install.packages("remotes")
devtools::install_github("jacobheng/cellwrangler")
table(is.na(mycds$Covidsplit))
mycds$Covidsplit[is.na(mycds$Covidsplit)] <- "use"
table(is.na(mycds$Covidsplit))
library("cellwrangler")
mycds.nonCov <- subset_cds("notuse",mycds,"Covidsplit")
str(mycds.nonCov)
dim(mycds.nonCov)
dim(mycds)
table(mycds.nonCov$Covidsplit)
mycds.nonCov <- detectGenes(mycds.nonCov, min_expr = 0.1)
mycds_expressed_genes_nonCov <-  row.names(subset(fData(mycds.nonCov),
                                               num_cells_expressed >= 5))
length(mycds_expressed_genes_nonCov)
mycds_filteredLow_nonCov <- mycds.nonCov[mycds_expressed_genes_nonCov,] 
disp_table <- dispersionTable(mycds_filteredLow_nonCov)
disp.genes <- subset(disp_table, 
                     mean_expression >= 0.1 
                     &dispersion_empirical >= 1 * dispersion_fit)
length(disp.genes$gene_id)
mycds_filteredLowDisp_nonCov <- mycds_filteredLow_nonCov[disp.genes$gene_id,] 
diff_test_res_nonCov <- differentialGeneTest(mycds_filteredLowDisp_nonCov,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
str(diff_test_res_nonCov)
write.csv(diff_test_res_nonCov,"~/diff_test_res_nonCov.csv")
library(dplyr)
diff_test_res_Cov[,c("gene_short_name", "pval", "qval")]%>%head()
sig.gene.nonCov <-row.names(subset(diff_test_res_nonCov, qval < 0.05))
length(sig.gene.nonCov)
write.table(sig.gene.nonCov,"~/sig.gene.nonCov.txt")
plot_genes_in_pseudotime(mycds_filteredLowDisp_nonCov[sig.gene.nonCov[1:10]],
                         color_by = 'State',ncol = 5)
plot_pseudotime_heatmap(mycds_filteredLowDisp_nonCov[sig.gene.nonCov,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap=TRUE)
saveRDS(mycds_filteredLowDisp_nonCov,"~/mycds_filteredLowDisp_nonCov.rds")
class(mycds_filteredLowDisp_nonCov[sig.gene.nonCov,]@assayData$exprs)
nonCov.siggene.matrix <-
  mycds_filteredLowDisp_nonCov[sig.gene.nonCov,]@assayData$exprs
dist <- dist(nonCov.siggene.matrix, method="euclidean")
clust <- hclust(dist,method="ward.D2")
clust$order
plot(clust, labels=nonCov.siggene.matrix[,2],cex=0.3)
a <- plot_pseudotime_heatmap(mycds_filteredLowDisp_nonCov[sig.gene.nonCov,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap=TRUE) 
k <- (cutree(a$tree_row,k=3)-2) %% 3 +1
levels(k)
length(k)
str(k)
hclu1 <- subset(k,k==1)
length(hclu1)
hclu1.genes <- names(hclu1)
hclu1.genes
hclu2 <- subset(k,k==2)
length(hclu2)
hclu2.genes <- names(hclu2)
hclu2.genes
hclu3 <- subset(k,k==3)
length(hclu3)
str(hclu3)
hclu3.genes <- names(hclu3)
hclu3.genes
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("GOplot")
install.packages("stringr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
hclu3.genes <- as.character(hclu3.genes)
hclu3.genes
hclu3.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu3.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu3.go_Ssiggenes <- enrichGO(gene = hclu3.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(hclu3.go_Ssiggenes,"~/hclu3.go_nonCsiggenes.rds")
hclu3.go_Ssiggenes_results <- as.data.frame(hclu3.go_Ssiggenes@result)
write.csv(hclu3.go_Ssiggenes_results,"~/hclu3.go_nonCsiggenes_results.csv")
colnames(hclu3.go_Ssiggenes_results)
hclu3.go_Ssiggenes_results_used <- hclu3.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu3.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu3.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu3.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu3.go_Ssiggenes)
str(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes_lowP <- 
  subset(hclu3.go_Ssiggenes,hclu3.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu3.go_Ssiggenes_lowP)
dim(hclu3.go_Ssiggenes_lowP)
heatplot(hclu3.go_Ssiggenes)
hclu3.go_Ssiggenes@result$Description
hclu3.go_Ssiggenes@result$pvalue
cnetplot(hclu3.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu3.kk <- enrichKEGG(hclu3.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu3.kk)
barplot(hclu3.kk,showCategory = 24)
hclu2.genes <- as.character(hclu2.genes)
hclu2.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu2.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu2.go_Ssiggenes <- enrichGO(gene = hclu2.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
saveRDS(hclu2.go_Ssiggenes,"~/hclu2.go_nonCsiggenes.rds")
hclu2.go_Ssiggenes_results <- as.data.frame(hclu2.go_Ssiggenes@result)
write.csv(hclu2.go_Ssiggenes_results,"~/hclu2.go_nonCsiggenes_results.csv")
colnames(hclu2.go_Ssiggenes_results)
hclu2.go_Ssiggenes_results_used <- hclu2.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu2.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu2.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu2.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu2.go_Ssiggenes)
str(hclu2.go_Ssiggenes)
hclu2.go_Ssiggenes_lowP <- 
  subset(hclu2.go_Ssiggenes,hclu2.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu2.go_Ssiggenes_lowP)
dim(hclu2.go_Ssiggenes_lowP)
heatplot(hclu2.go_Ssiggenes)
hclu2.go_Ssiggenes@result$Description
hclu2.go_Ssiggenes@result$pvalue
cnetplot(hclu2.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu2.kk <- enrichKEGG(hclu2.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu2.kk)
barplot(hclu2.kk,showCategory = 18)
hclu1.genes <- as.character(hclu1.genes)
hclu1.Ssig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=hclu1.genes,
                        keytype = "SYMBOL",
                        column="ENTREZID")
hclu1.go_Ssiggenes <- enrichGO(gene = hclu1.Ssig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
str(hclu1.go_Ssiggenes)
saveRDS(hclu1.go_Ssiggenes,"~/hclu1.go_nonCsiggenes.rds")
hclu1.go_Ssiggenes_results <- as.data.frame(hclu1.go_Ssiggenes@result)
write.csv(hclu1.go_Ssiggenes_results,"~/hclu1.go_nonCsiggenes_results.csv")
colnames(hclu1.go_Ssiggenes_results)
hclu1.go_Ssiggenes_results_used <- hclu1.go_Ssiggenes_results[,c(1,2,3,9,7)]
hclu1.go_Ssiggenes_results_used$geneID <- 
  str_replace(hclu1.go_Ssiggenes_results_used$geneID,"/",",")
names(hclu1.go_Ssiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(hclu1.go_Ssiggenes)
str(hclu1.go_Ssiggenes)
hclu1.go_Ssiggenes_lowP <- 
  subset(hclu1.go_Ssiggenes,hclu1.go_Ssiggenes@result$p.adjust < 0.05)
str(hclu1.go_Ssiggenes_lowP)
dim(hclu1.go_Ssiggenes_lowP)
heatplot(hclu1.go_Ssiggenes)
hclu1.go_Ssiggenes@result$Description
hclu1.go_Ssiggenes@result$pvalue
cnetplot(hclu1.go_Ssiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
hclu1.kk <- enrichKEGG(hclu1.Ssig_entrezId, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(hclu1.kk)
barplot(hclu1.kk,showCategory = 4)
mycds.cluI <- readRDS('Mt-Cleared-Outputs/Secondary_preprocess_codes_output/CT/Unsup_Secondary/method1/mycds_unsup.rds')
a <- plot_cell_clusters(mycds.cluI, 1, 2, color = "Cluster",
                                show_cell_names =T,
                                cell_name_size =3,
                                cell_size = 1.5)
a$layers[[1]]$aes_params$alpha <- 0.5
a
b <- plot_cell_clusters(mycds.cluI, 1, 2, color = "cond.idents",
                                show_cell_names =T,
                                cell_name_size =3,
                                cell_size = 1.5)
b$layers[[1]]$aes_params$alpha <- 0.5
b
mycds.cluI$RNA_snn_res.0.5
mycds.cluI$seurat_clusters
mycds.cluI$SRlabelsmain
mycds.cluI$Cluster
mycds_expressed_genes.cluI <-  row.names(subset(fData(mycds.cluI),
                                           num_cells_expressed >= 10))
mycds_filteredLow <- mycds.unsup[mycds_expressed_genes.cluI,]
set.seed(1)
clustering_DEG_genes_unsup.cluI <-
  differentialGeneTest(mycds.unsup[mycds_expressed_genes.cluI,],
                       fullModelFormulaStr = '~Cluster')
saveRDS(clustering_DEG_genes_unsup.cluI,"Mt-Cleared-Outputs/Secondary_preprocess_codes_output/CT/Unsup_Secondary/cclustering_DEG_genes_unsup.cluI.rds")
clustering_DEG_genes_unsup.cluI
sig.genes.cluI <- subset(clustering_DEG_genes_unsup.cluI,qval<0.1)
sig.genes.cluI.sortIndex <- sort(sig.genes.cluI$qval,index.return=T)
sig.genes.cluI.sortIndex$ix
sig.genes.cluI.sorted <- sig.genes.cluI[sig.genes.cluI.sortIndex$ix,]
sig.genes.cluI.sorted
sg.plot <- as.character(sig.genes.cluI.sorted$gene_short_name)
length(sg.plot)
blast.genes <- sg.plot[1:28]
cluIplot <- plot_genes_jitter(mycds.unsup[blast.genes,],
                  grouping = "Cluster",color_by="Cluster",nrow=4,ncol=7,cell_size = 0.08,relative_expr = TRUE)
cluIplot$layers[[1]]$aes_params$alpha <- 0.3
cluIplot
diff_test_res <- differentialGeneTest(mycds[mycds_expressed_genes,],
                                      fullModelFormulaStr = "~Cluster")
diff_test_res[1:4,1:4]
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
sessionInfo()
BiocManager::install(version = "3.14")
BiocManager::install(version = "3.18")
if(!require(monocle3))BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',                                             'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',                                             'SummarizedExperiment', 'batchelor', 'Matrix.utils',                                             'HDF5Array', 'terra', 'ggrastr'))
if(!require(devtools))install.packages("devtools")
if(!require(monocle3))devtools::install_github('cole-trapnell-lab/monocle3')
gene_M_df <- find_gene_modules(mycds.Cov[sig.gene.Cov,],resolution=c(10^seq(-6,-1)))
devtools::install_github("cole-trapnell-lab/cicero-release",ref="monocle3")
mycds <- readRDS("~/noNormmycdsDDRTree3.rds")
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$CTPseudotime <- mycds$Pseudotime
str(mycds$Pseudotime)
newPseu <- na.omit(mycds$Pseudotime)
length(newPseu)
table(is.na(mycds$Pseudotime))
table(is.na(mycds$PLT.idents))
newmycds = mycds[,mycds$PLT.idents %in% c("cluA","cluB","cluC","cluD","cluE","cluF",
                                          "cluG","cluH","cluI")]
table(mycds$PLT.idents)
table(is.na(newmycds$PLT.idents))
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$CTPseudotime <- newmycds$Pseudotime
BigMPLTintegrateduse.forWGCNAv2refpca@meta.data$CTState <- newmycds$State
saveRDS(BigMPLTintegrateduse.forWGCNAv2refpca,"~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
CXCL_Ssig <- c("")
CXCL_Allsig <- c("CXCR3","CXCR2","PF4","PF4V1","CXCL5","PPBP")
CXCL_Allsig_entrezId <- mapIds(x=org.Hs.eg.db,
                        keys=CXCL_Allsig,
                        keytype = "SYMBOL",
                        column="ENTREZID")
CXCL_Allsiggenes <- enrichGO(gene = CXCL_Allsig_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont="ALL",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
CXCL_Allsiggenes_results <- as.data.frame(CXCL_Allsiggenes@result)
colnames(CXCL_Allsiggenes_results)
CXCL_Allsiggenes_results_used <- CXCL_Allsiggenes_results[,c(1,2,3,9,7)]
CXCL_Allsiggenes_results_used$geneID <- 
  str_replace(CXCL_Allsiggenes_results_used$geneID,"/",",")
names(CXCL_Allsiggenes_results_used) <- c("Category","ID","term","Genes","adj_pval")
dim(CXCL_Allsiggenes)
str(CXCL_Allsiggenes)
str(CXCL_Allsiggenes)
str(CXCL_Allsiggenes[temp2[1:50],])
plotdata <- CXCL_Allsiggenes[temp2,]
str(plotdata)
CXCL_Allsiggenes@result$Description
CXCL_Allsiggenes@result$pvalue
cnetplot(CXCL_Allsiggenes,categorySize="pvalue",circular=TRUE,colorEdge=TRUE)
CXCL_Allsig.kk <- enrichKEGG(CXCL_Allsig_entrezId, organism = "hsa", 
                             keyType = "kegg", pvalueCutoff = 0.05,
                             pAdjustMethod = "BH",qvalueCutoff = 0.05)
dim(CXCL_Allsig.kk)
barplot(CXCL_Allsig.kk)
BigMergereU <- readRDS("~/BigMergereU.rds")
BigMPLTintegrateduse.forWGCNAv2refpca <- readRDS("~/BigMPLTintegrateduse.forWGCNAv2refpca.anotherday2.rds")
CovintegratedreU0 <- CovintegratedAnnorerun2SR
CovintegratedreU0@meta.data$finalgeneralumap
PLTCovMSH.cluster.final@meta.data$RNA_snn_res.1.5
CovintegratedreU0@meta.data$finalgeneralumap[1]
length(CovintegratedreU0@meta.data$finalgeneralumap)
CovintegratedreU0@meta.data$finalgeneralumapreU0 <- CovintegratedreU0@meta.data$finalgeneralumap
n = 1
for (i in (1:48333)){
  if (CovintegratedreU0@meta.data$finalgeneralumapreU0[i] == "Platelet"){
    CovintegratedreU0@meta.data$finalgeneralumapreU0[i] <-
      PLTCovMSH.cluster.final@meta.data$RNA_snn_res.1.5[n]
    n = n+1
  }
  i = i+1
}
levels(as.factor(CovintegratedreU0@meta.data$finalgeneralumapreU0))
levels(as.factor(PLTCovMSH.cluster.final@meta.data$RNA_snn_res.1.5))
saveRDS(CovintegratedreU0,'4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
BigMergereU <- readRDS("~/BigMergereU20.rds")
table(BigMergereU$CCCreU)
table(BigMergereU$predicted.celltypeBasicwithnewreU)
table(BigMergereU$predicted.celltypewithnew)
Idents(BigMergereU) <- BigMergereU$predicted.celltypeBasicwithnewreU
suppressMessages(if(!require(CellChat))devtools::install_github("sqjin/CellChat"))
suppressMessages(if(!require(ggplot2))install.packages('ggplot2'))
suppressMessages(if(!require(patchwork))install.packages('patchwork') )
suppressMessages(if(!require(ggalluvial))install.packages('ggalluvial'))
suppressMessages(if(!require(igraph))install.packages('igraph'))
suppressMessages(if(!require(dplyr))install.packages('dplyr'))
suppressMessages(options(stringsAsFactors = FALSE))
suppressMessages(options(futrue.globlas.Maxsize=2*1024**3))
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 8)))
library(CellChat)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
CellChatDB <- CellChatDB.human

BigMergereU <- readRDS("~/BigMergereU20.rds")
table(BigMergereU$predicted.celltypeBasicwithnewreU)
table(BigMergereU$CCCreU)
BigMergereU <- readRDS("~/BigMergereU.rds")
Idents(BigMergereU) <- BigMergereU$predicted.celltypeBasicwithnewreU
data.input <- BigMergereU[["RNA"]]@data
meta <- BigMergereU@meta.data
unique(meta$Covidsplit)
Covcell.use <- rownames(meta)[meta$Covidsplit == "use"]
data.input <- data.input[, Covcell.use]
meta = meta[Covcell.use, ]
unique(meta$predicted.celltypeBasicwithnewreU)
cellchatCov <- createCellChat(object = data.input, meta = meta, group.by = "predicted.celltypeBasicwithnewreU")
cellchatCov <- setIdent(cellchatCov, 
                        ident.use = "predicted.celltypeBasicwithnewreU")
levels(cellchatCov@idents)
groupSizeCov <- as.numeric(table(cellchatCov@idents))
groupSizeCov
saveRDS(cellchatCov,"~/cellchatCovpre.rds")

Idents(BigMergereU) <- BigMergereU$predicted.celltypeBasicwithnewreU
data.input <- BigMergereU[["RNA"]]@data
meta <- BigMergereU@meta.data
unique(meta$Covidsplit)
nonCovcell.use <- rownames(meta)[meta$Covidsplit == "notuse"]
data.input <- data.input[, nonCovcell.use]
meta = meta[nonCovcell.use, ]
unique(meta$predicted.celltypeBasicwithnewreU)
cellchatnonCov <- createCellChat(object = data.input, meta = meta, group.by = "predicted.celltypeBasicwithnewreU")
cellchatnonCov <- setIdent(cellchatnonCov, 
                        ident.use = "predicted.celltypeBasicwithnewreU")
levels(cellchatnonCov@idents)
groupSizenonCov <- as.numeric(table(cellchatnonCov@idents))
groupSizenonCov
saveRDS(cellchatnonCov,"~/cellchatnonCovpre.rds")
cellchatCov <- readRDS("~/cellchatCovprewithcluX.rds")
cellchatnonCov <- readRDS("~/cellchatnonCovprewithcluX.rds")
cellchat <- cellchatCov
cellchat <- cellchatnonCov
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
options(future.globals.maxSize=4000000000)
cellchat <- computeCommunProb(cellchat, raw.use = T)
 
cellchat <- filterCommunication(cellchat, min.cells = 5)
df.net <- subsetCommunication(cellchat)
class(df.net)
DT::datatable(df.net)
write.csv(df.net,'~/Cov.df.netwithcluX.csv')
write.csv(df.net,'~/nonCov.df.netwithcluX.csv')
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,"~/cellchatCovwithcluX.rds")
saveRDS(cellchat,"~/cellchatnonCovwithcluX.rds")
cellchatCov <- readRDS("~/cellchatCov.rds")
cellchat <- cellchatCov
cellchatnonCov <- readRDS("~/cellchatnonCov.rds")
cellchat <- cellchatnonCov
groupSize <- as.numeric(table(cellchat@idents))
groupSize
library(patchwork)
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength",
                 idents.use = 'Platelet',arrow.size = 0.001)
df.net <- read.csv("~/nonCov.df.net.csv")
pathways.show <- df.net$pathway_name
levels(as.factor(pathways.show))
pathways.show
par(mfrow = c(1,2), xpd=TRUE)
levels(cellchat@idents)
vertex.receiver = c(7,10,12,14)
netVisual_aggregate(cellchat, signaling =c("APP"),  
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("CCL"),  
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("CD34"),  
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("CD99"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("CLEC"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("COLLAGEN"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("CXCL"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("CXCL"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("ESAM"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("GP1BA"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("ICAM"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("ITGB2"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("JAM"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("MHC-I"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("RESISTIN"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("SELPLG"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling =c("THBS"),
                    layout = 'circle')
netVisual_aggregate(cellchat, signaling = pathways.show[1], layout = "circle")
netVisual_aggregate(cellchat, signaling = "IL1", layout = "chord")
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
names(group.cellType) <- levels(cellchat@idents)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_cell(cellchat, signaling = pathways.show[1],
                     group = group.cellType,
                     title.name = paste0(pathways.show[1], " signaling network"))
netVisual_chord_cell(cellchat, signaling = pathways.show[1],
                     title.name = paste0(pathways.show[1], " signaling network"))
p1 <- netAnalysis_contribution(cellchat, signaling = pathways.show[1],
                         title = pathways.show[1])
p2 <- netAnalysis_contribution(cellchat, signaling = pathways.show)
cowplot::plot_grid(p1, p2, align = "h",ncol=2)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show[1],
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,]
vertex.receiver = seq(1,4)
p1<-netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = LR.show, vertex.receiver = vertex.receiver,
                     layout = 'hierarchy')
vertex.receiver = seq(1,7)
p2<-netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = LR.show, vertex.receiver = vertex.receiver,
                     layout = 'hierarchy')
cowplot::plot_grid(p1, p2 ,align = "h",ncol=2)
netVisual_individual(cellchat, signaling = pathways.show[1], 
                     pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show, layout = "chord")
pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4)
dir.create('04.pathwat.show')
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0('04.pathwat.show/',pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
gg
levels(cellchat@idents)
par(mfrow = c(1,3), xpd=TRUE)
p1 <- netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c(5:11), remove.isolate = FALSE)
p2 <- netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c("cDC1","cDC2","LC","Inflam. DC","TC","Inflam. TC",  
                                 "CD40LG+ TC"), remove.isolate = FALSE)
p3 <-netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c("cDC1","cDC2"), remove.isolate = FALSE)
p4 <-netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c(5:11), signaling = c("CCL","CXCL"), 
                 remove.isolate = FALSE)
cowplot::plot_grid(p1, p2,p3, align = "h",ncol=4)
par(mfrow = c(1,3), xpd=TRUE)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
pairLR.use
p1 <-netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)
p2 <-netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)
p3 <-netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), 
                 signaling = c("CCL","CXCL"), remove.isolate = T)
cowplot::plot_grid(p1, p2,p3, align = "h",ncol=3)
cellchatCov <- readRDS("~/cellchatCov.rds")
cellchatnonCov <- readRDS("~/cellchatnonCov.rds")
groupSizeCov <- as.numeric(table(cellchatCov@idents))
groupSizeCov
groupSizenonCov <- as.numeric(table(cellchatnonCov@idents))
groupSizenonCov
df.net.ReUCov <- subsetCommunication(cellchatCov)
df.net.ReUnonCov <- subsetCommunication(cellchatnonCov)
pathways.show.ReUCov <- df.net.ReUCov$pathway_name
levels(as.factor(pathways.show.ReUCov))
pathways.show.ReUCovSpecific <- c("ALCAM","CCL","CD6",
                                  "CD86","COLLAGEN","SN")
pathways.show.ReUnonCov <- df.net.ReUnonCov$pathway_name
levels(as.factor(pathways.show.ReUnonCov))
pathways.show.ReUnonCovSpecific <- c("CADM","CD226","GRN",
                                     "IFN-II","MK","NECTIN",
                                     "PARs")
cellchatCov <- readRDS("~/cellchatCov.rds")
cellchatnonCov <- readRDS("~/cellchatnonCov.rds")
groupSizeCov <- as.numeric(table(cellchatCov@idents))
groupSizeCov
groupSizenonCov <- as.numeric(table(cellchatnonCov@idents))
groupSizenonCov
levels(as.factor(pathways.show.ReUCovSpecific))
new.cell.seq <- c("Neutrophil","NK cell",
                  "T cell","B cell","Monocyte",
                  "Erythroblast","PLT","PLT","PLT","PLT","PLT","PLT")
new.cell.seq <- factor(new.cell.seq,c("Neutrophil","NK cell",
                                      "T cell","B cell","Monocyte",
                                      "Erythroblast","PLT"))
names(new.cell.seq) <- levels(reUScellchat@idents)
pathways.show.ReUCovSpecific <- c("ALCAM","CCL","CD6",
                                  "CD86","COLLAGEN","SN")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_aggregate(cellchatCov, signaling ="CCL",
                    layout = 'chord', 
                    vertex.label.cex = 1,
                    title.space=1,pt.title=10)
pathways.show.ReUnonCovSpecific <- c("CADM","CD226","GRN",
                                     "IFN-II","MK","NECTIN",
                                     "PARs")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_aggregate(cellchatnonCov, signaling ="CD226",
                    layout = 'chord', 
                    vertex.label.cex = 1,
                    title.space=1,pt.title=10)
netVisual_aggregate(cellchatnonCov, signaling ="MK",
                    layout = 'chord', 
                    vertex.label.cex = 1,
                    title.space=1,pt.title=10)
netVisual_aggregate(cellchatnonCov, signaling ="NECTIN",
                    layout = 'chord', 
                    vertex.label.cex = 1,
                    title.space=1,pt.title=10)
netVisual_aggregate(cellchatnonCov, signaling ="PARs",
                    layout = 'chord', 
                    vertex.label.cex = 1,
                    title.space=1,pt.title=10)
netVisual_aggregate(reUScellchat, signaling ="CCL",
                    layout = 'chord', 
                    vertex.label.cex = 1,group = new.cell.seq,
                    title.space=5,big.gap = 20,
                    signaling.name = "CCL in Cov")
netVisual_aggregate(reUScellchat, signaling ="CD226",
                    layout = 'chord', 
                    vertex.label.cex = 1,group = new.cell.seq,
                    title.space=5,big.gap = 20,
                    signaling.name = "CCL in Cov")
netVisual_aggregate(reUScellchat, signaling ="NECTIN",
                    layout = 'chord', 
                    vertex.label.cex = 1,group = new.cell.seq,
                    title.space=5,big.gap = 20,
                    signaling.name = "CCL in Cov")
netVisual_aggregate(reUScellchat, signaling ="PARs",
                    layout = 'chord', 
                    vertex.label.cex = 1,group = new.cell.seq,
                    title.space=5,big.gap = 20,
                    signaling.name = "CCL in Cov")
netVisual_aggregate(reUScellchat, signaling ="MK",
                    layout = 'chord', 
                    vertex.label.cex = 1,group = new.cell.seq,
                    title.space=5,big.gap = 20,
                    signaling.name = "CCL in Cov")
pathways.show.ReUCov <- df.net.ReUCov$pathway_name
levels(as.factor(pathways.show.ReUCov))
pathways.show.ReUnonCov <- df.net.ReUnonCov$pathway_name
levels(as.factor(pathways.show.ReUnonCov))
netVisual_aggregate(cellchatCov, signaling ="THBS",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "THBS in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="THBS",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "THBS in nonCov")
netVisual_aggregate(cellchatCov, signaling ="SELPLG",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "SELPLG in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="SELPLG",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "SELPLG in nonCov")
netVisual_aggregate(cellchatCov, signaling ="RESISTIN",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "RESISTIN in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="RESISTIN",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "RESISTIN in nonCov")
netVisual_aggregate(cellchatCov, signaling ="PECAM1",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "JAM in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="PECAM1",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "JAM in nonCov")
netVisual_aggregate(cellchatCov, signaling ="MHC-I",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "JAM in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="MHC-I",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "JAM in nonCov")
netVisual_aggregate(cellchatCov, signaling ="JAM",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "JAM in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="JAM",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "JAM in nonCov")
netVisual_aggregate(cellchatCov, signaling ="ITGB2",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "GP1BA in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="ITGB2",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "GP1BA in nonCov")
netVisual_aggregate(cellchatCov, signaling ="ICAM",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "GP1BA in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="ICAM",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "GP1BA in nonCov")
netVisual_aggregate(cellchatCov, signaling ="GP1BA",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "GP1BA in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="GP1BA",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "GP1BA in nonCov")
netVisual_aggregate(cellchatCov, signaling ="ESAM",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CXCL in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="ESAM",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CXCL in nonCov")
netVisual_aggregate(cellchatCov, signaling ="CXCL",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CXCL in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="CXCL",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CXCL in nonCov")
netVisual_aggregate(cellchatCov, signaling ="CLEC",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CLEC in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="CLEC",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CLEC in nonCov")
netVisual_aggregate(cellchatCov, signaling ="CD99",
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CD99 in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="CD99",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CD99 in nonCov")
netVisual_aggregate(cellchatCov, signaling ="CD34",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CD34 in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="CD34",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CD34 in nonCov")
netVisual_aggregate(cellchatCov, signaling ="CD23",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CD23 in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="CD23",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "CD23 in nonCov")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_aggregate(cellchatCov, signaling ="APP",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "APP in Cov")
netVisual_aggregate(cellchatnonCov, signaling ="APP",  
                    layout = 'chord', vertex.label.cex = 1,
                    title.space=1,pt.title=10,
                    signaling.name = "APP in nonCov")
"THBS","SELPLG","RESISTIN","PECAM1","MHC-I","JAM","ITGB2","ICAM","GP1BA","ESAM","CXCL","CLEC","CD99","CD34","CD23","APP"
subsetCellChatnew <- function (
    object, cells.use = NULL, idents.use = NULL, group.by = NULL, 
    invert = FALSE, thresh = 0.05) 
{
    if (!is.null(idents.use)) {
        if (is.null(group.by)) {
            labels <- object@idents
            if (object@options$mode == "merged") {
                message("Use the joint cell labels from the merged CellChat object")
                labels <- object@idents$joint
            }
        }
        else {
            labels <- object@meta[[group.by]]
        }
        if (!is.factor(labels)) {
            labels <- factor(labels)
        }
        level.use0 <- levels(labels)
        level.use <- levels(labels)[levels(labels) %in% unique(labels)]
        if (invert) {
            level.use <- level.use[!(level.use %in% idents.use)]
        }
        else {
            level.use <- level.use[level.use %in% idents.use]
        }
        cells.use.index <- which(as.character(labels) %in% level.use)
        cells.use <- names(labels)[cells.use.index]
    }
    else if (!is.null(cells.use)) {
        labels <- object@idents
        if (object@options$mode == "merged") {
            message("Use the joint cell labels from the merged CellChat object")
            labels <- object@idents$joint
        }
        level.use0 <- levels(labels)
        level.use <- levels(labels)[levels(labels) %in% unique(as.character(labels[cells.use]))]
        cells.use.index <- which(names(labels) %in% cells.use)
    }
    else {
        stop("USER should define either `cells.use` or `idents.use`!")
    }
    cat("The subset of cell groups used for CellChat analysis are ", 
        level.use, "\n")
    if (nrow(object@data) > 0) {
        data.subset <- object@data[, cells.use.index]
    }
    else {
        data.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    if (nrow(object@data.project) > 0) {
        data.project.subset <- object@data.project[, cells.use.index]
    }
    else {
        data.project.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    data.signaling.subset <- object@data.signaling[, cells.use.index]
    meta.subset <- object@meta[cells.use.index, , drop = FALSE]
    if (object@options$mode == "merged") {
        idents <- object@idents[1:(length(object@idents) - 1)]
        group.existing <- level.use0[level.use0 %in% level.use]
        group.existing.index <- which(level.use0 %in% level.use)
        net.subset <- vector("list", length = length(object@net))
        netP.subset <- vector("list", length = length(object@netP))
        idents.subset <- vector("list", length = length(idents))
        names(net.subset) <- names(object@net)
        names(netP.subset) <- names(object@netP)
        names(idents.subset) <- names(object@idents[1:(length(object@idents) - 
            1)])
        images.subset <- vector("list", length = length(idents))
        names(images.subset) <- names(object@idents[1:(length(object@idents) - 
            1)])
        for (i in 1:length(idents)) {
            cat("Update slots object@images, object@net, object@netP, object@idents in dataset ", 
                names(object@idents)[i], "\n")
            images <- object@images[[i]]
            for (images.j in names(images)) {
                values <- images[[images.j]]
                if (images.j %in% c("coordinates")) {
                  values.new <- values[cells.use.index, ]
                  images[[images.j]] <- values.new
                }
                if (images.j %in% c("distance")) {
                  values.new <- values[group.existing.index, 
                    group.existing.index, drop = FALSE]
                  images[[images.j]] <- values.new
                }
            }
            images.subset[[i]] <- images
            net <- object@net
            for (net.j in names(net)) {
                values <- net[[net.j]]
                if (net.j %in% c("prob", "pval")) {
                  values.new <- values[group.existing.index, 
                    group.existing.index, ]
                  net[[net.j]] <- values.new
                }
                if (net.j %in% c("count", "sum", "weight")) {
                  values.new <- values[group.existing.index, 
                    group.existing.index]
                  net[[net.j]] <- values.new
                }
            }
            net.subset[[i]] <- net
            netP = computeCommunProbPathway(net = net.subset[[i]], 
                pairLR.use = object@LR[[i]]$LRsig, thresh = thresh)
            netP$centr = netAnalysis_computeCentrality(net = net.subset[[i]]$prob)
            netP.subset[[i]] <- netP
            idents.subset[[i]] <- idents[[i]][names(idents[[i]]) %in% 
                cells.use]
            idents.subset[[i]] <- factor(idents.subset[[i]], 
                levels = levels(idents[[i]])[levels(idents[[i]]) %in% 
                  level.use])
        }
        idents.subset$joint <- factor(object@idents$joint[cells.use.index], 
            levels = level.use)
    }
    else {
        cat("Update slots object@images, object@net, object@netP in a single dataset...", 
            "\n")
        group.existing <- level.use0[level.use0 %in% level.use]
        group.existing.index <- which(level.use0 %in% level.use)
        images <- object@images
        for (images.j in names(images)) {
            values <- images[[images.j]]
            if (images.j %in% c("coordinates")) {
                values.new <- values[cells.use.index, ]
                images[[images.j]] <- values.new
            }
            if (images.j %in% c("distance")) {
                values.new <- values[group.existing.index, group.existing.index, 
                  drop = FALSE]
                images[[images.j]] <- values.new
            }
        }
        images.subset <- images
        net <- object@net
        for (net.j in names(net)){
          values <- net[[net.j]]
          if (net.j %in% c("prob","pval")){
            values.new <- values[group.existing.index,group.existing.index, ,drop = F]
            net[[net.j]] <- values.new
            }
          if (net.j %in% c("count","sum","weight")){
            values.new <- values[group.existing.index, group.existing.index, drop=F]
            net[[net.j]] <- values.new
          }
        }
        net.subset <- net
        netP = computeCommunProbPathway(net = net.subset, pairLR.use = object@LR$LRsig, 
            thresh = thresh)
        netP$centr = netAnalysis_computeCentrality(net = net.subset$prob)
        netP.subset <- netP
        idents.subset <- object@idents[cells.use.index]
        idents.subset <- factor(idents.subset, levels = level.use)
    }
    object.subset <- methods::new(Class = "CellChat", data = data.subset, 
        data.signaling = data.signaling.subset, data.project = data.project.subset, 
        images = images.subset, net = net.subset, netP = netP.subset, 
        meta = meta.subset, idents = idents.subset, var.features = object@var.features, 
        LR = object@LR, DB = object@DB, options = object@options)
    return(object.subset)
}
p1 <- netAnalysis_contribution(cellchatCov, 
                               signaling = pathways.show.ReUCov,
                               title =  "Cov",signaling.name = T)
p2 <- netAnalysis_contribution(cellchatnonCov, 
                               signaling = pathways.show.ReUnonCov,
                               title =  "nonCov",signaling.name = T)
cowplot::plot_grid(p1, p2,align = "h",ncol=1)
levels(as.factor(cellchatCov@meta$predicted.celltypeBasicwithnewreU))
cellchatCovnoUnk <- subsetCellChatnew(
  cellchatCov,idents.use = c("B","CD14 Mono","CD16 Mono","CD4 T",
  "CD8 T","DC","Eryth","gdT","HSPC","MAIT","Neutrophil",
  "NK","Plasmablast","Platelet","T Prolif.","Treg"))
cellchatnonCovnoUnk <- subsetCellChatnew(
  cellchatnonCov,idents.use = c("B","CD14 Mono","CD16 Mono","CD4 T",
  "CD8 T","DC","Eryth","gdT","HSPC","MAIT","Neutrophil",
  "NK","Plasmablast","Platelet","T Prolif.","Treg"))
p1 <- netAnalysis_contribution(cellchatCovnoUnk, 
                               signaling = pathways.show.ReUCov,
                               title =  "Cov",signaling.name = T)
p2 <- netAnalysis_contribution(cellchatnonCovnoUnk, 
                               signaling = pathways.show.ReUnonCov,
                               title =  "nonCov",signaling.name = T)
cowplot::plot_grid(p1, p2,align = "h",ncol=1)
pairLR.CLEC <- extractEnrichedLR(cellchatCovnoUnk, 
                                signaling = "CLEC",
                                geneLR.return = T)
pairLR.CLEC$pairLR
LR.show <- pairLR.CLEC$pairLR[4,]
netVisual_individual(cellchatCov.common, signaling = "CLEC",  
                     pairLR.use = LR.show, 
                     layout = 'chord')
netVisual_individual(cellchatnonCovnoUnk, signaling = "CLEC",  
                     pairLR.use = LR.show, 
                     layout = 'chord')
LR.show <- pairLR.CXCL$pairLR[2,]
netVisual_individual(reUScellchat, signaling = "CXCL",  
                     pairLR.use = LR.show, 
                     layout = 'chord',group = new.cell.seq)
netVisual_individual(reUScellchat, signaling = "CXCL", 
                     pairLR.use = LR.show, layout = "circle")
levels(reUHcellchat@idents)
par(mfrow = c(1,3), xpd=TRUE)
p1 <- netVisual_bubble(cellchatCovnoUnk, sources.use = "Platelet", 
                 remove.isolate = FALSE)
pairLR.use <- extractEnrichedLR(cellchatCovnoUnk, signaling = c("CCL"))
pairLR.use
p1 <-netVisual_bubble(cellchatCovnoUnk, sources.use = "Platelet",
                      pairLR.use = pairLR.use, remove.isolate = TRUE)
p2 <-netVisual_bubble(cellchatCovnoUnk, sources.use = "Platelet",
                      signaling = c("CCL","CXCL"), remove.isolate = FALSE)
p3 <-netVisual_bubble(cellchatCovnoUnk, sources.use = "Platelet",
                      signaling = c("CCL","CXCL"), remove.isolate = F)
cowplot::plot_grid(p1, p2,p3, align = "h",ncol=3)
p1
p2
p3
netVisual_bubble(reUScellchat, 
                      signaling = c("THBS","CXCL"), remove.isolate = F)
pairLR.use.Cov <- extractEnrichedLR(cellchatCovnoUnk, signaling = c("CLEC"))
pairLR.use.Cov
pairLR.use.nonCov <- extractEnrichedLR(cellchatnonCovnoUnk, signaling = c("CLEC"))
pairLR.use.nonCov
p1 <- netVisual_bubble(cellchatCovnoUnk, 
                 pairLR.use = pairLR.use.Cov, remove.isolate = TRUE)
p2 <- netVisual_bubble(cellchatnonCovnoUnk, 
                 pairLR.use = pairLR.use.nonCov, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, align = "v",nrow=2)
netVisual_bubble(cellchatnonCovnoUnk, sources.use = "Platelet",
                 pairLR.use = pairLR.use.nonCov, remove.isolate = TRUE)
netVisual_bubble(cellchat,comparison = c(1,2),sources.use = "Platelet",
                 remove.isolate = T,thresh=0.05)
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
VlnPlot(CovintegratedreUall,features=c("CD99","CD99L2"),group.by = "CCClable1",pt.size=0,ncol=2,split.by="cond.idents")+
  theme(legend.position="bottom")
FeaturePlot(CovintegratedreUall,features=c("CD99","CD99L2"),pt.size=0,ncol=2,split.by="cond.idents")+
  theme(legend.position="bottom")
DotPlot(CovintegratedreUall,features=c("CD99"),group.by = "CCClable1",split.by="cond.idents")
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
DotPlot(PLTCovMSH.cluster.final,features=c("CD99"),group.by = "CCClable1")
DotPlot(PLTCovMSH.cluster.final,features=c("CD99"),group.by = "cond.idents")
FeaturePlot(PLTCovMSH.cluster.final,features=c("CD99"),pt.size=0,ncol=2,split.by="cond.idents")+
  theme(legend.position="bottom")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("CD99"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("APP"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("APP"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("APP"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("CCL"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("CCL"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("CCL"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("MHC-I"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("MHC-I"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("MHC-I"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("HLA-E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
  DotPlot(immPLT,features=c("HLA−E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(NFKBPLT,features=c("HLA−E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTLympho,features=c("HLA−E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTx,features=c("HLA−E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(proliferPLT,features=c("HLA−E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
  
DotPlot(aggrePLT,features=c("HLA−A"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(immPLT,features=c("HLA−A"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(NFKBPLT,features=c("HLA−A"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTLympho,features=c("HLA-A"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTx,features=c("HLA−A"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(proliferPLT,features=c("HLA−A"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
  
DotPlot(aggrePLT,features=c("HLA−B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(immPLT,features=c("HLA−B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(NFKBPLT,features=c("HLA−B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTLympho,features=c("HLA−B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTx,features=c("HLA−B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(proliferPLT,features=c("HLA−B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
  
DotPlot(aggrePLT,features=c("HLA−C"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(immPLT,features=c("HLA−C"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(NFKBPLT,features=c("HLA−C"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTLympho,features=c("HLA−C"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(PLTx,features=c("HLA−C"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
  DotPlot(proliferPLT,features=c("HLA−C"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(CovintegratedreUall,features=c("HLA−E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
NK = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("NK cell")]
DotPlot(NK,features=c("KLRD1","KLRC2"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
DotPlot(CovintegratedreUall,features=c("KLRD1","KLRC2"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
Tcell = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("T cell")]
DotPlot(Tcell,features=c("CD8A","CD8B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "MHC-I", split.by = "cond.idents", colors.ggplot = T)
levels(as.factor(CovintegratedreUall$finalgeneralumap))
PLT = CovintegratedreUall[,CovintegratedreUall$finalgeneralumap %in% c("Platelet")]
DotPlot(PLT,features=c("HLA-E"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
DotPlot(PLT,features=c("HLA-A"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
DotPlot(PLT,features=c("HLA-B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
DotPlot(PLT,features=c("HLA-C"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
library(Matrix)
expMtx <- Matrix(PLTCovMSH.cluster.final@assays$RNA@counts,sparse=T)
geneIDs <- data.frame(ID=expMtx@Dimnames[[1]],Name=expMtx@Dimnames[[1]],EX="Gene Expression")
write.table(x=geneIDs,file="4th/PartA/geneIDs.tsv",sep="\t",quote=F,col.names = F, row.names=F)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("CLEC"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("CLEC"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("CLEC"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("CLEC1B"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
NK = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("NK cell")]
DotPlot(NK,features=c("KLRB1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
Tcell = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("T cell")]
DotPlot(Tcell,features=c("KLRB1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "MHC-I", split.by = "cond.idents", colors.ggplot = T)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("GP1BA"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("GP1BA"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("GP1BA"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("GP1BA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
Mono = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("Monocyte")]
DotPlot(Mono,features=c("ITGAM","ITGB2"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
NK = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("NK cell")]
DotPlot(NK,features=c("ITGAM","ITGB2"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "MHC-I", split.by = "cond.idents", colors.ggplot = T)
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("JAM"))
pairLR.use.M
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p2
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("JAM3"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
Mono = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("Monocyte")]
DotPlot(Mono,features=c("ITGAM"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
DotPlot(Mono,features=c("ITGB2"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
NK = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("NK cell")]
DotPlot(NK,features=c("ITGAM"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
DotPlot(NK,features=c("ITGB2"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "JAM", split.by = "cond.idents", colors.ggplot = T)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("CD45"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("CD45"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("CD45"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("CXCL"))
pairLR.use.S
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
p3
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("PPBP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
Neu = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("Neutrophil")]
DotPlot(Neu,features=c("CXCR2"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "CXCL", split.by = "cond.idents", colors.ggplot = T)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("SELPLG"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("SELPLG"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("SELPLG"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("SELP"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
Mono = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("Monocyte")]
DotPlot(Mono,features=c("SELPLG"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "SELPLG", split.by = "cond.idents", colors.ggplot = T)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("RESISTIN"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("RESISTIN"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("RESISTIN"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("CAP1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
Mono = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("Monocyte")]
DotPlot(Mono,features=c("RETN"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "RESISTIN", split.by = "cond.idents", colors.ggplot = T)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("CD22"))
pairLR.use.H
pairLR.use.M <- extractEnrichedLR(reUMcellchat, signaling = c("CD22"))
pairLR.use.M
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("CD22"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p2 <- netVisual_bubble(reUMcellchat, 
                 pairLR.use = pairLR.use.M, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p2, p3, align = "v",nrow=3)
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("PTPRC"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
B = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("B cell")]
DotPlot(B,features=c("CD22"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
plotGeneExpression(reUAllcellchat, signaling = "CD22", split.by = "cond.idents", colors.ggplot = T)
pairLR.use.H <- extractEnrichedLR(reUHcellchat, signaling = c("THBS"))
pairLR.use.H
pairLR.use.S <- extractEnrichedLR(reUScellchat, signaling = c("THBS"))
pairLR.use.S
p1 <- netVisual_bubble(reUHcellchat, 
                 pairLR.use = pairLR.use.H, remove.isolate = TRUE)
p3 <- netVisual_bubble(reUScellchat, 
                 pairLR.use = pairLR.use.S, remove.isolate = TRUE)
cowplot::plot_grid(p1, p3, align = "v",nrow=2)
PLTCovMSH.cluster.final <- readRDS("PLTusingConventionalPC1-11/CT/PLTCovMSH.cluster.withWGCNAwithCT.rds")
levels(as.factor(PLTCovMSH.cluster.final$CCClable1))
aggrePLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("aggrePLT")]
immPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
NFKBPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("NFKBPLT")]
PLTLympho = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTLympho")]
PLTx = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("PLTx")]
proliferPLT = PLTCovMSH.cluster.final[,PLTCovMSH.cluster.final$CCClable1 %in% c("immPLT")]
DotPlot(aggrePLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(immPLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(NFKBPLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTLympho,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(PLTx,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()|
DotPlot(proliferPLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
DotPlot(aggrePLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(immPLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(NFKBPLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTLympho,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(PLTx,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") |
  DotPlot(proliferPLT,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral")
CovintegratedreUall <- readRDS('4th/PartA/steptwoAnno/clusModification/Covint.filtered.clusteredAllmoreres.anno.rds')
levels(as.factor(CovintegratedreUall$CCClable1))
Mono = CovintegratedreUall[,CovintegratedreUall$CCClable1 %in% c("Monocyte")]
DotPlot(Mono,features=c("CD47","THBS1"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") +
  NoLegend()
DotPlot(Mono,features=c("SIRPA"),group.by = "cond.idents")+
  scale_colour_distiller(palette = "Spectral") 
plotGeneExpression(reUAllcellchat, signaling = "THBS", split.by = "cond.idents", colors.ggplot = T)
getwd()
setwd('Mt-Cleared-Outputs/Secondary_preprocess_codes_output/CCC-1st')
getwd()
library(CellChat)
library(patchwork)
library(cowplot)
cellchatCov <- readRDS("~/cellchatCov.rds")
cellchatnonCov <- readRDS("~/cellchatnonCov.rds")
cellchatCov <- readRDS("~/cellchatCovwithcluX.rds")
cellchatnonCov <- readRDS("~/cellchatnonCovwithcluX.rds")
temp <- cbind(levels(cellchatCov@idents),levels(cellchatnonCov@idents))
temp
a <- levels(cellchatCov@idents)
str(a)
idents.use <- subset(a, a != "DC")
idents.use <- subset(idents.use, idents.use != "Neutrophil")
str(idents.use)
cellchatCov.common <- subsetCellChatnew(
  cellchatCov, 
  idents.use = idents.use)
cellchatCov.common@net$centr <-
  netAnalysis_computeCentrality(object=cellchatCov.common,
                                slot.name = "net")
cellchatnonCov@net$centr <-
  netAnalysis_computeCentrality(object=cellchatnonCov,
                                slot.name = "net")
cellchatCov.common <-
  netAnalysis_computeCentrality(object=cellchatCov.common,
                                slot.name = "netP")
cellchatnonCov <-
  netAnalysis_computeCentrality(object=cellchatnonCov,
                                slot.name = "netP")
object.list <- list(Cov = cellchatCov.common, 
                    nonCov = cellchatnonCov)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
saveRDS(cellchat,"~/cellchatAcrosswithcluX.rds")
cellchat <- readRDS("~/cellchatAcrosswithcluX.rds")
cellchatCov <- readRDS("~/cellchatCov.rds")
cellchatnonCov <- readRDS("~/cellchatnonCov.rds")
cellchatCovnoUnk <- subsetCellChatnew(
  cellchatCov,idents.use = c("B","CD14 Mono","CD16 Mono","CD4 T",
  "CD8 T","DC","Eryth","gdT","HSPC","MAIT","Neutrophil",
  "NK","Plasmablast","Platelet","T Prolif.","Treg"))
cellchatnonCovnoUnk <- subsetCellChatnew(
  cellchatnonCov,idents.use = c("B","CD14 Mono","CD16 Mono","CD4 T",
  "CD8 T","DC","Eryth","gdT","HSPC","MAIT","Neutrophil",
  "NK","Plasmablast","Platelet","T Prolif.","Treg"))
temp <- cbind(levels(cellchatCovnoUnk@idents),levels(cellchatnonCovnoUnk@idents))
temp
a <- levels(cellchatCovnoUnk@idents)
str(a)
idents.use <- subset(a, a != "DC")
idents.use <- subset(idents.use, idents.use != "Neutrophil")
idents.use
cellchatCov.common <- subsetCellChatnew(cellchatCovnoUnk, 
                                     idents.use = idents.use)
cellchatCov.common@net$centr <-
  netAnalysis_computeCentrality(object=cellchatCov.common,
                                slot.name = "net")
cellchatnonCovnoUnk@net$centr <-
  netAnalysis_computeCentrality(object=cellchatnonCovnoUnk,
                                slot.name = "net")
object.list <- list(Cov = cellchatCov.common, 
                    nonCov = cellchatnonCovnoUnk)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
saveRDS(cellchat,"~/cellchatAcrossnoUnk.rds")
cellchat <- readRDS("~/cellchatAcrossnoUnk.rds")
cellchatUnk <- cellchat
cellchat <- readRDS("~/cellchatAcross.rds")
cellchat <- readRDS("~/cellchatAcrosswithcluX.rds")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1, 2), stacked = T, do.stat = TRUE,do.flip = T,targets.use = "Platelet")
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1, 2), stacked = F, do.stat = TRUE,do.flip = T)
gg1 + gg2
gg2
par(mfrow = c(1,2), xpd=TRUE)
levels(cellchat@idents$joint)
netVisual_diffInteractionnew(cellchat,comparison = c(1, 2),
                             weight.scale = T
                             )
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteractionnew(cellchat,comparison = c(1, 2), 
                             weight.scale = T,
                             targets.use = "Platelet",arrow.size = 0,
                             vertex.size.max = 4)
netVisual_diffInteractionnew(cellchat, comparison = c(1, 2), 
                             weight.scale = T,
                             sources.use = "Platelet",arrow.size = 0,
                             vertex.size.max = 4
)
netVisual_diffInteractionnew(cellchat,comparison = c(1, 2), 
                          weight.scale = T,
                          targets.use = "Platelet",
                          measure = "weight"
                          )
netVisual_diffInteractionnew(cellchat, comparison = c(1, 2), 
                          weight.scale = T,
                          sources.use = "Platelet",
                          measure = "weight"
                          )
netVisual_diffInteractionnew(cellchat, 
                          comparison = c(1, 2),
                          weight.scale = T, 
                          measure = "weight")
netVisual_diffInteractionnew(cellchat, 
                          comparison = c(1, 2),
                          weight.scale = T, 
                          measure = "weight")
"THBS","SELPLG","RESISTIN","PECAM1","MHC-I","JAM","ITGB2","ICAM","GP1BA","ESAM","CXCL","CLEC","CD99","CD34","CD23","APP"
pathways.show <- c("THBS") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
pathways.show <- c("CXCL")
weight.max <- getMaxWeight(object.list2, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_aggregate(object.list2[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list2)[i]))
}
pathways.show <- c("THBS") 
par(mfrow = c(1,3), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pathways.show <- c("SELPLG") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
netVisual_heatmap(cellchat, signaling = "CLEC", measure = "weight",
                  font.size = 6,font.size.title = 10,comparison = c(2,1))
pathways.show <- c("CCL") 
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
netVisual_diffInteractionnew <- function (object, comparison = c(1, 2), measure = c("count", 
    "weight", "count.merged", "weight.merged"), color.use = NULL, 
    # color.edge = c("
    sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
    top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
    vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
    edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
    label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
    edge.curved = 0.2, shape = "circle", layout = in_circle(), 
    margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
{
    options(warn = -1)
    measure <- match.arg(measure)
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure %in% c("count", "count.merged")) {
        if (is.null(title.name)) {
            title.name = "Differential number of interactions"
        }
    }
    else if (measure %in% c("weight", "weight.merged")) {
        if (is.null(title.name)) {
            title.name = "Differential interaction strength"
        }
    }
    net <- net.diff
    if ((!is.null(sources.use)) | (!is.null(targets.use))) {
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- rownames(net.diff)[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- rownames(net.diff)[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        cells.level <- rownames(net.diff)
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0
        net <- tapply(df.net[["value"]], list(df.net[["source"]], 
            df.net[["target"]]), sum)
        net[is.na(net)] <- 0
    }
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)
        idx2 <- which(Matrix::colSums(net) == 0)
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
    }
    net[abs(net) < stats::quantile(abs(net), probs = 1 - top, 
        na.rm = T)] <- 0
    g <- graph_from_adjacency_matrix(net, mode = "directed", 
        weighted = T)
    edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
    coords <- layout_(g, layout)
    if (nrow(coords) != 1) {
        coords_scale = scale(coords)
    }
    else {
        coords_scale <- coords
    }
    if (is.null(color.use)) {
        color.use = scPalette(length(igraph::V(g)))
    }
    if (is.null(vertex.weight.max)) {
        vertex.weight.max <- max(vertex.weight)
    }
    vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
        5
    
    loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
        2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
        2]/coords_scale[igraph::V(g), 1]))
    if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle <- NA
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
    }
    igraph::V(g)$size <- vertex.weight
    igraph::V(g)$color <- color.use[igraph::V(g)]
    igraph::V(g)$frame.color <- color.use[igraph::V(g)]
    igraph::V(g)$label.color <- vertex.label.color
    igraph::V(g)$label.cex <- vertex.label.cex
    if (label.edge) {
        igraph::E(g)$label <- igraph::E(g)$weight
        igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
    }
    igraph::E(g)$arrow.width <- arrow.width
    igraph::E(g)$arrow.size <- arrow.size
    igraph::E(g)$label.color <- edge.label.color
    igraph::E(g)$label.cex <- edge.label.cex
    igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
        color.edge[2])
    igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
        alpha.edge)
    igraph::E(g)$weight <- abs(igraph::E(g)$weight)
    if (is.null(edge.weight.max)) {
        edge.weight.max <- max(igraph::E(g)$weight)
    }
    if (weight.scale == TRUE) {
        igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
            edge.width.max
    }
    else {
        igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    }
    if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
        igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
            1])] <- loop.angle[edge.start[which(edge.start[, 
            2] == edge.start[, 1]), 1]]
    }
    radian.rescale <- function(x, start = 0, direction = 1) {
        c.rotate <- function(x) (x + start)%%(2 * pi) * direction
        c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }
    label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
        direction = -1, start = 0)
    label.dist <- vertex.weight/max(vertex.weight) + 2
    plot(g, edge.curved = edge.curved, vertex.shape = shape, 
        layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
        vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
        edge.label.family = "Helvetica")
    if (!is.null(title.name)) {
        text(0, 1.5, title.name, cex = 1.1)
    }
    gg <- recordPlot()
    return(gg)
}
"THBS","SELPLG","RESISTIN","PECAM1","MHC-I","JAM","ITGB2",
"ICAM","GP1BA","ESAM","CXCL","CLEC","CD99","CD34","CD23","APP"
plotGeneExpression(cellchat, signaling = "THBS", split.by = "Covidsplit", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "SELPLG", split.by = "cond.idents", colors.ggplot = T)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) +
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}
gg
gg1 <- netAnalysis_signalingRole_scatter(cellchatCovnoUnk,title = "Cov",
                                         weight.MinMax = weight.MinMax,
                                         point.shape = 2)
gg2 <- netAnalysis_signalingRole_scatter(cellchatnonCovnoUnk,title = "nonCov",
                                         weight.MinMax = weight.MinMax)
patchwork::wrap_plots(plots = gg)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat,
                                            idents.use = "Platelet")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat3, idents.use = "Macrophage")
patchwork::wrap_plots(plots = list(gg1,gg2))
netAnalysis_signalingChanges_scatter(cellchat3,
                                            idents.use = "NK_cell")
netAnalysis_signalingChanges_scatter(cellchat3, idents.use = "M",
                                            signaling.exclude = c("MIF"))
suppressMessages(library(ComplexHeatmap))
i = 1
pathway.union <- union(object.list3[[i]]@netP$pathways, object.list3[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list3[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list3)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list3[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list3)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
BigMPLTintegrated <- readRDS("~/BigMPLTintegratedProjected.rds")
BigMPLTintegrated@meta.data$PLT.idents
VlnPlot(BigMPLTintegrated,features=c("KLRB1","CLEC1B","THBS1"),
        pt.size=0,ncol=2,split.by="MSH.idents",
        group.by = "PLT.idents")+
  theme(legend.position="bottom")
DotPlot(BigMPLTintegrated,features=c("KLRB1"),
        group.by = "MSH.idents")+
  theme(legend.position="bottom")
DefaultAssay(BigMPLTintegrated) <- "RNA"
BigMPLTintegrated@assays
tempdataC <- subset(BigMPLTintegrated,CLEC1B>0)
VlnPlot(tempdataC,features=c("CLEC1B"),
        pt.size=0,ncol=2,
        group.by = "MSH.idents")+
  theme(legend.position="bottom")
library(magrittr)
StackedVlnPlot(BigMPLTintegrated,group.by = c("MSH.idents"), features = c("KLRB1","CLEC1B","THBS1","CD47"),pt.size = 1)+
  theme(legend.position="bottom")+
  geom_point(position="jitter",col="black")
AllCovMSH.integrated@meta.data
FeaturePlot(BigMPLTintegrated,reduction = 'ref.umap',
            features = "KLRB1",pt.size=0.00000000001)
FeaturePlot(BigMPLTintegrated,reduction = 'ref.umap',
            features = "CLEC1B",pt.size=0.01)
FeaturePlot(BigMPLTintegrated,reduction = 'ref.umap',
            features = "THBS1",pt.size=0.01)
FeaturePlot(BigMPLTintegrated,reduction = 'ref.umap',
            features = "CD47",pt.size=0.01)
DotPlot(BigMPLTintegrated,
        features = c("KLRB1","CD47","CLEC1B","THBS1"),
        group.by = "PLT.idents")
PLTCovMSH.final <- readRDS("PLTusingConventionalPC1-11/F2a/PLTCovMSH.final.rds")
table(PLTCovMSH.final$CCClable1)
DotPlot(PLTCovMSH.final,
        features = c("KLRB1","CD47","CLEC1B","THBS1"),
        group.by = "CCClable1")
DotPlot(PLTCovMSH.final,
        features = c("KLRB1"),
        group.by = "CCClable1")
FeaturePlot(AllCovMSH.integrated,reduction = 'umap',features = "SIRPA",split.by = 'cond.idents',pt.size=0.00000000001)
VlnPlot(AllCovMSH.integrated,features=c("SIRPA","METTL3"),pt.size=0,ncol=2,split.by="cond.idents",group.by = c("SRlabelsmain"))+
  theme(legend.position="bottom")
