## Set your working directory
setwd("/Users/pgb13/repositories/crainbow/") ## setwd("your directory")

## Install BiocManager if not yet installed
if(!requireNamespace("BiocManager")) install.packages("BiocManager")

## Load the required R libraries
## If the libraries are not installed, use BiocManager
lib.list=c("Seurat","plyr","ggpubr","dplyr","ggthemes","data.table","GSVA", "biomaRt")


for(i in 1:length(lib.list)){
  if(any(installed.packages()[,1]==lib.list[i])){
    library(lib.list[i],character.only=T)}else{
      BiocManager::install(lib.list[i])
      library(lib.list[i],character.only=T)
    };
}

library(GSVA)





## CCA




## R scripts with utility functions
# source("~/Research/scripts/r_scripts/useful_functions.R")
# source("~/Research/scripts/r_scripts/plotfns.R")
# Assuming that these R scripts are in the current working dir.
source("useful_functions.R")
source("plotfns.R")

### Read the gene-set file
ach = read.gmt.file("acharya_genesets_mouse.gmt")

### Read all the possible .tsv files and identity select files
### Create a list object containing gene expression matrices of the select files
#files = list.files(path = "./data/", pattern=".tsv")
#files = files[files %in% c("MT1.tsv","TU1.tsv","TU2.tsv","TU3.tsv","TU4.tsv")]
files = c("./data/MT1.tsv", "./data/TU1.tsv", "./data/TU2.tsv", "./data/TU3.tsv", "./data/TU4.tsv")

data_list = llply(1:length(files),.progress="time",function(i) read.delim(files[i],header=T,row.names=1))

### Create a SEURAT object from each gene expression matrix within the list
### Each SEURAT object is preprocessed separately and highly variable genes are computed
seu_list = gene_list = vector(mode="list",length=length(data_list))
for(i in 1:length(seu_list)){
  print(i)
  seu = CreateSeuratObject(raw.data=data_list[[i]],min.cells=5)
  seu@meta.data$group = gsub(".tsv","",files[i])
  seu = FilterCells(seu, subset.names = "nGene", low.thresholds = round(quantile(seu@meta.data$nGene,0.1)), high.thresholds = Inf)
  seu = NormalizeData(seu,display.progress = F)
  seu = ScaleData(seu, display.progress = F,vars.to.regress = c("nUMI"))
  seu = FindVariableGenes(seu,display.progress = F, do.plot = F, mean.function = ExpMean, dispersion.function = LogVMR)
  genes_use = head(rownames(seu@hvg.info), 2000)
  seu_list[[i]] = seu
  gene_list[[i]] = genes_use
}

### Run multivariate Canonical Correlation Analysis (CCA) on the common genes across all the genesets
### Calculate the ratio of total variance explained by PPCA vs total variance explained by CCA, and filter cells based on these values
### CCA Align the matrices and compute the alignment metric score
### Map a t-SNE plot and find clusters

pdf("CCA_plots.pdf",width=10,height=7)
cca_out = RunMultiCCA(seu_list,genes.use=Reduce("intersect",gene_list),num.ccs = 20)
DimPlot(object = cca_out, reduction.use = "cca", group.by = "group", pt.size = 1, do.return = F)
VlnPlot(object = cca_out, features.plot = "CC1", group.by = "group", do.return = F)
MetageneBicorPlot(cca_out, grouping.var = "group", dims.eval = 1:20, display.progress = FALSE)
DimHeatmap(object = cca_out, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
cca_out = CalcVarExpRatio(cca_out,reduction.type = "pca", grouping.var = "group", dims.use = 1:10)
cca_out = SubsetData(cca_out, subset.name = "var.ratio.pca",accept.low = 0.5)
cca_out = AlignSubspace(cca_out, reduction.type = "cca", grouping.var = "group", dims.align = 1:10)
CalcAlignmentMetric(cca_out,reduction.use = "cca.aligned",dims.use = 1:10, grouping.var =  "group")
cca_out = RunTSNE(cca_out, reduction.use = "cca.aligned", dims.use = 1:10, do.fast = T)
TSNEPlot(object = cca_out, do.label = F,group.by="group")
cca_out = FindClusters(cca_out, reduction.type = "cca.aligned", resolution = c(0.6), dims.use = 1:10,print.output = F)
TSNEPlot(object = cca_out, do.label =T,no.legend = FALSE,dark.theme = FALSE,label.size=5)
FeaturePlot(cca_out, features.plot = c("Krt5","Krt14","Krt6a"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Krt8","Krt18","Krt19"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Il6","Il12b","Serpinb2","Cxcl3"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Arg1","Mrc1","Egr2","Cd83"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Ptprc","Foxp3","Cd3d","Cd4","Cd14"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

### Compute sample-set GSEA (ssGSEA) scores
## Plor


ssgsea_out = gsva(as.matrix(cca_out@data),gset.idx.list=ach$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
pam50 = scan("pam50_genes.txt",what="")
pam50_mm = convertHumanGeneList(pam50)
pam50_ssgsea = gsva(as.matrix(cca_out@data),gset.idx.list=list(pam50_mm),method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
ssgsea_out = rbind(ssgsea_out,pam50_ssgsea)
rownames(ssgsea_out)[49] = "PAM50"
FM = ssgsea_out
FM = as.matrix(Matrix::t(scale(Matrix::t(FM))))
FM = melt(data.frame(t(FM),"Cluster"=cca_out@ident,"Group"=cca_out@meta.data$group))
colnames(FM)[3:4]=c("GeneSets","EnrichmentScore")

### Heatmap of cluster-wise gene-set enrichment scores for every gene-set
ES_mean = FM %>% dplyr::group_by(Cluster,GeneSets,Group) %>% summarize_all(funs(mean))
ggplot(ES_mean, aes(Cluster, GeneSets)) + geom_tile(aes(fill = EnrichmentScore), colour = "white") + scale_fill_gradient2(low ="dark blue", high ="dark red", mid ="white",space = "Lab") +
      theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90)) + facet_grid(~Group)

### Heatmap of cluster-specific gene-set enrichment scores
ES1_mean = FM[,-1] %>% dplyr::group_by(GeneSets,Group) %>% summarize_all(funs(mean))
ggplot(ES1_mean, aes(Group, GeneSets)) + geom_tile(aes(fill = EnrichmentScore), colour = "white") + scale_fill_gradient2(low ="dark blue", high ="dark red", mid ="white",space = "Lab") +
      theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90))

### Gene-set feature plots
cca_out@meta.data = cbind(cca_out@meta.data,t(as.matrix(Matrix::t(scale(Matrix::t(ssgsea_out))))))
FeaturePlot(cca_out,features.plot=c("BASAL","LUMINAL","MASC","STROMAL","LUMINAL PROGENITORS","EMT"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("MACROPHAGE_ACTIVITY","Treg","INFLAMMATORY RESPONSE","NK_CELLS","TCELL_ACTIVATION","G2M CHECKPOINT"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("HYPOXIA","INVASIVENESS GENE SIGNATURE","ANGIOGENESIS","ECM","FIBROBLAST GROWTH FACTOR RECEPTOR","EPIGENETIC STEM CELL"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("WNT BETA CATENIN SIGNALING","TGF BETA SIGNALING","E2F TARGETS","P53 PATHWAY","HEDGEHOG SIGNALING","PI3K-AKT-MTOR SIGNALING"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("BASE EXCISION REPAIR","MISMATCH EXCISION REPAIR","NUCLEOTIDE EXCISION REPAIR","HOMOLOGOUS RECOMBINATION","NONHOMOLOGOUS END JOINING"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("PAM50"),reduction.use="tsne",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
###

###################################################################################################
# Finding differentially expressed genes (cluster biomarkers) and save the results #
# NOTE: Can play around with the min.pct param - make it high - so that canonical markers get an advantage (hopefully)!!!!
###################################################################################################
cca_out.markers <- FindAllMarkers(object = cca_out, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cca_out.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)




print("Entering marker identification")

cluster0.markers <- FindMarkers(object = cca_out, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = cca_out, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = cca_out, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = cca_out, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = cca_out, ident.1 = 4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = cca_out, ident.1 = 5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object = cca_out, ident.1 = 6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object = cca_out, ident.1 = 7, min.pct = 0.25)
cluster8.markers <- FindMarkers(object = cca_out, ident.1 = 8, min.pct = 0.25)
cluster9.markers <- FindMarkers(object = cca_out, ident.1 = 9, min.pct = 0.25)
cluster10.markers <- FindMarkers(object = cca_out, ident.1 = 10, min.pct = 0.25)
cluster11.markers <- FindMarkers(object = cca_out, ident.1 = 11, min.pct = 0.25)


#

write.table(cluster0.markers, "./out/cluster0.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster1.markers, "./out/cluster1.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster2.markers, "./out/cluster2.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster3.markers, "./out/cluster3.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster4.markers, "./out/cluster4.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster5.markers, "./out/cluster5.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster6.markers, "./out/cluster6.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster7.markers, "./out/cluster7.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster8.markers, "./out/cluster8.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster9.markers, "./out/cluster9.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster10.markers, "./out/cluster10.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(cluster11.markers, "./out/cluster11.tsv", sep="\t", col.names=T, row.names=T, quote=F)


print("Leaving marker identification")

cluster1.markers <- FindMarkers(object = cca_out, ident.1 = 1, ident.2 = 7, min.pct = 0.25, only.pos = TRUE)
print(x = head(x = cluster1.markers, n = 20))




dev.off()
