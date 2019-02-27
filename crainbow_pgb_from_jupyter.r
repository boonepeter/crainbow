
if(!requireNamespace("BiocManager")) install.packages("BiocManager")

lib.list=c("Seurat","plyr","ggpubr","dplyr","ggthemes","data.table","GSVA","biomaRt", "GSA")

for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			BiocManager::install(lib.list[i])
			library(lib.list[i],character.only=T)
    };
}

source("useful_functions.R")
source("plotfns.R")

## Read the list of mouse genesets
ach = read.gmt.file("acharya_genesets_mouse.gmt")

ach$geneset.names

pam50 = scan("pam50_genes.txt",what="")
pam50_mm = convertHumanGeneList(pam50)

ach$geneset.descriptions[51]="PAM50_mouse_homologs"
ach$geneset.names[51]="PAM50"
ach$genesets[[51]]=pam50_mm
names(ach$genesets)[51]="PAM50"

files = c("./data/MT1.tsv","./data/TU1.tsv","./data/TU2.tsv","./data/TU3.tsv","./data/TU4.tsv")

data_list = llply(1:length(files),.progress="time",function(i) read.delim(files[i],header=T,row.names=1))


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

cca_out = RunMultiCCA(seu_list,genes.use=Reduce("intersect",gene_list),num.ccs = 20)

DimPlot(object = cca_out, reduction.use = "cca", group.by = "group", pt.size = 1, do.return = F)

VlnPlot(object = cca_out, features.plot = "CC1", group.by = "group", do.return = F)

MetageneBicorPlot(cca_out, grouping.var = "group", dims.eval = 1:20, display.progress = TRUE)

cca_out = CalcVarExpRatio(cca_out,reduction.type = "pca", grouping.var = "group", dims.use = 1:10)
cca_out = SubsetData(cca_out, subset.name = "var.ratio.pca",accept.low = 0.5)

cca_out = AlignSubspace(cca_out, reduction.type = "cca", grouping.var = "group", dims.align = 1:10)

CalcAlignmentMetric(cca_out,reduction.use = "cca.aligned",dims.use = 1:10, grouping.var =  "group")

set.seed(1234)
cca_out = RunTSNE(cca_out, reduction.use = "cca.aligned", dims.use = 1:10, do.fast = T)
TSNEPlot(object = cca_out, do.label = F,group.by="group")

cca_out = FindClusters(cca_out, reduction.type = "cca.aligned", resolution = c(0.6), dims.use = 1:10,print.output = F);
table(cca_out@ident)


TSNEPlot(object = cca_out, do.label =T,no.legend = FALSE,dark.theme = FALSE,label.size=5)

FeaturePlot(cca_out, features.plot = c("Krt5","Krt14","Krt6a"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out, features.plot = c("Krt8","Krt18","Krt19"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out, features.plot = c("Il6","Il12b","Serpinb2","Cxcl3"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out, features.plot = c("Arg1","Mrc1","Egr2","Cd83"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out, features.plot = c("Ptprc","Foxp3","Cd3d","Cd4","Cd14"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

ssgsea_out = gsva(as.matrix(cca_out@data),gset.idx.list=ach$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)

FM = ssgsea_out
FM = as.matrix(Matrix::t(scale(Matrix::t(FM))))
FM = melt(data.frame(t(FM),"Cluster"=cca_out@ident,"Group"=cca_out@meta.data$group))
colnames(FM)[3:4]=c("GeneSets","EnrichmentScore")
head(FM)

ES_mean = FM %>% dplyr::group_by(Cluster,GeneSets,Group) %>% summarize_all(funs(mean))
ggplot(ES_mean, aes(Cluster, GeneSets)) + geom_tile(aes(fill = EnrichmentScore), colour = "white") + scale_fill_gradient2(low ="dark blue", high ="dark red", mid ="white",space = "Lab") +
      theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90)) + facet_grid(~Group)

ES1_mean = FM[,-1] %>% dplyr::group_by(GeneSets,Group) %>% summarize_all(funs(mean))
ggplot(ES1_mean, aes(Group, GeneSets)) + geom_tile(aes(fill = EnrichmentScore), colour = "white") + scale_fill_gradient2(low ="dark blue", high ="dark red", mid ="white",space = "Lab") +
      theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90))

cca_out@meta.data = cbind(cca_out@meta.data,t(as.matrix(Matrix::t(scale(Matrix::t(ssgsea_out))))))
colnames(cca_out@meta.data)

FeaturePlot(cca_out,features.plot=c("BASAL","LUMINAL","MASC","STROMAL","LUMINAL PROGENITORS","EMT"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out,features.plot=c("MACROPHAGE_ACTIVITY","Treg","INFLAMMATORY RESPONSE","NK_CELLS","TCELL_ACTIVATION","G2M CHECKPOINT"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out,features.plot=c("HYPOXIA","INVASIVENESS GENE SIGNATURE","ANGIOGENESIS","ECM","FIBROBLAST GROWTH FACTOR RECEPTOR","EPIGENETIC STEM CELL"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out,features.plot=c("WNT BETA CATENIN SIGNALING","TGF BETA SIGNALING","E2F TARGETS","P53 PATHWAY","HEDGEHOG SIGNALING","PI3K-AKT-MTOR SIGNALING"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out,features.plot=c("BASE EXCISION REPAIR","MISMATCH EXCISION REPAIR","NUCLEOTIDE EXCISION REPAIR","HOMOLOGOUS RECOMBINATION","NONHOMOLOGOUS END JOINING"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

FeaturePlot(cca_out,features.plot=c("PAM50"),reduction.use="tsne",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)


sessionInfo()

num_clusters <- 11
for(i in 0:(num_clusters-1)){
    markers <- FindMarkers(object = cca_out, ident.1 = i, min.pct = 0.25)
    filename <- paste("./out/cluster", i, ".tsv")
    write.table(markers, file = filename, sep = "\t", col.names = T, row.names = T, quote = F)
}
