
# 1. Head --------------------------------------------------------------------

# set working directory, amend as necessary
setwd("/media/data/Josh/NK_scRNA_code")
# set date/time, this is used as a prefix to the output file names
date_time <- gsub("-", "", Sys.Date())


# 2. Read files in -----------------------------------------------------------

# File matching ENS IDs with HUGO gene names
hugo_ENS_match <- read.table("HUGO_ENS_matching.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)

# GTF file for hg19 concatenated with ERCC92
hg19_ercc_gtf <- "hg19_ercc92.gtf"
library(rtracklayer)
hg19_ercc_gtf <- readGFF(hg19_ercc_gtf, filter = list(type = "gene"), tags = c("gene_id", "gene_name"))
gene.name.lookup <- hg19_ercc_gtf$gene_name
names(gene.name.lookup) <- hg19_ercc_gtf$gene_id

# file containing hg19 ENS IDs and gene lengths, needed for RPKM conversion
hg19_genelength <- read.table("hg19_genelengths.txt", header =  TRUE, sep = "\t", as.is = TRUE)

# read tdata
tdata <- read.table("annotations_and_STAR.txt", header =  TRUE, 
                    sep = "\t", row.names = 1, as.is = TRUE)
dim(tdata)
sample_order <- rownames(tdata)

# read counts
counts <- read.table("counts.txt", header =  TRUE, 
                    sep = "\t", row.names = 1, as.is = TRUE)
colnames(counts) <- gsub("^X", "", colnames(counts)) # remove "X" added to colnames while reading table in

# confirm that technical df and count matrix matches
identical(rownames(tdata), colnames(counts)) # TRUE

#   3. QC ----------------------------------------------------------------------

# distribution of cases
type_order <- c("PBMC", "NLN", "PRI", "MET")
tdata$type <- factor(tdata$type, levels = type_order)
table(tdata$pt_id, tdata$type)

pdf(file=paste(date_time, "_QC_plots.pdf", sep=""),height=4,width=5)

# STAR percentaage uniquely mapped
hist(tdata$uniquely_mapped_percent, breaks = 100, xlim=range(0,100), main="Percent of Reads Uniquely Mapped", xlab="Percent")
abline(v=40, col = "red", lwd = 2)

# Number of assigned reads
hist(tdata$Assigned, breaks = 100, main="Number of Assigned Reads", xlab = "No. of Reads")
abline(v=100000, col ="red", lwd = 2)

# ERCC fraction 
ERCC_names <- rownames(counts)[grep("ERCC", rownames(counts))] # 92 ERCC features
ENS_names <- rownames(counts)[grep("ENS", rownames(counts))] # 57820 ENS features

ERCC_total_m <- colSums(counts[ERCC_names,]) # library size for ERCC
identical(names(ERCC_total_m), rownames(tdata)) #TRUE
tdata$ERCC_total <- ERCC_total_m

ENS_total_m <- colSums(counts[ENS_names,]) # library size for hg19 transcriptome
identical(names(ENS_total_m), rownames(tdata)) #TRUE
tdata$ENS_total <- ENS_total_m

tdata$ERCC_fraction <- tdata$ERCC_total / tdata$Assigned #calculate ERCC fraction
hist(tdata$ERCC_fraction, breaks = 100, main="Fraction of Library Assigned to ERCC", xlab = "Fraction")
abline(v=0.8, col = "red", lwd=2)


# Number of ENS features 
tdata$ENS_features <- colSums(counts[ENS_names,] != 0)
hist(tdata$ENS_features, breaks=100, main = "Number of Features per Library", xlab="No. of Features")
abline(v=1200, col = "red", lwd=2)
mean(tdata$ENS_features) # mean before filtering
mean(tdata$ENS_features[tdata$ENS_features >= 1200]) # mean after filtering

# Number of HUGO genes
tdata$HUGO_features <- colSums(counts[ENS_names[(sub("\\..*", "", ENS_names) %in% hugo_ENS_match$Ensembl.ID.supplied.by.Ensembl.)],] != 0)
hist(tdata$HUGO_features, breaks=100, main = "Number of Features per Library", xlab="No. of Features")
abline(v=1200, col = "red", lwd=2)
mean(tdata$HUGO_features) # mean before filtering
mean(tdata$HUGO_features[tdata$HUGO_features >= 1200]) # mean after filtering

# Mitochondrial fraction
# mitochondrial ENS IDs from https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#expression-qc-reads
hemberg_mito_ids <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
                      "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
                      "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
                      "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
                      "ENSG00000198840")
mito_ids <- rownames(counts)[(sub("\\..*", "", rownames(counts)) %in% hemberg_mito_ids)]
mito_total_m <- colSums(counts[mito_ids,])
identical(names(mito_total_m), rownames(tdata)) #TRUE
tdata$mito_total <- mito_total_m
tdata$mito_fraction <- tdata$mito_total / tdata$Assigned
hist(tdata$mito_fraction, breaks = 100, main = "Fraction of Library Assigned to Mitochondrial Genes", xlab = "Fraction")
abline(v=0.1, col = "red", lwd=2)

dev.off()


#   4. Convert to RPKM --------------------------------------------------------------

library(edgeR)
table(colSums(counts)>0) # need to remove one case where library size is 0
counts_rpkm <- counts[ENS_names,colSums(counts)>0]
counts_rpkm <- DGEList(counts=counts_rpkm, genes=hg19_genelength[,c("GeneID","Length")])
counts_rpkm <- as.data.frame(rpkm(counts_rpkm,counts_rpkm$genes$Length))
dim(counts_rpkm) 


# 5. Export Artis & Colonna in RPKM to create signature matrix  --------

# get gene lengths by gene names, needed for Artis data
genelength <- hg19_annotation[,c("GeneID","Length")]
gene.id <- genelength$GeneID
gene.name <- gene.name.lookup[gene.id]
table(duplicated(gene.name))
genelength$GeneID[duplicated(gene.name) == FALSE] <- gene.name[duplicated(gene.name) == FALSE]

# Artis
artis_counts <- read.table("Artis_2019_GSE126107_counts.txt", header = TRUE, sep = "\t", quote="\"", row.names = 1, as.is = TRUE)
# get genelength for genes found in Artis matrix
common_genelist <- rownames(artis_counts)[rownames(artis_counts) %in% genelength$GeneID] #21957 genes
artis_counts <- artis_counts[common_genelist,]
artis_genelength <- genelength[genelength$GeneID %in% common_genelist,]
rownames(artis_genelength) <- artis_genelength$GeneID
artis_genelength <- artis_genelength[common_genelist,]
library(edgeR)
artis_rpkm <- DGEList(counts=artis_counts, genes=artis_genelength[,c("GeneID","Length")])
artis_rpkm <- as.data.frame(rpkm(artis_rpkm,artis_rpkm$genes$Length))
# inclue only features with HGNC symbols
artis_rpkm <- artis_rpkm[rownames(artis_rpkm) %in% hugo_ENS_match$Approved.symbol,]
artis_rpkm <- cbind(rownames(artis_rpkm),artis_rpkm)
colnames(artis_rpkm)[1] <- "Gene"
write.table(artis_rpkm, file="Artis_rpkm_hgnc.txt", 
            sep = "\t", quote = FALSE, row.names=FALSE, col.names=TRUE)


# Colonna 
colonna_matrix <- read.table("Colonna_2019_paper_log2rpkm.txt", header = T, sep = "\t", as.is = TRUE)
# inclue only features with HGNC symbols
colonna_matrix <- colonna_matrix[colonna_matrix$external_gene_source == "HGNC Symbol",]
colonna_matrix <- colonna_matrix[colonna_matrix$external_gene_name %in% hugo_ENS_match$Approved.symbol,] # for consistency with other datasets
colonna_matrix$external_gene_name[duplicated(colonna_matrix$external_gene_name)] # these are the genes that are duplicated
colonna_matrix <- colonna_matrix[!duplicated(colonna_matrix$external_gene_name),] # remove duplicated genes

colonna_matrix_temp <- colonna_matrix[,c(6:21)] # remove extra columns
colonna_matrix_temp <- 2^colonna_matrix_temp # convert back to rpkm instead of log2rpkm
rownames(colonna_matrix_temp) <- colonna_matrix$external_gene_name
colonna_matrix_rpkm <- as.data.frame(cbind(colonna_matrix$external_gene_name, colonna_matrix_temp))
write.table(colonna_matrix_rpkm, file="Colonna_rpkm_hgnc.txt", sep = "\t", quote = FALSE, row.names=FALSE)


# 6. Seurat analysis for all cells meeting filter criteria -----------------------------------------------------------

prefix <- paste(date_time, "_s_ALL_", sep="")

# filter criteria
filter_ALLcells <- tdata$uniquely_mapped_percent>40 & tdata$Assigned>100000 & 
  tdata$ERCC_fraction<0.8 & tdata$mito_fraction<0.1 & tdata$ENS_features>1200
table(filter_ALLcells) # no_of_cells = 1960

cases_filter_ALLcells <- rownames(tdata)[filter_ALLcells]
t_ALL <- tdata[cases_filter_ALLcells,]
m_ALL <- counts[ENS_names,cases_filter_ALLcells]
identical(rownames(t_ALL),colnames(m_ALL)) # TRUE

# replace ENS IDs with gene names when not duplicated
gene.id <- rownames(m_ALL)
gene.name <- gene.name.lookup[gene.id]
table(duplicated(gene.name))
rownames(m_ALL)[duplicated(gene.name) == FALSE] <- gene.name[duplicated(gene.name) == FALSE]
# write.table(m_ALL, file=paste(prefix, "matrix.txt", sep=""), 
#             sep = "\t", quote = FALSE, row.names=TRUE, col.names=NA)

# may need to restart R - maximum number of DLLs reached
library(Seurat)
library(mclust)
library(dplyr)
library(ggplot2)

s_ALL <- CreateSeuratObject(
  counts = m_ALL,
  min.cells = 3, 
  min.features = 200
)

s_ALL <- NormalizeData(
  object = s_ALL,
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

#   6.1 Find Variable Genes, PCA -----------------------------------------------------
pdf(file=paste(prefix, "findvargenes.pdf", sep=""),height=7.5,width=10)
s_ALL <- FindVariableFeatures(
  object = s_ALL,
  selection.method = "vst",
  nfeatures = 2000
)
dev.off()

# scale data
all.genes <- rownames(s_ALL)
s_ALL <- ScaleData(object = s_ALL, features = all.genes) 
# linear dimension reduction
s_ALL <- RunPCA(s_ALL, features = VariableFeatures(object= s_ALL))

pdf(file=paste(prefix, "PCAs.pdf", sep=""),height=6,width=8)
VizDimLoadings(s_ALL, dims = 1:2, reduction = "pca")
DimPlot(object = s_ALL, reduction = "pca")
dev.off()

pdf(file=paste(prefix, "PC_Heatmaps.pdf", sep=""),height=12,width=15)
DimHeatmap(s_ALL, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

s_ALL <- JackStraw(object = s_ALL, num.replicate = 100 )
s_ALL <- ScoreJackStraw(s_ALL, dims = 1:20)

pdf(file=paste(prefix, "Jackstraw.pdf", sep=""),height=12,width=15)
JackStrawPlot(object = s_ALL, dims = 1:20)
ElbowPlot(s_ALL)
dev.off()


#   6.2 Find clusters ------------------------------------------------
s_ALL_dim <- 1:12
s_ALL <- FindNeighbors(s_ALL, dims = s_ALL_dim)
s_ALL <- FindClusters(s_ALL, resolution = 0.9)

s_ALL <- RunUMAP(s_ALL, dims = 1:12)
pdf(file=paste(prefix, "umap_unlabelled.pdf", sep=""),height=5,width=6)
DimPlot(s_ALL, dims=c(1:2), reduction = "umap")
dev.off()

# find markers for each cluster
s_ALL_markers <- FindAllMarkers(s_ALL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
s_ALL_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top_markers <- s_ALL_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf(file=paste(prefix, "hm_markers.pdf", sep=""),height=12,width=15)
DoHeatmap(s_ALL, features = top_markers$gene) +
  theme(axis.text.y = element_text(size = 3))
dev.off()

# based on the top markers, clusters 5, 8 and 9 comprise CD3+/CD4+ T-cells and CD86+ antigen presenting cells
markers_NK_vs_nonNK <- FindMarkers(s_ALL, ident.1=c("0","1","2","3","4","6","7"), ident.2=c("5","8","9"))
write.table(markers_NK_vs_nonNK, file=paste(prefix, "markers_NK_vs_nonNK.txt", sep=""), 
            sep = "\t", quote = FALSE, row.names=TRUE, col.names=TRUE)

#   6.3 Show NK vs non-NK clusters --------------------------------------------

cluster_identity_all <- FetchData(s_ALL, vars = "ident")
dim(cluster_identity_all)
identical(rownames(cluster_identity_all), rownames(t_ALL)) #TRUE
cluster_identity_all$cell_type <- t_ALL$type

library(dplyr)
library(gdata)
cluster_identity_all$new_labels <- NA
cluster_identity_all$new_labels <- recode(cluster_identity_all$ident,
                                      '0'="C1",
                                      '1'="C6",
                                      '2'="C2",
                                      '3'="C5",
                                      '4'="C4",
                                      '5'="C8",
                                      '6'="C7",
                                      '7'="C3",
                                      '8'="C10",
                                      '9'="C9")

levels_10 <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10")
cluster_identity_all$new_labels <- factor(cluster_identity_all$new_labels, levels=levels_10)
cluster_identity_all$cell_type <- drop.levels(cluster_identity_all$cell_type)

# get sample IDs for those cells in the NK clusters only *** this is needed for section 6 where we only look at NK cells alone
samples_nkclusters_only <- row.names(cluster_identity_all)[!cluster_identity_all$new_labels %in% c("C8", "C9", "C10")]
samples_nonnkclusters <- row.names(cluster_identity_all)[cluster_identity_all$new_labels %in% c("C8", "C9", "C10")]


# plot for new labels
s_ALL_lab <- s_ALL
Idents(s_ALL_lab) <- cluster_identity_all$new_labels
s_ALL_lab <- AddMetaData(s_ALL_lab, cluster_identity_all$new_labels, col.name = "new_labels")
colors_10 <- c("thistle3",
               "thistle2",
               "thistle1",
               "lavender",
               "lavenderblush2",
               "lavenderblush3",
               "darkgrey",
               "burlywood2",
               "burlywood3",
               "burlywood4")

pdf(file=paste(prefix, "umap_labels.pdf", sep=""),height=4,width=4.5)
DimPlot(s_ALL_lab, reduction = "umap", cols = colors_10, label = TRUE, label.size = 4) +
  geom_hline(yintercept = 4.45, linetype="dashed") +
  annotate(geom="text", x=4, y=-6, label="NK clusters", color="black", size=4.5) +
  annotate(geom="text", x=3.5, y=10, label="non-NK clusters", color="black", size=4.5)
dev.off()


# plot by celltype
s_ALL_celltype <- s_ALL
Idents(s_ALL_celltype) <- cluster_identity_all$cell_type
colors4 <- c("azure4","deepskyblue3","plum3","deeppink3")
pdf(file=paste(prefix, "umap_celltype.pdf", sep=""),height=4,width=4.8)
DimPlot(s_ALL_celltype, reduction = "umap", cols = colors4)
dev.off()

# plot by FACS
s_ALL_facs <- s_ALL
Idents(s_ALL_facs) <- t_ALL$FACS_code
colors3 <- c("deepskyblue3","deeppink3", "plum3")
pdf(file=paste(prefix, "umap_facs.pdf", sep=""),height=4,width=5.5)
DimPlot(s_ALL_facs, reduction = "umap", cols = colors3)
dev.off()

# plot by pt ID
s_ALL_ptid <- s_ALL
Idents(s_ALL_ptid) <- t_ALL$pt_id
colors10 <- c("pink","red","green3",
              "orange","chocolate4","darkorchid", 
              "cadetblue2", "darkgreen","azure3","azure4")
pdf(file=paste(prefix, "umap_sampleid.pdf", sep=""),height=4,width=4.8)
DimPlot(s_ALL_ptid, reduction = "umap", cols = colors10)
dev.off()

# markers for new labels
# find markers for each cluster (repeating for new clusternames)
s_ALL_markers <- FindAllMarkers(s_ALL_lab, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
s_ALL_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(s_ALL_markers, 
            file=paste(prefix, "markers.txt", sep=""), 
            sep = "\t",
            quote = FALSE,
            row.names=TRUE,
            col.names=TRUE)

top_markers <- s_ALL_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf(file=paste(prefix, "hm_markers_labels.pdf", sep=""),height=12,width=15)
DoHeatmap(s_ALL_lab, features = top_markers$gene) +
  theme(axis.text.y = element_text(size = 3),
        axis.text.x = element_text(size = 8))
dev.off()

# featureplots
fea_plot <- c("CD3G", "CD3D", "CD4","FOXP3","CD86","HLA-DMB") 
pdf(file=paste(prefix, "featureplot.pdf", sep=""),height=7.5,width=6)
FeaturePlot(object = s_ALL_lab, features = fea_plot, pt.size = 0.3, 
            cols = c("gray92", "red"), sort.cell=F, max.cutoff = 3, ncol = 2) 
dev.off()

# 7. Seurat NK & ILC cells only -----------------------------------------------------------

prefix <- paste(date_time, "_s_NK_", sep="")
filter_NKcells <- tdata$uniquely_mapped_percent>40 & tdata$Assigned>100000 & 
  tdata$ERCC_fraction<0.8 & tdata$mito_fraction<0.1 & tdata$ENS_features>1200 & 
  rownames(tdata) %in% samples_nkclusters_only
table(filter_NKcells) # no_of_cells = 1683

# Assign the right filter
cases_filter_NKcells <- rownames(tdata)[filter_NKcells]
# filtered matrices
t_NK <- tdata[cases_filter_NKcells,]
dim(t_NK)
m_NK <- counts[ENS_names,cases_filter_NKcells]

# replace ENS IDs with gene names when not duplicated
gene.id <- rownames(m_NK)
gene.name <- gene.name.lookup[gene.id]
table(duplicated(gene.name))
rownames(m_NK)[duplicated(gene.name) == FALSE] <- gene.name[duplicated(gene.name) == FALSE]

library(Seurat)
library(mclust)
library(dplyr)
library(ggplot2)

s_NK <- CreateSeuratObject(
  counts = m_NK,
  min.cells = 3, 
  min.features = 200
)

s_NK <- NormalizeData(
  object = s_NK,
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

#   7.1 Find Variable Genes, PCA -----------------------------------------------------
pdf(file=paste(prefix, "findvargenes.pdf", sep=""),height=7.5,width=10)
s_NK <- FindVariableFeatures(
  object = s_NK,
  selection.method = "vst",
  nfeatures = 2000
)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(s_NK)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot2
dev.off()

# scale data
all.genes <- rownames(s_NK)
s_NK <- ScaleData(object = s_NK, features = all.genes)

# linear dimension reduction
s_NK <- RunPCA(s_NK, features = VariableFeatures(object= s_NK))
print(s_NK[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file=paste(prefix, "PCAs.pdf", sep=""),height=6,width=8)
VizDimLoadings(s_NK, dims = 1:2, reduction = "pca")
DimPlot(object = s_NK, reduction = "pca")
dev.off()

pdf(file=paste(prefix, "PC_Heatmaps.pdf", sep=""),height=12,width=15)
DimHeatmap(s_NK, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

s_NK <- JackStraw(object = s_NK, num.replicate = 100 )
s_NK <- ScoreJackStraw(s_NK, dims = 1:20)

pdf(file=paste(prefix, "Jackstraw.pdf", sep=""),height=12,width=15)
JackStrawPlot(object = s_NK, dims = 1:20)
ElbowPlot(s_NK)
dev.off()

#   7.2 Find clusters ------------------------------------------------
s_NK_dim <- 1:12
s_NK_resolution <- 0.8
s_NK <- FindNeighbors(s_NK, dims = s_NK_dim)
s_NK <- FindClusters(s_NK, resolution = s_NK_resolution)

s_NK <- RunUMAP(s_NK, dims = s_NK_dim)
pdf(file=paste(prefix, "umap_unlabelled.pdf", sep=""),height=5,width=6)
DimPlot(s_NK, dims=c(1:2), reduction = "umap")
dev.off()

# find markers for each cluster
s_NK_markers <- FindAllMarkers(s_NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- s_NK_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf(file=paste(prefix, "hm_markers.pdf", sep=""),height=12,width=15)
DoHeatmap(s_NK, features = top_markers$gene) +
  theme(axis.text.y = element_text(size = 3))
dev.off()

# cluster_identities
cluster_ident_NK <- FetchData(s_NK, vars = "ident")
cluster_ident_NK$cell_type <- t_NK$type
cluster_ident_NK$facs <- t_NK$FACS_code
cluster_ident_NK$pt_id <- t_NK$pt_id

# plot by celltype
s_NK_celltype <- s_NK
Idents(s_NK_celltype) <- cluster_ident_NK$cell_type
colors4 <- c("lightskyblue3","orange","mistyrose3","mistyrose4")
pdf(file=paste(prefix, "umap_others_celltype.pdf", sep=""),height=3,width=5.2)
DimPlot(s_NK_celltype, reduction = "umap", pt.size=0.6) +
  labs(colour="Cell Origin") +
  scale_color_manual(labels = c("Peripheral blood", "Normal lymph node", "Primary tumor", 
                                "Metastatic lymph node"),values=colors4)
dev.off()

# plot by FACS
FACS_levels <- c("CD56+ CD127-", "CD56+ CD127+", "CD56- CD127+")
cluster_ident_NK$facs <- factor(cluster_ident_NK$facs, levels=FACS_levels)

s_NK_facs <- s_NK
Idents(s_NK_facs) <- cluster_ident_NK$facs
colors3 <- c("deepskyblue3","deeppink3", "plum3")
pdf(file=paste(prefix, "umap_others_facs.pdf", sep=""),height=3,width=4.8)
DimPlot(s_NK_facs, reduction = "umap", cols = colors3, pt.size = 0.6) +
  labs(colour="FACS Indexing")
dev.off()

# plot by pt ID
ptid_levels <- c("7202", "7204", "7205", "7207", "7211", "7257", "7267", "7329", "N001", "N002")
cluster_ident_NK$pt_id <- factor(cluster_ident_NK$pt_id, levels = ptid_levels)

s_NK_ptid <- s_NK
Idents(s_NK_ptid) <- cluster_ident_NK$pt_id
colors10 <- c("red","green3",
              "pink","chocolate4","cadetblue2", 
              "mediumpurple1", "darkgreen","orange", "azure3","azure4")

pdf(file=paste(prefix, "umap_others_source.pdf", sep=""),height=3,width=4.4)
DimPlot(s_NK_ptid, reduction = "umap", cols = colors10, pt.size=0.4) +
  labs(colour="Source Patient")
dev.off()

#   7.3 Deconvolution with CIBERSORTx --------------------------

# export single-cell matrix as rpkm
rpkm_NK <- counts_rpkm[ENS_names,cases_filter_NKcells]
dim(rpkm_NK)
identical(rownames(t_NK),colnames(rpkm_NK)) # TRUE
# replace ENS IDs with Gene Names
gene.id <- rownames(rpkm_NK)
gene.name <- gene.name.lookup[gene.id]
table(duplicated(gene.name))
rownames(rpkm_NK)[duplicated(gene.name) == FALSE] <- gene.name[duplicated(gene.name) == FALSE]
# for consistency, only HUGO gene symbols are allowed
rpkm_NK_hgnc <- rpkm_NK[rownames(rpkm_NK) %in% hugo_ENS_match$Approved.symbol,]
write.table(rpkm_NK_hgnc, file=paste(prefix, "rpkm_hgnc.txt", sep=""), sep = "\t", quote = FALSE, row.names=TRUE, col.names=NA)


library(tidyr)
library(ggplot2)

# Deconvolution by signature from Yudanin et al. (Artis)
cib_artis_abs_NK <- read.table("deconvoluted_artis.txt", header =  TRUE, 
                     sep = "\t", row.names = 1, as.is = TRUE)
# get relative fractions
cib_artis_rel_NK <- cib_artis_abs_NK
cib_artis_rel_NK$NK <- cib_artis_rel_NK$NK / cib_artis_rel_NK$Absolute.score..sig.score.
cib_artis_rel_NK$ILC1 <- cib_artis_rel_NK$ILC1 / cib_artis_rel_NK$Absolute.score..sig.score.
cib_artis_rel_NK$ILC2 <- cib_artis_rel_NK$ILC2 / cib_artis_rel_NK$Absolute.score..sig.score.
cib_artis_rel_NK$ILC3 <- cib_artis_rel_NK$ILC3 / cib_artis_rel_NK$Absolute.score..sig.score.

cell_levels <- c("NK", "ILC1", "ILC2", "ILC3")
cell_colors <- c("lightpink2","lightskyblue2","lavenderblush2","palegreen2")
untidy <- as.data.frame(cib_artis_rel_NK[,cell_levels])
untidy$cluster_label <- cluster_ident_NK[row.names(untidy),]$ident
tidied <- untidy %>%
  gather(cell_levels, key="celltype", value="fraction")
tidied <- tidied[!is.na(tidied$cluster_label),]
tidied$celltype <- factor(tidied$celltype, levels=cell_levels)
p <- ggplot(tidied, aes(celltype,fraction,fill = cluster_label))
p + 
  geom_bar(position = 'dodge', stat = 'summary', fun.y = 'median') +
  labs(title = "Deconvolution (Yudanin et al. 2019)", x = "Yudanin et al. signature", y = "Fraction", fill = "") +
  # scale_fill_manual(labels = levels_NK, values = colors_NK) +
  geom_point(aes(x = celltype), shape = 21, cex=0.2, alpha=0.3, position = position_jitterdodge(jitter.width = 0.3,jitter.height=0, dodge.width=0.9)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=12), axis.text=element_text(size=12))


# Deconvolution by signature from Collins et al. (Colonna)
cell_levels <- c("CD57_pos", "CD57_neg", "CD56_bright", "ieILC1", "ILC3")
cell_colors <- c("lightpink2","lightskyblue2","palegreen2","antiquewhite2","antiquewhite3")

cib_colonna_abs_NK <- read.table("deconvoluted_colonna.txt", header =  TRUE, 
                               sep = "\t", row.names = 1, as.is = TRUE)

# the "ILC1s" in Colonna are actually "ieILC1s"
colnames(cib_colonna_abs_NK)[4] <- "ieILC1"
# get relative fractions
cib_colonna_rel_NK <- cib_colonna_abs_NK
cib_colonna_rel_NK$CD56_bright <- cib_colonna_rel_NK$CD56_bright / cib_colonna_rel_NK$Absolute.score..sig.score.
cib_colonna_rel_NK$ieILC1 <- cib_colonna_rel_NK$ieILC1 / cib_colonna_rel_NK$Absolute.score..sig.score.
cib_colonna_rel_NK$ILC3 <- cib_colonna_rel_NK$ILC3 / cib_colonna_rel_NK$Absolute.score..sig.score.
cib_colonna_rel_NK$CD57_neg <- cib_colonna_rel_NK$CD57_neg / cib_colonna_rel_NK$Absolute.score..sig.score.
cib_colonna_rel_NK$CD57_pos <- cib_colonna_rel_NK$CD57_pos / cib_colonna_rel_NK$Absolute.score..sig.score.

library(tidyr)
library(ggplot2)

untidy <- as.data.frame(cib_colonna_rel_NK[,cell_levels])
untidy$cluster_label <- cluster_ident_NK[row.names(untidy),]$ident
tidied <- untidy %>%
  gather(cell_levels, key="celltype", value="fraction")
tidied <- tidied[!is.na(tidied$cluster_label),]
tidied$celltype <- factor(tidied$celltype, levels=cell_levels)
p <- ggplot(tidied, aes(celltype,fraction,fill = cluster_label))
p + 
  geom_bar(position = 'dodge', stat = 'summary', fun.y = 'median') +
  labs(title = "Deconvolution (Collins et al. 2019)", x = "Collins et al. signature", y = "Fraction", fill = "") +
  # scale_fill_manual(labels = levels_NK, values = colors_NK) +
  geom_point(aes(x = celltype), shape = 21, cex=0.2, alpha=0.3, position = position_jitterdodge(jitter.width = 0.3,jitter.height=0, dodge.width=0.9)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=12), axis.text=element_text(size=12))


# 7.4 Label Clusters ------------------------------------------------------

# Assign new cluster names
library(dplyr)
cluster_ident_NK$new_labels <- NA
cluster_ident_NK$new_labels <- recode(cluster_ident_NK$ident,
                                      '0'="Peripheral",
                                      '1'="ILC1",
                                      '2'="ieILC1",
                                      '3'="NK-1",
                                      '4'="NK-2",
                                      '5'="NK-ILC1 intermediate",
                                      '6'="ieILC1 cycling",
                                      '7'="ILC3")

levels_NK <- c("Peripheral",
               "NK-1",
               "NK-2",
               "NK-ILC1 intermediate",
               "ILC1",
               "ILC3",
               "ieILC1",
               "ieILC1 cycling")
cluster_ident_NK$new_labels <- factor(cluster_ident_NK$new_labels, levels=levels_NK)

# plot for new labels
s_NK_lab <- s_NK
Idents(s_NK_lab) <- cluster_ident_NK$new_labels
s_NK_lab <- AddMetaData(s_NK_lab, cluster_ident_NK$new_labels, col.name = "new_labels")
colors_NK <- c("snow4",
               "orange",
               "green3",
               "cyan",
               "cornflowerblue",
               "blue",
               "hotpink2",
               "red")

pdf(file=paste(prefix, "umap_labels.pdf", sep=""),height=4,width=4.5)
DimPlot(s_NK_lab, reduction = "umap", cols = colors_NK)
dev.off()

# recode with new numbers assigned to each cluster
cluster_ident_NK$num_labels <- recode(cluster_ident_NK$new_labels,
                                      'Peripheral'='1: Peripheral',
                                      'ILC1'='5: ILC1',
                                      'ieILC1'='7: ieILC1',
                                      'NK-1'='2: NK-1',
                                      'NK-2'='3: NK-2',
                                      'NK-ILC1 intermediate'='4: NK-ILC1 intermediate',
                                      'ieILC1 cycling'='8: ieILC1 cycling',
                                      'ILC3'='6: ILC3')
levels_NK_num <- c('1: Peripheral', '2: NK-1', '3: NK-2', '4: NK-ILC1 intermediate', 
                   '5: ILC1', '6: ILC3', '7: ieILC1', '8: ieILC1 cycling')

cluster_ident_NK$num_labels <- factor(cluster_ident_NK$num_labels, levels=levels_NK_num)

s_NK_num <- s_NK_lab
Idents(s_NK_num) <- cluster_ident_NK$num_labels

pdf(file=paste(prefix, "umap_labels_numbered.pdf", sep=""),height=4,width=6)
DimPlot(s_NK_num, reduction = "umap", cols = colors_NK, label=FALSE, 
        label.size = 7, repel = FALSE)
dev.off()


# find markers for each cluster (repeating for new clusternames)
s_NK_markers <- FindAllMarkers(s_NK_lab, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
s_NK_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(s_NK_markers, 
            file=paste(prefix, "markers.txt", sep=""), 
            sep = "\t",
            quote = FALSE,
            row.names=TRUE,
            col.names=NA)

top_markers <- s_NK_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

pdf(file=paste(prefix, "hm_markers_labels.pdf", sep=""),height=4,width=8)
DoHeatmap(s_NK_lab, features = top_markers$gene, 
          group.bar = TRUE, group.colors = colors_NK, 
          group.bar.height = 0.03, label = FALSE,
          draw.lines=TRUE, lines.width = 6) +
  theme(axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 8))
dev.off()


# Replot deconvolution plots to show concordance between cluster identities and reference signatures

# Yudanin et al. (Artis)
cell_levels <- c("NK", "ILC1", "ILC2", "ILC3")
cell_colors <- c("lightpink2","lightskyblue2","lavenderblush2","palegreen2")

pdf(file=paste(prefix, "cib_Artis.pdf", sep=""),height=3,width=9)
untidy <- as.data.frame(cib_artis_rel_NK[,cell_levels])
untidy$cluster_label <- cluster_ident_NK[row.names(untidy),]$new_labels
tidied <- untidy %>%
  gather(cell_levels, key="celltype", value="fraction")
tidied <- tidied[!is.na(tidied$cluster_label),]
tidied$celltype <- factor(tidied$celltype, levels=cell_levels)
p <- ggplot(tidied, aes(celltype,fraction,fill = cluster_label))
p + 
  geom_bar(position = 'dodge', stat = 'summary', fun.y = 'median') +
  labs(title = "Deconvolution (Yudanin et al. 2019)", x = "Yudanin et al. signature", y = "Fraction", fill = "") +
  scale_fill_manual(labels = levels_NK, values = colors_NK) +
  geom_point(aes(x = celltype), shape = 21, cex=0.2, alpha=0.3, position = position_jitterdodge(jitter.width = 0.3,jitter.height=0, dodge.width=0.9)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=12), axis.text=element_text(size=12))
dev.off()

# Collins et al. (Colonna)
cell_levels <- c("CD57_pos", "CD57_neg", "CD56_bright", "ieILC1", "ILC3")
cell_colors <- c("lightpink2","lightskyblue2","palegreen2","antiquewhite2","antiquewhite3")

pdf(file=paste(prefix, "cib_Colonna.pdf", sep=""),height=3,width=9)
untidy <- as.data.frame(cib_colonna_rel_NK[,cell_levels])
untidy$cluster_label <- cluster_ident_NK[row.names(untidy),]$new_labels
tidied <- untidy %>%
  gather(cell_levels, key="celltype", value="fraction")
tidied <- tidied[!is.na(tidied$cluster_label),]
tidied$celltype <- factor(tidied$celltype, levels=cell_levels)
p <- ggplot(tidied, aes(celltype,fraction,fill = cluster_label))
p + 
  geom_bar(position = 'dodge', stat = 'summary', fun.y = 'median') +
  labs(title = "Deconvolution (Collins et al. 2019)", x = "Collins et al. signature", y = "Fraction", fill = "") +
  scale_fill_manual(labels = levels_NK, values = colors_NK) +
  geom_point(aes(x = celltype), shape = 21, cex=0.2, alpha=0.3, position = position_jitterdodge(jitter.width = 0.3,jitter.height=0, dodge.width=0.9)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=12), axis.text=element_text(size=12))
dev.off()

#   7.5 Stacked barplot of tumors --------------------------------------
clust_ident_NK_tum <- cluster_ident_NK[!cluster_ident_NK$new_labels=="Peripheral" & 
                                         !cluster_ident_NK$cell_type %in% c("PBMC", "NLN"),]
library(dplyr)
stack_NK <- clust_ident_NK_tum %>% 
  group_by(pt_id,new_labels) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

stack_NK$pt_id <- as.character(stack_NK$pt_id)

colors_stack_NK <- c(
  "orange",
  "green3",
  "cyan2",
  "cornflowerblue",
  "blue",
  "hotpink2",
  "red")
pdf(file=paste(prefix, "stacked_barplot.pdf", sep=""),height=4,width=6)
ggplot(stack_NK, aes(x = factor(pt_id), y = perc*100, fill = factor(new_labels))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(title = "Distribution of Cell Types by Tumor", x = "Patient", y = "Percent", fill = "Cell Type") +
  scale_fill_manual(breaks = levels_NK[2:length(levels_NK)], values= colors_stack_NK) +
  theme_minimal() +
  theme(axis.title.y=element_text(size=12), axis.title.x=element_text(size=12), axis.text=element_text(size=12)) 
dev.off()

#     7.6 Pseudotime  ----------------------------------------------------------
library(monocle)
seuratX <- s_NK_lab
# remove ILC1, ILC3 and NK-ILC1 intermediate
seuratX <- SubsetData(seuratX, ident.remove = c("ILC1", "ILC3", "NK-ILC1 intermediate"))

# extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(seuratX@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seuratX@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# construct monocle celldataset
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = uninormal() ) # since data is already normalized, thresholded and scalled in Seurat
pData(HSMM)
fData(HSMM)

# select top 50 genes from each cluster
top_markers_pseudotime <- s_NK_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

# run ordering algorithm
HSMM <- setOrderingFilter(HSMM, top_markers_pseudotime$gene)
# reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM <- reduceDimension(HSMM,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=2,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=0)

## order cells change colors and theta to match your plot
HSMM <- orderCells(HSMM)

levels_NK_new <- levels_NK[c(1,2,3,7,8)]
colors_NK_new <- colors_NK[c(1,2,3,7,8)]

pdf(file=paste(prefix, "pseudotime_trajectory.pdf", sep=""),height=2.2,width=3)
plot_cell_trajectory(HSMM, x = 1, y = 2,
                     color_by = "new_labels",
                     theta = -15,
                     show_branch_points = F,
                     show_tree = T,
                     cell_link_size = 0.5,
                     cell_size = 0.3) + 
  scale_color_manual(breaks = levels_NK_new, values= colors_NK_new) 
plot_cell_trajectory(HSMM, 
                     color_by = "State",
                     theta = -15,
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 0.8) 
dev.off()

branch_activation <- c("ITGAE", "ITGA1", "GZMA", "NR4A2")
pdf(file=paste(prefix, "pseudotime_branched.pdf", sep=""),height=5,width=4)
plot_genes_branched_pseudotime(HSMM[branch_activation,],
                               panel_order = branch_activation,
                               branch_point = 2,
                               #                               branch_labels = c("NK-2", "ieILC1"),
                               method = "loess",
                               color_by = "new_labels",
                               expression_curve_linetype_by = "new_labels",
                               cell_size=0.2,
                               relative_expr = FALSE,
                               ncol = 1) +
  scale_color_manual(values= c("hotpink2", "red", "orange","green3","snow4","grey65", "grey15")) +
  scale_fill_manual(values=c("white", "white"))
dev.off()


#     7.7 Log fold-change plots -------------------------------------------------------

# north
markers_ieILC1_vs_ILC3 <- FindMarkers(s_NK_lab, ident.1="ieILC1", ident.2="ILC3", 
                                      logfc.threshold = 0, min.pct=0, min.cells.group = 0)
colnames(markers_ieILC1_vs_ILC3) <- paste("north", colnames(markers_ieILC1_vs_ILC3), sep="_")
write.table(markers_ieILC1_vs_ILC3, file=paste(prefix, "markers_ieILC1ILC3_directcomparison.txt", sep=""), 
            sep = "\t", quote = FALSE, row.names=TRUE, col.names=NA)

# east
markers_NK2_vs_hILC1 <- FindMarkers(s_NK_lab, ident.1="NK-2", ident.2="ILC1",
                                    logfc.threshold = 0, min.pct=0, min.cells.group = 0)
colnames(markers_NK2_vs_hILC1) <- paste("east", colnames(markers_NK2_vs_hILC1), sep="_")

# include common genes
o <- rownames(markers_NK2_vs_hILC1)[rownames(markers_NK2_vs_hILC1) %in% rownames(markers_ieILC1_vs_ILC3)]
markers_NK2_vs_hILC1 <- markers_NK2_vs_hILC1[o,]
markers_ieILC1_vs_ILC3 <- markers_ieILC1_vs_ILC3[o,]
identical(rownames(markers_ieILC1_vs_ILC3), rownames(markers_NK2_vs_hILC1))

map_CD56_brightdim <- cbind(markers_ieILC1_vs_ILC3, markers_NK2_vs_hILC1)
map_CD56_brightdim$symbol <- rownames(map_CD56_brightdim)
write.table(map_CD56_brightdim, file=paste(prefix, "map_CD56_brightdim.txt", sep=""), 
            sep = "\t", quote = FALSE, row.names=TRUE, col.names=NA)
map_CD56_brightdim_fil <- map_CD56_brightdim[map_CD56_brightdim$north_p_val_adj < 0.05 | map_CD56_brightdim$east_p_val_adj < 0.05,]

pdf(file=paste(prefix, "logFC_plot_A.pdf", sep=""),height=5,width=6)

vp_first <- map_CD56_brightdim
vp <- map_CD56_brightdim_fil
library(calibrate)
logcri <- log10(3)

with(vp_first, plot(east_avg_logFC, north_avg_logFC, pch=20, 
                    main="", col=alpha("grey", 0.8),
                    ylab="Log FC ieILC1 / ILC3", xlab="Log FC NK-2 / ILC1",
                    xlim = c(-2.7, 3.6), ylim=c(-4.6,6)))
abline(v = logcri, col="grey", lwd=2, lty=2)
abline(v = -logcri, col="grey", lwd=2, lty=2)
abline(h = logcri, col="grey", lwd=2, lty=2)
abline(h = -logcri, col="grey", lwd=2, lty=2)
# groups
text(0, 5.9, "ieILC1", cex=1, font=2)
text(0, -4.5, "ILC3", cex=1, font=2)
text(3.3, 0, "NK-2", cex=1, font=2)
text(-2.5, 0, "ILC1", cex=1, font=2)

# ieILC1
with(subset(vp, north_avg_logFC > logcri & abs(east_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("hotpink2",0.4)))
# 
with(subset(vp, east_avg_logFC > logcri & north_avg_logFC > logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("chocolate4",0.4)))
# peripheral
with(subset(vp, east_avg_logFC > logcri & abs(north_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("chartreuse3",0.4)))
# 
with(subset(vp, east_avg_logFC > logcri & north_avg_logFC < -logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("darkolivegreen",0.4)))
# ILC3
with(subset(vp, north_avg_logFC < -logcri & abs(east_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("blue",0.4)))
# 
with(subset(vp, east_avg_logFC < -logcri & north_avg_logFC < -logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("navy",0.4)))
# ILC1
with(subset(vp, east_avg_logFC < -logcri & abs(north_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("cornflowerblue",0.4)))
# 
with(subset(vp, east_avg_logFC < -logcri & north_avg_logFC > logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("plum3",0.4)))
# selected genes
vp_targets <- c("GNLY", "GZMA", "GZMB", "GZMH", "NKG7", "PRF1", "CCL4",
                "ITGA1", "ITGAE", "IFNG",
                "KLRC1", "KIR2DL4",
                "NR4A1","NR4A2", "NR4A3", "CXCR4",
                "IL23R", "RORC",
                "KIT", "CXCR4", "SELL", "FOS",
                "CCR7", "SELL", "IL7R", "KIT")
vp_select <- vp[vp$symbol %in% vp_targets,]
library(basicPlotteR)
addTextLabels(vp_select$east_avg_logFC, vp_select$north_avg_logFC, vp_select$symbol, cex.label = 0.8, 
              col.label="black", col.background = "white")
dev.off()

#   7.8 LogFC plot peripheral NK-1 NK2 ieILC1  -------------------------------------------------------
# north
markers_peripheral_vs_NK2 <- FindMarkers(s_NK_lab, ident.1="Peripheral", ident.2="NK-2", 
                                         logfc.threshold = 0, min.pct=0, min.cells.group = 0)
colnames(markers_peripheral_vs_NK2) <- paste("north", colnames(markers_peripheral_vs_NK2), sep="_")

# east
markers_NK1_vs_ieILC1 <- FindMarkers(s_NK_lab, ident.1="NK-1", ident.2="ieILC1",
                                     logfc.threshold = 0, min.pct=0, min.cells.group = 0)
colnames(markers_NK1_vs_ieILC1) <- paste("east", colnames(markers_NK1_vs_ieILC1), sep="_")

# include common genes
o <- rownames(markers_NK1_vs_ieILC1)[rownames(markers_NK1_vs_ieILC1) %in% rownames(markers_peripheral_vs_NK2)]
markers_NK1_vs_ieILC1 <- markers_NK1_vs_ieILC1[o,]
markers_peripheral_vs_NK2 <- markers_peripheral_vs_NK2[o,]
identical(rownames(markers_peripheral_vs_NK2), rownames(markers_NK1_vs_ieILC1))

map_CD56_brightdim <- cbind(markers_peripheral_vs_NK2, markers_NK1_vs_ieILC1)
map_CD56_brightdim$symbol <- rownames(map_CD56_brightdim)
map_CD56_brightdim_fil <- map_CD56_brightdim[map_CD56_brightdim$north_p_val_adj < 0.05 | map_CD56_brightdim$east_p_val_adj < 0.05,]

# plot logFC map
pdf(file=paste(prefix, "logFC_plot_B", sep=""),height=6,width=7.5)

vp_first <- map_CD56_brightdim
vp <- map_CD56_brightdim_fil
library(calibrate)
logcri <- log10(3)

with(vp_first, plot(east_avg_logFC, north_avg_logFC, pch=20, 
                    main="", col=alpha("grey", 0.8),
                    ylab="Log FC peripheral / NK-2", xlab="Log FC NK-1 / ieILC1",
                    xlim = c(-2, 2.3), ylim=c(-3.4,2.6)))
abline(v = logcri, col="grey", lwd=2, lty=2)
abline(v = -logcri, col="grey", lwd=2, lty=2)
abline(h = logcri, col="grey", lwd=2, lty=2)
abline(h = -logcri, col="grey", lwd=2, lty=2)
# groups
text(0, 2.6, "peripheral", cex=1, font=2)
text(0, -3.4, "NK-2", cex=1, font=2)
text(2.2, 0, "NK-1", cex=1, font=2)
text(-1.9, 0, "ieILC1", cex=1, font=2)

# ieILC1
with(subset(vp, north_avg_logFC > logcri & abs(east_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("hotpink2",0.4)))
# 
with(subset(vp, east_avg_logFC > logcri & north_avg_logFC > logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("chocolate4",0.4)))
# peripheral
with(subset(vp, east_avg_logFC > logcri & abs(north_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("chartreuse3",0.4)))
# 
with(subset(vp, east_avg_logFC > logcri & north_avg_logFC < -logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("darkolivegreen",0.4)))
# ILC3
with(subset(vp, north_avg_logFC < -logcri & abs(east_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("blue",0.4)))
# 
with(subset(vp, east_avg_logFC < -logcri & north_avg_logFC < -logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("navy",0.4)))
# ILC1
with(subset(vp, east_avg_logFC < -logcri & abs(north_avg_logFC) < logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("cornflowerblue",0.4)))
# 
with(subset(vp, east_avg_logFC < -logcri & north_avg_logFC > logcri ), points(east_avg_logFC, north_avg_logFC, pch=20, col=alpha("plum3",0.4)))
# selected genes
vp_targets <- c("S100B", "CX3CR1", "ITGB2", "FGFBP2",
                "SELL", "IFIT1", "IFIT3", "STAT1",
                "IL7R", "CXCR4", 
                "NR4A2", "NFKB1", "REL", "FOS", "NFKBIA", "NFKBIZ",
                "NR4A1", "NR4A3",
                "FOSB", "JUNB", "ZEB2",
                "IFNG", "GZMA", "GZMB",
                "GNLY", "CCND2",
                "GZMK", "ITGA1", "ITGAE", "CSF1")
vp_select <- vp[vp$symbol %in% vp_targets,]
library(basicPlotteR)
addTextLabels(vp_select$east_avg_logFC, vp_select$north_avg_logFC, vp_select$symbol, 
              cex.label = 0.7, 
              col.label="black", col.background = "white")
dev.off()

#   7.9 Feature Plots -------------------------------------------------------

# Figure 2b
pdf(file=paste(prefix, "fp_fig2b.pdf", sep=""),height=7.2,width=8.8)
fea_plot <- c("ITGA1", "ITGAE", "CXCR6", "NCR2",
              "GZMA", "GZMB", "PRF1", "GNLY",
              "NR4A1", "NR4A2", "NR4A3", "CXCR4",
              "IL7R", "SELL", "KIT", "RORC")
p <- FeaturePlot(s_NK_lab, fea_plot, pt.size = 0.1, cols = c("green", "blue"), 
                 combine = FALSE, order = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()
}
cowplot::plot_grid(plotlist = p, ncol=4)
dev.off()

# Supplementary Fig 2
pdf(file=paste(prefix, "fp_supp_fig2.pdf", sep=""),height=7.2,width=8.8)
fea_plot <- c("STAT1", "NR4A2", "NFKB1", "REL", 
              "NFKBIZ", "NFKBIA", "ZEB2", "FOS", "FOSB",
              "JUNB", "AREG", "AHR", "EOMES", 
              "TBX21")
p <- FeaturePlot(s_NK_lab, fea_plot, pt.size = 0.1, cols = c("green", "blue"), 
                 combine = FALSE, order = TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()
}
cowplot::plot_grid(plotlist = p, ncol=4)
dev.off()
