# Script to compare PRAME's gene expression between ADRN/MES clustered cells using scRNA-Seq data
# to be compared only across tumor clusters
# setwd("~/KP/singleCellProjects/PRAME_ARDN_MES")
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

# load PMC (Princess Maxima Center) data -------------------------
load('~/KP/RShiny/single-cell-dts/nbPMC_data.RData')
str(nbPMC_data)

# get annotations
anno <- nbPMC_data@meta.data
# get cluster IDs
ident.info <- as.data.frame(Idents(nbPMC_data))
# merge metadata with indent.info
anno <- merge(anno, ident.info, by = 'row.names')


# to get cluster numbers for tumor clusters
unique(anno[grepl('Tumour', anno$Annotation),"Idents(nbPMC_data)"])
clusters.of.interest <- c('10','24','13')
# to visualize and make sure we are picking the right clusters
p <- DimPlot(nbPMC_data, reduction = "umap", label = TRUE, pt.size = 0.5) # default label = idents (cluster identity)
p <- DimPlot(nbPMC_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'Annotation')
ggsave(p, filename = 'temp.pdf', width = 10, height = 10)


# List of ADRN- and MES-signature genes from Groningen et.al, 2017 -------------
adrn.mes.list <- read.delim('ADRN_MES_signature_genes_list.txt', header = F)
adrn <- adrn.mes.list[adrn.mes.list$V2 == 'ADRN',]
adrn <- as.data.frame(adrn)
mes <- adrn.mes.list[adrn.mes.list$V2 == 'MES',]

# to label clusters as ADRN or MES ------------------
# find all markers for every cluster compared to all remaining cells
# Pre-filter features that have less than a two-fold change between the average expression in all clusters
nbPMC.markers <- FindAllMarkers(nbPMC_data, min.pct = 0.25, logfc.threshold = log(2))


# ...to classify clusters as ADRN/MES, get percentages of genes overlap in each cluster ----------

# get total genes in each cluster group
nbPMC.markers <- nbPMC.markers %>%
  group_by(cluster) %>%
  mutate(N=n())

# get number of overlapping ADRN genes in each group
freq.adrn <- nbPMC.markers %>%
  inner_join(., adrn, by=c("gene" = "V1")) %>% 
  group_by(cluster) %>%
  count()

# add freq.adrn to nbPMC.markers
nbPMC.markers <- merge(nbPMC.markers, freq.adrn, by = 'cluster')
names(nbPMC.markers)[9] <- 'n.adrn'

# to get percent ADRN for each cluster
pct.adrn <- nbPMC.markers %>%
  group_by(cluster) %>%
  summarise(percent_ADRN=n.adrn/N*100) %>%
  distinct()

# get number of overlapping MES genes in each group
freq.mes <- nbPMC.markers %>%
  inner_join(., mes, by=c("gene" = "V1")) %>% 
  group_by(cluster) %>%
  count()

# add freq.adrn to nbPMC.markers
nbPMC.markers <- merge(nbPMC.markers, freq.mes, by = 'cluster')
names(nbPMC.markers)[10] <- 'n.mes'

# to get percent ADRN for each cluster
pct.mes <- nbPMC.markers %>%
  group_by(cluster) %>%
  summarise(percent_MES=n.mes/N*100) %>%
  distinct()

# merge pct.adrn and pct.mes
subgroups <- merge(pct.adrn, pct.mes, by = 'cluster')
# get annotations for non-tumor clusters
anno.subset <- anno %>%
  dplyr::select(3,9) %>%
  distinct()
names(anno.subset) <- c('annotation','cluster')

# merge to get original annotations
subgroups <- merge(subgroups, anno.subset, by = 'cluster')
# replace Tumour Clusters by ADRN
subgroups$annotation <- gsub('Tumour cluster [0-9]','ADRN',subgroups$annotation)

# check cluster IDs (groups)
levels(nbPMC_data)

# assign cluster identities as ADRN/MES ------------------
subgroups <- subgroups[order(subgroups$cluster),]
new.cluster.ids <- subgroups$annotation
names(new.cluster.ids) <- levels(nbPMC_data)
nbPMC_data <- RenameIdents(nbPMC_data, new.cluster.ids)
DimPlot(nbPMC_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'Annotation')

# find all markers of ADRN level
# at least 10% of cells in cluster ADRN or all other cluster must express the gene for it to be included in Differential expression
adrn.de.markers <- FindMarkers(nbPMC_data, ident.1 = 'ADRN', min.pct = 0.25)
adrn.de.markers %>% 
  rownames_to_column(var = 'gene') %>%
  filter(gene == 'PRAME')


# Visualize expression of feature (i.e. PRAME) in ADRN and MES clusters
VlnPlot(nbGOSH_data, features = c("PRAME"), log = TRUE)
FeaturePlot(nbGOSH_data, features = c("PRAME"), label = TRUE)
DotPlot(object = nbGOSH_data, features = 'PRAME')
RidgePlot(nbGOSH_data, features = "PRAME", ncol = 2)
