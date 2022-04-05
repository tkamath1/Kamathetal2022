########################################
## Annotating individual cell types ####
########################################
library(Seurat)
library(rliger)
library(qs)
source("/home/tkamath/DA/PDpaper_reproducibility/genevalidatehelper.R")
source("/home/tkamath/DA/PDpaper_reproducibility/SeuratExtrafunctions.R")
source("/home/tkamath/DA/PDpaper_reproducibility/extrafuncs.R")
source('/home/tkamath/DA/PDpaper_reproducibility/prestowrapper.R')

############################################################
#### First merge all individually annotated datasets #######
############################################################
library(Seurat)
library(tidyverse)
library(NNLM)
library(stringr)
library(doParallel)
library(tidyr)
library(qs)

# First list all files of new objects
setwd('/home/tkamath/DA/individuals/')
list.individuals <- lapply(grep('qs',list.files(),value = T),function(x){qread(x)})
list.individuals <- lapply(list.individuals,function(z){
  z@meta.data$annotation <- z@ident
  return(z)})
names(list.individuals) <- readr::parse_number(grep('qs',list.files(),value = T))

metadata.all <- lapply(list.individuals,function(x){x@meta.data})
names(metadata.all) <- names(list.individuals)
metadata.all <- bind_rows(metadata.all,.id = 'names')
table(metadata.all$names,metadata.all$annotation)

# Loop through major cell types and subset eaach of the datasets
out.raw <- lapply(unique(metadata.all$names),function(x){
  out.raw <- lapply(c('NonDA','DA','Endofibro','Oligo','OPC','Astro','MG','PBMC'),function(y){
    idx.use <- names(list.individuals[[x]]@ident[which(list.individuals[[x]]@ident == y)] )
    return(list.individuals[[x]]@raw.data[,idx.use])
  })
  names(out.raw) <- c('NonDA','DA','Endofibro','Oligo','OPC','Astro','MG','PBMC')
  return(out.raw)
})

names(out.raw) <- unique(metadata.all$names)

sn.nonda <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 1))
qsave(sn.nonda,'../nurrintegration/sn_nonda.qs')

sn.da <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 2))
qsave(sn.da,'../nurrintegration/sn_da.qs')

sn.endo <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 3))
qsave(sn.endo,'../nurrintegration/sn_endo.qs')

sn.olig <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 4))
qsave(sn.olig,'../nurrintegration/sn_olig.qs')

sn.opc <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 5))
qsave(sn.opc,'../nurrintegration/sn_opc.qs')

sn.astro <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 6))
qsave(sn.astro,'../nurrintegration/sn_astro.qs')

sn.mg <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 7))
qsave(sn.mg,'../nurrintegration/sn_mg.qs')

sn.pbmc <- rliger::createLiger(raw.data = lapply(out.raw, `[[`, 8))
qsave(sn.pbmc,'../nurrintegration/sn_pbmc.qs')


####################
# DA neurons #######
####################
source("/home/tkamath/DA/PDpaper_reproducibility/SeuratExtrafunctions.R")
source("/home/tkamath/DA/PDpaper_reproducibility/extrafuncs.R")
source('/home/tkamath/DA/PDpaper_reproducibility/prestowrapper.R')

pd.meta <- read.csv('/home/tkamath/DA/PD_snrna_metadata.csv')
pd.meta <- pd.meta %>% mutate(disease = ifelse(Status %in% c('C','CA'),'Ctrl','Disease'))
pd.meta.use <-pd.meta[match(unique(pd.meta$Donor.ID),pd.meta$Donor.ID),]

# 
sn.da <- qread('/home/tkamath/DA/nurrintegration/sn_da.qs')
sn.da <- AddMito(sn.da,species = 'human')
sn.da@cell.data$status <- pd.meta.use[match(sn.da@cell.data$dataset,pd.meta.use$Donor.ID),]$disease

# Look at controls first
idx.use <- rownames(sn.da@cell.data[which(sn.da@cell.data$status == 'Ctrl'),])
sn.da.ctrl <- rliger::subsetLiger(object = sn.da, cells.use = idx.use)
sn.da.ctrl = rliger::normalize(sn.da.ctrl)
sn.da.ctrl <- addcelldata(sn.da.ctrl)

# Make an object from the rbinded matrices to detect HVG
tmp1 <- createLiger(raw.data = list('total' = MergeSparseDataAll(sn.da.ctrl@raw.data)))
tmp1 <- rliger::normalize(tmp1)
tmp1 <- selectGenes(tmp1, num.genes = 3000, do.plot = T)
genes.use <- lapply(sn.da.ctrl@raw.data, function(x){rownames(x)})
genes.use <- Reduce(intersect,genes.use)
var.genes.use <- intersect(tmp1@var.genes, genes.use)
sn.da.ctrl@var.genes <- var.genes.use

sn.da.ctrl = scaleNotCenter(sn.da.ctrl)
sn.da.ctrl = optimizeALS(sn.da.ctrl, k = 30,lambda = 10)
sn.da.ctrl <- quantile_norm(sn.da.ctrl)
sn.da.ctrl = runUMAP(sn.da.ctrl)


toremove <- c(2,11,21)
sn.da.ctrl <- rliger::subsetLiger(sn.da.ctrl,clusters.use = 
                            setdiff(levels(sn.da.ctrl@clusters),toremove))
sn.da.ctrl <- AddMito(sn.da.ctrl,species = 'human')
sn.da.ctrl <- rliger::subsetLiger(sn.da.ctrl,cells.use = 
                            rownames(sn.da.ctrl@cell.data[sn.da.ctrl@cell.data$percent.mito < 0.1,]))

sn.da.ctrl = optimizeALS(sn.da.ctrl, k = 25,lambda = 10)
sn.da.ctrl <- quantile_norm(sn.da.ctrl)
sn.da.ctrl <- louvainCluster(sn.da.ctrl,resolution = 0.6)
sn.da.ctrl = runUMAP(sn.da.ctrl)

# Remove 7 (seems like interneurons), 17,18
toremove <- c(7,17,18)
sn.da.ctrl <- rliger::subsetLiger(sn.da.ctrl, clusters.use = setdiff(levels(sn.da.ctrl@clusters),toremove))

#Remove 11,13 (oligos)
toremove <- c(11,13)
sn.da.ctrl <- subsetLiger(sn.da.ctrl, clusters.use = setdiff(levels(sn.da.ctrl@clusters),toremove))
sn.da.ctrl = optimizeALS(sn.da.ctrl, k = 25,lambda = 10)
sn.da.ctrl <- quantile_norm(sn.da.ctrl)
sn.da.ctrl <- louvainCluster(sn.da.ctrl,resolution = 0.8)
sn.da.ctrl = runUMAP(sn.da.ctrl)

sn.da.ctrl <- subsetLiger(sn.da.ctrl, clusters.use = setdiff(levels(sn.da.ctrl@clusters),c(6)))
sn.da.ctrl <- quantile_norm(sn.da.ctrl)
sn.da.ctrl <- louvainCluster(sn.da.ctrl,resolution = 0.8)
sn.da.ctrl = runUMAP(sn.da.ctrl)

# Remove 7,10,11,14,15,16,17,18
toremove <- c(7,10,11,14,15,16,17,18)
sn.da.ctrl <- subsetLiger(sn.da.ctrl, clusters.use = setdiff(levels(sn.da.ctrl@clusters),toremove))
sn.da.ctrl <- quantile_norm(sn.da.ctrl)
sn.da.ctrl <- louvainCluster(sn.da.ctrl,resolution = 0.8)
sn.da.ctrl = runUMAP(sn.da.ctrl)

sn.da.seurat <- ligerToSeurat(sn.da.ctrl)
markers.da <-prestowrapper(sn.da.seurat,all.clusters = T, one.sided = T)

# Remove 11,13
toremove <- c(11,13)
sn.da.ctrl <- subsetLiger(sn.da.ctrl, clusters.use = setdiff(levels(sn.da.ctrl@clusters),toremove))
sn.da.ctrl = optimizeALS(sn.da.ctrl, k = 20,lambda = 10)
sn.da.ctrl <- quantile_norm(sn.da.ctrl)
sn.da.ctrl <- louvainCluster(sn.da.ctrl,resolution = 0.8)
sn.da.ctrl = runUMAP(sn.da.ctrl)

levels(sn.da.ctrl@clusters) <- c('SOX6_PART1','SOX6_AGTR1','CALB1_CALCR','SOX6_DDT','SOX6_AGTR1',
                                 'SOX6_AGTR1','SOX6_AGTR1','CALB1_CRYM','CALB1_CD44','SOX6_GFRA2',
                                 'SOX6_AGTR1','CALB1_GEM','CALB1_TRHR','SOX6_AGTR1','SOX6_AGTR1',
                                 'SOX6_AGTR1')

sn.da.ctrl = runUMAP(sn.da.ctrl)
qsave(sn.da.ctrl,'sn_da_ctrl_021821.qs')


sn.da <- qread('/home/tkamath/DA/nurrintegration/sn_da.qs')
pd.meta <- read.csv('../PD_snrna_metadata.csv')
pd.meta <- pd.meta %>% mutate(disease = ifelse(Status %in% c('C','CA'),'Ctrl','Disease'))
pd.meta.use <-pd.meta[match(unique(pd.meta$Donor.ID),pd.meta$Donor.ID),]

sn.da@cell.data$status <- pd.meta.use[match(sn.da@cell.data$dataset,pd.meta.use$Donor.ID),]$disease

dataset.list <- unique(sn.da@cell.data[sn.da@cell.data$status == 'Disease',]$dataset)
for (x in c(1:length(unique(sn.da@cell.data[sn.da@cell.data$status == 'Disease',]$dataset)))  ){
  idx.use <- rownames(sn.da@cell.data[which(sn.da@cell.data$dataset == dataset.list[x]),])
  sn.individual <- subsetLiger(sn.da,cells.use = idx.use,remove.missing = F)
  sn.individual <- addcelldata(sn.individual)
  sn.individual <- rliger::normalize(sn.individual)
  sn.individual@var.genes <- sn.da.ctrl@var.genes
  genes.add <- setdiff(sn.da.ctrl@var.genes,rownames(sn.individual@raw.data[[1]]))
  matrix.add <- Matrix(matrix(0,nrow = length(genes.add),ncol = ncol(sn.individual@norm.data[[1]])),sparse = T)
  rownames(matrix.add) <- genes.add
  colnames(matrix.add) <- rownames(sn.individual@cell.data)
  sn.individual@norm.data <- list(rbind(sn.individual@norm.data[[1]],matrix.add))
  names(sn.individual@norm.data) <- as.character(x)
  sn.individual <- scaleNotCenter(sn.individual)
  sn.individual <- AddMito(sn.individual,species = 'human')
  lis.use <- list(sn.individual)
  names(lis.use) <- dataset.list[x]
  sn.da.ctrl = online_iNMF(sn.da.ctrl, X_new = lis.use,project = TRUE,miniBatch_size = 8000)
}

sn.da.use <- quantile_norm(sn.da.ctrl)
sn.da.use <- louvainCluster(sn.da.use, resolution = 0.4)
sn.da.use = runUMAP(sn.da.use)
plotByDatasetAndCluster(sn.da.use, axis.labels = c("UMAP1","UMAP2"))

sn.da.seurat <- ligerToSeurat(sn.da.use)
markers.da <- prestowrapper(sn.da.seurat,all.clusters = T,one.sided = T)

# Remove 11,12
sn.da.use <- subsetLiger(sn.da.use,clusters.use = setdiff(levels(sn.da.use@clusters),c(11,12)))
sn.da.use <- quantile_norm(sn.da.use)
sn.da.use <- louvainCluster(sn.da.use, resolution = 0.4)
sn.da.use = runUMAP(sn.da.use)
plotByDatasetAndCluster(sn.da.use, axis.labels = c("UMAP1","UMAP2"))
sn.da.use <- rliger::normalize(sn.da.use)

sn.da.use <- subsetLiger(sn.da.use,cells.use = rownames(sn.da.use@cell.data[which(
  sn.da.use@cell.data$percent.mito < 0.1),]))
sn.da.use <- louvainCluster(sn.da.use, resolution = 0.6)
sn.da.use = runUMAP(sn.da.use,rand.seed = 55)
plotByDatasetAndCluster(sn.da.annot, axis.labels = c("UMAP1","UMAP2"))

sn.da.seurat <- ligerToSeurat(sn.da.use)
markers.da <- prestowrapper(sn.da.seurat,all.clusters = T,one.sided = T)
View(markers.da[markers.da$group == 5,])

pdf('da_factorplot_final.pdf')
plotFactors_new(sn.da.use,plot.tsne = T)
dev.off()

pdf('da_tsne_final.pdf')
plotByDatasetAndCluster(sn.da.annot, axis.labels = c("UMAP1","UMAP2"))
rliger::plotGene(sn.da.annot,gene = 'CALB1',plot.by = 'none')
rliger::plotGene(sn.da.annot,gene = 'SOX6',plot.by = 'none')
dev.off()


sn.da.annot <- subsetLiger(sn.da.use,clusters.use = setdiff(levels(sn.da.use@clusters),c(14,15)))
sn.da.annot@cell.data$unannotatedclusters <- sn.da.annot@clusters
sn.da.seurat <- ligerToSeurat(sn.da.annot)
sn.da.annot <- subsetLiger(sn.da.annot,clusters.use = setdiff(levels(sn.da.annot@clusters),c(12,13)))
markers.da <- prestowrapper(sn.da.seurat,all.clusters = T,one.sided = T)
levels(sn.da.annot@clusters) <- c('SOX6_AGTR1','SOX6_PART1','SOX6_AGTR1','CALB1_CALCR','SOX6_DDT',
                                  'SOX6_GFRA2','SOX6_AGTR1','CALB1_PPP1R17','CALB1_CRYM_CCDC68',
                                  'CALB1_GEM','CALB1_RBP4','CALB1_TRHR')


addcelldata <- function(object){
  nUMI <- unlist(lapply(object@raw.data, function(x) {
    Matrix::colSums(x)
  }), use.names = F)
  nGene <- unlist(lapply(object@raw.data, function(x) {
    Matrix::colSums(x > 0)
  }), use.names = F)
  dataset <- unlist(lapply(seq_along(object@raw.data), function(i) {
    rep(names(object@raw.data)[i], ncol(object@raw.data[[i]]))
  }), use.names = F)
  object@cell.data <- data.frame(nUMI, nGene, dataset)
  rownames(object@cell.data) <- unlist(lapply(object@raw.data, function(x) {
    colnames(x)
  }), use.names = F)
  return(object)
}

########################
### Non-DA neurons #####
########################
#non da annotations
setwd('/home/tkamath/DA/nonda/')
options(mc.cores = 60)

data.use <- qread('/home/tkamath/DA/nonda/sn_nonda.qs')
pd.meta <- read.csv('/home/tkamath/DA/PD_snrna_metadata.csv')
pd.meta <- pd.meta %>% mutate(disease = ifelse(Status %in% c('C','CA'),'Ctrl','Disease'))
pd.meta.use <-pd.meta[match(unique(pd.meta$Donor.ID),pd.meta$Donor.ID),]
pd.meta.use$sex <- id.sex$sex

data.use <- rliger::convertOldLiger(data.use)
data.use@cell.data$dx <- pd.meta.use[match(data.use@cell.data$dataset,pd.meta.use$Donor.ID),]$disease

# First annotate the control datasets
idx.use <- rownames(data.use@cell.data[which(data.use@cell.data$dx == 'Ctrl'),])
sn.nonda.liger <- rliger::subsetLiger(data.use,cells.use = idx.use)
sn.nonda.liger@cell.data <- data.use@cell.data[idx.use,]

sn.nonda.liger <- AddMito(sn.nonda.liger,species = 'human')
idx.use <- rownames(sn.nonda.liger@cell.data[which(sn.nonda.liger@cell.data$percent.mito < 0.1),])
celldata <- sn.nonda.liger@cell.data
sn.nonda.liger <- rliger::subsetLiger(data.use,cells.use = idx.use)
sn.nonda.liger@cell.data <- celldata
sn.nonda.liger = rliger::normalize(sn.nonda.liger)

tmp1 <- rliger::createLiger(raw.data = list('total' = MergeSparseDataAll(sn.nonda.liger@raw.data)))
tmp1 <- rliger::normalize(tmp1)
tmp1 <- rliger::selectGenes(tmp1, do.plot = T,num.genes = 3000)
genes.use <- lapply(sn.nonda.liger@raw.data, function(x){rownames(x)})
genes.use <- Reduce(intersect,genes.use)
var.genes.use <- intersect(tmp1@var.genes, genes.use)
sn.nonda.liger@var.genes <- var.genes.use

sn.nonda.liger = rliger::scaleNotCenter(sn.nonda.liger)
sn.nonda.liger = rliger::optimizeALS(sn.nonda.liger, k = 25,lambda = 10)
sn.nonda.liger <- rliger::quantile_norm(sn.nonda.liger)
sn.nonda.liger = rliger::runUMAP(sn.nonda.liger)

# Remove 1,7,9,10,11,13,19,23,24
sn.nonda.liger2 <- rliger::subsetLiger(sn.nonda.liger,clusters.use = setdiff(levels(sn.nonda.liger@clusters),
                                                                             c(1,7,9,10,11,13,19,23,24)))
sn.nonda.liger2 <- rliger::quantile_norm(sn.nonda.liger2)
sn.nonda.liger2 <- rliger::louvainCluster(sn.nonda.liger2, resolution = 0.5)
sn.nonda.liger2 = rliger::runUMAP(sn.nonda.liger2)

# Part of 11, MBP
# Cluster 10 - mainly composed of two individuals
#17 - remove
# 18 - remove
# 16 - remove

sn.nonda.liger2 <-  rliger::subsetLiger(sn.nonda.liger2,clusters.use = setdiff(levels(sn.nonda.liger2@clusters),
                                                                               c(16,17,18)))
qsave(sn.nonda.liger2,'/home/tkamath/DA/nonda/nonda_liger_ctrl_0229.qs')

dataset.list <- unique(data.use@cell.data[data.use@cell.data$dx == 'Disease',]$dataset)
for (x in c(1:length(unique(data.use@cell.data[data.use@cell.data$status == 'Disease',]$dataset)))  ){
  idx.use <- rownames(data.use@cell.data[which(data.use@cell.data$dataset == dataset.list[x]),])
  sn.individual <- rliger::subsetLiger(data.use,cells.use = idx.use,remove.missing = F)
  sn.individual@cell.data <- data.use@cell.data[idx.use,]
  sn.individual <- rliger::normalize(sn.individual)
  sn.individual@var.genes <- sn.nonda.liger2@var.genes
  genes.add <- setdiff(sn.nonda.liger2@var.genes,rownames(sn.individual@raw.data[[1]]))
  matrix.add <- Matrix(matrix(0,nrow = length(genes.add),ncol = ncol(sn.individual@norm.data[[1]])),sparse = T)
  if (length(genes.add) > 0){
    rownames(matrix.add) <- genes.add
    colnames(matrix.add) <- rownames(sn.individual@cell.data)
    sn.individual@norm.data <- list(rbind(sn.individual@norm.data[[1]],matrix.add))
  }
  names(sn.individual@norm.data) <- as.character(x)
  sn.individual <- rliger::scaleNotCenter(sn.individual)
  sn.individual <- AddMito(sn.individual,species = 'human')
  idx.use2 <- rownames(sn.individual@cell.data[which(sn.individual@cell.data$percent.mito < 0.1),])
  sn.individual2 <- rliger::subsetLiger(sn.individual,cells.use = idx.use2)
  sn.individual2@cell.data <- sn.individual@cell.data[idx.use2,]
  lis.use <- list(sn.individual2)
  names(lis.use) <- dataset.list[x]
  sn.nonda.liger2 = rliger::online_iNMF(sn.nonda.liger2, X_new = lis.use,project = TRUE)
  sn.nonda.liger2 <- rliger::normalize(sn.nonda.liger2)
}
sn.nonda.liger2 <- rliger::quantile_norm(sn.nonda.liger2)
sn.nonda.liger2<- rliger::louvainCluster(sn.nonda.liger2)
sn.nonda.liger2 = rliger::runUMAP(sn.nonda.liger2)

qsave(sn.nonda.liger2,'nonda_liger_0229.qs')

sn.nonda.liger2 <- rliger::subsetLiger(sn.nonda.liger2,clusters.use = 
                                         setdiff(levels(sn.nonda.liger2@clusters),c(21,25,29,30,31)))
sn.nonda.liger2<- rliger::louvainCluster(sn.nonda.liger2)
sn.nonda.liger2 = rliger::runUMAP(sn.nonda.liger2)
sn.nonda.annot <- sn.nonda.liger2
sn.nonda.annot <- rliger::subsetLiger(sn.nonda.annot,
                                      clusters.use = setdiff(levels(sn.nonda.annot@clusters),c(4,30,26)))
levels(sn.nonda.annot@clusters) <- c('Inh_PRLR_RP11-384J4.2','Ex_EBF2_CTC-552D5.1','Inh_PRLR_RP11-384J4.2',
                                     'Ex_MYO5B','Inh_SIX3','Inh_IGFBP5','Ex_LAMP5_BAIAP3','Ex_VWA5B1_CALB1',
                                     'Ex_OPRD1','Ex_MYO5B','Ex_POSTN','Ex_LAMP5_NTNG2','Ex_SATB2',
                                     'Ex_CYP2J2','Ex_PPP1R1C','Inh_PAX5_CCBE1','Inh_PAX5_CCBE1',
                                     'Ex_LAMP5_BAIAP3','Inh_OTX2_CASR','Ex_LAMP5_NTNG2','Ex_LAMP5_NTNG2',
                                     'Inh_PAX5_VCAN','Inh_INHBA','Ex_VWA5B1_CALB1','Ex_LAMP5_NTNG2','Ex_PPP1R1C',
                                     'Inh_OTX2_CASR','Ex_MYO5B')

###########################################################
### Astrocytes: note these code were run in R-4.0 #########
## NB: liger was used back when it was called liger #######
###########################################################
rm(list =ls())
library(googleCloudStorageR)
library(liger)
library(qs)
source("source_code.R")
gcs_source("vgazesta/code/SN_alignment.R")

rm(list=ls())
runIndx=26 #26 #29, 27,28
sensitiveSearch=1
input_highly_var_genes=NULL
exNonMicCells=F
ncores=6
newRun=F
includeHuman=F
includeMouse=F
FinnishSbjBased=F
DE_supportingFractionThr=0.1
DE_n.adaptiveKernel=20
DE_nPropIter=3
uniformZscore=F
dist_zscore_gamma=F
dist_zscore_nbinom=F
regularize=F
geoMean=F
pseudocell_count=200
pval_thr=0.001
exclude_non_freq_pseudocells=F
saveDir="~/data/data/tmpBucket/results/SN_astro_indScaling"

#external_DE_path="~/mouse_concensus_HVGgenes.rda"
#external_DE_name="mouseDE"

#external_DE_path="~/hm_liger_var_genes_arranged.rda"
#external_DE_name="ligerVar"

external_DE_path=NULL
external_DE_name=NULL


dist_zscore_gamma=F
dist_zscore_norm=T
dist_zscore_nbinom=F

#prefix="MGonly_"
#prefix="PCprojection2"
#prefix="organismConserved3"
prefix="SN_astro_indScaling"
.ArgList=.myArgFn(runIndx=runIndx,exNonMicCells=F,ncores=ncores,sensitiveSearch=sensitiveSearch,includeHuman=includeHuman,includeMouse=includeMouse,FinnishSbjBased=FinnishSbjBased,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,newRun = newRun,inputDf=NULL,pseudocell_count=pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name,saveDir=saveDir)
.ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
.ArgList$input_highly_var_genes=input_highly_var_genes
.ArgList$saveDirGlobal="~/data/data/tmpBucket/results/sn_astro_indScaling_Global/"
if(!dir.exists(.ArgList$saveDirGlobal)){
  dir.create(.ArgList$saveDirGlobal,recursive = T)
}

.ArgList$HVG_count=6

.ArgList$covariates=NULL
.ArgList$includeHuman=F
.ArgList$includeMouse=F
.ArgList$min_ds_size=50

rm(list=ls())

library(qs)
data=qread(paste0("~/data/data/SN2020/","sn_astro","_arranged.qs"))

reArranged=list()
for(i in 1:length(data$data)){
  tmp=.mySplitObject(data$data[[i]],"batch_merging")
  reArranged=c(reArranged,tmp)
}
data_m=.mycBindFn(inputList=reArranged,batchNames=names(reArranged),verbose=F)

data=list(data=reArranged,data_m=data_m)

for(i in 1:length(data$data)){
  data$data[[i]]$batch_merging=paste0(data$data[[i]]$subject,"_",data$data[[i]]$batch_merging)
}
data$data_m$batch_merging=paste0(data$data_m$subject,"_",data$data_m$batch_merging)
.data=data


#dataorganism="Human";argList = .ArgList
data=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
.data2=data
length(data$varFeatures)

expData=data$data_m@assays$RNA@data
expData=expData[row.names(expData) %in% data$varFeatures,]
expData=as.matrix(expData)
nUMIdata=data$data_m$QC_Gene_unique_count

cordata=NULL
for(ik in 1:nrow(expData)){
  cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
}

print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
#data$varFeatures=as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])

tmp=data[c("varFeatures","allGenes" )]
save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))


tmp=data$data_m
save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))


cat('Scaling the data and PCA\n')
.myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)

#save(.ArgList,file="~/ArgList_IndScaling.rda")
save(.ArgList,file=paste0("~/data/data/SN2020/ArgList_","sn_astro_indScaling_HVG",.ArgList$HVG_count,".rda"))

#load("~/ArgList_IndScaling.rda")
load(paste0("~/data/data/SN2020/ArgList_","sn_astro_indScaling_HVG",.ArgList$HVG_count,".rda"))

.mySecondStepFn()
.myThirdStepFn()


###########################################################
### MG: note these code were run in R-4.0 #################
## NB: liger was used back when it was called liger #######
###########################################################

library(googleCloudStorageR)
library(liger)
library(qs)
source("source_code.R")
gcs_source("vgazesta/code/SN_alignment.R")


rm(list=ls())

runIndx=26 #26 #29, 27,28
sensitiveSearch=1
input_highly_var_genes=NULL
exNonMicCells=F
ncores=6
newRun=F
includeHuman=F
includeMouse=F
FinnishSbjBased=F
DE_supportingFractionThr=0.1
DE_n.adaptiveKernel=20
DE_nPropIter=3
uniformZscore=F
dist_zscore_gamma=F
dist_zscore_nbinom=F
regularize=F
geoMean=F
pseudocell_count=200
pval_thr=0.001
exclude_non_freq_pseudocells=F
saveDir="~/data/data/tmpBucket/results/SN_mg_indScaling"

#external_DE_path="~/mouse_concensus_HVGgenes.rda"
#external_DE_name="mouseDE"

#external_DE_path="~/hm_liger_var_genes_arranged.rda"
#external_DE_name="ligerVar"

external_DE_path=NULL
external_DE_name=NULL


dist_zscore_gamma=F
dist_zscore_norm=T
dist_zscore_nbinom=F


#prefix="MGonly_"
#prefix="PCprojection2"
#prefix="organismConserved3"
prefix="SN_mg_indScaling"
.ArgList=.myArgFn(runIndx=runIndx,exNonMicCells=F,ncores=ncores,sensitiveSearch=sensitiveSearch,includeHuman=includeHuman,includeMouse=includeMouse,FinnishSbjBased=FinnishSbjBased,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,newRun = newRun,inputDf=NULL,pseudocell_count=pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name,saveDir=saveDir)
.ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
.ArgList$input_highly_var_genes=input_highly_var_genes
.ArgList$saveDirGlobal="~/data/data/tmpBucket/results/mg_indScaling_Global"
if(!dir.exists(.ArgList$saveDirGlobal)){
  dir.create(.ArgList$saveDirGlobal,recursive = T)
}
.ArgList$HVG_count=6

.ArgList$covariates=NULL
.ArgList$includeHuman=F
.ArgList$includeMouse=F
.ArgList$min_ds_size=50

rm(list=ls())
data=qread(paste0("~/data/data/SN2020/","sn_mg","_arranged.qs"))

reArranged=list()
for(i in 1:length(data$data)){
  tmp=.mySplitObject(data$data[[i]],"batch_merging")
  reArranged=c(reArranged,tmp)
}
data_m=.mycBindFn(inputList=reArranged,batchNames=names(reArranged),verbose=F)

data=list(data=reArranged,data_m=data_m)

for(i in 1:length(data$data)){
  data$data[[i]]$batch_merging=paste0(data$data[[i]]$subject,"_",data$data[[i]]$batch_merging)
}
data$data_m$batch_merging=paste0(data$data_m$subject,"_",data$data_m$batch_merging)
.data=data


#dataorganism="Human";argList = .ArgList
data=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
.data2=data
length(data$varFeatures)

expData=data$data_m@assays$RNA@data
expData=expData[row.names(expData) %in% data$varFeatures,]
expData=as.matrix(expData)
nUMIdata=data$data_m$QC_Gene_unique_count

cordata=NULL
for(ik in 1:nrow(expData)){
  cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
}

print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
#data$varFeatures=as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])

tmp=data[c("varFeatures","allGenes" )]
save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))


tmp=data$data_m
save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))


cat('Scaling the data and PCA\n')
.myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)

#save(.ArgList,file="~/ArgList_IndScaling.rda")
save(.ArgList,file=paste0("~/data/data/SN2020/ArgList_","sn_mg_indScaling_HVG",.ArgList$HVG_count,".rda"))


.mySecondStepFn()
.myThirdStepFn()


###########################################################
### Endo: note these code were run in R-4.0 ###############
## NB: liger was used back when it was called liger #######
###########################################################
#SN Endo

library(googleCloudStorageR)
library(liger)
library(qs)
source("source_code.R")
gcs_source("vgazesta/code/SN_alignment.R")


rm(list=ls())
runIndx=26 #26 #29, 27,28
sensitiveSearch=1
input_highly_var_genes=NULL
exNonMicCells=F
ncores=6
newRun=F
includeHuman=F
includeMouse=F
FinnishSbjBased=F
DE_supportingFractionThr=0.1
DE_n.adaptiveKernel=20
DE_nPropIter=3
uniformZscore=F
dist_zscore_gamma=F
dist_zscore_nbinom=F
regularize=F
geoMean=F
pseudocell_count=200
pval_thr=0.001
exclude_non_freq_pseudocells=F
saveDir="~/data/data/tmpBucket/results/SN_endo_indScaling"

#external_DE_path="~/mouse_concensus_HVGgenes.rda"
#external_DE_name="mouseDE"

#external_DE_path="~/hm_liger_var_genes_arranged.rda"
#external_DE_name="ligerVar"

external_DE_path=NULL
external_DE_name=NULL


dist_zscore_gamma=F
dist_zscore_norm=T
dist_zscore_nbinom=F


#prefix="MGonly_"
#prefix="PCprojection2"
#prefix="organismConserved3"
prefix="SN_endo_indScaling"
.ArgList=.myArgFn(runIndx=runIndx,exNonMicCells=F,ncores=ncores,sensitiveSearch=sensitiveSearch,includeHuman=includeHuman,includeMouse=includeMouse,FinnishSbjBased=FinnishSbjBased,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,newRun = newRun,inputDf=NULL,pseudocell_count=pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name,saveDir=saveDir)
.ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
.ArgList$input_highly_var_genes=input_highly_var_genes
.ArgList$saveDirGlobal="~/data/data/tmpBucket/results/endo_indScaling_Global"
if(!dir.exists(.ArgList$saveDirGlobal)){
  dir.create(.ArgList$saveDirGlobal,recursive = T)
}
.ArgList$HVG_count=6

.ArgList$covariates=NULL
.ArgList$includeHuman=F
.ArgList$includeMouse=F
.ArgList$min_ds_size=50

rm(list=ls())
data=qread(paste0("~/data/data/SN2020/","sn_endo","_arranged.qs"))

reArranged=list()
for(i in 1:length(data$data)){
  tmp=.mySplitObject(data$data[[i]],"batch_merging")
  reArranged=c(reArranged,tmp)
}
data_m=.mycBindFn(inputList=reArranged,batchNames=names(reArranged),verbose=F)

data=list(data=reArranged,data_m=data_m)

for(i in 1:length(data$data)){
  data$data[[i]]$batch_merging=paste0(data$data[[i]]$subject,"_",data$data[[i]]$batch_merging)
}
data$data_m$batch_merging=paste0(data$data_m$subject,"_",data$data_m$batch_merging)
.data=data


#dataorganism="Human";argList = .ArgList
data=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
.data2=data
length(data$varFeatures)

expData=data$data_m@assays$RNA@data
expData=expData[row.names(expData) %in% data$varFeatures,]
expData=as.matrix(expData)
nUMIdata=data$data_m$QC_Gene_unique_count

cordata=NULL
for(ik in 1:nrow(expData)){
  cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
}

print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
#data$varFeatures=as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])

tmp=data[c("varFeatures","allGenes" )]
save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))


tmp=data$data_m
save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))


cat('Scaling the data and PCA\n')
.myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)

#save(.ArgList,file="~/ArgList_IndScaling.rda")
save(.ArgList,file=paste0("~/data/data/SN2020/ArgList_","sn_endo_indScaling_HVG",.ArgList$HVG_count,".rda"))


.mySecondStepFn()
.myThirdStepFn()

###########################################################
### Oligo: note these code were run in R-4.0 ##############
## NB: liger was used back when it was called liger #######
###########################################################
#SN Oligo

library(googleCloudStorageR)
library(liger)
library(qs)
source("source_code.R")
gcs_source("vgazesta/code/SN_alignment.R")


rm(list=ls())
runIndx=26 #26 #29, 27,28
sensitiveSearch=1
input_highly_var_genes=NULL
exNonMicCells=F
ncores=6
newRun=F
includeHuman=F
includeMouse=F
FinnishSbjBased=F
DE_supportingFractionThr=0.1
DE_n.adaptiveKernel=20
DE_nPropIter=3
uniformZscore=F
dist_zscore_gamma=F
dist_zscore_nbinom=F
regularize=F
geoMean=F
pseudocell_count=200
pval_thr=0.001
exclude_non_freq_pseudocells=F
saveDir="~/data/data/tmpBucket/results/SN_olig_indScaling"

#external_DE_path="~/mouse_concensus_HVGgenes.rda"
#external_DE_name="mouseDE"

#external_DE_path="~/hm_liger_var_genes_arranged.rda"
#external_DE_name="ligerVar"

external_DE_path=NULL
external_DE_name=NULL


dist_zscore_gamma=F
dist_zscore_norm=T
dist_zscore_nbinom=F


#prefix="MGonly_"
#prefix="PCprojection2"
#prefix="organismConserved3"
prefix="SN_olig_indScaling"
.ArgList=.myArgFn(runIndx=runIndx,exNonMicCells=F,ncores=ncores,sensitiveSearch=sensitiveSearch,includeHuman=includeHuman,includeMouse=includeMouse,FinnishSbjBased=FinnishSbjBased,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,newRun = newRun,inputDf=NULL,pseudocell_count=pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name,saveDir=saveDir)
.ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
.ArgList$input_highly_var_genes=input_highly_var_genes
.ArgList$saveDirGlobal="~/data/data/tmpBucket/results/olig_indScaling_Global"
if(!dir.exists(.ArgList$saveDirGlobal)){
  dir.create(.ArgList$saveDirGlobal,recursive = T)
}
.ArgList$HVG_count=6

.ArgList$covariates=NULL
.ArgList$includeHuman=F
.ArgList$includeMouse=F
.ArgList$min_ds_size=50

rm(list=ls())
data=qread(paste0("~/data/data/SN2020/","sn_olig","_arranged.qs"))

reArranged=list()
for(i in 1:length(data$data)){
  tmp=.mySplitObject(data$data[[i]],"batch_merging")
  reArranged=c(reArranged,tmp)
}
data_m=.mycBindFn(inputList=reArranged,batchNames=names(reArranged),verbose=F)

data=list(data=reArranged,data_m=data_m)

for(i in 1:length(data$data)){
  data$data[[i]]$batch_merging=paste0(data$data[[i]]$subject,"_",data$data[[i]]$batch_merging)
}
data$data_m$batch_merging=paste0(data$data_m$subject,"_",data$data_m$batch_merging)
.data=data


#dataorganism="Human";argList = .ArgList
data=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
.data2=data
length(data$varFeatures)

expData=data$data_m@assays$RNA@data
expData=expData[row.names(expData) %in% data$varFeatures,]
expData=as.matrix(expData)
nUMIdata=data$data_m$QC_Gene_unique_count

cordata=NULL
for(ik in 1:nrow(expData)){
  cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
}

print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
#data$varFeatures=as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])

tmp=data[c("varFeatures","allGenes" )]
save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))


tmp=data$data_m
save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))


cat('Scaling the data and PCA\n')
.myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)

#save(.ArgList,file="~/ArgList_IndScaling.rda")
save(.ArgList,file=paste0("~/data/data/SN2020/ArgList_","sn_olig_indScaling_HVG",.ArgList$HVG_count,".rda"))


.mySecondStepFn()
.myThirdStepFn()
.myForthStepFn()

###########################################################
### OPC: note these code were run in R-4.0 ################
## NB: liger was used back when it was called liger #######
###########################################################
#SN OPC

library(googleCloudStorageR)
library(liger)
library(qs)
source("source_code.R")
gcs_source("vgazesta/code/SN_alignment.R")


rm(list=ls())
runIndx=26 #26 #29, 27,28
sensitiveSearch=1
input_highly_var_genes=NULL
exNonMicCells=F
ncores=6
newRun=F
includeHuman=F
includeMouse=F
FinnishSbjBased=F
DE_supportingFractionThr=0.1
DE_n.adaptiveKernel=20
DE_nPropIter=3
uniformZscore=F
dist_zscore_gamma=F
dist_zscore_nbinom=F
regularize=F
geoMean=F
pseudocell_count=200
pval_thr=0.001
exclude_non_freq_pseudocells=F
saveDir="~/data/data/tmpBucket/results/SN_opc_indScaling"

#external_DE_path="~/mouse_concensus_HVGgenes.rda"
#external_DE_name="mouseDE"

#external_DE_path="~/hm_liger_var_genes_arranged.rda"
#external_DE_name="ligerVar"

external_DE_path=NULL
external_DE_name=NULL


dist_zscore_gamma=F
dist_zscore_norm=T
dist_zscore_nbinom=F


#prefix="MGonly_"
#prefix="PCprojection2"
#prefix="organismConserved3"
prefix="SN_opc_indScaling"
.ArgList=.myArgFn(runIndx=runIndx,exNonMicCells=F,ncores=ncores,sensitiveSearch=sensitiveSearch,includeHuman=includeHuman,includeMouse=includeMouse,FinnishSbjBased=FinnishSbjBased,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,newRun = newRun,inputDf=NULL,pseudocell_count=pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name,saveDir=saveDir)
.ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
.ArgList$input_highly_var_genes=input_highly_var_genes
.ArgList$saveDirGlobal="~/data/data/tmpBucket/results/opc_indScaling_Global"
if(!dir.exists(.ArgList$saveDirGlobal)){
  dir.create(.ArgList$saveDirGlobal,recursive = T)
}
.ArgList$HVG_count=4
.ArgList$allGenesFraction=0.5

.ArgList$covariates=NULL
.ArgList$includeHuman=F
.ArgList$includeMouse=F
.ArgList$min_ds_size=50

rm(list=ls())
data=qread(paste0("~/data/data/SN2020/","sn_opc","_arranged.qs"))

reArranged=list()
for(i in 1:length(data$data)){
  tmp=.mySplitObject(data$data[[i]],"batch_merging")
  reArranged=c(reArranged,tmp)
}
data_m=.mycBindFn(inputList=reArranged,batchNames=names(reArranged),verbose=F)

data=list(data=reArranged,data_m=data_m)

for(i in 1:length(data$data)){
  data$data[[i]]$batch_merging=paste0(data$data[[i]]$subject,"_",data$data[[i]]$batch_merging)
}
data$data_m$batch_merging=paste0(data$data_m$subject,"_",data$data_m$batch_merging)
.data=data


#dataorganism="Human";argList = .ArgList
data=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
.data2=data
length(data$varFeatures)

expData=data$data_m@assays$RNA@data
expData=expData[row.names(expData) %in% data$varFeatures,]
expData=as.matrix(expData)
nUMIdata=data$data_m$QC_Gene_unique_count

cordata=NULL
for(ik in 1:nrow(expData)){
  cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
}

print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
#data$varFeatures=as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])

tmp=data[c("varFeatures","allGenes" )]
save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))


tmp=data$data_m
save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))


cat('Scaling the data and PCA\n')
.myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)

#save(.ArgList,file="~/ArgList_IndScaling.rda")
save(.ArgList,file=paste0("~/data/data/SN2020/ArgList_","sn_opc_indScaling_HVG",.ArgList$HVG_count,".rda"))

.mySecondStepFn()
.myThirdStepFn()
.myForthStepFn()

