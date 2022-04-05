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

############################################
## Annotation of non-neuronal cell types ###
############################################
library(devtools)
#install_github("MacoskoLab/liger", ref = "online")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source("/home/tkamath/scripts/prestowrapper.R")
#library(Seurat)
library(rliger)
library(Matrix)
library(Seurat)
library("tidyverse")
library(qs)
paths = list.files(path="/home/tkamath/DA/SN2020_selected/",pattern = "*indScaling",full.names = T)

objs=lapply(paths,function(q){
  make.seurat(qread(paste0(q,"/data.qs")))
})
objs = lapply(objs,NormalizeData)
names(objs)= c("astro","endo","mg",'neurons',"olig","opc")

#### Astros
setwd('/home/tkamath/DA/astro/')
levels.use = levels(objs$astro@ident)[which(table(objs$astro@ident)>3)]
markers.astro <- prestowrapper(objs$astro,all.clusters = T,one.sided = T)
pdf('tsne.pdf')
DimPlot(objs$astro, pt.size = 0.5,reduction.use = 'umap',do.label = T,no.legend = T)
dev.off()
markers.astro <- split(markers.astro,f = markers.astro$group)
View(markers.astro$`3`)
FeaturePlot(objs$astro, features = "ZBBX",reduction.use = 'umap',pt.size = 0.5)

m = prestowrapper(objs$astro,ident.1 = 0,ident.2 = 2)
m2 = prestowrapper(objs$astro,ident.1 = 8,ident.2 = 9)
View(m2[m2$group == 'ident.1',])

# 0,2,3 - Astro_GJB6_OXTR
# 1 - Astro_CYP4F12
# 6 - Astro_GUCY1A2
# 7 - Astro_SIDT1
# 8 - Astro_GLYATL2
# 9 - Astro_VIM_TNFSRF12A
# 11 - Astro_SERPINA3
# 13 - Astro_GBP2_SPOCD1
# 14 - Astro_VIM_LHX2
# 15 - Ependyma_ZBBX

# 4 - REMOVE, oligo doublet
# 5- REMOVE, MT genes
# 12 - REMOVE mg doublet
# 10 -REMOVE neuron doublet
# 16 - REMOVE neuron doublet
# 17 - REMOVE OPC doublet

objs$astro <- SubsetData(objs$astro,ident.remove = c(4,5,12,10,16,17))
levels(objs$astro@ident) = c("Astro_GJB6_OXTR",'Astro_CYP4F12','Astro_GJB6_OXTR','Astro_GJB6_OXTR',
                             'Astro_GUCY1A2','Astro_SIDT1','Astro_GLYATL2','Astro_VIM_TNFSRF12A',
                             'Astro_SERPINA3','Astro_GBP2_SPOCD1','Astro_VIM_LHX2','Ependyma_ZBBX')
objs$astro <- RunUMAP(objs$astro,reduction.use = 'hpca',dims.use = c(1:30))

pdf('tsne_annot.pdf')
DimPlot(objs$astro, label = TRUE, pt.size = 0.5,reduction.use = 'umap',
        do.label = T,no.legend = T,label.size = 4)
dev.off()
astro.annotated <- objs$astro

qsave(astro.annotated,file = "/home/tkamath/DA/SN2020_selected/astro_full_annotated.qs")


##### Endos
setwd('/home/tkamath/DA/endofibro/')
levels.use = levels(objs$endo@ident)[which(table(objs$endo@ident)>3)]
markers.endo <- prestowrapper(objs$endo,all.clusters = T,one.sided = T)
View(markers.endo$`8`)

pdf('tsne.pdf')
DimPlot(objs$endo, pt.size = 0.5,reduction.use = 'umap',do.label = T,no.legend = T)
DimPlot(objs$endo, pt.size = 0.5,reduction.use = 'umap',group.by = 'subject')
dev.off()

FeaturePlot(objs$endo, features = "ABCC9",reduction.use = 'umap',pt.size = 0.5)

# Endos:
# 0,3,14 - Endo_DCN_ABCC9
# 1 - Endo_COL6A3
# 2,5,8 - Endo_IL27RA
# 4 - Endo_NOTCH3_PLK2
# 6 - Endo_SNTG2
# 7 - Endo_SLIT3
# 13 - Endo_MET

# 9 - Oligo doublet, REMOVE
# 10 - Neuron doublet, REMOVE
# 11 - MG doublet, REMOVE
# 12 - MG doublet, REMOVE
# 15 - Astro doublet, REMOVE
# 16 - Astro doublet, REMOVE

objs$endo <- SubsetData(objs$endo,ident.remove = c(9,10,11,12,15,16))
endo.annotated <- objs$endo
levels(endo.annotated@ident) <- c('Endo_DCN_ABCC9','Endo_COL6A3','Endo_IL27RA','Endo_DCN_ABCC9',
                                  'Endo_NOTCH3_PLK2','Endo_IL27RA','Endo_SNTG2','Endo_SLIT3',
                                  'Endo_IL27RA','Endo_MET','Endo_DCN_ABCC9')

endo.annotated <- RunUMAP(endo.annotated,reduction.use = 'hpca',dims.use = c(1:30))

pdf('tsne_annot.pdf')
DimPlot(endo.annotated, label = TRUE, pt.size = 0.5,reduction.use = 'umap',
        do.label = T,no.legend = T,label.size = 4)
dev.off()

qsave(endo.annotated,file = "/home/tkamath/DA/SN2020_selected/SN_endo_indScaling/endo_full_annotated.qs")

##### MG
setwd('/home/tkamath/DA/mg/')
levels.use = levels(objs$mg@ident)[which(table(objs$mg@ident)>3)]
markers.mg <- prestowrapper(objs$mg,all.clusters = T,one.sided = T)
markers.mg <- split(markers.mg,f = markers.mg$group)
View(markers.mg$`8`)

pdf('tsne.pdf')
DimPlot(objs$mg, pt.size = 0.5,reduction.use = 'umap',do.label = T,no.legend = T)
DimPlot(objs$mg, pt.size = 0.5,reduction.use = 'umap',group.by = 'dataset')
dev.off()

objs$mg <- AddMito(objs$mg,species = 'human')

FeaturePlot(objs$mg,features.plot = 'percent.mito',reduction.use = 'umap',
            pt.size = 0.5,no.legend = F)
FeaturePlot(mg.annotated,features.plot = 'IL15',reduction.use = 'umap',pt.size = 0.5)

mg.annotated@meta.data <- mg.annotated@meta.data %>% mutate(nurr = ifelse(grepl('DAPI',seq_batch),'DAPI','Nurr'))
DimPlot(mg.annotated, pt.size = 0.5,reduction.use = 'umap',do.label = T,no.legend = T)
DimPlot(mg.annotated, pt.size = 0.5,reduction.use = 'umap',group.by = 'subject')

m2 <- prestowrapper(mg.annotated,ident.1 = 6,ident.2 = c(1),one.sided = T)

m3 <- prestowrapper(mg.annotated,ident.1 = 11,ident.2 = c(20),one.sided = T)

# 0 - MG_CECR2_FGL1
# 1,7 - MG_TSPO_VIM
# 2 - MG_OPRM1
# 5 - MG_GPNMB_SULT1C2
# 6 - MG_GPNMB_SUSD1
# 8 - MG_SPON1... check this
# 9 - Macro_CD200R1
# 10 - MG_LPL
# 11 - MG_FOS
# 14 - RP genes
# 15 - MG_MGAM
# 16 - REMOVE
# 18 - REMOVE
# 19 - MG_CCL3
# 20 - MG_FOSL2
# 21 - MG_MKI67

# 3 - High MT, REEMOVE
# 4 - Oligo, REMOVE
# 12,13 - Neuron, REMOVE
# 17 - Astrocyte, REMOVE
# 22 - T-cell, REMOVE
# 23 - Endo, REMOVE

mg.annotated <- SubsetData(objs$mg,ident.remove = c(3,4,12,13,17,22,23))
mg.annotated <- RunUMAP(mg.annotated,reduction.use = 'hpca',dims.use = c(1:30))

DimPlot(mg.annotated, pt.size = 0.5,reduction.use = 'umap',do.label = T,no.legend = T)

mg.annotated <- SubsetData(mg.annotated,ident.remove = c(14,16,18))

levels(mg.annotated@ident) <- c('MG_CECR2_FGL1','MG_TSPO_VIM','MG_OPRM1','MG_GPNMB_SULT1C2','MG_GPNMB_SUSD1',
                                'MG_TSPO_VIM','MG_SPON1','Macro_CD200R1','MG_GPNMB_LPL','MG_FOSL2',
                                'MG_MGAM','MG_CCL3','MG_FOSL2','MG_MKI67')
mg.annotated <- RunUMAP(mg.annotated,reduction.use = 'hpca',dims.use = c(1:30))

DimPlot(mg.annotated, pt.size = 0.5,reduction.use = 'umap',do.label = T,no.legend = T)

qsave(mg.annotated,'/home/tkamath/DA/SN2020_selected/SN_mg_indScaling/mg_full_annotated.qs')

##### Neurons
sn.nonda<- qread('/home/tkamath/DA/SN2020_selected/solution1/NonDANEuron_AnnotatedSeurat.qs')
head(sn.nonda@meta.data)
ident.use <- sn.nonda@active.ident[which(sn.nonda@active.ident != 'REMOVE')]

sn.nonda.use<- qread('/home/tkamath/DA/SN2020_selected/solution1/data.qs')
sn.nonda.seurat <- CreateSeuratObject(raw.data = sn.nonda.use$countData,project = "SeuratProject",
                                      assay = "RNA",
                                      min.cells = 0,
                                      min.features = 0,
                                      names.field = 1,
                                      names.delim = "-",
                                      meta.data = sn.nonda.use$pdata)
sn.nonda.seurat = SetDimReduction(object=sn.nonda.seurat,reduction.type='umap',
                                  slot='cell.embeddings',new.data = as.matrix(sn.nonda@meta.data[,c('UMAP_1','UMAP_2')]))
sn.nonda.seurat = SetDimReduction(sn.nonda.seurat,reduction.type='umap',slot='key',new.data = 'UMAP_')

hpca.obj<-new(Class="dim.reduction",
              cell.embeddings=sn.nonda.use$embedding,key="HPCA",misc=list("raw" = sn.nonda.use$embedding))
colnames(hpca.obj@cell.embeddings) = paste0("Factor",1:ncol(sn.nonda.use$embedding))
sn.nonda.seurat@dr$hpca <- hpca.obj

sn.nonda.seurat <- SubsetData(sn.nonda.seurat,cells.use = names(ident.use))
sn.nonda.seurat@ident <- ident.use
sn.nonda.seurat <- SubsetData(sn.nonda.seurat,cells.use = names(ident.use),subset.raw = T)
#sn.nonda.seurat <-RunUMAP(sn.nonda.seurat,reduction.use = 'hpca',dims.use = c(1:30))

pdf('/home/tkamath/DA/nonda/tsne_annot.pdf')
DimPlot(sn.nonda.seurat,reduction.use = 'umap',do.label = T,no.legend = T)
dev.off()

qsave(sn.nonda.seurat,'/home/tkamath/DA/SN2020_selected/solution1/snnonda_full_annotated.qs')

#######Olig
setwd('/home/tkamath/DA/oligo/')
levels.use = levels(objs$olig@ident)[which(table(objs$olig@ident)>3)]
markers.olig <- prestowrapper(objs$olig,all.clusters = T,one.sided = T)

m1 <- prestowrapper(objs$olig,ident.1 = c(0,3,4),one.sided = T)
View(markers.olig[markers.olig$group == '14',])

pdf('tsne.pdf')
DimPlot(objs$olig, do.label = TRUE, pt.size = 0.5,reduction.use = 'umap')
dev.off()

FeaturePlot(objs$olig,features.plot = 'KCNAB1',reduction.use = 'umap',pt.size = 0.5)

objs$olig <- AddMito(objs$olig,species = 'human')
FeaturePlot(objs$olig,features.plot = 'PLXDC2',reduction.use = 'umap',pt.size = 0.5,no.legend = F)

objs$olig <- SubsetData(objs$olig,ident.remove =  c(10,11))
objs$olig <- RunUMAP(objs$olig,reduction.use = 'hpca',dims.use = c(1:30))

olig.annotated <- SubsetData(objs$olig,ident.remove =  c(14))
levels(olig.annotated@ident) <- c('Olig_PLXDC2','Olig_PLXDC2_SFRP1','Olig_PLXDC2_KCNAB1','Olig_PLXDC2_KCNK10','Olig_ENPP6_LUCAT1',
                                  'Olig_PLXDC2_KCNAB1','Olig_PLXDC2_SFRP1','Olig_ENPP6_LUCAT1','Olig_PLXDC2_SFRP1','Olig_ENPP6_EMILIN2',
                                  'Olig_ENPP6_ACTN2','Olig_ENPP6_EMILIN2')

olig.annotated <- RunUMAP(olig.annotated,reduction.use = 'hpca',dims.use = c(1:30))
pdf('tsne.pdf')
DimPlot(olig.annotated, pt.size = 0.5,reduction.use = 'umap',do.label = T,no.legend = T)
dev.off()

#OLIGOS
# 0 - PLXDC2
# 1,4,6 - PLXDC2_SFRP1
# 3,12 - PLXDC2_KCNAB1
# 13 - PLXDC2_KCNK10

# 2,5 - ENPP6_LUCAT1
# 8 - ENPP6_ACTN2
# 7,9 - ENPP6_EMILIN2

# 10 - astro doublet
# 11 - astro doublet
# 14 - REMOVE

qsave(olig.annotated,file = "/home/tkamath/DA/SN2020_selected/SN_olig_indScaling/olig_full_annotated.qs")


### OPC
levels.use = levels(objs$opc@ident)[which(table(objs$opc@ident)>3)]
markers.opc <- prestowrapper(objs$opc,all.clusters = T,one.sided = T)
View(markers.opc[markers.opc$group == 1,])

DimPlot(objs$opc, do.label = TRUE, pt.size = 0.5,reduction.use = 'umap')
FeaturePlot(opc.annotated, features ="CACNG4",reduction.use = 'umap')

# 0 -all cacng4
# 1 - 
# 2 - 
# 3- 
# 4 - 
# 5 - 
# 8 - OPC_HOXD3
# 9 - OPC_KIAA0040
# 10 - OPC_ADM
# 16 - OPC_MDFI

# 6 - REMOVE, mature oligo
# 7 - REMOVE
# 11 - REEMOVE
# 12 - REMOVE
# 13 - Neeuron doublet, remove
# 14 - MG doublet, remove
# 15- astro doublet remove
# 17 - REMOVE


opc.annotated <- SubsetData(objs$opc,ident.remove = c(6,7,11,12,13,14,15,17))
opc.annotated <- RunUMAP(opc.annotated,reduction.use = 'hpca',dims.use = c(1:30))
DimPlot(opc.annotated, do.label = TRUE, pt.size = 0.5,reduction.use = 'umap')

levels(opc.annotated@ident) <- c('OPC_CACNG4','OPC_CACNG4','OPC_CACNG4','OPC_CACNG4','OPC_CACNG4',
                                 'OPC_CACNG4','OPC_HOXD3','OPC_KIAA0040','OPC_ADM','OPC_MDFI')
qsave(opc.annotated,'/home/tkamath/DA/SN2020_selected/SN_opc_indScaling/opc_full_annotated.qs')


make.seurat<-function(data){
  nbt=Seurat::CreateSeuratObject(raw.data=data$countData,
                                 project = "SeuratProject",
                                 min.cells = 0,
                                 min.features = 0,
                                 names.field = 1,
                                 names.delim = "-",
                                 meta.data = data$pdata)
  nbt@ident <- data$pdata$anno_cluster_res
  names(nbt@ident) <- rownames(data$pdata)
  nbt = SetDimReduction(object=nbt,reduction.type='umap',
                        slot='cell.embeddings',new.data = as.matrix(data$pdata[,c('UMAP_1','UMAP_2')]))
  nbt = SetDimReduction(nbt,reduction.type='umap',slot='key',new.data = 'UMAP_')
  
  hpca.obj<-new(Class="dim.reduction",
                cell.embeddings=data$embedding,key="HPCA",misc=list("raw" = data$embedding))
  colnames(hpca.obj@cell.embeddings) = paste0("Factor",1:ncol(data$embedding))
  nbt@dr$hpca <- hpca.obj
  
  return(nbt)
}

