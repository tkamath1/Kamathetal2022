## Creating Nurr+ and Nurr- integrations
# nurr integration
library(Seurat)
library(rliger)
library(qs)
source("/home/tkamath/DA/PDpaper_reproducibility/genevalidatehelper.R")
source("/home/tkamath/DA/PDpaper_reproducibility/SeuratExtrafunctions.R")
source("/home/tkamath/DA/PDpaper_reproducibility/extrafuncs.R")
source('/home/tkamath/DA/PDpaper_reproducibility/prestowrapper.R')

#Generate cleaned objects and metadata
pd.meta.use <- qread('/home/tkamath/DA/PDpaper_reproducibility/pdmeta_use.qs')

# First combine all major cell type objcects
# Code to generate these objects can be found at: &&&& and objects can be downloaded from: ****
list.out <- lapply(grep('annotated.qs',list.files('.',recursive = T),value = T),function(x){qread(x)})
names(list.out) <- c('astro','endo','mg','olig', 'opc','nonda')

# Note: code to generate these objects can be found at: 'nonda_annotations.R' and 'daanalysis_final.R'
# Objects can be downloaded from: ****
list.out$nonda <- qread('/home/tkamath/DA/nonda/nonda_liger_annot_0229.qs')
list.out$da <- qread('/home/tkamath/DA/da/sn_da_annot_021821.qs')
list.out$da <- rliger::ligerToSeurat(list.out$da)
list.out$nonda <- rliger::ligerToSeurat(list.out$nonda)

#qsave(list.out,'/home/tkamath/DA/cleanedobjects.qs')

# Now grab list metadata
list.metadata <- lapply(list.out,function(x){
  x@meta.data$clusters <- x@ident
  return(x@meta.data)
})

list.metadata.cleaned <- lapply(list.metadata[c(1:5)],function(x){
  x <- x[,c('nGene','nUMI','subject','seq_batch','clusters')]
  x$sex <- as.factor(pd.meta.use[match(x$subject,pd.meta.use$Donor.ID),]$sex)
  x$status <- as.factor(pd.meta.use[match(x$subject,pd.meta.use$Donor.ID),]$Status)
  x$disease <- as.factor(pd.meta.use[match(x$subject,pd.meta.use$Donor.ID),]$disease)
  x$lib <- as.factor(paste0(x$subject,x$seq_batch))
  x$seq_batch <- NULL
  return(x)
})

list.add <- lapply(list.metadata[c(6,7)],function(x){
  colnames(x) <- c('nGene','nUMI','subject','clusters')
  x$clusters <- as.factor(x$clusters)
  x$sex <- as.factor(pd.meta.use[match(x$subject,pd.meta.use$Donor.ID),]$sex)
  x$status <- as.factor(pd.meta.use[match(x$subject,pd.meta.use$Donor.ID),]$Status)
  x$disease <- as.factor(pd.meta.use[match(x$subject,pd.meta.use$Donor.ID),]$disease)
  x$lib <- as.factor(unlist(sapply(strsplit(rownames(x),'_'),`[`,2)))
  return(x)
})
list.metadata.cleaned$nonda <- list.add$nonda
list.metadata.cleaned$da <- list.add$da

#qsave(list.metadata.cleaned,'/home/tkamath/DA/cleanedmetadata.qs')
all.obj <- qread('/home/tkamath/DA/cleanedobjects.qs')

# Note: code to generate cleaned object metadata is available at: 'makecleanedobjects.R'
cleaned.meta <- qread('/home/tkamath/DA/cleanedmetadata.qs')
cleaned.meta <- bind_rows(cleaned.meta,.id = 'names')
cleaned.meta.use <- cleaned.meta[which(cleaned.meta$disease == 'Ctrl' & 
                                         grepl('Nurr',cleaned.meta$lib)),]
idx.use <- rownames(cleaned.meta.use)

out.use <- MergeSparseDataAll(lapply(all.obj,function(x){
  x@raw.data[,intersect(colnames(x@raw.data),idx.use)]
}))
out.use <- lapply(unique(cleaned.meta.use$subject),function(x){
  out.use[,which(grepl(x,colnames(out.use)))]
})
names(out.use) <- unique(cleaned.meta.use$subject)

# Make nurr and non-nurr objects
lig.nurr <- createLiger(raw.data = out.use)
lig.nurr = rliger::normalize(lig.nurr)
lig.nurr = selectGenes(lig.nurr, num.genes = 1000)
lig.nurr = scaleNotCenter(lig.nurr)
lig.nurr = optimizeALS(lig.nurr, k = 40, lambda = 10)
lig.nurr <- quantile_norm(lig.nurr)
lig.nurr <- rliger::louvainCluster(lig.nurr,resolution = 0.4)
lig.nurr = runUMAP(lig.nurr)
levels.use <- rep('nonda',length(unique(lig.nurr@clusters)))
levels.use[c(4,10,12,13,8)+1] <- 'da'
levels(lig.nurr@clusters) <- levels.use

# Cluster DAPI separately
cleaned.meta.use <- cleaned.meta[which(cleaned.meta$disease == 'Ctrl' & 
                                         grepl('DAPI',cleaned.meta$lib)),]
idx.use <- rownames(cleaned.meta.use)

out.use <- MergeSparseDataAll(lapply(all.obj,function(x){
  x@raw.data[,intersect(colnames(x@raw.data),idx.use)]
}))
out.use <- lapply(unique(cleaned.meta.use$subject),function(x){
  out.use[,which(grepl(x,colnames(out.use)))]
})
names(out.use) <- unique(cleaned.meta.use$subject)
lig.neg <- createLiger(raw.data = out.use)
lig.neg = rliger::normalize(lig.neg)
lig.neg = selectGenes(lig.neg, num.genes = 1000)
lig.neg = scaleNotCenter(lig.neg)
lig.neg = optimizeALS(lig.neg, k = 25, lambda = 5)
lig.neg <- quantile_norm(lig.neg)
lig.neg <- rliger::louvainCluster(lig.neg,resolution = 1.2)
lig.neg = runUMAP(lig.neg,rand.seed = 102)





