###########################
### Species integraiton ###
###########################
library(Seurat)
library(rliger)
library(NNLM)
library(moments)
library(presto)
source('/home/tkamath/scripts/SeuratExtrafunctions.R')
source('/home/tkamath/scripts/extrafuncs.R')
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/species/')

sn.da.annot <- qread('/home/tkamath/DA/da/sn_da_annot_021821.qs')
sn.da.ctrl <- subsetLiger(sn.da.annot,cells.use = rownames(sn.da.annot@cell.data[sn.da.annot@cell.data$status == 'Ctrl',]),
                          remove.missing = F)

# Species projection final
a.rat <- readRDS('rat/rat_SN_TH_cleaned.rds')
a.macaque <- qread('newmacaque/newmacaque_th_cleaned_081621.qs')
a.tshrew <- readRDS('tree_shrew/treeshrew_SN_TH.rds')
a.mh <- readRDS("species/species_da.rds")
a.6173 <- readRDS('6173_seurat.rds')
a.6173THuse <- qread('../6173/6173_seurat_th.qs')

new.mouse <- qread('Midbrain_DAcells_liger.qs')
data.mouse.use <- MergeSparseDataAll(new.mouse@raw.data)
data.mouse.use <- data.mouse.use[,sample(colnames(data.mouse.use),size = 1500)]
mh.data <- mousehumanhomologs(data.human = MergeSparseDataAll(sn.da.ctrl@raw.data),
                              data.mouse = data.mouse.use)
rownames(mh.data$MOUSE) <- toupper(rownames(mh.data$MOUSE))
a.mhnew <- rliger::createLiger(list('mouse' = mh.data$MOUSE))
a.mhnew <- rliger::normalize(a.mhnew)
a.mhnew@var.genes <- sn.da.ctrl@var.genes
genes.add <- setdiff(a.mhnew@var.genes,rownames(a.mhnew@norm.data$mouse))
matrix.add <- matrix(data = 0,nrow = length(genes.add),ncol = length(colnames(a.mhnew@norm.data$mouse)))
colnames(matrix.add) <- rownames(a.mhnew@cell.data)
rownames(matrix.add) <-genes.add
a.mhnew@norm.data$mouse <- rbind(a.mhnew@norm.data$mouse,matrix.add)
a.mhnew@raw.data$mouse <- rbind(a.mhnew@raw.data$mouse,matrix.add)
#a.mhnew <- rliger::selectGenes(a.mhnew )
a.mhnew <- scaleNotCenter(a.mhnew)
sn.da.ctrl@cell.data <- sn.da.ctrl@cell.data[,c(1,2,3)]
species.proj2 = online_iNMF(sn.da.ctrl, X_new = list(mouse = a.mhnew),projection = T)

species.proj2 = quantile_norm(species.proj2,ref_dataset = '6173',quantiles = 30)
species.proj2 <- louvainCluster(species.proj2,resolution = 2.5)
species.proj2 = runUMAP(species.proj2)
plotByDatasetAndCluster(species.proj2,return.plots = T)
rliger::plotGene(species.proj2,gene = '',plot.by = 'none')

table(species.proj2@cell.data$dataset,species.proj2@clusters)
pdf('factorplot_0526.pdf')
plotFactors_new(species.proj2,plot.tsne = T)
dev.off()

liger::plotGene(a.mhnew, gene = 'Reln')
qsave(a.mhnew,'mousehuman_new.qs')
a.mhnew2 <- createLiger(raw.data = a.mhnew@raw.data)
a.mhnew2@var.genes <- a.mhnew@var.genes
a.mhnew2@H.norm <- a.mhnew@H.norm
a.mhnew2@H <- a.mhnew@H
a.mhnew2@W <- a.mhnew@W
a.mhnew2@V <- a.mhnew@V
a.mhnew2 <- rliger::normalize(a.mhnew2)
a.mhnew2 <- scaleNotCenter(a.mhnew2)
a.mhnew2@cell.data$lib <- as.factor(unlist(sapply(strsplit(rownames(a.mhnew2@cell.data),'_'),`[`,1)))
idx.remove <- rownames(a.mhnew2@cell.data[which(a.mhnew2@cell.data$lib == '6173DAPI'),])
a.mhnew3 <- rliger::subsetLiger(a.mhnew2, cells.use = setdiff(rownames(a.mhnew2@cell.data),idx.remove),
                                remove.missing = F)

a.mhnew2 <- louvainCluster(a.mhnew2, resolution = 0.4)
a.mhnew2 <- runTSNE(object = a.mhnew2)

# Add in other datasets
rownames(a.macaque@raw.data) <- firstupper(tolower(rownames(a.macaque@raw.data)))
a.macaque@meta.data$dataset <- 'macaque'
a.macaque <- seuratToLiger(a.macaque,combined.seurat = T,meta.var = 'dataset',remove.missing = F)
a.macaque@var.genes <- a.mhnew@var.genes
genes.add <- setdiff(a.mhnew2@var.genes,rownames(a.macaque@raw.data$macaque))

mat.add <- Matrix(data = 0, nrow = length(genes.add),ncol = nrow(a.macaque@cell.data))
rownames(mat.add) <- genes.add
colnames(mat.add) <- rownames(a.macaque@cell.data)
a.macaque@raw.data$macaque <- rbind(a.macaque@raw.data$macaque,mat.add)
a.macaque <- rliger::normalize(a.macaque)
a.macaque <- scaleNotCenter(a.macaque)

a.rat@meta.data$dataset <- 'rat'
a.rat <- seuratToLiger(a.rat,combined.seurat = T,meta.var = 'dataset',remove.missing = F)
a.rat@var.genes <- a.mhnew@var.genes
a.rat <- scaleNotCenter(a.rat)

a.mhnew2@cell.data$lib <- NULL
species.proj2 = online_iNMF(a.mhnew2, X_new = list(macaque = a.macaque),
                            projection = T)
species.proj2 = online_iNMF(species.proj2, X_new = list(rat = a.rat),
                            projection = T,miniBatch_size = 50)

species.proj2 = quantile_norm(species.proj2,ref_dataset = 'human',knn_k = 45)
species.proj2 = runUMAP(species.proj2)
species.proj2 <- louvainCluster(species.proj2,resolution = 0.3)
plotByDatasetAndCluster(species.proj2,return.plots = T)


rownames(a.tshrew@raw.data) <- firstupper(tolower(rownames(a.tshrew@raw.data)))
a.tshrew@meta.data$dataset <- 'treeshrew'
a.tshrew <- seuratToLiger(a.tshrew,combined.seurat = T,meta.var = 'dataset',remove.missing = F)
a.tshrew@var.genes <- a.mhnew2@var.genes
a.tshrew <- scaleNotCenter(a.tshrew)
species.proj2 = online_iNMF(species.proj2, X_new = list(treeshrew = a.tshrew),
                            projection = T)
species.proj2 = quantile_norm(species.proj2,ref_dataset = 'human')
species.proj2 = runUMAP(species.proj2)
species.proj2 <- louvainCluster(species.proj2,resolution = 0.5,eps = 0.01)
plotByDatasetAndCluster(species.proj2,return.plots = T)
pdf('tsne_updated_species.pdf',useDingbats = F)
p1 <- plotByDatasetAndCluster(species.proj2,return.plots = T,pt.size = 0.03)
p1[[1]]
p1[[2]]
t1 <- table(species.proj2@cell.data$dataset,species.proj2@clusters)
t2 <- as.data.frame(melt(t1/rowSums(t1)))
t2$Var1 <- as.factor(t2$Var1)
ggplot(t2, aes(fill=Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity")
liger::plotGene(species.proj2,gene = 'Fam83b',plot.by = 'none')
liger::plotGene(species.proj2,gene = 'Reln',plot.by = 'none')
dev.off()

liger::plotGene(species.proj2,gene = 'Fam83b',plot.by = 'none')

p1 <- rliger::plotByDatasetAndCluster(species.proj2,return.plots = T,pt.size = 0.5,do.shuffle = T)
p1[[1]]$data$Dataset <- as.factor(p1[[1]]$data$Dataset)
levels(p1[[1]]$data$Dataset) <- c('Primate','Primate','Non-primate','Non-primate','Non-primate')
p1[[1]] + scale_color_manual(values = c('Primate' = 'red','Non-primate' = 'blue'))
species.proj2 <- qread('/home/tkamath/DA/species/liger_species_012621.qs')
species.proj2 <- convertOldLiger(species.proj2,override.raw = T)

species.proj.use <- rliger::createLiger(raw.data = species.proj2@raw.data)
species.proj.use@norm.data <- species.proj2@norm.data
species.proj.use@cell.data <- species.proj2@cell.data
species.proj.use@H.norm <- species.proj2@H.norm
species.proj.use@H <- species.proj2@H
species.proj.use@W <- species.proj2@W
species.proj.use@V <- species.proj2@V
species.proj.use@var.genes <- species.proj2@var.genes
species.proj.use@tsne.coords <- species.proj2@tsne.coords
species.proj.use@clusters <- species.proj2@clusters
species.proj.use@scale.data <- species.proj2@scale.data

plotByDatasetAndCluster(species.proj.use)
rliger::plotGene(species.proj.use, gene = "Fam83b",plot.by = 'none')
qsave(species.proj.use,'/home/tkamath/DA/species/liger_species_0321.qs')

sn.da.ctrl <- qread('/home/tkamath/DA/da/sn_da_ctrl_021821.qs')
species.proj.use <- qread('/home/tkamath/DA/species/liger_species_091321.qs')

dataset.list <- unique(sn.da.ctrl@cell.data$dataset)[c(1,6,7)]
for (x in c(1:length(dataset.list))) {
  idx.use <- rownames(sn.da.ctrl@cell.data[which(sn.da.ctrl@cell.data$dataset == dataset.list[x]),])
  sn.individual <- rliger::subsetLiger(sn.da.ctrl,cells.use = idx.use,remove.missing = F)
  rownames(sn.individual@raw.data[[1]]) <- firstupper(tolower(rownames(sn.individual@raw.data[[1]])))
  sn.individual <- rliger::normalize(sn.individual)
  sn.individual@var.genes <- species.proj.use@var.genes
  genes.add <- setdiff(species.proj.use@var.genes,rownames(sn.individual@raw.data[[1]]))
  matrix.add <- Matrix(matrix(0,nrow = length(genes.add),ncol = ncol(sn.individual@norm.data[[1]])),sparse = T)
  rownames(matrix.add) <- genes.add
  colnames(matrix.add) <- rownames(sn.individual@cell.data)
  sn.individual@norm.data <- list(rbind(sn.individual@norm.data[[1]],matrix.add))
  names(sn.individual@norm.data) <- as.character(x)
  sn.individual <- scaleNotCenter(sn.individual)
  sn.individual@cell.data$percent.mito <- NULL
  lis.use <- list(sn.individual)
  names(lis.use) <- dataset.list[x]
  species.proj.use = online_iNMF(species.proj.use, X_new = lis.use,project = TRUE)
}
species.proj.use <- quantile_norm(species.proj.use,ref_dataset = 'human',eps = 0.1)
species.proj.use <- louvainCluster(species.proj.use)
species.proj.use = runUMAP(species.proj.use)


# Remove: 24, 23,22,20,19,16
rliger::plotGene(species.proj.use,gene = 'Calcr',plot.by = 'none')
species.proj.use2 <- subsetLiger(species.proj.use,clusters.use = setdiff(levels(
  species.proj.use@clusters),c(16,19,20,22,23,24)),remove.missing = F)
species.proj.use2 <- louvainCluster(species.proj.use2,resolution = 1.2)
species.proj.use2 = runUMAP(species.proj.use2)
p1 <- rliger::plotByDatasetAndCluster(species.proj.use2,return.plots = T,pt.size = 0.1,do.shuffle = T)
qsave(species.proj.use2,'liger_species_allfinal_091321.qs')
rliger::plotGene(species.proj.use2,gene = 'Sox6',plot.by = 'none')
species.proj.use2 <- subsetLiger(species.proj.use2,clusters.use = setdiff(levels(species.proj.use2@clusters),
                                                                          c(11)))

qsave(species.proj.use2,'liger_species_0914_allhumans.qs')

t2 <- species.proj.annot2@cell.data %>% mutate(primate = ifelse(dataset %in% c('mouse','rat','treeshrew'),
                                                                'non-primate','primate'))
species.proj.annot2@cell.data <- t2
t1 <- table(t2$primate,species.proj.annot2@clusters)
df.use <- melt(t1/rowSums(t1))
ggplot(df.use,aes(x = Var2, y = value, fill = Var1)) + geom_bar(position="fill", stat="identity")

qsave(species.proj.annot2,'/home/tkamath/DA/species/liger_species_0914_annot.qs')
