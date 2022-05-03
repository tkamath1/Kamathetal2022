library(Seurat)
library(rliger)
library(moments)
library(MAST)
library(dplyr)
library(lme4)
library(devtools)
source("/home/tkamath/scripts/genevalidatehelper.R")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/figures/suppfigures/')
options(mc.cores = 60)
base.dir <- c('/home/tkamath/DA/figures/suppfigures/')

######################################################
################# Supp Figure 3#######################
######################################################

liger.species <- qread('/home/tkamath/DA/species/liger_species_0914_annot.qs')
liger.species@cell.data$dataset <- as.factor(liger.species@cell.data$dataset)
levels(liger.species@cell.data$dataset) <- c('Human','Human','Human',"Human",
                                             'Macaque',"Mouse","Rat",'Tree shrew')
liger.species@cell.data$dataset <- factor(liger.species@cell.data$dataset,levels = c('Human',
                                                                                     'Macaque',"Tree shrew",
                                                                                     'Rat','Mouse'))

idx.use <- sample(rownames(liger.species@cell.data[which(liger.species@cell.data$dataset == 'Macaque'),]),
                  size = 2000)
idx.use <- c(idx.use,rownames(liger.species@cell.data[which(liger.species@cell.data$dataset != 'Macaque'),]))
liger.species2 <- subsetLiger(liger.species,cells.use = idx.use,remove.missing = F)
levels(liger.species2@clusters) <- c('DA1','DA5','DA2','DA4','DA6','DA3','DA7','DA8')

raw.human <- MergeSparseDataAll(liger.species3@raw.data[c(1,6,7,8)])
new.species <- createLiger(raw.data = list('Human' = raw.human, 'Macaque' = liger.species2@raw.data$macaque,
                                           'Tree shrew' = liger.species2@raw.data$treeshrew,
                                           'Rat' = liger.species2@raw.data$rat,
                                           'Mouse' = liger.species2@raw.data$mouse))
new.species <- rliger::normalize(new.species)
new.species@H.norm <- liger.species2@H.norm
new.species@tsne.coords <- liger.species2@tsne.coords[rownames(new.species@cell.data),]

p1 <- rliger::plotGene(new.species,gene = 'Slc18a2',return.plots = T,
                       set.dr.lims = T,pt.size = 0.4)
p1.use <- lapply(p1,function(z){
  z + theme(legend.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),axis.line = element_blank())
})

p1 <- rliger::plotGene(new.species,gene = 'Calb1',return.plots = T,set.dr.lims = T,pt.size = 0.4)
p2.use <- lapply(p1,function(z){
  z + theme(legend.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),axis.line = element_blank())
})

p1 <- rliger::plotGene(new.species,gene = 'Sox6',return.plots = T,set.dr.lims = T,pt.size = 0.4)
p3.use <- lapply(p1,function(z){
  z + theme(legend.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),axis.line = element_blank())
})

pdf('/home/tkamath/DA/figures/suppfigures/suppfig3/species_markers.pdf',
    useDingbats = F,height = 4,width = 24)
grid.arrange(p1.use[[1]],p1.use[[2]],p1.use[[3]],p1.use[[4]],p1.use[[5]],nrow = 1)
grid.arrange(p2.use[[1]],p2.use[[2]],p2.use[[3]],p2.use[[4]],p2.use[[5]],nrow = 1)
grid.arrange(p3.use[[1]],p3.use[[2]],p3.use[[3]],p3.use[[4]],p3.use[[5]],nrow = 1)
dev.off()

seurat.species <- ligerToSeurat(liger.species3)
p1 <- DotPlot(seurat.species,genes.plot = c('Sox6','Calb1','Gem','Fam83b','Nptx2'),do.return = T,plot.legend = T)
pdf('/home/tkamath/DA/figures/suppfigures/suppfig3/dotplot_markerspecies.pdf',useDingbats = F,width = 8, height = 4)
p1 + scale_color_gradient(low = 'grey',high = 'darkgreen') + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

t1 <- table(liger.species3@cell.data$humanannots,liger.species3@clusters)
df.use <- melt(t1/rowSums(t1))
df.use$Var1 <- factor(df.use$Var1,levels = c('SOX6_AGTR1','SOX6_DDT','SOX6_PART1','CALB1_CALCR','SOX6_GFRA2',
                                             'CALB1_PPP1R17','CALB1_CRYM_CCDC68','CALB1_RBP4','CALB1_TRHR','CALB1_GEM'))
df.use$Var2 <- factor(df.use$Var2,levels = c('DA1','DA2','DA3','DA6','DA4', 'DA5','DA7', 'DA8'))

pdf('/home/tkamath/DA/figures/suppfigures/suppfig3/confusionmat_species.pdf',useDingbats = F,width = 12)
ggplot(df.use,aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'black') + 
  scale_fill_gradient(low = 'white',high = 'darkgreen')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + xlab('Human cluster definitions') +
  ylab('Integrative species clusters')
dev.off()