## Figure 1 plots
library(Seurat)
library(rliger)
library(moments)
library(Matrix)
library(broom)
library(qs)
library(MAST)
library(dplyr)
library(lme4)
library(NNLM)
library(devtools)
library(ggrepel)
library(Rtsne)
source("/home/tkamath/scripts/genevalidatehelper.R")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/figures/')
options(mc.cores = 60)

######################################
### Fig 1B ###########################
######################################
lig.nurr <- qread('/home/tkamath/DA/nurrintegration/ligall_nurrpos.qs')
lig.neg <-qread('/home/tkamath/DA/nurrintegration/ligall_nurrneg.qs')

levels(lig.neg@clusters) <-  c('Non-DA','DA neuron')
p2 <- rliger::plotByDatasetAndCluster(lig.neg,return.plots = T,text.size = 0)
p2 <- p2[[2]] +  scale_color_manual(values = c("DA neuron" = "#1B9E77", "Non-DA" = "grey")) +
  theme(legend.position = 'none') 

levels(lig.nurr@clusters) <-  c('Non-DA','DA neuron')
p1 <- rliger::plotByDatasetAndCluster(lig.nurr,return.plots = T,text.size = 0)
p1 <- p1[[2]] +  scale_color_manual(values = c("DA neuron" = "#1B9E77", "Non-DA" = "grey")) +
  theme(legend.position = 'none')

# Nurr enrichment
pdf('/home/tkamath/DA/figures/nurrintegration_fig.pdf',useDingbats = F)
p1
p2
dev.off()

#######################################
###### Fig 1C #########################
#######################################
t1 <- table(lig.neg@cell.data$dataset,lig.neg@clusters)
t2 <- as.data.frame(melt(t1/rowSums(t1)))
t3 <- t2[t2$Var2 == 'DA neuron',]

t1 <- table(lig.nurr@cell.data$dataset,lig.nurr@clusters)
t2 <- as.data.frame(melt(t1/rowSums(t1)))
t4 <- t2[t2$Var2 == 'DA neuron',]

t5 <- merge(t3,t4, by = 'Var1')
t5$Var2.y <- NULL
colnames(t5) <- c('ID','celltype','nurrneg','nurrpos')
t5 <- melt(t5)
t5 <- t5[t5$variable != 'ID',]
t6 <- aggregate(value ~ variable, t5, mean)
pdf('proportiondaneurons.pdf',width = 3.5, height = 8,useDingbats = F)
ggplot(t5, aes(x = variable, y = value)) + geom_point() + geom_bar(data = t6, aes(x = variable, y = value,fill = variable),
                                                                   inherit.aes = F,stat = 'identity',alpha = 0.5) +
  scale_fill_manual(values= c('orange','purple'))+
  ylim(NA,1) + xlab('Library') + ylab('Proportion of DA neurons per replicate') + theme(legend.position = 'none')
dev.off()

##########################################
########## Fig 1D ########################
##########################################

## Dopamine neuron subtypes
sn.da.annot <- qread('/home/tkamath/DA/da/sn_da_annot_021821.qs')
# First plot the control datasets
idx.use <- rownames(sn.da.annot@cell.data[which(sn.da.annot@cell.data$status == 'Ctrl'),])
sn.da.ctrl <- rliger::subsetLiger(sn.da.annot,cells.use = idx.use)

pdf('/home/tkamath/DA/figures/daneurons.pdf',useDingbats = F)
p1 <- plotByDatasetAndCluster(sn.da.ctrl,return.plots = T,pt.size = 0.8,text.size = 0)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)
p1[[2]] + theme(legend.position ='none') + scale_color_manual(values = mycolors) + 
  geom_label_repel(data = p1[[2]]$layers[[2]]$data,aes(x = tsne1, y = tsne2,label = Cluster),size = 3,
                   label.size = 0.5,label.r = 0,inherit.aes = F)
rliger::plotGene(sn.da.ctrl,gene = 'SOX6',points.only = T,plot.by = 'none',pt.size = 1)
rliger::plotGene(sn.da.ctrl,gene = 'CALB1',points.only = T,plot.by = 'none',pt.size = 1)
rliger::plotGene(sn.da.ctrl,gene = 'AGTR1',points.only = T,plot.by = 'none',pt.size = 1)
rliger::plotGene(sn.da.ctrl,gene = 'TMEM200A',points.only = T,plot.by = 'none',pt.size = 1)
dev.off()

###############################
#### Fig 1 E ##################
###############################

#Plot dots
sn.da.seurat <- ligerToSeurat(sn.da.ctrl)
p2 <- DotPlot(sn.da.seurat,genes.plot = c('SOX6','CALB1',"ALDH1A1", 'TMEM200A','AGTR1', 'SYNPR','FAM19A4',
                                          'PART1','DDT','GFRA2','CALCR','CRYM','CCDC68',
                                          'PPP1R17','RBP4','GEM','TRHR'),x.lab.rot = T,
              do.return = T,plot.legend = T)
p2$data$id <- factor(p2$data$id,levels = c('SOX6_AGTR1','SOX6_PART1','SOX6_DDT','SOX6_GFRA2',
                                           'CALB1_CALCR','CALB1_CRYM_CCDC68','CALB1_PPP1R17','CALB1_RBP4',
                                           'CALB1_GEM','CALB1_TRHR'))
limit <- max(abs(p2$data$avg.exp.scale)) * c(-1, 1)
pdf('~/DA/figures/fig1/dotplot.pdf',width = 8,height = 5,useDingbats = F)
p2 + scale_color_gradient(low = 'grey',high = "darkgreen") + scale_shape_manual() + coord_flip() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2 + scale_color_scico(palette = 'cork',limit = limit)  + scale_shape_manual() + coord_flip() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()


#Species
species.da <- qread('/home/tkamath/DA/species/liger_species_0321_allhumans.qs')
species.da.use <- subsetLiger(species.da, clusters.use = setdiff(levels(species.da@clusters),c('REMOVE')),
                              remove.missing = F)
levels(species.da.use@cell.data$dataset) <- c('Human','Mouse','Macaque','Rat','Tree shrew',
                                              'Human','Human','Human')
idx.use <- sample(rownames(species.da.use@cell.data[which(species.da.use@cell.data$dataset == 'Human'),]),
                  size = 1400)
idx.use <- c(idx.use,rownames(species.da.use@cell.data[which(species.da.use@cell.data$dataset != 'Human'),]))
species.da.use2 <- subsetLiger(species.da.use,cells.use = idx.use,remove.missing = F)
species.da.use2 <- runTSNE(species.da.use2)
p2 <- rliger::plotByDatasetAndCluster(species.da.use2,
                                      return.plots = T,legend.size = 0,text.size = 0,pt.size = 1.2)
p2[[1]]$data$Dataset <- as.factor(p2[[1]]$data$Dataset)
levels(p2[[1]]$data$Dataset) <- c('Human','Human','Human','Human','Macaque',
                                  'Mouse','Rat','Tree shrew')
mycolors <- colorRampPalette(brewer.pal(5, "Dark2"))(5)
p.use1 <- p2[[1]] + scale_color_manual(values = mycolors) + theme(legend.position = 'none')

p2 <- rliger::plotByDatasetAndCluster(species.da.use2,return.plots = T,legend.size = 0,text.size = 0,
                                      pt.size = 1)
p2[[1]]$data$Dataset <- as.factor(p2[[1]]$data$Dataset)
levels(p2[[1]]$data$Dataset) <- c('Primate', 'Primate', 'Primate','Primate','Primate',
                                  'Non-primate','Non-primate','Non-primate')
p.use6 <- p2[[1]]  + 
  theme(legend.position = 'none')

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(8)
p.use2 <- p2[[2]] + scale_color_manual(values = mycolors)  + 
  geom_label_repel(data = p2[[2]]$layers[[2]]$data,aes(x = tsne1, y = tsne2,label = Cluster),size = 7,
                   label.size = 1.5,label.r = 0,inherit.aes = F)+ theme(legend.position = 'none')

t1 <- table(species.da.use2@clusters,species.da.use2@cell.data$dataset)
t2 <- as.data.frame(t1/rowSums(t1))
t2$Var2 <- as.factor(t2$Var2)
levels(t2$Var2) <- c('Primate','Non-primate','Primate','Non-primate','Non-primate')
t2$Var1 <- factor(t2$Var1,levels = c('DA7','DA5','DA2','DA3','DA6','DA4','DA1','DA8'))
p.use3 <- ggplot(t2, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust= 1)) + 
  xlab('Cluster') +ylab('Ratio per cluster')

p.use4 <- ggplot(t2, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity")  + 
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust= 1)) + 
  xlab('Cluster') +ylab('Percent makeup of cluster')

species.seurat <- ligerToSeurat(species.da.use)
levels(species.seurat@meta.data$orig.ident) <- c('Human','Mouse','Macaque','Rat','Tree shrew',
                                                 'Human','Human','Human')
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(8)
p.use5 <- DotPlot(species.seurat, genes.plot = rev(c('Aqp3','Gem','Reln','Calb1','Sox6')),
                  plot.legend = T,scale.by = 'size',do.return = T)
p.use5 <- p.use5 + scale_color_gradient(low = 'grey',high = "darkgreen")
pdf('/home/tkamath/DA/figures/figure1_species.pdf',useDingbats = F)
p.use1 + theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
               axis.title.x = element_blank(),axis.title.y = element_blank(),axis.line = element_blank())
p.use2+ theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
              axis.title.x = element_blank(),axis.title.y = element_blank(),axis.line = element_blank())
dev.off()

pdf('/home/tkamath/DA/figures/figure1_species2.pdf',useDingbats = F)
p.use3
p.use4
dev.off()

pdf('/home/tkamath/DA/figures/figure1_species3.pdf',useDingbats = F,width = 6, height = 3)
p.use5 + coord_flip()
dev.off()
pdf('/home/tkamath/DA/figures/figure1_species4.pdf',useDingbats = F)
p.use6+ theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
              axis.title.x = element_blank(),axis.title.y = element_blank(),axis.line = element_blank())
dev.off()

# Redo with new macaque data
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
idx.use <- sample(rownames(liger.species@cell.data[which(liger.species@cell.data$dataset == 'Human'),]),
                  size = 1400)
idx.use <- c(idx.use,rownames(liger.species@cell.data[which(liger.species@cell.data$dataset != 'Human'),]))
liger.species3 <- subsetLiger(liger.species2,cells.use = idx.use,remove.missing = F)


liger.species3 <- runUMAP(liger.species3)
levels(liger.species3@clusters) <- c('DA1','DA5','DA2','DA4','DA6','DA3','DA7','DA8')
p1 <- plotByDatasetAndCluster(liger.species3,return.plots = T,pt.size = 0.8,text.size = 8)
p1[[1]]$data$Dataset <- as.factor(p1[[1]]$data$Dataset)
levels(p1[[1]]$data$Dataset) <- c('Human','Human','Human','Human',
                                  'Macaque','Mouse','Rat','Tree shrew')
p1[[1]]$data$Dataset <- factor(p1[[1]]$data$Dataset,levels = c('Human',
                                                               'Macaque',"Tree shrew",
                                                               'Rat','Mouse'))

pdf('/home/tkamath/DA/species/tsne_0914_2.pdf',useDingbats = F)
p1[[1]] + scale_color_brewer(palette = 'Dark2') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),legend.position = 'none',
        axis.line = element_blank())
p1[[2]] + scale_color_brewer(palette = 'Dark2')+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),legend.position = 'none',
        axis.line = element_blank())
levels(p1[[1]]$data$Dataset) <- c('Primate','Primate','Non-primate','Non-primate','Non-primate')
p1[[1]] +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),legend.position = 'none',
        axis.line = element_blank())
dev.off()

liger.species3@cell.data <- liger.species3@cell.data %>% mutate(primate = ifelse(dataset %in% c('Mouse','Rat','Tree shrew'),'non-primate',
                                                                                 'primate'))
t1 <- table(liger.species3@clusters,liger.species3@cell.data$primate)
t2 <- melt(t1/rowSums(t1))
t2$Var2 <- factor(t2$Var2,levels = c('primate','non-primate'))

qsave(liger.species3,'~/DA/species/liger_species_0426_annot.qs')

pdf('/home/tkamath/DA/species/barplot_species_0914.pdf',useDingbats = F,width = 4)
ggplot(t2,aes(x = Var1, y = value, fill = Var2)) + geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + theme(legend.position = 'none')
dev.off()

##########################################
######### Fig 1I #########################
###########################################
# TF analysis
to.plot.df <- qread('/home/tkamath/DA/tfanalysis/tf_toplot.qs')
limit <- max(abs(to.plot.df$value)) * c(-1, 1)
pdf('/home/tkamath/DA/figures/fig1/tfs_toplot.pdf',height = 6,width = 6)
ggplot(to.plot.df,aes(y = DA_type, x = Regulon, fill = value)) + geom_tile() +
  scale_fill_scico(palette = 'cork',limit = limit) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + coord_flip()

ggplot(to.plot.df,aes(y = DA_type, x = Regulon, fill = value)) + geom_tile() +
  scale_fill_gradient2(low = 'white',high = 'darkgreen',midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + coord_flip()
dev.off()


