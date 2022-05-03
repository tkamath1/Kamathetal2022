library(Seurat)
library(rliger)
library(moments)
library(MAST)
library(dplyr)
library(lme4)
library(devtools)
library(scico)
library(rliger)
library(SCENIC)
source("/home/tkamath/scripts/genevalidatehelper.R")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/figures/suppfigures/')
options(mc.cores = 60)
base.dir <- c('/home/tkamath/DA/figures/suppfigures/')

######################################################
################# Supp Figure 7 ######################
######################################################
library(ggpubr)
cleaned.meta <- qread('/home/tkamath/DA/cleanedmetadata.qs')
pd.meta.use <- qread('/home/tkamath/DA/pdmeta_use.qs')
cleaned.meta <- bind_rows(cleaned.meta,.id  = 'names')
t2 <- aggregate(nUMI ~ subject, cleaned.meta, median)
t2$status <- as.factor(pd.meta.use[match(t2$subject,pd.meta.use$Donor.ID),]$disease)
p1 <- ggplot(t2,aes(x = status, y = nUMI,fill = status)) + geom_jitter(width = 0.2) + 
  geom_boxplot(width = 0.5,alpha= 0.2) +
  ylim(0,NA) + stat_compare_means(method = "t.test", aes(label = ..p.signif..), 
                                  label.x = 1.5, label.y = 12500, size = 12) + 
  ylab('Median number of UMIs') + theme(legend.position = 'none',
                                        axis.text.x = element_text(size = 20),
                                        axis.text.y = element_text(size = 20),
                                        axis.title = element_text(size = 20))

t1 <- aggregate(nGene ~ subject, cleaned.meta, median)
t1$status <- as.factor(pd.meta.use[match(t1$subject,pd.meta.use$Donor.ID),]$disease)
p2 <- ggplot(t1,aes(x = status, y = nGene,fill = status)) + geom_jitter(width = 0.2) + 
  geom_boxplot(width = 0.5,alpha= 0.2) +
  ylim(0,NA) + stat_compare_means(method = "t.test", aes(label = ..p.signif..), 
                                  label.x = 1.5, label.y = 4500,size = 12) +
  ylab('Median number of genes') + theme(legend.position = 'none',
                                         axis.text.x = element_text(size = 20),
                                         axis.text.y = element_text(size = 20),
                                         axis.title = element_text(size = 20))
pdf('suppfig7/qcmetrics.pdf',useDingbats = F,height = 4,width = 5)
p1
p2
dev.off()

# Age and sex
View(pd.meta.use)
pd.meta.use$age <- c(75,90,76,79,92,83,82,82,79,91,78,49,49,70,90,74,89,82)
pd.meta.use$pmi <- c(13.8,22,6.8,22,23.3,13,13.7,12.9,14,10,22.3,7.2,7,21,2,13,10,6)
qsave(pd.meta.use,'/home/tkamath/DA/pdmeta_use.qs')
summary(lm('age ~ disease',pd.meta.use))
p3 <- ggplot(pd.meta.use,aes(x = disease, y = age,fill = disease))+ geom_jitter(width = 0.2) + 
  geom_boxplot(width = 0.5,alpha= 0.2) +
  ylim(0,NA) + stat_compare_means(method = "t.test", aes(label = ..p.signif..), 
                                  label.x = 1.5, label.y = 95,size = 15) + xlab('status') +
  ylab('Age') + theme(legend.position = 'none',
                      axis.text.x = element_text(size = 20),
                      axis.text.y = element_text(size = 20),
                      axis.title = element_text(size = 20))

p5 <- ggplot(pd.meta.use,aes(x = disease, y = pmi,fill = disease))+ geom_jitter(width = 0.2) + 
  geom_boxplot(width = 0.5,alpha= 0.2) +
  ylim(0,NA) + stat_compare_means(method = "t.test", aes(label = ..p.signif..), 
                                  label.x = 1.5, label.y = 25,size = 15) +
  ylab('PMI') + theme(legend.position = 'none',
                      axis.text.x = element_text(size = 20),
                      axis.text.y = element_text(size = 20),
                      axis.title = element_text(size = 20))

t1 <- table(pd.meta.use$disease,pd.meta.use$sex)
t2 <- t1/rowSums(t1)
t2 <- melt(t2)
p4 <- ggplot(t2,aes(x = Var1, y = value,fill = Var2)) +
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90,size = 20),
                                                     axis.text.y = element_text(size = 20),
                                                     axis.title = element_text(size = 20))

pdf('suppfig7/agemetrics.pdf',useDingbats = F)
p3
dev.off()

pdf('suppfig7/sexmetrics.pdf',useDingbats = F)
p4
dev.off()

pdf('suppfig7/pmimetrics.pdf',useDingbats = F)
p5
dev.off()

#BY major cell type
cleaned.meta[which(grepl('Ex',cleaned.meta$clusters)),]$names <- 'Ex'
cleaned.meta[which(grepl('Inh',cleaned.meta$clusters)),]$names <- 'Inh'
cleaned.meta$names <- as.factor(cleaned.meta$names)
cleaned.meta$names <- factor(cleaned.meta$names,levels = c('Ex','Inh','da','astro','opc','olig','endo','mg'))

dodge <- position_dodge(width = 0.8)
cleaned.meta$celltype <- cleaned.meta$names
levels(cleaned.meta$celltype) <- c('Ex neuron','Inh neuron','DA neuron',
                                   'Astrocyte','OPC','Oligodendrocyte','Endothelial cell',
                                   'Microglia/macrophage')
pdf('suppfig7/qcmetrics_bydisease.pdf',width = 14,height = 8)
ggplot(cleaned.meta,aes(x = celltype, y = nGene,fill = disease)) + geom_violin(position = dodge) + 
  scale_y_log10() + geom_boxplot(width = 0.3,position = dodge,outlier.alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20)) + ylab('Number of genes')

ggplot(cleaned.meta,aes(x = celltype, y = nUMI,fill = disease)) + geom_violin(position = dodge) + 
  scale_y_log10() + geom_boxplot(width = 0.3,position = dodge,outlier.alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 20),        
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20)) + ylab('Number of UMIs')
dev.off()

# Alignment scores
alignment.df <- data.frame('Cell_type' = as.factor(c('DA neuron','Non-DA neuron','Astrocyte','Endothelial cell/pericyte',
                                                     'Microglia/macrophage','Oligodendrocyte','OPC')),
                           'Alignment' = c(1,0.99,0.74,0.61,0.76,0.76,0.87))
alignment.df$Cell_type <- factor(alignment.df$Cell_type,levels = rev(c('DA neuron','Non-DA neuron','Astrocyte',
                                                                       'Endothelial cell/pericyte','Microglia/macrophage',
                                                                       'Oligodendrocyte','OPC')))

pdf('suppfig7/alignmentscores_celltype.pdf',width = 10)
ggplot(alignment.df, aes(x = Cell_type, y= Alignment)) + geom_bar(stat = 'identity') +
  coord_flip() + ylab('Alignment score') + xlab('Cell type') +
  theme(axis.text.x = element_text(size = 20),        
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20))
dev.off()

lapply(unique(cleaned.meta$names),function(x){
  t1 <- aggregate(nGene ~ subject,cleaned.meta[cleaned.meta$names == x,],median)
  t1$status <- as.factor(pd.meta.use[match(t1$subject,pd.meta.use$Donor.ID),]$disease)
  summary(lm('nGene ~ status',t1))
  
  t1 <- aggregate(nUMI ~ subject,cleaned.meta[cleaned.meta$names == x,],median)
  t1$status <- as.factor(pd.meta.use[match(t1$subject,pd.meta.use$Donor.ID),]$disease)
  summary(lm('nUMI ~ status',t1))
})

all.obj <- qread('/home/tkamath/DA/cleanedobjects.qs')
## Astro
p1 <- DimPlot(all.obj$astro,reduction.use = 'umap',pt.size = 0.3,do.label = T,no.legend = T,do.return = T)
p1.use <- DimPlot(all.obj$astro,reduction.use = 'umap',pt.size = 0.3,do.label = F,no.legend = T,
                  do.return = T)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)
p1.astro <- p1.use + scale_color_manual(values = mycolors) + 
  geom_text_repel(data = p1$layers[[2]]$data,aes(x = x, y = y,label = ident),size = 6,inherit.aes = F)+ 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

p2 <- DimPlot(all.obj$astro,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'subject')
p2.astro <- p2 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank(),
        legend.text = element_text(size = 15))

all.obj$astro@meta.data$status <- as.factor(pd.meta.use[match(all.obj$astro@meta.data$subject,
                                                              pd.meta.use$Donor.ID),]$disease)
p3 <- DimPlot(all.obj$astro,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'status')
p3.astro <- p3 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

pdf('suppfig7/astrocyte_tsne.pdf',useDingbats = F,width = 15)
grid.arrange(p1.astro,p2.astro,nrow = 1)
dev.off()

all.obj$astro@dr$cca <- all.obj$astro@dr$hpca
astro.liger <- seuratToLiger(all.obj$astro,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
astro.liger@H.norm <- all.obj$astro@dr$hpca@cell.embeddings
calcAlignment(astro.liger) # 0.74


## Endo
p1 <- DimPlot(all.obj$endo,reduction.use = 'umap',pt.size = 0.3,do.label = T,no.legend = T,do.return = T)
p1.use <- DimPlot(all.obj$endo,reduction.use = 'umap',pt.size = 0.3,do.label = F,no.legend = T,
                  do.return = T)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)
p1.endo <- p1.use + scale_color_manual(values = mycolors) + 
  geom_text_repel(data = p1$layers[[2]]$data,aes(x = x, y = y,label = ident),size = 6,inherit.aes = F)+ 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

p2 <- DimPlot(all.obj$endo,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'subject')
p2.endo <- p2 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank(),legend.text = element_text(size = 15))

all.obj$endo@meta.data$status <- as.factor(pd.meta.use[match(all.obj$endo@meta.data$subject,
                                                             pd.meta.use$Donor.ID),]$disease)
p3 <- DimPlot(all.obj$endo,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'status')
p3.endo <- p3 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

pdf('endo_tsne.pdf',useDingbats = F,width = 15)
grid.arrange(p1.endo,p2.endo,nrow = 1)
dev.off()

all.obj$endo@dr$cca <- all.obj$endo@dr$hpca
endo.liger <- rliger::seuratToLiger(all.obj$endo,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
endo.liger@H.norm <- all.obj$endo@dr$hpca@cell.embeddings
calcAlignment(endo.liger) # 0.61


#### mg
p1 <- DimPlot(all.obj$mg,reduction.use = 'umap',pt.size = 0.3,do.label = T,no.legend = T,do.return = T)
p1.use <- DimPlot(all.obj$mg,reduction.use = 'umap',pt.size = 0.3,do.label = F,no.legend = T,
                  do.return = T)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(12)
p1.mg <- p1.use + scale_color_manual(values = mycolors) + 
  geom_text_repel(data = p1$layers[[2]]$data,aes(x = x, y = y,label = ident),size = 6,inherit.aes = F)+ 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

p2 <- DimPlot(all.obj$mg,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'subject')
p2.mg <- p2 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank(),legend.text = element_text(size = 15))

all.obj$mg@meta.data$status <- as.factor(pd.meta.use[match(all.obj$mg@meta.data$subject,
                                                           pd.meta.use$Donor.ID),]$disease)
p3 <- DimPlot(all.obj$mg,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'status')
p3.mg <- p3 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

pdf('mg_tsne.pdf',useDingbats = F,width = 15)
grid.arrange(p1.mg,p2.mg,nrow = 1)
dev.off()

all.obj$mg@dr$cca <- all.obj$mg@dr$hpca
mg.liger <- seuratToLiger(all.obj$mg,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
mg.liger@H.norm <- all.obj$mg@dr$hpca@cell.embeddings
calcAlignment(mg.liger) # 0.76

#### Oligo
p1 <- DimPlot(all.obj$olig,reduction.use = 'umap',pt.size = 0.3,do.label = T,no.legend = T,do.return = T)
p1.use <- DimPlot(all.obj$olig,reduction.use = 'umap',pt.size = 0.3,do.label = F,no.legend = T,
                  do.return = T)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(12)
p1.olig <- p1.use + scale_color_manual(values = mycolors) + 
  geom_text_repel(data = p1$layers[[2]]$data,aes(x = x, y = y,label = ident),size = 6,inherit.aes = F)+ 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

p2 <- DimPlot(all.obj$olig,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'subject')
p2.olig <- p2 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank(),legend.text = element_text(size = 15))

all.obj$olig@meta.data$status <- as.factor(pd.meta.use[match(all.obj$olig@meta.data$subject,
                                                             pd.meta.use$Donor.ID),]$disease)
p3 <- DimPlot(all.obj$olig,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'status')
p3.olig <- p3 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

pdf('olig_tsne.pdf',useDingbats = F,width = 15)
grid.arrange(p1.olig,p2.olig,nrow = 1)
dev.off()

all.obj$olig@dr$cca <- all.obj$olig@dr$hpca
olig.liger <- seuratToLiger(all.obj$olig,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
olig.liger@H.norm <- all.obj$olig@dr$hpca@cell.embeddings
calcAlignment(olig.liger) # 0.76

## OPC
p1 <- DimPlot(all.obj$opc,reduction.use = 'umap',pt.size = 0.3,do.label = T,no.legend = T,do.return = T)
p1.use <- DimPlot(all.obj$opc,reduction.use = 'umap',pt.size = 0.3,do.label = F,no.legend = T,
                  do.return = T)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(12)
p1.opc <- p1.use + scale_color_manual(values = mycolors) + 
  geom_text_repel(data = p1$layers[[2]]$data,aes(x = x, y = y,label = ident),size = 6,inherit.aes = F) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

p2 <- DimPlot(all.obj$opc,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'subject')
p2.opc <- p2 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank(),legend.text = element_text(size = 15))

all.obj$opc@meta.data$status <- as.factor(pd.meta.use[match(all.obj$opc@meta.data$subject,
                                                            pd.meta.use$Donor.ID),]$disease)
p3 <- DimPlot(all.obj$opc,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'status')
p3.opc <- p3 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

pdf('opc_tsne.pdf',useDingbats = F,width = 15)
grid.arrange(p1.opc,p2.opc,nrow = 1)
dev.off()

all.obj$opc@dr$cca <- all.obj$opc@dr$hpca
opc.liger <- rliger::seuratToLiger(all.obj$opc,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
opc.liger@H.norm <- all.obj$opc@dr$hpca@cell.embeddings
calcAlignment(opc.liger)

# Non-DA neuron
names(all.obj$nonda@dr)[1] <- 'umap'
colnames(all.obj$nonda@dr$umap@cell.embeddings) <- c('UMAP_1','UMAP_2')
all.obj$nonda@dr$umap@key <- 'UMAP_'
rownames(all.obj$nonda@dr$umap@cell.embeddings) <- rownames(all.obj$nonda@meta.data)

p1 <- DimPlot(all.obj$nonda,pt.size = 0.3,reduction.use = 'umap',
              do.label = T,no.legend = T,do.return = T)
p1.use <- DimPlot(all.obj$nonda,reduction.use = 'umap',pt.size = 0.3,do.label = F,no.legend = T,
                  do.return = T)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(17)
p1.nonda <- p1.use + scale_color_manual(values = mycolors) + 
  geom_text_repel(data = p1$layers[[2]]$data,aes(x = x, y = y,label = ident),size = 4,inherit.aes = F)+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

p2 <- DimPlot(all.obj$nonda,reduction.use = 'umap',pt.size = 0.01,
              do.label = F,no.legend = F,do.return = T,group.by = 'orig.ident')
p2.nonda <- p2 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank(),legend.text = element_text(size= 15))

all.obj$nonda@meta.data$status <- as.factor(pd.meta.use[match(all.obj$nonda@meta.data$orig.ident,
                                                              pd.meta.use$Donor.ID),]$disease)
p3 <- DimPlot(all.obj$nonda,reduction.use = 'umap',pt.size = 0.1,
              do.label = F,no.legend = F,do.return = T,group.by = 'status')
p3.nonda <- p3 + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())

pdf('nonda_tsne.pdf',useDingbats = F,width = 15)
grid.arrange(p1.nonda,p2.nonda,nrow = 1)
dev.off()

all.obj$nonda@dr$cca <- all.obj$nonda@dr$hpca
nonda.liger <- rliger::seuratToLiger(all.obj$nonda,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
nonda.liger@H.norm <- all.obj$nonda@dr$hpca@cell.embeddings
calcAlignment(nonda.liger)

nonda.liger <- qread('/home/tkamath/DA/nonda/nonda_liger_annot_0229.qs')
p1 <- plotByDatasetAndCluster(nonda.liger,return.plots = T,pt.size = 0.8,text.size = 0)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(20)
p1.nonda <- p1[[2]] + theme(legend.position ='none') + scale_color_manual(values = mycolors) + 
  geom_label_repel(data = p1[[2]]$layers[[2]]$data,aes(x = tsne1, y = tsne2,label = Cluster),size = 3,
                   label.size = 0.5,label.r = 0,inherit.aes = F)
