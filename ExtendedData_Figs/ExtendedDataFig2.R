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

##################################################
################## SUPP FIGURE 2 #################
##################################################
###################
## 2A #############
###################
# Relative total number of DA neurons
sn.ligerpaper <- readRDS('/home/tkamath/DA/ligerpaper/SN_7ids_analogizer.rds')
t1 <- table(sn.ligerpaper@clusters)
cleaned.meta$subject <- droplevels(cleaned.meta$subject)
t2 <- table(cleaned.meta$names)
thisstudy.vec <- as.vector(t2['DA neuron'])
welch.vec <- t1['NEUROdop']
smajic.vec <- sum(c(22,0,0,2,5,3))
agrawal.vec <- sum(c(0,0,0,7,18,22,2,23))
lamanno.vec <- sum(c(47,38,37))
list.use <- list('This study' = thisstudy.vec,'Welch et al 2019' = welch.vec,
                 'Smajic et al 2021' = smajic.vec, 'Agrawal et al 2020' = agrawal.vec,
                 'La Manno et al 2016' = lamanno.vec)
list.df <- melt(list.use)
list.df$L1 <- as.factor(list.df$L1)
list.df$L1 <- reorder(list.df$L1,-list.df$value)
pdf('suppfig2/numberdaneurons.pdf',width = 8,height = 1.29*8,useDingbats = F)
ggplot(list.df, aes(x = L1, y = value,fill = L1)) + 
  geom_bar(stat = 'identity') + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.position = 'none') + xlab('Study') + 
  ylab('Total number of DA neurons profiled') + scale_fill_brewer(palette = "Dark2")
dev.off()



# Supp figure sox6/calb1 ratios
sn.da.ctrl@cell.data$sox6 <- rep('sox6',nrow(sn.da.ctrl@cell.data))
sn.da.ctrl@cell.data[which(grepl('CALB1',sn.da.ctrl@clusters)),]$sox6 <- 'calb1'
t1 <- table(sn.da.ctrl@cell.data$dataset,sn.da.ctrl@clusters)
t2 <- as.data.frame(melt(t1/rowSums(t1)))
t2$Var1 <- as.factor(t2$Var1)
ggplot(t2, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="fill", stat="identity")


sn.da.annot <- qread('/home/tkamath/DA/da/sn_da_annot_021821.qs')
idx.ctrl <- rownames(sn.da.annot@cell.data[which(sn.da.annot@cell.data$status == 'Ctrl'),])
sn.da.ctrl <- subsetLiger(sn.da.annot,cells.use = idx.ctrl,remove.missing = F)
p1 <-plotByDatasetAndCluster(sn.da.ctrl,return.plots = T,do.shuffle = T,pt.size = 0.5)
p1[[1]]$labels$colour <- 'Donor ID'
setwd('~/DA/figures/suppfigures/suppfig2/')
pdf('tsne_da_byindiv.pdf',useDingbats = F)
p1[[1]]+ geom_point_rast() + theme(axis.text = element_blank(),
                                   axis.title = element_blank(),
                                   axis.ticks = element_blank(),
                                   axis.line = element_blank()) + theme(legend.position = c(0.8, 0.8))
dev.off()
pdf_convert('tsne_da_byindiv.pdf',format = 'png',dpi = 500)

# Plot markers 
marker.genes <- c('PBX1','PITX3','BNC2','SLC6A3','KCNJ6','LMO3','ALDH1A1')
p.use <- lapply(marker.genes,function(x){
  p1 <- VlnPlot(sn.da.ctrl.seurat,features.plot = x,x.lab.rot = T,point.size.use = 0,do.return = T)
  p1$data$ident <- factor(p1$data$ident,levels = c('SOX6_AGTR1','SOX6_PART1','SOX6_DDT','SOX6_GFRA2',
                                                   'CALB1_CALCR','CALB1_CRYM_CCDC68','CALB1_PPP1R17',
                                                   'CALB1_RBP4','CALB1_GEM','CALB1_TRHR'))
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)
  p1 + scale_fill_manual(values = mycolors) + theme(axis.title.x = element_blank(),axis.title.y =
                                                      element_text(color = 'black',face = 'italic')) +
    ylab(as.character(x)) +ggtitle('')
})
pdf('~/DA/figures/suppfigures/suppfig2/violinplots_DAothermarkers.pdf',
    useDingbats = F,width = 25,height = 15)
do.call('grid.arrange',p.use)
dev.off()


p1 <- rliger::plotGene(sn.da.ctrl,gene = 'TH',return.plots = T,pt.size = 1,
                       set.dr.lims = T,points.only = T,axis.labels = c('UMAP1','UMAP2'))
p1 <- lapply(p1,function(x){x + theme(plot.title = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.line = element_blank(),legend.title = element_blank())})
library(gridExtra)
png('markergenes_dath.png',res = 300,width = 8000,height = 1000)
grid.arrange(p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[1]],nrow = 1)
dev.off()

p1 <- rliger::plotGene(sn.da.ctrl,gene = 'SLC6A3',return.plots = T,pt.size = 1,
                       set.dr.lims = T,points.only = T,axis.labels = c('UMAP1','UMAP2'))
p1 <- lapply(p1,function(x){x + theme(plot.title = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.line = element_blank(),legend.title = element_blank())})
library(gridExtra)
png('markergenes_da6a3.png',res = 300,width = 8000,height = 1000)
grid.arrange(p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[1]],nrow = 1)
dev.off()

p1 <- rliger::plotGene(sn.da.ctrl,gene = 'SLC18A2',return.plots = T,pt.size = 1,
                       set.dr.lims = T,points.only = T,axis.labels = c('UMAP1','UMAP2'))
p1 <- lapply(p1,function(x){x + theme(plot.title = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.line = element_blank(),legend.title = element_blank())})
png('markergenes_da18a2.png',res = 300,width = 8000,height = 1000)
grid.arrange(p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[1]],nrow = 1)
dev.off()

p1 <- rliger::plotGene(sn.da.ctrl,gene = 'SOX6',return.plots = T,pt.size = 1,
                       set.dr.lims = T,points.only = T,axis.labels = c('UMAP1','UMAP2'))
p1 <- lapply(p1,function(x){x + theme(plot.title = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.line = element_blank(),legend.title = element_blank())})
png('markergenes_dasox6.png',res = 300,width = 8000,height = 1000)
grid.arrange(p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[1]],nrow = 1)
dev.off()

p1 <- rliger::plotGene(sn.da.ctrl,gene = 'CALB1',return.plots = T,pt.size = 1,
                       set.dr.lims = T,points.only = T,axis.labels = c('UMAP1','UMAP2'))
p1 <- lapply(p1,function(x){x + theme(plot.title = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.line = element_blank(),legend.title = element_blank())})
png('markergenes_dacalb1.png',res = 300,width = 8000,height = 1000)
grid.arrange(p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[1]],nrow = 1)
dev.off()

sn.da.ctrl@cell.data$sox6 <- rep('sox6',nrow(sn.da.ctrl@cell.data))
sn.da.ctrl@cell.data[which(grepl('CALB1',sn.da.ctrl@clusters)),]$sox6 <- 'calb1'
t1 <- table(sn.da.ctrl@cell.data$dataset,sn.da.ctrl@cell.data$sox6)
t2 <- as.data.frame(melt(t1/rowSums(t1)))
t2$Var1 <- as.factor(t2$Var1)
levels(t2$Var2) <- c('CALB1+','SOX6+')
pdf('suppfig2/sox6proportions.pdf',useDingbats = F)
ggplot(t2, aes(fill=Var2, y=value, x=Var2)) + 
  geom_boxplot(alpha = 0.2) + geom_jitter(width = 0.2) +
  theme(legend.position = 'none') + xlab('Broad DA subtype') + ylab('Proportion per individual') 
dev.off()

Hs_norm <- sn.da.ctrl@H.norm
W_mat <- sn.da.ctrl@W
factors.use <- list('CALB1/TRHR' = 1,
                    'SOX6/DDT' = 2,
                    'CALB1/RBP4' = 6,
                    'CALB1/PPP1R17' = 8,
                    'CALB1/CRYM/CCDC68' = 9,
                    'SOX6/PART1' = 10,
                    'SOX6/GFRA2' = 12,
                    'CALB1/GEM' = 15,
                    'SOX6/AGTR1' = 16,
                    'CALB1/CALCR' = 19)
factors.use <- factors.use[c(9,2,6,7,5,10,3,4,1,8)]
tsne.colors <- c('lemonchiffon', 'red')
p.list <- lapply(names(factors.use),function(z){
  idx.use <- which(names(factors.use) == z)
  tsne_df <- data.frame(Hs_norm[, factors.use[[z]] ], sn.da.ctrl@tsne.coords)
  factorlab <- paste0("Factor", idx.use)
  colnames(tsne_df) <- c(factorlab, "UMAP1", "UMAP2")
  genes.plot <- sn.da.ctrl@var.genes[head(order(W_mat[factors.use[[z]],],decreasing = T),5)]
  lab.use <- paste0('Factor', idx.use,': ',z,' -- \n',paste(genes.plot[1:2],collapse = ', '),'\n',
                    paste(genes.plot[3:5],collapse = ', '))
  #lab.use <- paste0(lab.use,paste(genes.plot[4:7],collapse = ', '),'\n')
  #lab.use <- paste0(lab.use,paste(genes.plot[8:10],collapse = ', '))
  p1 <- ggplot(tsne_df, aes_string(x = "UMAP1", y = "UMAP2", color = factorlab)) + 
    geom_point_rast(size = 0.8) +
    scale_color_gradient(low = tsne.colors[1], high = tsne.colors[2]) + 
    ggtitle(label = lab.use)+ 
    guides(color=guide_legend(title="Factor loading")) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),plot.title = element_text(hjust = 0.3,size = 20),
          legend.title = element_blank())
  return(p1)
})

pdf('/home/tkamath/DA/figures/suppfigures/suppfig2/inmffactors.pdf',useDingbats = F,
    width = 10*128/50,height = 10)
grid.arrange(p.list[[1]],p.list[[2]],p.list[[3]],p.list[[4]],p.list[[5]],
             p.list[[6]],p.list[[7]],p.list[[8]],p.list[[9]],p.list[[10]],nrow = 2)
dev.off()

pdf_convert('inmffactors.pdf',format = 'png',dpi = 500)
