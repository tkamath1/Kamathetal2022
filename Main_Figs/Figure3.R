#######################
## Figure 3 plots #####
#######################
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


#############################################
##### Figure 3A #############################
#############################################

pd.meta.use <- qread('/home/tkamath/DA/pdmeta_use.qs')

out.masc <- qread('/home/tkamath/DA/ratio_analysis/mascresults_simplemodel.qs')
out.masc <- bind_rows(out.masc, .id = 'names')
out.ratios.use <- out.masc[out.masc$term =='diseaseDisease', ]
out.ratios.use$color <- rep('black',nrow(out.ratios.use))
out.ratios.use$fdr <- p.adjust(out.ratios.use$p.value,method = 'fdr')
out.ratios.use <- out.ratios.use %>% mutate(color = ifelse(fdr < 0.05, 'red','black'))

mycols <- brewer.pal(n = 8, 'Dark2')
out.ratios.use$name <- as.factor(out.ratios.use$name)
out.ratios.use$name <- factor(out.ratios.use$name,levels = c('da','olig','opc','astro','endo','mg','nonda'))
out.ratios.use$cluster <- substring(out.ratios.use$cluster,first = 8)
#levels(out.ratios.use$name) <- c('DA neuron','Oligodendrocyte','OPC','Astrocyte','Endothelial cell/fibroblast',
#                                 'Microglia/macrophage','Ex/Inh Neuron')
out.ratios.use$estimate <- out.ratios.use$estimate/log(2)
out.ratios.use$std.error <- out.ratios.use$std.error/log(2)
pdf('/home/tkamath/DA/figures/volcanoplot_fig2.pdf',useDingbats = F,width = 7, height = 7)
ggplot(out.ratios.use,aes(x = estimate, y = -log10(fdr), color = name)) + 
  geom_point() + theme(legend.position = 'right') +
  scale_color_brewer(palette = 'Dark2') + xlab('log(Odds Ratio)') +
  ylab(expression(~-log[10]~(FDR))) + geom_label_repel(data = out.ratios.use[out.ratios.use$fdr < 0.05,],
                                                       aes(x = estimate, y = -log10(fdr),label = cluster),
                                                       label.size = 1.5,label.r = 0,size = 5) +
  theme(legend.position = 'none',axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))
dev.off()


############################################
############### Fig 3B #####################
############################################

sn.da.annot <- qread('/home/tkamath/DA/da/sn_da_annot_021821.qs')
idx.ctrl <- rownames(sn.da.annot@cell.data[which(sn.da.annot@cell.data$status == 'Ctrl'),])
sn.da.ctrl <- subsetLiger(sn.da.annot, cells.use = idx.ctrl, remove.missing = F)
p1 <- plotByDatasetAndCluster(sn.da.ctrl,return.plots = T)
data.use <- p1[[2]]$data
data.text <- p1[[2]]$layers[[2]]$data
p3 <- ggplot(data.use, aes(x=tsne1, y=tsne2) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F) +
  scale_fill_viridis_c(option = "magma")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  )

idx.disease <- rownames(sn.da.annot@cell.data[which(sn.da.annot@cell.data$status == 'Disease'),])
sn.da.disease <- subsetLiger(sn.da.annot, cells.use = idx.disease, remove.missing = F)
p1 <- plotByDatasetAndCluster(sn.da.disease,return.plots = T,do.legend = F)
data.use2 <- p1[[2]]$data
data.text <- p1[[2]]$layers[[2]]$data
p4 <- ggplot(data.use2, aes(x=tsne1, y=tsne2) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F) +
  scale_fill_viridis_c(option = "magma")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  )

sn.celldata <- sn.da.annot@cell.data
sn.celldata$tsne1 <- sn.da.annot@tsne.coords[,1]
sn.celldata$tsne2 <- sn.da.annot@tsne.coords[,2]

t1 <- table(sn.celldata$status)

tsne_df <- data.frame('tsne1' = sn.celldata$tsne1, 'tsne2' = sn.celldata$tsne2)
rownames(tsne_df) <- rownames(sn.celldata)

p1 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2)) +
  geom_point(size = 0.3) +
  guides(color = guide_legend(override.aes = list(size = 2)))+ geom_bin2d(bins = 40) +
  scale_fill_continuous(type = "viridis") + ggtitle('iNPH UMAP density')
p1
ggp1 <- ggplot_build(p1)
densitydf <- ggp1$data[[2]]

lis.bx <- lapply(c(1:dim(densitydf)[1]),function(x){
  rownames(tsne_df)[which( (tsne_df$tsne1 < densitydf$xmax[x] & tsne_df$tsne1 > densitydf$xmin[x]) & 
                             (tsne_df$tsne2 < densitydf$ymax[x] & tsne_df$tsne2 > densitydf$ymin[x]) ) ]
})

ratios.density <- do.call(rbind,lapply(lis.bx, function(x){
  t1 <- table(sn.celldata[x,]$dataset)
  t2 <- t1/sum(t1)
  return(t1)
}))

m.density <- melt(scale(ratios.density))

pd.meta.use <- qread('/home/tkamath/DA/pdmeta_use.qs')
m.density$status <- pd.meta.use[match(m.density$Var2,pd.meta.use$Donor.ID),]$disease
m2 <- aggregate(value ~ Var1 + status, m.density,mean)

mdisease <- m2[which(m2$status == 'Disease'),]
mctrl <- m2[which(m2$status == 'Ctrl'),]


mid =0
df.plotcommon <- data.frame('x' = densitydf$x, 'y' = densitydf$y,
                            'diseasez' =  scale(mdisease$value - mctrl$value ))
sc <- scale_fill_gradient2(low= 'slateblue4', high = 'red',mid = 'white')
p1 <- ggplot(df.plotcommon, aes(x=x, y=y, fill = diseasez)) + 
  geom_point(size = 3.5,shape = 21,color = 'white',stroke = 0.05) + sc + scale_shape_manual()
pdf('/home/tkamath/DA/figures/densityplot_Daneurons.pdf',useDingbats = F,height = 5, width = 6)
p1
dev.off()

p1 <- plotByDatasetAndCluster(sn.da.annot,return.plots = T,do.legend = F,text.size = 2)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)
pdf('tsne_Daneurons.pdf',useDingbats = F,height = 3.5, width = 3.5)
p1[[2]] + scale_color_manual(values = mycolors)
dev.off()

out.ratios.da <- out.ratios.use[out.ratios.use$name == 'da',]
out.ratios.da$cluster <- reorder(out.ratios.da$cluster,out.ratios.da$estimate)


############################################
############### Fig 3C #####################
############################################

mycolors <- colorRampPalette(brewer.pal(10, "RdBu"))(13)
pdf('oddsplot.pdf',useDingbats = F,height = 6,width = 5)
ggplot(out.ratios.da,aes(x=cluster,ymin=estimate-2.5*std.error, ymax=estimate+2.5*std.error)) +
  geom_errorbar(inherit.aes = T,width=0.4, alpha=0.9, size=0.3, aes(colour = cluster) ) + 
  geom_point(aes(x = cluster, y = estimate),fill = 'white',size = 3,shape = 21,stroke = 0.3,color = 'grey') +
  coord_flip() +
  theme(legend.position = 'none') + scale_color_manual(values = rev(mycolors)) + 
  geom_hline(yintercept = 0, linetype="dotted", 
             color = "grey", size=0.3)
dev.off()

