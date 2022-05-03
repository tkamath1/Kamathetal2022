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
################# Supp Figure 8  #####################
######################################################


### Panel A ########
t2 <- qread('/home/tkamath/DA/ratio_analysis/ratios_allcells.qs')
pdf('/home/tkamath/DA/figures/suppfigures/suppfig5/normalized_ratio.pdf',useDingbats = F,
    width = 10,height = 5)
ggplot(t2,aes(x = Var2, y = log2(z),fill = status2 )) + geom_jitter(position=position_jitterdodge(),
                                                                    aes(group=status2),color = 'black',size = 1) + 
  geom_boxplot(aes(fill = status2),alpha = 0.2,outlier.size = NULL,outlier.color = 'grey') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab('Cell type') + ylab('log2(Normalized ratio)')
dev.off()

### Panel B ########
sn.mg.seurat <- qread('/home/tkamath/DA/SN2020_selected/SN_mg_indScaling/mg_full_annotated.qs')
sn.mg.seurat@raw.data <- sn.mg.seurat@raw.data[,rownames(sn.mg.seurat@meta.data)]
sn.mg.seurat@dr$cca <- sn.mg.seurat@dr$hpca
FeaturePlot(sn.mg.seurat,reduction.use = 'umap',features.plot = '')

sn.mg <- rliger::seuratToLiger(sn.mg.seurat,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
sn.mg@tsne.coords <- sn.mg.seurat@dr$umap@cell.embeddings
plotByDatasetAndCluster(sn.mg)

sn.mg <- rliger::normalize(sn.mg)
sn.mg <- scaleNotCenter(sn.mg)
sn.mg@cell.data <- sn.mg@cell.data[names(getGeneValues(sn.mg@raw.data,gene = 'GPNMB')),]
sn.mg@tsne.coords <- sn.mg@tsne.coords[names(getGeneValues(sn.mg@raw.data,gene = 'GPNMB')),]

pdf('/home/tkamath/DA/figures/suppfigures/suppfig5/plotgene_mg_GPNMB.pdf',useDingbats = F)
rliger::plotGene(sn.mg,gene = 'GPNMB',plot.by = 'none') + geom_point_rast()
dev.off()
pdf('/home/tkamath/DA/figures/suppfigures/suppfig5/plotgene_mg_SULT1C2.pdf',useDingbats = F)
rliger::plotGene(sn.mg,gene = 'SULT1C2',plot.by = 'none') + geom_point_rast()
dev.off()

sn.astro.seurat <- qread('/home/tkamath/DA/SN2020_selected/SN_astro_indScaling/astro_full_annotated.qs')
sn.astro.seurat@raw.data <- sn.astro.seurat@raw.data[,rownames(sn.astro.seurat@meta.data)]
sn.astro.seurat@dr$cca <- sn.astro.seurat@dr$hpca
sn.astro <- rliger::seuratToLiger(sn.astro.seurat,combined.seurat = T,meta.var = 'subject',cca.to.H = T)
sn.astro@tsne.coords <- sn.astro.seurat@dr$umap@cell.embeddings
plotByDatasetAndCluster(sn.astro)
sn.astro <- rliger::normalize(sn.astro)
sn.astro@cell.data <- sn.astro@cell.data[names(getGeneValues(sn.astro@raw.data,gene = 'VIM')),]
sn.astro@tsne.coords <- sn.astro@tsne.coords[names(getGeneValues(sn.astro@raw.data,gene = 'VIM')),]

pdf('/home/tkamath/DA/figures/suppfigures/suppfig5/plotgene_astro_VIM.pdf',useDingbats = F,
    width = 7)
rliger::plotGene(sn.astro,gene = 'VIM',plot.by = 'none') + geom_point_rast()
dev.off()
pdf('/home/tkamath/DA/figures/suppfigures/suppfig5/plotgene_astro_LHX2.pdf',useDingbats = F)
rliger::plotGene(sn.astro,gene = 'LHX2',plot.by = 'none') + geom_point_rast()
dev.off()

### Panel C ########
# Downsampling analysis for astro and microglia
out.vals2 <- qread('/home/tkamath/DA/ratio_analysis/downsample_analysis_astro.qs')
pdf('~/DA/figures/suppfigures/suppfig7/downsampling_astro.pdf',useDingbats = F,width = 4*57/26,height = 4)
ggplot(out.vals2, aes(x = id, y = -log10(value), fill = variable,shape=21)) + geom_point(size = 4) +
  scale_shape_identity() +  guides(fill=guide_legend(override.aes=list(shape=21)))+
  geom_hline(yintercept = 1.3,col="red", lwd=0.8, lty=2) + scale_x_log10() + ylim(0,NA) +
  theme(legend.position = c(0.6,0.3))
dev.off()

out.vals2 <- qread('/home/tkamath/DA/ratio_analysis/downsample_analysis_mg.qs')
pdf('~/DA/figures/suppfigures/suppfig7/downsampling_mg.pdf',useDingbats = F,width = 4*57/26,height = 4)
ggplot(out.vals2, aes(x = id, y = -log10(value), fill = variable,shape=21)) + geom_point(size = 4) +
  scale_shape_identity() +  guides(fill=guide_legend(override.aes=list(shape=21)))+
  geom_hline(yintercept = 1.3,col="red", lwd=0.8, lty=2) + ylim(0,NA) +
  theme(legend.position = c(0.1,0.9))
dev.off()
