library(rliger)
library(moments)
library(Matrix)
library(qs)
library(MAST)
library(dplyr)
library(lme4)
library(org.Hs.eg.db)
library(NNLM)
library(devtools)
library(Rtsne)
source("/home/tkamath/scripts/genevalidatehelper.R")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/striatum/')
options(mc.cores = 60)

######################################################
################# Supp Figure 10  #####################
######################################################

##############################
#### Ex Data 10A/B ############
##############################

ad.cn <- read.table('/home/tkamath/ldscore/ldsc/AD_CNcelltypes.cell_type_results.txt',header = T)
ad.sn <- read.table('/home/tkamath/ldscore/ldsc/AD_SNcelltypes.cell_type_results.txt',header = T)

pd.cn <- read.table('/home/tkamath/ldscore/ldsc/PD_CNcelltypes.cell_type_results.txt',header = T)
pd.sn <- read.table('/home/tkamath/ldscore/ldsc/PD_maincelltypes.cell_type_results.txt',header = T)

df.ldsc <- data.frame('pval' = c(ad.cn$Coefficient_P_value,pd.cn$Coefficient_P_value,
                                 ad.sn$Coefficient_P_value,pd.sn$Coefficient_P_value),
                      'Cell_types' = c(as.vector(ad.cn$Name),as.vector(pd.cn$Name),
                                       as.vector(ad.sn$Name),as.vector(pd.sn$Name)),
                      'region' = c(rep('CN',length = 16),rep('SN',length = 16)),
                      'Disease' = c(rep('AD',length = 8),rep('PD',length = 8),
                                    rep('AD',length = 8),rep('PD',length = 8) )) 

# MAGMA major cell types
df.cn <- qread('/home/tkamath/DA/striatum/gwas_analysis.qs')
colnames(df.cn) <- c('Cell_types','Disease','pval')
df.cn <- df.cn[,c('pval', 'Cell_types','Disease')]
levels(df.cn$Disease) <- c('AD','PD')
df.sn <- qread('/home/tkamath/DA/gwas/snmagma_majorcelltypes_01.qs')
colnames(df.sn) <- c('Cell_types','Disease','pval')
levels(df.sn$Disease) <- c('AD','PD')

df.magma <- bind_rows(list('CN' = df.cn,'SN' = df.sn),.id = 'region')

df.heritability <- bind_rows(list('MAGMA' = df.magma, 'sLDSC' = df.ldsc),.id = 'method')
df.heritability$Cell_types <- as.factor(df.heritability$Cell_types)
levels(df.heritability$Cell_types) <- c('Astrocyte','Astrocyte','DA neuron',
                                        'dSPN','iSPN','Endothelial/pericyte','Endothelial/pericyte',
                                        'Excitatory neuron','Inhibitory neuron','Inhibitory neuron',
                                        'Microglia/macrophage','Oligodendrocyte','Oligodendrocyte',
                                        'OPC','OPC')
df.heritability$region <- as.factor(df.heritability$region)
df.heritability$region <- factor(df.heritability$region,levels = c('SN','CN'))
df.heritability$Cell_types <- factor(df.heritability$Cell_types,levels = rev(c('DA neuron','Excitatory neuron',
                                                                               'Inhibitory neuron','dSPN','iSPN',
                                                                               'Microglia/macrophage','Astrocyte','Endothelial/pericyte',
                                                                               'OPC','Oligodendrocyte')))
pdf('/home/tkamath/DA/figures/suppfigures/suppfig10/sldsc_supp.pdf',width = 10,useDingbats = F)
ggplot(df.heritability[df.heritability$method == 'sLDSC',],aes(x = Cell_types,-log10(pval))) + 
  geom_bar(stat = 'identity',aes(fill = Cell_types) ) +
  scale_fill_manual(values = c('#E6AB02', '#1B9E77',"#7570B3","#E7298A","#66A61E")) +
  facet_wrap(~region + Disease,scales = 'free_y') + coord_flip() + gghighlight::gghighlight(pval < 0.05/8,
                                                                                            calculate_per_facet = T) +
  ylim(NA,6) + ggtitle(label = 'sLDSC')
dev.off()

df.ldsc <- read.table('/home/tkamath/ldscore/ldsc/PD_subtypes.cell_type_results.txt',header = T)
df.ldsc$magmap <- -log10(df.ldsc$Coefficient_P_value)
df.ldsc <- df.ldsc[,c('Name','magmap')]
colnames(df.ldsc)[1] <- c('celltype')
df.ldsc$cluster <- as.character(unlist(sapply(strsplit(as.character(df.ldsc$celltype),'_'),`[`,1)))
df.ldsc[which(grepl('SOX6|CALB1',df.ldsc$cluster)),]$cluster <- 'DA'

df.ldsc$cluster <- as.factor(df.ldsc$cluster)

df.ldsc[which(df.ldsc$cluster == 'Ependyma'),]$cluster <- "Astro"
df.ldsc[which(df.ldsc$cluster == 'Macro'),]$cluster <- "MG"
df.ldsc$cluster <- as.factor(df.ldsc$cluster)
df.ldsc$cluster <- droplevels(df.ldsc$cluster)
levels(df.ldsc$cluster) <- c('Astrocyte','DA neuron','Endothelial cell/pericyte',
                             "Ex neuron",'Inh neuron','Microglia/macrophage','Oligodendrocyte','OPC')
df.ldsc$cluster <- factor(df.ldsc$cluster,levels = rev(c('DA neuron','Ex neuron','Inh neuron',
                                                         'Microglia/macrophage','Endothelial cell/pericyte','Oligodendrocyte','OPC' ,'Astrocyte')))
mycolors <- colorRampPalette(brewer.pal(n = 8, "Dark2"))(8)
pdf('/home/tkamath/DA/figures/suppfigures/suppfig10/ldsc_subtypes_supp.pdf',width =6.5,height = 5,useDingbats = F)
ggplot() + geom_jitter(data = df.ldsc, aes(x = cluster, y = magmap,fill = cluster,shape = 21),
                       size = 7,stroke = 0.5,width = 0.2) +
  geom_hline(yintercept=-log10(0.05/68), linetype="dashed", color = "red") + 
  scale_fill_manual(values = rev(mycolors) ) + scale_shape_identity() + 
  geom_label_repel(data = df.ldsc[df.ldsc$magmap > -log10(0.05/68),],aes(x = cluster, y = magmap, 
                                                                         label = celltype)) +
  xlab('Celltype') + ylab(expression(~-log[10]~(p-value))) + theme(axis.text.x = 
                                                                     element_text(angle = 45,hjust = 1,size = 18),
                                                                   legend.position = 'none') +
  coord_flip() + theme(panel.spacing = unit(2, "lines"))
dev.off()

##############################
#### Ex Data 10C #############
##############################

# P-value/Nuclei number
sn.da.annot <- qread('/home/tkamath/DA/da/sn_da_annot_021821.qs')
sn.da.ctrl <- subsetLiger(sn.da.annot, cells.use = rownames(sn.da.annot@cell.data[
  which(sn.da.annot@cell.data$status == "Ctrl"),]),remove.missing = F)
t1 <- table(sn.da.ctrl@clusters)

a3 <- lapply(names(t1),function(x){
  a1 <- read.table(paste0('/home/tkamath/magma/pd',x,'mast.gsa.out'),header = T)
  return(a1$NGENES)
})

a2 <- lapply(names(t1),function(x){
  a1 <- read.table(paste0('/home/tkamath/magma/pd',x,'mast.gsa.out'),header = T)
  return(a1$P)
})
names(a2) <- names(t1)
df.plot <- melt(t1)
df.plot$value2 <- a2
df.plot$value2 <- as.vector(unlist(df.plot$value2))
summary(lm(value2 ~ value, data = df.plot))

df.plot.use <- df.plot[df.plot$Var1 != 'SOX6_AGTR1',]

pdf('/home/tkamath/DA/figures/suppfigures/suppfig10/pvalue_nucleisize.pdf',useDingbats = F)
ggplot(data = df.plot.use,aes(x = value, y = -log10(value2),label = Var1 )) + 
  geom_point() + geom_smooth(method = 'lm',se = F) + ylim(0,NA) + geom_label_repel() +
  xlab('Number of nuclei sampled') + ylab('-log10(MAGMA p-value)') + 
  stat_fit_glance(method = 'lm',geom = 'text',label.x.npc = 0.9,
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 2), sep = "")))
dev.off()

pd.subtypes <- read.table('/home/tkamath/ldscore/ldsc/PD_subtypes.cell_type_results.txt',header = T)
pd.subtypes.use <- pd.subtypes[grep('CALB1.|SOX6',pd.subtypes$Name),]
pd.subtypes.use$value2 <- df.plot[match(pd.subtypes.use$Name,df.plot$Var1),]$value

pd.subtypes.use <- pd.subtypes.use[pd.subtypes.use$Name != 'SOX6_AGTR1',]

summary(lm(-log10(Coefficient_P_value) ~ value2, data = pd.subtypes.use))


pdf('/home/tkamath/DA/figures/suppfigures/suppfig10/pvalue_nucleisize_sldsc.pdf',useDingbats = F)
ggplot(data = pd.subtypes.use,aes(x = value2, y = -log10(Coefficient_P_value),label = Name )) + 
  geom_point() + geom_smooth(method = 'lm',se = F) + ylim(0,NA) + geom_label_repel() +
  xlab('Number of nuclei sampled') + ylab('-log10(s-LDSC p-value)') + 
  stat_fit_glance(method = 'lm',geom = 'text',label.x.npc = 0.9,
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 2), sep = "")))
dev.off()

##############################
#### Ex Data 10D #############
##############################
df.robustness <- qread('/home/tkamath/DA/gwas/dfplot_magma_robustness_differentdiseases.qs')
pdf('/home/tkamath/DA/figures/suppfigures/suppfig10/poweranalysis.pdf',useDingbats = F,width = 10)
ggplot(df.robustness,aes( x= genesize, y = -log10(value),fill = disease,shape = 21,group = disease)) + 
  geom_point(size = 4) +
  scale_shape_identity()  +geom_hline(yintercept = -log10(0.05/8),linetype = 'dashed',color = 'red') + xlim(1000,4000)
dev.off()

df.robustness2 <- qread('~/DA/gwas/df_plot_magma_robustness.qs')
df.robustness2 <- df.robustness2[df.robustness2$FDR_thresh == 0.05,]
df.robustness2 <- df.robustness2[grep('CALB1.|SOX6',df.robustness2$Cell_subtype),]
df.robustness2$Cell_subtype <- droplevels(df.robustness2$Cell_subtype)
df.robustness2$Cell_subtype <- factor(df.robustness2$Cell_subtype,levels = c('SOX6_AGTR1','SOX6_PART1','CALB1_CALCR',
                                                                             'SOX6_DDT','SOX6_GFRA2','CALB1_PPP1R17','CALB1_CRYM_CCDC68',
                                                                             'CALB1_GEM','CALB1_RBP4','CALB1_TRHR'))

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)
pdf('/home/tkamath/DA/figures/suppfigures/suppfig10/robustness_dasubtypes.pdf',width = 5*126/46,height = 5,useDingbats = F)
ggplot(df.robustness2,aes(x = Gene_set_size, y = -log10(`p-val`),fill = Cell_subtype,shape = 21 )) + 
  geom_jitter(size = 5,width = 0.2)+
  scale_shape_identity() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  scale_fill_manual(values = mycolors) + ylim(0.0000000001,NA) + guides(fill=guide_legend(override.aes=list(shape=21)))
dev.off()

##############################
#### Ex Data 10E #############
##############################
# Load in magma genes from Nalls et al PD GWAS
magma.out <- read.csv('/home/tkamath/magma/PDFUMAresults/magma.genes.out',header = T,sep = '\t')

files.use <- grep('mastDE_.*(CALB1|SOX6)',list.files('/home/tkamath/DA/gwas/genesets/subtypes/'),
                  value = T)
files.use <- files.use[c(c(1:6,9:12))]
out.mast <- lapply(files.use,function(x){
  qread(paste0('/home/tkamath/DA/gwas/genesets/subtypes/',x))
})

names(out.mast) <- as.vector(unlist(lapply(lapply(strsplit(files.use,'_'),
                                                  function(x){paste0(x[[2]],"_",x[[3]])}),function(x){
                                                    sapply(strsplit(x,'\\.'),`[`,1)})))
da.mast <- CombineMAST(out.mast,component.pvalue.use = 'H',component.coef.use = 'D',
                       correction = 'fdr',contrast.use = 'clusteruse1')
da.mast.use <- lapply(da.mast,function(x){x[x$fdr < 0.01 & x$coef > 0,]})
genes.mast.use <- lapply(da.mast.use,function(x){
  x[order(x$z, decreasing = T)[c(1:min(nrow(x),3500)) ],]$primerid
})


magma.df <- lapply(names(genes.mast.use),function(x){
  genes.agtr1 <- genes.mast.use[[x]]
  genes.nonagtr1 <- unique(as.vector(unlist(genes.mast.use[
    setdiff(c(1:10),grep(x,names(genes.mast.use) )) ])))
  genes.plot <- setdiff(genes.agtr1,genes.nonagtr1)
  
  magma.out$agtr1zscore <- -log10(da.mast[[x]][match(
    magma.out$SYMBOL,da.mast[[x]]$primerid),]$`Pr(>Chisq)`)
  magma.out.use <- magma.out[!is.na(magma.out$agtr1zscore),]
  
  magma.out.use$agtr1bin <- rep(0,length = nrow(magma.out.use))
  magma.out.use[magma.out.use$SYMBOL %in% genes.plot,]$agtr1bin <- 1
  
  magma.out.use$agtr1bin <- as.factor(magma.out.use$agtr1bin)
  magma.out.use <- magma.out.use[magma.out.use$agtr1bin == 1,]
  magma.out.use <- magma.out.use[order(magma.out.use$agtr1zscore,decreasing = T),]
  magma.out.use$orderz <- rev(c(1:nrow(magma.out.use)))
  magma.out.use$orderz <- magma.out.use$orderz/max(magma.out.use$orderz)
  return(magma.out.use)
})
names(magma.df)<- names(genes.mast.use)

magma.df <- bind_rows(magma.df,.id= 'celltypes')

magma.df$celltypes <- as.factor(magma.df$celltypes)

# Make sure to see that the LEVELS MATCH UP!!!
levels(magma.df$celltypes) <- c('CALB1_CALCR','CALB1_CRYM_CCDC68','CALB1_GEM',
                                'CALB1_PPP1R17','CALB1_RBP4','CALB1_TRHR',
                                'SOX6_AGTR1','SOX6_DDT','SOX6_GFRA2','SOX6_PART1')
magma.df$celltypes <- reorder(magma.df$celltypes,new.order = c('SOX6_AGTR1','SOX6_DDT',
                                                               "SOX6_PART1",'CALB1_CALCR',
                                                               'CALB1_CRYM_CCDC68','SOX6_GFRA2','CALB1_PPP1R17','CALB1_RBP4',
                                                               'CALB1_TRHR','CALB1_GEM'))

pdf('/home/tkamath/DA/gwas/dotplot_zscores.pdf',width = 24,height = 10,useDingbats = F)
ggplot(magma.df,aes(x = orderz, y = ZSTAT)) + 
  geom_point(alpha = 0.4,aes(size = ZSTAT)) + 
  geom_label_repel(data = magma.df[magma.df$ZSTAT > 4.568,],
                   aes(x = orderz, y = ZSTAT, label = SYMBOL)) + 
  facet_wrap(~celltypes,nrow = 2,scales = 'free_x') + xlab('Specific gene sets') +
  ylab('MAGMA Z scores')
dev.off()



