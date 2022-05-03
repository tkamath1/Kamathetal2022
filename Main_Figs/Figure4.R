## Figure 3 plots
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

#########################################
################ Fig 4 A/B ##############
#########################################

# Known PD mendelian genes
# first plot heatmap
out.mast.rare <- qread('/home/tkamath/DA/gwas/rarevariants_mast.qs')
out.mast.rare$Var1 <- reorder(out.mast.rare$Var1,new.order = c(20:26,1:19))
out.mast.rare$Var2 <- factor(out.mast.rare$Var2,levels = c('Astrocyte','OPC','Oligodendrocyte',
                                                           'Endothelial cell/pericyte',
                                                           'Microglia/macrophage',
                                                           'Inhibitory neuron',
                                                           'Excitatory neuron',
                                                           'DA neuron'))
pdf('/home/tkamath/DA/figures/fig4/rarevariantsgwas.pdf',width = 9, height = 5.5,useDingbats = F)
ggplot(out.mast.rare, 
       aes(x = Var2,y = Var1,fill = value)) + geom_tile(color = 'black')+
  scale_fill_distiller(palette= 'RdYlBu',limits = c(0,1)*max(abs(out.mast.rare$value)))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11),axis.text.y = element_text(size = 13)) +
  scale_shape_identity() + coord_flip()
dev.off()

rare.pvals <- qread('/home/tkamath/DA/gwas/rarevariant_enrichmentpvalues.qs')
rare.pvals$celltype <- factor(rare.pvals$celltype,levels = c('astro','opc','olig','endo',
                                                             'Inh','Ex','da','mg'))
levels(rare.pvals$celltype) <-  c('Astrocyte','OPC','Oligodendrocyte','Endothelial cell',
                                  'Inhibitory neuron','Excitatory neuron','DA neuron',
                                  'Microglia/macrophage')
colnames(rare.pvals) <- c('celltype','pval','odds','intersection')
pdf('/home/tkamath/DA/figures/fig4/rarevariants_pvalues.pdf',width = 3, height = 6.5)
ggplot(rare.pvals, aes(x = celltype, y = -log10(pval))) + geom_bar(stat = 'identity',width = 0.75) + coord_flip() +
  gghighlight(pval < 0.05/8) +
  theme(axis.text.x = element_text(size =16),axis.text.y = element_blank(),
        axis.title.x = element_text(size = 20)) + ylab(expression(~-log[10]~(p-value))) +
  geom_hline(yintercept = -log10(0.05/8),linetype = 'dashed',color = 'red')
dev.off()


#######################################
################ Fig 4C ##############
#######################################

celltypes <- c('astro','da','endo','Ex','Inh','mg','olig','opc')
out.pvals <- lapply(celltypes,function(x){
  pd.use <-read.table(paste0('/home/tkamath/magma/pd',as.character(x),'mast.gsa.out'),header = T)
  ad.use <-read.table(paste0('/home/tkamath/magma/ad',as.character(x),'mast.gsa.out'),header = T)
  out.pvals <- do.call(rbind,list('PD' = pd.use$P,'AD' = ad.use$P))
  return(out.pvals)})
out.pvals <- do.call(rbind,out.pvals)
df.use <- data.frame('magmap' = out.pvals,'disease' = rownames(out.pvals))
df.use$celltype <- rep(celltypes,each = 4)
df.use <- df.use %>% mutate(magmapbin = ifelse(magmap < 0.05/8, 0.05/(8), magmap))
df.use$celltype <- as.factor(df.use$celltype)
levels(df.use$celltype) <- c('Astrocyte','DA neuron','Endothelial cell','Excitatory neuron',
                             'Inhibitory neuron','Microglia/macrophage','Oligodendrocyte','OPC')
df.use$disease <- factor(df.use$disease,levels = c('AD','PD','ADHD','SCZ'))
df.use$celltype <- factor(df.use$celltype,levels = rev(c('DA neuron','Excitatory neuron','Inhibitory neuron',
                                                         'Microglia/macrophage','Endothelial cell',
                                                         'Oligodendrocyte','OPC','Astrocyte')))

striatum.df <- qread('/home/tkamath/DA/striatum/gwas_analysis.qs')
colnames(striatum.df) <- c('celltype','disease','magmap')
striatum.df$disease <- as.factor(striatum.df$disease)
levels(striatum.df$disease) <- c('SCZ','AD','PD','ADHD')
striatum.df$disease <- factor(striatum.df$disease,levels = c('AD','PD','ADHD','SCZ'))
striatum.df <- striatum.df %>% mutate(magmapbin = ifelse(magmap < 0.05/8, 0.05/(8), magmap))
striatum.df$region <- rep('striatum',nrow(striatum.df))
df.use$region <- rep('sn',nrow(df.use))
levels(striatum.df$celltype) <- c('Astrocyte','dSPN','iSPN','Endothelial cell','Inhibitory neuron',
                                  'Microglia/macrophage','Oligodendrocyte','OPC')
striatum.df$celltype <- factor(striatum.df$celltype,
                               levels = c('iSPN','dSPN','Inhibitory neuron','Microglia/macrophage','Endothelial cell',
                                          'Oligodendrocyte','OPC','Astrocyte'))
df.use2 <- rbind(df.use,striatum.df)

pdf('/home/tkamath/DA/figures/fig4/magmanerichment_alldiseases2.pdf',width = 7,height = 5,useDingbats = F)
ggplot(df.use2[df.use2$disease %in% c('AD','PD'),]) + 
  geom_bar(stat = 'identity',aes(x = -log10(magmap), y = celltype, fill = celltype)) +
  scale_fill_manual(values = c("#E6AB02","#1B9E77","#7570B3","#E7298A","#66A61E")) +
  facet_wrap(~region + disease,
             scales = 'free_y') + gghighlight::gghighlight(magmap < 0.05/8,
                                                           calculate_per_facet = T)

dev.off()


#######################################
################ Fig 4D ##############
#######################################

subtypes.use <- qread('/home/tkamath/DA/gwas/genesets/subtypes/subtypes_test.qs')

magmap <- as.vector(unlist(lapply(subtypes.use,function(x){
  if (file.exists(paste0('/home/tkamath/magma/pd',as.character(x),'mast.gsa.out'))){
    magma.use <- read.table(paste0('/home/tkamath/magma/pd',as.character(x),'mast.gsa.out'),header = T)
    return(-log10(magma.use$P))
  }
  else
    return(0)
})))
df.use <- data.frame('celltype' = subtypes.use, 'magmap' = magmap)

df.use$celltype <- reorder(df.use$celltype, df.use$magmap)
df.use$cluster <- as.character(unlist(sapply(strsplit(as.character(df.use$celltype),'_'),`[`,1)))
df.use[which(grepl('SOX6|CALB1',df.use$cluster)),]$cluster <- 'DA'
mycolors <- colorRampPalette(brewer.pal(n = 8, "Dark2"))(10)
df.use[which(df.use$cluster == 'Ependyma'),]$cluster <- 'Astro'
df.use[which(df.use$cluster == 'Macro'),]$cluster <- 'MG'
df.use$cluster <- as.factor(df.use$cluster)
df.use$cluster <- factor(df.use$cluster,levels = rev(c('DA','Ex','Inh','MG', 'Endo','Olig','OPC','Astro')))
levels(df.use$cluster) <- rev(c('DA neuron','Excitatory neuron','Inhibitory neuron',
                                'Microglia/macrophage','Endothelial cell','Oligodendrocyte',
                                'OPC','Astrocyte'))

pdf('/home/tkamath/DA/figures/fig4/PDmagmaenrichment_mastsubtypesmarkers.pdf',width =5.5,height = 5,useDingbats = F)
ggplot() + geom_jitter(data = df.use, aes(x = cluster, y = magmap,fill = cluster,shape = 21),size = 7,stroke = 0.5,width = 0.2) +
  geom_hline(yintercept=-log10(0.05/68), linetype="dashed", color = "red") + 
  scale_fill_manual(values = rev(c("#1B9E77","#7570B3","#E7298A","#E6AB02","#66A61E",mycolors[7],
                                   mycolors[8],mycolors[9]))) + scale_shape_identity() + 
  geom_label_repel(data = df.use[df.use$magmap > -log10(0.05/68),],aes(x = cluster, y = magmap, 
                                                                       label = celltype)) +
  xlab('Celltype') + ylab(expression(~-log[10]~(p-value))) + theme(axis.text.x = 
                                                                     element_text(angle = 45,hjust = 1,size = 18),
                                                                   legend.position = 'none') +
  coord_flip() + ylim(0,NA)
dev.off()

##########################
## Fig 4E/F###############
##########################
magma.out.agtr1 <- qread('~/DA/gwas/gwas_pseudomanhattanmainfig.qs')

pdf('/home/tkamath/DA/figures/fig4/dotplot_zscores_agtr1_rd2.pdf',useDingbats = F,width = 12)
ggplot(magma.out.agtr1,aes(x = ZORDER2, y = ZSTAT,fill = spec,size = -log10(P),shape = shape  )) + 
  geom_point(aes(alpha = ZSTAT)) +
  geom_label_repel(data = magma.out.agtr1[magma.out.agtr1$ZSTAT > 4.568,],
                   aes(x = ZORDER2, y = ZSTAT,label = SYMBOL),fontface = 'italic') + 
  scale_fill_manual(values =c('grey', 'yellow','orange','red')) + scale_shape_manual(values = c(21,23)) +
  xlab('Scaled -log10(p-value)') + ylab('MAGMA Z-score')
dev.off()

magma.out.agtr1$spec <- NULL
magma.out.agtr1 <- magma.out.agtr1 %>% mutate(spec = ifelse(SYMBOL %notin% genes.nonagtr1,2,1))
magma.out.agtr1 <- magma.out.agtr1 %>% mutate(spec = ifelse(SYMBOL %in% genes.nonsox6,0,spec))

genes.use <- magma.out.agtr1[which(magma.out.agtr1$spec %in% c(1,2) & magma.out.agtr1$ZSTAT > 4.568  ),]$SYMBOL

library(enrichR)
dbs <- listEnrichrDbs()$libraryName
out.enrich <- enrichR::enrichr(genes = genes.use, databases = grep('GO_Biological_Process_2017',dbs,value = T))
out.enrich.use <- bind_rows(out.enrich,.id = 'names')

out.enrich.use2 <- out.enrich.use[order(out.enrich.use$Adjusted.P.value,decreasing = F)[c(1:10)],]
out.enrich.use2$Term <- as.factor(out.enrich.use2$Term)
out.enrich.use2$Term <-reorder(out.enrich.use2$Term,-log10(out.enrich.use2$P.value))
pdf("/home/tkamath/DA/gwas/gwas_goterms.pdf",useDingbats = F,height = 14,width = 14)
ggplot(out.enrich.use2,
       aes(x = Term,y = -log10(P.value))) + geom_bar(stat = 'identity',fill = 'grey',color = 'black') + coord_flip()
dev.off()

############################
## Fig 4G/H ################
############################
out.results.use <- qread('~/DA/diffexp/outresults_TFDE.qs')

out.results.use2 <- out.results.use[setdiff(names(out.results.use),
                                            names(which(sapply(out.results.use,nrow) == 0)))]
out.results.use4 <- bind_rows(out.results.use2,.id = 'names2')

df.out <- lapply(unique(out.results.use4$pathway),function(x){
  if (nrow(out.results.use4[out.results.use4$pathway == x & 
                            out.results.use4$names2 != 'SOX6_AGTR1_nonurr',]) > 0 &
      nrow(out.results.use4[out.results.use4$pathway == x & 
                            out.results.use4$names2 == 'SOX6_AGTR1_nonurr',])){
    a1 <- min(out.results.use4[out.results.use4$pathway == x & 
                                 out.results.use4$names2 != 'SOX6_AGTR1_nonurr',]$padj)
    a2 <- min(out.results.use4[out.results.use4$pathway == x & 
                                 out.results.use4$names2 == 'SOX6_AGTR1_nonurr',]$padj)
    return(data.frame('pvalout' = a1,'pvalin' = a2))
  }
  else
    return(NULL)
})
names(df.out) <- unique(out.results.use4$pathway)
df.out <- bind_rows(df.out,.id = 'names')

out.results <- mclapply(df.test$names2,function(z){
  out.plots <- lapply(unique(names(mast.da)),function(x){
    gene.set <- mast.da[[x]][which(mast.da[[x]]$fdr <= 0.99),]
    gene.set.use <- (sign(gene.set[order(gene.set$fdr,decreasing = F),]$coef))*(-log10(gene.set[
      order(gene.set$fdr,decreasing = F),]$`Pr(>Chisq)`))
    names(gene.set.use) <- gene.set[order(gene.set$fdr,decreasing = F),]$primerid
    fgseaRes <- fgsea(pathways = list('p53' = unique(genes.db2[[z]] ) ),stats = gene.set.use,
                      minSize  = 1,maxSize  = 500,nperm = 1000)
    p1 <- plotEnrichment(pathway = genes.db2[[z]] ,stats = gene.set.use)
    return(p1$data)
  })
  names(out.plots)<- unique(names(mast.da))
  bind_rows(out.plots,.id = 'names')
},mc.cores = 20)
names(out.results) <- df.test$names2
out.results.plot <- bind_rows(out.results,.id = 'names2')

masc.results <- qread('/home/tkamath/DA/ratio_analysis/mascresults_simplemodel.qs')
masc.da <- masc.results[masc.results$term == 'diseaseDisease' & masc.results$name == 'da',]
masc.da$cluster <- substring(masc.da$cluster,first = 8)
masc.da$cluster <- as.vector(unlist(lapply(strsplit(masc.da$cluster,'_'),
                                           function(x){paste0(x[[1]],'_',x[[2]])})))

out.results.plot$names <- as.vector(unlist(lapply(strsplit(as.character(out.results.plot$names),
                                                           '_'),function(x){paste0(x[[1]],'_',x[[2]])})))
df.test$names <- as.vector(unlist(lapply(strsplit(as.character(df.test$names),
                                                  '_'),function(x){paste0(x[[1]],'_',x[[2]])})))
out.results.plot$oddsratio <- masc.da[match(out.results.plot$names,masc.da$cluster),]$estimate
masc.da$gsea <- df.test[match(masc.da$cluster,df.test$names),]$ES
values.use <- colorRampPalette(brewer.pal(11,"RdBu"))(12)
out.results.plot$names <- reorder(out.results.plot$names,-out.results.plot$oddsratio)

qsave(out.results.plot,'~/DA/tfanalysis/tf_regulons_toplot.qs')

values.use <- colorRampPalette(brewer.pal(11,"RdBu"))(12)
pdf('~/DA/tfanalysis/tf_deregulons_supp.pdf',useDingbats = F,width = 14)
ggplot(out.results.plot[out.results.plot$names2 %notin% c('LMX1A','TP53','NR2F2'),], 
       aes(x = x, y = y, color = names)) + geom_point(size = 0.2) + 
  facet_wrap(facets=~names2,scales = 'free_x',nrow = 2) + scale_color_manual(values = values.use[c(2:11)]) +
  geom_hline(yintercept = 0)
dev.off()




