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

################# Supp Figure 1#######################
######################################################
library(DropSeq.util)
# Grab these from dropviz.org
dge.path <- '/home/tkamath/DA/dropviz/F_GRCm38.81.P60SubstantiaNigra.raw.dge.txt.gz'
dge <- loadSparseDge(dge.path)
cluster.sn <- readRDS('/home/tkamath/DA/dropviz/F_GRCm38.81.P60SubstantiaNigra.cluster.assign.RDS')
dropviz.sn <- CreateSeuratObject(dge)
dropviz.sn@ident <- cluster.sn
dropviz.sn <- SubsetData(dropviz.sn,cells.use = names(dropviz.sn@ident),subset.raw = T)
dropviz.sn <- SubsetData(dropviz.sn, cells.use = names(dropviz.sn@ident[dropviz.sn@ident %in% c(1:14)]),
                         subset.raw = T)
levels(dropviz.sn@ident) <- c('Neuron_CA3_C1ql3','Neuron_Rora','Neuron_Gad2','Neuron_Th', 
                              'Polydendrocyte_Tnr','Polydendrocyte_Cacng4','Astrocyte_Gja1','Ependyma',
                              'Microglia_Macrophage_C1qb','	Oligodendrocyte_Tfr','Oligodendrocyte_Mbp',
                              '	Endothelial_Flt1','Mural_Rgs5Acta2','	Fibroblast-Like_Dcn')
p1 <- DotPlot(dropviz.sn, genes.plot = c('Th','Nr4a2'),cols.use = c('orange','red'),
              x.lab.rot = T,plot.legend = T,do.return = T)
p1$data$id <- reorder(p1$data$id,p1$data$pct.exp)
pdf('suppfig1/nurr_mousebrain.pdf',width = 5, height = 5,useDingbats = F)
p1
dev.off()

# Make stacked violin plot instead
p1 <- StackedVlnPlot(dropviz.sn,features = c('Th','Nr4a2'))
pdf('suppfig1/nurr_mousebrain_violinplot.pdf',width = 5, height = 5,useDingbats = F)
p1 + theme(axis.text.x = element_text(angle = 90))
dev.off()

markers.th <- prestowrapper(dropviz.sn, log.FC = log(1.2),ident.1 = c('Neuron_Th'),one.sided = T)
tf.mouse <- read.csv('/home/tkamath/DA/ligerpaper/tf.name',sep = '\t',header = F)
View(markers.th[markers.th$feature %in% tf.mouse$V2,])

# nurr integrations
lig.neg <- qread('/home/tkamath/DA/nurrintegration/ligall_nurrneg.qs')
lig.neg <- louvainCluster(lig.neg)
plotByDatasetAndCluster(lig.neg)
lig.neg.annot <- lig.neg
levels(lig.neg.annot@clusters) <- c('Olig','Astro','Olig','Olig','MG',
                                    'Olig','Olig','Olig','MG','OPC',
                                    'Olig','Astro','Olig','Olig','Olig',
                                    'Olig','Olig','Olig','Endo/pericyte','Non-DA neuron',
                                    'DA neuron','Non-DA neuron','Non-DA neuron','MG',
                                    'Endo/pericyte','Endo/pericyte','Olig','Olig')
lig.neg.annot@clusters <- factor(lig.neg.annot@clusters,levels = c('DA neuron','Non-DA neuron','Astro','MG',
                                                                   'OPC','Olig','Endo/pericyte'))
p2 <- plotByDatasetAndCluster(lig.neg.annot,return.plots = T,do.legend = T,text.size = 0)
p2[[2]] <- p2[[2]] + theme(axis.text = element_blank(),
                           axis.title = element_blank(),
                           axis.ticks = element_blank(),
                           axis.line = element_blank(),legend.position = 'none') + scale_color_brewer(palette = 'Dark2')
pdf('tsne_ligneg.pdf',useDingbats = F)
p2[[2]]
dev.off()

pdf_convert('tsne_ligneg.pdf',dpi = 400,format = 'png')

raw.neg <- MergeSparseDataAll(lig.neg.annot@raw.data)
lig.neg.seurat <- CreateSeuratObject(raw.neg)
lig.neg.seurat@ident <- lig.neg.annot@clusters

lig.neg.seurat@ident <- factor(lig.neg.seurat@ident,levels = c('DA neuron','Non-DA neuron','Astro',
                                                               'MG','OPC','Olig','Endo/pericyte'))
p1 <- DotPlot(lig.neg.seurat,genes.plot = c('TH','SLC17A6','GAD1','GJA1','CX3CR1',
                                            'CSPG4','OLIG1','CLDN5'),
              x.lab.rot = T,do.return = T,cols.use = c('orange','red'),plot.legend = T)
pdf('suppfig1/dotplot_ligneg.pdf',useDingbats = F)
p1
dev.off()


p1 <- StackedVlnPlot(lig.neg.seurat,features = rev(c('CLDN5','OLIG1','VCAN','CX3CR1',
                                                     'AQP4','RBFOX3','TH',
                                                     'SLC6A3','SLC18A2')) )
pdf('stackvln_ligneg.pdf',useDingbats = F,height = 14, width = 14)
p1 + theme(axis.text.x = element_text(angle = 90))
dev.off()


lig.pos.annot <- qread('/home/tkamath/DA/nurrintegration/ligall_nurrpos.qs')
levels(lig.pos.annot@clusters)[c(1,2)] <- c('Non-DA neuron','Non-DA neuron')
lig.pos.annot@clusters <- factor(lig.pos.annot@clusters,levels = c('DA neuron','Non-DA neuron','Astrocyte',
                                                                   'Microglia/macrophage','OPC','Oligodendrocyte',
                                                                   'Endothelial cell/pericyte'))

p1 <- plotByDatasetAndCluster(lig.pos.annot,return.plots = T,do.legend = T,text.size = 0)
setwd('suppfig1/')
pdf('tsne_ligpos.pdf')
p1[[2]]+ theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_blank(),legend.position = 'none') + scale_color_brewer(palette = 'Dark2')
dev.off()
library(pdftools)
pdf_convert('tsne_ligpos.pdf',format = 'png',dpi = 400)

raw.pos <- MergeSparseDataAll(lig.pos.annot@raw.data)
lig.pos.seurat <- CreateSeuratObject(raw.pos)
lig.pos.seurat@ident <- lig.pos.annot@clusters

p1 <- StackedVlnPlot(lig.pos.seurat,features = rev(c('CLDN5','OLIG1','VCAN','CX3CR1',
                                                     'AQP4','SYT1','TH',
                                                     'SLC6A3','SLC18A2')))
pdf('stackvln_ligpos.pdf',useDingbats = F,height = 14, width = 14)
p1 + theme(axis.text.x = element_text(angle = 90))
dev.off()

p1 <- rliger::plotGene(lig.pos, gene = 'TH',plot.by = 'none',return.plots = T,do.legend = F)
png(paste0(base.dir,'suppfig1/thexpression_nurrpos.png'),res = 300,width = 1000, height = 1000)
p1 + theme(axis.text = element_blank(),
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           axis.line = element_blank())
dev.off()

p1 <- rliger::plotGene(lig.neg, gene = 'TH',plot.by = 'none',return.plots = T,do.legend = F)
png(paste0(base.dir,'suppfig1/thexpression_nurrneg.png'),res = 300,width = 1000, height = 1000)
p1 + theme(axis.text = element_blank(),
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           axis.line = element_blank())
dev.off()

cleaned.meta <- qread('/home/tkamath/DA/cleanedmetadata.qs')
cleaned.meta <- bind_rows(cleaned.meta,.id = 'names')
cleaned.meta <- cleaned.meta[cleaned.meta$disease == 'Ctrl',]
#cleaned.meta[which(grepl('Ex',cleaned.meta$clusters)),]$names <- 'Ex'
#cleaned.meta[which(grepl('Inh',cleaned.meta$clusters)),]$names <- 'Inh'
cleaned.meta$names <- as.factor(cleaned.meta$names)
cleaned.meta$names <- factor(cleaned.meta$names,levels = 
                               c('da','nonda','opc','astro', 'olig','endo','mg'))
levels(cleaned.meta$names) <- c('DA neuron','Non-DA neuron','OPC','Astrocyte',
                                'Oligodendrocyte','Endothelial cell','MG')
cleaned.meta <- cleaned.meta %>% mutate(nurr = ifelse(grepl('DAPI',lib),'dapi','nurr' ))

t1 <- aggregate(nGene ~ subject,cleaned.meta,median)
t1$column <- rep(1,nrow(t1))
cleaned.meta$subject <- droplevels(cleaned.meta$subject)
cleaned.meta$subject <- factor(cleaned.meta$subject,
                               levels = c('3346','5610','3482','6173','4956','3345','3322','3298'))

pdf(paste0(base.dir,'suppfig1/qcmetrics2.pdf'),useDingbats = F)
ggplot(data  = cleaned.meta,aes(x = subject,y = nGene,fill = subject)) + geom_violin() +
  geom_boxplot(width = 0.2, fill = 'white') + scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20),legend.position = 'none',
        axis.text.y = element_text(size = 20)) + ylab('Number of genes')

ggplot(data  = cleaned.meta,aes(x = subject,y = nUMI,fill = subject)) + geom_violin() +
  geom_boxplot(width = 0.2,fill = 'white') + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20),legend.position = 'none',
        axis.text.y = element_text(size = 20)) + ylab('Number of UMIs')
dev.off()

mycolors = RColorBrewer::brewer.pal(8,'Dark2')
pdf(paste0(base.dir,'suppfig1/qcmetrics.pdf'),useDingbats = F,height = 8, width = 8*44/26)
ggplot(data = cleaned.meta,aes(x = names, y = log10(nGene),fill = names)) + geom_violin() +
  geom_boxplot(width = 0.2,fill = 'white') +
  theme(axis.text.x = element_blank(),legend.position = 'none',
        axis.text.y = element_text(size = 20)) + ylab('log10(number of genes)') +
  scale_fill_manual(values = c(mycolors[1],mycolors[2],mycolors[5],
                               mycolors[3],mycolors[6],mycolors[7],mycolors[4]))
dev.off()

pdf(paste0(base.dir,'suppfig1/qcmetrics_3.pdf'),useDingbats = F,height = 8, width = 8*44/37)
ggplot(data = cleaned.meta,aes(x = names, y = log10(nUMI),fill = names)) + geom_violin() +
  geom_boxplot(width = 0.2,fill = 'white') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20),legend.position = 'none',
        axis.text.y = element_text(size = 20)) + ylab('log10(number of UMIs)')+
  scale_fill_manual(values = c(mycolors[1],mycolors[2],mycolors[5],
                               mycolors[3],mycolors[6],mycolors[7],mycolors[4]))
dev.off()

