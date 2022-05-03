library(Seurat)
library(qs)
library(MAST)
library(dplyr)
library(lme4)
library(org.Hs.eg.db)
source("/home/tkamath/scripts/genevalidatehelper.R")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/gwas/')
options(mc.cores = 20)

pd.meta.use <- qread('/home/tkamath/DA/pdmeta_use.qs')

setwd('/home/tkamath/DA/diffexp/da/')
source('/home/tkamath/scripts/genevalidatehelper.R')

p.thresh <- 0.05
ngene.thresh <- 3500
setwd('/home/tkamath/DA/gwas/genesets/')
source('/home/tkamath/scripts/genevalidatehelper.R')
out.mast <- lapply(grep('mastDE',list.files(),value = T),function(x){qread(x)})
out.mast <- CombineMAST(out.mast,component.pvalue.use = 'H',component.coef.use = 'D',
                        correction = 'fdr',contrast.use = 'clusteruse1')
out.mast <- lapply(out.mast,function(x){x[x$fdr < p.thresh & x$coef > 0,]})
names(out.mast) <- sapply(strsplit(sapply(strsplit(grep('mastDE',list.files(),value = T),'_'),`[`,2),'\\.'),`[`,1)

lapply(names(out.mast),function(z){
  out.mast.use <- out.mast[[z]]
  out.mast.use$gene.use <- mapIds(org.Hs.eg.db, keys = as.character(out.mast.use$primerid),
                                  keytype = 'SYMBOL',column = 'ENSEMBL')
  out.mast.use <- out.mast.use[which(!is.na(out.mast.use$gene.use)),]
  df.print <- data.frame('gene' = out.mast.use[order(out.mast.use$z,decreasing = T),]$gene.use[1:min(nrow(out.mast.use),ngene.thresh)],
                         'set' = rep(as.character(z), min(nrow(out.mast.use),ngene.thresh) ))
  df.print$gene <- mapIds(org.Hs.eg.db, keys = as.character(df.print$gene),keytype = 'SYMBOL',column = 'ENSEMBL')
  write.table(df.print,paste0('/home/tkamath/DA/gwas/genesets/',ngene.thresh,'_',p.thresh,'_',as.character(z),
                              '_mast_geneset.csv'),quote = F,row.names = F,col.names = F,sep = '\t')
})  

out.pvals <- mclapply(names(out.mast),function(z){
  setwd('/home/tkamath/magma/')
  system(paste0('./magma --gene-results PDFUMAresults/magma.genes.raw --set-annot /home/tkamath/DA/gwas/genesets/',ngene.thresh,'_',
                p.thresh,'_',as.character(z),'_mast_geneset.csv gene-col=1 set-col=2 --out pdsn',ngene.thresh,'_',
                p.thresh,'_',as.character(z),'mast'))
},mc.cores = 15)

celltypes <- names(out.mast)
qsave(celltypes,'/home/tkamath/DA/gwas/celltypes_test.qs')
df.use <- data.frame('celltype' = celltypes, 'magmap' = as.vector(unlist(lapply(celltypes,
                                                                                function(x){magma.use <- read.table(paste0('/home/tkamath/magma/pdsn',ngene.thresh,'_',
                                                                                                                           p.thresh,'_',as.character(x),'mast.gsa.out'),header = T)
                                                                                return(-log10(magma.use$P))}))))

df.use$celltype <- factor(df.use$celltype,levels = rev(c('da','olig','Ex','Inh','opc','endo','mg','astro')))
levels(df.use$celltype) <- rev(c('DA neuron','Oligodendrocyte','Excitatory neuron','Inhibitory neuron','OPC',
                                 'Endothelial cell','Microglia/macrophage','Astrocyte'))
qsave(df.use,'/home/tkamath/DA/gwas/pdmagma_mastmarkers.qs')

# MAST AD
p.thresh <- 0.05
ngene.thresh <- 3500
setwd('/home/tkamath/DA/gwas/genesets/')
source('/home/tkamath/scripts/genevalidatehelper.R')
out.mast <- lapply(grep('mastDE',list.files(),value = T),function(x){qread(x)})
out.mast <- CombineMAST(out.mast,component.pvalue.use = 'H',component.coef.use = 'D',
                        correction = 'fdr',contrast.use = 'clusteruse1')
out.mast <- lapply(out.mast,function(x){x[x$fdr < p.thresh & x$coef > 0,]})
names(out.mast) <- sapply(strsplit(sapply(strsplit(grep('mastDE',list.files(),value = T),'_'),`[`,2),'\\.'),`[`,1)

lapply(names(out.mast),function(z){
  out.mast.use <- out.mast[[z]]
  df.print <- data.frame('gene' = out.mast.use[order(out.mast.use$z,decreasing = T),]$primerid[1:ngene.thresh],
                         'set' = rep(as.character(z),ngene.thresh))
  df.print$gene <- mapIds(org.Hs.eg.db, keys = as.character(df.print$gene),keytype = 'SYMBOL',column = 'ENTREZID')
  write.table(df.print,paste0('/home/tkamath/DA/gwas/genesets/ad_',p.thresh,'_',ngene.thresh,'_',as.character(z),'_mast_geneset.csv'),
              quote = F,row.names = F,col.names = F,sep = '\t')
})  

out.pvals <- mclapply(names(out.mast),function(z){
  setwd('/home/tkamath/magma/')
  system(paste0('./magma --gene-results AD2021/AD2021.genes.raw --set-annot /home/tkamath/DA/gwas/genesets/ad_',
                p.thresh,'_',ngene.thresh,'_',as.character(z),
                '_mast_geneset.csv gene-col=1 set-col=2 --out adsn',ngene.thresh,'_',
                p.thresh,'_',as.character(z),'mast'))
},mc.cores = 15)

# Combine SN into one data frame
out.df <- lapply(c('ad','pd'),function(z){
  out.pvals2 <- lapply(celltypes,function(x){
    magma.use <-read.table(paste0('/home/tkamath/magma/',as.character(z),'sn3500_0.01_',as.character(x),
                                  'mast.gsa.out'),header = T)
    magma.use$P})
  names(out.pvals2) <- celltypes
  out.pvals2 <- do.call(rbind,out.pvals2)
  colnames(out.pvals2) <- z
  return(out.pvals2)
})
names(out.df) <- c('ad','pd')
out.df <- do.call(cbind,out.df)
out.df <- melt(out.df)
qsave(out.df,'/home/tkamath/DA/gwas/snmagma_majorcelltypes_01.qs')

##################
### Striatum #####
##################
setwd('/home/tkamath/DA/striatum/')
source("/home/tkamath/scripts/genevalidatehelper.R")
out.mast <- lapply(grep('mastDE',list.files(),value = T),function(x){qread(x)})
out.mast <- CombineMAST(out.mast,component.pvalue.use = 'H',component.coef.use = 'D',
                        correction = 'fdr',contrast.use = 'clusteruse1')
out.mast <- lapply(out.mast,function(x){x[x$fdr < 0.05 & x$coef > 0,]})
names(out.mast) <- sapply(strsplit(sapply(strsplit(grep('mastDE',list.files(),value = T),'_'),`[`,2),'\\.'),`[`,1)

lapply(names(out.mast),function(z){
  out.mast.use <- out.mast[[z]]
  df.print <- data.frame('gene' = out.mast.use[order(out.mast.use$z,decreasing = T),]$primerid[1:3500],
                         'set' = rep(as.character(z),3500 ))
  df.print$gene <- mapIds(org.Hs.eg.db, keys = as.character(df.print$gene),keytype = 'SYMBOL',column = 'ENSEMBL')
  write.table(df.print,paste0('/home/tkamath/DA/striatum/pdcn_3500_0.05_',as.character(z),'_mast_geneset.csv'),quote = F,
              row.names = F,col.names = F,sep = '\t')
})  

lapply(names(out.mast),function(z){
  out.mast.use <- out.mast[[z]]
  df.print <- data.frame('gene' = out.mast.use[order(out.mast.use$z,decreasing = T),]$primerid[1:3500],
                         'set' = rep(as.character(z),3500 ))
  df.print$gene <- mapIds(org.Hs.eg.db, keys = as.character(df.print$gene),keytype = 'SYMBOL',column = 'ENTREZID')
  write.table(df.print,paste0('/home/tkamath/DA/striatum/pdcn_3500_0.05_',as.character(z),'_mast_geneset2.csv'),quote = F,
              row.names = F,col.names = F,sep = '\t')
})  

out.pvals <- mclapply(names(out.mast),function(z){
  setwd('/home/tkamath/magma/')
  system(paste0('./magma --gene-results PDFUMAresults/magma.genes.raw --set-annot /home/tkamath/DA/striatum/pdcn_3500_0.05_',
                as.character(z),'_mast_geneset.csv gene-col=1 set-col=2 --out pdstriatum',as.character(z),'mast'))
},mc.cores = 15)

out.pvals <- mclapply(names(out.mast),function(z){
  setwd('/home/tkamath/magma/')
  system(paste0('./magma --gene-results AD2021/AD2021.genes.raw --set-annot /home/tkamath/DA/striatum/pdcn_3000_0.05_',
                as.character(z),'_mast_geneset2.csv gene-col=1 set-col=2 --out adstriatum',as.character(z),'mast'))
},mc.cores = 15)

celltypes <- names(out.mast)
out.df <- lapply(c('ad','pd'),function(z){
  out.pvals2 <- lapply(celltypes,function(x){
    magma.use <-read.table(paste0('/home/tkamath/magma/',as.character(z),'striatum',as.character(x),'mast.gsa.out'),
                           header = T)
    magma.use$P
  })
  names(out.pvals2) <- celltypes
  out.pvals2 <- do.call(rbind,out.pvals2)
  colnames(out.pvals2) <- z
  return(out.pvals2)
})
names(out.df) <- c('ad','pd')
out.df <- do.call(cbind,out.df)
out.df <- melt(out.df)
qsave(out.df,'/home/tkamath/DA/striatum/gwas_analysis.qs')

### Subtypes
setwd('/home/tkamath/DA/gwas/genesets/subtypes/')
out.mast <- lapply(grep('mastDE',list.files(),value = T),function(x){qread(x)})
out.mast <- CombineMAST(out.mast,component.pvalue.use = 'H',component.coef.use = 'D',
                        correction = 'fdr',contrast.use = 'clusteruse1')
out.mast <- lapply(out.mast,function(x){
  x$geneuse <- mapIds(org.Hs.eg.db, keys = as.character(x$primerid),keytype = 'SYMBOL',column = 'ENSEMBL')
  x <- x[which(!is.na(x$geneuse)),]
  return(x)
})
out.mast <- lapply(out.mast,function(x){x[ x$z > 0 & x$fdr < 0.05,]})
names(out.mast) <- sapply(strsplit(substring(grep('mastDE',list.files(),value = T),first = 8),'\\.'),`[`,1)

qsave(names(out.mast),'/home/tkamath/DA/gwas/genesets/subtypes/subtypes_test.qs')

subtypes.test <- names(out.mast)

lapply(subtypes.test,function(z){
  out.mast.use <- out.mast[[z]]
  out.mast.use$gene.use <- mapIds(org.Hs.eg.db, keys = as.character(out.mast.use$primerid),
                                  keytype = 'SYMBOL',column = 'ENSEMBL')
  out.mast.use <- out.mast.use[which(!is.na(out.mast.use$gene.use)),]
  
  df.print <- data.frame('gene' = out.mast.use[order(out.mast.use$z,decreasing = T),]$gene.use[1:min(nrow(out.mast.use),3500) ],
                         'set' = rep(as.character(z),min(nrow(out.mast.use),3500) ))
  #df.print$gene <- mapIds(org.Hs.eg.db, keys = as.character(df.print$gene),keytype = 'SYMBOL',column = 'ENSEMBL')
  write.table(df.print,paste0('/home/tkamath/DA/gwas/genesets/subtypes/magmagenerated/',as.character(z),'_mast_geneset_subtypes.csv'),
              quote = F,row.names = F,col.names = F,sep = '\t')
})

out.pvals <- mclapply(subtypes.test,function(z){
  setwd('/home/tkamath/magma/')
  system(paste0('./magma --gene-results PDFUMAresults/magma.genes.raw --set-annot /home/tkamath/DA/gwas/genesets/subtypes/magmagenerated/',
                as.character(z),'_mast_geneset_subtypes.csv gene-col=1 set-col=2 --out a1pd',as.character(z),'mast'))
},mc.cores = 40)

celltypes <- names(out.mast)
df.use <- data.frame('celltype' = celltypes, 'magmap' = as.vector(unlist(lapply(celltypes,
                                                                                function(x){magma.use <- read.table(paste0('/home/tkamath/magma/a1pd',as.character(x),'mast.gsa.out'),header = T)
                                                                                return(-log10(magma.use$P))}))))

df.use$celltype <- reorder(df.use$celltype, df.use$magmap)
df.use$cluster <- as.character(unlist(sapply(strsplit(as.character(df.use$celltype),'_'),`[`,1)))
df.use[which(grepl('SOX6|CALB1',df.use$cluster)),]$cluster <- 'DA'
pdf('/home/tkamath/DA/gwas/PDmagmaenrichment_mastsubtypesmarkers.pdf',width = 10,height = 10,useDingbats = F)
ggplot() + geom_bar(data = df.use,aes( x = celltype, y = magmap,fill = cluster),stat = 'identity') + 
  coord_flip() +
  ylab('-log10(p-value, MAGMA)') + xlab('Cell type') + 
  geom_hline(yintercept=-log10(0.05/68), linetype="dashed", color = "red")+
  gghighlight(cluster == 'DA') + scale_fill_brewer(palette = "Dark2")
dev.off()

########################
## Robustness testing###
########################
gene.size <- c(1000,1500,2000,2500,3000,3500,4000)
pthresh <- c(0.01,0.05)
cell.subtypes <- qread('/home/tkamath/DA/gwas/genesets/subtypes/subtypes_test.qs')
subtypes.test <- qread('/home/tkamath/DA/gwas/genesets/subtypes/subtypes_test_above1000_0.05.qs')
da.subtypes <- cell.subtypes[c(10:15,65:68)]

out.vec3 <- lapply(pthresh,function(y){
  out.vec2 <- lapply(gene.size,function(x){
    out.vec1 <- lapply(subtypes.test,function(z){
      if (file.exists(paste0('/home/tkamath/magma/pd',as.character(y),'_',as.character(x),'_',as.character(z),'mast.gsa.out'))){
        magma.use <- read.table(paste0('/home/tkamath/magma/pd',as.character(y),'_',as.character(x),'_',as.character(z),'mast.gsa.out'),header = T)
        return(magma.use$P)
      }
      else
        return(1)
    })
    names(out.vec1) <- subtypes.test
    return(bind_rows(out.vec1,.id = 'names') )
  })
  names(out.vec2) <- gene.size
  return(bind_rows(out.vec2,.id = 'Gene_size'))
})
names(out.vec3) <- pthresh
df.plot <- melt(as.data.frame(bind_rows(out.vec3,.id = 'names2')))
colnames(df.plot) <- c('FDR_thresh', 'Gene_set_size','Cell_subtype','p-val')
library(RColorBrewer)
df.plot$Cell_type <- as.vector(unlist(sapply(strsplit(as.character(df.plot$Cell_subtype),'_'),`[`,1)))
df.plot[grep('CALB1|SOX6',df.plot$Cell_type),]$Cell_type <- 'DA'
df.plot$Cell_type <- as.factor(df.plot$Cell_type)
qsave(df.plot,'/home/tkamath/DA/gwas/df_plot_magma_robustness.qs')


rare.variants <- read.csv('/home/tkamath/DA/gwas/supp15_IPDrarevariants.csv',header = T)
rare.variants <- rare.variants$Known.PD.gene.in.EOPD
rare.variants <- as.factor(rare.variants)
levels(rare.variants)[17] <- 'PARK2'

markers.wilcox <- qread('/home/tkamath/DA/gwas/markergenes_wilcox.qs')
markers.wilcox$pctdiff <- markers.wilcox$pct_in - markers.wilcox$pct_out
# 
setwd('/home/tkamath/DA/gwas/genesets/')
out.mast <- lapply(grep('mastDE',list.files(),value = T),function(x){qread(x)})
out.mast <- CombineMAST(out.mast,component.pvalue.use = 'H',component.coef.use = 'D',
                        correction = 'fdr',contrast.use = 'clusteruse1')
names(out.mast) <- sapply(strsplit(substring(grep('mastDE',list.files(),value = T),first = 8),'\\.'),`[`,1)
out.mast.tot <- bind_rows(out.mast,.id = 'names')
out.mast.tot.use <- out.mast.tot

out.mast.tot.use <- merge(out.mast.tot.use,markers.wilcox,by.x = c('names','primerid'),by.y = c('group','feature'))

out.test <- lapply(unique(out.mast.tot.use$names),function(x){
  out.mast.test <- out.mast.tot.use[out.mast.tot.use$names == x,]
  y <- floor(nrow(out.mast.test)/10)
  print(y)
  out.mast.test <- out.mast.test[order(out.mast.test$auc,decreasing = T)[c(1:y )],]
  geneA <- out.mast.test$primerid
  geneB <- rare.variants
  go.obj <- newGeneOverlap(geneB,geneA)
  go.obj <- testGeneOverlap(go.obj)
  return(list(go.obj@pval,go.obj@odds.ratio,go.obj@intersection))
})
df.use <- data.frame('celltype' = unique(out.mast.tot.use$names),'hpval' = as.vector(unlist(sapply(out.test,`[`,1))),
                     'odds' = as.vector(unlist(sapply(out.test,`[`,2))),
                     'intersection' = length(sapply(out.test,`[`,3)))
df.use$celltype <- reorder(df.use$celltype,df.use$hpval)
qsave(df.use,'/home/tkamath/DA/gwas/rarevariant_enrichmentpvalues.qs')

matrix.out <- matrix(nrow= length(rare.variants),ncol =length(unique(out.mast.tot.use$names)))
for (i in 1:length(rare.variants)){
  for (j in 1:length(unique(out.mast.tot$names))){
    coef.use <- out.mast.tot.use[out.mast.tot.use$primerid == rare.variants[i] &
                                   out.mast.tot.use$names == unique(out.mast.tot.use$names)[j],]$`Pr(>Chisq)`
    z.use <- out.mast.tot.use[out.mast.tot.use$primerid == rare.variants[i] &
                                out.mast.tot.use$names == unique(out.mast.tot.use$names)[j],]$auc
    if( length(coef.use != 0)){
      matrix.out[i,j] <- z.use
    }
    else{
      matrix.out[i,j] <- 0
    }
    
  }
}
rownames(matrix.out) <- rare.variants
colnames(matrix.out) <-unique(out.mast.tot.use$names)

matrix.out[which(!is.finite(matrix.out))] <- max(matrix.out[which(is.finite(matrix.out))])

out.mast.rare <- as.data.frame(melt(matrix.out))
p1 <- ComplexHeatmap::Heatmap(t(matrix.out))
order.cols <- ComplexHeatmap::column_order(p1)
order.rows <- ComplexHeatmap::row_order(p1)
out.mast.rare$Var1 <- reorder(out.mast.rare$Var1,new.order = order.cols)
out.mast.rare$Var2 <- reorder(out.mast.rare$Var2,new.order = rev(as.vector(unlist(order.rows))))

levels(out.mast.rare$Var2) <- c('DA neuron','Inhibitory neuron','Excitatory neuron',
                                'Endothelial cell/pericyte','OPC', 'Oligodendrocyte','Microglia/macrophage','Astrocyte')

out.mast.rare$Var1 <- reorder(out.mast.rare$Var1,new.order = c(13:26,1:12))
qsave(out.mast.rare,'/home/tkamath/DA/gwas/rarevariants_mast.qs')


#  Interpret GWAS (for response to reviewers)
library(Seurat)
library(qs)
library(MAST)
library(dplyr)
library(lme4)
library(org.Hs.eg.db)
source("/home/tkamath/scripts/genevalidatehelper.R")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/gwas/')
options(mc.cores = 60)

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

library(UpSetR)
pdf('/home/tkamath/DA/figures/revisions_plot/upsetplot_DAmarkers.pdf',useDingbats = F)
upset(fromList(genes.mast.use),nsets = 10,nintersects = 10)
dev.off()


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

magma.out$mastz <- da.mast.use$SOX6_AGTR1[match(magma.out$SYMBOL,da.mast.use$SOX6_AGTR1$primerid),]$z
magma.out$logp <- -log10(da.mast.use$SOX6_AGTR1[match(magma.out$SYMBOL,da.mast.use$SOX6_AGTR1$primerid),]$`Pr(>Chisq)`)
magma.out.agtr1 <- magma.out[which(!is.na(magma.out$mastz)),]
magma.out.agtr1 <- magma.out[which(!is.na(magma.out$logp)),]
setorderv(x = magma.out.agtr1,cols = 'mastz',order = -1)
setorderv(x = magma.out.agtr1,cols = 'logp',order = -1)
magma.out.agtr1$ZORDER <- rev(c(1:nrow(magma.out.agtr1)))

genes.nonagtr1 <- unique(as.vector(unlist(genes.mast.use[
  setdiff(c(1:10),grep('SOX6_AGTR1',names(genes.mast.use) )) ])))
genes.nonsox6 <- unique(as.vector(unlist(genes.mast.use[
  setdiff(c(1:10),grep('SOX6',names(genes.mast.use) )) ])))

magma.out.agtr1$spec <- NULL
magma.out.agtr1 <- magma.out.agtr1 %>% mutate(spec = ifelse(SYMBOL %notin% genes.nonagtr1,2,1))
magma.out.agtr1 <- magma.out.agtr1 %>% mutate(spec = ifelse(SYMBOL %in% genes.nonsox6,0,spec))
magma.out.agtr1 <- magma.out.agtr1 %>% mutate(spec = ifelse(ZSTAT > 4.568,spec,-1))
magma.out.agtr1$spec <- as.factor(magma.out.agtr1$spec)

magma.out.agtr1$ZORDER2 <- magma.out.agtr1$ZORDER/max(magma.out.agtr1$ZORDER)

magma.out.agtr1$shape <- rep(21,nrow(magma.out.agtr1))
#magma.out.agtr1[magma.out.agtr1$SYMBOL %in% c('LRRK2','GBA'),]$shape <- '23'
magma.out.agtr1$shape <- as.factor(magma.out.agtr1$shape)
qsave(magma.out.agtr1,'~/DA/gwas/gwas_pseudomanhattanmainfig.qs')

pdf('/home/tkamath/DA/gwas/dotplot_zscores_agtr1_rd2.pdf',useDingbats = F,width = 12)
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


