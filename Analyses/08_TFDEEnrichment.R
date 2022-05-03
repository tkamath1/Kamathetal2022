library(Seurat)
library(rliger)
library(qs)
library(MAST)
library(dplyr)
library(lme4)
library(NNLM)
library(devtools)
library(enrichR)
library(Rtsne)
library(UpSetR)
source("/home/tkamath/scripts/genevalidatehelper.R")
source("/home/tkamath/scripts/SeuratExtrafunctions.R")
source("/home/tkamath/scripts/extrafuncs.R")
source('/home/tkamath/scripts/prestowrapper.R')
setwd('/home/tkamath/DA/diffexp/')

setwd('/home/tkamath/DA/diffexp/da/')
source('/home/tkamath/scripts/genevalidatehelper.R')
mast.da <- lapply(grep('nonurr',list.files(),value = T),function(x){qread(x)})
mast.da <- CombineMAST(mast.da,component.pvalue.use = 'H',component.coef.use = 'D',
                       correction = 'fdr',contrast.use = 'statusDisease')
names(mast.da) <- as.factor(unlist(sapply(strsplit(substring(grep('nonurr',list.files(),value = T),first = 8),
                                                   '\\.'),`[`,1)))
mast.da.DE <- lapply(mast.da,function(x){x[x$fdr < 0.01 & x$coef > 0,]})

genes.in <- mast.da.DE$SOX6_AGTR1_nonurr$primerid
genes.out <- unique(as.vector(unlist(lapply(mast.da.DE[which(names(mast.da.DE) != 'SOX6_AGTR1_nonurr')],function(x){
  x$primerid
}))))
genes.use <- setdiff(genes.in,genes.out)
tf.table <- read.table('/home/tkamath/TF_names_v_1.01.txt')
tf.unique <- intersect(tf.table$V1,genes.use)

# Do it for all TFs...
trrust.db <- read.csv('/home/tkamath/DA/diffexp/TRRUST_Transcription_Factors_2019.txt',sep = '\t',header = F)
encode.db <- read.csv('/home/tkamath/DA/diffexp/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt',sep = '\t',header = F)
archs.db <- read.csv('/home/tkamath/DA/diffexp/ARCHS4_TFs_Coexp.txt',sep = '\t',header = F)
list.databases <- list('TRRUST' = trrust.db,'ENCODE' = encode.db,'ARCHS4' = archs.db)

mast.da.DEtot <- bind_rows(mast.da.DE,.id = 'names2')
tf.table <- read.csv('/home/tkamath/TF_names_v_1.01.txt',header = F)
tfs.test <- mast.da$SOX6_AGTR1_nonurr[mast.da$SOX6_AGTR1_nonurr$primerid %in% 
                                        tf.table$V1,]$primerid

genes.db <- lapply(tfs.test,function(x){
  c2 <- lapply(list.databases,function(z){
    c1 <- as.vector(unlist(z[grep(paste0('\\b',x,'\\b'),z$V1),]))
    c1 <- unique(c1)
    return(c1)
  })
  as.vector(unlist(c2))
})
names(genes.db) <- tfs.test
genes.db.use <- genes.db[which(sapply(genes.db,length) > 0)]

genes.db2 <- lapply(genes.db.use,function(z){
  a1 <- as.vector(unlist(mapIds(org.Hs.eg.db, keys = as.character(z),keytype = 'SYMBOL',column = 'ENSEMBL')))
  a2 <- a1[which(!is.na(a1))]
  a2 <- as.vector(unlist(mapIds(org.Hs.eg.db, keys = as.character(a2),keytype = 'ENSEMBL',column = 'SYMBOL')))
})
qsave(genes.db2,'~/DA/tfanalysis/genesdb2.qs')

out.plots <- lapply(unique(names(mast.da)),function(x){
  gene.set <- mast.da[[x]][which(mast.da[[x]]$fdr <= 0.99),]
  gene.set.use <- (sign(gene.set[order(gene.set$fdr,decreasing = F),]$coef))*(-log10(gene.set[
    order(gene.set$fdr,decreasing = F),]$`Pr(>Chisq)`))
  names(gene.set.use) <- gene.set[order(gene.set$fdr,decreasing = F),]$primerid
  fgseaRes <- fgsea(pathways = genes.db2,nperm = 10000,
                    stats = gene.set.use,minSize = 1,maxSize = 500)
  return(fgseaRes)
})

names(out.plots) <- unique(names(mast.da))

out.plots.use <- lapply(out.plots,function(x){x[x$padj < 0.05,]$pathway})

setdiff(out.plots.use$SOX6_AGTR1_nonurr,as.vector(unlist(out.plots.use[c(1:6,8:10)])))

out.results.use <- lapply(names(out.results),function(z){
  names(out.results[[z]]) <- unique(names(mast.da))
  out.return <- bind_rows(out.results[[z]],.id = 'names')
  return(out.return)
})

qsave(out.results.use,'~/DA/diffexp/outresults_TFDE.qs')
