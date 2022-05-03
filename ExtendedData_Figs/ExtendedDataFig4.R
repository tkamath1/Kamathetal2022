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
################# Supp Figure 4#######################
######################################################

# Load in scenic data
library(scico)
tf.to.plot <- qread('/home/tkamath/DA/tfanalysis/tf_forsupp.qs')
limits = c(-1,1)*max(tf.to.plot$value)
out.use <-unique(as.vector(unlist(lapply(unique(tf.to.plot$DA_type),function(x){
  tmp1 <- tf.to.plot[tf.to.plot$DA_type == x,]
  as.vector(tmp1[order(tmp1$value,decreasing = T)[1:10],]$Regulon)
}))))


tf.use <- tf.to.plot[tf.to.plot$Regulon %in% out.use,]
tf.use$Regulon <- droplevels(tf.use$Regulon)
out.vals <- reshape(tf.use, idvar = "Regulon", timevar = "DA_type", direction = "wide")
rownames(out.vals) <- out.vals$Regulon
out.vals$Regulon <- NULL
p1 <- pheatmap::pheatmap(out.vals)
pdf('~/DA/tmp2.pdf')
p1
dev.off()

as.vector(unlist(sapply(strsplit(colnames(out.vals),'\\.'),`[`,2)))

tf.use$Regulon <- reorder(tf.use$Regulon,new.order = rownames(out.vals))
tf.use$Regulon <- reorder(tf.use$Regulon,new.order = p1[[1]]$order)
#tf.use$DA_type <- reorder(tf.use$DA_type,new.order = p1$tree_col$order)

pdf('/home/tkamath/DA/figures/suppfigures/tfs_forsupp.pdf',useDingbats = F,width = 28)
ggplot(tf.use,aes(y = DA_type,x = Regulon, fill = value)) + geom_tile() +
  scale_fill_scico(palette = 'cork',limits= limits) + 
  theme(axis.text.x = element_text(angle = 90,size = 25,vjust = 1),
        axis.text.y = element_text(size = 20))
dev.off()

tf.human <- qread('/home/tkamath/DA/tfanalysis/tf_toplot.qs')
tf.macaque <- qread('/home/tkamath/DA/tfanalysis/tf_forsupp_macaque.qs')


scenicOptions <- qread("~/DA/tfanalysis/macaque_SCENIC/scenicoptions_macaque.qs") 
setwd('~/DA/tfanalysis/')
regulonAUC<-loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC<-regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
order.da <- rev(c('SOX6_AGTR1','SOX6_PART1',
                  'SOX6_DDT','SOX6_GFRA2',
                  'CALB1_CALCR','CALB1_CRYM_CCDC68','CALB1_PPP1R17',
                  'CALB1_RBP4','CALB1_GEM','CALB1_TRHR'))
ligerex<-qread("/home/tkamath/DA/species/newmacaque/newmacaque_human.qs")
macaque.idx <- rownames(ligerex@cell.data[which(
  ligerex@cell.data$dataset == 'macaque'),])
ligerex <- subsetLiger(ligerex,cells.use =macaque.idx,remove.missing = F )
ligerex <- ligerToSeurat(ligerex)
exprMat<-as.array(GetAssayData(ligerex, slot = 'data')) #Note this is normalized data 
cellInfo<-data.frame(ligerex@ident)#Pull the cluster information from the Seurat object (barcode and clustera ids) 
cellInfo$bc<-rownames(cellInfo)
num_clusters<-length(unique(cellInfo$ligerex.ident))

cellInfo$ligerex.ident <- factor(cellInfo$ligerex.ident,levels = rev(order.da))
regulonActivity_byCellType<-sapply(split(rownames(cellInfo), cellInfo$ligerex.ident),
                                   function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
minVal<-0.05
regulonActivity_byCellType2<-regulonActivity_byCellType #Create a new object to not disrupt downstream analysis
regulonActivity_byCellType2<-regulonActivity_byCellType2[which(rowSums(regulonActivity_byCellType2>minVal)>0),]

regulonActivity_byCellType_Scaled<-t(scale(t(regulonActivity_byCellType), center = T, scale=T))

macaque.regulons <- melt(regulonActivity_byCellType_Scaled)
macaque.regulons$Var1 <- sapply(lapply(as.character(macaque.regulons$Var1),
                                       function(x){strsplit(x," ")[[1]]}),`[`,1)
macaque.regulons$Var1 <- as.factor(unlist(sapply(strsplit(macaque.regulons$Var1,'_'),`[`,1)))

both.regulons <- merge(y = macaque.regulons,x = tf.human, by.y = c('Var1','Var2'), by.x = c('Regulon','DA_type'))

# Get Unique lines in the data table
both.regulons <- both.regulons %>% distinct(Regulon, DA_type, .keep_all = TRUE)

both.regulons$meanboth <- (both.regulons$value.x + both.regulons$value.y)/2
lapply(unique(both.regulons$DA_type),function(x){
  tmp1 <- both.regulons[which(both.regulons$DA_type == x),]
  tmp1[order(tmp1$meanboth,decreasing = T)[1:10],]$Regulon
})

pdf('~/DA/figures/suppfigures/suppfig4/macaque_human_TFscatter.pdf',useDingbats = F,
    width = 12,height = 4)
ggplot(both.regulons,aes(x = value.x, y = value.y,fill = DA_type,shape = 21)) + 
  geom_point(size = 2.5) +scale_shape_identity() + facet_wrap(~DA_type,nrow = 2) +
  theme(legend.position = 'none') + xlab('Macaque scaled regulon scores') + 
  ylab('Human scaled regulon scores') +theme(
    strip.text.x = element_text(
      size = 12, face = "bold"
    ),
    strip.text.y = element_text(
      size = 12, face = "bold"
    )
  )
dev.off()
