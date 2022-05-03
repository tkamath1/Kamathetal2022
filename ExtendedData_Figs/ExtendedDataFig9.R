# Striatum analysis
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

#Human striatum

source("source_code.R")
rm(list=ls())
.ArgList=.myArgConstructor("Allhuman")
.ArgList$HVG_list=3

hfiles=dir("~/data/striatum/human/rawData_arranged/")
dataList=list()
for(ifile in hfiles){
  tmp=qread(file.path("~/data/striatum/human/rawData_arranged",ifile))
  tmp$batch_merging=gsub("\\.qs","",ifile)
  tmp$ds_batch=tmp$batch_merging
  tmp$organism="Human"
  tmp=.myIndFileReaderFn(inputData=tmp,organism="Human",argList = .ArgList,convertToHuman=F,rankedBased=F)
  dataList=c(dataList,list(tmp))
  names(dataList)[length(data)]=gsub("\\.qs","",ifile)
  rm(tmp)
  gc()
}


data=list(data=dataList,data_m=NULL)
.data=data
data=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
.data2=data
dim(data$varFeatures)


tmp=data[c("varFeatures","allGenes" )]
save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))

cat('Scaling the data and PCA\n')
#argList = .ArgList;UMI_cor_thr=.ArgList$UMI_cor_thr
.myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)
rm(data)
gc()

res=parallel::mclapply(.ArgList$HVG_list,.myHarmonyFn,mc.cores = 7)
if(class(res)[1]=="matrix"){
  res=list(res)
}
names(res)=paste0("HVG_",.ArgList$HVG_list)
library(qs)

qsave(res,file=paste0("~/data/striatum/human_harmony/all_gropus/resHarmony_IndScaled_HVGlist","2_4","_geneFraction",.ArgList$allGenesFraction,"_MGIsymbol_correct_nocovariate.qs"))
qsave(.ArgList,file=paste0("~/data/striatum/human_harmony/all_gropus/","ArgList_HVGlist","_IndScaling_Feb16_MGIsymbol_correct_nocovariate.qs"))
#.ArgList=qread(paste0("~/data/striatum/human_harmony/all_gropus/","ArgList_HVGlist","_IndScaling_Feb16_MGIsymbol_correct_nocovariate.qs"))


#Figures
res=qread(paste0("~/data/striatum/human_harmony/all_gropus/resHarmony_IndScaled_HVGlist","2_4","_geneFraction",.ArgList$allGenesFraction,"_MGIsymbol_correct_nocovariate.qs"))


offsetCount=0
newRun=T
res_graph=parallel::mclapply(1:length(res),.AnnoFn_graph,res=res,newRun=newRun,.ArgList=.ArgList,offsetCount=offsetCount,changeSaveDir=F,mc.cores = length(res))
gc()

{
  .ArgList$HVG_count=3
  load(.myFilePathMakerFn("pca_anno",argList=.ArgList))
  
  harmony_res=res[[i]]
  #j=3;res=harmony_res;data_m=data$data;ncores=5;doStriatum=T;changeSaveDir=F
  tmp=parallel::mclapply(1:length(res),.AnnoFn,res=harmony_res,data_m=data$data,pd=pd,.ArgList=.ArgList,ncores=5,offsetCount=offsetCount,doStriatum=T,changeSaveDir=F,mc.cores = 2)
  gc()
  
}

#############
#arranging the high performing solution as a seurat object
#HVG3-v1

.ArgList=qread(paste0("~/data/striatum/human_harmony/all_gropus/","ArgList_HVGlist","_IndScaling_Feb16_MGIsymbol_correct_nocovariate.qs"))
res=qread(paste0("~/data/striatum/human_harmony/all_gropus/resHarmony_IndScaled_HVGlist","2_4","_geneFraction",.ArgList$allGenesFraction,"_MGIsymbol_correct_nocovariate.qs"))
j=1

hfiles=dir("~/data/striatum/human/rawData_arranged/")
dataList=list()
for(ifile in hfiles){
  tmp=qread(file.path("~/data/striatum/human/rawData_arranged",ifile))
  tmp$batch_merging=gsub("\\.qs","",ifile)
  tmp$ds_batch="human_postmortem"
  tmp$organism="Human"
  tmp=.myIndFileReaderFn(inputData=tmp,organism="Human",argList = .ArgList,convertToHuman=F,rankedBased=F)
  dataList=c(dataList,list(tmp))
  names(dataList)[length(data)]=gsub("\\.qs","",ifile)
  rm(tmp)
  gc()
}

data_m=.mycBindFn(inputList = dataList,batchNames = unlist(lapply(dataList,function(x) as.character(x$batch_merging[1])))) #44858 347459
tmp_name=unlist(lapply(strsplit(colnames(data_m),"_"),function(x) paste(x[(length(x)-2):length(x)],collapse = "_")))
colnames(data_m)=tmp_name

.ArgList$HVG_count=8
load(.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))

.ArgList$HVG_count=3
varFeatures=tmp$varFeatures[tmp$varFeatures[,2]>(.ArgList$HVG_count-1),1]
load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony1"),argList=.ArgList,pdf=F))

embeddings=res[[paste0("HVG_",.ArgList$HVG_count)]][[j]]

data_m=.extraExport2SeuratFn(inputData=data_m,project="RNA")
data_m=NormalizeData(data_m)
data_m@assays$RNA@var.features=varFeatures
data_m=ScaleData(data_m,features = varFeatures)

pd=pd[match(colnames(data_m),row.names(pd)),]

data_m[["umap"]]= CreateDimReducObject(embeddings = as.matrix(pd[,c("UMAP_1","UMAP_2")]), 
                                       key = "UMAP_", assay = "RNA", global = TRUE)


colnames(embeddings)=gsub("PC_","",colnames(embeddings))
data_m[["pca"]] <- CreateDimReducObject(embeddings = embeddings, 
                                        assay = "RNA",  
                                        key = "PC_")


data_m$anno_cluster_res=pd$anno_cluster_res
qsave(data_m,file="~/myBucket/human_striatum.qs")


##############
data_m=qread("~/myBucket/human_striatum.qs")
data_m$UMAP_1=data_m@reductions$umap@cell.embeddings[,1]
data_m$UMAP_2=data_m@reductions$umap@cell.embeddings[,2]

pd=as.data.frame(data_m@meta.data)
markerGenes=read.table("~/myBucket/markergenes/SPN_markers_set2",sep="\t",header = F,stringsAsFactors = F)
markerGenes=toupper(as.character(markerGenes[,1]))
markerGenes=markerGenes[!duplicated(markerGenes)]
p=.myFeaturePlot(inputSeurat = data_m,inputDimData = pd,inputGenes = markerGenes,combine_figs = T,order=F)
ggsave(plot=p,file="~/myBucket/plot.png",width=49,height=49)

markerGenes=c("GAD1","GAD2","SLC17A7","OLIG1","OLIG2","PPP1R1B","DRD1","DRD2","P2RY12","C1QB","AQP4","GABRB2","SST","PVALB","VIP","LAMP5","KCNJ8","KIT","DCN","CSPG4","FLT1","CD8","TOP2A","CNKSR3","BRINP3")
markerGenes=toupper(as.character(markerGenes))
markerGenes=markerGenes[!duplicated(markerGenes)]
p=.myFeaturePlot(inputSeurat = data_m,inputDimData = pd,inputGenes = markerGenes,combine_figs = T,order =F)
ggsave(plot=p,file="~/myBucket/plot2.png",width=19,height=19)

Idents(data_m)=as.character(data_m$anno_cluster_res)
p=.myDotPlot(object=data_m, features=markerGenes, cols = c("blue","red"))
ggsave(plot=p,file="~/myBucket/plot3.png",width=19,height=10)


data_m$anno_cluster_label="unassigned"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("34","29"))]="doublet_Olig-PPP1R1B"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("30","6","27"))]="Olig_CSPG4"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("0","12","7"))]="DRD1"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("22","2","9","16","39"))]="DRD2"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("26"))]="doublet_MG-Olig"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("21"))]="SST"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("10","13"))]="KIT"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("5"))]="MG"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("1","17","18","23","20"))]="Astro"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("14"))]="doublet_DRD-Olig"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("24"))]="doublet_Astro-GAD"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("25"))]="doublet_Astro-Olig"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("38","41"))]="doublet_DRD-Astro"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("32"))]="doublet_MG-Astro"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("40"))]="doublet_MG-Olig"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("28"))]="doublet_MG-DRD"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("8","35","31"))]="doublet_DRD-KIT"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("11","3","4"))]="Oligo"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("19","33"))]="DCN"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("15"))]="Flt1"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("29"))]="doublet_DRD-FLT-Olig"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("34"))]="doublet_DRD-Astro-Olig"
data_m$anno_cluster_label[which(data_m$anno_cluster_res %in% c("36","37"))]="doublet"

Idents(data_m)=data_m$anno_cluster_label

p=DimPlot(data_m, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(plot=p,file="~/myBucket/plot4.png",width=10,height=10)


data=list(countData=data_m@assays$RNA@counts,pdata=as.data.frame(data_m@meta.data),umapdata=data_m@reductions$umap@cell.embeddings,embedding=data_m@reductions$pca@cell.embeddings,fdata=as.data.frame(data_m@assays$RNA@meta.features))

qsave(data,file="~/myBucket/striatum_human/overall_classes.qs")

#Analysis of SPN clusters

rm(list=ls())
data_m=qread("~/myBucket/striatum_human/seurat4_obj.qs")
data_m$ds_batch=data_m$batch_merging
data=data_m[,which(Idents(data_m) %in% c("DRD1","DRD2"))]

data=SingleCellExperiment(assays = list(counts = data@assays$RNA@counts),colData = as.data.frame(data@meta.data),rowData=as.data.frame(data@assays$RNA@meta.features))
data=list(data=.mySplitObject(data,'batch_merging'),data_m=NULL)
.myAnalysisFn(cellClass="SPNhuman",inputData=data)


#Analysis of Oligo clusters

rm(list=ls())
data_m=qread("~/myBucket/striatum_human/seurat4_obj.qs")
data_m$ds_batch=data_m$batch_merging

table(Idents(data_m))

data=data_m[,which(Idents(data_m) %in% c("Oligo"))]

data=SingleCellExperiment(assays = list(counts = data@assays$RNA@counts),colData = as.data.frame(data@meta.data),rowData=as.data.frame(data@assays$RNA@meta.features))
data=list(data=.mySplitObject(data,'batch_merging'),data_m=NULL)
#cellClass="Oligohuman";inputData=data
.myAnalysisFn(cellClass="Oligohuman",inputData=data)


#Analysis of MG clusters
rm(list=ls())
data_m=qread("~/myBucket/striatum_human/seurat4_obj.qs")
data_m$ds_batch=data_m$batch_merging

table(Idents(data_m))

data=data_m[,which(Idents(data_m) %in% c("MG"))]

data=SingleCellExperiment(assays = list(counts = data@assays$RNA@counts),colData = as.data.frame(data@meta.data),rowData=as.data.frame(data@assays$RNA@meta.features))
data=list(data=.mySplitObject(data,'batch_merging'),data_m=NULL)
#cellClass="MGhuman";inputData=data
.myAnalysisFn(cellClass="MGhuman",inputData=data)



#Analysis of Astro clusters
rm(list=ls())
data_m=qread("~/myBucket/striatum_human/seurat4_obj.qs")
data_m$ds_batch=data_m$batch_merging

table(Idents(data_m))

data=data_m[,which(Idents(data_m) %in% c("Astro"))]

data=SingleCellExperiment(assays = list(counts = data@assays$RNA@counts),colData = as.data.frame(data@meta.data),rowData=as.data.frame(data@assays$RNA@meta.features))
data=list(data=.mySplitObject(data,'batch_merging'),data_m=NULL)
#cellClass="Astrohuman";inputData=data
.myAnalysisFn(cellClass="Astrohuman",inputData=data)


#Analysis of Interneuron clusters
rm(list=ls())
data_m=qread("~/myBucket/striatum_human/seurat4_obj.qs")
data_m$ds_batch=data_m$batch_merging

table(Idents(data_m))

data=data_m[,which(Idents(data_m) %in% c("SST","KIT"))]

data=SingleCellExperiment(assays = list(counts = data@assays$RNA@counts),colData = as.data.frame(data@meta.data),rowData=as.data.frame(data@assays$RNA@meta.features))
data=list(data=.mySplitObject(data,'batch_merging'),data_m=NULL)
#cellClass="Interneuronhuman";inputData=data
.myAnalysisFn(cellClass="Interneuronhuman",inputData=data)

####################
#Arranging the data
rm(list=ls())
files=c("astro_anno.qs","Interneurons_anno.qs","MG_anno.qs","Oligo_anno.qs","SPN_anno.qs")
anno=NULL
for(ifile in files){
  tmp=qread(file.path("~/myBucket/striatum_human/",ifile))
  anno=rbind(anno,tmp)
}

data_m=qread("~/myBucket/striatum_human/seurat4_obj.qs")


anno=anno[,c("anno_cluster_res","UMAP_1","UMAP_2","anno_cluster_label")]
colnames(anno)=paste0("subcluster_",colnames(anno))
anno$sample=row.names(anno)

pd=as.data.frame(data_m@meta.data)
pd$sample=row.names(pd)
anno=merge(pd,anno,by="sample",all=T)

anno$anno_cluster_label[!is.na(anno$subcluster_anno_cluster_label)]=anno$subcluster_anno_cluster_label[!is.na(anno$subcluster_anno_cluster_label)]
anno=anno[match(row.names(pd),anno$sample),]

data_m@meta.data=anno

data=list(countData=data_m@assays$RNA@counts,pdata=as.data.frame(data_m@meta.data),umapdata=data_m@reductions$umap@cell.embeddings,embedding=data_m@reductions$pca@cell.embeddings,fdata=as.data.frame(data_m@assays$RNA@meta.features))
qsave(data,file="~/myBucket/striatum_human/final_annotated_obj.qs")


t1 <- aggregate(nUMI ~ dataset,orig.striatum.use@cell.data,median)
t2 <- aggregate(nGene ~ dataset,orig.striatum.use@cell.data,median)
qsave(list(t1,t2),'/home/tkamath/DA/striatum/striatum_qc.qs')


# Cutoffs for libraries
libs.use <- grep('CN',list.files('/home/tkamath/data/'),value = T)
libs.use <-paste0('/home/tkamath/data/',libs.use)
out.cn <- rliger::read10X(sample.dirs = libs.use,min.umis = c(1000,1000,1000,1000,1000,1000,1000,1000),
                          sample.names = c('5730DAPI','5730NeuN','3839DAPIA','3839DAPIB',
                                           '3898DAPIA','3898DAPIB','4340DAPIA','4340DAPIB'))
tmp1 <- CreateSeuratObject(out.cn)
tmp1 <- AddMito(tmp1,species = 'human')
orig.striatum.use@cell.data$percent.mito <- tmp1@meta.data[match(rownames(orig.striatum.use@cell.data),
                                                                 rownames(tmp1@meta.data) ),]$percent.mito

striatum.obj <- qread('../final_annotated_obj.qs')

striatum.seurat <- make.seurat(striatum.obj)

striatum.seurat@ident <- as.factor(striatum.obj$pdata$anno_cluster_label)
names(striatum.seurat@ident) <- rownames(striatum.seurat@meta.data)

idx.use <- names(striatum.seurat@ident[which(!grepl('doublet',as.character(striatum.seurat@ident)))])
striatum.seurat <- SubsetData(striatum.seurat, cells.use = idx.use,subset.raw = T)

striatum.liger <- rliger::seuratToLiger(striatum.seurat,combined.seurat = T,meta.var = 'batch_merging')
striatum.liger@tsne.coords <- as.matrix(striatum.seurat@dr$umap@cell.embeddings)
striatum.liger@H.norm <- striatum.seurat@dr$hpca@cell.embeddings
striatum.liger@W <- as.matrix(striatum.seurat@dr$hpca@gene.loadings)
striatum.liger@H <- lapply(names(striatum.liger@raw.data),function(x){
  as.matrix(striatum.liger@H.norm[which(grepl(x,rownames(striatum.liger@H.norm))),])
})
names(striatum.liger@H) <- names(striatum.liger@raw.data)
qsave(striatum.liger,'striatum_object_cleaned.qs')
levels(striatum.liger@clusters) <- c('Astro','Astro','Astro','Endofibro','DRD1','DRD1',
                                     'DRD1','DRD1','DRD1','DRD1','DRD2','DRD2',
                                     'DRD2','DRD2','DRD2','DRD2','DRD2','Endofibro',
                                     'Interneuron','Interneuron','Interneuron','Interneuron','Interneuron','Interneuron',
                                     'Interneuron','Interneuron','Interneuron','Interneuron','Microglia/macrophage','Microglia/macrophage',
                                     'Microglia/macrophage','Microglia/macrophage','Microglia/macrophage','Microglia/macrophage','Microglia/macrophage',
                                     'OPC','Olig','Olig','Olig','Olig','Olig','Olig','Olig','Olig', 'Interneuron')
# Generate wilcoxon cluster markers
striatum.seurat@var.genes <- striatum.obj$fdata
striatum.liger@var.genes <- striatum.seurat@var.genes
striatum.seurat@ident <- striatum.liger@clusters

rownames(striatum.seurat@data) <- striatum.obj$fdata$symbol

orig.striatum <- readRDS('allstriatum.liger.rds')
orig.raw <- orig.striatum@raw.data
orig.striatum <- createLiger(orig.raw)
idx.use <- as.character(unlist(lapply(strsplit(rownames(striatum.seurat@meta.data),'_'),function(x){paste0(x[[2]],'_',x[[3]])})))
orig.striatum.use <- rliger::subsetLiger(orig.striatum, cells.use = idx.use,remove.missing = F)
orig.striatum.use@cell.data <- orig.striatum@cell.data[idx.use,]
orig.striatum.use <- AddMito(orig.striatum.use, species = 'human')
names(striatum.seurat@ident) <- as.character(unlist(lapply(strsplit(names(striatum.seurat@ident),'_'),
                                                           function(x){paste0(x[[2]],"_",x[[3]])})))
orig.striatum.use@clusters <- striatum.seurat@ident
orig.striatum.use@cell.data$clusters <- orig.striatum.use@clusters

rownames(striatum.seurat@meta.data) <- idx.use
orig.striatum.use@clusters <- striatum.seurat@ident[match(rownames(orig.striatum.use@cell.data),
                                                          rownames(striatum.seurat@meta.data))]
names(orig.striatum.use@clusters) <- rownames(orig.striatum.use@cell.data)
orig.striatum.use@tsne.coords <- as.matrix(striatum.seurat@dr$umap@cell.embeddings)
orig.striatum.use@H.norm <- striatum.seurat@dr$hpca@cell.embeddings
orig.striatum.use@W <- as.matrix(striatum.seurat@dr$hpca@gene.loadings)
orig.striatum.use@H <- lapply(names(orig.striatum.use@raw.data),function(x){
  as.matrix(orig.striatum.use@H.norm[which(grepl(x,rownames(orig.striatum.use@H.norm))),])
})
names(orig.striatum.use@H) <- names(orig.striatum.use@raw.data)
orig.striatum.use@cell.data$clusters <- orig.striatum.use@clusters
qsave(orig.striatum.use,'/home/tkamath/DA/striatum/striatum_ligeruse.qs')

qsave(striatum.seurat,'/home/tkamath/DA/striatum/striatum_seurat.qs')
striatum.markers <- prestowrapper(striatum.seurat,log.FC = log(1.2),all.clusters = T,one.sided = T,use.raw = F)

lapply(unique(striatum.markers$group),function(z){
  da.wilcox <- striatum.markers[striatum.markers$group == z,]
  df.print <- data.frame('gene' = da.wilcox[order(da.wilcox$auc,decreasing = T),]$feature[1:2000],
                         'set' = rep('da',2000))
  df.print$gene <- as.character(mapIds(org.Hs.eg.db, keys = as.character(df.print$gene),keytype = 'SYMBOL',column = 'ENSEMBL'))
  if (z == 'Microglia/macrophage'){
    write.table(df.print,paste0('/home/tkamath/DA/striatum/mggeneset.csv'),
                quote = F,row.names = F,col.names = F,sep = '\t')
  }
  else
    write.table(df.print,paste0('/home/tkamath/DA/striatum/',as.character(z),'geneset.csv'),
                quote = F,row.names = F,col.names = F,sep = '\t')
})  

celltypes <- unique(striatum.markers$group)
celltypes[6] <- 'mg'
out.pvals <- mclapply(unique(striatum.markers$group),function(z){
  setwd('/home/tkamath/magma/')
  system(paste0('./magma --gene-results PDFUMAresults/magma.genes.raw --set-annot /home/tkamath/DA/striatum/',
                as.character(z),'geneset.csv gene-col=1 set-col=2 --out pd',as.character(z)))
},mc.cores = 60)

df.use <- data.frame('celltype' = celltypes, 'magmap' = as.vector(unlist(lapply(celltypes,
                                                                                function(x){
                                                                                  magma.use <- read.table(paste0('/home/tkamath/magma/pd',as.character(x),'.gsa.out'),header = T)
                                                                                  return(-log10(magma.use$P))
                                                                                }))))

# Now find markers using MAST
orig.striatum.use <- qread('/home/tkamath/DA/striatum/striatum_ligeruse.qs')
object.use <- orig.striatum.use
data.use <- MergeSparseDataAll(object.use@raw.data)
celldata.use <- orig.striatum.use@cell.data
data.use <- data.use[,rownames(celldata.use)]
data.use <- log1p(10000*sweep(data.use,MARGIN=2,FUN="/",STATS=Matrix::colSums(data.use)))
celltypes <- unique(celldata.use$clusters)
celldata.use$lognumi <- log(celldata.use$nUMI)
celldata.use$lib <- as.factor(unlist(sapply(strsplit(rownames(celldata.use),'_'),`[`,1)))

lapply(celltypes,function(z){
  data.use2 <- data.use[which(Matrix::rowMeans(data.use[,rownames(celldata.use[
    celldata.use$cluster == z,])] > 0) > 0.03),]
  sca.test3 <- FromMatrix(as.matrix(data.use2 ))   
  colData(sca.test3) <- cbind(colData(sca.test3),celldata.use)
  colData(sca.test3)$clusteruse <- as.factor(ifelse(colData(sca.test3)$clusters == z,'1','0'))
  system.time(summary.test1 <- zlm(sca = sca.test3,formula = as.formula
                                   ('~ lognumi + clusteruse + percent.mito + (1|lib)'),
                                   method  = 'glmer',ebayes = F,parallel = T,fitArgsD = list(nAGQ = 0)))
  summaryCond <- summary(summary.test1, doLRT= 'clusteruse1') 
  summaryDt <- summaryCond$datatable
  if (z == 'Microglia/macrophage'){
    qsave(summaryDt,paste0('/home/tkamath/DA/striatum/mastDE_mg.qs'))
  }
  else
    qsave(summaryDt,paste0('/home/tkamath/DA/striatum/mastDE_',as.character(z),'.qs'))
})

make.seurat<-function(data){
  nbt=Seurat::CreateSeuratObject(raw.data=data$countData,
                                 project = "SeuratProject",
                                 min.cells = 0,
                                 min.features = 0,
                                 names.field = 1,
                                 names.delim = "-",
                                 meta.data = data$pdata)
  nbt@ident <- data$pdata$anno_cluster_res
  names(nbt@ident) <- rownames(data$pdata)
  nbt = SetDimReduction(object=nbt,reduction.type='umap',
                        slot='cell.embeddings',new.data = as.matrix(data$pdata[,c('UMAP_1','UMAP_2')]))
  nbt = SetDimReduction(nbt,reduction.type='umap',slot='key',new.data = 'UMAP_')
  
  hpca.obj<-new(Class="dim.reduction",
                cell.embeddings=data$embedding,key="HPCA",misc=list("raw" = data$embedding))
  colnames(hpca.obj@cell.embeddings) = paste0("Factor",1:ncol(data$embedding))
  nbt@dr$hpca <- hpca.obj
  
  return(nbt)
}


list.out <- qread('/home/tkamath/DA/striatum/striatum_qc.qs')
t1 <- list.out[[1]]
t2 <- list.out[[2]]
pdf('/home/tkamath/DA/figures/suppfigures/suppfig7/qcbydataset_striatum.pdf',useDingbats = F,
    width = 8,height = 5)
ggplot(t1, aes(x = dataset, y = nUMI,fill = dataset)) + geom_bar(stat = 'identity') +
  theme(legend.position = 'none')

ggplot(t2, aes(x = dataset, y = nGene,fill = dataset)) + geom_bar(stat = 'identity') +
  theme(legend.position = 'none')
dev.off()

orig.striatum.use <- qread('/home/tkamath/DA/striatum/striatum_ligeruse.qs')
orig.striatum.use@cell.data$clusters <- as.factor(orig.striatum.use@cell.data$clusters)
orig.striatum.use@cell.data$clusters <- factor(orig.striatum.use@cell.data$clusters,
                                               levels = c('DRD1','DRD2','Interneuron',
                                                          'Astro','OPC','Microglia/macrophage',
                                                          'Endofibro','Olig'))

pdf('/home/tkamath/DA/figures/suppfigures/suppfig7/qcbycell_striatum.pdf',useDingbats = F,
    width = 10,height = 5)
ggplot(orig.striatum.use@cell.data,aes(x = clusters, y = nUMI,fill = clusters)) + geom_violin() +
  scale_y_log10() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  geom_boxplot(width = 0.2,fill = 'white') +
  theme(legend.position = 'none') + scale_fill_manual(values = mycolors)

ggplot(orig.striatum.use@cell.data,aes(x = clusters, y = nGene,fill = clusters)) + geom_violin() +
  scale_y_log10() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  geom_boxplot(width = 0.2,fill = 'white') +
  theme(legend.position = 'none') + scale_fill_manual(values = mycolors)
dev.off()

t1 <- table(orig.striatum.use@cell.data$dataset, orig.striatum.use@clusters)
t1 <- as.data.frame(melt(t1/rowSums(t1)))
t1$Var1 <- as.factor(t1$Var1)

t1$Var2 <- factor(t1$Var2,levels =c('DRD1','DRD2','Interneuron',
                                    'Astro','OPC','Microglia/macrophage',
                                    'Endofibro','Olig'))

pdf('/home/tkamath/DA/figures/suppfigures/suppfig7/clusterbyindividual.pdf',useDingbats = F,
    height = 8,width = 7)
ggplot(t1, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

orig.striatum.use@clusters <- factor(orig.striatum.use@clusters,levels = c('DRD1','DRD2','Interneuron',
                                                                           'Astro','OPC','Microglia/macrophage',
                                                                           'Endofibro','Olig'))
p2 <- plotByDatasetAndCluster(orig.striatum.use,return.plots = T,do.legend = F)
p1[[2]]
p1 <- plotByDatasetAndCluster(orig.striatum.use,return.plots = T,do.legend = F,text.size = 0)
p2[[2]]

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)
p1.striatum <- p1[[2]] + scale_color_manual(values = mycolors) + 
  geom_label_repel(data = p2[[2]]$layers[[2]]$data,aes(x = tsne1, y = tsne2,label = Cluster),size = 5,
                   label.size = 0.5,label.r = 0,inherit.aes = F)+ 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.border = element_blank())
pdf('/home/tkamath/DA/figures/suppfigures/suppfig7/tsne_striatum.pdf',useDingbats = F)
p1.striatum
dev.off()

orig.striatum.seurat <- CreateSeuratObject(raw.data = MergeSparseDataAll(orig.striatum.use@raw.data))
orig.striatum.seurat@ident <- orig.striatum.use@clusters[colnames(orig.striatum.seurat@data)]

markers.striatum <- prestowrapper(orig.striatum.seurat,all.clusters = T)
orig.striatum.seurat@ident <- factor(orig.striatum.seurat@ident,levels = c('Olig','Endofibro',
                                                                           'Microglia/macrophage',
                                                                           'OPC','Astro','Interneuron','DRD2','DRD1'))
p1 <- DotPlot(orig.striatum.seurat, genes.plot = c('DRD1','DRD2' ,'GAD1','GJA1','CSPG4', 'CX3CR1',
                                                   'DCN', 'MBP'),do.return = T,x.lab.rot = T,plot.legend = T)
pdf('/home/tkamath/DA/figures/suppfigures/suppfig7/markergenes_striatum.pdf',useDingbats = F)
p1 + coord_flip()
dev.off()
