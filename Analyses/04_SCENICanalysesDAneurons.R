############################################
## Human and macaque SCENIC analyses #######
############################################

#Running SCENIC on pooled cells (20) as per Suo et al. (2018)
#Load the required libraries
library(rlang) #Must be loaded first 
library(BiocManager)
library(devtools)
library(dplyr)
library(rliger)
library(Seurat)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCENIC)
library(shiny)
library(pheatmap)
library(qs)
library(data.table)
library(reshape2)
library(stringr)
#Set seed for randomizations 
set.seed(1)
#Load in the dataset, convert to Seurat object and pull the normalized data 
ligerex<-qread("/home/tkamath/DA/da/sn_da_annot_021821.qs")
idx.use <- rownames(ligerex@cell.data[which(ligerex@cell.data$status == 'Ctrl'),])
ligerex.ctrl <- subsetLiger(ligerex,cells.use = idx.use,remove.missing = F)
object.2<-ligerToSeurat(ligerex.ctrl) 
exprMat<-as.array(GetAssayData(object.2, slot = 'data')) #Note this is normalized data 
cellInfo<-data.frame(object.2@ident)#Pull the cluster information from the Seurat object (barcode and cluster ids) 
cellInfo$bc<-rownames(cellInfo)
num_clusters<-length(unique(cellInfo$object.2.ident)) #Determine number of cluters 
#Subset a matrix to test this with
#exprMat<-exprMat[1:1000,1:1000, drop=F]
#cellInfo<-cellInfo[1:1000,]
#Create the pooled experimental matrix for regulon determination
exprMat.pooled<-matrix(ncol=1, nrow=nrow(exprMat)) #Need to set number of rows 
#Cycle through each cluster and perform the pooling 
for(c in 1:num_clusters) {
  cluster<-unique(cellInfo$object.2.ident)[c] #Set the cluster for analysis 
  num_cells<-nrow(filter(cellInfo, object.2.ident==cluster))
  bcs_to_pool<-filter(cellInfo, object.2.ident==cluster)$bc #Pull just the barcodes 
  submat<-exprMat[,bcs_to_pool] #Pull cells from each cluster from exprMatrix
  num_iterations<-floor(num_cells/10)  #Determine number of pooling iterations to run 
  random.submat<-submat[,sample(ncol(submat), replace = FALSE)] #Shuffle the columns of the subsetted matrix (submat)
  #Initialize a new expression matrix for this cluster, and cycle through columns pooling the cells. 
  cluster.mat<-matrix(nrow = nrow(exprMat), ncol=num_iterations)
  colnames(cluster.mat)<-rep(as.character(cluster), num_iterations) #Need to name the cells something to initilize
  for(i in 1:num_iterations) { 
    stop<-10*i
    start<-stop-9
    avg<-as.matrix(rowMeans(random.submat[,start:stop]))  #Calculate the rowmeans for pooled cells 
    avg.name<-paste(colnames(random.submat[,start:stop]), collapse = "Z") #Paste together the names  
    cluster.mat[,i]<-avg
    rownames(cluster.mat)<-rownames(avg) #Keep the gene names
    colnames(cluster.mat)[i]<-paste0(cluster, "_", avg.name) #Include the original cluster ID in the pooled cell name also
  }
  exprMat.pooled<-cbind(exprMat.pooled, cluster.mat)
  print(cluster)
}
exprMat.pooled<-exprMat.pooled[,-1] #Remove the column needed to initialize the matrix and check numbers again 
#Save this pooled matrix 
qsave(exprMat.pooled, "/home/tkamath/DA/tfanalysis/pooledexpression_da.qs")

#Proceed with the SCENIC pipeline using this pooled matrix to create the regulons 
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

## File names after download
downloadedDBFiles <- c("cisTarget_databases/hg19-500bp-upstream-7species.mc9nr.feather",
                       "cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.feather")

## The databases are downloaded only in the case when the files do not exist.
if( !(file.exists(downloadedDBFiles[1]) & file.exists(downloadedDBFiles[2])) ){
  
  ## Create and move to the new folder
  dir.create("cisTarget_databases");
  setwd("cisTarget_databases") 
  for(featherURL in dbFiles)
  {
    download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
    descrURL <- gsub(".feather$", ".descr", featherURL)
    if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
  }
  setwd("..")
  
}


## Config object
org="hgnc" ## Set to Homo sapiens
dbDir="databases" # RcisTarget databases location
scenicOptions <- initializeScenic(org = org,datasetTitle = 'test1',nCores = 50,
                                  dbDir = '/home/tkamath/DA/tfanalysis/cisTarget_databases/')

#Filter (already done in Seurat so it's set to 0) and compare with Rcistarget databases 
genesKept <- geneFiltering(exprMat.pooled, scenicOptions=scenicOptions
                           ,minCountsPerGene=0
                           ,minSamples=0)

interestingGenes <- c("LMX1A",'NR4A2','SOX6','OTX2')
interestingGenes[which(!interestingGenes %in% genesKept)]
#Filter the expression matrix to contain only these genes
exprMat_filtered <- exprMat.pooled[genesKept, ]
dim(exprMat_filtered)
#Before running GENIE3, run the correlation to determine the coexpression sign
runCorrelation(exprMat_filtered, scenicOptions)
#Run GENIE3 (normalize your data if necessary), takes ~3 days
dim(exprMat_filtered)
runGenie3(exprMat_filtered, scenicOptions)
#Run SCENIC on pooled cells using the wrappers and examine the output of each step (should no longer scale with num cells)
runSCENIC_1_coexNetwork2modules(scenicOptions)
qsave(scenicOptions,"scenicoptions.qs") 
#Check the tf modules (and the threshold/method) that come out before pruning with Rcistarget. 
tfmodules_preStep2<-readRDS(file = "int/1.6_tfModules_asDF.Rds")
View(tfmodules_preStep2)
#Next run step 2 to prune the tf modules using RcisTarget
runSCENIC_2_createRegulons(scenicOptions) 
#Score all the cells with AUCell (not just pooled cells)
runSCENIC_3_scoreCells(scenicOptions, exprMat)
#Heatmap of all the cells, regulon -x- liger cluster 
regulonAUC<-loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC<-regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulon.presto <- wilcoxauc(as.matrix(regulonAUC@assays$data@listData$AUC),y = cellInfo$object.2.ident)
regulon.supp <- regulon.presto[regulon.presto$padj < 0.01,]$feature


order.da <- rev(c('SOX6_AGTR1','SOX6_PART1',
                  'SOX6_DDT','SOX6_GFRA2',
                  'CALB1_CALCR','CALB1_CRYM_CCDC68','CALB1_PPP1R17',
                  'CALB1_RBP4','CALB1_GEM','CALB1_TRHR'))
cellInfo$object.2.ident <- factor(cellInfo$object.2.ident,levels = rev(order.da))
regulonActivity_byCellType<-sapply(split(rownames(cellInfo), cellInfo$object.2.ident),
                                   function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

regulon.presto <- split(regulon.presto,regulon.presto$group)
regulons.plot <- rev(unique(as.vector(unlist(lapply(rev(order.da),function(x){
  regulon.presto[[x]][order(regulon.presto[[x]]$auc,decreasing = T)[1:3],]$feature
})))))

#Filter and remove lowly expressing regulons 
minVal<-0.05
regulonActivity_byCellType2<-regulonActivity_byCellType #Create a new object to not disrupt downstream analysis
regulonActivity_byCellType2<-regulonActivity_byCellType2[which(rowSums(regulonActivity_byCellType2>minVal)>0),]

regulonActivity_byCellType_Scaled<-t(scale(t(regulonActivity_byCellType), center = T, scale=T))
to.plot.df <- melt(to.plot)
to.plot <- regulonActivity_byCellType_Scaled[regulons.plot,]
colnames(to.plot.df) <- c('Regulon','DA_type','value')
to.plot.df$Regulon <- as.factor(to.plot.df$Regulon)
to.plot.df$Regulonshow <- sapply(lapply(as.character(to.plot.df$Regulon),
                                        function(x){strsplit(x," ")[[1]]}),`[`,1)
to.plot.df$Regulonshow <- as.factor(unlist(sapply(strsplit(to.plot.df$Regulonshow,'_'),`[`,1)))
to.plot.df$Regulonshow <- as.factor(to.plot.df$Regulonshow)
qsave(to.plot.df,'/home/tkamath/DA/tfanalysis/tf_toplot.qs')

to.plot <- regulonActivity_byCellType_Scaled[regulon.supp,]
rownames(to.plot) <- sapply(lapply(as.character(rownames(to.plot)),
                                   function(x){strsplit(x," ")[[1]]}),`[`,1)
rownames(to.plot) <- as.factor(unlist(sapply(strsplit(rownames(to.plot),'_'),`[`,1)))
p1 <- pheatmap::pheatmap(to.plot)
ord <- p1$tree_row$order
to.plot.df <- melt(to.plot)
colnames(to.plot.df) <- c('Regulon','DA_type','value')
to.plot.df$Regulon <- reorder(to.plot.df$Regulon,new.order = ord)
qsave(to.plot.df,'/home/tkamath/DA/tfanalysis/tf_forsupp.qs')

#Make unfiltered dataframe 
regulonActivity_byCellType_Scaled_df<-as.data.frame(regulonActivity_byCellType_Scaled)
regulonActivity_byCellType_Scaled_df$regulon<-rownames(regulonActivity_byCellType_Scaled_df)
#Plot same heatmap by filtering out negative or lowly expressed regulons (note the value is scaled)
minVal<-2
regulonActivity_byCellType_Scaled<-regulonActivity_byCellType_Scaled[which(rowSums(regulonActivity_byCellType_Scaled>=minVal)>0),]
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,
                   fontsize_row=3)
#Future plans to filter these regulons 
#1) Filter out weakly expressing regulons (unscaled)
#2) Filter out extended motifs based on motif similarity (weak evidence)
#2) Filter out small or useless clusters 
#3) Filter by TF expression counts 
#4) Filter by number of target genes 
#In preparation for binarizaition: Now adjust the AUC threshold based on data.  
aucellApp<-plotTsne_AUCellApp(scenicOptions, exprMat)
savedSelections<-shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#Binarize the regulons per cell based on an AUC threshold
runSCENIC_4_aucell_binarize(scenicOptions)
#Create heatmap based on binarized regulons,plotting percentage of cells in that cluster with regulon active (1)
#minPerc let's you cut off weakly expressing regulons (although thresholding might have set regulons of interest to 0)
minPerc <- 0
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$object.2.ident), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5, 
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,
                   fontsize_row=3)
#Save the SCENIC options file 
qsave(scenicOptions, file="int/scenicOptions.qs") 


### Subtype analysis with macaque dataset
#Load in the dataset, convert to Seurat object and pull the normalized data 
setwd('~/DA/tfanalysis/macaque_SCENIC/')
ligerex<-qread("/home/tkamath/DA/species/newmacaque/newmacaque_human.qs")
macaque.idx <- rownames(ligerex@cell.data[which(
  ligerex@cell.data$dataset == 'macaque'),])
ligerex <- subsetLiger(ligerex,cells.use =macaque.idx,remove.missing = F )
ligerex <- ligerToSeurat(ligerex)
exprMat<-as.array(GetAssayData(ligerex, slot = 'data')) #Note this is normalized data 
cellInfo<-data.frame(ligerex@ident)#Pull the cluster information from the Seurat object (barcode and clustera ids) 
cellInfo$bc<-rownames(cellInfo)
num_clusters<-length(unique(cellInfo$ligerex.ident)) #Determine number of cluters 
exprMat.pooled<-matrix(ncol=1, nrow=nrow(exprMat)) #Need to set number of rows 
#Cycle through each cluster and perform the pooling 
for(c in 1:num_clusters) {
  cluster<-unique(cellInfo$ligerex.ident)[c] #Set the cluster for analysis 
  num_cells<-nrow(dplyr::filter(cellInfo, ligerex.ident==cluster))
  bcs_to_pool<-dplyr::filter(cellInfo, ligerex.ident==cluster)$bc #Pull just the barcodes 
  submat<-exprMat[,bcs_to_pool] #Pull cells from each cluster from exprMatrix
  num_iterations<-floor(num_cells/10)  #Determine number of pooling iterations to run 
  random.submat<-submat[,sample(ncol(submat), replace = FALSE)] #Shuffle the columns of the subsetted matrix (submat)
  #Initialize a new expression matrix for this cluster, and cycle through columns pooling the cells. 
  cluster.mat<-matrix(nrow = nrow(exprMat), ncol=num_iterations)
  colnames(cluster.mat)<-rep(as.character(cluster), num_iterations) #Need to name the cells something to initilize
  for(i in 1:num_iterations) { 
    stop<-10*i
    start<-stop-9
    avg<-as.matrix(rowMeans(random.submat[,start:stop]))  #Calculate the rowmeans for pooled cells 
    avg.name<-paste(colnames(random.submat[,start:stop]), collapse = "Z") #Paste together the names  
    cluster.mat[,i]<-avg
    rownames(cluster.mat)<-rownames(avg) #Keep the gene names
    colnames(cluster.mat)[i]<-paste0(cluster, "_", avg.name) #Include the original cluster ID in the pooled cell name also
  }
  exprMat.pooled<-cbind(exprMat.pooled, cluster.mat)
  print(cluster)
}
exprMat.pooled<-exprMat.pooled[,-1] #Remove the column needed to initialize the matrix and check numbers again 
#Save this pooled matrix 
qsave(exprMat.pooled, "/home/tkamath/DA/tfanalysis/pooledexpression_da_macaque.qs")
qsave(exprMat, "/home/tkamath/DA/tfanalysis/expression_da_macaque.qs")

#Proceed with the SCENIC pipeline using this pooled matrix to create the regulons 
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

## File names after download
downloadedDBFiles <- c("cisTarget_databases/hg19-500bp-upstream-7species.mc9nr.feather",
                       "cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.feather")

## The databases are downloaded only in the case when the files do not exist.
if( !(file.exists(downloadedDBFiles[1]) & file.exists(downloadedDBFiles[2])) ){
  
  ## Create and move to the new folder
  dir.create("cisTarget_databases");
  setwd("cisTarget_databases") 
  for(featherURL in dbFiles)
  {
    download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
    descrURL <- gsub(".feather$", ".descr", featherURL)
    if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
  }
  setwd("..")
  
}


## Config object
org="hgnc" ## Set to Homo sapiens
dbDir="databases" # RcisTarget databases location
scenicOptions <- initializeScenic(org = org,datasetTitle = 'scenice_macaque',nCores = 15,
                                  dbDir = '/home/tkamath/DA/tfanalysis/cisTarget_databases/')

#Filter (already done in Seurat so it's set to 0) and compare with Rcistarget databases 
genesKept <- geneFiltering(exprMat.pooled, scenicOptions=scenicOptions
                           ,minCountsPerGene=0
                           ,minSamples=0)

interestingGenes <- c("LMX1A",'NR4A2','SOX6','OTX2')
interestingGenes[which(interestingGenes %in% genesKept)]
#Filter the expression matrix to contain only these genes
exprMat_filtered <- exprMat.pooled[genesKept, ]
dim(exprMat_filtered)
#Before running GENIE3, run the correlation to determine the coexpression sign
runCorrelation(exprMat_filtered, scenicOptions)
#Run GENIE3 (normalize your data if necessary), takes ~3 days
dim(exprMat_filtered)
runGenie3(exprMat_filtered, scenicOptions)
#Run SCENIC on pooled cells using the wrappers and examine the output of each step (should no longer scale with num cells)
runSCENIC_1_coexNetwork2modules(scenicOptions)
qsave(scenicOptions,"scenicoptions_macaque.qs") 
#Check the tf modules (and the threshold/method) that come out before pruning with Rcistarget. 
tfmodules_preStep2<-readRDS(file = "~/DA/tfanalysis/int/1.6_tfModules_asDF.Rds")
View(tfmodules_preStep2)
#Next run step 2 to prune the tf modules using RcisTarget
setwd('../')
runSCENIC_2_createRegulons(scenicOptions) 
#Score all the cells with AUCell (not just pooled cells)
runSCENIC_3_scoreCells(scenicOptions, exprMat)
#Heatmap of all the cells, regulon -x- liger cluster 
regulonAUC<-loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC<-regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

qsave(regulonAUC,'~/DA/tfanalysis/macaque_SCENIC/macaque_AUCdata.qs')

regulon.presto <- wilcoxauc(as.matrix(regulonAUC@assays$data@listData$AUC),y = cellInfo$ligerex.ident)
regulon.supp <- regulon.presto[regulon.presto$padj < 0.01,]
qsave(regulon.supp,'~/DA/tfanalysis/tf_forsupp_macaque.qs')

order.da <- rev(c('SOX6_AGTR1','SOX6_PART1',
                  'SOX6_DDT','SOX6_GFRA2',
                  'CALB1_CALCR','CALB1_CRYM_CCDC68','CALB1_PPP1R17',
                  'CALB1_RBP4','CALB1_GEM','CALB1_TRHR'))
cellInfo$ligerex.ident <- factor(cellInfo$ligerex.ident,levels = rev(order.da))
regulonActivity_byCellType<-sapply(split(rownames(cellInfo), cellInfo$ligerex.ident),
                                   function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

regulon.presto <- split(regulon.presto,regulon.presto$group)
regulons.plot <- rev(unique(as.vector(unlist(lapply(rev(order.da),function(x){
  regulon.presto[[x]][order(regulon.presto[[x]]$auc,decreasing = T)[1:3],]$feature
})))))

#Filter and remove lowly expressing regulons 
minVal<-0.05
regulonActivity_byCellType2<-regulonActivity_byCellType #Create a new object to not disrupt downstream analysis
regulonActivity_byCellType2<-regulonActivity_byCellType2[which(rowSums(regulonActivity_byCellType2>minVal)>0),]

regulonActivity_byCellType_Scaled<-t(scale(t(regulonActivity_byCellType), center = T, scale=T))
to.plot <- regulonActivity_byCellType_Scaled[regulons.plot,]
to.plot.df <- melt(to.plot)
colnames(to.plot.df) <- c('Regulon','DA_type','value')
to.plot.df$Regulon <- as.factor(to.plot.df$Regulon)
to.plot.df$Regulonshow <- sapply(lapply(as.character(to.plot.df$Regulon),
                                        function(x){strsplit(x," ")[[1]]}),`[`,1)
to.plot.df$Regulonshow <- as.factor(unlist(sapply(strsplit(to.plot.df$Regulonshow,'_'),`[`,1)))
to.plot.df$Regulonshow <- as.factor(to.plot.df$Regulonshow)
qsave(to.plot.df,'/home/tkamath/DA/tfanalysis/tf_toplot_macaque.qs')

to.plot <- regulonActivity_byCellType_Scaled[regulon.supp,]
rownames(to.plot) <- sapply(lapply(as.character(rownames(to.plot)),
                                   function(x){strsplit(x," ")[[1]]}),`[`,1)
rownames(to.plot) <- as.factor(unlist(sapply(strsplit(rownames(to.plot),'_'),`[`,1)))
p1 <- pheatmap::pheatmap(to.plot)
ord <- p1$tree_row$order
to.plot.df <- melt(to.plot)
colnames(to.plot.df) <- c('Regulon','DA_type','value')
to.plot.df$Regulon <- reorder(to.plot.df$Regulon,new.order = ord)
qsave(to.plot.df,'/home/tkamath/DA/tfanalysis/tf_forsupp.qs')

#Make unfiltered dataframe 
regulonActivity_byCellType_Scaled_df<-as.data.frame(regulonActivity_byCellType_Scaled)
regulonActivity_byCellType_Scaled_df$regulon<-rownames(regulonActivity_byCellType_Scaled_df)
#Plot same heatmap by filtering out negative or lowly expressed regulons (note the value is scaled)
minVal<-2
regulonActivity_byCellType_Scaled<-regulonActivity_byCellType_Scaled[which(rowSums(regulonActivity_byCellType_Scaled>=minVal)>0),]
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,
                   fontsize_row=3)
#Future plans to filter these regulons 
#1) Filter out weakly expressing regulons (unscaled)
#2) Filter out extended motifs based on motif similarity (weak evidence)
#2) Filter out small or useless clusters 
#3) Filter by TF expression counts 
#4) Filter by number of target genes 
#In preparation for binarizaition: Now adjust the AUC threshold based on data.  
aucellApp<-plotTsne_AUCellApp(scenicOptions, exprMat)
savedSelections<-shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#Binarize the regulons per cell based on an AUC threshold
runSCENIC_4_aucell_binarize(scenicOptions)
#Create heatmap based on binarized regulons,plotting percentage of cells in that cluster with regulon active (1)
#minPerc let's you cut off weakly expressing regulons (although thresholding might have set regulons of interest to 0)
minPerc <- 0
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$object.2.ident), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5, 
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,
                   fontsize_row=3)
#Save the SCENIC options file 
qsave(scenicOptions, file="macaque_SCENIC/scenicoptions_macaque.qs") 

