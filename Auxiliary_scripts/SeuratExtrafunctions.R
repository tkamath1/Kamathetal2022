

#A simple function to convert the weird 10x format into a sparse matrix for loading into Seurat.
Convert10xtoMatrix <-function(input){
  
  raw.data = exprs(input)
  genes = fData(input)
  idx.unique =which(table(genes[,2])==1)
  genes.use = names(table(genes[,2]))[idx.unique]
  idx = match(genes.use,genes[,2])
  raw.data = raw.data[idx,]
  rownames(raw.data) = genes.use
  return(raw.data)  
}

scalar1 <- function(x) {x / sqrt(sum(x^2))}


#I prefer how we select variable genes to Seurat's approach, which has some weird binning and normalization
mean.var.auto.seurat = function(object, alphathresh=0.99,varthresh=0.1,cex.use=0.3) {
  dge = object@raw.data
  xm = sweep(dge,2,colSums(dge),"/")
  trx_per_cell <- colSums(dge)
  
  gene_expr_mean <- rowMeans(xm)  # Each gene's mean expression level (across all cells)
  gene_expr_var  <- apply( xm, 1, var  )  # Each gene's expression variance (across all cells)
  nolan_constant <- mean((1 / trx_per_cell )) 
  alphathresh.corrected=alphathresh/dim(dge)[1]
  genemeanupper <- gene_expr_mean+qnorm(1-alphathresh.corrected/2)*sqrt(gene_expr_mean*nolan_constant/dim(dge)[2])
  genes.use=names(gene_expr_var)[which(gene_expr_var/nolan_constant> genemeanupper & log10(gene_expr_var) > log10(gene_expr_mean)+(log10(nolan_constant)+varthresh))]
  plot( log10(gene_expr_mean), log10(gene_expr_var), cex=cex.use
  )
  
  points(log10(gene_expr_mean[genes.use]),log10(gene_expr_var[genes.use]),cex=cex.use,col='green')
  abline(log10(nolan_constant),1,col='purple')
  
  legend("bottomright",paste0("Selected genes: ",length(genes.use)),pch=20,col='green')
  return(genes.use)
}


environment(mean.var.auto.seurat)<-asNamespace('Seurat')


#The main curation code for ICA.  Plots gene loadings and cell embeddings for each IC, and displays highest loading genes.  Helps determine whether an IC relates to a doublet class, an artifact, or is more likely to be a real biological signal.

#Intention is to run this twice for each ICA: run once with curation.data = NULL, which displays all ICs and draws a tSNE from all ICs; then, use the generated PDF to curate which ICs to retain, and feed in that curation data into the next function call.

CurateICA.seurat<-function(object,make.tsne = T, feature.plot=TRUE,kurtosis_threshold = 6,skewness_threshold = 1,curation.data=NULL,cluster.doublets = T,doubletalpha=0.01) {
  
  x = list(object@dr$ica@gene.loadings,object@dr$ica@cell.embeddings)
  
  n_ICs = dim(x[[1]])[2]
  # Calculate skews of cell scores and gene loadings
  gene_skews_initial = apply( x[[1]], 2, skewness )
  cell_skews_initial = apply( x[[2]], 2, skewness )
  
  # Flip certain ICs so that all ICs have positive skew for cell scores
  wh_ICs_flip = which( cell_skews_initial < 0 )
  x[[1]][,wh_ICs_flip] = -1 * x[[1]][,wh_ICs_flip]
  x[[2]][,wh_ICs_flip] = -1 * x[[2]][,wh_ICs_flip]
  
  # Reorder the ICs by their skew (from best to worst)
  IC_reorder = order( abs(cell_skews_initial), decreasing=T )
  x[[1]] = x[[1]][,IC_reorder]
  x[[2]] = x[[2]][,IC_reorder]
  
  # Calculate cell-kurtosis(k), gene-skew(gs), and cell-skew(cs) for the reordered ICs
  gk = apply( x[[1]], 2, kurtosis )
  ck = apply( x[[2]], 2, kurtosis )
  gs = apply( x[[1]], 2, skewness )
  cs = apply( x[[2]], 2, skewness )
  if (!is.null(curation.data)) {
    curation.names = curation.data$status
    
    keep_IC=(curation.names == "Real")
    
    idx.keep=which(keep_IC)
    
  } else {
    keep_IC = rep(TRUE,times = ncol(x[[2]]))
    idx.keep = which(keep_IC)
    cluster.doublets = F
    curation.names = rep("Uncurated",times = ncol(x[[2]]))
  }
  
  
  
  # Routine for 1-D clustering to remove Doublet and Outlier cells.
  if(cluster.doublets) {
    
    
    idx.doublet=which(grepl(curation.names,pattern="Doublet",ignore.case=T) | grepl(curation.names,pattern="Outlier",ignore.case=T))
    
    #Obtain the mode for centering the 1-D gaussian using the density() function, and the turnpoints function from pastecs
    
    doublets=lapply(idx.doublet,function(i){
      q=x[[2]][,i]
      d = density(q)
      ts_y = ts(d$y)
      tp = turnpoints(ts_y)
      idx = which.max(d$y[tp$tppos])
      center=d$x[tp$tppos[idx]]
      pvalsUpper=pnorm(q,mean=center,sd= sd(q),lower.tail=F)
      pvalsUpper=p.adjust(pvalsUpper,"fdr")
      idx=which(pvalsUpper <= doubletalpha)
      return(idx)
    })
    
    
  } else {
    idx.doublet=integer(0)
  }
  
  #Store reoriented ICs
  colnames(x[[2]])<-paste0("IC",1:ncol(x[[2]]))
  colnames(x[[1]])<-paste0("IC",1:ncol(x[[2]]))
  
  
  object@dr$ica@cell.embeddings = x[[2]]
  object@dr$ica@gene.loadings = x[[1]]
  object@dr$ica@misc = curation.names
  
  #If plotting on new tSNE, make new embedding with only retained ICs.
  if (make.tsne)
    object = RunTSNE(object = object,reduction.use = 'ica',dims.use = which(keep_IC),do.fast=T,pca=F,check_duplicates=F)
  
  
  # Skew-kurtosis plot of the ICs
  par(mfrow=c(1,1))
  plot(  log(ck,2), cs , col="white",
         xlab = "log2( Cell-score kurtosis)", ylab="Cell-score skew",
         main = "Cell-score skew and kurtosis (after component reordering)") 
  rect( 0, 0, log2(kurtosis_threshold), skewness_threshold, col="lightgray", border=NA )
  
  
  text( log(ck, 2), cs, 1:length(cs)  , col="black" )
  
  abline(v=0); abline(h=0)
  
  #  Review the cell-score and gene-loading distributions for each component
  
  for(i in 1:n_ICs) {
    par(mfrow=c(2,1))
    print(i)
    top_genes = row.names( x[[1]] )[ order(x[[1]][,i], decreasing=T )[1:7] ]
    top_genes_textstring = paste(top_genes, collapse=", ")
    stats_textstring1 = paste( "k=", round(ck[i],1), ", s=", round(cs[i],1))
    stats_textstring2 = paste( "k=", round(gk[i],1), ", s=", round(gs[i],1))
    component_textstring = paste( "Component", i, "  (", curation.names[i],")" )
    if( !keep_IC[i] ) component_textstring=paste(component_textstring, "(REMOVE)" )
    plot_title1 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring1 )
    plot_title2 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring2 )
    
    # Plot cell scores
    plot( x[[2]][,i], cex=0.2, main=plot_title1, xlab="Cell", ylab="Score"  )
    
    # Highlight the doublet-flagged cells.
    if (i %in% idx.doublet) {
      idx.use = which(idx.doublet==i)
      
      points(doublets[[idx.use]],x[[2]][doublets[[idx.use]],i],cex=0.5,col='red')
    }
    
    # Plot gene loadings
    plot( x[[1]][,i], cex=0.2, main=plot_title2, xlab="Gene", ylab="Loading"  )
    if (feature.plot) { 
      tsrots = object@dr$tsne@cell.embeddings
      par(mfrow=c(1,1))
      fplot(tsrots[rownames(x[[2]]),],x[[2]][,i],title=component_textstring)
    }
  }
  
  # Skew-skew plot of the ICs
  par(mfrow=c(1,1))
  plot(  gs, cs , col="white",
         xlab = "Gene-loading skew", ylab="Cell-score skew",
         main = "Cell-score and gene-loading skew (after component reordering)") 
  text( gs, cs,1:length(cs)   , col="black" )
  abline(v=0); abline(h=0)
  
  
  if(cluster.doublets) {
    object = SubsetData(object = object,cells.use = setdiff(rownames(x[[2]]),rownames(x[[2]])[unlist(doublets)]),do.scale = F)
    object@raw.data=object@raw.data[,colnames(object@data)]
    object@raw.data=object@raw.data[which(Matrix::rowSums(object@raw.data)>0),]
  }
  return(object)
}

fplot<-function(tsne,factor,title,cols.use=heat.colors(10),pt.size=0.7,pch.use=20) {
  c=intersect(rownames(tsne),names(factor))
  data.use=factor[c]
  data.cut=as.numeric(as.factor(cut(as.numeric(data.use),breaks=length(cols.use))))
  data.col=rev(cols.use)[data.cut]
  plot(tsne[c,1],tsne[c,2],col=data.col,cex=pt.size,pch=pch.use,main=title)
  
}

environment(CurateICA.seurat)<-asNamespace('Seurat')

StoreSpatialData<-function(object,spatial_coord,bead.size =10) {
  cells.use = intersect(names(object@ident),rownames(spatial_coord))
  spatial.obj = new(Class="dim.reduction",gene.loadings=matrix(0),cell.embeddings=as.matrix(spatial_coord[cells.use,]),key="SPATIAL_")
  colnames(spatial.obj@cell.embeddings) = paste0("SPATIAL_",1:ncol(spatial.obj@cell.embeddings))
  object = SubsetData(object,cells.use = cells.use)
  object@raw.data = object@raw.data[,colnames(object@data)]
  object@dr$spatial=spatial.obj
  return(object)
}
environment(StoreSpatialData)<-asNamespace('Seurat')


AtlasToSeurat<-function(data_folder_path,subcluster=NULL,keep.removed = F) {
  data_folder_path=sub(x=data_folder_path,pattern = "/$",replacement = "")
  
  split=strsplit(data_folder_path,split="/")[[1]]
  samplename=split[length(split)]
  idents=readRDS(paste0(data_folder_path,"/assign/",samplename,".subcluster.assign.RDS"))
  idents.level1 = readRDS(paste0(data_folder_path,"/assign/",samplename,".cluster.assign.RDS"))
  if (is.null(subcluster)) {
    rawdata=readRDS(paste0(data_folder_path,"/dge/",samplename,".filtered.raw.dge.RDS"))
    if (class(rawdata)[1] == "dgTMatrix") {
      temp = summary(rawdata)
      col.names = colnames(rawdata)
      row.names = rownames(rawdata)
      rawdata = sparseMatrix(i=temp[,1],j=temp[,2],x=temp[,3])
      rownames(rawdata) = row.names
      colnames(rawdata) = col.names
    }
    temp.factor = factor(sapply(as.character(idents),function(i){strsplit(i,split="-")[[1]][1]}))
    names(temp.factor) = names(idents)
    cells.keep = unlist(lapply(levels(idents.level1),function(i){
      if (i %in% levels(temp.factor)) {
        x = names(temp.factor)[which(temp.factor == i)]
        y = names(idents.level1)[which(idents.level1==i)]
        return(intersect(x,y))
      } else {
        return(names(idents.level1)[which(idents.level1==i)])
      }
      
    }))
    new.idents = rep("NA",times=length(cells.keep))
    names(new.idents) = cells.keep
    cells.sub = intersect(cells.keep,names(idents))
    new.idents[cells.sub] = as.character(idents[cells.sub])
    cells.top = setdiff(cells.keep,cells.sub)
    new.idents[cells.top] = paste0(idents.level1[cells.top],"-1")
    idents.use = factor(new.idents)
    names(idents.use) = names(new.idents)
    normalized.dge=readRDS(paste0(data_folder_path,"/dge/",samplename,".filtered.scaled.dge.RDS"))
    orig.dge=rawdata[,names(idents.use)]
    tsne=readRDS(paste0(data_folder_path,"/tSNE/",samplename,"_tSNExy.RDS"))[names(idents.use),]
    ics=readRDS(paste0(data_folder_path,"/components/",samplename,".ica.RDS"))
    my.seurat=CreateSeuratObject(raw.data=orig.dge)
    my.seurat@data=normalized.dge[,names(idents.use)]
    
    my.seurat@var.genes = rownames(ics$gene_loadings)
    my.seurat@scale.data=t(scale(t(normalized.dge[my.seurat@var.genes,]),center=T,scale=T))
    ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations[names(idents.use),],key="IC")
    colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))    
    colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
    tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
    colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
    my.seurat@dr$ica=ica.obj
    my.seurat@dr$tsne=tsne.obj
    my.seurat@ident=idents.use
    class.path = list.files(path=paste0(data_folder_path,"/curation_sheets/"),pattern="*cluster_class.csv",full.names=T)
    class.data = read.csv(class.path)
    return(my.seurat)
    
    
    
  } else {
    rawdata=readRDS(paste0(data_folder_path,"/cluster",subcluster,"/",samplename,".subcluster_inputs.RDS"))
    orig.dge=rawdata$raw_dge
    scaled.dge=rawdata$scaled_gene_selected_dge
    
    tsne=readRDS(paste0(data_folder_path,"/tSNE/",samplename,".cluster",subcluster,".CURATEDtSNE.RDS"))
    if(is.null(tsne)) {
      tsne=readRDS(paste0(data_folder_path,"/cluster",subcluster,"/",samplename,".cluster",subcluster,".auto.tSNExy.RDS"))
      
    }
    icfiles=list.files(paste0(data_folder_path,"/components/"))
    icfile.use=icfiles[grep(icfiles,pattern=paste0("cluster",subcluster,"\\."))]
    ics=readRDS(paste0(data_folder_path,"/components/",icfile.use))
    curationsheets=list.files(paste0(data_folder_path,"/curation_sheets/"))
    curationfile=curationsheets[grep(curationsheets,pattern=paste0("Subcluster_",subcluster,"[_\\.]"))]
    ic.annotations = read.table(paste0(data_folder_path,"/curation_sheets/",curationfile), stringsAsFactors=F, sep=",", header=T, quote="\"")
    my.seurat=CreateSeuratObject(raw.data=orig.dge)
    
    my.seurat@data = sweep(my.seurat@raw.data,2,Matrix::colSums(my.seurat@raw.data),"/")
    my.seurat@scale.data=scaled.dge
    my.seurat@var.genes = rownames(ics$gene_loadings)
    ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations,key="IC",misc=ic.annotations)
    colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))
    colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
    tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
    colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
    cells.int=intersect(rownames(tsne),names(idents))
    cells.removed=setdiff(rownames(tsne),cells.int)
    idents.new=factor(x = c(as.character(idents[cells.int]),rep("REMOVED",times=length(cells.removed))))
    names(idents.new) = c(cells.int,cells.removed)
    idents.new = idents.new[rownames(tsne)]
    my.seurat@dr$ica=ica.obj
    my.seurat@dr$tsne=tsne.obj
    my.seurat@ident=idents.new
    return(my.seurat)
    
  }
  
}
environment(AtlasToSeurat)<-asNamespace('Seurat')

AtlasToSeuratWindows<-function(data_folder_path,subcluster=NULL,keep.removed = F) {
  data_folder_path=sub(x=data_folder_path,pattern = "\\$",replacement = "")
  
  split=strsplit(data_folder_path,split="\\\\")[[1]]
  samplename=split[length(split)]
  idents=readRDS(paste0(data_folder_path,"\\assign\\",samplename,".subcluster.assign.RDS"))
  idents.level1 = readRDS(paste0(data_folder_path,"\\assign\\",samplename,".cluster.assign.RDS"))
  if (is.null(subcluster)) {
    rawdata=readRDS(paste0(data_folder_path,"\\dge\\",samplename,".filtered.raw.dge.RDS"))
    if (class(rawdata)[1] == "dgTMatrix") {
      temp = summary(rawdata)
      col.names = colnames(rawdata)
      row.names = rownames(rawdata)
      rawdata = sparseMatrix(i=temp[,1],j=temp[,2],x=temp[,3])
      rownames(rawdata) = row.names
      colnames(rawdata) = col.names
    }
    temp.factor = factor(sapply(as.character(idents),function(i){strsplit(i,split="-")[[1]][1]}))
    names(temp.factor) = names(idents)
    cells.keep = unlist(lapply(levels(idents.level1),function(i){
      if (i %in% levels(temp.factor)) {
        x = names(temp.factor)[which(temp.factor == i)]
        y = names(idents.level1)[which(idents.level1==i)]
        return(intersect(x,y))
      } else {
        return(names(idents.level1)[which(idents.level1==i)])
      }
      
    }))
    new.idents = rep("NA",times=length(cells.keep))
    names(new.idents) = cells.keep
    cells.sub = intersect(cells.keep,names(idents))
    new.idents[cells.sub] = as.character(idents[cells.sub])
    cells.top = setdiff(cells.keep,cells.sub)
    new.idents[cells.top] = paste0(idents.level1[cells.top],"-1")
    idents.use = factor(new.idents)
    names(idents.use) = names(new.idents)
    normalized.dge=readRDS(paste0(data_folder_path,"\\dge\\",samplename,".filtered.scaled.dge.RDS"))
    orig.dge=rawdata[,names(idents.use)]
    tsne=readRDS(paste0(data_folder_path,"\\tSNE\\",samplename,"_tSNExy.RDS"))[names(idents.use),]
    ics=readRDS(paste0(data_folder_path,"\\components\\",samplename,".ica.RDS"))
    my.seurat=CreateSeuratObject(raw.data=orig.dge)
    my.seurat@data=normalized.dge[,names(idents.use)]
    
    my.seurat@var.genes = rownames(ics$gene_loadings)
    my.seurat@scale.data=t(scale(t(normalized.dge[my.seurat@var.genes,]),center=T,scale=T))
    ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations[names(idents.use),],key="IC")
    colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))    
    colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
    tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
    colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
    my.seurat@dr$ica=ica.obj
    my.seurat@dr$tsne=tsne.obj
    my.seurat@ident=idents.use
    class.path = list.files(path=paste0(data_folder_path,"\\curation_sheets\\"),pattern="*cluster_class.csv",full.names=T)
    class.data = read.csv(class.path)
    return(my.seurat)
    
    
    
  } else {
    rawdata=readRDS(paste0(data_folder_path,"\\cluster",subcluster,"\\",samplename,".subcluster_inputs.RDS"))
    orig.dge=rawdata$raw_dge
    scaled.dge=rawdata$scaled_gene_selected_dge
    
    tsne=readRDS(paste0(data_folder_path,"\\tSNE\\",samplename,".cluster",subcluster,".CURATEDtSNE.RDS"))
    if(is.null(tsne)) {
      tsne=readRDS(paste0(data_folder_path,"\\cluster",subcluster,"\\",samplename,".cluster",subcluster,".auto.tSNExy.RDS"))
      
    }
    icfiles=list.files(paste0(data_folder_path,"\\components\\"))
    icfile.use=icfiles[grep(icfiles,pattern=paste0("cluster",subcluster,"\\."))]
    ics=readRDS(paste0(data_folder_path,"\\components\\",icfile.use))
    curationsheets=list.files(paste0(data_folder_path,"\\curation_sheets\\"))
    curationfile=curationsheets[grep(curationsheets,pattern=paste0("Subcluster_",subcluster,"[_\\.]"))]
    ic.annotations = read.table(paste0(data_folder_path,"\\curation_sheets\\",curationfile), stringsAsFactors=F, sep=",", header=T, quote="\"")
    my.seurat=CreateSeuratObject(raw.data=orig.dge)
    
    my.seurat@data = sweep(my.seurat@raw.data,2,Matrix::colSums(my.seurat@raw.data),"\\")
    my.seurat@scale.data=scaled.dge
    my.seurat@var.genes = rownames(ics$gene_loadings)
    ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations,key="IC",misc=ic.annotations)
    colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))
    colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
    tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
    colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
    cells.int=intersect(rownames(tsne),names(idents))
    cells.removed=setdiff(rownames(tsne),cells.int)
    idents.new=factor(x = c(as.character(idents[cells.int]),rep("REMOVED",times=length(cells.removed))))
    names(idents.new) = c(cells.int,cells.removed)
    idents.new = idents.new[rownames(tsne)]
    my.seurat@dr$ica=ica.obj
    my.seurat@dr$tsne=tsne.obj
    my.seurat@ident=idents.new
    return(my.seurat)
    
  }
  
}
environment(AtlasToSeuratWindows)<-asNamespace('Seurat')





RunNMF<-function(object, factors.compute = 50, log.norm = FALSE, print.results = TRUE, 
                 factors.print = 1:factors.compute, genes.print = 50, seed.use = 1, ...) {
  if (log.norm) {
    norm=log(Matrix.column_norm(object@raw.data)*10000+1)
  } else {
    norm=Matrix.column_norm(object@raw.data)
  }
  uncentered.scaled = scale(t(norm[object@var.genes,]),center=F,scale=T)
  object@scale.data = uncentered.scaled
  print(dim(uncentered.scaled))
  # uncentered.scaled = apply(norm[object@var.genes,],1,scalar1)
  factors.compute=min(c(factors.compute,dim(uncentered.scaled)))
  set.seed(seed=seed.use)
  nmf.results = nnmf(uncentered.scaled,k=factors.compute)
  nmf.obj<-new(Class="dim.reduction",gene.loadings=t(nmf.results$H),cell.embeddings=nmf.results$W,key="NMF",misc=list("raw" = nmf.results))
  colnames(nmf.obj@cell.embeddings) = paste0("Factor",1:ncol(nmf.obj@cell.embeddings))
  colnames(nmf.obj@gene.loadings) = paste0("Factor",1:ncol(nmf.obj@gene.loadings))
  object@dr$nmf <- nmf.obj
  if(print.results) {
    for (i in factors.print) {
      code <- paste0(GetDimReduction(object = object, reduction.type = "nmf",
                                     slot = "key"), i)
      sx <- DimTopGenes(object = object, dim.use = i, reduction.type = 'nmf',
                        num.genes = genes.print*2, use.full = F,do.balanced = FALSE)
      print(code)
      #print((sx[1:genes.print]))
      #print("")
      print(rev(x = (sx[(length(x = sx) - genes.print + 1):length(x = sx)])))
      print("")
      print("")
    }
    
    
  }
  return(object)
}

environment(RunNMF)<-asNamespace('Seurat')



CurateNMF.seurat<-function(object,make.tsne = T, feature.plot=TRUE,reduction.embedding = 'tsne',kurtosis_threshold = 6,skewness_threshold = 1,curation.data=NULL,cluster.doublets = T,doubletalpha=0.01,pt.size = 0.7,dark.theme = F,cols.use = c("yellow","red"),reduction.use = 'nmf',do.reorder = T,counterstain = NULL) {
  x = list(object@dr[[reduction.use]]@gene.loadings,object@dr[[reduction.use]]@cell.embeddings)
  
  n_ICs = dim(x[[1]])[2]
  # Calculate skews of cell scores and gene loadings
  gene_skews_initial = apply( x[[1]], 2, skewness )
  cell_skews_initial = apply( x[[2]], 2, skewness )
  
  # Flip certain ICs so that all ICs have positive skew for cell scores
  wh_ICs_flip = which( cell_skews_initial < 0 )
  x[[1]][,wh_ICs_flip] = -1 * x[[1]][,wh_ICs_flip]
  x[[2]][,wh_ICs_flip] = -1 * x[[2]][,wh_ICs_flip]
  
  # Reorder the ICs by their skew (from best to worst)
  if (do.reorder) {
    IC_reorder = order( abs(cell_skews_initial), decreasing=T )
    x[[1]] = x[[1]][,IC_reorder]
    x[[2]] = x[[2]][,IC_reorder]
  } else {
    IC_reorder = 1:ncol(x[[1]])
  }
  # Calculate cell-kurtosis(k), gene-skew(gs), and cell-skew(cs) for the reordered ICs
  gk = apply( x[[1]], 2, kurtosis )
  ck = apply( x[[2]], 2, kurtosis )
  gs = apply( x[[1]], 2, skewness )
  cs = apply( x[[2]], 2, skewness )
  if (!is.null(curation.data)) {
    curation.names = curation.data$status
    
    keep_IC=(curation.names == "Real")
    
    idx.keep=which(keep_IC)
    
  } else {
    keep_IC = rep(TRUE,times = ncol(x[[2]]))
    idx.keep = which(keep_IC)
    cluster.doublets = F
    curation.names = rep("Uncurated",times = ncol(x[[2]]))
  }
  
  
  
  # Routine for 1-D clustering to remove Doublet and Outlier cells.
  if(cluster.doublets) {
    
    
    idx.doublet=which(grepl(curation.names,pattern="Doublet",ignore.case=T) | grepl(curation.names,pattern="Outlier",ignore.case=T))
    
    #Obtain the mode for centering the 1-D gaussian using the density() function, and the turnpoints function from pastecs
    
    doublets=lapply(idx.doublet,function(i){
      q=x[[2]][,i]
      d = density(q)
      ts_y = ts(d$y)
      tp = turnpoints(ts_y)
      idx = which.max(d$y[tp$tppos])
      center=d$x[tp$tppos[idx]]
      pvalsUpper=pnorm(q,mean=center,sd= sd(q),lower.tail=F)
      pvalsUpper=p.adjust(pvalsUpper,"fdr")
      idx=which(pvalsUpper <= doubletalpha)
      return(idx)
    })
    
    
  } else {
    idx.doublet=integer(0)
  }
  
  #Store reoriented ICs
  colnames(x[[2]])<-paste0("Factor",1:ncol(x[[2]]))
  colnames(x[[1]])<-paste0("Factor",1:ncol(x[[2]]))
  object@dr[[reduction.use]]@cell.embeddings = x[[2]]
  object@dr[[reduction.use]]@gene.loadings = x[[1]]
  object@dr[[reduction.use]]@misc$curation = curation.names
  object@dr[[reduction.use]]@misc$raw$W = object@dr[[reduction.use]]@misc$raw$W[,IC_reorder]
  object@dr[[reduction.use]]@misc$raw$H = object@dr[[reduction.use]]@misc$raw$H[IC_reorder,]
  #If plotting on new tSNE, make new embedding with only retained ICs.
  if (make.tsne)
    object = RunTSNE(object = object,reduction.use = reduction.use,dims.use = which(keep_IC))
  
  
  # Skew-kurtosis plot of the ICs
  par(mfrow=c(1,1))
  plot(  log(ck,2), cs , col="white",
         xlab = "log2( Cell-score kurtosis)", ylab="Cell-score skew",
         main = "Cell-score skew and kurtosis (after component reordering)") 
  rect( 0, 0, log2(kurtosis_threshold), skewness_threshold, col="lightgray", border=NA )
  
  
  text( log(ck, 2), cs, 1:length(cs)  , col="black" )
  
  abline(v=0); abline(h=0)
  
  #  Review the cell-score and gene-loading distributions for each component
  
  for(i in 1:n_ICs) {
    par(mfrow=c(2,1))
    print(i)
    top_genes = row.names( x[[1]] )[ order(x[[1]][,i], decreasing=T )[1:7] ]
    top_genes_textstring = paste(top_genes, collapse=", ")
    stats_textstring1 = paste( "k=", round(ck[i],1), ", s=", round(cs[i],1))
    stats_textstring2 = paste( "k=", round(gk[i],1), ", s=", round(gs[i],1))
    component_textstring = paste( "Component", i, "  (", curation.names[i],")" )
    if( !keep_IC[i] ) component_textstring=paste(component_textstring, "(REMOVE)" )
    plot_title1 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring1 )
    plot_title2 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring2 )
    
    # Plot cell scores
    plot( x[[2]][,i], cex=0.2, main=plot_title1, xlab="Cell", ylab="Score"  )
    
    # Highlight the doublet-flagged cells.
    if (i %in% idx.doublet) {
      idx.use = which(idx.doublet==i)
      
      points(doublets[[idx.use]],x[[2]][doublets[[idx.use]],i],cex=0.5,col='red')
    }
    
    # Plot gene loadings
    plot( x[[1]][,i], cex=0.2, main=plot_title2, xlab="Gene", ylab="Loading"  )
    if (feature.plot) { 
      par(mfrow=c(1,1))
      FeaturePlot(object,pt.size = pt.size,reduction.use = reduction.embedding,features.plot = paste0(object@dr[[reduction.use]]@key,i),dark.theme = dark.theme,cols.use = cols.use)
      
    }
  }
  
  # Skew-skew plot of the ICs
  par(mfrow=c(1,1))
  plot(  gs, cs , col="white",
         xlab = "Gene-loading skew", ylab="Cell-score skew",
         main = "Cell-score and gene-loading skew (after component reordering)") 
  text( gs, cs,1:length(cs)   , col="black" )
  abline(v=0); abline(h=0)
  
  
  if(cluster.doublets) {
    object = SubsetData(object = object,cells.use = setdiff(rownames(x[[2]]),unlist(doublets)),do.scale = F)
  }
  return(object)
}
environment(CurateNMF.seurat)<-asNamespace('Seurat')


RunNMFreg<-function(object,refmat,umi.weighted = F) {
  genes.use = intersect(rownames(object@data),colnames(refmat$H))
  sc = scale(t(object@data[genes.use,]),center=F,scale=T)
  refmat$H = refmat$H[,genes.use]
  reg = predict(refmat,newdata=sc,which="W")
  if (umi.weighted) {
    reg.obj  <- new(Class = "dim.reduction",gene.loadings = t(refmat$H),cell.embeddings = sweep(reg$coefficients,1,colSums(object@raw.data)[colnames(object@data)],"*"),key="NMFreg",misc = list("raw" = reg))
  }else {
    reg.obj  <- new(Class = "dim.reduction",gene.loadings = t(refmat$H),cell.embeddings = reg$coefficients,key="NMFreg",misc = list("raw" = reg)) 
  }
  rownames(reg.obj@gene.loadings) = genes.use
  object@dr$nmfreg <- reg.obj
  return(object)
  
}


fplot<-function(tsne,factor,title,cols.use=heat.colors(10),pt.size=0.7,pch.use=20) {
  c=intersect(rownames(tsne),names(factor))
  data.use=factor[c]
  data.cut=as.numeric(as.factor(cut(as.numeric(data.use),breaks=length(cols.use))))
  data.col=rev(cols.use)[data.cut]
  plot(tsne[c,1],tsne[c,2],col=data.col,cex=pt.size,pch=pch.use,main=title)
  
}

correlate.factors<-function(object1,object2) {
  H1=object1@dr$nmf@gene.loadings
  H1=H1[,which(object1@dr$nmf@misc$curation=="Real")]
  H2=object2@dr$nmf@gene.loadings
  H2=H2[,which(object2@dr$nmf@misc$curation=="Real")]
  shared.genes = intersect(rownames(H1),rownames(H2))
  cors=cor(H1[shared.genes,],H2[shared.genes,])
  if(ncol(H1)==ncol(H2))
    cors=matrix.sort(cors)
  return(cors)
  
}

correlate.ICs<-function(object1,object2) {
  H1=object1@dr$ica@gene.loadings
  H2=object2@dr$ica@gene.loadings
  shared.genes = intersect(rownames(H1),rownames(H2))
  cors=cor(H1[shared.genes,],H2[shared.genes,])
  if(ncol(H1)==ncol(H2))
    cors=matrix.sort(cors)
  return(cors)
  
}



matrix.sort <- function(matrix, require_square=TRUE) {
  
  if (require_square && nrow(matrix) != ncol(matrix)) stop("Not diagonal")
  if(is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
  
  row.max <- apply(matrix,1,which.max)
  if(all(table(row.max) != 1)) stop("Ties cannot be resolved")
  
  matrix[names(sort(row.max)),]
}


BuildRFClassifierSCALED<- function (object, training.genes = NULL, training.classes = NULL, 
                                    verbose = TRUE, ...) 
{
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(x = training.genes, default = rownames(x = object@data))
  training.data <- as.data.frame(x = as.matrix(x = scale(t(x = object@data[training.genes, 
                                                                           ])),center=T,scale=T))
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    print("Training Classifier ...")
  }
  classifier <- ranger(data = training.data, dependent.variable.name = "class", 
                       classification = TRUE, write.forest = TRUE, ...)
  return(classifier)
}
environment(BuildRFClassifierSCALED)<-asNamespace('Seurat')


downsample.matrix <- function(mat, samplerate=0.8,seed=1) {
  set.seed(seed)
  new <- matrix(0, nrow(mat), ncol(mat))
  colnames(new) <- colnames(mat)
  rownames(new) <- rownames(mat)
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      new[i,j] <- sum(runif(mat[i,j], 0, 1) < samplerate)
    }
  }
  return(new)
}




#Function takes in a list of DGEs, with gene rownames and cell colnames, and merges them into a single DGE.

generate.merged.dge = function(base.path,libraries,library.names,min.umis = 800) {
  data.list = lapply(1:length(libraries),function(j) {
    dat = Convert10xtoMatrix(load_cellranger_matrix(paste0(base.path,"/",libraries[j]),barcode_filtered = F))
    umis.per.cell = colSums(dat)
    idx = which (umis.per.cell > min.umis)
    dat = dat[,idx]
    
    dat = dat[which(rowSums(dat)>0),]
    return(dat)
  })
  MergeSparseDataAll(data.list,library.names = library.names)
  
}

Matrix.column_norm <- function(A){
  if (class(A)[1] == "dgTMatrix") {
    temp = summary(A)
    col.names = colnames(A)
    row.names = rownames(A)
    A = sparseMatrix(i=temp[,1],j=temp[,2],x=temp[,3])
    rownames(A) = row.names
    colnames(A) = col.names
  }
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  return(A)
}


Matrix.tapply <- function(X, INDEX, FUN=NULL, ..., simplify=TRUE) {
  ## Matrix tapply
  ## X: matrix with n rows; INDEX: vector or list of vectors of length n
  ## FUN: function to operate on submatrices of x by INDEX
  ## ...: arguments to FUN; simplify: see sapply
  
  idx.list <- tapply(seq(ncol(X)), INDEX, c)
  sapply(idx.list, function(idx,x,fun,...) fun(x[,idx,drop=FALSE],...),
         x=X, fun=FUN, ..., simplify=simplify)
}




suggest.centers <- function(ICAcell.matrix, peakthresh = 0.15, peakpos = 0.9, centeralpha=0.1, min.cells.ic = 5, idx.keep,manual=F,threshes) {
  cell.list=lapply(1:ncol(ICAcell.matrix),function(q){
    ic = ICAcell.matrix[,q]
    if (!manual) {
      d = density(ic)
      
      ts_y = ts(d$y)
      tp = turnpoints(ts_y)
      
      if (tp$firstispeak) {
        peak.pos = tp$tppos[((1:length(tp$tppos) %% 2) ==1)]
        #pit.pos =  tp$tppos[((1:length(tp$tppos) %% 2) ==0)]
      } else {
        peak.pos = tp$tppos[((1:length(tp$tppos) %% 2) ==0)]
        #pit.pos =  tp$tppos[((1:length(tp$tppos) %% 2) ==1)]
      }
      
      idx = min(which(d$y[peak.pos] > peakthresh))
      center.use =d$x[peak.pos[idx]]
      
      #For each center, use the left trough of each peak to identify a reasonable range of the distribution from which to calculate the variance.
      low.idx = max(which(d$y[tp$tppos] < (peakpos*d$y[peak.pos[idx]]) & d$x[tp$tppos] < center.use))
      if(!is.finite(low.idx)) {
        low.idx = min(ic)
      } else {
        low.idx=d$x[tp$tppos[low.idx]]
      }
      
      #abline(v = low.idx,lty = 2,col ='red')
      ic.use = ic[which(ic < center.use & ic > low.idx)]
      ic.reflected = 2*center.use - ic.use
      #  parms=JohnsonFit(t=c(center.use,sd(c(ic.use,ic.reflected)),skewness(c(ic.use,ic.reflected)),kurtosis(c(ic.use,ic.reflected))),moment = 'use')
      pvalsUpper=pnorm(ic,mean=center.use,sd= sd(c(ic.use,ic.reflected)),lower.tail=F)
      # pvalsUpper=pJohnson(ic,parms,lower.tail = F)
      pvalsUpper=p.adjust(pvalsUpper,"fdr")
      
      
      idx=which(pvalsUpper <= centeralpha)
    } else {
      idx=which(ic > threshes[q])
    }
    if (length(idx) > min.cells.ic) {
      return((rownames(ICAcell.matrix)[idx]))
    } else {
      return(NULL)
    }
  })
  
  # In many instances, there will be a large population of "null" or "default" cells that do not have their own positive IC loading.   In these instances, add an extra cluster center to accommodate them.
  
  return(cell.list)
}


skewsbystate = function(ics,ic.names,threshes,state1="FC",state2="PC",iterations=10,alpha = 0.1,cols=c("#F90606","#7A0303"),johnson.mode = 'quant',manual.threshes,do.manual=F,min.cells=50) {
  
  rowMeans(sapply(1:iterations,function(x){
    cells.1 = grep(rownames(ics),pattern=state1)
    cells.2 = grep(rownames(ics),pattern=state2)
    s1 = length(cells.1)
    s2 = length(cells.2)
    ds = min(c(s1,s2))
    set.seed(seed=x)
    cells.1 = sample(cells.1,size=ds)
    cells.2=sample(cells.2,size= ds)
    ics = rbind(ics[cells.1,],ics[cells.2,])
    
    pos.cells = suggest.centers(ics,idx.keep = ic.names,centeralpha = alpha,manual = do.manual,threshes=manual.threshes)
    
    #  dev.off()
    skews=unlist(lapply(1:length(pos.cells),function(i){
      q = pos.cells[[i]]
      if (length(q) > min.cells) {
        ic = ics[,i]
        
        names(ic)<-rownames(ics)
        idx.state1 = q[grep(q,pattern=state1)]
        idx.state2 = q[grep(q,pattern=state2)]
        #skew.sum = (sum(ic[idx.state1])/(sum(ic[idx.state1])+sum(ic[idx.state2])))
        #return(skew.sum)
        return(length(idx.state1)/(length(idx.state1)+length(idx.state2)))
      } else {
        return(NA)
      }
    }))  
    
    if (x ==1) {
      for (i in 1:ncol(ics)) {
        plot(1:nrow(ics),ics[,i],pch=20,cex=0.5,col=gray(0.5))
        tissue.id = vector(length = length(pos.cells[[i]]))
        idx.1 = grep(pos.cells[[i]],pattern=state1)
        idx.2 = grep(pos.cells[[i]],pattern=state2)
        tissue.id[idx.1] = cols[1]
        tissue.id[idx.2] = cols[2]
        points(which(rownames(ics) %in% pos.cells[[i]]),ics[pos.cells[[i]],i],pch=20,cex=0.5,col=tissue.id)
        legend('topright',legend = c(state1,state2),col = c(cols[1],cols[2]),pch=20)
        abline(h=min(ics[pos.cells[[i]],i]),lty=2,col='black')
        title(main=paste0("IC ",ic.names[i]," skew = ",round(skews[i],3)))
      }
      
    }
    return(skews)
  }))
  
  
}



align_quantiles = function(Ps,clusters,quantiles=50)
{
  dims = ncol(Ps[[1]])
  num_clusters = dims
  for (k in 1:length(Ps))
  {
    for (i in 1:dims)
    {
      for (j in 1:num_clusters)
      {
        if (sum(clusters[[1]]==j)==0 | sum(clusters[[k]]==j)==0){next}
        if (sum(clusters[[k]]==j)==1){
          Ps[[k]][clusters[[k]]==j,i] = mean(Ps[[1]][clusters[[1]]==j,i])
          next
        }
        q2 = quantile(Ps[[k]][clusters[[k]]==j,i],seq(0,1,by=1/quantiles))
        q1 = quantile(Ps[[1]][clusters[[1]]==j,i],seq(0,1,by=1/quantiles))
        if (sum(q1)==0 | sum(q2)==0)
        {
          new_vals = rep(0,sum(clusters[[k]]==j))
        }
        else
        {
          warp_func = approxfun(q2,q1)
          new_vals = warp_func(Ps[[k]][clusters[[k]]==j,i])
        }
        
        Ps[[k]][clusters[[k]]==j,i] = new_vals
      }
    }
  }
  return(Reduce(rbind,Ps))
}

store.iNMF<-function(object,iNMF.H, iNMF.W,iNMF.V) {
  object@dr$iNMF=new(Class="dim.reduction",gene.loadings = as.matrix(iNMF.W), misc= iNMF.V,cell.embeddings=as.matrix(iNMF.H),key="iNMF")
  return(object)
}
environment(store.iNMF) <- asNamespace('Seurat')



selectGenes_sparse = function (object, alphathresh = 0.99, varthresh = 0.1, cex.use = 0.3,
                               combine = "union", keep.unique = F, capitalize = F)
{
  genes.use = c()
  for (i in 1:length(object@raw.data)) {
    if (capitalize) {
      rownames(object@raw.data[[i]]) = toupper(rownames(object@raw.data[[i]]))
      rownames(object@norm.data[[i]]) = toupper(rownames(object@norm.data[[i]]))
    }
    trx_per_cell <- Matrix::colSums(object@raw.data[[i]])
    gene_expr_mean <- Matrix::rowMeans(object@norm.data[[i]])
    gene_expr_var <- sparse_var(object@norm.data[[i]])
    names(gene_expr_var) = rownames(object@norm.data[[i]])
    nolan_constant <- mean((1/trx_per_cell))
    alphathresh.corrected = alphathresh/dim(object@raw.data[[i]])[1]
    genemeanupper <- gene_expr_mean + qnorm(1 - alphathresh.corrected/2) *
      sqrt(gene_expr_mean * nolan_constant/dim(object@raw.data[[i]])[2])
    genes.new = names(gene_expr_var)[which(gene_expr_var/nolan_constant >
                                             genemeanupper & log10(gene_expr_var) > log10(gene_expr_mean) +
                                             (log10(nolan_constant) + varthresh))]
    plot(log10(gene_expr_mean), log10(gene_expr_var), cex = cex.use)
    points(log10(gene_expr_mean[genes.new]), log10(gene_expr_var[genes.new]),
           cex = cex.use, col = "green")
    abline(log10(nolan_constant), 1, col = "purple")
    legend("bottomright", paste0("Selected genes: ", length(genes.new)),
           pch = 20, col = "green")
    if (combine == "union") {
      genes.use = union(genes.use, genes.new)
    }
    if (combine == "intersection") {
      genes.use = intersect(genes.use, genes.new)
    }
  }
  if (!keep.unique) {
    for (i in 1:length(object@raw.data)) {
      genes.use = genes.use[genes.use %in% rownames(object@raw.data[[i]])]
    }
  }
  object@var.genes = genes.use
  return(object)
}



selectGenesSeurat_sparse = function (object, alphathresh = 0.99, varthresh = 0.1, cex.use = 0.3)
{
  trx_per_cell <- Matrix::colSums(object@raw.data)
  norm = Matrix.column_norm(object@raw.data)
  gene_expr_mean <- Matrix::rowMeans(norm)
  gene_expr_var <- Matrix::rowMeans(norm^2) - gene_expr_mean^2
  names(gene_expr_var) = rownames(norm)
  nolan_constant <- mean((1/trx_per_cell))
  alphathresh.corrected = alphathresh/dim(object@raw.data)[1]
  genemeanupper <- gene_expr_mean + qnorm(1 - alphathresh.corrected/2) *
    sqrt(gene_expr_mean * nolan_constant/dim(object@raw.data)[2])
  genes.use = names(gene_expr_var)[which(gene_expr_var/nolan_constant >
                                           genemeanupper & log10(gene_expr_var) > log10(gene_expr_mean) +
                                           (log10(nolan_constant) + varthresh))]
  #plot(log10(gene_expr_mean), log10(gene_expr_var), cex = cex.use)
  #points(log10(gene_expr_mean[genes.new]), log10(gene_expr_var[genes.new]),
  #cex = cex.use, col = "green")
  #abline(log10(nolan_constant), 1, col = "purple")
  #legend("bottomright", paste0("Selected genes: ", length(genes.new)),
  #pch = 20, col = "green")
  object@var.genes = genes.use
  return(object)
}

sparse_var = function(x){
  rms = Matrix::rowMeans(x)
  Matrix::rowSums((x-rms)^2)/(dim(x)[2]-1)
}

Sparse_transpose = function(x){
  h = summary(x)
  sparseMatrix(i = h[,2],j=h[,1],x=h[,3])
  
}

plot_gene<-function (object, gene,pt.size=0.5,by.dataset=T){
  gene_vals = c()
  for (i in 1:length(object@norm.data)) {
    gene_vals = c(gene_vals, object@norm.data[[i]][gene,
                                                   ])
  }
  gene_vals = log(10000 * gene_vals + 1)
  gene_df = data.frame(object@tsne.coords)
  rownames(gene_df) = names(object@clusters)
  gene_df$Gene = gene_vals[rownames(gene_df)]
  colnames(gene_df) = c("tSNE1", "tSNE2", gene)
  gene_plots = list()
  if (by.dataset) {
    for (i in 1:length(object@norm.data)) {
      plot_i = (ggplot(gene_df[rownames(object@scale.data[[i]]),
                               ], aes_string(x = "tSNE1", y = "tSNE2", color = gene)) +
                  geom_point(size=pt.size) + scale_color_gradient2(low = "yellow",
                                                                   mid = "red", high = "black", midpoint = (max(gene_vals) -
                                                                                                              min(gene_vals))/2, limits = c(min(gene_vals),
                                                                                                                                            max(gene_vals))) + ggtitle(names(object@scale.data)[i]))
      gene_plots[[i]] = plot_i
    }
  } else {
    
    plot_i = (ggplot(gene_df, aes_string(x = "tSNE1", y = "tSNE2", color = gene)) +
                geom_point(size=pt.size) + scale_color_gradient2(low = "yellow",
                                                                 mid = "red", high = "black", midpoint = (max(gene_vals) -
                                                                                                            min(gene_vals))/2, limits = c(min(gene_vals),
                                                                                                                                          max(gene_vals))) + ggtitle("All Data"))
    gene_plots = list(plot_i)
    
  }
  print(plot_grid(plotlist = gene_plots, ncol = 1))
}

RunUMAP<-function(object, rand.seed = 42,use.raw = F,k=2,distance = 'euclidean',reduction.use="pca",dims.use = 1:5,reduction.name = 'umap',reduction.key = 'UMAP_') {
  
  UMAP<-import("umap")
  umapper = UMAP$UMAP(n_components=as.integer(k),metric = distance)
  Rumap = umapper$fit_transform
  data.use = GetDimReduction(object=object,reduction.type=reduction.use,slot="cell.embeddings")[,dims.use]
  umap.out = Rumap(data.use)
  colnames(x = umap.out) = paste0("UMAP",1:ncol(umap.out))
  rownames(umap.out) = rownames(data.use)
  object = SetDimReduction(object=object,reduction.type=reduction.name,slot='cell.embeddings',new.data = umap.out)
  object = SetDimReduction(object=object,reduction.type=reduction.name,slot='key',new.data = reduction.key)
  return(object)
  
  
}
environment(RunUMAP)<-asNamespace('Seurat')



UMAPPlot<-
  function (object, do.label = FALSE, pt.size = 1, label.size = 4,
            cells.use = NULL, colors.use = NULL, ...)
  {
    return(DimPlot(object = object, reduction.use = "umap", cells.use = cells.use,
                   pt.size = pt.size, do.label = do.label, label.size = label.size,
                   cols.use = colors.use, ...))
  }

environment(UMAPPlot)<-asNamespace('Seurat')

load.multi<-function(base.path = "/broad/macosko/data/libraries/",lib.list,num.nuclei=NULL,lib.names,min.umis = 0) {
  dge.list = lapply(1:length(lib.list),function(i){
    temp = Convert10xtoMatrix(load_cellranger_matrix(paste0(base.path,lib.list[i]),barcode_filtered =F))
    cs = colSums(temp)
    if (!is.null(num.nuclei)) {
      temp = temp[,order(cs,decreasing=T)[1:num.nuclei[i]]]
    }
    temp = temp[,which(colSums(temp)>min.umis)]
    return(temp)
  })
  if (length(dge.list)>1) {
    MergeSparseDataAll(dge.list,lib.names)
  } else {
    return(dge.list[[1]])
  }
}
load.multi.v3<-function(base.path = "/broad/macosko/data/libraries/",lib.list,num.nuclei=NULL,lib.names,min.umis = 0,use.filtered = F, merge = F) {
  dge.list = lapply(1:length(lib.list),function(i){
    print(paste0('Reading in library from: ' ,lib.names[i]))
    if (use.filtered) {
      rawdata = readMM(paste0(base.path,lib.list[i],"/outs/filtered_feature_bc_matrix/matrix.mtx.gz"))
      bcs = read.table(paste0(base.path,lib.list[i],"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),stringsAsFactors =F)[,1]
      gs = read.table(paste0(base.path,lib.list[i],"/outs/filtered_feature_bc_matrix/features.tsv.gz"),stringsAsFactors = F)[,2]
      gs = factor(gs,levels = unique(gs))
      colnames(rawdata) = bcs
    } else {
      rawdata = readMM(paste0(base.path,lib.list[i],"/outs/raw_feature_bc_matrix/matrix.mtx.gz"))
      bcs = read.table(paste0(base.path,lib.list[i],"/outs/raw_feature_bc_matrix/barcodes.tsv.gz"),stringsAsFactors =F)[,1]
      gs = read.table(paste0(base.path,lib.list[i],"/outs/raw_feature_bc_matrix/features.tsv.gz"),stringsAsFactors = F)[,2]
      gs = factor(gs,levels = unique(gs))
      colnames(rawdata) = bcs
    }
    rawdata = rawdata[,which(Matrix::colSums(rawdata)>min.umis)]
    s = summary(rawdata)
    mapping = cbind(1:length(gs),as.numeric(gs))
    i = mapvalues(s[,1],mapping[,1],mapping[,2],warn_missing=F)
    raw.new = sparseMatrix(i= i, j = s[,2],x = s[,3], dims = c(length(unique(gs)),ncol(rawdata)))
    rownames(raw.new) = unique(gs)
    colnames(raw.new) = colnames(rawdata)	
    cs = Matrix::colSums(raw.new)
    if (!is.null(num.nuclei)) {
      if (num.nuclei[i] < ncol(rawdata)) {
        print(warning(paste0("You selected more nuclei than are in output matrix.  Returning all nuclei in matrix.")))
        num.nuclei[i] = ncol(raw.new)
      }
      raw.new = raw.new[,order(cs,decreasing=T)[1:num.nuclei[i]]]
    }
    
    
    return(raw.new)
  })
  if (length(dge.list)>1 & merge == T) {
    MergeSparseDataAll(dge.list,lib.names)
  } else {
    return(dge.list)
  }
}


fastRead<-function(inFile,skip=3,data.table=T,sep="\t",make.sparse =T,header=T,data.type="DGE") {
  if (length(grep(".gz", inFile))==1) {
    t=tempfile()
    cmd=paste("gunzip -c", inFile, ">", t)
    system(cmd,intern=F)
    inFileFinal=t
  } else {
    inFileFinal=inFile
  }
  
  
  a=fread(inFileFinal, data.table=data.table,skip = skip,header = header,sep=sep)
  if (data.type == "DGE") {
    b=a[,-c(1,2)]
    rownames(b) = a$GENE
    a = b
  }
  if (make.sparse) {
    return(DataTableToMatrix(a))
  } else {
    return (a)
  }
}

DataTableToMatrix<-function(d_table) {
  i_list <- lapply(d_table, function(x) which(x != 0))
  counts <- unlist(lapply(i_list, length), use.names = F)
  
  sparseMatrix(
    i = unlist(i_list, use.names = F),
    j = rep(1:ncol(d_table), counts),
    x = unlist(lapply(d_table, function(x) x[x != 0]), use.names = F),
    dims = dim(d_table),
    dimnames = list(rownames(d_table), names(d_table)))
}
mousehumanhomologs<-function(data.mouse,data.human,homolog.file = 
                               "/home/tkamath/scripts/HOM_MouseHumanSequence.rpt",mouse.names = TRUE) {
  homologs = read.table(homolog.file,header=T,stringsAsFactors=F,sep="\t")
  #Go from mouse homologs to human.
  mouse.genes.use = intersect(homologs$Symbol,rownames(data.mouse))
  homolog.ids.mouse = homologs$HomoloGene.ID[which((homologs$Symbol %in% mouse.genes.use) & homologs$Common.Organism.Name!="human")]
  h =homologs$HomoloGene.ID[which(homologs$HomoloGene.ID %in% homolog.ids.mouse)]
  t = table(h)
  #Use only 1:1 matches
  homolog.ids.use = names(t)[which(t==2)]
  h.human = homologs$HomoloGene.ID[intersect(which(homologs$HomoloGene.ID %in% as.numeric(homolog.ids.use)),which(homologs$Common.Organism.Name=="human"))]
  h.mouse =homologs$HomoloGene.ID[intersect(which(homologs$HomoloGene.ID %in% as.numeric(homolog.ids.use)),which(homologs$Common.Organism.Name!="human"))]
  homolog.ids.use = intersect(h.human,h.mouse)
  idx = intersect(which(homologs$HomoloGene.ID %in% as.numeric(homolog.ids.use)),which(homologs$Common.Organism.Name=="human"))
  human.genes.use = homologs$Symbol[idx]
  idx = intersect(which(homologs$HomoloGene.ID %in% as.numeric(homolog.ids.use)),which(homologs$Common.Organism.Name!="human"))
  mouse.genes.use = homologs$Symbol[idx]
  gene.table = cbind(human.genes.use,mouse.genes.use)
  idx = which(gene.table[,1] %in% rownames(data.human))
  gene.table=gene.table[idx,]
  
  mouse.data.use = data.mouse[gene.table[,2],]
  
  human.data.use= data.human[gene.table[,1],]
  if (mouse.names) {
    rownames(human.data.use) = gene.table[,2]
  } else {
    rownames(mouse.data.use) = genes.table[,1]
  }
  data.all = list(MOUSE=mouse.data.use,HUMAN = human.data.use)
  return(data.all)
}

convertGenes <- function(data.1 = data.1, data.2 = data.2, species = c('human','rat'),
                         species.file = '/home/tkamath/scripts/HOM_AllOrganism.rpt'){
  homologs = read.table(species.file,header=T,stringsAsFactors=F,sep="\t", fill = T)
  id.use.1 <- grep(paste0('^',as.character(species[1])),homologs$Common.Organism.Name)
  id.use.2 <- grep(paste0('^',as.character(species[2])),homologs$Common.Organism.Name)
  h1 = homologs[c(id.use.1,id.use.2),]
  h1 <- h1[order(h1$HomoloGene.ID),]
  a1 <- aggregate(nchar(Common.Organism.Name)~HomoloGene.ID,h1,mean)
  gene.ids <- a1[which(a1$`nchar(Common.Organism.Name)` == mean(nchar(unique(h1$Common.Organism.Name)))),]$HomoloGene.ID
  h1 <- h1[h1$HomoloGene.ID %in% gene.ids,]
  gene.table <- data.frame(cbind(h1[grep(paste0('^',as.character(species[1])),h1$Common.Organism.Name),]$Symbol,
                                 h1[grep(paste0('^',as.character(species[2])),h1$Common.Organism.Name),]$Symbol))
  idx.1 <- which(gene.table[,1] %in% rownames(data.1))
  idx.2 <- which(gene.table[,2] %in% rownames(data.2))
  genes.use <- gene.table[intersect(idx.1,idx.2),]
  data.1[match(genes.use[,2], rownames(data.1))]
}


subsetAnalogizer<-function(object,idents.use = NULL,cells.use = NULL) {
  if (!is.null(idents.use)){
    cells.use = names(object@clusters)[which(object@clusters %in% idents.use)]
    
  }
  
  raw.data = lapply(object@raw.data,function(q){
    cells = intersect(cells.use,colnames(q))
    if (length(cells) > 1) {
      q[,cells]
    } else {
      return(NULL)
    }
  })
  raw.data = raw.data[!sapply(raw.data,is.null)]
  nms = names(raw.data)
  a = Analogizer(raw.data)
  
  a@norm.data = lapply(names(a@raw.data),function(i){
    object@norm.data[[i]][,colnames(a@raw.data[[i]])]
    
  })
  a@scale.data = lapply(names(a@raw.data),function(i){
    object@scale.data[[i]][colnames(a@raw.data[[i]]),]
    
  })
  a@H = lapply(names(a@raw.data),function(i){
    object@H[[i]][colnames(a@raw.data[[i]]),]
  })
  a@clusters = object@clusters[unlist(lapply(a@H,rownames))]
  a@clusters = droplevels(a@clusters)
  a@tsne.coords = object@tsne.coords[names(a@clusters),]
  a@H.norm = object@H.norm[names(a@clusters),]
  a@W = object@W
  names(a@scale.data) = names(a@raw.data) = names(a@norm.data) = names(a@H) = nms
  return(a)
}

AnalogizerToSeurat<-
  function (object, need.sparse = T,add.names=F)
  {
    if (add.names) {
      nms = names(object@H)
    } else {
      nms = NULL
    }
    if (need.sparse) {
      object@raw.data = lapply(object@raw.data, function(x) {
        Matrix(as.matrix(x), sparse = T)
      })
      object@norm.data = lapply(object@norm.data, function(x) {
        Matrix(as.matrix(x), sparse = T)
      })
    }
    raw.data = MergeSparseDataAll(object@raw.data, nms)
    norm.data = MergeSparseDataAll(object@norm.data, nms)
    scale.data = do.call(rbind, object@scale.data)
    rownames(scale.data) = colnames(norm.data)
    inmf.obj = new(Class = "dim.reduction", gene.loadings = t(object@W),
                   cell.embeddings = object@H.norm, key = "iNMF")
    tsne.obj = new(Class = "dim.reduction", cell.embeddings = object@tsne.coords,
                   key = "tSNE_")
    rownames(tsne.obj@cell.embeddings) = rownames(scale.data)
    rownames(inmf.obj@cell.embeddings) = rownames(scale.data)
    colnames(tsne.obj@cell.embeddings) = paste0("tSNE_", 1:2)
    new.seurat = CreateSeuratObject(raw.data,min.genes = 0,min.cells=0)
    new.seurat = NormalizeData(new.seurat)
    new.seurat@scale.data = t(scale.data)
    new.seurat@dr$tsne = tsne.obj
    new.seurat@dr$inmf = inmf.obj
    new.seurat = SetIdent(new.seurat, ident.use = as.character(object@clusters))
    return(new.seurat)
  }

AssignByAtlas<-function(obj,ids.ref,knn_k=10) {
  assigned.cells = intersect(names(ids.ref),names(obj@clusters))
  unassigned.cells = setdiff(names(obj@clusters),assigned.cells)
  ids.use = ids.ref[assigned.cells]
  ids.use = droplevels(ids.use)
  nn.10 = get.knnx(obj@H.norm[assigned.cells,],obj@H.norm[unassigned.cells,],k=knn_k)$nn.index
  assignments = apply(nn.10,1,function(n){
    ids = ids.use[n]
    return(ids[which.max(ids)])
    
  })
  ids.all = rep("NA",length(obj@clusters))
  names(ids.all) = names(obj@clusters)
  ids.all[assigned.cells] = as.character(ids.use)
  ids.all[unassigned.cells] = as.character(assignments)
  obj@clusters = factor(ids.all)
  names(obj@clusters) = names(ids.all)
  return(obj)
}

# Function takes in a list of DGEs, with gene rownames and cell colnames, 
# and merges them into a single DGE.
# Also adds library.names to cell.names if expected to be overlap (common with 10X barcodes)
MergeSparseDataAll <- function(datalist, library.names = NULL) {
  
  # Use summary to convert the sparse matrices into three-column indexes where i are the
  # row numbers, j are the column numbers, and x are the nonzero entries
  col_offset <- 0
  allGenes <- unique(unlist(lapply(datalist, rownames)))
  allCells <- c()
  for (i in 1:length(datalist)) {
    print(paste0('Merging matrix: ',as.character(i) ))
    curr <- datalist[[i]]
    curr_s <- summary(curr)
    
    # Now, alter the indexes so that the two 3-column matrices can be properly merged.
    # First, make the current and full column numbers non-overlapping.
    curr_s[, 2] <- curr_s[, 2] + col_offset
    
    # Update full cell names
    if (!is.null(library.names)) {
      cellnames <- paste0(library.names[i], "_", colnames(curr))
    } else {
      cellnames <- colnames(curr)
    }
    allCells <- c(allCells, cellnames)
    
    # Next, change the row (gene) indexes so that they index on the union of the gene sets,
    # so that proper merging can occur.
    idx <- match(rownames(curr), allGenes)
    newgenescurr <- idx[curr_s[, 1]]
    curr_s[, 1] <- newgenescurr
    
    # Now bind the altered 3-column matrices together, and convert into a single sparse matrix.
    if (!exists("full_mat")) {
      full_mat <- curr_s
    } else {
      full_mat <- rbind(full_mat, curr_s)
    }
    col_offset <- length(allCells)
  }
  M <- sparseMatrix(
    i = full_mat[, 1],
    j = full_mat[, 2],
    x = full_mat[, 3],
    dims = c(
      length(allGenes),
      length(allCells)
    ),
    dimnames = list(
      allGenes,
      allCells
    )
  )
  return(M)
}

DotPlot_scaled =
  function (object, genes.plot, cols.use = c("lightgrey", "blue"),
            col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
            scale.by = "radius", scale.min = NA, scale.max = NA, group.by,
            plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE)
  {
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                         radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    if (!missing(x = group.by)) {
      object <- SetAllIdent(object = object, id = group.by)
    }
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    colnames(x = data.to.plot) <- genes.plot
    data.to.plot$cell <- rownames(x = data.to.plot)
    data.to.plot$id <- object@ident
    data.to.plot$dataset <- object@meta.data$orig.ident
    data.to.plot <- data.to.plot %>% gather(key = genes.plot,
                                            value = expression, -c(cell, id,dataset))
    data.to.plot <- data.to.plot %>% group_by(id, genes.plot,dataset) %>%
      summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression,
                                                                              threshold = 0))
    data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot,dataset) %>%
      mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale,
                                                                                   max = col.max, min = col.min))
    data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot,
                                      levels = rev(x = genes.plot))
    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot,
                                                   y = id)) + geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
      scale.func(range = c(0, dot.scale), limits = c(scale.min,
                                                     scale.max)) + theme(axis.title.x = element_blank(),
                                                                         axis.title.y = element_blank())
    if (length(x = cols.use) == 1) {
      p <- p + scale_color_distiller(palette = cols.use)
    }
    else {
      p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
    }
    if (!plot.legend) {
      p <- p + theme(legend.position = "none")
    }
    if (x.lab.rot) {
      p <- p + theme(axis.text.x = element_text(angle = 90,
                                                vjust = 0.5))
    }
    suppressWarnings(print(p))
    if (do.return) {
      return(p)
    }
  }


plotGene = function (object, gene, methylation.indices = NULL, pt.size = 0.1, 
                     min.clip = 0, max.clip = 1, points.only = F, low.col = "yellow", 
                     high.col = "red", return.plots = F,option="plasma",zero_color="#F5F5F5",by.dataset=T) 
{
  gene_vals <- c()
  for (i in 1:length(object@norm.data)) {
    if (i %in% methylation.indices) {
      tmp <- object@norm.data[[i]][gene, ]
      max_v <- quantile(tmp, probs = max.clip, na.rm = T)
      min_v <- quantile(tmp, probs = min.clip, na.rm = T)
      tmp[tmp < min_v & !is.na(tmp)] <- min_v
      tmp[tmp > max_v & !is.na(tmp)] <- max_v
      gene_vals <- c(gene_vals, tmp)
    }
    else {
      if (gene %in% rownames(object@norm.data[[i]])) {
        gene_vals_int <- log2(10000 * object@norm.data[[i]][gene,] + 1)
        gene_vals_int[gene_vals_int==0]=NA
      }
      else {
        gene_vals_int <- rep(list(NA), ncol(object@norm.data[[i]]))
        names(gene_vals_int) <- colnames(object@norm.data[[i]])
      }
      gene_vals <- c(gene_vals, gene_vals_int)
    }
  }
  gene_df <- data.frame(object@tsne.coords)
  rownames(gene_df) <- names(object@clusters)
  gene_df$gene <- as.numeric(gene_vals[rownames(gene_df)])
  colnames(gene_df) <- c("tSNE1", "tSNE2", "gene")
  gene_plots <- list()
  if(by.dataset){
    for (i in 1:length(object@norm.data)) {
      gene_df.sub <- gene_df[rownames(object@scale.data[[i]]),]
      max_v <- max(gene_df.sub["gene"], na.rm = T)
      min_v <- min(gene_df.sub["gene"], na.rm = T)
      midpoint <- (max_v - min_v)/2
      plot_i <- (ggplot(gene_df.sub, aes_string(x = "tSNE1", 
                                                y = "tSNE2", color = "gene")) + geom_point(size = pt.size) + 
                   scale_color_viridis(option=option,direction=-1,na.value=zero_color) + labs(col = gene) + ggtitle(names(object@scale.data)[i]))
      gene_plots[[i]] <- plot_i
    }
  }
  else
  {
    plot_i <- ggplot(gene_df, aes_string(x = "tSNE1", 
                                         y = "tSNE2", color = "gene")) + geom_point(size = pt.size) + 
      scale_color_viridis(option=option,direction=-1,na.value=zero_color) + labs(col = gene) + ggtitle("Full dataset")
    gene_plots[[1]] <- plot_i
  }
  if (points.only) {
    for (i in 1:length(gene_plots)) {
      gene_plots[[i]] <- gene_plots[[i]] + theme(axis.line = element_blank(), 
                                                 axis.text.x = element_blank(), axis.text.y = element_blank(), 
                                                 axis.ticks = element_blank(), axis.title.x = element_blank(), 
                                                 axis.title.y = element_blank(), legend.position = "none", 
                                                 panel.background = element_blank(), panel.border = element_blank(), 
                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                 plot.background = element_blank(), plot.title = element_blank())
    }
  }
  if (return.plots) {
    return(gene_plots)
  }
  else {
    for (i in 1:length(gene_plots)) {
      print(gene_plots[[i]])
    }
  }
}

environment(DotPlot_scaled)<-asNamespace('Seurat')

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


selectGenesLowMem<-function (object, alpha.thresh = 0.99, var.thresh = 0.1, combine = "union", 
                             keep.unique = F, capitalize = F, do.plot = T, cex.use = 0.3) 
{
  if (length(var.thresh) == 1) {
    var.thresh <- rep(var.thresh, length(object@raw.data))
  }
  genes.use <- c()
  for (i in 1:length(object@raw.data)) {
    if (capitalize) {
      rownames(object@raw.data[[i]]) <- toupper(rownames(object@raw.data[[i]]))
      rownames(object@norm.data[[i]]) <- toupper(rownames(object@norm.data[[i]]))
    }
    trx_per_cell <- colSums(object@raw.data[[i]])
    gene_expr_mean <- rowMeans(object@norm.data[[i]])
    chunks = chunk2(1:nrow(object@norm.data[[i]]),10)
    gene_expr_var <- unlist(lapply(chunks,function(q){
      sparse.var(object@norm.data[[i]][q,])
    }))
    names(gene_expr_var) = rownames(object@raw.data[[i]])
    # gene_expr_var <- sparse.var(object@norm.data[[i]])
    nolan_constant <- mean((1/trx_per_cell))
    alphathresh.corrected <- alpha.thresh/nrow(object@raw.data[[i]])
    genemeanupper <- gene_expr_mean + qnorm(1 - alphathresh.corrected/2) * 
      sqrt(gene_expr_mean * nolan_constant/ncol(object@raw.data[[i]]))
    genes.new <- names(gene_expr_var)[which(gene_expr_var/nolan_constant > 
                                              genemeanupper & log10(gene_expr_var) > log10(gene_expr_mean) + 
                                              (log10(nolan_constant) + var.thresh[i]))]
    if (do.plot) {
      plot(log10(gene_expr_mean), log10(gene_expr_var), 
           cex = cex.use, xlab = "Gene Expression Mean (log10)", 
           ylab = "Gene Expression Variance (log10)")
      points(log10(gene_expr_mean[genes.new]), log10(gene_expr_var[genes.new]), 
             cex = cex.use, col = "green")
      abline(log10(nolan_constant), 1, col = "purple")
      legend("bottomright", paste0("Selected genes: ", 
                                   length(genes.new)), pch = 20, col = "green")
      title(main = names(object@raw.data)[i])
    }
    if (combine == "union") {
      genes.use <- union(genes.use, genes.new)
    }
    if (combine == "intersection") {
      genes.use <- intersect(genes.use, genes.new)
    }
  }
  if (!keep.unique) {
    for (i in 1:length(object@raw.data)) {
      genes.use <- genes.use[genes.use %in% rownames(object@raw.data[[i]])]
    }
  }
  object@var.genes <- genes.use
  return(object)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


read10X <- function(sample.dirs, sample.names, merge = T, num.cells = NULL, min.umis = 0,
                    use.filtered = F, reference = NULL) {
  datalist <- list()
  datatypes <- c("Gene Expression")
  
  if (length(num.cells) == 1) {
    num.cells <- rep(num.cells, length(sample.dirs))
  }
  for (i in seq_along(sample.dirs)) {
    print(paste0("Processing sample ", sample.names[i]))
    sample.dir <- sample.dirs[[i]]
    inner1 <- paste0(sample.dir, "/outs")
    if (dir.exists(inner1)) {
      sample.dir <- inner1
      is_v3 <- dir.exists(paste0(sample.dir, "/filtered_feature_bc_matrix"))
      matrix.prefix <- ifelse(use.filtered, "filtered", "raw")
      if (is_v3) {
        sample.dir <- paste0(sample.dir, "/", matrix.prefix, "_feature_bc_matrix")
      } else {
        if (is.null(reference)) {
          references <- list.dirs(paste0(sample.dir, "/raw_gene_bc_matrices"),
                                  full.names = F,
                                  recursive = F
          )
          if (length(references) > 1) {
            stop("Multiple reference genomes found. Please specify a single one.")
          } else {
            reference <- references[1]
          }
        }
        sample.dir <- paste0(sample.dir, "/", matrix.prefix, "_gene_bc_matrices/", reference)
      }
    } else {
      is_v3 <- file.exists(paste0(sample.dir, "/features.tsv.gz"))
    }
    suffix <- ifelse(is_v3, ".gz", "")
    features.file <- ifelse(is_v3, paste0(sample.dir, "/features.tsv.gz"),
                            paste0(sample.dir, "/genes.tsv")
    )
    matrix.file <- paste0(sample.dir, "/matrix.mtx", suffix)
    barcodes.file <- paste0(sample.dir, "/barcodes.tsv", suffix)
    
    rawdata <- readMM(matrix.file)
    # convert to dgc matrix
    if (class(rawdata)[1] == "dgTMatrix") {
      rawdata <- as(rawdata, "CsparseMatrix")
    }
    
    # filter for UMIs first to increase speed
    umi.pass <- which(Matrix::colSums(rawdata) > min.umis)
    if (length(umi.pass) == 0) {
      print("No cells pass UMI cutoff. Please lower it.")
    }
    rawdata <- rawdata[, umi.pass, drop = F]
    
    barcodes <- readLines(barcodes.file)[umi.pass]
    # Remove -1 tag from barcodes
    if (all(grepl(barcodes, pattern = "\\-1$"))) {
      barcodes <- as.vector(sapply(barcodes, function(x) {
        strsplit(x, "-")[[1]][1]
      }))
    }
    
    features <- read.delim(features.file, header = F, stringsAsFactors = F)
    # since some genes are only differentiated by ENSMBL
    rownames(rawdata) <- make.unique(features[, 2])
    colnames(rawdata) <- barcodes
    
    # split based on 10X datatype -- V3 has Gene Expression, Antibody Capture, CRISPR, CUSTOM
    # V2 has only Gene Expression by default and just two columns
    if (ncol(features) < 3) {
      samplelist <- list(rawdata)
      names(samplelist) <- c("Gene Expression")
    } else {
      sam.datatypes <- features[, 3]
      sam.datatypes.unique <- unique(sam.datatypes)
      # keep track of all unique datatypes
      datatypes <- union(datatypes, unique(sam.datatypes))
      samplelist <- lapply(sam.datatypes.unique, function(x) {
        rawdata[which(sam.datatypes == x), ]
      })
      names(samplelist) <- sam.datatypes.unique
    }
    
    # num.cells filter only for gene expression data
    if (!is.null(num.cells)) {
      cs <- Matrix::colSums(samplelist[["Gene Expression"]])
      limit <- ncol(samplelist[["Gene Expression"]])
      if (num.cells[i] > limit) {
        print(paste0(
          "You selected more cells than are in matrix ", i,
          ". Returning all ", limit, " cells."
        ))
        num.cells[i] <- limit
      }
      samplelist[["Gene Expression"]] <- samplelist[["Gene Expression"]][, order(cs, decreasing = T)
                                                                         [1:num.cells[i]]]
    }
    
    datalist[[i]] <- samplelist
  }
  if (merge) {
    print('Merging samples')
    return_dges <- lapply(datatypes, function(x) {
      mergelist <- lapply(datalist, function(d) {
        d[[x]]
      })
      mergelist <- mergelist[!sapply(mergelist, is.null)]
      sample.names.x <- sample.names[!sapply(mergelist, is.null)]
      MergeSparseDataAll(mergelist, sample.names)
    })
    names(return_dges) <- datatypes
    
    # if only one type of data present
    if (length(return_dges) == 1) {
      print(paste0("Returning ", datatypes, " data matrix"))
      return(return_dges[[1]])
    }
    return(return_dges)
  } else {
    names(datalist) <- sample.names
    return(datalist)
  }
}




