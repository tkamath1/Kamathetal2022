
prestowrapper <- function(object, log.FC = log(1.2), all.clusters = F, ident.1 = NA,ident.2 = NA, 
                          metadata.column = NA,one.sided = F, subcluster.use = NA, use.raw = F,use.oldlfc = T){
  if (class(object) == 'seurat'){
    if (use.raw == T){
      X_mat <- object@raw.data
    }
    else {
      X_mat <- object@data
    }
    if (all.clusters == T && is.na(c(ident.1,ident.2) )){
      y_vec <- object@ident
      presto.return <- runpresto(X_mat, y_vec, log.foldchange = log.FC)
      return(presto.return)
    }
    else
      if (all.clusters == F && is.na(metadata.column)){
        if (is.na(ident.2) && (ident.1 %in% object@ident) ){
          y_vec <- rep('ident.2', length = length(colnames(X_mat) ))
          y_vec[which(object@ident %in% ident.1)] <- 'ident.1'
          presto.return <- runpresto(X_mat, y_vec, log.foldchange = log.FC)
          if (one.sided == T){
            presto.return <- presto.return[which(presto.return$group == 'ident.1'),]
          }
          return(presto.return)
        }
        else
          if (all(ident.1 %in% object@ident) && all(ident.2 %in% object@ident) ){
            y1 <- object@ident[object@ident %in% c(ident.1,ident.2)]
            X_mat <- X_mat[,names(y1)]
            y_vec <- rep('ident.2', length = length(colnames(X_mat) ))
            y_vec[which(y1 %in% ident.1)] <- 'ident.1'
            presto.return <- runpresto(X_mat, y_vec, log.foldchange = log.FC)
            if (one.sided == T){
              presto.return <- presto.return[which(presto.return$group == 'ident.1'),]
            }
            return(presto.return)
          }
        else
          stop('Some clusters are not in Seurat ident vector')
      }
    else
      if (metadata.column %in% colnames(object@meta.data)){
        if (all(c(ident.1,ident.2) %in% object@meta.data[[metadata.column]] )){
          if (!is.na(subcluster.use)){
            cells.use <- which(object@ident == subcluster.use)
          }
          else
            cells.use <- c(1:dim(object@meta.data)[1])
          df.use <- object@meta.data[cells.use,]
          #print(df.use)
          X_mat <- X_mat[,rownames(df.use)]
          y_vec <- rep('ident.2', length = length(colnames(X_mat) ))
          y_vec[which(df.use[[metadata.column]] %in% ident.1)] <- 'ident.1'
          presto.return <- runpresto(X_mat, y_vec, log.foldchange = log.FC)
          if (one.sided == T){
            presto.return <- presto.return[which(presto.return$group == 'ident.1'),]
          }
          return(presto.return)
        }
      }
    else
      stop('Metadata column specified not present in Seurat meta.data')
  }
  if (class(object) == 'liger'){
    #log.norm <- lapply(object@norm.data, function(x){
    #  return(log1p(x*10000 + 1))
    #})
    #X_dat <- do.call(rbind,log.norm)
    return('Not implemented yet')
    if (all.clusters == T && is.na(c(ident.1,ident.2) )){
      
      presto.return <- runpresto(X_mat, y_vec, log.foldchange = log.FC)
      return(presto.return)
    }
  }
}

runpresto <- function(X_mat = X_mat, y_vec = y_vec, log.foldchange){
  presto.obj <- wilcoxauc(X_mat, y_vec)
  presto.obj <- presto.obj[which(presto.obj$logFC > log.foldchange),]
  presto.obj <- presto.obj[order(presto.obj[,'group'],presto.obj[,'pval']),]
  presto.obj$statistic <- NULL
  presto.obj$pct.diff <- presto.obj$pct_in - presto.obj$pct_out
  return(presto.obj)
}

