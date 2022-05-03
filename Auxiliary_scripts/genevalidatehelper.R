library(rliger)
library(qs)

plot.by.dataset <-function(lig.obj,genes,sample.sheet = NULL) {
  for (i in genes) {
    x = rliger::plotGene(lig.obj,gene=i,return.plots = T,plot.by = 'dataset',
                         points.only = T,scale.by = 'none')
    if (!is.null(sample.sheet)){
      x = lapply(1:length(x),function(q){
        x[[q]]+theme(axis.line = element_blank(), axis.text.x = element_blank(),
                     axis.text.y = element_blank(), axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(), legend.position = "none",
                     panel.background = element_blank(), panel.border = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     plot.background = element_blank()) + 
          ggtitle(paste0(sample.sheet[
            match(names(lig.obj@norm.data[q]),sample.sheet$Code1),]$Status,'_',names(lig.obj@norm.data)[q]))
      })
      names(x) <- names(lig.obj@norm.data)
      
      x <- x[order(sample.sheet[match(names(x),sample.sheet$Code1),]$Status)]
    }
    else
      x = lapply(1:length(x),function(q){
        x[[q]]+theme(axis.line = element_blank(), axis.text.x = element_blank(),
                     axis.text.y = element_blank(), axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(), legend.position = "none",
                     panel.background = element_blank(), panel.border = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     plot.background = element_blank()) + 
          ggtitle(names(lig.obj@norm.data)[q])
      })
    grid.arrange(grobs = x,top=i)
  }
}
multi.plot<-function(object,genes.plot) {
  l= lapply(genes.plot,function(q){
    p = rliger::plotGene(object,q,points.only = T,return.plots = T)
    p[[1]] + theme(plot.title = element_text()) + ggtitle(q)
  })
  grid.arrange(grobs = l)
}

plotAvgViolin <- function(lig.obj, gene.use, celltype.use,sample.sheet,label.size,pt.size,return = F){
  out.gene <- lapply(lig.obj@norm.data,function(x,gene = gene.use,celltype = celltype.use,
                                                object.use = lig.obj){
    if (!is.na(celltype)){
      idx.use <- intersect(colnames(x),names(object.use@clusters[which(object.use@clusters == celltype)]))
      if (length(idx.use) == 0){
        return(0)
      }
    }
    else
      idx.use <- colnames(x)
    if(gene %in% rownames(x) )
      return(mean(x[gene,idx.use] > 0))
    else
      return(0)
  })
  df.plot <- data.frame('gene' = as.vector(unlist(out.gene)),'dataset' = names(out.gene), 
                        'status' = sample.sheet[match(names(out.gene),sample.sheet$Code1),]$Status)
  df.plot$status <- factor(df.plot$status,levels = c('Ctrl','Abeta','AbetaTau'))
  g1 <- ggplot(df.plot,aes(x = status, y =  gene,fill = status)) + geom_point(size = pt.size) +
    theme(legend.position = 'none') + geom_boxplot(alpha = 0.2) + 
    geom_text_repel(data= df.plot, aes(x = status, y =  gene, label = dataset,size = label.size)) +
    ggtitle(gene.use)
  if (return == F){
    return(g1)
  }
  else
    return(df.plot)
}

CombineMAST <- function(mast.results,component.pvalue.use, component.coef.use,correction,contrast.use){
  mast.results <- lapply(mast.results,function(x){
    if (contrast.use %in% x[['contrast']]){
      x.out <- merge(x[contrast==contrast.use & component==component.pvalue.use,.(primerid, `Pr(>Chisq)`)], 
                     x[contrast==contrast.use & component==component.coef.use, .(primerid, coef, ci.hi, ci.lo,z)],
                     by='primerid')
      x.out[,fdr:=p.adjust(`Pr(>Chisq)`, correction)]
      return(x.out)
    }
    else{
      return(NULL)
    }
  })
  return(mast.results)
}
