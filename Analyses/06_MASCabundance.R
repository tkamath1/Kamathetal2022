#########################################################
## MASC and differential abundance analysis #############
#########################################################
rm(list = ls())
library(Seurat)
library(rliger)
library(NNLM)
library(dplyr)
library(stringr)
library(qs)
library(metap)
library(lme4)
library(plyr)

setwd('/home/tkamath/DA/ratio_analysis/')
source('/home/tkamath/scripts/extrafuncs.R')
source('/home/tkamath/scripts/SeuratExtrafunctions.R')

# Available at: &&&&
lis.celldata <- qread('/home/tkamath/DA/cleanedmetadata.qs')
pd.meta.use <- qread('/home/tkamath/DA/pdmeta_use.qs')

lis.celldata <- bind_rows(lis.celldata,.id = 'celltype')

# Convert library to nurr
lis.celldata <- lis.celldata %>% mutate(nurr = ifelse(grepl('DAPI',lib),'DAPI','Nurr'))
lis.celldata$nurr <- as.factor(lis.celldata$nurr)

# Replace non-DA cell type with excitatory and inhibitory 
lis.celldata[which(grepl('Ex',lis.celldata$clusters)),]$celltype <- 'Ex'
lis.celldata[which(grepl('Inh',lis.celldata$clusters)),]$celltype <- 'Inh'


t1 <- table(lis.celldata$subject,lis.celldata$celltype)
t2 <- as.data.frame(melt(t1/rowSums(t1)))
t2$status <- as.factor(pd.meta.use[match(t2$Var1,pd.meta.use$Donor.ID),]$disease)
t2 <- t2 %>% mutate(status2 = ifelse(grepl('C',status),'C',as.character(status) ))

setDT(t2)[, 'z' := lapply(.SD[,grep("^val", names(.SD)), with=FALSE], 
                          function(x) x/median(x[status=="Ctrl"])), Var2]

t2$Var2 <- factor(t2$Var2,levels = c('da','Ex','Inh', 'olig','astro','endo','opc','mg'))
levels(t2$Var2) <- c('DA neuron','Ex neuron','Inh neuron','Oligodendrocyte',
                     'Astrocyte','Endothelial cell/pericyte','OPC','Microglia/macrophage')
qsave(t2,'/home/tkamath/DA/ratio_analysis/ratios_allcells.qs')

lapply(unique(t2$Var2),function(x){
  p1 <- wilcox.test(t2[t2$Var2 == x & t2$status2 == 'PD',]$value,
                    t2[t2$Var2 == x & t2$status2 == 'C',]$value)$p.value
  p2 <- wilcox.test(t2[t2$Var2 == x & t2$status2 == 'LBD',]$value,
                    t2[t2$Var2 == x & t2$status2 == 'C',]$value)$p.value
  return(list(p1,p2))
})

# Run MASC
out.ratios <- lapply(unique(lis.celldata$celltype),function(x){
  out.masc <- MASC(dataset = lis.celldata[lis.celldata$celltype ==x,],
                   cluster = lis.celldata[lis.celldata$celltype ==x,]$clusters,
                   contrast = 'disease',random_effects = c('lib'),fixed_effects = c('sex','nurr'),
                   verbose = T)
  
  out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
  colnames(out.masc2)[1] <- 'cluster'
  return(out.masc2)
})
names(out.ratios) <- unique(lis.celldata$celltype)
out.ratios2 <- bind_rows(out.ratios,.id='name')
qsave(out.ratios2,'/home/tkamath/DA/ratio_analysis/mascresults_simplemodel.qs')
out.ratios.use <- out.ratios2[out.ratios2$term =='diseaseDisease', ]
out.ratios.use$color <- rep('black',nrow(out.ratios.use))
out.ratios.use$fdr <- p.adjust(out.ratios.use$p.value,method = 'fdr')
out.ratios.use <- out.ratios.use %>% mutate(color = ifelse(fdr < 0.05, 'red','black'))


# Run MASC for DA subtypes, PD versus LBD
out.ratios.pd <- lapply(unique(lis.celldata$celltype)[8],function(x){
  lis.celldata <- lis.celldata[lis.celldata$status != 'LBD',]
  out.masc <- MASC(dataset = lis.celldata[lis.celldata$celltype ==x,],
                   cluster = lis.celldata[lis.celldata$celltype ==x,]$clusters,
                   contrast = 'disease',random_effects = c('lib'),fixed_effects = c('sex','nurr'),
                   verbose = T)
  
  out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
  colnames(out.masc2)[1] <- 'cluster'
  return(out.masc2)
})
out.ratios.pd <- out.ratios.pd[[1]]
out.ratios.pd <- out.ratios.pd[out.ratios.pd$term == 'diseaseDisease',]
out.ratios.pd$fdr <- p.adjust(out.ratios.pd$p.value,method = 'fdr')
out.ratios.pd$cluster <- substring(out.ratios.pd$cluster,first = 8)
qsave(out.ratios.pd,'/home/tkamath/DA/ratio_analysis/mascresults_PDonly.qs')

out.ratios.lbd <- lapply(unique(lis.celldata$celltype)[8],function(x){
  lis.celldata <- lis.celldata[lis.celldata$status != 'PD',]
  out.masc <- MASC(dataset = lis.celldata[lis.celldata$celltype ==x,],
                   cluster = lis.celldata[lis.celldata$celltype ==x,]$clusters,
                   contrast = 'disease',random_effects = c('lib'),fixed_effects = c('sex','nurr'),
                   verbose = T)
  
  out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
  colnames(out.masc2)[1] <- 'cluster'
  return(out.masc2)
})

out.ratios.lbd <- out.ratios.lbd[[1]]
out.ratios.lbd <- out.ratios.lbd[out.ratios.lbd$term == 'diseaseDisease',]
qsave(out.ratios.lbd,'/home/tkamath/DA/ratio_analysis/mascresults_LBDonly.qs')

# Run without low quality libraries
lis.celldata <- qread('/home/tkamath/DA/cleanedmetadata.qs')
lis.celldata <- bind_rows(lis.celldata,.id = 'celltype')
lis.celldata <- lis.celldata %>% mutate(nurr = ifelse(grepl('DAPI',lib),'DAPI','Nurr'))
lis.celldata$nurr <- as.factor(lis.celldata$nurr)
lis.celldata[which(grepl('Ex',lis.celldata$clusters)),]$celltype <- 'Ex'
lis.celldata[which(grepl('Inh',lis.celldata$clusters)),]$celltype <- 'Inh'

sn.da.annot <- qread('/home/tkamath/DA/da/sn_da_annot_021821.qs')
pd.meta.use <- qread('/home/tkamath/DA/pdmeta_use.qs')
t1 <- aggregate(nUMI ~ lib,sn.da.annot@cell.data,mean)
t1$dataset <- sn.da.annot@cell.data[match(t1$lib,sn.da.annot@cell.data$lib),]$dataset
t1$disease <- pd.meta.use[match(t1$dataset,pd.meta.use$Donor.ID),]$disease

qsave(t1,'/home/tkamath/DA/ratio_analysis/numibydisease_perlib.qs')

`%notin%` <- Negate(`%in%`)
out.ratios2 <- lapply(unique(lis.celldata$celltype)[8],function(x){
  lis.celldata <- lis.celldata[lis.celldata$lib %notin% c('5610DAPI','3298DAPIC','4568DAPIB'),]
  out.masc <- MASC(dataset = lis.celldata[lis.celldata$celltype ==x,],
                   cluster = lis.celldata[lis.celldata$celltype ==x,]$clusters,
                   contrast = 'disease',random_effects = c('lib'),fixed_effects = c('sex','nurr'),
                   verbose = T)
  
  out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
  colnames(out.masc2)[1] <- 'cluster'
  return(out.masc2)
})
out.ratios.use <- bind_rows(out.ratios2,.id = 'names')
out.ratios.use <- out.ratios.use[out.ratios.use$term == 'diseaseDisease',]
out.ratios.use$fdr <- p.adjust(out.ratios.use$p.value,method = 'fdr',n = 10)
out.ratios.use$cluster <- substring(out.ratios.use$cluster,first = 8)

# Without DDT
out.ratios.noddt <- lapply(unique(lis.celldata$celltype)[8],function(x){
  lis.celldata <- lis.celldata[lis.celldata$clusters != 'SOX6_DDT',]
  out.masc <- MASC(dataset = lis.celldata[lis.celldata$celltype ==x,],
                   cluster = lis.celldata[lis.celldata$celltype ==x,]$clusters,
                   contrast = 'disease',random_effects = c('lib'),fixed_effects = c('sex','nurr'),
                   verbose = T)
  
  out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
  colnames(out.masc2)[1] <- 'cluster'
  return(out.masc2)
})
out.ratios.noddt <- bind_rows(out.ratios.noddt,.id = 'names')
out.ratios.noddt <- out.ratios.noddt[out.ratios.noddt$term == 'diseaseDisease',]
out.ratios.noddt$fdr <- p.adjust(out.ratios.noddt$p.value,method= 'fdr')

# Downsample
lis.celldata.da <- lis.celldata[lis.celldata$celltype == 'da',]

out.vals<-lapply(c(1.5,2,3,3.5,4,4.5,5,5.5,6),function(z){
  seed.values<-c(555,745,896,1021,2055)
  p.vec<-c()
  p.vec2 <- c()
  p.vec3 <- c()
  p.vec4 <- c()
  for (i in c(1:length(seed.values)) ){
    set.seed(seed.values[i])
    lis.celldata.ds <- lis.celldata.da[sample(rownames(lis.celldata.da),size = nrow(lis.celldata.da)/z),]
    out.ratios <- lapply(unique(lis.celldata.ds$celltype),function(x){
      out.masc <- MASC(dataset = lis.celldata.ds[lis.celldata.ds$celltype == x,],
                       cluster = lis.celldata.ds[lis.celldata.ds$celltype == x,]$clusters,
                       contrast = 'disease',random_effects = c('lib'),fixed_effects = c('nurr','sex'),
                       verbose = T)
      out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
      colnames(out.masc2)[1] <- 'cluster'
      return(out.masc2)
    })
    names(out.ratios) <- unique(lis.celldata.ds$celltype)
    out.ratios2 <- bind_rows(out.ratios,.id='name')
    out.ratios.ds <- out.ratios2[out.ratios2$term == 'diseaseDisease',]
    out.ratios.ds$fdr<-p.adjust(out.ratios.ds$p.value,method='fdr')
    p.vec[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterSOX6_AGTR1',]$fdr)
    p.vec2[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterSOX6_PART1',]$fdr)
    p.vec3[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterSOX6_GFRA2',]$fdr)
    p.vec4[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterSOX6_DDT',]$fdr)
  }
  return(list('p1' = 1-pnorm(mean(qnorm(1-p.vec))),'p2' = 1-pnorm(mean(qnorm(1-p.vec2))),
              'p3' = 1-pnorm(mean(qnorm(1-p.vec3))),'p4'=1-pnorm(mean(qnorm(1-p.vec4)))  ) )
})

out.vals2 <- as.data.frame(do.call(rbind,out.vals))
out.vals2$p1 <- as.vector(unlist(out.vals2$p1))
out.vals2$p2 <- as.vector(unlist(out.vals2$p2))
out.vals2$p3 <- as.vector(unlist(out.vals2$p3))
out.vals2$p4 <- as.vector(unlist(out.vals2$p4))
vec.use <- c(8464/c(1.5,2,3,3.5,4,4.5,5,5.5,6))
out.vals2 <- melt(out.vals2)
out.vals2$id <- rep(vec.use,4)
levels(out.vals2$variable) <- c('SOX6_AGTR1','SOX6_PART1','SOX6_GFRA2','SOX6_DDT')
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(4)
qsave(out.vals2,'/home/tkamath/DA/ratio_analysis/downsample_analysis.qs')


## For astrocytes
lis.celldata.astro <- lis.celldata[lis.celldata$celltype == 'astro',]

out.vals<-lapply(c(1.5,2,3,3.5,4,4.5,5,5.5,6),function(z){
  seed.values<-c(555,745,896,1021,2055,1,2,3,4,5)
  p.vec<-c()
  p.vec2 <- c()
  p.vec3 <- c()
  p.vec4 <- c()
  for (i in c(1:length(seed.values)) ){
    set.seed(seed.values[i])
    lis.celldata.ds <- lis.celldata.astro[sample(rownames(lis.celldata.astro),size = nrow(lis.celldata.astro)/z),]
    out.ratios <- lapply(unique(lis.celldata.ds$celltype),function(x){
      out.masc <- MASC(dataset = lis.celldata.ds[lis.celldata.ds$celltype == x,],
                       cluster = lis.celldata.ds[lis.celldata.ds$celltype == x,]$clusters,
                       contrast = 'disease',random_effects = c('lib'),fixed_effects = c('nurr','sex'),
                       verbose = T)
      out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
      colnames(out.masc2)[1] <- 'cluster'
      return(out.masc2)
    })
    names(out.ratios) <- unique(lis.celldata.ds$celltype)
    out.ratios2 <- bind_rows(out.ratios,.id='name')
    out.ratios.ds <- out.ratios2[out.ratios2$term == 'diseaseDisease',]
    out.ratios.ds$fdr<-p.adjust(out.ratios.ds$p.value,method='fdr')
    p.vec[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterAstro_VIM_LHX2',]$fdr)
    p.vec2[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterAstro_GJB6_OXTR',]$fdr)
    p.vec3[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterAstro_CYP4F12',]$fdr)
  }
  return(list('p1' = 1-pnorm(mean(qnorm(1-p.vec))),'p2' = 1-pnorm(mean(qnorm(1-p.vec2))),
              'p3' = 1-pnorm(mean(qnorm(1-p.vec3)))  ) )
})

out.vals2 <- as.data.frame(do.call(rbind,out.vals))
out.vals2$p1 <- as.vector(unlist(out.vals2$p1))
out.vals2$p2 <- as.vector(unlist(out.vals2$p2))
out.vals2$p3 <- as.vector(unlist(out.vals2$p3))
vec.use <- c(nrow(lis.celldata.astro)/c(1.5,2,3,3.5,4,4.5,5,5.5,6))
out.vals2 <- melt(out.vals2)
out.vals2$id <- rep(vec.use,3)
levels(out.vals2$variable) <- c('Astro_VIM_LHX2','Astro_GJB6_OXTR','Astro_CYP4F12')
qsave(out.vals2,'/home/tkamath/DA/ratio_analysis/downsample_analysis_astro.qs')


#microglia
lis.celldata.mg <- lis.celldata[lis.celldata$celltype == 'mg',]

out.vals<-lapply(c(1.5,2,3,3.5,4,4.5,5,5.5,6),function(z){
  seed.values<-c(555,745,896,1021,2055,1,2,3,4,5)
  p.vec<-c()
  p.vec2 <- c()
  for (i in c(1:length(seed.values)) ){
    set.seed(seed.values[i])
    lis.celldata.ds <- lis.celldata.mg[sample(rownames(lis.celldata.mg),size = nrow(lis.celldata.mg)/z),]
    out.ratios <- lapply(unique(lis.celldata.ds$celltype),function(x){
      out.masc <- MASC(dataset = lis.celldata.ds[lis.celldata.ds$celltype == x,],
                       cluster = lis.celldata.ds[lis.celldata.ds$celltype == x,]$clusters,
                       contrast = 'disease',random_effects = c('lib'),fixed_effects = c('nurr','sex'),
                       verbose = T)
      out.masc2 <- bind_rows(lapply(out.masc,function(x){x[[1]]}),.id='name')
      colnames(out.masc2)[1] <- 'cluster'
      return(out.masc2)
    })
    names(out.ratios) <- unique(lis.celldata.ds$celltype)
    out.ratios2 <- bind_rows(out.ratios,.id='name')
    out.ratios.ds <- out.ratios2[out.ratios2$term == 'diseaseDisease',]
    out.ratios.ds$fdr<-p.adjust(out.ratios.ds$p.value,method='fdr')
    p.vec[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterMG_GPNMB_SULT1C2',]$fdr)
    p.vec2[i] <- as.numeric(out.ratios.ds[out.ratios.ds$cluster=='clusterMG_MGAM',]$fdr)
  }
  return(list('p1' = 1-pnorm(mean(qnorm(1-p.vec))),
              'p2' = 1-pnorm(mean(qnorm(1-p.vec2))) ) )
})

out.vals2 <- as.data.frame(do.call(rbind,out.vals))
out.vals2$p1 <- as.vector(unlist(out.vals2$p1))
out.vals2$p2 <- as.vector(unlist(out.vals2$p2))
vec.use <- c(nrow(lis.celldata.mg)/c(1.5,2,3,3.5,4,4.5,5,5.5,6))
out.vals2 <- reshape2::melt(out.vals2)
out.vals2$id <- rep(vec.use,2)
levels(out.vals2$variable) <- c('MG_GPNMB_SULT1C2','MG_MGAM')
qsave(out.vals2,'/home/tkamath/DA/ratio_analysis/downsample_analysis_mg.qs')
