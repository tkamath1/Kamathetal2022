## Loading in and cleaning individual datasets
# Last run on 02/15/21
library(Seurat)
library(rliger)
library(qs)
setwd('~/DA/PDpaper_reproducibility/')
source('SeuratExtrafunctions.R')
source('extrafuncs.R')

## Note: these scripts assume you have all libraries loaded into the relative path: ~/data/libraries/
## Note 2: due to random initialization of t-SNE/UMAP, 
###the resulting data can look slightly different depending on seed

########################
###########3298#########
########################
pd.data <- read.table('~/data/3298.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1,
                 data.names = c('3298DAPIA','3298DAPIB',"3298DAPIC","3298NurrA","3298NurrB"),
                 min.umis = c(900,900,900,900,900))
C3298 = MergeSparseDataAll(t1)
C3298 = CreateSeuratObject(C3298)
C3298 = NormalizeData((C3298))
C3298 = FindVariableGenes(C3298)
print(length(C3298@var.genes))
C3298= ScaleData(C3298)
C3298= RunPCA(C3298, pcs.compute = 30)
C3298= RunTSNE(C3298, dims.use = 1:30, reduction.use= 'pca')
C3298 <- FindClusters(object = C3298, reduction.type = 'pca', dims.use = 1:30,
                   resolution = 0.9, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(C3298@ident) <- c('Oligo','Oligo','NonDA','Astro','NonDA','OPC',
                     'MG','Oligo','Astro','NonDA','Endofibro','Astro','Oligo',
                     'DA','NonDA','Endofibro','REMOVE','NonDA','Endofibro','PBMC',
                     'NonDA','PBMC','Endofibro','NonDA')
# Save file
#qsave(C3298,'3298_new.qs')

########################
###########3322#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/C3322/3322.txt')
t1 <- loadsclibs(lib.directory = '/home/aabdulra/JilongLibs/',data.libs = pd.data$V1[1:11],
                 data.names = c('3322DAPIA','3322DAPIB',"3322DAPIC","3322NurrA","3322NurrB","3322NurrC",
                                "3322NurrD","3322NurrE","3322NurrF","3322NurrG", "3322NurrH"),
                 min.umis = c(1000,1000,1000,1000,1000,1000,1000,1000, 1000,1000,1000))
C3322 = MergeSparseDataAll(t1)
C3322 = CreateSeuratObject(C3322)
C3322 = NormalizeData(C3322)
C3322 = FindVariableGenes(C3322)
print(length(C3322@var.genes))
qsave(C3322,'C3322_new.qs')
C3322 = ScaleData(C3322)
C3322 = RunPCA(C3322, pcs.compute = 25)
C3322 = RunTSNE(C3322,dims.use = 1:25,reduction.use='pca')
C3322  <- FindClusters(object = C3322, reduction.type = 'pca', dims.use = 1:25,
                       print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
C3322@scale.data <- NULL

qsave(C3322,'C3322_new.qs')

levels(C3322@ident) <- c('NonDA','NonDA','NonDA','Oligo','NonDA',
                         'NonDA','Oligo','NonDA','NonDA','Oligo',
                         'NonDA','Oligo','MG','MG','Endofibro',
                         'Oligo','Astro','Oligo','NonDA','NonDA',
                         'OPC','NonDA','REMOVE','Endofibro','Astro',
                         'Endofibro','DA','Oligo','NonDA','NonDA',
                         'OPC','PBMC','PBMC','NonDA','Astro','MG')

#qsave(C3322,'C3322_new.qs')

########################
###########3345#########
########################
pd.data <- read.table('/home/aabdulra/SN2020/IndividualAnalysisMay2020/C3345/C3345.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:4], 
                 data.names = c("3345DAPIA","3345DAPIB","3345DAPIC", "3345NurrA"),
                 min.umis = c(1000,1000,1000, 1000))
C3345 = MergeSparseDataAll(t1)
C3345 = CreateSeuratObject(C3345)
C3345 = NormalizeData((C3345))
C3345 = FindVariableGenes((C3345))
print(length(C3345@var.genes))
C3345 = ScaleData((C3345))
C3345 = RunPCA(C3345, pcs.compute = 30)
C3345 = RunTSNE(C3345,dims.use = 1:30,reduction.use='pca')
C3345  <- FindClusters(object = C3345, reduction.type = 'pca', dims.use = 1:30,
                       resolution = .7, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
C3345@scale.data <- NULL
levels(C3345@ident) <- c('Oligo','Oligo','Oligo','MG','Astro',
                         'DA','NonDA','OPC','DA','NonDA',
                         'NonDA','Endofibro','MG','NonDA','DA',
                         'NonDA','MG','PBMC','REMOVE','PBMC',
                         'REMOVE','REMOVE')
#qsave(C3345,'C3345_new.qs')

########################
###########3346#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/C3346/C3346.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:4],
                 data.names = c("3346DAPIB","3346NurrA","3346NurrB","3346NurrC"),
                 min.umis = c(2200,2700,2500,1800))

C3346 = MergeSparseDataAll(t1)
C3346 = CreateSeuratObject(C3346)
C3346 = NormalizeData((C3346))
C3346 = FindVariableGenes((C3346))
print(length(C3346@var.genes))
C3346 = ScaleData((C3346))
C3346 = RunPCA(C3346, pcs.compute = 30)
C3346 = RunTSNE(C3346,dims.use = 1:30,reduction.use='pca')
C3346  <- FindClusters(object = C3346, reduction.type = 'pca', dims.use = 1:30,
                       resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

levels(C3346@ident) <- c('MG','NonDA','NonDA','NonDA','Oligo',
                         'Oligo','Oligo','Oligo','DA','NonDA',
                         'NonDA','Astro','DA','NonDA','Endofibro',
                         'NonDA','NonDA','NonDA','Oligo','NonDA',
                         'NonDA','NonDA','OPC','Endofibro','NonDA',
                         'Oligo','NonDA','NonDA','NonDA','Endofibro',
                         'MG','MG','PBMC','Astro')
#qsave(C3346,'C3346_new.qs')

########################
###########3482#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/C3482/C3482.txt')
t1 <- loadsclibs(lib.directory = "~/JilongLibs/",data.libs = pd.data$V1[1:4],
                 data.names = c("3482DAPIA","3482DAPIB","3482NurrA","3482NurrB"),
                 min.umis = c(1200,1200,2600,3000))
C3482 = MergeSparseDataAll(t1)
C3482 = CreateSeuratObject(C3482)
C3482 = NormalizeData((C3482))
C3482 = FindVariableGenes((C3482))
C3482 = ScaleData((C3482))
C3482 = RunPCA(C3482, pcs.compute = 30)
C3482 = RunTSNE(C3482,dims.use = 1:30,reduction.use='pca')
C3482  <- FindClusters(object = C3482, reduction.type = 'pca', dims.use = 1:30,
                       resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

levels(C3482@ident) <- c('Oligo','Oligo','Oligo','Astro','MG',
                         'NonDA','DA','OPC','DA','Oligo',
                         'Endofibro','Oligo','NonDA','NonDA','NonDA',
                         'DA','NonDA','NonDA','Oligo','NonDA',
                         'Endofibro','Endofibro','Oligo','DA', 'NonDA',
                         'NonDA','MG','PBMC','Astro','NonDA',
                         'NonDA','NonDA','NonDA','REMOVE')
#qsave(C3482,'C3482_new.qs')

########################
###########4956#########
########################
pd.data <- read.table('/home/aabdulra/SN2020/IndividualAnalysisMay2020/C4956/C4956.txt')
t2 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:4],
                 data.names = c("4956DAPIA","4956DAPIB","4956NurrA","4956NurrB"),
                 min.umis = c(1000,1000,1000,1000))
C4956 = MergeSparseDataAll(t1)
C4956 = CreateSeuratObject(C4956)
C4956 = NormalizeData((C4956))
C4956 = FindVariableGenes((C4956))
print(length(C4956@var.genes))
C4956 = ScaleData((C4956))
C4956 = RunPCA(C4956, pcs.compute = 30)
C4956 = RunTSNE(C4956,dims.use = 1:30,reduction.use='pca')
C4956  <- FindClusters(object = C4956, reduction.type = 'pca', dims.use = 1:30,
                       print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(C4956@ident) <- c('Oligo','Astro','DA','DA','NonDA',
                         'MG','Oligo','DA','MG','Oligo',
                         'NonDA','Oligo','Endofibro','NonDA','Astro',
                         'DA','DA','NonDA','DA','NonDA',
                         'NonDA','Oligo','DA','MG','NonDA',
                         'MG','NonDA','NonDA','OPC','MG',
                         'PBMC','MG','NonDA')
#qsave(C4956,'C4956_new.qs')

########################
###########5610#########
########################
pd.data <- read.table('/home/aabdulra/SN2020/IndividualAnalysisMay2020/C5610/C5610.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:2],
                 data.names = c("5610DAPI","5610Nurr"),min.umis = c(700,800))

C5610 = MergeSparseDataAll(t1)
C5610 = CreateSeuratObject(C5610)
C5610 = NormalizeData((C5610))
C5610 = FindVariableGenes((C5610))
print(length(C5610@var.genes))
C5610 = ScaleData((C5610))
C5610 = RunPCA(C5610, pcs.compute = 30)
C5610 = RunTSNE(C5610,dims.use = 1:30,reduction.use='pca')
C5610  <- FindClusters(object = C5610, reduction.type = 'pca', dims.use = 1:30,
                       resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(C5610@ident) <- c('DA','DA','DA','REMOVE','Oligo',
                         'NonDA','Oligo','Astro','Oligo','REMOVE',
                         'MG','DA','NonDA','DA','NonDA',
                         'OPC','Endofibro','Astro','Astro','NonDA',
                         'Endofibro','PBMC','NonDA','Oligo','NonDA')
#qsave(C5610,'C5610_new.qs')

########################
###########6173#########
########################
pd.data <- read.table('/home/aabdulra/SN2020/IndividualAnalysisMay2020/C6173/C6173.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:3],
                 data.names = c("6173DAPI","6173NurrA","6173NurrB"),
                 min.umis = c(700,600,600))

C6173 = MergeSparseDataAll(t1)
C6173 = CreateSeuratObject(C6173)
C6173 = NormalizeData((C6173))
C6173 = FindVariableGenes((C6173))
print(length(C6173@var.genes))
C6173 = ScaleData((C6173))
C6173 = RunPCA(C6173, pcs.compute = 30)
C6173 = RunTSNE(C6173,dims.use = 1:30,reduction.use='pca')
C6173 <- FindClusters(object = C6173, reduction.type = 'pca', dims.use = 1:30,
                      resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

levels(C6173@ident) <- c('DA','Oligo','Oligo','Oligo','MG',
                         'Astro','NonDA','Oligo','NonDA','Oligo',
                         'NonDA','Oligo','OPC','NonDA','Oligo',
                         'DA','Endofibro','MG','Endofibro','NonDA',
                         'NonDA','PBMC','NonDA')
#qsave(C6173,'C6173_new.qs')

########################
###########2544#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/LB2544/LB2544.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:5],
                 data.names = c("2544DAPIA","2544DAPIB","2544NurrA","2544NurrB","2544NurrC"),
                 min.umis = c(600,1100,2200,2000,1000))
LB2544 = MergeSparseDataAll(t1)
LB2544 = CreateSeuratObject(LB2544)
LB2544 = NormalizeData((LB2544))
LB2544 = FindVariableGenes((LB2544))
print(length(LB2544@var.genes))
LB2544 = ScaleData((LB2544))
LB2544 = RunPCA(LB2544, pcs.compute = 30)
LB2544 = RunTSNE(LB2544,dims.use = 1:30,reduction.use='pca')
LB2544  <- FindClusters(object = LB2544, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(LB2544@ident) <- c('Oligo','Oligo','Oligo','NonDA','MG',
                          'DA','Astro','Astro','MG','DA',
                          'NonDA','OPC','NonDA','NonDA','Oligo',
                          'Oligo','Oligo','NonDA','Endofibro','Endofibro',
                          'NonDA','NonDA','Endofibro','PBMC','MG')
#qsave(LB2544,'LB2544_new.qs')

########################
###########2561#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/LB2561/LB2561.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:5],
                 data.names = c("LB2561DAPIA","LB2561DAPIB","LB2561NurrA", "LB2561NurrB","LB2561NurrC"),
                 min.umis = c(900,900,1000,1000, 1100))
LB2561 = MergeSparseDataAll(t1)
LB2561 = CreateSeuratObject(LB2561)
LB2561 = NormalizeData((LB2561))
LB2561 = FindVariableGenes(LB2561,varthresh = 0.3)
print(length(LB2561@var.genes))
LB2561 = ScaleData(LB2561)
LB2561 = RunPCA(LB2561, pcs.compute = 30)
LB2561 = RunTSNE(LB2561, dims.use = 1:30, reduction.use='pca')
LB2561  <- FindClusters(object = LB2561, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(LB2561@ident) <- c('NonDA','NonDA','NonDA','Oligo','NonDA',
                          'NonDA','Oligo','Astro','NonDA','NonDA',
                          'Oligo','NonDA','Oligo','OPC','NonDA',
                          'MG','MG','Oligo','Endofibro','NonDA',
                          'NonDA','Endofibro','NonDA','NonDA','NonDA',
                          'DA','PBMC','NonDA')
#qsave(LB2561,'LB2561_new.qs')

########################
###########2569#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/LB2569/LB2569.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:3],
                 data.names = c("2569DAPIA","2569DAPIB","2569Nurr"),
                 min.umis = c(800,800,800))
LB2569 = MergeSparseDataAll(t1)
LB2569 = CreateSeuratObject(LB2569)
LB2569 = NormalizeData((LB2569))
LB2569 = FindVariableGenes((LB2569))
print(length(LB2569@var.genes))
LB2569 = ScaleData((LB2569))
LB2569 = RunPCA(LB2569, pcs.compute = 30)
LB2569 = RunTSNE(LB2569,dims.use = 1:30,reduction.use='pca')
LB2569  <- FindClusters(object = LB2569, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(LB2569@ident) <- c('Oligo','Astro','MG','Oligo','NonDA',
                          'DA','Oligo','Oligo','Oligo','DA',
                          'OPC','MG','Astro','Endofibro','NonDA',
                          'Endofibro','DA','NonDA','NonDA','PBMC',
                          'Endofibro','Astro','Oligo','NonDA','REMOVE')
#qsave(LB2569,'LB2569_new.qs')

########################
###########4775#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/LB4775/LB4775.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:5],
                 data.names = c("4775DAPIA","4775DAPIB", "4775DAPIC", "4775NurrA","4775NurrB"),
                 min.umis = c(1000,1000,1000, 1400, 2200))
LB4775 = MergeSparseDataAll(t1)
LB4775 = CreateSeuratObject(LB4775)
LB4775 = NormalizeData((LB4775))
LB4775 = FindVariableGenes((LB4775))
print(length(LB4775@var.genes))
LB4775= ScaleData((LB4775))
LB4775= RunPCA(LB4775, pcs.compute = 30)
LB4775= RunTSNE(LB4775,dims.use = 1:30,reduction.use='pca')
LB4775 <- FindClusters(object = LB4775, reduction.type = 'pca', dims.use = 1:30,
                       resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(LB4775@ident) <- c('Oligo','Oligo','Oligo','Oligo','Oligo',
                          'MG','Astro','NonDA','NonDA','OPC',
                          'NonDA','MG','Astro','NonDA','Astro',
                          'Astro','MG','Astro','NonDA','Astro',
                          'NonDA','Endofibro','Endofibro','DA','NonDA',
                          'MG','Endofibro','NonDA','PBMC','Endofibro',
                          'MG','PBMC')
#qsave(LB4775,'LB4775_new.qs')

########################
###########1963#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/PD1963/PD1963.txt')
pd.1963 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:8],
                      data.names = c("1963DAPIA","1963DAPIB","1963DAPIC","1963NurrA","1963NurrB",
                                     "1963NurrC","1963NurrD","1963NurrE"),
                      min.umis = c(2100,2100,2100,2700,2700, 2700, 2700, 2700))
PD1963 = MergeSparseDataAll(pd.1963)
PD1963 = CreateSeuratObject(PD1963)
PD1963 = NormalizeData((PD1963))
PD1963 = FindVariableGenes((PD1963))
print(length(PD1963@var.genes))
PD1963= ScaleData((PD1963))
PD1963= RunPCA(PD1963, pcs.compute = 30)
PD1963= RunTSNE(PD1963,dims.use = 1:30,reduction.use='pca')
PD1963 <- FindClusters(object = PD1963, reduction.type = 'pca', dims.use = 1:30,
                       resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(PD1963@ident) <- c('NonDA','Oligo','Oligo','Oligo','Astro',
                          'MG','NonDA','NonDA','NonDA','Oligo',
                          'NonDA','Astro','Astro','OPC','NonDA',
                          'NonDA','DA','NonDA','Endofibro','Endofibro',
                          'NonDA','Astro','DA','Astro','Endofibro',
                          'Oligo','PBMC','MG','NonDA')
#qsave(PD1963,'PD1963_new.qs')

########################
###########2142#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/PD2142/PD2142.txt')
pd.2142 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:3], 
                      data.names = c("2142DAPIA","2142DAPIB","2142NurrA"),
                      min.umis = c(600,600,2000))
PD2142 = MergeSparseDataAll(pd.2142)
PD2142 <- CreateSeuratObject(PD2142)
PD2142 = NormalizeData((PD2142))
PD2142 = FindVariableGenes((PD2142))
print(length(PD2142@var.genes))
PD2142 = ScaleData((PD2142))
PD2142 = RunPCA(PD2142, pcs.compute = 30)
PD2142 = RunTSNE(PD2142,dims.use = 1:30,reduction.use='pca')
PD2142  <- FindClusters(object = PD2142, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(PD2142@ident) <- c('Oligo','Oligo','NonDA','NonDA','Oligo',
                          'NonDA','NonDA','MG','Astro','Oligo',
                          'Oligo','OPC','DA','Endofibro','NonDA')
#qsave(PD2142,'PD2142_new.qs')

########################
###########3873#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/PD3873/PD3873.txt')
pd.3873 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:8], 
                      data.names = c("3873DAPIA","3873DAPIB","3873NurrA","3873NurrB","3873NurrC",
                                     "3873NurrD","3873NurrE","3873NurrF"),
                      min.umis = c(700,900,1200,1000,1000, 1000, 1300, 1100))
PD3873 = MergeSparseDataAll(pd.3873)
PD3873 = CreateSeuratObject(PD3873)
PD3873 = NormalizeData(PD3873)
PD3873 = FindVariableGenes(PD3873)
print(length(PD3873@var.genes))
PD3873 = ScaleData(PD3873)
PD3873 = RunPCA(PD3873, pcs.compute = 30)
PD3873 = RunTSNE(PD3873,dims.use = 1:30,reduction.use='pca')
PD3873  <- FindClusters(object = PD3873, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(PD3873@ident) <- c('REMOVE','REMOVE','NonDA','NonDA','NonDA',
                          'Oligo','Oligo','Oligo','NonDA','NonDA',
                          'NonDA','Astro','NonDA','NonDA','MG',
                          'DA','Oligo','Astro','NonDA','OPC',
                          'NonDA','Oligo','Endofibro','Endofibro','REMOVE',
                          'NonDA','MG','PBMC','OPC')
#qsave(PD3873,'PD3873_new.qs')

########################
###########3887#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/PD3887/PD3887.txt')
pd.3887 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:4], 
                      data.names = c("PD3887DAPIA","PD3887DAPIB","PD3887NurrA","PD3887NurrB"),
                      min.umis = c(1700,2000,1000,1000))
PD3887 = MergeSparseDataAll(pd.3887)
PD3887 = CreateSeuratObject(PD3887)
PD3887 = NormalizeData((PD3887))
PD3887 = FindVariableGenes((PD3887))
print(length(PD3887@var.genes))
PD3887 = ScaleData((PD3887))
PD3887 = RunPCA(PD3887, pcs.compute = 30)
PD3887 = RunTSNE(PD3887,dims.use = 1:30,reduction.use='pca')
PD3887  <- FindClusters(object = PD3887, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(PD3887@ident) <- c('NonDA','NonDA','Oligo','NonDA','Oligo',
                          'MG','NonDA','NonDA','NonDA','Oligo',
                          'Astro','Oligo','OPC','NonDA','Astro',
                          'NonDA','Endofibro','MG','REMOVE','Endofibro',
                          'PBMC','NonDA','DA','MG')
#qsave(PD3887,'PD3887_new.qs')

########################
###########4560#########
########################
pd.data <- read.table('~/SN2020/IndividualAnalysisMay2020/PD4560/PD4560.txt')
t1 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:7], 
                 data.names = c("4560DAPIA","4560DAPIB","4560NurrA","4560NurrB","4560NurrC",
                                "4560NurrD", "4560NurrE"),
                 min.umis = c(1100,1100,1000,1000,1000, 900, 1000))
PD4560 = MergeSparseDataAll(t1)
PD4560 = CreateSeuratObject(PD4560)
PD4560 = NormalizeData((PD4560))
PD4560 = FindVariableGenes((PD4560))
print(length(PD4560@var.genes))
PD4560 = ScaleData((PD4560))
PD4560 = RunPCA(PD4560, pcs.compute = 30)
PD4560 = RunTSNE(PD4560,dims.use = 1:30,reduction.use='pca')
PD4560  <- FindClusters(object = PD4560, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(PD4560@ident) <- c('REMOVE','Oligo','NonDA','Oligo','NonDA',
                          'Oligo','Oligo','Astro','NonDA','DA',
                          'Oligo','MG','NonDA','Endofibro','MG',
                          'NonDA','OPC','NonDA','NonDA','NonDA',
                          'NonDA','REMOVE','NonDA','NonDA','Astro',
                          'NonDA','NonDA','NonDA','NonDA','Astro',
                          'Astro','PBMC','PBMC','MG')
#qsave(PD4560,'PD4560_new.qs')

########################
###########4568#########
########################
pd.data <- read.table('/home/aabdulra/SN2020/IndividualAnalysisMay2020/PD4568/PD4568.txt')
pd.4568 <- loadsclibs(lib.directory = "/home/aabdulra/JilongLibs/",data.libs = pd.data$V1[1:4], 
                      data.names = c("4568DAPIA","4568DAPIB","4568NurrA","4568NurrB"),
                      min.umis = c(700,900,1000,1000))
PD4568 = MergeSparseDataAll(pd.4568)
PD4568 = CreateSeuratObject(PD4568)
PD4568 = NormalizeData((PD4568))
PD4568 = FindVariableGenes(PD4568,varthresh = 0.25)
print(length(PD4568@var.genes))
PD4568 = ScaleData((PD4568))
PD4568 = RunPCA(PD4568, pcs.compute = 30)
PD4568 = RunTSNE(PD4568,dims.use = 1:30,reduction.use='pca')
PD4568  <- FindClusters(object = PD4568, reduction.type = 'pca', dims.use = 1:30,
                        resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
levels(PD4568@ident) <- c('Oligo','Oligo','Astro','Oligo','NonDA',
                          'Oligo','MG','NonDA','Oligo','OPC',
                          'NonDA','Endofibro','Oligo','NonDA','Endofibro',
                          'DA','NonDA','NonDA','NonDA','NonDA',
                          'Endofibro','Endofibro','NonDA','Endofibro','REMOVE',
                          'Astro','PBMC','NonDA')
#qsave(PD4568,'PD4568_new.qs')
