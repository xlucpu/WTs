#------------------------#
# PROJECT: WILMS PROJECY #
#    START: 2020/4/21    #

workdir <- "E:/IGBMC/myproject/Wilms"; setwd(workdir)
tumor.path <- workdir; tumorname <- "Wilms"
fig.path <- file.path(workdir,"Figures")
res.path <- file.path(workdir,"Results")
data.path <- file.path(workdir,"InputData")

if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }
if (!file.exists(data.path)) { dir.create(data.path) }

GinfoFile <- "overlapTable.txt"  ## Gene annotation file created by Jack
comRFun.path <- file.path(tumor.path,"commonFun") ## universal functions
shareFun.path <- file.path(tumor.path,"SharedScripts") ## functions shared by this type of analysis only
script.path <- file.path(tumor.path,"Scripts")
comAnn.path <- file.path(tumor.path,"Annotation") ## path containning annotation files created by Jack
Ginfo <- read.table(file.path(comAnn.path,GinfoFile),sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
createGinfoFileFlag <- F ## to create Gene annotation file FALSE by default
options(warn =0)
tailrows <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique") ## tail rows not related to gene names in read counts input file
lowcut.mRNA <- 1
lowcut.lnc <- 0.25 # FPKM cutoff for low expressed genes
lowpct <- 0.9 # if a gene is below lowcut in 90% of samples, then it will be not be included for the analysis

source(file.path(shareFun.path,"LoadCommonFunLib.R")) 
source(file.path(script.path,"createList.R"))

# set colors
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")

# customized function
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

pca2batch <- function(indata, batch, batchvar = colnames(batch), fig.dir, PCA.fig.title, pos1 = "bottomright", pos2 = "topright", xy=c(1,2), pch.orginal=16, cols=NULL, showID=FALSE, cex=1, showLegend=T, batch1move=0.3, batch2move=0.4, width=5, height=5, withoutgrid=TRUE) {
  
  library(ClassDiscovery)
  
  N.batch1 = length(unique(batch[,batchvar[1]]))    
  N.batch2 = length(unique(batch[,batchvar[2]]))    
  
  if (is.null(cols)) { 
    cols <- rainbow(N.batch1) 
  }else{
    if (length(cols) != N.batch1) {stop("cols length not equal to batch length")} 
  }           
  
  indata=na.omit(indata)
  pca <- SamplePCA(indata, usecor=F, center=F)
  pct1 <- round (pca@variances[xy[1]]/sum(pca@variances), digits=3)*100
  pct2 <- round (pca@variances[xy[2]]/sum(pca@variances), digits=3)*100 
  xlab.text = paste("PC", xy[1], ": ", as.character(pct1), "% variance", sep="")
  ylab.text = paste("PC", xy[2], ": ", as.character(pct2), "% variance", sep="")
  
  if(withoutgrid) {
    outfile = file.path(fig.dir, paste(PCA.fig.title, ".withoutgrid.pdf",sep="")) 
    
    pdf(file=outfile, width = width, height = height)
    par(mar = par()$mar + c(0,0,0,6)) 
    plot(pca@scores[,xy[1]], pca@scores[,xy[2]], 
         cex=cex, xlab=xlab.text, ylab=ylab.text,
         col=cols[factor(batch[,batchvar[1]])], 
         pch=(pch.orginal:(pch.orginal-1+N.batch2))[factor(batch[,batchvar[2]])])
    abline(h=0, v=0, col="brown", lty=2)
  } else {
    outfile = file.path(fig.dir, paste(PCA.fig.title, ".withgrid.pdf",sep=""))
    
    pdf(file=outfile, width = width, height = height)
    par(mar = par()$mar + c(0,0,0,6))
    plot(pca@scores[,xy[1]], pca@scores[,xy[2]], 
         cex=cex, xlab=xlab.text, ylab=ylab.text,
         col=cols[factor(batch[,batchvar[1]])], 
         pch=(pch.orginal:(pch.orginal-1+N.batch2))[factor(batch[,batchvar[2]])])
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray80")
    
    par(new=TRUE)
    plot(pca@scores[,xy[1]], pca@scores[,xy[2]], 
         cex=cex, xlab=xlab.text, ylab=ylab.text,
         col=cols[factor(batch[,batchvar[1]])], 
         pch=(pch.orginal:(pch.orginal-1+N.batch2))[factor(batch[,batchvar[2]])])
    grid (lty = 1, col = "white",lwd=2)
  }
  
  if (showID) { 
    text(pca@scores[,xy[1]], pca@scores[,xy[2]], colnames(indata), lwd=1, cex=cex)
  }
  if(showLegend){ 
    par(xpd = TRUE) #all plotting is clipped to the figure region
    legend(pos1,fill = cols,
           legend=names(table(factor(batch[,batchvar[1]]))),
           inset=c(-batch1move,0),
           border = NA,bty = "n")
    legend(pos2,legend=names(table(factor(batch[,batchvar[2]]))),
           border = NA,bty = "n",
           inset=c(-batch2move,0), 
           pch = (pch.orginal:(pch.orginal-1+N.batch2)))
  }
  invisible(dev.off())
}

#-----------------------------#
# generate meta file for gaby #
meta <- read.csv(file.path(data.path,"Gaby_Wilm_raw data RNAseq.csv"),row.names = 1,check.names = F,stringsAsFactors = F,header = T)
countsTab <- read.table(file.path(data.path,"Wilm_TARGET_Project_rawCounts_V15.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

meta$has_exp <- "No"
meta[colnames(countsTab),"has_exp"] <- "Yes"
write.csv(meta,file.path(res.path,"meta_data_for_Gaby.csv"),row.names = T,quote = F)

#-------------------------#
# RNA expression analysis #
library(pheatmap)
library(dendsort)

# 1. load the entire RNA expression matrix and convert to TPM
HERVKInFile <- NULL
exp.file <- "Wilm_TARGET_Project_rawCounts_V15.txt"
SinfoFile <- "matched_sample_info395.txt"
res <- loadData(data.path, shareFun.path, script.path, comAnn.path, createGinfoFileFlag, SinfoFile, GinfoFile, HERVKInFile, GENEInFile = exp.file, tailrows)
Sinfo=res$Sinfo; Ginfo<-res$Ginfo; Vids<-res$Vids; Lids<-res$Lids; Mids<-res$Mids; countsTable<-res$countsTable

res <-normData(countsTable, tailrows, res.path, tumorname)
countsNorm <- res$countsNorm; vsd <- res$vsd
FPKM <- calFPKM(countsTable, Vids, tailrows, Ginfo, outfile=file.path(res.path, paste(tumorname, ".FPKM.all.features.txt", sep="")))
FPKM[1:3,1:3]
TPM <- as.data.frame(round(apply(FPKM,2,fpkmToTpm),2))

Mids <- intersect(Mids,rownames(TPM))
Lids <- intersect(Lids,rownames(TPM))

# alter Sinfo
Sinfo[which(Sinfo$Subtype == "CDC"),"Subtype"] <- "Bellini"
Sinfo[which(Sinfo$Subtype == "Clear-cell RCC"),"Subtype"] <- "gKIRC"

# 2. clustering analysis by using all samples
rcc.sam <- Sinfo[which(Sinfo$SampleType != "Normal" & Sinfo$to.add.for.analysis.for.WT.project.06042020 %in% c("","yes")),]
rcc.sam <- rownames(rcc.sam[which(rcc.sam$Subtype %in% c("Diffuse anaplastic WT","Focal anaplastic WT","KICH","KIRC","KIRP","gKIRC")),]) # remove bellini cancers
target.sam <- rownames(Sinfo[which(Sinfo$Subtype == "TARGET_Wilm"),])
all.sam <- c(rcc.sam,target.sam)
annCol.all.gaby <- Sinfo[all.sam,1:2]
annCol.all.gaby$Cohort <- ifelse(annCol.all.gaby$Subtype %in% c("BLCA","KICH","KIRC","KIRP","TRCC"),"TCGA",
                                 ifelse(annCol.all.gaby$Subtype == "TARGET_Wilm","TARGET","GABY"))
mycol <- RColorBrewer::brewer.pal(n = 6, name = 'Paired')
annColors.all.gaby <- list("Subtype" = c("gKIRC" = mycol[1],
                                         "KIRC" = mycol[1],
                                         "Diffuse anaplastic WT" = mycol[2],
                                         "Focal anaplastic WT" = mycol[3],
                                         "TARGET_Wilm" = mycol[4],
                                         "KIRP" = mycol[5],
                                         "KICH" = mycol[6]),
                           "Cohort" = c("GABY" = darkblue,
                                        "TARGET" = seagreen,
                                        "TCGA" = yellow))

set.seed(20000112)
target.fhwt.sel <- sample(rownames(annCol.target.fhwt),15)

# generate expression data for all selected tumor samples
all.exp <- TPM[,c(rcc.sam,target.dawt.sel,target.fhwt.sel)]
all.logexp <- log2(all.exp + 1)
indata <- all.logexp[rowSums(all.logexp) > 0,]

# pca to view the batches
batchPCA(indata = indata[,c(rcc.sam,target.dawt.sel,target.fhwt.sel)],
         batch = annCol.all.gaby[c(rcc.sam,target.dawt.sel,target.fhwt.sel),"Cohort"],
         fig.dir = fig.path,
         PCA.fig.title = "PCA.all.gaby3",
         cols = annColors.all.gaby$Cohort,
         showID = F,
         cex = 0.9,
         showLegend = T)

tmp <- annCol.all.gaby
tmp <- rbind.data.frame(tmp[tmp$SampleType == "Tumor",],tmp[c(target.dawt.sel,target.fhwt.sel),])
tmp <- Sinfo[rownames(tmp),]
write.table(tmp,file.path(res.path,"randomly selected 95 kidney tumor3.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# remove batch effect by combat
library(sva)
modcombat = model.matrix(~1, data=annCol.all.gaby[c(rcc.sam,target.dawt.sel,target.fhwt.sel),])
all.logexp.combat = ComBat(dat=as.matrix(indata), batch=annCol.all.gaby[c(rcc.sam,target.dawt.sel,target.fhwt.sel),"Cohort"], mod=modcombat)
batchPCA(indata = all.logexp.combat[,c(rcc.sam,target.dawt.sel,target.fhwt.sel)],
         batch = annCol.all.gaby[c(rcc.sam,target.dawt.sel,target.fhwt.sel),"Cohort"],
         fig.dir = fig.path,
         PCA.fig.title = "PCA.all.gaby.combat3",
         cols = annColors.all.gaby$Cohort,
         showID = F,
         cex = 0.9,
         showLegend = T)

pca2batch(indata = all.logexp.combat[,c(rcc.sam,target.dawt.sel,target.fhwt.sel)],
          batch = annCol.all.gaby[c(rcc.sam,target.dawt.sel,target.fhwt.sel),],
          batchvar = c("Subtype","Cohort"),
          fig.dir = fig.path,
          PCA.fig.title = "PCA.all.gaby.combat.subtype3",
          cols = as.character(annColors.all.gaby$Subtype[sort(names(annColors.all.gaby$Subtype))]),
          showID = F,
          cex = 0.8,
          showLegend = T,
          width=6, height=5, 
          batch1move = 0.6,batch2move = 0,
          pos1="bottomright", pos2="topright",
          withoutgrid = T) 

# kirc sample batches
tmp <- annCol.all.gaby[c(rcc.sam,target.dawt.sel,target.fhwt.sel),]
tmp <- tmp[which(tmp$Subtype %in% c("KIRC","gKIRC")),]
batchPCA(indata = all.logexp.combat[,rownames(tmp)],
         batch = annCol.all.gaby[rownames(tmp),"Subtype"],
         fig.dir = fig.path,
         PCA.fig.title = "PCA.kirc.tumors.combat3",
         cols = c(red,blue),
         showID = T,
         cex = 0.6,
         showLegend = T)

# calculate estimate score
tmp <- as.data.frame(all.logexp.combat[,c(rcc.sam,target.dawt.sel,target.fhwt.sel)])
tmp$symbol <- Ginfo[rownames(tmp),"genename"]
tmp <- apply(tmp[,setdiff(colnames(tmp), "symbol")], 2, function(x) tapply(x, INDEX=factor(tmp$symbol), FUN=median, na.rm=TRUE))
write.table(tmp, file.path(res.path,"all_with_random_30sam_of_TARGET_estimate.input2.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f=file.path(res.path, "all_with_random_30sam_of_TARGET_estimate.input2.txt") , output.f=file.path(res.path,"all_with_random_30sam_of_TARGET_ESTIMATE2.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"all_with_random_30sam_of_TARGET_ESTIMATE2.txt"), file.path(res.path,"all_with_random_30sam_of_TARGET_estimate_score2.txt"), platform="affymetrix")
est <- read.table(file = file.path(res.path,"all_with_random_30sam_of_TARGET_estimate_score2.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est) <- est[,2]; colnames(est) <- est[1,]; est <- est[-1,c(-1,-2)];
est <- sapply(est, as.numeric); rownames(est) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")
colnames(est) <- colnames(indata); est.backup <- est
est <- annTrackScale(indata = est, halfwidth = 2, poolsd = F); est <- as.data.frame(t(est))

annCol.all.gaby <- Sinfo[all.sam,1:2]
annCol.all.gaby$Cohort <- ifelse(annCol.all.gaby$Subtype %in% c("BLCA","KICH","KIRC","KIRP","TRCC"),"TCGA",
                                 ifelse(annCol.all.gaby$Subtype == "TARGET_Wilm","TARGET","GABY"))
annCol.all.gaby <- cbind.data.frame(annCol.all.gaby, est[rownames(annCol.all.gaby),c("ImmuneScore","StromalScore")])
annCol.all.gaby <- cbind.data.frame(annCol.all.gaby, as.data.frame(t(est.backup))[rownames(annCol.all.gaby),"TumorPurity",drop = F])

annColors.all.gaby[["ImmuneScore"]] <- greenred(64)
annColors.all.gaby[["StromalScore"]] <- greenred(64)

# identify optimal cluster number
indata <- all.logexp.combat[,c(rcc.sam,target.dawt.sel,target.fhwt.sel)]
var <- apply(indata, 1, mad)
sel_gene <- var[var > quantile(var,probs = seq(0,1,0.1))[10]]
sel_gene2 <- var[order(var,decreasing = T)][1:500]
indata <- indata[names(sel_gene),]

annCol.random <- annCol.all.gaby[colnames(indata),c(2,3,4,5)]
annCol.random[which(annCol.random$Subtype == "Diffuse anaplastic WT"),"Subtype"] <- "DAWT"
annCol.random[which(annCol.random$Subtype == "Focal anaplastic WT"),"Subtype"] <- "FAWT"
annCol.random[which(annCol.random$Subtype == "gKIRC"),"Subtype"] <- "KIRC"
annCol.random[target.dawt.sel,"Subtype"] <- "DAWT"
annCol.random[target.fhwt.sel,"Subtype"] <- "FHWT"

annColors.random <- list("Subtype" = c(
  "KIRC" = mycol[1],
  "DAWT" = mycol[2],
  "FAWT" = mycol[3],
  "FHWT" = mycol[4],
  "KIRP" = mycol[5],
  "KICH" = mycol[6]),
  "Cohort" = c("GABY" = darkblue,
               "TARGET" = seagreen,
               "TCGA" = yellow))
annColors.random[["ImmuneScore"]] <- greenred(64)
annColors.random[["StromalScore"]] <- greenred(64)

# use consensus hierachical clustering
N.cluster <- 5
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.9*nrow(indata))
N.sample.per.bootstrap <- round(0.9*ncol(indata))
map.res.path <- file.path(tumor.path, paste("commonCluster3_dawt_fhwt_combat_log2TPM_ClusterNum",N.cluster,sep = ""))
featType <- "all_with_random_30sam_of_TARGET3"

cluster.all.ans <- plot.common.cluster(indata, 
                                       tumorname=tumorname, 
                                       N.cluster=N.cluster, 
                                       N.bootstrap=N.bootstrap, 
                                       N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                       N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                       map.res.path=map.res.path, fig.path=fig.path, 
                                       featType=featType,
                                       annCol=annCol.random[colnames(indata),], annColors=annColors.random, 
                                       seed=20000112, 
                                       dist0="euclidean", link0 = "ward.D",
                                       dist="euclidean", link = "ward.D",
                                       clstCol = mycol[1:5], 
                                       height = 7, dendsort = T)

pdf(file.path(fig.path, "consensus heatmap of all tumors with random 30 TARGET3.pdf"), height=8,width = 10)
pheatmap(as.matrix(cluster.all.ans$sum),
         cluster_rows = as.hclust(cluster.all.ans$dendro),
         cluster_cols = as.hclust(cluster.all.ans$dendro),
         border_color = NA,
         annotation_col = cluster.all.ans$annCol,
         annotation_colors = cluster.all.ans$annColors,
         # color = NMF:::ccRamp(c("#2165AA","#F2F1F3","#B91F33"),64),
         color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 4)
invisible(dev.off())

hcg <- hclust(distanceMatrix(as.matrix(t(indata[names(sel_gene2),])), "euclidean"), "ward.D")
plotdata <- standarize.fun(indata[names(sel_gene2),],halfwidth = 2)
use_raster <- F
p <- pheatmap(as.matrix(plotdata),
              cluster_rows = dendsort(hcg),
              cluster_cols = as.hclust(cluster.all.ans$dendro),
              border_color = NA,
              annotation_col = cluster.all.ans$annCol[,c(1,2,5)],
              annotation_colors = cluster.all.ans$annColors,
              color = NMF:::ccRamp(c("#2165AA","#F2F1F3","#B91F33"),64),
              # color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              show_colnames = F,
              show_rownames = F,
              cellwidth = 2,
              cellheight = 200/nrow(plotdata),
              fontsize_col = 4)
pdf(file.path(fig.path, "expression heatmap of all tumors with random 30 TARGET3.pdf"), height=8,width = 10)
draw(p, annotation_legend_side = "left", use_raster = F)
invisible(dev.off())

# calculate immune enrichment score and re-cluster
indata <- read.table(file.path(res.path,"all_with_random_30sam_of_TARGET_estimate.input2.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
indata <- gsva(expr = as.matrix(indata),
               gset.idx.list = immune.sig.ccr,
               method = "gsva")
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D"); hcs.random30 <- hcs
plotdata <- standarize.fun(indata,halfwidth = 2)
p <- pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
              cluster_rows = F,
              cluster_cols = hcs,
              border_color = NA,
              annotation_col = cluster.all.ans$annCol,
              annotation_colors = cluster.all.ans$annColors,
              #color = NMF:::ccRamp(c("#2165AA","#F2F1F3","#B91F33"),64),
              color = NMF:::ccRamp(heatmap.BlBkRd,64),
              #color = inferno(64),
              gaps_row = c(14,22),
              treeheight_row = 15,
              treeheight_col = 15,
              show_colnames = F,
              show_rownames = T,
              cellwidth = 2,
              cellheight = 8,
              fontsize_col = 4,
              fontsize_row = 8)
pdf(file = file.path(fig.path,"CCR immune signature heatmap of all with random 30sam of target in 2 clusters3.pdf"), width = 10,height = 8)
draw(p, annotation_legend_side = "left", use_raster = F)
invisible(dev.off())

table(cluster.all.ans$annCol$Subtype,cluster.all.ans$annCol$Clustering)
# cluster1 cluster2 cluster3 cluster4 cluster5
# DAWT       14        0        8        0        0
# FAWT        0        0        3        0        0
# FHWT       10        1        4        0        0
# KICH        0        0        0        0       15
# KIRC        2       22        0        0        1
# KIRP        0        1        2       12        0
fisher.test(matrix(c(3,0,12,14+10+1), byrow = T,nrow = 2))

# 3. clustering analysis by using only diffuse and focal WT
wt.gaby.sam <- rownames(Sinfo[which(Sinfo$Subtype %in% c("Diffuse anaplastic WT","Focal anaplastic WT") & Sinfo$to.add.for.analysis.for.WT.project.06042020 == "yes"),])
annCol.wt.gaby <- Sinfo[wt.gaby.sam,"Subtype",drop = F]
annColors.wt.gaby <- list("Subtype" = c("Diffuse anaplastic WT" = darkblue,"Focal anaplastic WT" = sun))

# get most variable genes (median absolute deviation)
indata <- log2(TPM[,wt.gaby.sam] + 1)
indata <- indata[rowSums(indata) > 0,]
var <- apply(indata, 1, mad)
sel_gene <- var[var > quantile(var,probs = seq(0,1,0.1))[10]]
indata <- indata[names(sel_gene),]
batchPCA(indata = indata,
         batch = annCol.wt.gaby$Subtype,
         fig.dir = fig.path,
         PCA.fig.title = "PCA.wt.gaby",
         cols = c(darkblue,sun),
         showID = F,
         cex = 0.9,
         showLegend = T)

# use consensus hierachical clustering
N.cluster <- 2
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.9*nrow(indata))
N.sample.per.bootstrap <- round(1*ncol(indata))
map.res.path <- file.path(tumor.path, paste("commonCluster_log2TPM_mad10_ClusterNum",N.cluster,sep = ""))
featType <- "wt.gaby"

cluster2.wt.gaby.ans <- plot.common.cluster(indata, 
                                            tumorname=tumorname, 
                                            N.cluster=N.cluster, 
                                            N.bootstrap=N.bootstrap, 
                                            N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                            N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                            map.res.path=map.res.path, fig.path=fig.path, 
                                            featType=featType,
                                            annCol=annCol.wt.gaby, annColors=annColors.wt.gaby, 
                                            seed=2020415, 
                                            dist0="pearson", link0 = "ward.D",
                                            dist="pearson", link = "ward.D",
                                            clstCol=c(jco[2],jco[1]), 
                                            height = 7, dendsort = T)

pdf(file.path(fig.path, "consensus heatmap of wilms tumor in gaby cohort.pdf"), height=5,width = 6)
pheatmap(as.matrix(cluster2.wt.gaby.ans$sum),
         cluster_rows = as.hclust(cluster2.wt.gaby.ans$dendro),
         cluster_cols = as.hclust(cluster2.wt.gaby.ans$dendro),
         border_color = NA,
         annotation_col = cluster2.wt.gaby.ans$annCol,
         annotation_colors = cluster2.wt.gaby.ans$annColors,
         color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         show_colnames = T,
         show_rownames = F)
invisible(dev.off())

hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward.D2")
pdf(file.path(fig.path, "expression heatmap of wilms tumors in gaby cohort.pdf"), height=4,width = 5.5)
plotdata <- standarize.fun(indata,halfwidth = 2)
pheatmap(as.matrix(plotdata),
         cluster_rows = dendsort(hcg),
         cluster_cols = as.hclust(cluster2.wt.gaby.ans$dendro),
         border_color = NA,
         annotation_col = cluster2.wt.gaby.ans$annCol,
         annotation_colors = cluster2.wt.gaby.ans$annColors,
         color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 8)
invisible(dev.off())

# 4. immunity evaluation by MCPcounter, estimate, CIBERSORT, and immune signatures from Nature

indata <- log2(TPM[Mids,wt.gaby.sam] + 1)
indata <- indata[rowSums(indata) > 0,]
indata$symbol <- Ginfo[rownames(indata),"genename"]
indata <- apply(indata[,setdiff(colnames(indata), "symbol")], 2, function(x) tapply(x, INDEX=factor(indata$symbol), FUN=median, na.rm=TRUE))
write.table(indata,file = file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# estimate
library(estimate)
filterCommonGenes(input.f=file.path(res.path, "Wilms_tumor_Gaby_log2TPM_hugo.txt") , output.f=file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo_ESTIMATE.txt"), file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo_estimate_score.txt"), platform="affymetrix")
est <- read.table(file = file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est) <- est[,2]; colnames(est) <- est[1,]; est <- est[-1,c(-1,-2)];
est <- sapply(est, as.numeric); rownames(est) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.backup = as.data.frame(est); colnames(est.backup) <- colnames(indata)
est <- annTrackScale(indata = est, halfwidth = 2, poolsd = F); est <- as.data.frame(t(est))
rownames(est) <- colnames(indata)

# MCPcounter
library(MCPcounter)
MCPscore <- MCPcounter.estimate(expression = indata,featuresType = "HUGO_symbols")
MCPscore <- as.data.frame(MCPscore); MCPscore.raw.backup <- MCPscore
MCPscore <- annTrackScale(indata = MCPscore, halfwidth = 2, poolsd = F); MCPscore <- as.data.frame(t(MCPscore))

MCPscore.pure <- as.data.frame(MCPscore.raw.backup)
for (sam in colnames(MCPscore.pure)) {
  for (cell in rownames(MCPscore.pure)) {
    # MCPscore.pure[cell,sam] <- MCPscore.raw.backup[cell,sam]/(1-est.backup[4,sam])
    MCPscore.pure[cell,sam] <- MCPscore.raw.backup[cell,sam]*est.backup[4,sam]
  }
}
MCPscore.pure <- annTrackScale(indata = MCPscore.pure, halfwidth = 2, poolsd = F); MCPscore.pure <- as.data.frame(t(MCPscore.pure))

# cibersort
# source(file.path(script.path,"CIBERSORT.R"))
# ciber.res.r <- CIBERSORT(file.path(script.path,"ref.txt"), file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo.txt"), perm=1000, QN=TRUE)
ciber.abs.res <- read.table(file.path(res.path,"CIBERSORT.Output.absolute.Wilms_tumor_Gaby.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
ciber.abs.res <- ciber.abs.res[,1:22]
ciber.abs.res <- ciber.abs.res[,colSums(ciber.abs.res)>0]
ciber.abs.res$m12_ratio <- ciber.abs.res$`Macrophages M1`/ciber.abs.res$`Macrophages M2`
ciber.abs.res <- annTrackScale(indata = t(ciber.abs.res), halfwidth = 2, poolsd = F); ciber.abs.res <- as.data.frame(t(ciber.abs.res))

ciber.rel.res <- read.table(file.path(res.path,"CIBERSORT.Output.relative.Wilms_tumor_Gaby.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
ciber.rel.res <- ciber.rel.res[,1:22]
ciber.rel.res <- ciber.rel.res[,colSums(ciber.rel.res)>0]
ciber.rel.res$m12_ratio <- ciber.rel.res$`Macrophages M1`/ciber.rel.res$`Macrophages M2`
ciber.rel.res <- annTrackScale(indata = t(ciber.rel.res), halfwidth = 2, poolsd = F); ciber.rel.res <- as.data.frame(t(ciber.rel.res))

# nature immune signature
library(GSVA)
indata <- read.table(file.path(res.path, "Wilms_tumor_Gaby_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
immune.signature <- read.table(file.path(comAnn.path,"Nature_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$`Cell Type`)
immune.sig.nature <- list()
for (i in cell.type) {
  immune.sig.nature[[i]] <- intersect(toupper(immune.signature[which(immune.signature$`Cell Type` == i),"Gene"]),rownames(indata))
}
imm.score.gaby.wt <- gsva(as.matrix(indata),immune.sig.nature,method="ssgsea")

# cell immune signature
immune.signature <- read.table(file.path(comAnn.path,"Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.cell <- list()
for (i in cell.type) {
  immune.sig.cell[[i]] <- intersect(toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"]),rownames(indata))
}
immune.sig.cell[["Angiogenesis"]] <- intersect(c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP"),rownames(indata))
immune.sig.cell <- immune.sig.cell[setdiff(names(immune.sig.cell),"pDC")]
imm.score.gaby.wt2 <- gsva(as.matrix(indata),immune.sig.cell,method="ssgsea")

# CCR immune signature
immune.signature <- read.table(file.path(comAnn.path,"CCR_Curated_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- intersect(toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"]),rownames(indata))
}
imm.score.gaby.wt3 <- gsva(as.matrix(indata),immune.sig.ccr,method="gsva",parallel.sz = 1)
imm.score.gaby.wt3.pure <- imm.score.gaby.wt3
for (sam in colnames(imm.score.gaby.wt3)) {
  for (cell in rownames(imm.score.gaby.wt3)) {
    imm.score.gaby.wt3.pure[cell,sam] <- imm.score.gaby.wt3[cell,sam]/(1-est.backup[4,sam])
  }
}
immune.sig.ccr.order <- c("T.cells.CD8",
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")

# me immune signature
indata <- read.table(file = file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
immune.signature <- read.table(file.path(comAnn.path,"Gaby_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.me <- list()
for (i in cell.type) {
  immune.sig.me[[i]] <- intersect(toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"]),rownames(indata))
}

imm.score.gaby.wt4 <- gsva(as.matrix(indata),immune.sig.me,method="gsva",parallel.sz = 1)
imm.score.gaby.wt4.pure <- imm.score.gaby.wt4
for (sam in colnames(imm.score.gaby.wt4)) {
  for (cell in rownames(imm.score.gaby.wt4)) {
    imm.score.gaby.wt4.pure[cell,sam] <- imm.score.gaby.wt4[cell,sam]/(1-est.backup[4,sam])
  }
}

# add immunity to annotation  
annCol.wt.gaby <- cbind.data.frame(annCol.wt.gaby[,1,drop = F],MCPscore[rownames(annCol.wt.gaby),])
annCol.wt.gaby <- cbind.data.frame(annCol.wt.gaby,est[rownames(annCol.wt.gaby),1:2])

annColors.wt.gaby[["T cells"]] <- annColors.wt.gaby[["CD8 T cells"]] <- greenred(64)
annColors.wt.gaby[["Cytotoxic lymphocytes"]] <- annColors.wt.gaby[["NK cells"]] <-annColors.wt.gaby[["B lineage"]] <- greenred(64)
annColors.wt.gaby[["Monocytic lineage"]] <- annColors.wt.gaby[["Myeloid dendritic cells"]] <- greenred(64)
annColors.wt.gaby[["Neutrophils"]] <- annColors.wt.gaby[["Endothelial cells"]] <- annColors.wt.gaby[["Fibroblasts"]] <- greenred(64)
annColors.wt.gaby[["ImmuneScore"]] <- blueyellow(64)
annColors.wt.gaby[["StromalScore"]] <- blueyellow(64)

# 5. use mRNA and lncRNA seperately for gaby WT clustering to see consistency
# mRNA
indata <- log2(TPM[Mids,wt.gaby.sam] + 1)
indata <- indata[rowSums(indata) > 0,]
var <- apply(indata, 1, mad)
sel_gene <- var[var > quantile(var,probs = seq(0,1,0.1))[10]]
indata <- indata[names(sel_gene),]

N.cluster <- 2
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.9*nrow(indata))
N.sample.per.bootstrap <- round(1*ncol(indata))
map.res.path <- file.path(tumor.path, paste("commonCluster_mRNA_log2TPM_mad10_ClusterNum",N.cluster,sep = ""))
featType <- "wt.gaby.mRNA"

cluster2.wt.gaby.mRNA.ans <- plot.common.cluster(indata, 
                                                 tumorname=tumorname, 
                                                 N.cluster=N.cluster, 
                                                 N.bootstrap=N.bootstrap, 
                                                 N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                 N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                 map.res.path=map.res.path, fig.path=fig.path, 
                                                 featType=featType,
                                                 annCol=annCol.wt.gaby, annColors=annColors.wt.gaby, 
                                                 seed=2020415, 
                                                 dist0="pearson", link0 = "ward.D",
                                                 dist="pearson", link = "ward.D",
                                                 clstCol=c(jco[2],jco[1]), 
                                                 height = 7, dendsort = T)

pdf(file.path(fig.path, "consensus heatmap of wilms tumor in gaby cohort using mRNA.pdf"), height=8,width = 10)
pheatmap(as.matrix(cluster2.wt.gaby.mRNA.ans$sum),
         cluster_rows = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = NA,
         annotation_col = annCol.wt.gaby,
         annotation_colors = annColors.wt.gaby,
         # color = NMF:::ccRamp(c("#2165AA","#F2F1F3","#B91F33"),64),
         # color = NMF:::ccRamp(c("#2165AA","#F2F1F3","#B91F33"),64),
         color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 15,
         show_colnames = T,
         show_rownames = F)
invisible(dev.off())

hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward.D2")
pdf(file.path(fig.path, "expression heatmap of wilms tumors in gaby cohort using mRNA.pdf"), height=8,width = 10)
plotdata <- standarize.fun(indata,halfwidth = 2)
pheatmap(as.matrix(plotdata),
         cluster_rows = dendsort(hcg),
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = NA,
         annotation_col = cluster2.wt.gaby.mRNA.ans$annCol,
         annotation_colors = cluster2.wt.gaby.mRNA.ans$annColors,
         # color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
         # color = viridisLite::inferno(64),
         # color = viridisLite::viridis(64),
         color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 0.1,
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 8)
invisible(dev.off())

# 6. heatmap showing cibersort, immune signatures using mRNA dendrogram
plotdata <- t(ciber.abs.res[rownames(cluster2.wt.gaby.mRNA.ans$annCol),])
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D2")
pdf(file.path(fig.path, "cibersort heatmap of wilms tumors in gaby cohort.pdf"), height=8,width = 10)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = NA,
         annotation_col = cluster2.wt.gaby.mRNA.ans$annCol,
         annotation_colors = cluster2.wt.gaby.mRNA.ans$annColors,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

plotdata <- standarize.fun(imm.score.gaby.wt[,rownames(cluster2.wt.gaby.mRNA.ans$annCol)],halfwidth = 1)
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D2")
pdf(file.path(fig.path, "nature immune signature heatmap of wilms tumors in gaby cohort.pdf"), height=8,width = 10)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = NA,
         annotation_col = cluster2.wt.gaby.mRNA.ans$annCol,
         annotation_colors = cluster2.wt.gaby.mRNA.ans$annColors,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

plotdata <- standarize.fun(imm.score.gaby.wt2[,rownames(cluster2.wt.gaby.mRNA.ans$annCol)],halfwidth = 1)
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D2")
pdf(file.path(fig.path, "cell immune signature heatmap of wilms tumors in gaby cohort.pdf"), height=8,width = 10)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = NA,
         annotation_col = cluster2.wt.gaby.mRNA.ans$annCol,
         annotation_colors = cluster2.wt.gaby.mRNA.ans$annColors,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

plotdata <- standarize.fun(imm.score.gaby.wt3.pure[,rownames(cluster2.wt.gaby.mRNA.ans$annCol)],halfwidth = 3)
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D2")
pdf(file.path(fig.path, "CCR immune signature heatmap of wilms tumors in gaby cohort.pdf"), height=8,width = 10)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = NA,
         annotation_col = cluster2.wt.gaby.mRNA.ans$annCol,
         annotation_colors = cluster2.wt.gaby.mRNA.ans$annColors,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         # color = viridisLite::inferno(64),
         # color = viridisLite::viridis(64),
         # color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

# 7. differential expression analysis between two RNA clusters of Wilms tumor
library(edgeR)
source(file.path(shareFun.path,"twoclassedgeR.R"))
annCol.wt.gaby$RNA_Clust <- cluster2.wt.gaby.ans$group[rownames(annCol.wt.gaby)]
annCol.wt.gaby$RNA_Clust <- ifelse(annCol.wt.gaby$RNA_Clust == "cluster1","TP53_Wild","TP53_Mutated")
annColors.wt.gaby[["RNA_Clust"]] <- c("TP53_Wild" = jco[2],"TP53_Mutated" = jco[1])

# mRNA
tmp <- annCol.wt.gaby$RNA_Clust; names(tmp) <- rownames(annCol.wt.gaby)
Groupinfo <- annCol.wt.gaby
colnames(Groupinfo)[14] <- "group"
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.RNAclust(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="Gaby_WT_mRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

# lncRNA
tmp <- annCol.wt.gaby$RNA_Clust; names(tmp) <- rownames(annCol.wt.gaby)
Groupinfo <- annCol.wt.gaby
colnames(Groupinfo)[14] <- "group"
PASSFlag <- rep(TRUE,length(Lids)); names(PASSFlag) <- Lids
complist <- createList.RNAclust(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Lids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Lids, featType="Gaby_WT_lncRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

# 8. heatmap using differentially expressed genes
fdr.cutoff <- 0.05
logfc.cutoff <- 2
deres <- read.table(file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.TP53_Mutated_vs_TP53_Wild.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
up.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC > logfc.cutoff),]
dn.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC < -logfc.cutoff),]
pdf(file.path(fig.path, "DEGs heatmap of wilms tumors in gaby cohort using mRNA.pdf"), height=8,width = 10)
indata <- log2(TPM[c(rownames(up.gene),rownames(dn.gene)),wt.gaby.sam] + 1)
plotdata <- standarize.fun(indata,halfwidth = 2)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.lncRNA.ans$dendro),
         border_color = NA,
         annotation_col = annCol.wt.gaby,
         annotation_colors = annColors.wt.gaby,
         # color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
         # color = viridisLite::inferno(64),
         # color = viridisLite::viridis(64),
         color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 0.1,
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 8)
invisible(dev.off())

deres <- read.table(file.path(res.path,"Gaby_WT_lncRNA_edgeR_test_result.TP53_Mutated_vs_TP53_Wild.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
up.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC > logfc.cutoff),]
dn.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC < -logfc.cutoff),]
pdf(file.path(fig.path, "DEGs heatmap of wilms tumors in gaby cohort using lncRNA.pdf"), height=8,width = 10)
indata <- log2(TPM[c(rownames(up.gene),rownames(dn.gene)),wt.gaby.sam] + 1)
plotdata <- standarize.fun(indata,halfwidth = 2)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.lncRNA.ans$dendro),
         border_color = NA,
         annotation_col = annCol.wt.gaby,
         annotation_colors = annColors.wt.gaby,
         # color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
         # color = viridisLite::inferno(64),
         # color = viridisLite::viridis(64),
         color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 18,
         cellheight = 1,
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 8)
invisible(dev.off())

# 9. GSEA for mRNA
library(clusterProfiler)
library(enrichplot)
MSigDB=read.gmt(file.path(comAnn.path,"msigdb.v7.1.symbols.gmt"))
MSigDB.hallmark=read.gmt(file.path(comAnn.path,"h.all.v7.1.symbols.gmt"))
MSigDB.oncogenic=read.gmt(file.path(comAnn.path,"c6.all.v7.1.symbols.gmt"))
MSigDB.hm.rt <- rbind.data.frame(MSigDB.hallmark,MSigDB[grep("REACTOME_",MSigDB$ont),])

deres <- read.table(file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.TP53_Mutated_vs_TP53_Wild.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- deres$logFC
names(geneList) <- deres$id
geneList <- sort(geneList,decreasing = T)
gsea_gaby_wt_tp53_mutvswild <- GSEA(geneList = geneList,
                                    TERM2GENE=MSigDB,
                                    pvalueCutoff = 0.25,
                                    nPerm = 10000,
                                    seed = T,
                                    verbose=F)
res <- data.frame(gsea_gaby_wt_tp53_mutvswild)
write.table(as.data.frame(gsea_gaby_wt_tp53_mutvswild),file.path(res.path,"GSEA results for Gaby WT mRNA TP53 MutatedvsWild DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

gsea_gaby_wt_tp53_mutvswild_hm.rt <- GSEA(geneList = geneList,
                                    TERM2GENE=MSigDB.hm.rt,
                                    pvalueCutoff = 0.25,
                                    nPerm = 10000,
                                    seed = T,
                                    verbose=F)
res <- data.frame(gsea_gaby_wt_tp53_mutvswild_hm.rt)
write.table(as.data.frame(gsea_gaby_wt_tp53_mutvswild_hm.rt),file.path(res.path,"GSEA hallmark and reactome results for Gaby WT mRNA TP53 MutatedvsWild DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

gseaplot2(x = gsea_gaby_wt_tp53_mutvswild_hm.rt,
          geneSetID = c("HALLMARK_P53_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
          pvalue_table = F,
          color = c(alpha(darkblue,0.3),alpha(darkblue,0.6),darkblue))
ggsave(filename = file.path(fig.path,"GSEA downregulated pathway in gaby cohort.pdf"), width = 5.5,height = 4)

gseaplot2(x = gsea_gaby_wt_tp53_mutvswild_hm.rt,
          #geneSetID = c("REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA","REACTOME_HDACS_DEACETYLATE_HISTONES","REACTOME_HOMOLOGY_DIRECTED_REPAIR"),
          geneSetID = c("REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_CHROMATIN_MODIFYING_ENZYMES","REACTOME_HDACS_DEACETYLATE_HISTONES"),
          pvalue_table = F,
          color = c(alpha(darkred,0.3),alpha(darkred,0.6),darkred))
ggsave(filename = file.path(fig.path,"GSEA upregulated pathway in gaby cohort.pdf"), width = 5.5,height = 4)

gsea_hallmark_gaby_wt_tp53_mutvswild <- GSEA(geneList = geneList,
                                        TERM2GENE=MSigDB.hallmark,
                                        pvalueCutoff = 0.25,
                                        nPerm = 10000,
                                        seed = T,
                                        verbose=F)
res <- data.frame(gsea_hallmark_gaby_wt_tp53_mutvswild)
write.table(as.data.frame(gsea_hallmark_gaby_wt_tp53_mutvswild),file.path(res.path,"GSEA hallmark results for Gaby WT mRNA TP53 MutatedvsWild DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

gseaplot2(x = gsea_hallmark_gaby_wt_tp53_mutvswild,
          geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          pvalue_table = F,
          color = seagreen)
ggsave(filename = file.path(fig.path,"GSEA inteferon gamma.pdf"), width = 5.5,height = 4)

# gseaplot(x = gsea_hallmark_gaby_wt_tp53_mutvswild,
#           geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#           pvalue_table = T)

gsea_oncogenic_gaby_wt_tp53_mutvswild <- GSEA(geneList = geneList,
                                              TERM2GENE=MSigDB.oncogenic,
                                              pvalueCutoff = 0.25,
                                              nPerm = 10000,
                                              seed = T,
                                              verbose=F)
res <- data.frame(gsea_oncogenic_gaby_wt_tp53_mutvswild)
write.table(as.data.frame(gsea_oncogenic_gaby_wt_tp53_mutvswild),file.path(res.path,"GSEA oncogenic results for Gaby WT mRNA TP53 MutatedvsWild DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

# 10. validate TP53 signatures in TARGET cohort
deres <- read.table(file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
up.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC > logfc.cutoff),]
dn.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC < -logfc.cutoff),]

# generate template
n.top <- 300
C1.signature <- up.gene[order(up.gene$logFC,decreasing = T),][1:n.top,] # tp53 wild and immune infiltrated
C2.signature <- dn.gene[order(dn.gene$logFC,decreasing = F),][1:n.top,] # tp53 mutated and immune deplete
tp53.templates <- data.frame(probe=c(rownames(C1.signature),rownames(C2.signature)),Description="na",class=rep(c("immune_infiltrated","immune_deplete"),c(n.top,n.top)))

Sinfo.target <- read.table(file.path(data.path,"curate annotation of TARGET tumors with both expression and mutation add clinical info.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
Sinfo.target <- Sinfo.target[which(Sinfo.target$TP53_MUT != "N/A"),]
Sinfo.target$TP53_ALTER <- "No"
Sinfo.target$TP53_ALTER <- ifelse(Sinfo.target$TP53_MUT == "Yes" | Sinfo.target$TP53_LOSS == "Yes","Yes","No")

TPM.TARGET <- TPM[,rownames(Sinfo.target)]
annCol.target <- data.frame(TP53_MUT = Sinfo.target$TP53_MUT,
                            TP53_LOSS = Sinfo.target$TP53_LOSS,
                            TP53_ALTER = Sinfo.target$TP53_ALTER,
                            HISTOLOGY = Sinfo.target$`Histologic Classification of Primary Tumor`,
                            row.names = rownames(Sinfo.target),
                            stringsAsFactors = F)
annColors.target <- list()
annColors.target[["TP53_MUT"]] <- c("Yes" = "black","No" = "white")
annColors.target[["TP53_LOSS"]] <- c("Yes" = "black","No" = "white")
annColors.target[["TP53_ALTER"]] <- c("Yes" = "black","No" = "white")
annColors.target[["HISTOLOGY"]] <- c("DAWT" = darkblue,"FHWT" = sun)

# nearest template prediction
library(CMScaller)
indata <- log2(TPM.TARGET + 1)
indata <- indata[rowSums(indata)>0,]
indata <- t(scale(t(indata),center = T,scale = T))
set.seed(2020415)
ntp.target <- ntp(indata, tp53.templates, nPerm=1000,distance = "cosine")
ntp.target$pred <- ifelse(ntp.target$FDR < 0.05,as.character(ntp.target$prediction),NA)

annCol.target$NTP_PRED <- ntp.target[rownames(annCol.target),"prediction"]
table(annCol.target$NTP_PRED,annCol.target$TP53_MUT)
table(annCol.target$NTP_PRED,annCol.target$TP53_LOSS)
table(annCol.target$NTP_PRED,annCol.target$TP53_ALTER)
fisher.test(table(annCol.target$NTP_PRED,annCol.target$TP53_ALTER))

# supervised clustering
deres <- read.table(file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
up.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC > logfc.cutoff),]
dn.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC < -logfc.cutoff),]

indata <- log2(TPM.TARGET[c(rownames(up.gene),rownames(dn.gene)),] + 1)
target.mRNA.nmf <- nmf(indata,2,seed = 2020415,method = "nsNMF",nrun = 30)
group <- predict(target.mRNA.nmf)
tmp <- data.frame(NMF = group,
                  TP53_MUT = annCol.target[names(group),"TP53_MUT"],
                  TP53_ALTER = annCol.target[names(group),"TP53_ALTER"])
fisher.test(table(tmp$NMF,tmp$TP53_MUT))

N.cluster <- 2
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.9*nrow(indata))
N.sample.per.bootstrap <- round(0.9*ncol(indata))
map.res.path <- file.path(tumor.path, paste("commonCluster_TARGET_mRNA_log2TPM_ClusterNum",N.cluster,sep = ""))
featType <- "TARGET.mRNA"

cluster2.wt.target.mRNA.ans <- plot.common.cluster(indata, 
                                                   tumorname=tumorname, 
                                                   N.cluster=N.cluster, 
                                                   N.bootstrap=N.bootstrap, 
                                                   N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                   N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                   map.res.path=map.res.path, fig.path=fig.path, 
                                                   featType=featType,
                                                   annCol=annCol.target, annColors=annColors.target, 
                                                   seed=2020415, 
                                                   dist0="euclidean", link0 = "ward.D",
                                                   dist="euclidean", link = "ward.D",
                                                   clstCol=c(jco[2],jco[1]), 
                                                   height = 7, dendsort = T)

deres <- read.table(file.path(res.path,"Gaby_WT_lncRNA_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
up.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC > logfc.cutoff),]
dn.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC < -logfc.cutoff),]

indata <- log2(TPM.TARGET[c(rownames(up.gene),rownames(dn.gene)),] + 1)
N.cluster <- 2
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.9*nrow(indata))
N.sample.per.bootstrap <- round(0.9*ncol(indata))
map.res.path <- file.path(tumor.path, paste("commonCluster_TARGET_lncRNA_log2TPM_ClusterNum",N.cluster,sep = ""))
featType <- "TARGET.lncRNA"

cluster2.wt.target.lncRNA.ans <- plot.common.cluster(indata, 
                                                     tumorname=tumorname, 
                                                     N.cluster=N.cluster, 
                                                     N.bootstrap=N.bootstrap, 
                                                     N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                     N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                     map.res.path=map.res.path, fig.path=fig.path, 
                                                     featType=featType,
                                                     annCol=annCol.target, annColors=annColors.target, 
                                                     seed=2020415, 
                                                     dist0="euclidean", link0 = "ward.D",
                                                     dist="euclidean", link = "ward.D",
                                                     clstCol=c(jco[2],jco[1]), 
                                                     height = 7, dendsort = T)

# 11. differential expression analysis in TARGET cohort regarding TP53 mutation
# mRNA
tmp <- annCol.target$TP53_ALTER; names(tmp) <- rownames(annCol.target)
Groupinfo <- annCol.target[,c("TP53_MUT","TP53_LOSS","TP53_ALTER")]
colnames(Groupinfo)[3] <- "group"
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.TP53(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="TARGET_WT_mRNA_TP53", complist,PASSFlag=PASSFlag, overwt=TRUE)

# lncRNA
tmp <- annCol.target$TP53_ALTER; names(tmp) <- rownames(annCol.target)
Groupinfo <- annCol.target[,c("TP53_MUT","TP53_LOSS","TP53_ALTER")]
colnames(Groupinfo)[3] <- "group"
PASSFlag <- rep(TRUE,length(Lids)); names(PASSFlag) <- Lids
complist <- createList.TP53(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Lids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Lids, featType="TARGET_WT_lncRNA_TP53", complist,PASSFlag=PASSFlag, overwt=TRUE)

# GSEA
deres <- read.table(file.path(res.path,"TARGET_WT_mRNA_TP53_edgeR_test_result.Yes_vs_No.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- deres$logFC
names(geneList) <- deres$id
geneList <- sort(geneList,decreasing = T)
gsea_target_tp53_mutvswild <- GSEA(geneList = geneList,
                              TERM2GENE=MSigDB,
                              pvalueCutoff = 0.25,
                              nPerm = 10000,
                              seed = T,
                              verbose=F)
res <- data.frame(gsea_target_tp53_mutvswild)
write.table(as.data.frame(gsea_target_tp53_mutvswild),file.path(res.path,"GSEA results for TARGET WT mRNA TP53 Mutated vs Wild DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

# 12. immunity in TARGET cohort regarding TP53
indata <- log2(TPM[Mids,rownames(annCol.target)] + 1)
indata <- indata[rowSums(indata) > 0,]
indata$symbol <- Ginfo[rownames(indata),"genename"]
indata <- apply(indata[,setdiff(colnames(indata), "symbol")], 2, function(x) tapply(x, INDEX=factor(indata$symbol), FUN=median, na.rm=TRUE))
write.table(indata,file = file.path(res.path,"Wilms_tumor_TARGET_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

indata <- read.table(file.path(res.path,"Wilms_tumor_TARGET_log2TPM_hugo.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
filterCommonGenes(input.f=file.path(res.path, "Wilms_tumor_TARGET_log2TPM_hugo.txt") , output.f=file.path(res.path,"Wilms_tumor_target_log2TPM_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"Wilms_tumor_target_log2TPM_hugo_ESTIMATE.txt"), file.path(res.path,"Wilms_tumor_target_log2TPM_hugo_estimate_score.txt"), platform="affymetrix")
est.target <- read.table(file = file.path(res.path,"Wilms_tumor_target_log2TPM_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.target) <- est.target[,2]; colnames(est.target) <- est.target[1,]; est.target <- est.target[-1,c(-1,-2)];
est.target <- sapply(est.target, as.numeric); rownames(est.target) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.target.backup = as.data.frame(est.target); colnames(est.target.backup) <- colnames(indata)
est.target <- annTrackScale(indata = est.target, halfwidth = 2, poolsd = F); est.target <- as.data.frame(t(est.target))
rownames(est.target) <- colnames(indata)

# MCPcounter
MCPscore.target <- MCPcounter.estimate(expression = indata,featuresType = "HUGO_symbols")
MCPscore.target <- as.data.frame(MCPscore.target); MCPscore.target.raw.backup <- MCPscore.target
MCPscore.target <- annTrackScale(indata = MCPscore.target, halfwidth = 2, poolsd = F); MCPscore.target <- as.data.frame(t(MCPscore.target))

# nature immune signature
indata <- read.table(file.path(res.path, "Wilms_tumor_TARGET_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
immune.signature <- read.table(file.path(comAnn.path,"Nature_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$`Cell Type`)
immune.sig.nature <- list()
for (i in cell.type) {
  immune.sig.nature[[i]] <- intersect(toupper(immune.signature[which(immune.signature$`Cell Type` == i),"Gene"]),rownames(indata))
}
imm.score.target.wt <- gsva(as.matrix(indata),immune.sig.nature,method="ssgsea")

# cell immune signature
immune.signature <- read.table(file.path(comAnn.path,"Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.cell <- list()
for (i in cell.type) {
  immune.sig.cell[[i]] <- intersect(toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"]),rownames(indata))
}
immune.sig.cell[["Angiogenesis"]] <- intersect(c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP"),rownames(indata))
immune.sig.cell <- immune.sig.cell[setdiff(names(immune.sig.cell),"pDC")]
imm.score.target.wt2 <- gsva(as.matrix(indata),immune.sig.cell,method="ssgsea")

# CCR immune signature
immune.signature <- read.table(file.path(comAnn.path,"CCR_Curated_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- intersect(toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"]),rownames(indata))
}
imm.score.target.wt3 <- gsva(as.matrix(indata),immune.sig.ccr,method="gsva",parallel.sz = 1)

# heatmap
pdf(file.path(fig.path,"MCPcounter heatmap of TARGET data regarding TP53 status.pdf"),width = 6,height = 3)
plotdata <- t(MCPscore.target[,1:8])
hcs <- hclust(distanceMatrix(as.matrix(plotdata), "euclidean"), "ward.D2")
group <- cutree(hcs,2)
annCol.target$MCPcounter_Clust <- paste0("C",group[rownames(annCol.target)])
annColors.target[["MCPcounter_Clust"]] <- c("C1" = jco[2],"C2" = jco[1])
fisher.test(table(annCol.target$TP53_ALTER,annCol.target$MCPcounter_Clust))
pheatmap(plotdata,
         cluster_rows = F,
         cluster_cols = hcs,
         annotation_col = annCol.target,
         annotation_colors = annColors.target,
         show_rownames = T,
         show_colnames = F,
         color = greenred(64))
invisible(dev.off())

plotdata <- standarize.fun(imm.score.target.wt3[,rownames(annCol.target)],halfwidth = 1)
hcs <- hclust(distanceMatrix(as.matrix(plotdata), "euclidean"), "ward.D2")
pdf(file.path(fig.path,"CCR immune signature heatmap of TARGET data regarding TP53 status.pdf"),width = 6,height = 3)
pheatmap(plotdata,
         cluster_rows = F,
         cluster_cols = hcs,
         annotation_col = annCol.target,
         annotation_colors = annColors.target,
         show_rownames = T,
         show_colnames = F,
         color = greenred(64))
invisible(dev.off())

# 13. differential expression analysis between Gaby WT and other kidney tumors (consider batches)
library(edgeR)
Groupinfo <- annCol.all.gaby[which(annCol.all.gaby$Subtype != "TARGET_Wilm"),]
Groupinfo$Subtype <- gsub("gKIRC","KIRC",Groupinfo$Subtype)
Groupinfo[grepl("WT",Groupinfo$Subtype),"Subtype"] <- "WT"
Groupinfo[!grepl("WT",Groupinfo$Subtype),"Subtype"] <- "OtherRCC"
colnames(Groupinfo) <- c("Tissue","group","batch")
tmp <- Groupinfo$group; names(tmp) <- rownames(Groupinfo)

# total RNA
PASSFlag <- rep(TRUE,length(Mids) + length(Lids)); names(PASSFlag) <- c(Mids,Lids)
complist <- createList.WT(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[c(Mids,Lids), names(tmp)], tailrows, Groupinfo=Groupinfo, features=c(Mids,Lids), featType="Gaby_WT_totalRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

# mRNA
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.WT(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="Gaby_WT_mRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

# lncRNA
PASSFlag <- rep(TRUE,length(Lids)); names(PASSFlag) <- Lids
complist <- createList.WT(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Lids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Lids, featType="Gaby_WT_lncRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

# GSEA
deres <- read.table(file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.WT_vs_OtherRCC.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- deres$logFC
names(geneList) <- deres$id
geneList <- sort(geneList,decreasing = T)
gsea_gabyWTvsOtherRCC <- GSEA(geneList = geneList,
                              TERM2GENE=MSigDB,
                              pvalueCutoff = 0.25,
                              nPerm = 10000,
                              seed = T,
                              verbose=F)
res <- data.frame(gsea_gabyWTvsOtherRCC)
write.table(as.data.frame(gsea_gabyWTvsOtherRCC),file.path(res.path,"GSEA results for Gaby WT vs Other RCCs DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

# heatmap of DEGs
deres <- read.table(file.path(res.path,"Gaby_WT_totalRNA_edgeR_test_result.WT_vs_OtherRCC.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
up.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC > logfc.cutoff),]
dn.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC < -logfc.cutoff),]
Groupinfo <- annCol.all.gaby[which(annCol.all.gaby$Subtype != "TARGET_Wilm"),]
Groupinfo$Subtype <- gsub("gKIRC","KIRC",Groupinfo$Subtype)
Groupinfo[grepl("WT",Groupinfo$Subtype),"Subtype"] <- "WT"

pdf(file.path(fig.path, "DEGs heatmap of wilms tumors in gaby cohort vs Other RCCs using totalRNA.pdf"), height=6,width = 7)
indata <- all.logexp.combat[c(rownames(up.gene),rownames(dn.gene)),rownames(Groupinfo)]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(plotdata), "euclidean"), "ward.D2")
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = hcs,
         border_color = NA,
         annotation_col = Groupinfo[,2,drop = F],
         annotation_colors = list(Subtype = c("WT" = "#33A02C",
                                              "KIRC" = "#1F78B4",
                                              "Bellini" = "#A6CEE3",
                                              "KIRP" = "#FB9A99",
                                              "KICH" = "#E31A1C")),
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         # color = viridisLite::inferno(64),
         # color = viridisLite::viridis(64),
         # color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 5,
         cellheight = 0.09,
         show_colnames = F,
         show_rownames = F,
         fontsize_col = 8)
invisible(dev.off())

# 14. comparsion of each element of immunity regarding TP53 alteration
# Gaby cohort
RNA_Clust1 <- rownames(annCol.wt.gaby[which(annCol.wt.gaby$RNA_Clust == "C1"),]) # TP53 wild and immune infiltrated
RNA_Clust2 <- rownames(annCol.wt.gaby[which(annCol.wt.gaby$RNA_Clust == "C2"),]) # TP53 mutated and immune deplete

mcp.comp.wt.gaby <- NULL
for (i in rownames(MCPscore.raw.backup)) {
  tmp1 <- as.numeric(MCPscore.raw.backup[i,RNA_Clust1])
  tmp2 <- as.numeric(MCPscore.raw.backup[i,RNA_Clust2])
  wilcox.p <- wilcox.test(tmp1,tmp2)$p.value
  ttes.p   <- t.test(tmp1,tmp2)$p.value
  mcp.comp.wt.gaby <- rbind.data.frame(mcp.comp.wt.gaby,
                                       data.frame(cell = i,
                                                  avg.tp53.wild = mean(tmp1),
                                                  avg.tp53.mut  = mean(tmp2),
                                                  wilcox.p = round(wilcox.p,4),
                                                  ttest.p  = round(ttes.p,4),
                                                  stringsAsFactors = F),
                                       stringsAsFactors = F)
}
write.table(mcp.comp.wt.gaby,file.path(res.path,"comparision of MCPcounter cells in gaby cohort regarding TP53 status.txt"),sep = "\t",row.names = F,quote = F)

nature.imm.comp.wt.gaby <- NULL
tmp <- standarize.fun(imm.score.gaby.wt,halfwidth = 1)
for (i in rownames(tmp)) {
  tmp1 <- as.numeric(tmp[i,RNA_Clust1])
  tmp2 <- as.numeric(tmp[i,RNA_Clust2])
  wilcox.p <- wilcox.test(tmp1,tmp2)$p.value
  ttes.p   <- t.test(tmp1,tmp2)$p.value
  nature.imm.comp.wt.gaby <- rbind.data.frame(nature.imm.comp.wt.gaby,
                                       data.frame(cell = i,
                                                  avg.tp53.wild = mean(tmp1),
                                                  avg.tp53.mut  = mean(tmp2),
                                                  wilcox.p = round(wilcox.p,4),
                                                  ttest.p  = round(ttes.p,4),
                                                  stringsAsFactors = F),
                                       stringsAsFactors = F)
}
write.table(nature.imm.comp.wt.gaby,file.path(res.path,"comparision of nature immune cells in gaby cohort regarding TP53 status.txt"),sep = "\t",row.names = F,quote = F)

cell.imm.comp.wt.gaby <- NULL
tmp <- standarize.fun(imm.score.gaby.wt2,halfwidth = 1)
for (i in rownames(tmp)) {
  tmp1 <- as.numeric(tmp[i,RNA_Clust1])
  tmp2 <- as.numeric(tmp[i,RNA_Clust2])
  wilcox.p <- wilcox.test(tmp1,tmp2)$p.value
  ttes.p   <- t.test(tmp1,tmp2)$p.value
  cell.imm.comp.wt.gaby <- rbind.data.frame(cell.imm.comp.wt.gaby,
                                              data.frame(cell = i,
                                                         avg.tp53.wild = mean(tmp1),
                                                         avg.tp53.mut  = mean(tmp2),
                                                         wilcox.p = round(wilcox.p,4),
                                                         ttest.p  = round(ttes.p,4),
                                                         stringsAsFactors = F),
                                              stringsAsFactors = F)
}
write.table(cell.imm.comp.wt.gaby,file.path(res.path,"comparision of cell immune cells in gaby cohort regarding TP53 status.txt"),sep = "\t",row.names = F,quote = F)

# target cohort
# annCol.target.mod <- annCol.target[which(annCol.target$HISTOLOGY == "DAWT"),]
tp53_wild <- rownames(annCol.target[which(annCol.target$TP53_ALTER == "No"),])
tp53_altered <- rownames(annCol.target[which(annCol.target$TP53_ALTER == "Yes"),])

mcp.comp.wt.target <- NULL
for (i in rownames(MCPscore.target.raw.backup)) {
  tmp1 <- as.numeric(MCPscore.target.raw.backup[i,tp53_wild])
  tmp2 <- as.numeric(MCPscore.target.raw.backup[i,tp53_altered])
  wilcox.p <- wilcox.test(tmp1,tmp2)$p.value
  ttes.p   <- t.test(tmp1,tmp2)$p.value
  mcp.comp.wt.target <- rbind.data.frame(mcp.comp.wt.target,
                                         data.frame(cell = i,
                                                    avg.tp53.wild = mean(tmp1),
                                                    avg.tp53.mut  = mean(tmp2),
                                                    wilcox.p = round(wilcox.p,4),
                                                    ttest.p  = round(ttes.p,4),
                                                    stringsAsFactors = F),
                                         stringsAsFactors = F)
}
write.table(mcp.comp.wt.target,file.path(res.path,"comparision of MCPcounter cells in TARGET cohort regarding TP53 status.txt"),sep = "\t",row.names = F,quote = F)

nature.imm.comp.wt.target <- NULL
tmp <- standarize.fun(imm.score.target.wt,halfwidth = 1)
for (i in rownames(tmp)) {
  tmp1 <- as.numeric(tmp[i,tp53_wild])
  tmp2 <- as.numeric(tmp[i,tp53_altered])
  wilcox.p <- wilcox.test(tmp1,tmp2)$p.value
  ttes.p   <- t.test(tmp1,tmp2)$p.value
  nature.imm.comp.wt.target <- rbind.data.frame(nature.imm.comp.wt.target,
                                              data.frame(cell = i,
                                                         avg.tp53.wild = mean(tmp1),
                                                         avg.tp53.mut  = mean(tmp2),
                                                         wilcox.p = round(wilcox.p,4),
                                                         ttest.p  = round(ttes.p,4),
                                                         stringsAsFactors = F),
                                              stringsAsFactors = F)
}
write.table(nature.imm.comp.wt.target,file.path(res.path,"comparision of nature immune cells in target cohort regarding TP53 status.txt"),sep = "\t",row.names = F,quote = F)

cell.imm.comp.wt.target <- NULL
tmp <- standarize.fun(imm.score.target.wt2,halfwidth = 1)
for (i in rownames(tmp)) {
  tmp1 <- as.numeric(tmp[i,tp53_wild])
  tmp2 <- as.numeric(tmp[i,tp53_altered])
  wilcox.p <- wilcox.test(tmp1,tmp2)$p.value
  ttes.p   <- t.test(tmp1,tmp2)$p.value
  cell.imm.comp.wt.target <- rbind.data.frame(cell.imm.comp.wt.target,
                                            data.frame(cell = i,
                                                       avg.tp53.wild = mean(tmp1),
                                                       avg.tp53.mut  = mean(tmp2),
                                                       wilcox.p = round(wilcox.p,4),
                                                       ttest.p  = round(ttes.p,4),
                                                       stringsAsFactors = F),
                                            stringsAsFactors = F)
}
write.table(cell.imm.comp.wt.target,file.path(res.path,"comparision of cell immune cells in target cohort regarding TP53 status.txt"),sep = "\t",row.names = F,quote = F)

pdf(file.path(fig.path,"MCPcounter heatmap of TARGET data regarding TP53 status without clustering.pdf"),width = 6,height = 3)
pheatmap(t(MCPscore.target[c(tp53_wild,tp53_altered),]),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = annCol.target[c(tp53_wild,tp53_altered),],
         annotation_colors = annColors.target,
         show_rownames = T,
         show_colnames = F,
         color = greenred(64))
invisible(dev.off())

plotdata <- standarize.fun(imm.score.target.wt3[,rownames(annCol.target)],halfwidth = 1)
pdf(file.path(fig.path,"CCR immune signature heatmap of TARGET data regarding TP53 status without clustering.pdf"),width = 6,height = 3)
hcs <- hclust(distanceMatrix(as.matrix(plotdata), "pearson"), "average")
group <- cutree(hcs,2)
annCol.target$ccr_Clust <- paste0("C",group[rownames(annCol.target)])
fisher.test(table(annCol.target$ccr_Clust,annCol.target$TP53_ALTER))
annColors.target[["ccr_Clust"]] <- c("C1" = jco[1],"C2" = jco[2])
pheatmap(plotdata,
         cluster_rows = F,
         cluster_cols = dendsort(hcs),
         annotation_col = annCol.target[,c(1:4,6)],
         annotation_colors = annColors.target,
         show_rownames = T,
         show_colnames = F,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64))
         #color = greenred(64))
invisible(dev.off())

library(survival)
library(survminer)
library(ggplot2)
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  ccr_Clust = annCol.target[rownames(Sinfo.target),"ccr_Clust"],
                  stringsAsFactors = F)

fitd <- survdiff(Surv(OS.time, OS) ~ ccr_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ ccr_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jama",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

# 15. comparsion of expression of TMEM173
indata <- read.table(file.path(res.path, "Wilms_tumor_Gaby_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
tmp1 <- as.numeric(indata["TMEM173",RNA_Clust1])
tmp2 <- as.numeric(indata["TMEM173",RNA_Clust2])
wilcox.test(tmp1,tmp2) #0.01587
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Wild","Mutated"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("Wild","Mutated"))
pdf(file.path(fig.path,"boxplot for TMEM173 expression in gaby WT cohort regarding TP53 status.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "TP53 status",
        ylab = "log2(TPM) in Gaby",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,3.5,"P = 0.016",cex = 1)
invisible(dev.off())

indata <- read.table(file.path(res.path, "Wilms_tumor_TARGET_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
tmp1 <- as.numeric(indata["TMEM173",tp53_wild])
tmp2 <- as.numeric(indata["TMEM173",tp53_altered])
wilcox.test(tmp1,tmp2) #0.0111
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Wild","Mutated"),c(length(tp53_wild),length(tp53_altered))))
tmp$mut <- factor(tmp$mut,levels = c("Wild","Mutated"))
pdf(file.path(fig.path,"boxplot for TMEM173 expression in TARGET WT cohort regarding TP53 status.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "TP53 status",
        ylab = "log2(TPM) in TARGET",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,6,"P = 0.011",cex = 1)
invisible(dev.off())

tmp1 <- as.numeric(indata["TMEM173",RNA_Clust1.target])
tmp2 <- as.numeric(indata["TMEM173",RNA_Clust2.target])
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("C1","C2"),c(length(RNA_Clust1.target),length(RNA_Clust2.target))))
pdf(file.path(fig.path,"boxplot for TMEM173 expression in TARGET WT cohort regarding 2 clusters.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunologic Cluster",
        ylab = "log2(TPM) in TARGET",
        outline = F,
        col = ggplot2::alpha(c(jco[2:1],seagreen),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1],seagreen))
text(1.5,6.5,"P < 0.001",cex = 1)
invisible(dev.off())

# TMEM173 is not associated with target survival
indata <- read.table(file.path(res.path, "Wilms_tumor_TARGET_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
tmp <- data.frame(exp = as.numeric(indata["TMEM173",]/as.numeric(est.target.backup[4,])),
                  OS.time = Sinfo.target[colnames(indata),"Overall Survival Time in Days"],
                  OS = ifelse(Sinfo.target[colnames(indata),"Vital Status"] == "Dead" , 1, 0),
                  stringsAsFactors = F)
coxph(Surv(OS.time,OS)~exp,data = tmp)
fit <- coxph(Surv(OS.time, OS)~rcspline.eval(tmp$exp,nk=4,inclx = T),data=tmp, x=TRUE)
hr1 <- smoothHR(data=tmp, coxfit=fit) 
print(hr1)

# normal sample of kidney
nor.gaby.sam <- rownames(Sinfo[which(Sinfo$SampleType == "Normal" & Sinfo$Subtype %in% c("Bellini","Pediatric RCC")),])
indata <- read.table(file.path(res.path, "Wilms_tumor_Gaby_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

nor.gaby.sam <- rownames(Sinfo[which(Sinfo$SampleType == "Normal" & Sinfo$Subtype %in% c("Bellini","Pediatric RCC")),])
indata <- read.table(file.path(res.path, "Wilms_tumor_Gaby_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
sting.tum.gaby <- as.numeric(TPM["ENSG00000184584.8",RNA_Clust1])
sting.nor.gaby <- as.numeric(TPM["ENSG00000184584.8",nor.gaby.sam])
wilcox.test(sting.tum.gaby,sting.nor.gaby)
boxplot(sting.tum.gaby,sting.nor.gaby,main = "GABY")

nor.target.sam <- rownames(Sinfo[which(Sinfo$SampleType == "Normal" & Sinfo$Subtype %in% "TARGET_Wilm"),])
sting.tum.target <- as.numeric(TPM["ENSG00000184584.8",rownames(Sinfo.target)])
sting.nor.target <- as.numeric(TPM["ENSG00000184584.8",nor.target.sam])
boxplot(sting.tum.target,sting.nor.target)

# 16. comparsion of immune genes
indata <- read.table(file.path(res.path, "Wilms_tumor_Gaby_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
tmp1 <- as.numeric(indata["CD4",RNA_Clust1])/mean(as.numeric(indata["CD8A",RNA_Clust1])+as.numeric(indata["CD8B",RNA_Clust1]))
tmp2 <- as.numeric(indata["CD4",RNA_Clust2])/mean(as.numeric(indata["CD8A",RNA_Clust2])+as.numeric(indata["CD8B",RNA_Clust2]))
t.test(tmp1,tmp2)

#--------------------------------------------------------------------------------------------------------#
# Reanalysis in TARGET data by using 39 (37 expression matched) samples with detailed genomic alteration #
Sinfo.ccr39 <- read.table(file.path(data.path,"ccr39 genomic and clinical information.txt"),sep = "\t", check.names = F,stringsAsFactors = F,header = T,row.names = 1)
annCol.ccr39 <- Sinfo.ccr39[,c(11,12,13,14,17,20,25,28,31,32:44)]
# annCol.ccr39$Histology <- annCol.ccr39$`Histology Classification in Primary Tumor`
annCol.ccr39$Stage <- ifelse(annCol.ccr39$Stage %in% c("I","II"),"I_II",">III")
annCol.ccr39$Age2 <- ifelse(annCol.ccr39$`Age at Diagnosis in Days` > 2*365,">2","<=2")
annCol.ccr39$Age3 <- ifelse(annCol.ccr39$`Age at Diagnosis in Days` > 3*365,">3","<=3")
annCol.ccr39$Agem <- ifelse(annCol.ccr39$`Age at Diagnosis in Days` > median(annCol.ccr39$`Age at Diagnosis in Days`),">5","<=5")

annColors.ccr39 <- list()
annColors.ccr39[["TP53_MUT"]] <- annColors.ccr39[["TP53_LOSS"]] <- annColors.ccr39[["TP53_ALTER"]] <- annColors.ccr39[["WT1 Mut"]] <-
  annColors.ccr39[["WT1 copy number loss"]] <- annColors.ccr39[["CTNNB1 Mut"]] <- annColors.ccr39[["WTX Mut"]] <-
  annColors.ccr39[["WTX copy number loss"]] <- annColors.ccr39[["DGCR8 Mut"]] <- annColors.ccr39[["DROSHA Mut"]] <-
  annColors.ccr39[["SIX1/2 Mut"]] <- annColors.ccr39[["4q loss"]] <- annColors.ccr39[["1q gain"]] <-
  annColors.ccr39[["TP53_MUT"]] <- annColors.ccr39[["16q loss"]] <- annColors.ccr39[["14q loss"]] <-
  annColors.ccr39[["22 loss"]] <- c("Yes" = "black","No" = "white")
annColors.ccr39[["Vital Status"]] <- c("Dead" = "black","Alive" = "white")
annColors.ccr39[["TMB"]] <- blueyellow(64)
annColors.ccr39[["Gender"]] <- c("Female" = lightred, "Male" = cyan)
annColors.ccr39[["Stage"]] <- c("I_II" = darkblue, ">III" = sun)
annColors.ccr39[["Age2"]] <- c(">2" = nake, "<=2" = seagreen)
annColors.ccr39[["Age3"]] <- c(">3" = nake, "<=3" = seagreen)
annColors.ccr39[["Agem"]] <- c(">5" = nake, "<=5" = seagreen)

# 1. differential expression analysis in TARGET cohort regarding TP53 mutation
# mRNA
tmp <- Sinfo.ccr39$TP53_ALTER; names(tmp) <- rownames(Sinfo.ccr39)
Groupinfo <- Sinfo.ccr39[,c("TP53_MUT","TP53_LOSS","TP53_ALTER")]
colnames(Groupinfo)[3] <- "group"
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.TP53(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="ccr39_WT_mRNA_TP53", complist,PASSFlag=PASSFlag, overwt=TRUE)

# lncRNA
tmp <- Sinfo.ccr39$TP53_ALTER; names(tmp) <- rownames(Sinfo.ccr39)
Groupinfo <- Sinfo.ccr39[,c("TP53_MUT","TP53_LOSS","TP53_ALTER")]
colnames(Groupinfo)[3] <- "group"
PASSFlag <- rep(TRUE,length(Lids)); names(PASSFlag) <- Lids
complist <- createList.TP53(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Lids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Lids, featType="ccr39_WT_lncRNA_TP53", complist,PASSFlag=PASSFlag, overwt=TRUE)

# GSEA
deres <- read.table(file.path(res.path,"ccr39_WT_mRNA_TP53_edgeR_test_result.Yes_vs_No.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- deres$logFC
names(geneList) <- deres$id
geneList <- sort(geneList,decreasing = T)
gsea_ccr39_tp53_mutvswild <- GSEA(geneList = geneList,
                                  TERM2GENE=MSigDB,
                                  pvalueCutoff = 0.25,
                                  nPerm = 10000,
                                  seed = T,
                                  verbose=F)
res <- data.frame(gsea_ccr39_tp53_mutvswild)
write.table(as.data.frame(gsea_ccr39_tp53_mutvswild),file.path(res.path,"GSEA results for ccr39 WT mRNA TP53 Mutated vs Wild DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

# 2. immunity in TARGET cohort regarding TP53
indata <- log2(TPM[Mids,rownames(annCol.ccr39)] + 1)
indata <- indata[rowSums(indata) > 0,]
indata$symbol <- Ginfo[rownames(indata),"genename"]
indata <- apply(indata[,setdiff(colnames(indata), "symbol")], 2, function(x) tapply(x, INDEX=factor(indata$symbol), FUN=median, na.rm=TRUE))
write.table(indata,file = file.path(res.path,"Wilms_tumor_ccr39_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# estimate
library(estimate)
filterCommonGenes(input.f=file.path(res.path, "Wilms_tumor_ccr39_log2TPM_hugo.txt") , output.f=file.path(res.path,"Wilms_tumor_ccr39_log2TPM_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"Wilms_tumor_ccr39_log2TPM_hugo_ESTIMATE.txt"), file.path(res.path,"Wilms_tumor_ccr39_log2TPM_hugo_estimate_score.txt"), platform="affymetrix")
est.ccr39 <- read.table(file = file.path(res.path,"Wilms_tumor_ccr39_log2TPM_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.ccr39) <- est.ccr39[,2]; colnames(est.ccr39) <- est.ccr39[1,]; est.ccr39 <- est.ccr39[-1,c(-1,-2)];
est.ccr39 <- sapply(est.ccr39, as.numeric); rownames(est.ccr39) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.ccr39.backup = as.data.frame(est.ccr39); colnames(est.ccr39.backup) <- colnames(indata)
est.ccr39 <- annTrackScale(indata = est.ccr39, halfwidth = 2, poolsd = F); est.ccr39 <- as.data.frame(t(est.ccr39))
rownames(est.ccr39) <- colnames(indata)

# MCPcounter
MCPscore.ccr39 <- MCPcounter.estimate(expression = indata,featuresType = "HUGO_symbols")
MCPscore.ccr39 <- as.data.frame(MCPscore.ccr39); MCPscore.ccr39.raw.backup <- MCPscore.ccr39
MCPscore.ccr39 <- annTrackScale(indata = MCPscore.ccr39, halfwidth = 2, poolsd = F); MCPscore.ccr39 <- as.data.frame(t(MCPscore.ccr39))
MCPscore.ccr39.pure <- as.data.frame(MCPscore.ccr39.raw.backup)
for (sam in colnames(MCPscore.ccr39.pure)) {
  for (cell in rownames(MCPscore.ccr39.pure)) {
    # MCPscore.ccr39.pure[cell,sam] <- MCPscore.ccr39.raw.backup[cell,sam]/(1-est.ccr39.backup[4,sam])
    MCPscore.ccr39.pure[cell,sam] <- MCPscore.ccr39.raw.backup[cell,sam]*est.ccr39.backup[4,sam]
  }
}

# ccr immune signature
library(GSVA)
indata <- read.table(file.path(res.path, "Wilms_tumor_ccr39_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
immune.signature <- read.table(file.path(comAnn.path,"CCR_Curated_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- intersect(toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"]),rownames(indata))
}
imm.score.ccr39.wt <- gsva(as.matrix(indata),immune.sig.ccr,method="gsva",parallel.sz = 1)
imm.score.ccr39.wt.pure <- imm.score.ccr39.wt
for (sam in colnames(imm.score.ccr39.wt)) {
  for (cell in rownames(imm.score.ccr39.wt)) {
    imm.score.ccr39.wt.pure[cell,sam] <- imm.score.ccr39.wt[cell,sam]*est.ccr39.backup[4,sam]
  }
}

# heatmap
pdf(file.path(fig.path,"MCPcounter heatmap of ccr39 data regarding TP53 status.pdf"),width = 6,height = 3)
plotdata <- standarize.fun(MCPscore.ccr39.pure,halfwidth = 2)
plotdata <- annTrackScale(MCPscore.ccr39.pure,halfwidth = 2,poolsd = F)

plotdata <- t(MCPscore.ccr39)
pheatmap(plotdata[1:8,],
         border_color = NA,
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annCol.ccr39,
         annotation_colors = annColors.ccr39,
         show_rownames = T,
         show_colnames = F,
         color = greenred(64))
invisible(dev.off())

pdf(file.path(fig.path,"CCR immune signature heatmap of ccr39 data in TARGET cohort.pdf"),width = 10,height = 20)
plotdata <- standarize.fun(imm.score.ccr39.wt.pure,halfwidth = 1)
hcs <- hclust(distanceMatrix(as.matrix(plotdata), "euclidean"), "ward.D")
group <- cutree(hcs,2)
annCol.ccr39$ccr_Clust <- paste0("C",group[rownames(annCol.ccr39)])
annColors.ccr39[["ccr_Clust"]] <- c("C1" = jco[1],"C2" = jco[2])
# plotdata <- annTrackScale(imm.score.ccr39.wt,halfwidth = 2,poolsd = F)
pheatmap(plotdata,
         border_color = NA,
         cluster_rows = F,
         cluster_cols = dendsort(hcs),
         annotation_col = annCol.ccr39,
         annotation_colors = annColors.ccr39,
         cellwidth = 10,
         cellheight = 10,
         show_rownames = T,
         show_colnames = F,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64))
         #color = greenred(64))
invisible(dev.off())

tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  ccr_Clust = annCol.ccr39[rownames(Sinfo.target),"ccr_Clust"],
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ ccr_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ ccr_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jama",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by ccr39 data in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

#-------------------------------------------------------------------------------------------------------------------------#
# Reanalysis in TARGET data (cbioportal) by using 101 samples with frequent mutations alteration and clinical information #
library(tableone)

annCol.cbio <- read.delim(file.path(data.path,"mutation cbioportal target data 101 cases simple.txt"),sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
annColors.cbio <- list()
for (i in colnames(annCol.cbio)[1:17]) {
  annColors.cbio[[i]] <- c("Yes" = "black","No" = "white")
}
indata <- log2(TPM[Mids,rownames(annCol.cbio)] + 1)
indata <- indata[rowSums(indata) > 0,]
indata$symbol <- Ginfo[rownames(indata),"genename"]
indata <- apply(indata[,setdiff(colnames(indata), "symbol")], 2, function(x) tapply(x, INDEX=factor(indata$symbol), FUN=median, na.rm=TRUE))
write.table(indata,file = file.path(res.path,"Wilms_tumor_cbio_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# estimate
library(estimate)
filterCommonGenes(input.f=file.path(res.path, "Wilms_tumor_cbio_log2TPM_hugo.txt") , output.f=file.path(res.path,"Wilms_tumor_cbio_log2TPM_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"Wilms_tumor_cbio_log2TPM_hugo_ESTIMATE.txt"), file.path(res.path,"Wilms_tumor_cbio_log2TPM_hugo_estimate_score.txt"), platform="affymetrix")
est.cbio <- read.table(file = file.path(res.path,"Wilms_tumor_cbio_log2TPM_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.cbio) <- est.cbio[,2]; colnames(est.cbio) <- est.cbio[1,]; est.cbio <- est.cbio[-1,c(-1,-2)];
est.cbio <- sapply(est.cbio, as.numeric); rownames(est.cbio) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.cbio.backup = as.data.frame(est.cbio); colnames(est.cbio.backup) <- colnames(indata)
est.cbio <- annTrackScale(indata = est.cbio, halfwidth = 2, poolsd = F); est.cbio <- as.data.frame(t(est.cbio))
rownames(est.cbio) <- colnames(indata)

# MCPcounter
MCPscore.cbio <- MCPcounter.estimate(expression = indata,featuresType = "HUGO_symbols")
MCPscore.cbio <- as.data.frame(MCPscore.cbio); MCPscore.cbio.raw.backup <- MCPscore.cbio
MCPscore.cbio <- annTrackScale(indata = MCPscore.cbio, halfwidth = 2, poolsd = F); MCPscore.cbio <- as.data.frame(t(MCPscore.cbio))
MCPscore.cbio.pure <- as.data.frame(MCPscore.cbio.raw.backup)
for (sam in colnames(MCPscore.cbio.pure)) {
  for (cell in rownames(MCPscore.cbio.pure)) {
    # MCPscore.cbio.pure[cell,sam] <- MCPscore.cbio.raw.backup[cell,sam]/(1-est.cbio.backup[4,sam])
    MCPscore.cbio.pure[cell,sam] <- MCPscore.cbio.raw.backup[cell,sam]/(1-est.cbio.backup[4,sam])
  }
}

# ccr immune signature
library(GSVA)
indata <- read.table(file.path(res.path, "Wilms_tumor_cbio_log2TPM_hugo.txt"), sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
immune.signature <- read.table(file.path(comAnn.path,"CCR_Curated_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- intersect(toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"]),rownames(indata))
}
imm.score.cbio.wt <- gsva(as.matrix(indata),immune.sig.ccr,method="gsva",parallel.sz = 1)
imm.score.cbio.wt.pure <- imm.score.cbio.wt
for (sam in colnames(imm.score.cbio.wt)) {
  for (cell in rownames(imm.score.cbio.wt)) {
    imm.score.cbio.wt.pure[cell,sam] <- imm.score.cbio.wt[cell,sam]/(1-est.cbio.backup[4,sam])
  }
}

plotdata <- standarize.fun(imm.score.cbio.wt.pure,halfwidth = 1)
hcs <- hclust(distanceMatrix(as.matrix(imm.score.cbio.wt.pure), "euclidean"), "ward.D")
group <- cutree(hcs,2)
annCol.cbio$ccr_Clust <- paste0("C",group[rownames(annCol.cbio)])
stabl <- CreateTableOne(vars=colnames(annCol.cbio),strata="ccr_Clust",
                        data=annCol.cbio,factorVars=colnames(annCol.cbio))
print(stabl,showAllLevels = TRUE)

pdf(file.path(fig.path, "CCR immune signature heatmap of cbioportal 102 wilms tumors in target cohort.pdf"), height=12,width = 10)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = dendsort(hcs),
         border_color = NA,
         annotation_col = annCol.cbio,
         annotation_colors = annColors.cbio,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         # color = viridisLite::inferno(64),
         # color = viridisLite::viridis(64),
         # color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 4,
         cellheight = 12,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  ccr_Clust = annCol.cbio[rownames(Sinfo.target),"ccr_Clust"],
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ ccr_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ ccr_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jama",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

#-----------------------------------------------#
# Analysis Genomic alteration of Gaby WT cohort #

gaby.wt.genomic.alter <- read.table(file.path(data.path,"Gaby_WT_CNA_data.txt"),sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
gaby.wt.mlpa <- read.table(file.path(data.path,"Gaby_WT_MLPA_data.txt"),sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)

# 1. comprehensive heatmap for immune signature and genomic alteration
annCol.wt.gaby.genomic.alter <- annCol.wt.gaby[,c(1,14),drop = F]
annCol.wt.gaby.genomic.alter$TLS <- annTrackScale(TLS.gaby[rownames(annCol.wt.gaby.genomic.alter)],halfwidth = 2)
annCol.wt.gaby.genomic.alter$OS <- gaby.wt.mlpa[rownames(annCol.wt.gaby.genomic.alter),"Status"]
annCol.wt.gaby.genomic.alter <- cbind.data.frame(annCol.wt.gaby.genomic.alter,gaby.wt.genomic.alter[rownames(annCol.wt.gaby.genomic.alter),])
annCol.wt.gaby.genomic.alter <- cbind.data.frame(annCol.wt.gaby.genomic.alter,gaby.wt.mlpa[rownames(annCol.wt.gaby.genomic.alter),7:20])

annColors.wt.gaby.genomic.alter <- annColors.wt.gaby[c("Subtype","RNA_Clust")]
annColors.wt.gaby.genomic.alter[["1q gain"]] <- annColors.wt.gaby.genomic.alter[["4q loss"]] <-
  annColors.wt.gaby.genomic.alter[["14q loss"]] <- annColors.wt.gaby.genomic.alter[["16q loss"]] <-
  annColors.wt.gaby.genomic.alter[["22 loss"]] <- annColors.wt.gaby.genomic.alter[["TP53_Mut"]] <- 
  annColors.wt.gaby.genomic.alter[["TP53 mutations"]] <- annColors.wt.gaby.genomic.alter[["DROSA"]] <- c("Yes" = "black","No" = "white")
annColors.wt.gaby.genomic.alter[["1p STATUS"]] <- annColors.wt.gaby.genomic.alter[["1q STATUS"]] <-
  annColors.wt.gaby.genomic.alter[["MYCN STATUS"]] <- annColors.wt.gaby.genomic.alter[["2p flanking probe status"]] <-
  annColors.wt.gaby.genomic.alter[["2q flanking probe status"]] <- annColors.wt.gaby.genomic.alter[["FBXW7 STATUS"]] <-
  annColors.wt.gaby.genomic.alter[["4q flanking probe status"]] <- annColors.wt.gaby.genomic.alter[["WT1 STATUS"]] <-
  annColors.wt.gaby.genomic.alter[["IGF2 flanking probe status"]] <- annColors.wt.gaby.genomic.alter[["16p STATUS"]] <-
  annColors.wt.gaby.genomic.alter[["16q STATUS"]] <- annColors.wt.gaby.genomic.alter[["TP53 STATUS"]] <-
  annColors.wt.gaby.genomic.alter[["17p flanking probe status"]] <- c("GAIN" = sun,"LOSS" = darkblue,"NORMAL" = "white")
annColors.wt.gaby.genomic.alter[["OS"]] <- c("alive" = "grey90","dead" = brown)
annColors.wt.gaby.genomic.alter[["TLS"]] <- blueyellow(64)

tmp <- gaby.wt.sinfo[colnames(mygene)[1:(ncol(mygene)-1)],c("Age at diagnostic (years)","Sexe","Risk (review/local)","Histology type (review/local)","Overall Stage (review/local)")]
colnames(tmp) <- c("Age","Gender","Risk","Histology","Stage")
tmp$Age <- round(tmp$Age,2)
tmp$Age <- ifelse(tmp$Age > 4.85, ">4.85", "<=4.85")
tmp$Histology <- gsub(" (311)", "", tmp$Histology, fixed = T)
tmp$Histology <- gsub(" (312)", "", tmp$Histology, fixed = T)
tmp$Gender <- ifelse(tmp$Gender == "F","Female","Male")
tmp$Risk <- gsub(" Risk","",tmp$Risk)
annCol.wt.gaby.genomic.alter <- cbind.data.frame(annCol.wt.gaby.genomic.alter, tmp[rownames(annCol.wt.gaby.genomic.alter),])
annCol.wt.gaby.genomic.alter$Purity <- as.numeric(est.backup[4,rownames(annCol.wt.gaby.genomic.alter)])

annColors.wt.gaby.genomic.alter[["Age"]] = c(">4.85" = orange, "<=4.85" = alpha(orange,0.6))
annColors.wt.gaby.genomic.alter[["Gender"]] = c("Female" = jco[2], "Male" = jco[1])
annColors.wt.gaby.genomic.alter[["Risk"]] = c("High" = darkred, "Intermediate" = alpha(darkred,0.6))
annColors.wt.gaby.genomic.alter[["Histology"]] = c("ANA diffuse" = darkblue,"Focal anaplasia" = darkred)
annColors.wt.gaby.genomic.alter[["Stage"]] = c("II" = alpha(seagreen,0.3), "III" = alpha(seagreen,0.6), "IV" = seagreen)
annColors.wt.gaby.genomic.alter[["Purity"]] <- NMF:::ccRamp(c("white","black"),n=64)
plotdata <- standarize.fun(imm.score.gaby.wt3.pure[,rownames(cluster2.wt.gaby.mRNA.ans$annCol)],halfwidth = 2)
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D2")
pdf(file.path(fig.path, "CCR immune signature heatmap of wilms tumors in gaby cohort (purified by tp).pdf"), height=26,width = 10)
pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = "black",
         annotation_col = annCol.wt.gaby.genomic.alter[,c("RNA_Clust","Histology","OS","Age","Gender","Risk","Stage","TP53 mutations","TP53 STATUS","FBXW7 STATUS","MYCN STATUS")],
         annotation_colors = annColors.wt.gaby.genomic.alter[c("RNA_Clust","Histology","OS","Age","Gender","Risk","Stage","TP53 mutations","TP53 STATUS","FBXW7 STATUS","MYCN STATUS")],
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         # color = viridisLite::inferno(64),
         # color = viridisLite::viridis(64),
         # color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 15,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         cutree_cols = 2,
         gaps_row = c(14,22),
         fontsize_col = 8)
invisible(dev.off())

indata <- gaby.wt[TLS,rownames(cluster2.wt.gaby.mRNA.ans$annCol)]
plotdata <- standarize.fun(indata,halfwidth = 2)
hm1 <- pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = "black",
         annotation_col = annCol.wt.gaby.genomic.alter[,c("RNA_Clust","TLS")],
         annotation_colors = annColors.wt.gaby.genomic.alter[c("RNA_Clust","TLS")],
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         #color = viridisLite::inferno(64),
         #color = viridisLite::viridis(64),
         # color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 15,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         cutree_cols = 2,
         fontsize_col = 8)
indata <- gaby.wt[c("TMEM173","MB21D1"),rownames(cluster2.wt.gaby.mRNA.ans$annCol)]
plotdata <- standarize.fun(indata,halfwidth = 2)
hm2 <- pheatmap(as.matrix(plotdata),
                cluster_rows = F,
                cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
                border_color = "black",
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                #color = viridisLite::inferno(64),
                #color = viridisLite::viridis(64),
                # color = greenred(64),
                treeheight_row = 15,
                treeheight_col = 15,
                cellwidth = 15,
                cellheight = 10,
                show_colnames = T,
                show_rownames = T,
                cutree_cols = 2,
                fontsize_col = 8)
plotdata <- standarize.fun(gaby.cgassting,halfwidth = 2)
hm3 <- pheatmap(as.matrix(plotdata),
                cluster_rows = F,
                cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
                border_color = "black",
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                #color = viridisLite::inferno(64),
                #color = viridisLite::viridis(64),
                # color = greenred(64),
                treeheight_row = 15,
                treeheight_col = 15,
                cellwidth = 15,
                cellheight = 10,
                show_colnames = T,
                show_rownames = T,
                cutree_cols = 2,
                fontsize_col = 8)
pdf(file.path(fig.path, "TLS heatmap of wilms tumors in gaby cohort.pdf"), height=5,width = 9)
draw(hm1 %v% hm2 %v% hm3)
invisible(dev.off())

# 2. comparison for TMB
tmp1 <- as.numeric(gaby.wt.mlpa[RNA_Clust1,"TMB"])
tmp2 <- as.numeric(gaby.wt.mlpa[RNA_Clust2,"TMB"])
wilcox.test(tmp1,tmp2) #0.6714

# 3. comparsion of CD4 memory activated
tmp1 <- imm.score.gaby.wt3.pure["T.cells.CD4.memory.activated",RNA_Clust1]
tmp2 <- imm.score.gaby.wt3.pure["T.cells.CD4.memory.activated",RNA_Clust2]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Infiltrated","Desert"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for CD4 memory activated in gaby WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Activated memory CD4+ T cells",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,3,"P = 0.032",cex = 1)
invisible(dev.off())

outTab <- NULL
for(cell in rownames(imm.score.gaby.wt3.pure)) {
  tmp1 <- as.numeric(imm.score.gaby.wt3.pure[cell,RNA_Clust1])
  tmp2 <- as.numeric(imm.score.gaby.wt3.pure[cell,RNA_Clust2])
  wilcox.p <- wilcox.test(tmp1,tmp2)$p.value
  outTab <- rbind.data.frame(outTab,data.frame(cell = cell,
                                               avg.C1 = mean(tmp1),
                                               avg.C2 = mean(tmp2),
                                               p.wilcox = wilcox.p,
                                               stringsAsFactors = F),
                             stringsAsFactors = F)
}
outTab$diff1.2 <- outTab$avg.C1-outTab$avg.C2
outTab$FDR <- p.adjust(outTab$p.wilcox,method = "BH")
write.table(outTab,file.path(res.path,"comparision of CCR immune signature between two gaby clusters.txt"),sep = "\t",row.names = F,quote = F)

#------------------------------------------------------------------------#
# Reanalysis immunity enrichment with consideration of tumor purity (TP) #

# 1. distribution of TP in Gaby and TARGET (n=114) cohorts
pdf(file.path(fig.path,"density plot of tumor purity.pdf"),width = 4,height = 3.5)
par(bty="o", mgp = c(2,0.3,0), mar = c(3.1,3.1,2.1,3.1),tcl=-.25, font.main=3)

dens.gaby <- density(as.numeric(est.backup[4,]))
plot(dens.gaby$x,dens.gaby$y, col=ggplot2::alpha(jco[1],0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n",xlim = c(0,1),ylim = c(0,13))
polygon(dens.gaby$x,dens.gaby$y,col = ggplot2::alpha(jco[1],0.5),border = ggplot2::alpha(jco[1],0.5)) 
axis(side=1, at = c(0,0.5,0.7,0.9,1),label = c(0,0.5,0.7,0.9,1))
rug(as.numeric(est.backup[4,]),side = 1,lwd = 1.5,col = jco[1])
mtext("Distribution of tumor purity", side=1, line=1.5)
mtext("Density", side=2, line=1.5)

par(new = T)
dens.target <- density(as.numeric(est.target.backup[4,]))
plot(dens.target$x,dens.target$y, col=ggplot2::alpha(jco[2],0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n",xlim = c(0,1),ylim = c(0,13))
polygon(dens.target$x,dens.target$y,col = ggplot2::alpha(jco[2],0.5),border = ggplot2::alpha(nake,0.5)) 
axis(side=2, at = pretty(range(dens.target$y))[-length(pretty(range(dens.target$y)))],las = 1)
rug(as.numeric(est.target.backup[4,]),side = 1,lwd = 1.5,col = jco[2])
legend("topleft",
       legend = c("GABY","TARGET"),
       fill = jco[1:2],
       x.intersp = 0.3,
       y.intersp = 0.8,
       border = "white",
       bty = "n")
invisible(dev.off())

# 2. remove effect of tumor purity
tune.tp <- function(tp.vec = NULL, tp.cutoff = 0.98, exp.df = NULL, sig.list = NULL, gsva.method = "gsva", core.use = 1) {
  
  # tp.vec is a named vector storing the numerical tumor purity ranging from 0-1 (estimated by estimate R package)
  # tp.cutoff is a threshold to decide a sample who has a tumor purity over than this threshold will be removed for microenviroment inference
  # exp.df is the expression data frame that used to perform enrichment analysis
  # sig.list is a list storing signatures used to infer enrichment score
  # gsva.method is the method that used to perform enrichment analysis
  
  library(GSVA)
  
  if(!all(names(tp.vec) == colnames(exp.df))) {stop("Mismatched in sample ID!\n")}
  rm.sam <- names(tp.vec)[tp.vec > tp.cutoff]
  mean.tp.rm.sam <- mean(tp.vec[rm.sam])
  cat(length(rm.sam),"samples has been removed!\n")
  cat("The mean tumor purity for removed samples is",mean.tp.rm.sam,"\n")
  
  keep.sam <- names(tp.vec)[tp.vec < tp.cutoff]
  
  exp <- exp.df[,keep.sam]
  es.raw <- gsva(as.matrix(exp),sig.list,method=gsva.method,parallel.sz = core.use)
  es.pure <- es.raw
  for (sam in colnames(es.raw)) {
    for (cell in rownames(es.raw)) {
      es.pure[cell,sam] <- es.pure[cell,sam]/(1-tp.vec[sam])
    }
  }
  cat("The range of raw enrichment score matrix is ",paste0(round(range(es.raw),3),collapse = " to "),"\n")
  cat("The range of puried enrichment score matrix is ",paste0(round(range(es.pure),3),collapse = " to "),"\n")
  
  return(list(keep.sam = keep.sam,
              es.raw = es.raw,
              es.pure = es.pure,
              range.es.raw = round(range(es.raw),3),
              range.es.pure = round(range(es.pure),3)))
}

# tuning target data (cbioportal 101)
TPM.TARGET.HUGO <- read.table(file.path(res.path,"Wilms_tumor_TARGET_log2TPM_hugo.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tp.target <- as.numeric(est.target.backup[4,]); names(tp.target) <- colnames(est.target.backup)

target.pure <- tune.tp(tp.vec = tp.target,
                       exp.df = TPM.TARGET.HUGO,
                       gsva.method = "gsva",
                       sig.list = immune.sig.ccr,
                       tp.cutoff = 0.98)
target.pure.me <- tune.tp(tp.vec = tp.target,
                       exp.df = TPM.TARGET.HUGO,
                       gsva.method = "gsva",
                       sig.list = immune.sig.me,
                       tp.cutoff = 0.98)
target.pure95 <- tune.tp(tp.vec = tp.target,
                         exp.df = TPM.TARGET.HUGO,
                         gsva.method = "gsva",
                         sig.list = immune.sig.ccr,
                         tp.cutoff = 0.95)

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.cbio))
indata <- target.pure$es.pure[,com_sam]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D2")
group <- cutree(hcs,2)
annCol.cbio$tune_Clust <- NA
annCol.cbio[com_sam,"tune_Clust"] <- paste0("C",group[com_sam])
annCol.cbio$Tumor_purity <- as.numeric(est.target.backup[4,rownames(annCol.cbio)])
annCol.cbio$Histology <- annCol.cbio$`Histology Classification in Primary Tumor`
annCol.cbio$Stage <- annCol.cbio$`Neoplasm American Joint Committee on Cancer Clinical Group Stage`
annCol.cbio$Stage <- ifelse(annCol.cbio$Stage == "I","I",
                            ifelse(annCol.cbio$Stage == "II","II",
                            ifelse(annCol.cbio$Stage %in% c("III","IIIB"),"III","IV")))
annCol.cbio$Age2 <- ifelse(annCol.cbio$`Diagnosis Age` > 2,">2","<=2")
annCol.cbio$Age3 <- ifelse(annCol.cbio$`Diagnosis Age` > 3,">3","<=3")
annCol.cbio$Agem <- ifelse(annCol.cbio$`Diagnosis Age` > median(annCol.cbio$`Diagnosis Age`),">5","<=5")

annColors.cbio[["Tumor_purity"]] <- NMF:::ccRamp(c("white","black"),n=64)
annColors.cbio[["tune_Clust"]] <- c("C1" = jco[1],"C2" = jco[2], "C3" = jco[4])
annColors.cbio[["Fraction Genome Altered"]] <- blueyellow(64)
annColors.cbio[["Mutation Count"]] <- blueyellow(64)
annColors.cbio[["Sex"]] <- c("Female" = lightred, "Male" = cyan)
annColors.cbio[["Histology"]] <- c("FHWT" = sun, "DAWT" = darkblue)
annColors.cbio[["Stage"]] <- c("I" = darkblue, "II" = cyan, "III" = lightred, "IV" = sun)
annColors.cbio[["Age2"]] <- c(">2" = nake, "<=2" = seagreen)
annColors.cbio[["Age3"]] <- c(">3" = nake, "<=3" = seagreen)
annColors.cbio[["Agem"]] <- c(">5" = nake, "<=5" = seagreen)

tmp <- annCol.cbio[,c(1:17,28,32,41:48)]
stabl <- CreateTableOne(vars=colnames(tmp),strata="tune_Clust",
                        data=annCol.cbio,factorVars=colnames(tmp)[c(1:17,20,23,24,25,26,27)])
print(stabl,showAllLevels = TRUE)
write.table(print(stabl,showAllLevels = TRUE),file.path(res.path,"statistic between CCR immune signatures groups of cbioportal 101 wilms tumors.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

pdf(file.path(fig.path, "CCR immune signature heatmap of cbioportal 101 wilms tumors in target cohort (purified by tp).pdf"), height=12,width = 10)
pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
         cluster_rows = F,
         cluster_cols = dendsort(hcs),
         border_color = NA,
         annotation_col = annCol.cbio[colnames(plotdata),c(42:48,1:17)],
         annotation_colors = annColors.cbio,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 4,
         cellheight = 12,
         gaps_row = c(14,22),
         cutree_cols = 3,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

# fraction genome altered comparison
tmp1 <- annCol.cbio[which(annCol.cbio$tune_Clust == "C1"),"Fraction Genome Altered"]
tmp2 <- annCol.cbio[which(annCol.cbio$tune_Clust == "C2"),"Fraction Genome Altered"]
wilcox.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Desert","Infiltrated"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for Fraction Genome Altered in cbioportal 101 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Fraction Genome Altered",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.7,"P = 0.047",cex = 1)
invisible(dev.off())

# comparsion of TMB
tmp1 <- annCol.cbio[which(annCol.cbio$tune_Clust == "C1"),"Mutation Count"]
tmp2 <- annCol.cbio[which(annCol.cbio$tune_Clust == "C2"),"Mutation Count"]
wilcox.test(tmp1,tmp2)

# comparsion of CD4 memory activated
tmp1 <- rownames(annCol.cbio[which(annCol.cbio$tune_Clust == "C1"),])
tmp2 <- rownames(annCol.cbio[which(annCol.cbio$tune_Clust == "C2"),])
tmp1 <- indata["T.cells.CD4.memory.activated",tmp1]
tmp2 <- indata["T.cells.CD4.memory.activated",tmp2]
wilcox.test(tmp1,tmp2)
t.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Desert","Infiltrated"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for CD4 memory activated in cbioportal 101 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Activated memory CD4+ T cells",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,18,"P = 0.025",cex = 1)
invisible(dev.off())

# survival comparsion
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.cbio[rownames(Sinfo.target),"tune_Clust"],
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jco",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by cbioportal 101 (purified by tp) data in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

# correlation between tumor purity and survival
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.cbio[rownames(Sinfo.target),"tune_Clust"],
                  tp = annCol.cbio[rownames(Sinfo.target),"Tumor_purity"],
                  stringsAsFactors = F)
tmp <- na.omit(tmp)
coxph(Surv(OS.time,OS)~tp,data = tmp) # tumor purity is not associated with survival

# integrative clustering for cbioportal target samples
library(iClusterPlus)
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.cbio))
indata <- target.pure$es.pure[,com_sam]

iMut <- annCol.cbio[com_sam,1:17]
iMut[iMut == "Yes"] <- 1
iMut[iMut == "No"] <- 0
iMut <- sapply(iMut, as.numeric); rownames(iMut) <- com_sam
iExp <- as.matrix(t(indata))
bayfit.cbio <- tune.iClusterBayes(cpus=1,dt1=iMut,dt2=iExp,
                                  type=c("binomial","gaussian"),
                                  K=1:3,n.burnin=18000,n.draw=12000,
                                  prior.gamma=c(0.5,0.5,0.5),sdev=0.05,thin=3)
allBIC = NULL
devratio = NULL
nK = length(bayfit.cbio$fit)
for(i in 1:nK){
  allBIC = c(allBIC,bayfit.cbio$fit[[i]]$BIC)
  devratio = c(devratio,bayfit.cbio$fit[[i]]$dev.ratio)
}

plotHMBayes(fit=bayfit.cbio$fit[[1]],datasets=list(iMut,iExp),
              type=c("binomial","gaussian"), 
              threshold=c(0.5,0.5),row.order=c(F,F),
              sparse=c(T,T),cap=c(F,T))

iMut.prob <- bayfit.cbio$fit[[1]]$beta.pp[[1]]; names(iMut.prob) <- colnames(iMut)
iExp.prob <- bayfit.cbio$fit[[1]]$beta.pp[[2]]; names(iExp.prob) <- colnames(iExp)

annCol.cbio$iClust <- NA
annCol.cbio[com_sam,"iClust"] <- paste0("iC",bayfit.cbio$fit[[1]]$clusters)
annColors.cbio[["iClust"]] <- c("iC1" = darkblue,"iC2" = sun)
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  iClust = annCol.cbio[rownames(Sinfo.target),"iClust"],
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ iClust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 

iClust1 <- rownames(annCol.cbio[which(annCol.cbio$iClust == "iC1"),])
iClust2 <- rownames(annCol.cbio[which(annCol.cbio$iClust == "iC2"),])

indata <- target.pure$es.pure[,com_sam]
plotdata <- standarize.fun(indata,halfwidth = 2)
pdf(file.path(fig.path,"iCluster CCR immune signature heatmap of cbioportal 102 wilms tumors in target cohort (purified by tp).pdf"), height=12,width = 10)
pheatmap(as.matrix(plotdata[immune.sig.ccr.order,c(iClust1,iClust2)]),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         annotation_col = annCol.cbio[colnames(plotdata),c(20,18,19,1:17)],
         annotation_colors = annColors.cbio,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 4,
         cellheight = 12,
         gaps_row = c(14,22),
         gaps_col = length(iClust1),
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

# tuning target data (ccr 39)
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.ccr39))
indata <- target.pure$es.pure[,com_sam]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D")
group <- cutree(hcs,2)
annCol.ccr39$tune_Clust <- NA
annCol.ccr39[com_sam,"tune_Clust"] <- paste0("C",group[com_sam])
annCol.ccr39$Tumor_purity <- as.numeric(est.target.backup[4,rownames(annCol.ccr39)])
annColors.ccr39[["tune_Clust"]] <- c("C1" = jco[1],"C2" = jco[2])
annColors.ccr39[["Tumor_purity"]] <- NMF:::ccRamp(c("white","black"),n=64)

tmp <- annCol.ccr39[,c(1:4,6,7,10:27)]
stabl <- CreateTableOne(vars=colnames(tmp),strata="tune_Clust",
                        data=tmp,factorVars=colnames(tmp)[c(1:22)])
print(stabl,showAllLevels = TRUE)
write.table(print(stabl,showAllLevels = TRUE),file.path(res.path,"statistic between CCR immune signatures groups of ccr 39 wilms tumors.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

pdf(file.path(fig.path, "CCR immune signature heatmap of ccr 39 wilms tumors in target cohort (purified by tp).pdf"), height=25,width = 10)
pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
         cluster_rows = F,
         cluster_cols = dendsort(hcs),
         border_color = "black",
         annotation_col = annCol.ccr39[colnames(plotdata),c(26,27,23:25,7,4,1,2,3,10:22)],
         annotation_colors = annColors.ccr39,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 8,
         cellheight = 10,
         gaps_row = c(14,22),
         cutree_cols = 2,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

# survival analysis
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.ccr39[rownames(Sinfo.target),"tune_Clust"],
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jco",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by ccr 39 (purified by tp) data in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.ccr39[rownames(Sinfo.target),"tune_Clust"],
                  tp = annCol.ccr39[rownames(Sinfo.target),"Tumor_purity"],
                  stringsAsFactors = F)
tmp <- na.omit(tmp)
coxph(Surv(OS.time,OS)~tp,data = tmp) # tumor purity is not associated with survival

# comparsion of CD4 memory activated
tmp1 <- rownames(annCol.ccr39[which(annCol.ccr39$tune_Clust == "C1"),])
tmp2 <- rownames(annCol.ccr39[which(annCol.ccr39$tune_Clust == "C2"),])
tmp1 <- indata["T.cells.CD4.memory.activated",tmp1]
tmp2 <- indata["T.cells.CD4.memory.activated",tmp2]
wilcox.test(tmp1,tmp2)
t.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Desert","Infiltrated"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for CD4 memory activated in ccr 39 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Activated memory CD4+ T cells",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,18,"P = 0.237",cex = 1)
invisible(dev.off())

# tuning all target data 114
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target))
indata <- target.pure$es.pure[,com_sam]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D")
group <- cutree(hcs,2)
annCol.target$tune_Clust <- NA
annCol.target[com_sam,"tune_Clust"] <- paste0("C",group[com_sam])
annCol.target$Tumor_purity <- as.numeric(est.target.backup[4,rownames(annCol.target)])
annColors.target[["tune_Clust"]] <- c("C1" = jco[1],"C2" = jco[2],"C3" = jco[3])
annColors.target[["Tumor_purity"]] <- NMF:::ccRamp(c("white","black"),n=64)

fisher.test(table(annCol.target$tune_Clust,annCol.target$HISTOLOGY))

pdf(file.path(fig.path, "CCR immune signature heatmap of 114 wilms tumors in target cohort (purified by tp).pdf"), height=8,width = 10)
pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
         cluster_rows = F,
         cluster_cols = dendsort(hcs),
         border_color = NA,
         annotation_col = annCol.target[colnames(plotdata),c(7,8,4,1:3)],
         annotation_colors = annColors.target,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 3,
         cellheight = 10,
         gaps_row = c(14,22),
         cutree_cols = 2,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

# comparsion of CD4 memory activated
tmp1 <- rownames(annCol.target[which(annCol.target$tune_Clust == "C1"),])
tmp2 <- rownames(annCol.target[which(annCol.target$tune_Clust == "C2"),])
tmp1 <- indata["T.cells.CD4.memory.activated",tmp1]
tmp2 <- indata["T.cells.CD4.memory.activated",tmp2]
wilcox.test(tmp1,tmp2)
t.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Desert","Infiltrated"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for CD4 memory activated in target 114 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Activated memory CD4+ T cells",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,18,"P = 0.024",cex = 1)
invisible(dev.off())

# survival analysis
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.target[rownames(Sinfo.target),"tune_Clust"],
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jco",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

# sensitivity analysis by using 95% cutoff
com_sam <- intersect(target.pure95$keep.sam,rownames(annCol.target))
indata <- target.pure95$es.pure[,com_sam]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D")
group <- cutree(hcs,2)
annCol.target$tune_Clust <- NA
annCol.target[com_sam,"tune_Clust"] <- paste0("C",group[com_sam])
annCol.target$Tumor_purity <- as.numeric(est.target.backup[4,rownames(annCol.target)])
annColors.target[["tune_Clust"]] <- c("C1" = jco[1],"C2" = jco[2],"C3" = jco[3])
annColors.target[["Tumor_purity"]] <- NMF:::ccRamp(c("white","black"),n=64)

fisher.test(table(annCol.target$tune_Clust,annCol.target$HISTOLOGY)) 
fisher.test(table(annCol.target$tune_Clust,annCol.target$TP53_ALTER))
fisher.test(table(annCol.target$tune_Clust,annCol.target$TP53_MUT))

pdf(file.path(fig.path, "CCR immune signature heatmap of 114 wilms tumors in target cohort (purified by tp95).pdf"), height=8,width = 10)
pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
         cluster_rows = F,
         #cluster_cols = dendsort(hcs),
         cluster_cols = hcs,
         border_color = NA,
         annotation_col = annCol.target[colnames(plotdata),c(7,8,4,1:3)],
         annotation_colors = annColors.target,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 3,
         cellheight = 10,
         gaps_row = c(14,22),
         cutree_cols = 2,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

# comparsion of CD4 memory activated
tmp1 <- rownames(annCol.target[which(annCol.target$tune_Clust == "C1"),])
tmp2 <- rownames(annCol.target[which(annCol.target$tune_Clust == "C2"),])
tmp1 <- indata["T.cells.CD4.memory.activated",tmp1]
tmp2 <- indata["T.cells.CD4.memory.activated",tmp2]
wilcox.test(tmp1,tmp2)
t.test(tmp1,tmp2,alternative = "greater")

# survival analysis
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.target[rownames(Sinfo.target),"tune_Clust"],
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jco",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by target 114 (purified by tp95) data in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

#--------------------------------------------------------------------------------------------#
# Reanalysis TARGET 114 samples (purified by tumor purity) with as much annotation as we can #
tmp1 <- annCol.target[,c(1:4,8)]; tmp1$sample <- rownames(tmp1)
tmp2 <- annCol.cbio[,c(2:17,28,32)]; tmp2$sample <- rownames(tmp2)
tmp3 <- annCol.ccr39[,c(9:22)]; tmp3$sample <- rownames(tmp3)

tmp <- merge(tmp1,tmp2,by = "sample",all.x = T)
tmp <- merge(tmp,tmp3,by = "sample",all.x = T)
rownames(tmp) <- tmp$sample; tmp <- tmp[,-1]
tmp <- cbind.data.frame(tmp,Sinfo.target[rownames(tmp),c(13,16,19,20)])
write.table(tmp,file = file.path(data.path,"curated annotation for 114 target data combined with ccr and cbioportal.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

annCol.target114 <- read.table(file.path(data.path,"curated annotation for 114 target data combined with ccr and cbioportal.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
annCol.target114 <- annCol.target114[,setdiff(colnames(annCol.target114),c("CTNNB1","DROSHA","SIX1","WT1","DGCR8"))]
annCol.target114$Stage <- ifelse(Sinfo.target$Stage %in% c("I","II"),"I_II",">III")
annCol.target114$Age2 <- ifelse(annCol.target114$Age > 2*365,">2","<=2")
annCol.target114$Age3 <- ifelse(annCol.target114$Age > 3*365,">3","<=3")
annCol.target114$Agem <- ifelse(annCol.target114$Age > median(annCol.target114$Age),">4.4","<=4.4")

annColors.target114 <- list()
annColors.target114[["TP53_MUT"]] <- c("Yes" = "black","No" = "grey90")
annColors.target114[["TP53_LOSS"]] <- c("Yes" = darkblue,"No" = "grey90")
#annColors.target114[["TP53_ALTER"]] <- c("Yes" = "black","No" = "grey70")
annColors.target114[["HISTOLOGY"]] <- c("DAWT" = darkblue,"FHWT" = darkred)
annColors.target114[["Stage"]] <- c("I_II" = alpha(seagreen,0.3),">III" = alpha(seagreen,0.6))
#annColors.target114[["Tumor_purity"]] <- NMF:::ccRamp(c("white","black"),n=64)
#annColors.target114 <- append(annColors.target114,annColors.ccr39)
#annColors.target114 <- append(annColors.target114,annColors.cbio)
annColors.target114[["Agem"]] <- c(">4.4" = orange,"<=4.4" = alpha(orange,0.6))
#annColors.target114[["tune_Clust"]] <- c("C1" = jco[2],"C2" = jco[1],"C3" = seagreen)
annColors.target114[["tune_Clust"]] <- c("C1" = jco[2],"C2" = jco[1])
annColors.target114[["OS"]] <- c("Alive" = "grey90","Dead" = brown)
annColors.target114[["Gender"]] <- c("Female" = jco[2], "Male" = jco[1])

# tuning all target data 114
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- target.pure$es.pure[,com_sam]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D"); hcs.target114 <- hcs
group <- cutree(hcs,3)
annCol.target114$tune_Clust <- NA
annCol.target114[com_sam,"tune_Clust"] <- paste0("C",group[com_sam])
annCol.target114$tune_Clust <- ifelse(annCol.target114$tune_Clust == "C2","C1",
                                      ifelse(annCol.target114$tune_Clust == "C3", "C2","C3"))
annCol.target114$tune_Clust <- ifelse(annCol.target114$tune_Clust == "C1","C1","C2")
#annCol.target114$Tumor_purity <- as.numeric(est.target.backup[4,rownames(annCol.target)])
annCol.target114$TLS <- annTrackScale(TLS.target[rownames(annCol.target114)],halfwidth = 2)
annColors.target114[["TLS"]] <- blueyellow(64)
fisher.test(table(annCol.target114$tune_Clust,annCol.target114$HISTOLOGY))
fisher.test(table(annCol.target114$tune_Clust,annCol.target114$TP53_ALTER))
fisher.test(table(annCol.target114$tune_Clust,annCol.target114$TP53_MUT))

mut6 <- read.table(file.path(data.path,"cbioportal_6_mutations_target.tsv"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mut6[mut6 == ""] <- "Wild"
mut6[mut6 != "Wild"] <- "Mutated"
mut6[is.na(mut6)] <- "Wild"
colnames(mut6) <- paste0(colnames(mut6),"A")
tmp <- intersect(colnames(mut6),rownames(matchID))
mut6 <- mut6[,tmp]
colnames(mut6) <- matchID[tmp,"SRA_RUN"]
mut6 <- mut6[,rownames(annCol.target114)]
annCol.target114 <- cbind.data.frame(annCol.target114,t(mut6))
annColors.target114[["CTNNB1"]] <- c("Mutated" = "black", "Wild" = "grey90")
annColors.target114[["DROSHA"]] <- c("Mutated" = "black", "Wild" = "grey90")
annColors.target114[["SIX1"]] <- c("Mutated" = "black", "Wild" = "grey90")
annColors.target114[["SIX2"]] <- c("Mutated" = "black", "Wild" = "grey90")
annColors.target114[["WT1"]] <- c("Mutated" = "black", "Wild" = "grey90")
annColors.target114[["DGCR8"]] <- c("Mutated" = "black", "Wild" = "grey90")

p <- pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
         cluster_rows = F,
         cluster_cols = dendsort(hcs.target114),
         border_color = NA,
         annotation_col = annCol.target114[colnames(plotdata),c(45,4,41,37,44,39,2,1,47:52)],
         annotation_colors = annColors.target114,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 5,
         cellheight = 13,
         gaps_row = c(14,22),
         cutree_cols = 2,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
pdf(file.path(fig.path, "CCR immune signature heatmap of 114 wilms tumors in target cohort (purified by tp) 2 clusters.pdf"), height=10,width = 15)
draw(p, annotation_legend_side = "bottom")
invisible(dev.off())

indata <- target.wt[TLS,colnames(plotdata)]
plotdata <- standarize.fun(indata,halfwidth = 2)
p <- pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = dendsort(hcs.target114),
         border_color = NA,
         annotation_col = annCol.target114[colnames(plotdata),c(45,46)],
         annotation_colors = annColors.target114[c("tune_Clust","TLS")],
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         #color = viridisLite::inferno(64),
         #color = viridisLite::viridis(64),
         # color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 5,
         cellheight = 13,
         show_colnames = F,
         show_rownames = T,
         cutree_cols = 2,
         fontsize_col = 8)
pdf(file.path(fig.path, "TLS heatmap of wilms tumors in target cohort.pdf"), height=10,width = 15)
draw(p, annotation_legend_side = "bottom")
invisible(dev.off())

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- target.pure$es.pure[,com_sam]
outTab <- NULL
for(cell in rownames(indata)) {
  tmp1 <- as.numeric(indata[cell,RNA_Clust1.target])
  tmp2 <- as.numeric(indata[cell,RNA_Clust2.target])
  #tmp3 <- as.numeric(indata[cell,RNA_Clust3.target])
  
  #kruskal.p <- kruskal.test(list(tmp1,tmp2,tmp3))$p.value
  wilcox.p12 <- wilcox.test(tmp1,tmp2)$p.value
  #wilcox.p13 <- wilcox.test(tmp1,tmp3)$p.value
  #wilcox.p23 <- wilcox.test(tmp2,tmp3)$p.value
  #wilcox.p1.23 <- wilcox.test(tmp1,c(tmp2,tmp3))$p.value
  
  outTab <- rbind.data.frame(outTab,data.frame(cell = cell,
                                               avg.C1 = mean(tmp1),
                                               avg.C2 = mean(tmp2),
                                               #avg.C3 = mean(tmp3),
                                               #p.kruskal = kruskal.p,
                                               p12.wilcox = wilcox.p12,
                                               #p13.wilcox = wilcox.p13,
                                               #p23.wilcox = wilcox.p23,
                                               #p1.23.wilcox = wilcox.p1.23,
                                               stringsAsFactors = F),
                             stringsAsFactors = F)
}
#outTab$diff1.23 <- outTab$avg.C1-(outTab$avg.C2+outTab$avg.C3)
outTab$diff <- outTab$avg.C1-outTab$avg.C2
#outTab$FDR <- p.adjust(outTab$p1.23.wilcox,method = "BH")
outTab$FDR <- p.adjust(outTab$p12.wilcox,method = "BH")
write.table(outTab,file.path(res.path,"comparision of CCR immune signature between 2 target clusters.txt"),sep = "\t",row.names = F,quote = F)

# differential expression of mRNA
# C1 vs Others
tmp <- annCol.target114$tune_Clust; names(tmp) <- rownames(annCol.target114)
tmp <- na.omit(tmp)
tmp <- ifelse(tmp == "C1","C1","Others")
Groupinfo <- annCol.target114[names(tmp),]
colnames(Groupinfo)[40] <- "group"
Groupinfo$group <- ifelse(Groupinfo$group == "C1","C1","Others")
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.RNAclust.targetC1(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="TARGET_WT_mRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

# C2 vs Others
tmp <- annCol.target114$tune_Clust; names(tmp) <- rownames(annCol.target114)
tmp <- na.omit(tmp)
tmp <- ifelse(tmp == "C2","C2","Others")
Groupinfo <- annCol.target114[names(tmp),]
colnames(Groupinfo)[40] <- "group"
Groupinfo$group <- ifelse(Groupinfo$group == "C2","C2","Others")
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.RNAclust.targetC2(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="TARGET_WT_mRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

deres <- read.table(file.path(res.path,"TARGET_WT_mRNA_edgeR_test_result.C2_vs_Others.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- deres$logFC
names(geneList) <- deres$id
geneList <- sort(geneList,decreasing = T)
gsea_target_c2vsc1 <- GSEA(geneList = geneList,
                           TERM2GENE=MSigDB.hallmark,
                           pvalueCutoff = 0.25,
                           nPerm = 10000,
                           seed = T,
                           verbose=F)
res <- data.frame(gsea_target_c2vsc1)
write.table(as.data.frame(gsea_target_c2vsc1),file.path(res.path,"GSEA hallmark results for TARGET WT mRNA TP53 Mutated C2 vs Wild C1 DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

gsea_target_c2vsc1_hm.rt <- GSEA(geneList = geneList,
                           TERM2GENE=MSigDB.hm.rt,
                           pvalueCutoff = 0.25,
                           nPerm = 10000,
                           seed = T,
                           verbose=F)
res <- data.frame(gsea_target_c2vsc1_hm.rt)
write.table(as.data.frame(gsea_target_c2vsc1_hm.rt),file.path(res.path,"GSEA hallmark and reactome results for TARGET WT mRNA TP53 Mutated C2 vs Wild C1 DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

gseaplot2(x = gsea_target_c2vsc1_hm.rt,
          geneSetID = c("HALLMARK_P53_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
          pvalue_table = F,
          color = c(alpha(darkblue,0.3),alpha(darkblue,0.6),darkblue))
ggsave(filename = file.path(fig.path,"GSEA downregulated pathway in target cohort.pdf"), width = 5.5,height = 4)

gseaplot2(x = gsea_target_c2vsc1_hm.rt,
          geneSetID = c("REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA","REACTOME_HDACS_DEACETYLATE_HISTONES","REACTOME_HOMOLOGY_DIRECTED_REPAIR"),
          pvalue_table = F,
          color = c(alpha(darkred,0.3),alpha(darkred,0.6),darkred))
ggsave(filename = file.path(fig.path,"GSEA upregulated pathway in target cohort.pdf"), width = 5.5,height = 4)

gseaplot2(x = gsea_target_c2vsc1,
          geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          pvalue_table = F,
          color = seagreen)
ggsave(filename = file.path(fig.path,"GSEA inteferon gamma in target.pdf"), width = 5.5,height = 4)

gsea_target_c2vsc1_all <- GSEA(geneList = geneList,
                           TERM2GENE=MSigDB,
                           pvalueCutoff = 0.25,
                           nPerm = 10000,
                           seed = T,
                           verbose=F)
res <- data.frame(gsea_target_c2vsc1_all)
write.table(as.data.frame(gsea_target_c2vsc1_all),file.path(res.path,"GSEA all results for TARGET WT mRNA TP53 Mutated C2 vs Wild C1 DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

# survival analysis
tmp <- data.frame(OS.time = annCol.target114$OS.time/365,
                  OS = ifelse(annCol.target114$OS == "Dead",1,0),
                  tune_Clust = annCol.target114$tune_Clust2,
                  Age = annCol.target114$Age/365,
                  Stage = annCol.target114$Stage,
                  stringsAsFactors = F)
tmp$Stage <- factor(tmp$Stage,levels = c("I_II",">III"))
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1],seagreen),
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Years)",
                ylab = "Overall Survival (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort 2 clusters.pdf"),width = 5,height = 5)
p
invisible(dev.off())

# pairwise_survdiff(Surv(OS.time, OS)~ tune_Clust, data=tmp)
# C1      C2     
# C2 0.14606 -      
# C3 0.00049 0.22074

# univariate cox regression to see if the cluster is independent predictors
ucox1 <- coxph(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude) # 0.001
ucox2 <- coxph(Surv(OS.time, OS) ~ Age, data=tmp, na.action=na.exclude) # 0.289
ucox3 <- coxph(Surv(OS.time, OS) ~ Stage, data=tmp, na.action=na.exclude) # 0.826

# multivariate cox regression to see if the cluster is independent predictors
mcox <- coxph(Surv(OS.time, OS) ~ tune_Clust + Age + Stage, data=tmp, na.action=na.exclude)
# coef exp(coef) se(coef)      z       p
# tune_ClustC2  1.08366   2.95549  0.33288  3.255 0.00113
# Age          -0.06889   0.93343  0.06641 -1.037 0.29956
# Stage>III    -0.03462   0.96597  0.33490 -0.103 0.91766
# 
# Likelihood ratio test=11.15  on 3 df, p=0.01092
# n= 96, number of events= 37 

# fraction genome altered comparison
tmp1 <- annCol.target114[which(annCol.target114$tune_Clust == "C1"),"Fraction Genome Altered"]
tmp2 <- annCol.target114[which(annCol.target114$tune_Clust == "C2"),"Fraction Genome Altered"]
#tmp3 <- annCol.target114[which(annCol.target114$tune_Clust == "C3"),"Fraction Genome Altered"]

wilcox.test(tmp2,tmp1,alternative = "greater")
tmp <- data.frame(exp = c(tmp2,tmp1),
                  mut = rep(c("Desert","Infiltrated"),c(length(tmp2),length(tmp1))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for Fraction Genome Altered in cbioportal 101 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Fraction Genome Altered",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.7,"P = 0.047",cex = 1)
invisible(dev.off())

# comparsion of TMB
tmp1 <- annCol.target114[which(annCol.target114$tune_Clust == "C2"),"TMB"]
tmp2 <- annCol.target114[which(annCol.target114$tune_Clust == "C1"),"TMB"]

wilcox.test(tmp1,tmp2)

# comparsion of CD4 memory activated
tmp1 <- rownames(annCol.target114[which(annCol.target114$tune_Clust == "C1"),])
tmp2 <- rownames(annCol.target114[which(annCol.target114$tune_Clust == "C2"),])

tmp1 <- indata["T.cells.CD4.memory.activated",tmp1]
tmp2 <- indata["T.cells.CD4.memory.activated",tmp2]

t.test(tmp1,tmp2,alternative = "less")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Infiltrated","Desert"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for CD4 memory activated in cbioportal 101 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Activated memory CD4+ T cells",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,18,"P = 0.025",cex = 1)
invisible(dev.off())

# immune cell comparsion among three groups
outTab <- NULL
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- target.pure$es.pure[,com_sam]

RNA_Clust1.target <- rownames(annCol.target114[which(annCol.target114$tune_Clust == "C1"),])
RNA_Clust2.target <- rownames(annCol.target114[which(annCol.target114$tune_Clust == "C2"),])

for (cell in rownames(indata)) {
  tmp1 <- as.numeric(indata[cell,RNA_Clust1.target])
  tmp2 <- as.numeric(indata[cell,RNA_Clust2.target])

  wp <- wilcox.test(tmp1,tmp2)$p.value
  outTab <- rbind.data.frame(outTab,
                             data.frame(cell = cell,
                                        avg.clust1 = mean(tmp1),
                                        avg.clust2 = mean(tmp2),
                                        wilcox.p = round(wp,4),
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
  
}
write.csv(outTab,file.path(res.path,"comparsion of CCR immune signatures among 2 TARGET clusters.csv"),row.names = F,quote = F)

library(tidyr)
dd <- as.data.frame(t(indata))
dd$RNA_Clust <- annCol.target114[rownames(dd),"tune_Clust"]
dd$sample <- rownames(dd)
d2 <- gather(dd, cell, expr, 1:24)

pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(expr ~ RNA_Clust, data = subset(d2, cell == x))
  return(res$p.value)
})
pv <- data.frame(gene = d2$cell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

ggplot(d2, aes(cell, expr, fill=RNA_Clust)) + scale_fill_manual(values = c(jco[2],jco[1])) +
  geom_boxplot() + 
  geom_text(aes(gene, y=max(d2$expr) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("Enrichment score") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(fig.path,"boxplot for all CCR immune signatures among 2 clusters in 114 target cohort (purified by tp).pdf"),width = 10,height = 6)

# clinical features comparison
tmp <- annCol.target114
stabl <- CreateTableOne(vars=colnames(tmp),strata="tune_Clust",
                        data=tmp,factorVars=colnames(tmp)[c(1:4,6:21,24:37,39,41:44)])
print(stabl,showAllLevels = TRUE)
write.table(print(stabl,showAllLevels = TRUE),file.path(res.path,"statistic between CCR immune signatures 2 groups of 114 target wilms tumors.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# Gene ontology for Gaby wt cohort mRNA
fdr.cutoff <- 0.05
logfc.cutoff <- 2
deres <- read.table(file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.TP53_Mutated_vs_TP53_Wild.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
up.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC > logfc.cutoff),]
dn.gene <- deres[which(deres$FDR < fdr.cutoff & deres$logFC < -logfc.cutoff),]
tmp <- bitr(as.character(up.gene$id), fromType = 'SYMBOL', toType = 'ENSEMBL', OrgDb = 'org.Hs.eg.db')
ego_up <-     enrichGO(gene          = as.character(tmp$ENSEMBL),
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENSEMBL",
                       ont           = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff  = 0.25,
                       qvalueCutoff  = 0.25,
                       minGSSize     = 10,
                       readable      = T)
write.table(data.frame(summary(ego_up)),file.path(res.path,"Gaby WT TP53_mutated group upregualted genes GO BP.txt"),sep = "\t",row.names = F,quote = F)
ego_up.sim <- simplify(ego_up)
dotplot(ego_up.sim,showCategory = 30)

tmp <- bitr(as.character(dn.gene$id), fromType = 'SYMBOL', toType = 'ENSEMBL', OrgDb = 'org.Hs.eg.db')
ego_dn <-     enrichGO(gene          = as.character(tmp$ENSEMBL),
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENSEMBL",
                       ont           = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff  = 0.25,
                       qvalueCutoff  = 0.25,
                       minGSSize     = 10,
                       readable      = T)
write.table(data.frame(summary(ego_dn)),file.path(res.path,"Gaby WT TP53_wild group upregualted genes GO BP.txt"),sep = "\t",row.names = F,quote = F)

#----------------------------#
# immunity in normal samples #
TPM.NORMAL <- TPM[,c(nor.gaby.sam,nor.target.sam)]
indata <- log2(TPM.NORMAL[Mids,] + 1)
indata$symbol <- Ginfo[rownames(indata),"genename"]
indata <- apply(indata[,setdiff(colnames(indata), "symbol")], 2, function(x) tapply(x, INDEX=factor(indata$symbol), FUN=median, na.rm=TRUE))
write.table(indata,file = file.path(res.path,"kidney_normal_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

filterCommonGenes(input.f=file.path(res.path, "kidney_normal_log2TPM_hugo.txt") , output.f=file.path(res.path,"kidney_normal_log2TPM_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"kidney_normal_log2TPM_hugo_ESTIMATE.txt"), file.path(res.path,"kidney_normal_log2TPM_hugo_estimate_score.txt"), platform="affymetrix")
est.normal <- read.table(file = file.path(res.path,"kidney_normal_log2TPM_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.normal) <- est.normal[,2]; colnames(est.normal) <- est.normal[1,]; est.normal <- est.normal[-1,c(-1,-2)];
est.normal <- sapply(est.normal, as.numeric); rownames(est.normal) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.normal.backup = as.data.frame(est.normal); colnames(est.normal.backup) <- colnames(indata)
est.normal <- annTrackScale(indata = est.normal, halfwidth = 2, poolsd = F); est.normal <- as.data.frame(t(est.normal))
rownames(est.normal) <- colnames(indata)

# submap
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

gaby.wt <- read.table(file.path(res.path,"Wilms_tumor_Gaby_log2TPM_hugo.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
var.gaby <- apply(gaby.wt, 1, mad)
sel_gene.gaby <- names(var.gaby[var.gaby > quantile(var.gaby)[4]])
gaby.wt.info <- annCol.wt.gaby
gaby.wt.info$rank <- ifelse(gaby.wt.info$RNA_Clust == "TP53_Wild",1,2)

target.wt <- read.table(file.path(res.path,"Wilms_tumor_TARGET_log2TPM_hugo.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
target.wt.info <- as.data.frame(annCol.target114[which(annCol.target114$tune_Clust %in% c("C1","C2","C3")),])
target.wt.info$rank <- ifelse(target.wt.info$tune_Clust == "C1",1,2)

target.wt <- target.wt[,rownames(target.wt.info)]
var.target <- apply(target.wt, 1, mad)
sel_gene.target <- names(var.target[var.target > quantile(var.target)[4]])
GENELIST <- intersect(rownames(gaby.wt),rownames(target.wt))

sam_info <- gaby.wt.info
in_gct <- gaby.wt[GENELIST,rownames(gaby.wt.info)]

gct_file <- file.path(res.path,"gaby.wt.for.SubMap.gct")
cls_file <- file.path(res.path,"gaby.wt.for.SubMap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

sam_info <- target.wt.info
in_gct <- target.wt[GENELIST,rownames(target.wt.info)]

gct_file <- file.path(res.path,"target.wt.for.SubMap.gct")
cls_file <- file.path(res.path,"target.wt.for.SubMap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

#-------------------#
# Mutation analysis #
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(deconstructSigs)
library(ClassDiscovery)
library(shape)

# 1. processing mutation data
exome.dat <- read.table(file.path(data.path,"Gaby WT exome data 168.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = NULL,header = T)
maf <- exome.dat %>% 
  gather( sample, sample_status, colnames(exome.dat)[8:19]) %>%
  filter( sample_status == 1 )

library(maftools)
tmp <- maf[,c(2:10,12)]
tmp$End_Position <- tmp$Start + nchar(tmp$AlleleA)
tmp$Variant_Type <- ifelse(tmp$Status %in% c("frameshift deletion","nonframeshift deletion"),"DEL",
                           ifelse(tmp$Status %in% c("frameshift insertion","nonframeshift insertion"),"INS",
                                  "SNP"))
colnames(tmp) <- c("Hugo_Symbol","Chromosome","Start_Position","cDNA_shift","peptide_shift","SNV_Diff","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Tumor_Sample_Barcode","End_Position","Variant_Type")

tmp$Variant_Classification <- gsub("nonframeshift insertion","In_Frame_Ins",tmp$Variant_Classification)
tmp$Variant_Classification <- gsub("nonframeshift deletion","In_Frame_Del",tmp$Variant_Classification)
tmp$Variant_Classification <- gsub("frameshift insertion","Frame_Shift_Ins",tmp$Variant_Classification)
tmp$Variant_Classification <- gsub("frameshift deletion","Frame_Shift_Del",tmp$Variant_Classification)
tmp$Variant_Classification <- gsub("nonsynonymous SNV","Missense_Mutation",tmp$Variant_Classification)
tmp$Variant_Classification <- gsub("splicing","Splice_Site",tmp$Variant_Classification)
tmp$Variant_Classification <- gsub("stopgain SNV","Nonsense_Mutation",tmp$Variant_Classification)
write.csv(tmp,file.path(res.path,"Gaby WT exome data long format 168.csv"),row.names = F,quote = F)

# 2. perform deconstructSigs
maf <- tmp
maf$Chr <- paste0("chr",maf$Chromosome)
unique(maf$Variant_Classification)
# [1] "Nonsense_Mutation" "Missense_Mutation" "Splice_Site"       "Frame_Shift_Del"  
# [5] "In_Frame_Ins"      "Frame_Shift_Ins"   "In_Frame_Del" 
sigs.input.sbs <- mut.to.sigs.input(mut.ref = maf, 
                                    sample.id = "Tumor_Sample_Barcode", 
                                    chr = "Chr", 
                                    pos = "Start_Position", 
                                    ref = "Reference_Allele", 
                                    alt = "Tumor_Seq_Allele2",
                                    sig.type = "SBS")
write.csv(sigs.input.sbs,file.path(res.path,"mutation.sig.input.bydeconstructSigs.csv"),row.names = T,quote = F)

cut.off <- 0.06
mut.wt.cosmic2013 <- mut.wt.cosmic2019 <- data.frame()
sigs.out.cosmic2013.list <- sigs.out.cosmic2019.list <- list()
for (sample in rownames(sigs.input.sbs)) {
  tmp <- whichSignatures(tumor.ref = sigs.input.sbs, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome',
                         signature.cutoff = cut.off)
  
  sigs.out.cosmic2013.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt.cosmic2013 <- rbind.data.frame(mut.wt.cosmic2013,tmp)
  
  tmp <- whichSignatures(tumor.ref = sigs.input.sbs, 
                         signatures.ref = signatures.exome.cosmic.v3.may2019, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome',
                         signature.cutoff = cut.off)
  
  sigs.out.cosmic2019.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt.cosmic2019 <- rbind.data.frame(mut.wt.cosmic2019,tmp)
}
write.csv(mut.wt.cosmic2013,file.path(res.path,"mutation.snp.signature.weightMatrix.cosmic2013.csv"),row.names = T,quote = F)
write.csv(mut.wt.cosmic2019,file.path(res.path,"mutation.snp.signature.weightMatrix.cosmic2019.csv"),row.names = T,quote = F)

#---------------------#
# COSMIC2013 analysis #
MB <- 1
mut.wt.cosmic2013.trim <- mut.wt.cosmic2013[,1:30]
mut.wt.cosmic2013.trim <- mut.wt.cosmic2013.trim[,colSums(mut.wt.cosmic2013.trim) > 0] #17 signature for snp
mut.wt.cosmic2013.trim.backup <- mut.wt.cosmic2013.trim

mut.wt.cosmic2013.trim <- mut.wt.cosmic2013.trim * MB
mut.wt.cosmic2013.trim$Subject <- rownames(mut.wt.cosmic2013.trim)

# 2. signature2013 comparison regarding systematic therapy
comp.cosmic2013 <- NULL
for (i in 1:ncol(mut.wt.cosmic2013.trim.backup)) {
  tmp1 <- mut.wt.cosmic2013.trim.backup[RNA_Clust1,i]
  tmp2 <- mut.wt.cosmic2013.trim.backup[RNA_Clust2,i]
  
  comp.cosmic2013 <- rbind.data.frame(comp.cosmic2013,data.frame(avg.RNAClust1=mean(tmp1),
                                                                 avg.RNAClust2=mean(tmp2),
                                                                 pvalue=wilcox.test(tmp1,tmp2)$p.value,
                                                                 stringsAsFactors = F,
                                                                 row.names = colnames(mut.wt.cosmic2013.trim.backup)[i]),
                                      stringsAsFactors = F)
}

# 3. NMF
keepsig <- colnames(mut.wt.cosmic2013.trim.backup)[colSums(mut.wt.cosmic2013.trim.backup > 0.10) >= 3]

nmf.input <- t(mut.wt.cosmic2013.trim.backup[,keepsig])
mut.nmf <- nmf(nmf.input + 0.0001,3,seed = 2020525,method = "nsNMF",nrun = 50)

group <- predict(mut.nmf)
annCol.mutsig <- gaby.wt.genomic.alter[,1:5]
annCol.mutsig <- cbind.data.frame(annCol.mutsig,gaby.wt.mlpa[rownames(annCol.mutsig),7:20])
annCol.mutsig$OS <- gaby.wt.mlpa[rownames(annCol.mutsig),"Status"]
annCol.mutsig$COSMIC2015 <- paste0("C",group[rownames(annCol.mutsig)])
annCol.mutsig$RNA_Clust <- "N/A"
annCol.mutsig[RNA_Clust1,"RNA_Clust"] <- "C1"
annCol.mutsig[RNA_Clust2,"RNA_Clust"] <- "C2"
annColors.mutsig <- annColors.wt.gaby.genomic.alter
annColors.mutsig[["COSMIC2015"]] <- c("C1" = red, "C2" = blue, "C3" = seagreen)
annColors.mutsig[["RNA_Clust"]] <- c("C1" = jco[2], "C2" = jco[1],"N/A" = "grey90")

index <- extractFeatures(mut.nmf,"max")
sig.order <- unlist(index)
sample.order <- names(group[order(group)])

pdf(file.path(fig.path,"heatmap of NMF clusters in mutation signatures2015 darkcolor.pdf"),width = 8,height = 40)
pheatmap(
  nmf.input[,sample.order],
  border_color = "black",
  cluster_cols = F,
  cluster_rows = F,
  show_rownames = T,
  show_colnames = T,
  annotation_col = annCol.mutsig[sample.order,],
  annotation_colors = annColors.mutsig,
  color = viridisLite::inferno(64),
  gaps_col = c(5,8),
  fontsize = 10,
  cellwidth = 12,
  cellheight = 15)
invisible(dev.off())

# 4. oncoprint with wilms specific mutations
library(ComplexHeatmap)
onco_show_gene <- c("DROSHA","TP53",
                    "WT1","CTNNB1","AMER1","IGF2","MYCN","SIX1","SIX2",
                    "SMARCA4","MLLT1","BCORL1","COL6A3","NF1","BCOR","NONO",
                    "ARID1A","MAP3K4","MAX","ASXL1","BRD7","FGFR1","HDAC4",
                    "CHD4","ACTB","DICER1","DGCR8","XPO5","DIS3L2","TARBP2")
onco_show_gene <- intersect(onco_show_gene,maf$Hugo_Symbol)
mut.matrix <- matrix("", # initialize a empty matrix
                     ncol = 12,
                     nrow = length(onco_show_gene),
                     dimnames = list(onco_show_gene,unique(maf$Tumor_Sample_Barcode)))

for (sam in colnames(mut.matrix)) {
  cat(sam,"\n")
  for (gene in onco_show_gene) {
    tmp <- maf[which(maf$Tumor_Sample_Barcode == sam),]
    
    if(is.element(gene,tmp$Hugo_Symbol)) {
      tmp <- tmp[which(tmp$Hugo_Symbol == gene),]
      # put the largest number of mutation type into the matrix
      mut.matrix[gene,sam] <- names(which.max(table(tmp$Variant_Classification))) 
    } else {next()}
  }
}
mut.matrix["TP53","PEDrna_13T"] <- "Germline_Mutation"
# calculate variant type for print
type <- c()
for (i in 1:nrow(mut.matrix)) {
  tmp <- as.character(mut.matrix[i,])
  type <- unique(c(type,tmp))
}
cat(paste0("You have a total of ",length(setdiff(type,""))," types of mutation!\nPlease check as the following:\n"))
print(setdiff(type,""))# this is all the mutation type now you have and you have to set the color one by one

# set colors for oncoprint
mycol <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')[c(2,4,6,8)]
lightgrey <- "#dcddde" # background color was set to grey
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = lightgrey, col = lightgrey))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[1], col = mycol[1])) # e.g., the first color in mycol vector is for Missense_Mutation
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[2], col = mycol[2])) 
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[3], col = mycol[3])) 
  },
  Germline_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[4], col = mycol[4])) 
  }
)

# set color again for mutation annotation which should exactely the same with above
col = c("Missense_Mutation" = mycol[1], 
        "Nonsense_Mutation" = mycol[2], 
        "Frame_Shift_Del" = mycol[3], 
        "Germline_Mutation" = mycol[4])
my_ann <- annCol.mutsig[colnames(mut.matrix),c(8,11,13,17)]
my_annotation = HeatmapAnnotation(df = data.frame(my_ann,check.names = F), # annotation data frame
                                  col = list("MYCN STATUS" = c("NORMAL" = "white", "GAIN" = darkred),
                                             "FBXW7 STATUS" = c("NORMAL" = "white", "LOSS" = darkblue),
                                             "WT1 STATUS" = c("NORMAL" = "white", "GAIN" = darkred),
                                             "TP53 STATUS" = c("NORMAL" = "white", "LOSS" = darkblue)))

p <- oncoPrint(mut.matrix[,rownames(my_ann)], # this is the detailed mutation matrix
               alter_fun = alter_fun,  # this is the alteration function we already set
               col = col, # this is the color for mutation
               bottom_annotation = my_annotation, # this is the annotation bar which was put bottom
               # column_order = rownames(my_ann), # sort your sample as your subtype order
               # row_order = mut.order, # sort mutation as mutation frequency
               show_pct = T, #show percentage in left
               column_title = "", # no title shown
               show_heatmap_legend=T, # show legend in the oncoprint
               top_annotation = NULL,
               show_column_names = T,
               # column_split = my_ann$Subtype,
               # some detailed size below and you may not have to change it
               column_title_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8)) 

# output this oncoprint to pdf in the local working path
pdf(file.path(fig.path,"oncoprint of gaby wilms with specific mutations.pdf"),width = 7,height = 3.5)
draw(p, annotation_legend_side  = "right")
invisible(dev.off())

# mutual exclusivity of TP53 and DROSHA
fisher.test(matrix(c(1,3,5,3),nrow = 2,byrow = T))

tmp <- read.table(file.path(data.path,"mutual exclusivity input.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp <- t(tmp)
res <- NULL
for(i in 1:ncol(tmp)){
  for(j in i:ncol(tmp)){
    tab <- table(tmp[,i], tmp[,j])
    if(i!=j){
      f <- fisher.test(tab) 
      # deal with zero in 2x2 cell
      if(f$estimate == 0) {
        f <- fisher.test(tab,alternative = "less")
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
      } else {
        if(is.infinite(f$estimate)) {
          f <- fisher.test(tab + 1)
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab + 1,alternative = "greater")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
          } else{
            f <- fisher.test(tab + 1,alternative = "less")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
          }
          
        } else {
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab,alternative = "greater")
          } else{f <- fisher.test(tab,alternative = "less")}
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
        }   
      }
    }
  }
}
# some formatting
res <- as.data.frame(res)
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
res$geneA <- factor(res$geneA,levels=colnames(tmp))
res$geneB <- factor(res$geneB,levels=colnames(tmp))
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$log2OR <- log2(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))
# plot with ggplot 2
write.table(res,file.path(res.path,"Mutual_exclusivity.txt"),sep = "\t",row.names = F)

p <- ggplot(res, aes(geneA, geneB)) + geom_tile(aes(fill = log2OR),colour = "white") + scale_fill_gradient2(low = "darkblue",mid = lightgrey,high = "darkred",midpoint=0) + 
  geom_text(aes(label=stars), color="white", size=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank())
p
ggsave(file.path(fig.path,"Mutual_exclusivity.pdf"),width = 4,height = 3.5)

# calculate TMB
require(dplyr)
mutect.dataframe <- function(x){
  somatic_sum <- x %>% group_by(Tumor_Sample_Barcode) %>% summarise(TMB = n())
}
variants_per_sample <- as.data.frame(mutect.dataframe(maf))
rownames(variants_per_sample) <- variants_per_sample$Tumor_Sample_Barcode
head(variants_per_sample)

tmp1 <- as.numeric(variants_per_sample[RNA_Clust1,"TMB"])
tmp2 <- as.numeric(variants_per_sample[RNA_Clust2,"TMB"])
wilcox.test(tmp1,tmp2) #0.9163

diffuse.sam <- rownames(gaby.wt.mlpa[which(gaby.wt.mlpa$Histology_type == "ANA diffuse"),])
focal.sam <- rownames(gaby.wt.mlpa[which(gaby.wt.mlpa$Histology_type == "Focal anaplasia"),])
tmp1 <- as.numeric(variants_per_sample[diffuse.sam,"TMB"])
tmp2 <- as.numeric(variants_per_sample[focal.sam,"TMB"])
wilcox.test(tmp1,tmp2) #0.9258

#----------------------#
# copy number analysis #

# load package
library(ggplot2)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationHub)
library(data.table)
# library(knitr)
library(biscuiteer)
library(copynumber)
library(Rmisc)

Sys.setenv(LANGUAGE = "en") # Display error message in English
options(stringsAsFactors = FALSE) # Forbid chr convert to factor
options(warn = -1)

# TARGET
seg <- read.table(file.path(data.path,"wt_target_2018_pub/data_cna_hg19.seg"),sep="\t",header=T,stringsAsFactors = F)
colnames(seg)<- c("ID","Chr","Start","End","Probes","Log2Ratio")
head(seg)

seg$ID <- paste0(seg$ID,"A")
com_sam <- intersect(unique(seg$ID),rownames(matchID))
seg <- seg[which(seg$ID %in% com_sam),]
seg$SRA_RUN <- matchID[match(seg$ID,rownames(matchID)),"SRA_RUN"]

tmp1 <- seg[which(seg$SRA_RUN %in% RNA_Clust1.target),] # 67
tmp2 <- seg[which(seg$SRA_RUN %in% RNA_Clust2.target),] # 27
#tmp3 <- seg[which(seg$SRA_RUN %in% RNA_Clust3.target),] # 15

tmp1$ID <- tmp1$SRA_RUN
tmp2$ID <- tmp2$SRA_RUN
#tmp3$ID <- tmp3$SRA_RUN

seg <- rbind.data.frame(tmp1,tmp2)
seg <- seg[,1:6]

# extract hg19 chromosome info
data(hg19.chromArm, package="biscuiteer")
arm <- rep("n",dim(seg)[1])
chromosomes <- rep(0,dim(seg)[1])
chromArm <- as.data.frame(hg19.chromArm[1:48,])

chromArm$seqnames2 <- rownames(chromArm)
chromArm$chr.length.sum <- cumsum(as.numeric(chromArm$width))
chromArm$chr.length.cumsum <- c(0, chromArm$chr.length.sum[-nrow(chromArm)])
chromArm$middle.chr <- round(diff(c(0, chromArm$chr.length.sum)) /2)
chromArm$middle.chr.genome <- chromArm$middle.chr + chromArm$chr.length.cumsum
head(chromArm)

chromArm$seqnames2 <- rownames(chromArm)
chromArm$chr.length.sum <- cumsum(as.numeric(chromArm$width))
chromArm$chr.length.cumsum <- c(0, chromArm$chr.length.sum[-nrow(chromArm)])
chromArm$middle.chr <- round(diff(c(0, chromArm$chr.length.sum)) /2)
chromArm$middle.chr.genome <- chromArm$middle.chr + chromArm$chr.length.cumsum
head(chromArm)

# convert segment input for copynumber package 
seg2 <- seg
seg2$Chr <- paste0("chr",seg2$Chr)
arm <- rep("n",dim(seg2)[1])
chromosomes <- rep(0,dim(seg2)[1])
chromArm <- as.data.frame(hg19.chromArm[1:48,])
chrom <- c(1:22,"X","Y")
names(chrom)<- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
for (i in names(chrom)){
  chromosomes[seg2$Chr== paste0("chr",i)] = chrom[i]
  nums=chromArm[chromArm$seqnames== paste0("chr",chrom[i]),"end"][1]
  arm[seg2$Chr== paste0("chr",i) & (seg2$Start < (nums+1) & seg2$End < (nums+1))] = "p"
  arm[seg2$Chr== paste0("chr",i) & !(seg2$Start < (nums+1) & seg2$End < (nums+1))] = "q"
}
seg2$arm <- arm
seg2$chromosome <- chromosomes
seg2 <- seg2[,c(1,2,7,3,4,5,6)]
colnames(seg2)<- c("sampleID","chrom","arm","start.pos","end.pos","n.probes","mean")
seg2$chrom <- as.character(lapply(seg2$chrom, function(x) {strsplit(x,"r")[[1]][2]}))
seg2[seg2$chrom == "X","chrom"]= "23"
seg2[seg2$chrom == "Y","chrom"]= "24"
seg2=seg2[!(seg2$chrom == "24" | seg2$chrom == "23"),]
length(unique(seg2$sampleID))
head(seg2)

# Frequence
pdf(file.path(fig.path,"target copy number percentage of samples.pdf"),width = 10,height = 3)
plotFreq(segments=seg2,thres.gain=0.2,thres.loss=-0.2)
invisible(dev.off())

# heatmap
pdf(file.path(fig.path,"target copy number segmation ordered by CCR immune cluster.pdf"),width = 13,height = 12)
plotAberration(segments=seg2,thres.gain=0.2,mar = c(3,5,3,1))
invisible(dev.off())

# GABY
seg <- read.table(file.path(data.path,"aWT-CNVSum-GM-01-Seg.txt"),sep="\t",header=T,stringsAsFactors = F)
colnames(seg)<- c("ID","Chr","Start","End","Probes","Log2Ratio")
head(seg)
seg$ID <- gsub("PED","PEDrna",seg$ID)

tmp1 <- seg[which(seg$ID %in% RNA_Clust1),] # 5
tmp2 <- seg[which(seg$ID %in% RNA_Clust2),] # 5

seg <- rbind.data.frame(tmp1,tmp2)

# extract hg19 chromosome info
data(hg19.chromArm, package="biscuiteer")
arm <- rep("n",dim(seg)[1])
chromosomes <- rep(0,dim(seg)[1])
chromArm <- as.data.frame(hg19.chromArm[1:48,])

chromArm$seqnames2 <- rownames(chromArm)
chromArm$chr.length.sum <- cumsum(as.numeric(chromArm$width))
chromArm$chr.length.cumsum <- c(0, chromArm$chr.length.sum[-nrow(chromArm)])
chromArm$middle.chr <- round(diff(c(0, chromArm$chr.length.sum)) /2)
chromArm$middle.chr.genome <- chromArm$middle.chr + chromArm$chr.length.cumsum
head(chromArm)

chromArm$seqnames2 <- rownames(chromArm)
chromArm$chr.length.sum <- cumsum(as.numeric(chromArm$width))
chromArm$chr.length.cumsum <- c(0, chromArm$chr.length.sum[-nrow(chromArm)])
chromArm$middle.chr <- round(diff(c(0, chromArm$chr.length.sum)) /2)
chromArm$middle.chr.genome <- chromArm$middle.chr + chromArm$chr.length.cumsum
head(chromArm)

# convert segment input for copynumber package 
seg2 <- seg
seg2$Chr <- paste0("chr",seg2$Chr)
arm <- rep("n",dim(seg2)[1])
chromosomes <- rep(0,dim(seg2)[1])
chromArm <- as.data.frame(hg19.chromArm[1:48,])
chrom <- c(1:22,"X","Y")
names(chrom)<- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
for (i in names(chrom)){
  chromosomes[seg2$Chr== paste0("chr",i)] = chrom[i]
  nums=chromArm[chromArm$seqnames== paste0("chr",chrom[i]),"end"][1]
  arm[seg2$Chr== paste0("chr",i) & (seg2$Start < (nums+1) & seg2$End < (nums+1))] = "p"
  arm[seg2$Chr== paste0("chr",i) & !(seg2$Start < (nums+1) & seg2$End < (nums+1))] = "q"
}
seg2$arm <- arm
seg2$chromosome <- chromosomes
seg2 <- seg2[,c(1,2,7,3,4,5,6)]
colnames(seg2)<- c("sampleID","chrom","arm","start.pos","end.pos","n.probes","mean")
seg2$chrom <- as.character(lapply(seg2$chrom, function(x) {strsplit(x,"r")[[1]][2]}))
seg2[seg2$chrom == "X","chrom"]= "23"
seg2[seg2$chrom == "Y","chrom"]= "24"
seg2=seg2[!(seg2$chrom == "24" | seg2$chrom == "23"),]
length(unique(seg2$sampleID))
head(seg2)

# Frequence
pdf(file.path(fig.path,"Gaby copy number percentage of samples.pdf"),width = 10,height = 3)
plotFreq(segments=seg2,thres.gain=0.2,thres.loss=-0.2)
invisible(dev.off())

# heatmap
pdf(file.path(fig.path,"Gaby copy number segmation ordered by CCR immune cluster.pdf"),width = 13,height = 3)
plotAberration(segments=seg2,thres.gain=0.2,mar = c(3,5,3,1))
invisible(dev.off())

#--------------------------------#
# FGA calculation in GABY cohort #
library(ggplot2)
library(magrittr)
library(patchwork)
library(Rmisc)
library(Cairo)
library(stringr)
Segment <- read.table(file.path(data.path,"aWT-CNVSum-GM-01-Seg.txt"),sep="\t",header=T,stringsAsFactors = F)
colnames(Segment)<- c("sample","chr","start","end","probes","value")

Segment$bases=Segment$end-Segment$start

data=data.frame()
for (i in 1:length(table(Segment$sample))) {
  tmp=Segment[Segment$sample==names(table(Segment$sample))[i],]
  FGA=sum(tmp[abs(tmp$value)>0.2,"bases"])/ sum(tmp[,"bases"])
  FGG=sum(tmp[tmp$value<(-0.2),"bases"])/sum(tmp[,"bases"])
  FGL=sum(tmp[tmp$value>0.2,"bases"])/sum(tmp[,"bases"])
  tmp=data.frame(Patient=names(table(Segment$sample))[i],FGA=FGA,FGG=FGG,FGL=FGL)
  data=rbind(data,tmp)
}
data$Patient <- gsub("PED","PEDrna",data$Patient)
fga.gaby <- data
rownames(fga.gaby) <- fga.gaby$Patient
rm(data)

tmp1 <- fga.gaby[RNA_Clust1,"FGA"]
tmp2 <- fga.gaby[RNA_Clust2,"FGA"]
wilcox.test(tmp1,tmp2,alternative = "less") # 0.0754

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Infiltrated","Desert"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for Fraction Genome Altered in Gaby WT regarding CCR immune signature group.pdf"),width = 1.8,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Fraction Genome Altered",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.075",cex = 1)
invisible(dev.off())

tmp1 <- fga.gaby[RNA_Clust1,"FGG"]
tmp2 <- fga.gaby[RNA_Clust2,"FGG"]
wilcox.test(tmp1,tmp2,alternative = "less")

tmp1 <- fga.gaby[RNA_Clust1,"FGL"]
tmp2 <- fga.gaby[RNA_Clust2,"FGL"]
wilcox.test(tmp1,tmp2,alternative = "less")

#-----------------#
# GISTIC analysis #

# Gaby cohort
Segment <- read.table(file.path(data.path,"aWT-CNVSum-GM-01-Seg.txt"),sep="\t",header=T,stringsAsFactors = F)
colnames(Segment)<- c("Sample","Chrom","Start","Stop","Mark","Seg.CN")
Segment$Sample <- gsub("PED","PEDrna",Segment$Sample)
write.table(Segment,file.path(res.path,"gistic input of gaby wilms cohort.txt"),sep = "\t",row.names = F,quote = F)

marker <- Segment[,1:4]
head(marker)

# create marker file for GISTIC2.0
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chrom"],2))
  c <- c(c,marker[i,"Start"],marker[i,"Stop"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"Gaby_wilms_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

# load gistic results with qvalue 0.25, cutoff 0.2, confident interval 95
gistic.gaby <- read.table(file.path(data.path,"236403/Gaby_wilms_0.2_95.all_thresholded.by_genes.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cnacomp.gaby <- NULL
for (i in 1:nrow(gistic.gaby)) {
  display.progress(index = i,totalN = nrow(gistic.gaby))
  
  # loss
  tmp <- gistic.gaby[i,]
  g <- rownames(tmp)
  loci <- tmp$Cytoband
  tmp <- data.frame(cna = as.numeric(tmp[,c(RNA_Clust1,RNA_Clust2)]),
                    group = rep(c("C1","C2"),c(5,5)),
                    stringsAsFactors = F)
  loss <- amp <- tmp
  loss$cna <- ifelse(loss$cna < 0,"LOSS","Others")
  amp$cna <- ifelse(amp$cna > 0,"AMP","Others")
  
  loss.dt <- as.data.frame.array(table(loss$cna,loss$group))
  amp.dt <- as.data.frame.array(table(amp$cna,amp$group))
  
  if(!is.element("AMP",rownames(amp.dt))) {
    amp.pct <- c(0,0)
    amp.p <- NA
  } else {
    amp.pct <- as.numeric(amp.dt[1,]/colSums(amp.dt))
    amp.p <- fisher.test(amp.dt)$p.value
  }
  
  if(!is.element("LOSS",rownames(loss.dt))) {
    loss.pct <- c(0,0)
    loss.p <- NA
  } else {
    loss.pct <- as.numeric(loss.dt[1,]/colSums(loss.dt))
    loss.p <- fisher.test(loss.dt)$p.value
  }
  
  cnacomp.gaby <- rbind.data.frame(cnacomp.gaby,
                                   data.frame(gene = g,
                                              loci = loci,
                                              loss.C1 = loss.pct[1],
                                              loss.C2 = loss.pct[2],
                                              p.loss = loss.p,
                                              amp.C1 = amp.pct[1],
                                              amp.C2 = amp.pct[2],
                                              p.amp = amp.p,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
  
}
rownames(cnacomp.gaby) <- cnacomp.gaby$gene
write.table(cnacomp.gaby,file.path(res.path,"comparison of gene level cna in gaby cohort between two immune clusters.txt"),sep = "\t",row.names = F,quote = F)

# copy number heatmap of microenvironment factors
micro.gene <- read.table(file.path(data.path,"microenviroment factors.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
rownames(micro.gene) <- micro.gene$Symbol
com_gene <- intersect(rownames(micro.gene),rownames(gaby.wt))
com_gene <- setdiff(com_gene,"HLA-DQB2")
com_cna <- intersect(com_gene,cnacomp.gaby$gene)
annRow.cna <- annRow[com_cna,,drop = F]

plotdata <- cnacomp.gaby[com_cna,c("loss.C1","loss.C2","amp.C1","amp.C2")]
blank <- "    " 
p.value <- cnacomp.gaby[com_cna,"p.loss"]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*","ns"))))
p.label <- formatC(p.value, 
                   format = "e",
                   digits = 2)
library(stringr)
add.label <- str_pad(paste0(rownames(plotdata)," ",sig.label),
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")
pdf(file.path(fig.path,"cna loss percentage of microenviroment factors in Gaby clusters.pdf"),width = 8,height = 22)
pheatmap(plotdata[,1:2],
         cluster_rows = F,cluster_cols = F,
         border_color = NA,
         color = bluered(64),
         annotation_row = annRow.cna[rownames(plotdata),,drop = F],
         annotation_colors = annColors.micro,
         gaps_row = c(46,87,93,116,120,126,128,142),
         #gaps_col = 2,
         labels_row = add.label)
invisible(dev.off())

plotdata <- cnacomp.gaby[com_cna,c("loss.C1","loss.C2","amp.C1","amp.C2")]
blank <- "    "
p.value <- cnacomp.gaby[com_cna,"p.amp"]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*","ns"))))
p.label <- formatC(p.value, 
                   format = "e",
                   digits = 2)
library(stringr)
add.label <- str_pad(paste0(rownames(plotdata)," ",sig.label),
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")
pdf(file.path(fig.path,"cna gain percentage of microenviroment factors in Gaby clusters.pdf"),width = 8,height = 22)
pheatmap(plotdata[,3:4],
         cluster_rows = F,cluster_cols = F,
         border_color = NA,
         color = bluered(64),
         annotation_row = annRow.cna[rownames(plotdata),,drop = F],
         annotation_colors = annColors.micro,
         gaps_row = c(46,87,93,116,120,126,128,142),
         #gaps_col = 2,
         labels_row = add.label)
invisible(dev.off())

# correlation between 4 CNV region and gene expression
loci <- c("11q14.3","16q22.1","16q24.3","21q11.2")
gistic.gaby.loci <- gistic.gaby[which(gistic.gaby$Cytoband %in% loci),]
gistic.gaby.loci <- gistic.gaby.loci[intersect(rownames(gistic.gaby.loci),rownames(gaby.wt)),]

cnacomp.gaby <- NULL
for (i in 1:nrow(gistic.gaby.loci)) {
  display.progress(index = i,totalN = nrow(gistic.gaby))
  
  tmp <- gistic.gaby.loci[i,]
  g <- rownames(tmp)
  loci <- tmp$Cytoband
  tmp <- data.frame(cna = as.numeric(tmp[,c(RNA_Clust1,RNA_Clust2)]),
                    exp = as.numeric(gaby.wt[g,c(RNA_Clust1,RNA_Clust2)]),
                    stringsAsFactors = F)
  
  if(loci %in% c("11q14.3","16q22.1","16q24.3")) {
    class <- "loss"
    tmp$class <- ifelse(tmp$cna < 0,"loss","others")
    tmp1 <- tmp[which(tmp$class == "loss"),"exp"]
    tmp2 <- tmp[which(tmp$class == "others"),"exp"]
    wt <- wilcox.test(tmp1,tmp2,alternative = "less")
    cnacomp.gaby <- rbind.data.frame(cnacomp.gaby,
                                     data.frame(gene = g,
                                                loci = loci,
                                                category = class,
                                                avg.cna = mean(tmp1),
                                                avg.others = mean(tmp2),
                                                p.value = wt$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)    
  } else {
    class <- "amp"
    tmp$class <- ifelse(tmp$cna < 0,"amp","others")
    tmp1 <- tmp[which(tmp$class == "amp"),"exp"]
    tmp2 <- tmp[which(tmp$class == "others"),"exp"]
    wt <- wilcox.test(tmp1,tmp2,alternative = "greater")
    cnacomp.gaby <- rbind.data.frame(cnacomp.gaby,
                                     data.frame(gene = g,
                                                loci = loci,
                                                category = class,
                                                avg.cna = mean(tmp1),
                                                avg.others = mean(tmp2),
                                                p.value = wt$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)  
  }
}
rownames(cnacomp.gaby) <- cnacomp.gaby$gene
write.table(cnacomp.gaby,file.path(res.path,"correlation between cna and expression in 4 loci of gaby cohort.txt"),sep = "\t",row.names = F,quote = F)

# TARGET cohort
Segment <- read.table(file.path(data.path,"wt_target_2018_pub/data_cna_hg19.seg"),sep="\t",header=T,stringsAsFactors = F)
colnames(Segment)<- c("Sample","Chrom","Start","Stop","Mark","Seg.CN")
write.table(Segment,file.path(res.path,"gistic input of target wilms cohort.txt"),sep = "\t",row.names = F,quote = F)

marker <- Segment[,1:4]
head(marker)

# create marker file for GISTIC2.0
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chrom"],2))
  c <- c(c,marker[i,"Start"],marker[i,"Stop"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"Target_wilms_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

# load gistic results with qvalue 0.25, cutoff 0.2, confident interval 95
gistic.target <- read.table(file.path(data.path,"236749/Target_wilms_0.2_95.all_thresholded.by_genes.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp1 <- gistic.target[,1:2]
tmp2 <- gistic.target[,3:ncol(gistic.target)]

colnames(tmp2) <- paste0(colnames(tmp2),"A")
com_sam <- intersect(colnames(tmp2),rownames(matchID))
tmp2 <- tmp2[,com_sam]
colnames(tmp2) <- matchID[com_sam,"SRA_RUN"]
gistic.target <- cbind.data.frame(tmp1,tmp2)

cnacomp.target <- NULL
for (i in 1:nrow(gistic.target)) {
  display.progress(index = i,totalN = nrow(gistic.target))
  
  # loss
  tmp <- gistic.target[i,]
  g <- rownames(tmp)
  loci <- tmp$Cytoband
  tmp <- data.frame(cna = as.numeric(tmp[,c(RNA_Clust1.target.cna,RNA_Clust2.target.cna)]),
                    group = rep(c("C1","C2"),c(67,27)),
                    stringsAsFactors = F)
  loss <- amp <- tmp
  loss$cna <- ifelse(loss$cna < 0,"LOSS","Others")
  amp$cna <- ifelse(amp$cna > 0,"AMP","Others")
  
  loss.dt <- as.data.frame.array(table(loss$cna,loss$group))
  amp.dt <- as.data.frame.array(table(amp$cna,amp$group))
  
  if(!is.element("AMP",rownames(amp.dt))) {
    amp.pct <- c(0,0)
    amp.p <- NA
  } else {
    amp.pct <- as.numeric(amp.dt[1,]/colSums(amp.dt))
    amp.p <- fisher.test(amp.dt)$p.value
  }
  
  if(!is.element("LOSS",rownames(loss.dt))) {
    loss.pct <- c(0,0)
    loss.p <- NA
  } else {
    loss.pct <- as.numeric(loss.dt[1,]/colSums(loss.dt))
    loss.p <- fisher.test(loss.dt)$p.value
  }
  
  cnacomp.target <- rbind.data.frame(cnacomp.target,
                                   data.frame(gene = g,
                                              loci = loci,
                                              loss.C1 = loss.pct[1],
                                              loss.C2 = loss.pct[2],
                                              p.loss = loss.p,
                                              amp.C1 = amp.pct[1],
                                              amp.C2 = amp.pct[2],
                                              p.amp = amp.p,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
  
}
rownames(cnacomp.target) <- cnacomp.target$gene
write.table(cnacomp.target,file.path(res.path,"comparison of gene level cna in target cohort between two immune clusters.txt"),sep = "\t",row.names = F,quote = F)

# copy number heatmap of microenvironment factors
micro.gene <- read.table(file.path(data.path,"microenviroment factors.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
rownames(micro.gene) <- micro.gene$Symbol
com_gene <- intersect(rownames(micro.gene),rownames(target.wt))
com_gene <- setdiff(com_gene,"HLA-DQB2")
com_cna <- intersect(com_gene,cnacomp.target$gene)
annRow.cna <- annRow[com_cna,,drop = F]

plotdata <- cnacomp.target[com_cna,c("loss.C1","loss.C2","amp.C1","amp.C2")]
blank <- "    " 
p.value <- cnacomp.target[com_cna,"p.loss"]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*","ns"))))
p.label <- formatC(p.value,
                   format = "e",
                   digits = 2)
library(stringr)
add.label <- str_pad(paste0(rownames(plotdata)," ",sig.label),
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")
pdf(file.path(fig.path,"cna loss percentage of microenviroment factors in target clusters.pdf"),width = 8,height = 22)
pheatmap(plotdata[,1:2],
         cluster_rows = F,cluster_cols = F,
         border_color = NA,
         color = bluered(64),
         annotation_row = annRow.cna[rownames(plotdata),,drop = F],
         annotation_colors = annColors.micro,
         gaps_row = c(47,89,97,121,125,131,134,149),
         #gaps_col = 2,
         labels_row = add.label)
invisible(dev.off())

plotdata <- cnacomp.target[com_cna,c("loss.C1","loss.C2","amp.C1","amp.C2")]
blank <- "    "
p.value <- cnacomp.target[com_cna,"p.amp"]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*","ns"))))
p.label <- formatC(p.value, 
                   format = "e",
                   digits = 2)
library(stringr)
add.label <- str_pad(paste0(rownames(plotdata)," ",sig.label),
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")
pdf(file.path(fig.path,"cna gain percentage of microenviroment factors in target clusters.pdf"),width = 8,height = 22)
pheatmap(plotdata[,3:4],
         cluster_rows = F,cluster_cols = F,
         border_color = NA,
         color = bluered(64),
         annotation_row = annRow.cna[rownames(plotdata),,drop = F],
         annotation_colors = annColors.micro,
         gaps_row = c(47,89,97,121,125,131,134,149),
         #gaps_col = 2,
         labels_row = add.label)
invisible(dev.off())

# redo target cohort using cbioportal data
micro.gene <- read.table(file.path(data.path,"microenviroment factors.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
rownames(micro.gene) <- micro.gene$Symbol
com_gene <- intersect(rownames(micro.gene),rownames(target.wt))
com_gene <- setdiff(com_gene,"HLA-DQB2")
micro.gene <- micro.gene[com_gene,]
com_cna <- intersect(com_gene,rownames(target.cna))

cnacomp.target.cbio <- NULL
for (i in 1:nrow(target.cna)) {
  display.progress(index = i,totalN = nrow(target.cna))
  
  # loss
  tmp <- target.cna[i,]
  g <- rownames(tmp)
  tmp <- data.frame(cna = as.numeric(tmp[,c(RNA_Clust1.target.cna,RNA_Clust2.target.cna,RNA_Clust3.target.cna)]),
                    group = rep(c("C1","C2"),c(67,27)),
                    stringsAsFactors = F)
  loss <- amp <- tmp
  loss$cna <- ifelse(loss$cna < 0,"LOSS","Others")
  amp$cna <- ifelse(amp$cna > 0,"AMP","Others")
  
  loss.dt <- as.data.frame.array(table(loss$cna,loss$group))
  amp.dt <- as.data.frame.array(table(amp$cna,amp$group))
  
  if(!is.element("AMP",rownames(amp.dt))) {
    amp.pct <- c(0,0)
    amp.p <- NA
  } else {
    amp.pct <- as.numeric(amp.dt[1,]/colSums(amp.dt))
    amp.p <- fisher.test(amp.dt)$p.value
  }
  
  if(!is.element("LOSS",rownames(loss.dt))) {
    loss.pct <- c(0,0)
    loss.p <- NA
  } else {
    loss.pct <- as.numeric(loss.dt[1,]/colSums(loss.dt))
    loss.p <- fisher.test(loss.dt)$p.value
  }
  
  cnacomp.target.cbio <- rbind.data.frame(cnacomp.target.cbio,
                                     data.frame(gene = g,
                                                loss.C1 = loss.pct[1],
                                                loss.C2 = loss.pct[2],
                                                p.loss = loss.p,
                                                amp.C1 = amp.pct[1],
                                                amp.C2 = amp.pct[2],
                                                p.amp = amp.p,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
  
}
rownames(cnacomp.target.cbio) <- cnacomp.target.cbio$gene
write.table(cnacomp.target.cbio,file.path(res.path,"comparison of gene level cna in target cohort between two immune clusters (cbioportal data).txt"),sep = "\t",row.names = F,quote = F)

micro.gene <- read.table(file.path(data.path,"microenviroment factors.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
rownames(micro.gene) <- micro.gene$Symbol
com_gene <- intersect(rownames(micro.gene),rownames(target.wt))
com_gene <- setdiff(com_gene,"HLA-DQB2")
com_cna <- intersect(com_gene,cnacomp.target.cbio$gene)
annRow.cna <- annRow[com_cna,,drop = F]

plotdata <- cnacomp.target.cbio[com_cna,c("loss.C1","loss.C2","amp.C1","amp.C2")]
blank <- "    "
p.value <- cnacomp.target.cbio[com_cna,"p.loss"]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*","ns"))))
p.label <- formatC(p.value, 
                   format = "e",
                   digits = 2)
library(stringr)
add.label <- str_pad(paste0(rownames(plotdata)," ",sig.label), 
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")
pdf(file.path(fig.path,"cna loss percentage of microenviroment factors in target clusters (cbioportal).pdf"),width = 8,height = 22)
pheatmap(plotdata[,1:2],
         cluster_rows = F,cluster_cols = F,
         border_color = NA,
         color = bluered(64),
         annotation_row = annRow.cna[rownames(plotdata),,drop = F],
         annotation_colors = annColors.micro,
         gaps_row = c(47,86,94,118,122,128,131,146),
         labels_row = add.label)
invisible(dev.off())

plotdata <- cnacomp.target.cbio[com_cna,c("loss.C1","loss.C2","amp.C1","amp.C2")]
blank <- "    "
p.value <- cnacomp.target.cbio[com_cna,"p.amp"]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*","ns"))))
p.label <- formatC(p.value,
                   format = "e",
                   digits = 2)
library(stringr)
add.label <- str_pad(paste0(rownames(plotdata)," ",sig.label),
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")
pdf(file.path(fig.path,"cna gain percentage of microenviroment factors in target clusters (cbioportal).pdf"),width = 8,height = 22)
pheatmap(plotdata[,3:4],
         cluster_rows = F,cluster_cols = F,
         border_color = NA,
         color = bluered(64),
         annotation_row = annRow.cna[rownames(plotdata),,drop = F],
         annotation_colors = annColors.micro,
         gaps_row = c(47,86,94,118,122,128,131,146),
         labels_row = add.label)
invisible(dev.off())

# common loss and gain in both cohort
com_loss <- intersect(cnacomp.gaby[which(cnacomp.gaby$p.loss < 0.05),"gene"],
                      cnacomp.target.cbio[which(cnacomp.target.cbio$p.loss < 0.05),"gene"])
com_loss <- intersect(micro.gene$Symbol,com_loss)

com_gain <- intersect(cnacomp.gaby[which(cnacomp.gaby$p.amp < 0.05),"gene"],
                      cnacomp.target.cbio[which(cnacomp.target.cbio$p.amp < 0.05),"gene"])
com_gain <- intersect(micro.gene$Symbol,com_gain)

#---------------------------------------------------------------#
# map full list of cibersort CD8 and CD4 to original CCR cluter #
cd8 <- immune.sig.me$T.cells.CD8
cd4.native <- immune.sig.me$T.cells.CD4.naive
cd4.activated <- immune.sig.me$T.cells.CD4.memory.activated
cd4.resting <- immune.sig.me$T.cells.CD4.memory.resting

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
plotdata <- target.wt[cd8,com_sam]
plotdata <- standarize.fun(plotdata,halfwidth = 2)
pdf(file.path(fig.path, "heatmap of cd8 in 114 wilms tumors in target cohort (purified by tp) 2 clusters.pdf"), height=20,width = 10)
pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = dendsort(hcs.target114),
         border_color = NA,
         annotation_col = annCol.target114[colnames(plotdata),c(45),drop = F],
         annotation_colors = annColors.target114,
         color = greenred(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 5,
         cellheight = 9,
         #gaps_row = c(14,22),
         cutree_cols = 2,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

#-------------------------------------------------------#
# prognostic value of CD4/CD8 compared with other cells #
rms2curve <- function(surv.df = NULL, immune.matrix = NULL, main.marker = NULL, tau = 60, n.grid = 100, fig.path = NULL,prefix = NULL,width = 4.5,height = 4.5) {
  
  require(survival)
  library(survcomp)
  library(survRM2)
  
  cat("Please make sure the survival time is based on month!\n")
  com_sam <- intersect(rownames(surv.df),colnames(immune.matrix))
  surv.df <- surv.df[com_sam,]
  #immune.matrix <- exp(immune.matrix[,com_sam])
  immune.matrix <- immune.matrix[,com_sam]
  
  surv.time <- surv.df$OS.time
  surv.event <- surv.df$OS
  marker1 <- as.numeric(main.marker[rownames(surv.df)])
  cat("1")
  for (cell in rownames(immune.matrix)) {
    marker2 <- as.numeric(immune.matrix[cell,rownames(surv.df)])
    
    # Set up data using age as the marker
    surv <- Surv(surv.time,surv.event)
    marker.pp <- seq(from=0,to=1,length=n.grid)
    marker1.qq <- quantile(marker1,marker.pp)
    marker2.qq <- quantile(marker2,marker.pp)
    
    fitdat.df1 <- data.frame(marker1=marker1)
    newdat.df1 <- data.frame(marker1=marker1.qq)
    
    fitdat.df2 <- data.frame(marker2=marker2)
    newdat.df2 <- data.frame(marker2=marker2.qq)
    # Calculations
    cat("2")
    cox.model1 <- coxph(surv~marker1,data=fitdat.df1)
    rms.calc1 <- summary(survfit(cox.model1,
                                 newdata=newdat.df1),rmean=tau)
    rms.mean1 <- rms.calc1$table[,"*rmean"]
    
    cox.model2 <- coxph(surv~marker2,data=fitdat.df2)
    rms.calc2 <- summary(survfit(cox.model2,newdata=newdat.df2),rmean=tau)
    rms.mean2 <- rms.calc2$table[,"*rmean"]
    
    tmp <- data.frame(OS.time = surv.df$OS.time,
                      OS = surv.df$OS,
                      main.marker = marker1,
                      minor.marker = marker2,
                      Cluster = surv.df$Cluster,
                      stringsAsFactors = F)
    fit1 <- coxph(Surv(OS.time,OS) ~ main.marker,data = tmp)
    cindex1 <- concordance.index(predict(fit1),surv.time = tmp$OS.time,surv.event = tmp$OS,method = "noether")
    fit2 <- coxph(Surv(OS.time,OS) ~ minor.marker,data = tmp)
    cindex2 <- concordance.index(predict(fit2),surv.time = tmp$OS.time,surv.event = tmp$OS,method = "noether")
    ccomp <- cindex.comp(cindex1, cindex2)
    
    # RMS Curve
    fig.name <- paste0(prefix,"_",cell,".pdf")
    pdf(file.path(fig.path,fig.name),width = width,height = height)
    par(bty="o", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
    plot(marker.pp,rms.mean1,type="p",pch = 19,col = jco[2],ylim = c(0,tau),
         xlab="Percentile of ES",ylab = "RMS",cex = 0.8)
    points(marker.pp,rms.mean2,pch = 19,col = jco[1],cex = 0.8)
    legend("bottomleft", 
           legend = c(paste0(cell, ": C-index = ",round(ccomp$cindex2,2)),paste0("CD4/CD8: C-index = ",round(ccomp$cindex1,2)),
                      paste0("P ",ifelse(ccomp$p.value < 0.001, "< 0.001", paste0("= ",round(ccomp$p.value,3))))),
           fill = c(jco[1],jco[2],NA),cex=0.8, border=NA, y.intersp=1, x.intersp=0.2,bty = "n")
    invisible(dev.off())
    
    cat(cell,"\n")
    print(list(fit1 = fit1,fit2 = fit2,cindex.compare = ccomp))
  }
}

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- as.data.frame(target.pure$es.pure[,com_sam])
indata <- as.data.frame(standarize.fun(indata,halfwidth = 2))
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`/30.5,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.target114[rownames(Sinfo.target),"tune_Clust"],
                  row.names = rownames(Sinfo.target),
                  stringsAsFactors = F)
tmp <- tmp[colnames(indata),]
tmp$Cluster <- ifelse(tmp$tune_Clust == "C1","hot","cold")
coxph(Surv(tmp$OS.time,tmp$OS)~as.numeric(indata["T.cells.CD4.memory.activated",]))
coxph(Surv(tmp$OS.time,tmp$OS)~as.numeric(indata["T.cells.CD4.memory.resting",]))
coxph(Surv(tmp$OS.time,tmp$OS)~as.numeric(indata["T.cells.CD4.naive",]))
coxph(Surv(tmp$OS.time,tmp$OS)~as.numeric(indata["T.cells.CD8",]))
tmp1 <- abs(as.numeric(indata["T.cells.CD4.memory.activated",])-as.numeric(indata["T.cells.CD8",]))/max(abs(as.numeric(indata["T.cells.CD4.memory.activated",])),abs(as.numeric(indata["T.cells.CD8",])))

coxph(Surv(tmp$OS.time,tmp$OS)~tmp1)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

cd48.activated <- indata["T.cells.CD4.memory.activated",]-indata["T.cells.CD8",]; names(cd48.activated) <- colnames(indata)
cd48.naive <- indata["T.cells.CD4.naive",]-indata["T.cells.CD8",]; names(cd48.naive) <- colnames(indata)
cd48.resting <- indata["T.cells.CD4.memory.resting",]-indata["T.cells.CD8",]; names(cd48.resting) <- colnames(indata)

indata <- as.data.frame(t(apply(indata, 1, range01)))
cd48.activated <- range01(cd48.activated)

rms2curve(surv.df = tmp,
          immune.matrix = indata,
          main.marker = cd48.activated,
          fig.path = fig.path,
          prefix = "CD48.activated")

tmp1 <- as.numeric(cd48.activated[RNA_Clust1.target])
tmp2 <- as.numeric(cd48.activated[RNA_Clust2.target])
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- as.numeric(indata["T.cells.CD4.memory.activated",RNA_Clust1.target])
tmp2 <- as.numeric(indata["T.cells.CD4.memory.activated",RNA_Clust2.target])
t.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- as.numeric(indata["T.cells.CD8",RNA_Clust1.target])
tmp2 <- as.numeric(indata["T.cells.CD8",RNA_Clust2.target])
tmp3 <- as.numeric(indata["T.cells.CD8",RNA_Clust3.target])
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

#-------------------------------------------#
# KM curve for CD4 and CD8 with best cutoff #
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- as.data.frame(target.pure$es.pure[,com_sam])
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  tune_Clust = annCol.target114[rownames(Sinfo.target),"tune_Clust"],
                  row.names = rownames(Sinfo.target),
                  stringsAsFactors = F)
tmp <- tmp[colnames(indata),]
tmp$Cluster <- ifelse(tmp$tune_Clust == "C1","hot","cold")
tmp$cd4 <- as.numeric(indata["T.cells.CD4.memory.activated",])
tmp$cd8 <- as.numeric(indata["T.cells.CD8",])
tmp$cd48 <- abs(tmp$cd4-tmp$cd8)/max(abs(tmp$cd4),abs(tmp$cd8))

cd48.cut <- surv_cutpoint(tmp, time = "OS.time", 
                         event = "OS", 
                         variables = "cd48", 
                         minprop = 0.2) 
cd48.cat <- surv_categorize(cd48.cut)
fitd <- survdiff(Surv(OS.time, OS) ~ cd48, data=cd48.cat, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ cd48, data=cd48.cat, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jama",
                pval = T,
                data=cd48.cat,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for CCR ratio of CD4 and CD8 with best cutoff in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

cd8.cut <- surv_cutpoint(tmp, time = "OS.time", 
                          event = "OS", 
                          variables = "cd8", 
                          minprop = 0.2) 
cd8.cat <- surv_categorize(cd8.cut)
fitd <- survdiff(Surv(OS.time, OS) ~ cd8, data=cd8.cat, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ cd8, data=cd8.cat, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jama",
                pval = T,
                data=cd8.cat,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for CCR CD8 with best cutoff in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

cd4.cut <- surv_cutpoint(tmp, time = "OS.time", 
                         event = "OS", 
                         variables = "cd4", 
                         minprop = 0.2) 
cd4.cat <- surv_categorize(cd4.cut)
fitd <- survdiff(Surv(OS.time, OS) ~ cd4, data=cd4.cat, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ cd4, data=cd4.cat, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jama",
                pval = T,
                data=cd4.cat,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for CCR CD4 with best cutoff in TARGET cohort.pdf"),width = 5,height = 6)
p
invisible(dev.off())

tmp1 <- tmp[which(tmp$Cluster == "hot"),"cd4"]
tmp2 <- tmp[which(tmp$Cluster == "cold"),"cd4"]
wilcox.test(tmp1,tmp2)
box.df <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Infiltrated","Desert"),c(length(tmp1),length(tmp2))))
box.df$mut <- factor(box.df$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for cd4 114 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = box.df,
        xlab = "Immunity status",
        ylab = "CD4+",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = box.df, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,17,"P = 0.161",cex = 1)
invisible(dev.off())

tmp1 <- tmp[which(tmp$Cluster == "hot"),"cd8"]
tmp2 <- tmp[which(tmp$Cluster == "cold"),"cd8"]
wilcox.test(tmp1,tmp2)
box.df <- data.frame(exp = c(tmp1,tmp2),
                     mut = rep(c("Infiltrated","Desert"),c(length(tmp1),length(tmp2))))
box.df$mut <- factor(box.df$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for cd8 114 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = box.df,
        xlab = "Immunity status",
        ylab = "CD8+",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = box.df, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,17,"P < 0.001",cex = 1)
invisible(dev.off())

tmp1 <- tmp[which(tmp$Cluster == "hot"),"cd48"]
tmp2 <- tmp[which(tmp$Cluster == "cold"),"cd48"]
wilcox.test(tmp1,tmp2)
box.df <- data.frame(exp = c(tmp1,tmp2),
                     mut = rep(c("Infiltrated","Desert"),c(length(tmp1),length(tmp2))))
box.df$mut <- factor(box.df$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for ratio of cd4 and cd8 114 WT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = box.df,
        xlab = "Immunity status",
        ylab = "Ratio of CD4+/CD8+",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = box.df, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,1.3,"P < 0.001",cex = 1)
invisible(dev.off())

# use cibersort results
tmp <- read.table(file.path(res.path,"CIBERSORT.Output.relative.Wilms_tumor_TARGET.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
tmp <- tmp[,1:22]
indata <- as.data.frame(t(tmp))
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  row.names = rownames(Sinfo.target),
                  stringsAsFactors = F)
tmp <- tmp[colnames(indata),]
tmp$cd4 <- as.numeric(indata["T cells CD4 memory activated",])
tmp$cd8 <- as.numeric(indata["T cells CD8",])
tmp$cd48 <- tmp$cd4-tmp$cd8

# univariate cox of cd4 and cd8 genes using full list of cibersort
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(viridis)
library(scales)

realdata <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                       OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                       row.names = rownames(Sinfo.target),
                       stringsAsFactors = F)
realdata <- cbind.data.frame(realdata,t(target.wt[cd4.activated,rownames(realdata)]))
realdata <- realdata[target.pure$keep.sam,]
Coxoutput=data.frame()
for(i in colnames(realdata[,3:ncol(realdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ realdata[,i], data = realdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
write.table(Coxoutput,file.path(res.path,"univariate cox of CD4 genes in target cohort.txt"),sep = "\t",row.names = F,quote = F)

realdata <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`,
                       OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                       row.names = rownames(Sinfo.target),
                       stringsAsFactors = F)
realdata <- cbind.data.frame(realdata,t(target.wt[cd8,rownames(realdata)]))
realdata <- realdata[target.pure$keep.sam,]
Coxoutput=data.frame()
for(i in colnames(realdata[,3:ncol(realdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ realdata[,i], data = realdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
write.table(Coxoutput,file.path(res.path,"univariate cox of CD8 genes in target cohort.txt"),sep = "\t",row.names = F,quote = F)

#-----------------------#
# Other TME calculation #
library(psych)
immune.checkpoint <- c("CD279","CD274","PDCD1LG2","CTLA4","HAVCR2","LAG3")
immune.checkpoint <- intersect(immune.checkpoint,rownames(target.wt))

immunosuppression <- c("CXCL12","TGFB1","TGFB3","LGALS1")
t.cell.activation <- c("CXCL9","CXCL10","CXCL16","IFNG","IL15")
t.cell.survival <- c("CD70","CD27")
regulatory.t.cell <- c("FOXP3","TNFRSF18")
MHC.I <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","B2M")
MCC <- c("CCL2")
TLS <- c("CCL21","CCL19","CXCL13","CXCL11","CCL8","CXCL10","CXCL9","CCL2","CCL3","CCL18","CCL5","CCL4")

geometric.mean.rm0 <- function (x, na.rm = TRUE) 
{
  x <- x[x>0]
  if(sum(x) == 0) {
    return(0)
  } else if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = TRUE))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}

immunosuppression.target <- apply(target.wt[immunosuppression,], 2, geometric.mean.rm0)
t.cell.activation.target <- apply(target.wt[t.cell.activation,], 2, geometric.mean.rm0)
t.cell.survival.target <- apply(target.wt[t.cell.survival,], 2, geometric.mean.rm0)
regulatory.t.cell.target <- apply(target.wt[regulatory.t.cell,], 2, geometric.mean.rm0)
MHC.I.target <- apply(target.wt[MHC.I,], 2, geometric.mean.rm0)
TLS.target <- apply(target.wt[TLS,], 2, geometric.mean.rm0)

indata <- rbind.data.frame(immunosuppression.target,
                           t.cell.activation.target,
                           t.cell.survival.target,
                           regulatory.t.cell.target,
                           MHC.I.target,
                           target.wt[MCC,],
                           #target.wt[TLS,],
                           TLS.target)
rownames(indata) <- c("Immunosuppression",
                      "T cell activation",
                      "T cell survival",
                      "Regulatory T cell",
                      "Class MHC I",
                      "myeloid cell chemotaxis",
                      "tertiary lymphoid structures")
colnames(indata) <- colnames(target.wt)
TME.raw.dat <- indata

#-----------------------------------------------------#
# Recluster to see if immune group could be seperated #
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
group <- cutree(hcs.target114,4)
annCol.target114$tune_Clust2 <- NA
annCol.target114[com_sam,"tune_Clust2"] <- paste0("C",group[com_sam])
annCol.target114$tune_Clust2 <- ifelse(annCol.target114$tune_Clust2 == "C1","C4",
                                      ifelse(annCol.target114$tune_Clust2 == "C3", "C3",
                                             ifelse(annCol.target114$tune_Clust2 == "C2","C2","C1")))
annCol.target114$tune_Clust3 <- ifelse(annCol.target114$tune_Clust2 == "C1","C1",
                         ifelse(annCol.target114$tune_Clust2 == "C2","C2",
                                ifelse(annCol.target114$tune_Clust2 %in% c("C3","C4"),"C3",NA)))
annColors.target114[["tune_Clust2"]] <- c("C1" = jco[2],"C2" = jco[1],"C3" = seagreen,"C4" = sun)
annColors.target114[["tune_Clust3"]] <- c("C1" = jco[2],"C2" = jco[1],"C3" = seagreen)

pdf(file.path(fig.path, "comprehensive immune signature heatmap of 114 wilms tumors in target cohort (purified by tp) 4 clusters.pdf"), height=62,width = 15)
pheatmap(as.matrix(TME.dat[setdiff(row.order,c("Regulatory T cell","myeloid cell chemotaxis")),]),
         cluster_rows = F,
         cluster_cols = dendsort(hcs.target114),
         border_color = NA,
         annotation_col = annCol.target114[colnames(plotdata),c(45,46,47,5,4,41,37,42,43,44,1:3,6:21,24:36)],
         annotation_colors = annColors.target114,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 5,
         cellheight = 10,
         gaps_row = c(14,22,24,29),
         #cutree_cols = 4,
         show_colnames = F,
         show_rownames = T,
         fontsize_col = 8)
invisible(dev.off())

tmp <- data.frame(OS.time = annCol.target114$OS.time,
                  OS = ifelse(annCol.target114$OS == "Dead",1,0),
                  tune_Clust = annCol.target114$tune_Clust,
                  Age = annCol.target114$Age/365,
                  Stage = annCol.target114$Stage,
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1],seagreen,sun),
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort 4 clusters.pdf"),width = 5,height = 6)
p
invisible(dev.off())

pairwise_survdiff(Surv(OS.time, OS)~ tune_Clust, data=tmp)
# C1      C2      C3     
# C2 0.26489 -       -      
# C3 0.40697 0.09169 -      
# C4 0.04582 0.00075 0.26489

tmp <- data.frame(OS.time = annCol.target114$OS.time,
                  OS = ifelse(annCol.target114$OS == "Dead",1,0),
                  tune_Clust = annCol.target114$tune_Clust3,
                  Age = annCol.target114$Age/365,
                  Stage = annCol.target114$Stage,
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1],seagreen),
                pval = T,
                data=tmp,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf(file.path(fig.path,"KM for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort modified 3 clusters.pdf"),width = 5,height = 6)
p
invisible(dev.off())

pairwise_survdiff(Surv(OS.time, OS)~ tune_Clust, data=tmp, p.adjust.method = "none")
# C1     C2    
# C2 0.1993 -     
# C3 0.0736 0.0014

dd <- as.data.frame(t(TME.raw.dat[c(1,2,3,5,7),target.pure$keep.sam]))
dd$RNA_Clust <- annCol.target114[rownames(dd),"tune_Clust3"]
dd$sample <- rownames(dd)
d2 <- gather(dd, cell, expr, 1:5)
d2$cell <- factor(d2$cell,levels = rownames(indata))

pvalues <- sapply(d2$cell, function(x) {
  res <- kruskal.test(expr ~ RNA_Clust, data = subset(d2, cell == x))
  return(res$p.value)
})
pv <- data.frame(gene = d2$cell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

ggplot(d2, aes(cell, expr, fill=RNA_Clust)) + scale_fill_manual(values = c(jco[2],jco[1],seagreen)) +
  geom_boxplot() + 
  geom_text(aes(gene, y=max(d2$expr) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("Enrichment score") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(fig.path,"boxplot for other TMEs among modified three clusters in 114 target cohort.pdf"),width = 6,height = 5)

dd <- as.data.frame(t(target.wt[immune.checkpoint,target.pure$keep.sam]))
dd$RNA_Clust <- annCol.target114[rownames(dd),"tune_Clust3"]
dd$sample <- rownames(dd)
d2 <- gather(dd, cell, expr, 1:5)
d2$cell <- factor(d2$cell,levels = immune.checkpoint)

pvalues <- sapply(d2$cell, function(x) {
  res <- kruskal.test(expr ~ RNA_Clust, data = subset(d2, cell == x))
  return(res$p.value) 
})
pv <- data.frame(gene = d2$cell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

ggplot(d2, aes(cell, expr, fill=RNA_Clust)) + scale_fill_manual(values = c(jco[2],jco[1],seagreen)) +
  geom_boxplot() + 
  geom_text(aes(gene, y=max(d2$expr) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("Enrichment score") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(fig.path,"boxplot for immune checkpoint among modified three clusters in 114 target cohort.pdf"),width = 6,height = 5)

TME_Clust1.target <- rownames(annCol.target114[which(annCol.target114$tune_Clust3 == "C1"),])
TME_Clust2.target <- rownames(annCol.target114[which(annCol.target114$tune_Clust3 == "C2"),])
TME_Clust3.target <- rownames(annCol.target114[which(annCol.target114$tune_Clust3 == "C3"),])

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- TME.raw.dat[c(1,2,3,5,7),com_sam]
outTab <- NULL
for(cell in rownames(indata)) {
  tmp1 <- as.numeric(indata[cell,TME_Clust1.target])
  tmp2 <- as.numeric(indata[cell,TME_Clust2.target])
  tmp3 <- as.numeric(indata[cell,TME_Clust3.target])
  
  kruskal.p <- kruskal.test(list(tmp1,tmp2,tmp3))$p.value
  wilcox.p12 <- wilcox.test(tmp1,tmp2)$p.value
  wilcox.p13 <- wilcox.test(tmp1,tmp3)$p.value
  wilcox.p23 <- wilcox.test(tmp2,tmp3)$p.value
  wilcox.p12.3 <- wilcox.test(tmp3,c(tmp1,tmp2))$p.value
  
  outTab <- rbind.data.frame(outTab,data.frame(cell = cell,
                                               avg.C1 = mean(tmp1),
                                               avg.C2 = mean(tmp2),
                                               avg.C3 = mean(tmp3),
                                               p.kruskal = kruskal.p,
                                               p12.wilcox = wilcox.p12,
                                               p13.wilcox = wilcox.p13,
                                               p23.wilcox = wilcox.p23,
                                               p12.3.wilcox = wilcox.p12.3,
                                               stringsAsFactors = F),
                             stringsAsFactors = F)
}
write.table(outTab,file.path(res.path,"comparision of TME between three modified target clusters.txt"),sep = "\t",row.names = F,quote = F)

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- target.wt[immune.checkpoint,com_sam]
outTab <- NULL
for(cell in rownames(indata)) {
  tmp1 <- as.numeric(indata[cell,TME_Clust1.target])
  tmp2 <- as.numeric(indata[cell,TME_Clust2.target])
  tmp3 <- as.numeric(indata[cell,TME_Clust3.target])
  
  kruskal.p <- kruskal.test(list(tmp1,tmp2,tmp3))$p.value
  wilcox.p12 <- wilcox.test(tmp1,tmp2)$p.value
  wilcox.p13 <- wilcox.test(tmp1,tmp3)$p.value
  wilcox.p23 <- wilcox.test(tmp2,tmp3)$p.value
  wilcox.p12.3 <- wilcox.test(tmp3,c(tmp1,tmp2))$p.value
  
  outTab <- rbind.data.frame(outTab,data.frame(cell = cell,
                                               avg.C1 = mean(tmp1),
                                               avg.C2 = mean(tmp2),
                                               avg.C3 = mean(tmp3),
                                               p.kruskal = kruskal.p,
                                               p12.wilcox = wilcox.p12,
                                               p13.wilcox = wilcox.p13,
                                               p23.wilcox = wilcox.p23,
                                               p12.3.wilcox = wilcox.p12.3,
                                               stringsAsFactors = F),
                             stringsAsFactors = F)
}
write.table(outTab,file.path(res.path,"comparision of immune checkpoint between three modified target clusters.txt"),sep = "\t",row.names = F,quote = F)

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- as.data.frame(TME.raw.dat[,com_sam])
tmp <- data.frame(OS.time = Sinfo.target$`Overall Survival Time in Days`/30.5,
                  OS = ifelse(Sinfo.target$`Vital Status` == "Dead",1,0),
                  row.names = rownames(Sinfo.target),
                  stringsAsFactors = F)
tmp <- tmp[colnames(indata),]
coxph(Surv(tmp$OS.time,tmp$OS)~as.numeric(indata["Immunosuppression",]))
coxph(Surv(tmp$OS.time,tmp$OS)~as.numeric(indata["tertiary lymphoid structures",]))

# differentially anaysis between modified C1 vs C2
# mRNA
tmp <- annCol.target114$tune_Clust3; names(tmp) <- rownames(annCol.target114)
tmp <- na.omit(tmp[tmp!="C3"])
Groupinfo <- annCol.target114[which(annCol.target114$tune_Clust3 != "C3"),]
colnames(Groupinfo)[47] <- "group"
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.TME.target(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="TARGET_WT_mRNA_modified_3_clusters", complist,PASSFlag=PASSFlag, overwt=TRUE)

# GSEA
deres <- read.table(file.path(res.path,"TARGET_WT_mRNA_modified_3_clusters_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- deres$logFC
names(geneList) <- deres$id
geneList <- sort(geneList,decreasing = T)
gsea_target_modified_C1vsC2 <- GSEA(geneList = geneList,
                                    TERM2GENE=MSigDB,
                                    pvalueCutoff = 0.25,
                                    nPerm = 10000,
                                    seed = T,
                                    verbose=F)
res <- data.frame(gsea_target_modified_C1vsC2)
write.table(as.data.frame(gsea_target_modified_C1vsC2),file.path(res.path,"GSEA results for TARGET WT mRNA modified 3 clusters of C1 vs C2 DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

#--------------------#
# TME in Gaby cohort #
immunosuppression.gaby <- apply(gaby.wt[immunosuppression,], 2, geometric.mean.rm0)
t.cell.activation.gaby <- apply(gaby.wt[t.cell.activation,], 2, geometric.mean.rm0)
t.cell.survival.gaby <- apply(gaby.wt[t.cell.survival,], 2, geometric.mean.rm0)
regulatory.t.cell.gaby <- apply(gaby.wt[regulatory.t.cell,], 2, geometric.mean.rm0)
MHC.I.gaby <- apply(gaby.wt[setdiff(MHC.I,"HLA-G"),], 2, geometric.mean.rm0)
TLS.gaby <- apply(gaby.wt[TLS,], 2, geometric.mean.rm0)

indata <- rbind.data.frame(immunosuppression.gaby,
                           t.cell.activation.gaby,
                           t.cell.survival.gaby,
                           regulatory.t.cell.gaby,
                           MHC.I.gaby,
                           gaby.wt[MCC,],
                           TLS.gaby)
rownames(indata) <- c("Immunosuppression",
                      "T cell activation",
                      "T cell survival",
                      "Regulatory T cell",
                      "Class MHC I",
                      "myeloid cell chemotaxis",
                      "tertiary lymphoid structures")
colnames(indata) <- colnames(gaby.wt)

plotdata <- imm.score.gaby.wt3.pure[,rownames(cluster2.wt.gaby.mRNA.ans$annCol)]
plotdata <- rbind.data.frame(plotdata,indata[,rownames(cluster2.wt.gaby.mRNA.ans$annCol)])
plotdata <- rbind.data.frame(plotdata,gaby.wt[immune.checkpoint,rownames(cluster2.wt.gaby.mRNA.ans$annCol)])
plotdata <- standarize.fun(plotdata,halfwidth = 2)

pdf(file.path(fig.path, "comprehensive immune signature heatmap of wilms tumors in gaby cohort (purified by tp).pdf"), height=26,width = 10)
pheatmap(as.matrix(plotdata[setdiff(row.order,c("Regulatory T cell","myeloid cell chemotaxis")),]),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = "black",
         annotation_col = annCol.wt.gaby.genomic.alter[,setdiff(colnames(annCol.wt.gaby.genomic.alter),c("TP53 mutations","1q gain","16p STATUS"))],
         annotation_colors = annColors.wt.gaby.genomic.alter,
         color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 15,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         cutree_cols = 2,
         gaps_row = c(14,22,24,29),
         fontsize_col = 8)
invisible(dev.off())

#--------------------------------#
# check EZH2 and CBX2 expression #
ezh2.gaby.c1 <- as.numeric(gaby.wt["EZH2",RNA_Clust1])
ezh2.gaby.c2 <- as.numeric(gaby.wt["EZH2",RNA_Clust2])
wilcox.test(ezh2.gaby.c1,ezh2.gaby.c2) #0.007937
tmp <- data.frame(exp = c(ezh2.gaby.c1,ezh2.gaby.c2),
                  mut = rep(c("C1","C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for EZH2 expression in gaby WT cohort regarding RNA cluster.pdf"),width = 1.8,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "RNA Cluster",
        ylab = "log2(TPM) in Gaby",
        ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,5,"P = 0.008",cex = 1)
invisible(dev.off())

# EZH2 pathway
tmp <- MSigDB[grep("KAMMINGA_EZH2_TARGETS",MSigDB$ont),]
ezh2.path.gaby <- gsva(expr = as.matrix(gaby.wt),
                       gset.idx.list = list("KAMMINGA_EZH2_TARGETS" = tmp$gene),
                       method = "ssgsea")
ezh2.path.gaby.c1 <- as.numeric(ezh2.path.gaby[1,RNA_Clust1])
ezh2.path.gaby.c2 <- as.numeric(ezh2.path.gaby[1,RNA_Clust2])
wilcox.test(ezh2.path.gaby.c1,ezh2.path.gaby.c2) #0.007937
tmp <- data.frame(exp = c(ezh2.path.gaby.c1,ezh2.path.gaby.c2),
                  mut = rep(c("C1","C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for EZH2 target in gaby WT cohort regarding RNA cluster.pdf"),width = 1.8,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "RNA Cluster",
        ylab = "ssGSEA score",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,1,"P = 0.008",cex = 1)
invisible(dev.off())

ezh2.target.c1 <- as.numeric(target.wt["EZH2",RNA_Clust1.target])
ezh2.target.c2 <- as.numeric(target.wt["EZH2",RNA_Clust2.target])
wilcox.test(ezh2.target.c1,ezh2.target.c2) #9.02e-05
tmp <- data.frame(exp = c(ezh2.target.c1,ezh2.target.c2),
                  mut = rep(c("C1","C2"),c(69,27)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for EZH2 expression in target WT cohort regarding RNA cluster.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "RNA Cluster",
        ylab = "log2(TPM) in target",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,7.5,"P < 0.001",cex = 1)
invisible(dev.off())

# EZH2 pathway
tmp <- MSigDB[grep("KAMMINGA_EZH2_TARGETS",MSigDB$ont),]
ezh2.path.target <- gsva(expr = as.matrix(target.wt),
                       gset.idx.list = list("KAMMINGA_EZH2_TARGETS" = tmp$gene),
                       method = "ssgsea")
ezh2.path.target.c1 <- as.numeric(ezh2.path.target[1,RNA_Clust1.target])
ezh2.path.target.c2 <- as.numeric(ezh2.path.target[1,RNA_Clust2.target])
wilcox.test(ezh2.path.target.c1,ezh2.path.target.c2) #6.754e-05
tmp <- data.frame(exp = c(ezh2.path.target.c1,ezh2.path.target.c2),
                  mut = rep(c("C1","C2"),c(69,27)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for EZH2 target in target WT cohort regarding RNA cluster.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "RNA Cluster",
        ylab = "ssGSEA score",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,1.6,"P < 0.001",cex = 1)
invisible(dev.off())

#----------------------------------------------------------------------#
# create oncoprint with mutated features, cna and clinical information #
library(ComplexHeatmap)
library(ggplot2)
mygene  <- read.table(file.path(data.path,"Gaby_WT_oncoprint_input.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mypathway <- data.frame(PathwayID = rep(c("Mutation","CNV"),c(3,5)),
                        row.names = rownames(mygene))
mygene$pathway <- mypathway[as.character(rownames(mygene)),]

gaby.wt.sinfo <- read.csv(file.path(data.path,"all clinical data exome 23122020.csv"),row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
gaby.wt.sinfo <- gaby.wt.sinfo[which(gaby.wt.sinfo$samID != ""),]
rownames(gaby.wt.sinfo) <- gaby.wt.sinfo$samID

tmp <- gaby.wt.sinfo[colnames(mygene)[1:(ncol(mygene)-1)],c("Age at diagnostic (years)","Sexe","Risk (review/local)","Histology type (review/local)","Overall Stage (review/local)")]
colnames(tmp) <- c("Age","Gender","Risk","Histology","Stage")
tmp$Age <- round(tmp$Age,2)
tmp$Age <- ifelse(tmp$Age > 4.85, ">4.85", "<=4.85")
tmp$Histology <- gsub(" (311)", "", tmp$Histology, fixed = T)
tmp$Histology <- gsub(" (312)", "", tmp$Histology, fixed = T)
tmp$Gender <- ifelse(tmp$Gender == "F","Female","Male")
tmp$Risk <- gsub(" Risk","",tmp$Risk)
my_annotation = HeatmapAnnotation(df = tmp,
                                  col = list(Age = c(">4.85" = orange, "<=4.85" = alpha(orange,0.6)),
                                             Gender = c("Female" = jco[2], "Male" = jco[1]),
                                             Risk = c("High" = darkred, "Intermediate" = alpha(darkred,0.6)),
                                             Histology = c("ANA diffuse" = darkblue,"Focal anaplasia" = darkred),
                                             Stage = c("II" = alpha(seagreen,0.3), "III" = alpha(seagreen,0.6), "IV" = seagreen)),
                                  border = T,
                                  height = unit(30,"mm"))
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "grey90", col = NA))
  },
  SM = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = darkred, col = NA)) 
  },
  GM = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = lightred, col = NA)) 
  },
  Yes = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = darkblue, col = NA)) 
  }
)

col = c("SM" = darkred, 
        "GM" = lightred, 
        "Yes" = darkblue)

top_anno <- HeatmapAnnotation(TMB = anno_barplot(c(variants_per_sample[colnames(mygene[,-ncol(mygene)]),"TMB"]),
                                                 border = T,
                                                 gp = gpar(fill = "black"),
                                                 height = unit(20,"mm")))
p <- oncoPrint(mygene[1:(ncol(mygene)-1)],
          row_split = mygene$pathway,
          get_type = function(x) x,
          alter_fun = alter_fun, col = col,
          remove_empty_columns = FALSE,
          column_order = colnames(mygene)[1:12],
          show_pct = T,
          show_column_names = T,
          right_annotation = NULL,
          top_annotation = top_anno,
          bottom_annotation = my_annotation,
          show_heatmap_legend = T)
pdf(file.path(fig.path,"somatic muation and cna with clinical information in Gaby WT cohort.pdf"), width = 6, height = 5.8)
draw(p, annotation_legend_side = "bottom")
invisible(dev.off())

#-----------------------------#
# forest plot for gaby cohort #
library(survival)
library(forestplot)
library(survminer)
tmp <- gaby.wt.sinfo[colnames(mygene)[1:(ncol(mygene)-1)],c("Age at diagnostic (years)","Sexe","Risk (review/local)","Histology type (review/local)","Overall Stage (review/local)","OS","OS.time","PFS","PFS.time")]
colnames(tmp) <- c("Age","Gender","Risk","Histology","Stage","OS","OS.time","PFS","PFS.time")
tmp$Age <- round(tmp$Age,2)
tmp$Age <- ifelse(tmp$Age > 4.85, ">4.85", "<=4.85")
tmp$Histology <- gsub(" (311)", "", tmp$Histology, fixed = T)
tmp$Histology <- gsub(" (312)", "", tmp$Histology, fixed = T)
tmp$Gender <- ifelse(tmp$Gender == "F","Female","Male")
tmp$Risk <- gsub(" Risk","",tmp$Risk)
tmp$RNA_clust <- NA
tmp[RNA_Clust1,"RNA_clust"] <- "C1"
tmp[RNA_Clust2,"RNA_clust"] <- "C2"
tmp$TP53 <- as.character(mygene["TP53",rownames(tmp)])
tmp$TP53 <- ifelse(tmp$TP53 != "",1, 0)
fisher.test(table(tmp$TP53,tmp$RNA_clust))

fitd <- survdiff(Surv(OS.time, OS) ~ RNA_clust,
                 data      = tmp,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ RNA_clust,
               data      = tmp,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("RNA_clust=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = TRUE,
                risk.table.col    = "strata",
                palette           = jco[2:1],
                data              = tmp,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                surv.median.line  = "hv",
                xlab              = "Time (Years)",
                ylab              = "Survival probability (%)",
                risk.table.y.text = FALSE)
p$plot <- p$plot + scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0,100,25))
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.55,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km in 2 RNA_clust in Gaby wt cohort.pdf"), width = 4.5, height = 5)
print(p)
dev.off()

data <- read_xlsx("E:/IGBMC/myproject/Wilms/Manuscript/New version/Submission/Nat Commun/Supplementary Tables-R1.xlsx", sheet = 1, skip = 1)
data <- as.data.frame(data)
tmp <- data[which(data$Has_RNAseq == "Yes"),]
tmp$RNA_clust <- rep(c("dWT","iWT"), c(5,5))
tmp$RNA_clust <- factor(tmp$RNA_clust, levels = c("iWT","dWT"))
tmp$RFS.time <- tmp$RFS.time/12
fitd <- survdiff(Surv(RFS.time, RFS) ~ RNA_clust,
                 data      = tmp,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(RFS.time, RFS)~ RNA_clust,
               data      = tmp,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("RNA_clust=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = TRUE,
                risk.table.col    = "strata",
                palette           = jco[2:1],
                data              = tmp,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                surv.median.line  = "hv",
                xlab              = "Time (Years)",
                ylab              = "Recurrence-free survival (%)",
                risk.table.y.text = FALSE)
p$plot <- p$plot + scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0,100,25))
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 2, y = 0.55,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km recurrence in 2 RNA_clust in Gaby wt cohort.pdf"), width = 4.5, height = 5)
print(p)
dev.off()

#----------------------------------------#
# compare CD count ragarding TP53 status #
cdcount <- read.table(file.path(data.path,"CDcount_TP53.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tp53.yes <- cdcount[which(cdcount$`TP53 gene (NM_001276697) mutations` %in% c("Yes","Germline mutation yes")),]
tp53.no <- cdcount[which(cdcount$`TP53 gene (NM_001276697) mutations` %in% c("NO")),]

tmp1 <- na.omit(cdcount[rownames(tp53.yes),"CD3"])
tmp2 <- na.omit(cdcount[rownames(tp53.no),"CD3"])
wilcox.test(tmp1,tmp2) # 0.5303

tmp1 <- na.omit(cdcount[rownames(tp53.yes),"CD4"])
tmp2 <- na.omit(cdcount[rownames(tp53.no),"CD4"])
wilcox.test(tmp1,tmp2) # 0.7436

tmp1 <- na.omit(cdcount[rownames(tp53.yes),"CD8"])
tmp2 <- na.omit(cdcount[rownames(tp53.no),"CD8"])
wilcox.test(tmp1,tmp2) # 0.07323

tmp1 <- na.omit(cdcount[rownames(tp53.yes),"CD4:CD8"])
tmp2 <- na.omit(cdcount[rownames(tp53.no),"CD4:CD8"])
wilcox.test(tmp1,tmp2) # 0.002525

tmp <- cdcount[c(rownames(tp53.yes),rownames(tp53.no)),c("TP53 gene (NM_001276697) mutations","CD3","CD4","CD8","CD4:CD8")]
tmp <- as.data.frame(na.omit(tmp))
colnames(tmp)[1] <- "TP53"
tmp$TP53 <- ifelse(tmp$TP53 == "NO","Wild","Mutated")
tmp$TP53 <- factor(tmp$TP53,levels = c("Mutated","Wild"))
tmp$`CD4:CD8` <- tmp$`CD4:CD8` * 10
library(data.table)
long <- as.data.frame(melt(tmp, id.vars = c("TP53"), variable.name = "Lymph"))


std <- function(x) sd(x)/sqrt(length(x))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      std = std(x[[col]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(long, varname="value", 
                    groupnames=c("TP53", "Lymph"))
# Convert dose to a factor variable
df2$Lymph=as.factor(df2$Lymph)
head(df2)

df2.1 <- df2[which(df2$Lymph %in% c("CD3","CD4","CD8")),]
df2.2 <- df2[-which(df2$Lymph %in% c("CD3","CD4","CD8")),]

# Default bar plot
p<- ggplot(df2, aes(x=Lymph, y=value, fill=TP53)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-std, ymax=value+std), width=.2,
                position=position_dodge(.9)) +labs(title="", x="Lymphocytes ", y = "Count")+
  theme_classic() +
  theme(legend.position = "top") + 
  scale_fill_manual(values=c(heatmap.L.BlYlRd[5],heatmap.L.BlYlRd[1]))
p
ggsave(filename = file.path(fig.path,"barplot for CD count regarding TP53 mutations.pdf"), width = 3,height = 4)

#-----------------#
# TIDE prediction #
TIDE <- round(sweep(gaby.wt,2, apply(gaby.wt, 2, median)),2)
TIDE <- round(sweep(gaby.wt,1, apply(gaby.wt, 1, median)),2)
write.table(TIDE,file.path(res.path,"Gaby_WT_TIDE_input.self_subtract"),sep = "\t",row.names = T,col.names = NA,quote = F)

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
TIDE <- round(sweep(target.wt[,com_sam],1, apply(target.wt[,com_sam], 1, mean)),2)
write.table(TIDE,file.path(res.path,"Target_WT_TIDE_input.self_subtract"),sep = "\t",row.names = T,col.names = NA,quote = F)

tide.res <- read.csv(file.path(res.path,"target_tide_output.csv"),row.names = 1,check.names = F,stringsAsFactors = F,header = T)
annCol.target114[rownames(tide.res),"TIDE"] <- tide.res$Responder
fisher.test(table(annCol.target114$TIDE,annCol.target114$tune_Clust))

# compre TLS
tmp1 <- TLS.gaby[RNA_Clust1]
tmp2 <- TLS.gaby[RNA_Clust2]
wilcox.test(tmp1,tmp2)

# compare CNA and TP53
tmp  <- read.table(file.path(data.path,"Gaby_WT_oncoprint_input.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp <- tmp[c(1,5:8),]
tmp[tmp!=""] <- 1
tmp <- as.data.frame(t(tmp))
fisher.test(table(tmp$TP53,tmp$`4q loss`)) # 0.002165
fisher.test(table(tmp$TP53,tmp$`14q loss`)) # 0.1818
fisher.test(table(tmp$TP53,tmp$`16q loss`)) # 0.2424
fisher.test(table(tmp$TP53,tmp$`22 loss`)) # 0.2424

# clinical information for 55 wts
wt3 <- read.table(file.path(data.path,"clinical data for 55 wts.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp <- data.frame(exp = c(wt3$`CD3 count`,wt3$`CD4 count`,wt3$`CD8 count`),
                  mut = rep(c("CD3","CD4","CD8"),c(55,55,55)),
                  sam = rep(rownames(wt3),3))

pdf(file.path(fig.path,"boxplot for lymphcytes of 55 wts.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data=tmp,
        xlab = "",
        ylab = "Number of cells/10 HPF",
        names = c("CD3","CD4","CD8"),
        outline = F,
        ylim = c(0,800),
        col = c(alpha("#5CC169",0.8),alpha("#410155",0.8),alpha(jco[2],0.8)))
stripchart(exp~mut,data= tmp, vertical = TRUE, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col =  c(alpha("#5CC169",0.8),alpha("#410155",0.8),alpha(jco[2],0.8)))
invisible(dev.off())

par(bty="o", mgp = c(1.5,.33,0), mar=c(3,3,1,0.1), las=1, tcl=-.25,las = 1)
barplot(wt3$`CD4 count`,col = alpha("#5CC169",0.8), border = NA)
barplot(wt3$`CD8 count`,add = T, col = alpha("#410155",0.8), border = NA)
tmp2 <- tmp
samorder <- rownames(wt3[order(wt3$`CD3 count`,decreasing = F),])
tmp2$sam <- factor(tmp2$sam,levels = samorder)
ggplot (tmp2, aes(x=sam, y=exp, fill=mut)) + 
  geom_bar (stat="identity", position = position_dodge(width = 0.8)) +   
  scale_fill_manual(values = c(alpha("#5CC169",0.8),alpha("#410155",0.8),alpha(jco[2],0.8))) +
theme_classic() + ylab("Number of cells/10 HPF") + 
  theme(axis.line.y = element_line(size = 0.8),
        axis.ticks.y = element_line(size = 0.2),
        axis.text.x = element_text(size = 8, color = "black", angle = 90),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(
    breaks = seq(0,1000,100),
    labels = seq(0,1000,100)
  )
ggsave(file.path(fig.path,"distribution of lymphcytes in 55 wts.pdf"),width = 6,height = 6)

# survival analysis
tmp <- wt3[,c("CD4:CD8 ratio","OS","OS.time","PFS","PFS.time")]
colnames(tmp)[1] <- c("ratio")
tmp$group <- ifelse(tmp$ratio > median(tmp$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = F,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0.5,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.55,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os in 55 wt cohort regarding cd4_8 ratio.pdf"), width = 4, height = 3)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = F,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0.5,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.55,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs in 55 wt cohort regarding cd4_8 ratio.pdf"), width = 4, height = 3)
print(p)
dev.off()

# multicox
tmp <- wt3[,c("CD4:CD8 ratio","OS","OS.time","PFS","PFS.time","Overall Stage (review/local)","Risk (review/local)")]
colnames(tmp)[1] <- c("ratio")
tmp$ratio <- tmp$ratio/10
tmp$group <- factor(ifelse(tmp$ratio > median(tmp$ratio),"High","Low"),levels = c("Low","High"))
tmp$stage <- ifelse(tmp$`Overall Stage (review/local)` %in% c("I","II"),"I_II","III_IV")
tmp$risk <- factor(ifelse(tmp$`Risk (review/local)` %in% c("High risk","High Risk"),"HRisk","LRisk"),levels = c("LRisk","HRisk"))
cox1 <- summary(coxph(Surv(PFS.time,PFS) ~ ratio + risk + stage, data = tmp))
cox2 <- summary(coxph(Surv(OS.time,OS) ~ ratio + risk + stage, data = tmp))

summary(coxph(Surv(PFS.time,PFS) ~ stage, data = tmp))
summary(coxph(Surv(OS.time,OS) ~ risk, data = tmp))

# forest plot
mulcox <- data.frame(variable = rownames(cox1$conf.int),
                         HR = cox1$conf.int[,1],
                         lower.95CI = cox1$conf.int[,3],
                         upper.95CI = cox1$conf.int[,4],
                         p = cox1$coefficients[,5],
                         stringsAsFactors = F)

hrtable <- rbind(c("RFS",NA,NA,NA,NA),
                 mulcox)
tabletext <- cbind(c("Variable",hrtable$variable),
                   c("HR (95% CI)",paste0(round(as.numeric(hrtable$HR),2)," (",round(as.numeric(hrtable$lower.95CI),2),"-",round(as.numeric(hrtable$upper.95CI),2),")")),
                   c("pvalue",round(as.numeric(hrtable$p),3)))
tabletext[2,2:3] <- NA

pdf(file.path(fig.path,"forestplot of risk table of 55 WTs regarding RFS.pdf"), width = 8, height = 6)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))),
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=3,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-2,-1,0,1,2,3,4),
           lwd.xaxis=2,
           xlab="log(Hazard Ratio)",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "6" = gpar(lwd=1, col="black", lty=1)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1.2,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())

mulcox <- data.frame(variable = rownames(cox2$conf.int),
                     HR = cox2$conf.int[,1],
                     lower.95CI = cox2$conf.int[,3],
                     upper.95CI = cox2$conf.int[,4],
                     p = cox2$coefficients[,5],
                     stringsAsFactors = F)

hrtable <- rbind(c("OS",NA,NA,NA,NA),
                 mulcox)
tabletext <- cbind(c("Variable",hrtable$variable),
                   c("HR (95% CI)",paste0(round(as.numeric(hrtable$HR),2)," (",round(as.numeric(hrtable$lower.95CI),2),"-",round(as.numeric(hrtable$upper.95CI),2),")")),
                   c("pvalue",round(as.numeric(hrtable$p),3)))
tabletext[2,2:3] <- NA

pdf(file.path(fig.path,"forestplot of risk table of 55 WTs regarding OS.pdf"), width = 8, height = 6)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))),
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=3,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-2,-1,0,1,2,3,4,5,6),
           lwd.xaxis=2,
           xlab="log(Hazard Ratio)",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "6" = gpar(lwd=1, col="black", lty=1)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1.2,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())

# high replication stress signature
gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}
HRS.signature <- gmt2list(file.path(comAnn.path,"high_replication_stress.gmt"))
HRS.signature <- sapply(HRS.signature, function(x) setdiff(x,""))

HRS.score <- gsva(as.matrix(gaby.wt),
                  HRS.signature,
                  method = "ssgsea")

hcs <- hclust(distanceMatrix(as.matrix(HRS.score[,rownames(annCol.wt.gaby.genomic.alter)]), "euclidean"), "ward.D2")
HRS.group.gaby <- cutree(hcs, k =2)
annCol.wt.gaby.genomic.alter$HRS <- ifelse(HRS.group.gaby[rownames(annCol.wt.gaby.genomic.alter)] == 1,"RS-Low","RS-High")
annColors.wt.gaby.genomic.alter[["HRS"]] <- c("RS-Low" = "#4558A3","RS-High" = "#C86C0D")
plotdata <- standarize.fun(HRS.score[,rownames(annCol.wt.gaby.genomic.alter)],halfwidth = 1)
p <- pheatmap(as.matrix(plotdata),
         cluster_rows = F,
         cluster_cols = hcs,
         border_color = "black",
         annotation_col = annCol.wt.gaby.genomic.alter[,c("HRS","RNA_Clust","Histology","OS","TP53 mutations","TP53 STATUS")],
         annotation_colors = annColors.wt.gaby.genomic.alter[c("HRS","RNA_Clust","Histology","OS","TP53 mutations","TP53 STATUS")],
         color = viridisLite::viridis(64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 15,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         cutree_cols = 2,
         fontsize_col = 8)
pdf(file.path(fig.path, "HRS signature heatmap of wilms tumors in gaby cohort.pdf"), height=8,width = 15)
draw(p, annotation_legend_side = "bottom", heatmap_legend_side = "left")
invisible(dev.off())
fisher.test(table(annCol.wt.gaby.genomic.alter$HRS,annCol.wt.gaby.genomic.alter$`TP53 mutations`),alternative = "less") # 0.119
fisher.test(matrix(c(4,0,1,5),byrow = T,ncol = 2)) # 0.04762
fisher.test(table(annCol.wt.gaby.genomic.alter$HRS,annCol.wt.gaby.genomic.alter$RNA_Clust)) # 0.04762

HRS.score.target <- gsva(as.matrix(target.wt),
                    HRS.signature,
                    method = "ssgsea")

hcs <- hclust(distanceMatrix(as.matrix(HRS.score.target), "euclidean"), "ward.D2")
HRS.group.target <- cutree(hcs, k =2)
annCol.target114[colnames(HRS.score.target),"HRS"] <- ifelse(HRS.group.target[colnames(HRS.score.target)] == 1,"RS-High","RS-Low")
annColors.target114[["HRS"]] <- c("RS-Low" = "#4558A3","RS-High" = "#C86C0D")

plotdata <- standarize.fun(HRS.score.target,halfwidth = 1)
p <- pheatmap(as.matrix(plotdata),
              cluster_rows = F,
              cluster_cols = dendsort(hcs),
              border_color = NA,
              annotation_col = annCol.target114[colnames(plotdata),c(48,40,4,34,1)],
              annotation_colors = annColors.target114[c("HRS","tune_Clust", "HISTOLOGY", "OS","TP53_MUT")],
              color = viridisLite::viridis(64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 5,
              cellheight = 13,
              cutree_cols = 2,
              show_colnames = F,
              show_rownames = T,
              fontsize_col = 8)
pdf(file.path(fig.path, "HRS signature heatmap of wilms tumors in target cohort.pdf"), height=10,width = 20)
draw(p, annotation_legend_side = "bottom", heatmap_legend_side = "left")
invisible(dev.off())

fisher.test(table(annCol.target114$HRS,annCol.target114$TP53_MUT),alternative = "less") # 0.04707
fisher.test(table(annCol.target114$HRS,annCol.target114$tune_Clust)) # 7.497e-05

fisher.test(table(annCol.target114$HRS,annCol.target114$HISTOLOGY)) # 0.003
#         DAWT FHWT
# RS-High   31   34
# RS-Low     5   26

fisher.test(table(annCol.wt.gaby.genomic.alter$HRS,annCol.wt.gaby.genomic.alter$Histology)) # 0.2
#         ANA diffuse Focal anaplasia
# RS-High           4               0
# RS-Low            3               3

# drug sensitivity analysis
library(pRRophetic)
library(impute)
library(SimDesign)
gdsc1 <- read.table(file.path(data.path,"GDSC1_fitted_dose_response_25Feb20.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
gdsc2 <- read.table(file.path(data.path,"GDSC2_fitted_dose_response_25Feb20.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
gdsc <- rbind.data.frame(gdsc1, gdsc2); gdsc.bk <- gdsc
gdsc <- gdsc2; gdsc.bk <- gdsc
gdsc$index <- paste(gdsc$CELL_LINE_NAME, gdsc$DRUG_NAME)
gdsc <- gdsc[!duplicated(gdsc$index),c("CELL_LINE_NAME","DRUG_NAME","LN_IC50")]
gdsc.auc <- reshape(gdsc, idvar = "CELL_LINE_NAME", timevar = "DRUG_NAME", direction = "wide")
gdsc.auc$CELL_LINE_NAME <- gsub("-","",gdsc.auc$CELL_LINE_NAME,fixed = T)
gdsc.auc <- gdsc.auc[!duplicated(gdsc.auc$CELL_LINE_NAME),]
rownames(gdsc.auc) <- gdsc.auc$CELL_LINE_NAME
gdsc.auc <- gdsc.auc[,-1]
colnames(gdsc.auc) <- gsub("LN_IC50.","",colnames(gdsc.auc),fixed = T)

load(file.path(data.path,"GDSC_EXPR.RData"))
gdsc.expr <- expr
comccl <- intersect(colnames(gdsc.expr),rownames(gdsc.auc))
gdsc.expr <- gdsc.expr[,comccl]
gdsc.auc <- gdsc.auc[comccl,]

gdsc.auc <- gdsc.auc[,apply(gdsc.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(gdsc.auc)]
gdsc.auc.knn <- as.data.frame(impute.knn(as.matrix(gdsc.auc))$data)

AT1.inhibitor <- unique(gdsc.bk[which(gdsc.bk$PUTATIVE_TARGET == "ATR"),"DRUG_NAME"])
WEE1.inhibitor <- unique(gdsc.bk[which(gdsc.bk$PUTATIVE_TARGET %in% c("WEE1, CHEK1","WEE1, PLK1")),"DRUG_NAME"])
EZH2.inhibitor <- unique(gdsc.bk[which(gdsc.bk$PUTATIVE_TARGET == "EZH2"),"DRUG_NAME"])
DHAC.inhibitor <- unique(gdsc.bk[grep("HDAC",gdsc.bk$PUTATIVE_TARGET),"DRUG_NAME"])
gdsc.auc.knn <- gdsc.auc.knn[,c(AT1.inhibitor,WEE1.inhibitor,EZH2.inhibitor,DHAC.inhibitor)]

## GDSC
trainExpr <- gdsc.expr
trainPtype <- as.data.frame(gdsc.auc.knn)
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

testExpr <- gaby.wt
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) {
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- as.vector(trainPtype[,d])
  
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = T,
                                  selection = 1))
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
gdsc.pred.auc.gaby <- as.data.frame(t(outTab))
gdsc.pred.auc.gaby$HRS <- annCol.wt.gaby.genomic.alter[rownames(gdsc.pred.auc.gaby),"HRS"]
gdsc.pred.auc.gaby$TME <- annCol.wt.gaby.genomic.alter[rownames(gdsc.pred.auc.gaby),"RNA_Clust"]
gdsc.pred.auc.gaby$TP53 <- annCol.wt.gaby.genomic.alter[rownames(gdsc.pred.auc.gaby),"TP53 mutations"]
write.table(gdsc.pred.auc.gaby,file = file.path(res.path,"predicted ln_ic50 using GDSC dataset for gaby cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-High"),5]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-Low"),5]
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TP53 == "Yes"),5]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TP53 == "No"),5]
wilcox.test(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),5]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),5]
wilcox.test(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-High"),2]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-Low"),2]
wilcox.test(tmp1,tmp2, alternative = "less")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(4,6)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD6739 ATR inhibitor in gaby WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "AZD6739\nEstimated ln(IC50)",
        main = "ATR inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2.8,"P = 0.057",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-High"),5]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-Low"),5]
wilcox.test(tmp1,tmp2, alternative = "less")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(4,6)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD1775 WEE1 inhibitor in gaby WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "AZD1775\nEstimated ln(IC50)",
        main = "WEE1 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,1.18,"P = 0.005",cex = 1)
invisible(dev.off())

## EZH2 inhibitor
tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-High"),"GSK343"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-Low"),"GSK343"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(4,6)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of GSK343 EZH2 inhibitor in gaby WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "GSK343\nEstimated ln(IC50)",
        main = "EZH2 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2.8,"P = 0.038",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),"GSK343"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),"GSK343"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of GSK343 EZH2 inhibitor in gaby WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "GSK343\nEstimated ln(IC50)",
        main = "EZH2 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,2.8,"P = 0.0008",cex = 1)
invisible(dev.off())

## HDAC inhibitor
tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-High"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-Low"),"Vorinostat"]
wilcox.test(tmp1,tmp2,alternative = "less")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(4,6)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Vorinostat HDAC inhibitor in gaby WT cohort regarding replication stress group.pdf"),width = 4,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "Vorinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,1.65,"P = 0.086",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),"Vorinostat"]
wilcox.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Vorinostat HDAC inhibitor in gaby WT cohort regarding TME.pdf"),width = 4,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "Vorinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.5,"P = 0.075",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-High"),"Entinostat"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-Low"),"Entinostat"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(4,6)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Entinostat HDAC inhibitor in gaby WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "Entinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2.6,"P = 0.61",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),"Entinostat"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),"Entinostat"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Entinostat HDAC inhibitor in gaby WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "Entinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.5,"P = 0.548",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-High"),"PCI-34051"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$HRS == "RS-Low"),"PCI-34051"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(4,6)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of PCI-34051 HDAC inhibitor in gaby WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "PCI-34051\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2.6,"P = 0.476",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),"PCI-34051"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),"PCI-34051"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of PCI-34051 HDAC inhibitor in gaby WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "PCI-34051\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.5,"P = 1",cex = 1)
invisible(dev.off())


## TARGET
testExpr <- target.wt
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) {
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- as.vector(trainPtype[,d])
  
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = T,
                                  selection = 1))
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
gdsc.pred.auc.target <- as.data.frame(t(outTab))
gdsc.pred.auc.target$HRS <- annCol.target114[rownames(gdsc.pred.auc.target),"HRS"]
gdsc.pred.auc.target$TME <- annCol.target114[rownames(gdsc.pred.auc.target),"tune_Clust"]
gdsc.pred.auc.target$TP53 <- annCol.target114[rownames(gdsc.pred.auc.target),"TP53_MUT"]
write.table(gdsc.pred.auc.target,file = file.path(res.path,"predicted ln_ic50 using GDSC dataset for target cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-High"),2]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-Low"),2]
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TP53 == "Yes"),2]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TP53 == "No"),2]
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C2"),2]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C1"),2]
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-High"),2]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-Low"),2]
wilcox.test(tmp1,tmp2, alternative = "less")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(65,31)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD6739 ATR inhibitor in target WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "AZD6739\nEstimated ln(IC50)",
        main = "ATR inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,3.5,"P = 0.003",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-High"),5]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-Low"),5]
wilcox.test(tmp1,tmp2, alternative = "less")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(65,31)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD1775 WEE1 inhibitor in target WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "AZD1775\nEstimated ln(IC50)",
        main = "WEE1 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2,"P < 0.001",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-High"),"GSK343"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-Low"),"GSK343"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(65,31)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of GSK343 EZH2 inhibitor in target WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "GSK343\nEstimated ln(IC50)",
        main = "EZH2 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,3.2,"P = 0.471",cex = 1)
invisible(dev.off())

## EZH2 inhibitor
tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C1"),"GSK343"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C2"),"GSK343"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("tTME-C1","tTME-C2"),c(69,27)))
tmp$mut <- factor(tmp$mut,levels = c("tTME-C1","tTME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of GSK343 EZH2 inhibitor in target WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "GSK343\nEstimated ln(IC50)",
        main = "EZH2 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,3.2,"P = 0.864",cex = 1)
invisible(dev.off())

## HDAC inhibitor
tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-High"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-Low"),"Vorinostat"]
wilcox.test(tmp1,tmp2,alternative = "less")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(65,31)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Vorinostat HDAC inhibitor in target WT cohort regarding replication stress group.pdf"),width = 4,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "Vorinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,1.8,"P = 0.002",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C1"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C2"),"Vorinostat"]
wilcox.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("tTME-C1","tTME-C2"),c(69,27)))
tmp$mut <- factor(tmp$mut,levels = c("tTME-C1","tTME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Vorinostat HDAC inhibitor in target WT cohort regarding TME.pdf"),width = 4,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "Vorinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.8,"P < 0.001",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-High"),"Entinostat"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-Low"),"Entinostat"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(65,31)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Entinostat HDAC inhibitor in target WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "Entinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2.6,"P = 0.51",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C1"),"Entinostat"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C2"),"Entinostat"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("tTME-C1","tTME-C2"),c(69,27)))
tmp$mut <- factor(tmp$mut,levels = c("tTME-C1","tTME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Entinostat HDAC inhibitor in target WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "Entinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.5,"P = 0.563",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-High"),"PCI-34051"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$HRS == "RS-Low"),"PCI-34051"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(65,31)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of PCI-34051 HDAC inhibitor in target WT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "PCI-34051\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,5,"P = 0.012",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C1"),"PCI-34051"]
tmp2 <- gdsc.pred.auc.target[which(gdsc.pred.auc.target$TME == "C2"),"PCI-34051"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("tTME-C1","tTME-C2"),c(69,27)))
tmp$mut <- factor(tmp$mut,levels = c("tTME-C1","tTME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of PCI-34051 HDAC inhibitor in target WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "PCI-34051\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,5,"P = 0.018",cex = 1)
invisible(dev.off())

#---------------------------#
# cell lines from preprints #
gse156065.path <- file.path(data.path,"GSE156065_RAW")
tmp1 <- read.delim(file.path(gse156065.path,"GSM4721824_COG-W-408.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp2 <- read.delim(file.path(gse156065.path,"GSM4721825_17_94.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp3 <- read.delim(file.path(gse156065.path,"GSM4721826_PCB-00007.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp4 <- read.delim(file.path(gse156065.path,"GSM4721827_CF-00108.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp5 <- read.delim(file.path(gse156065.path,"GSM4721828_CF-00136A.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp6 <- read.delim(file.path(gse156065.path,"GSM4721829_CF-00333.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp7 <- read.delim(file.path(gse156065.path,"GSM4721830_IM-WT-10.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp8 <- read.delim(file.path(gse156065.path,"GSM4721831_IM-WT-1.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp9 <- read.delim(file.path(gse156065.path,"GSM4721832_IM-WT-6.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp10 <- read.delim(file.path(gse156065.path,"GSM4721833_KT13_PDX.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp11 <- read.delim(file.path(gse156065.path,"GSM4721834_KT18_PDX.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
tmp12 <- read.delim(file.path(gse156065.path,"GSM4721835_Wit49.genes.resort.results.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)

tmp1$gene <- sapply(strsplit(tmp1$gene_id,"_",fixed = T),"[",2)
tmp1 <- as.data.frame(apply(tmp1[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp1$gene), FUN=max, na.rm=TRUE)))

tmp2$gene <- sapply(strsplit(tmp2$gene_id,"_",fixed = T),"[",2)
tmp2 <- as.data.frame(apply(tmp2[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp2$gene), FUN=max, na.rm=TRUE)))

tmp3$gene <- sapply(strsplit(tmp3$gene_id,"_",fixed = T),"[",2)
tmp3 <- as.data.frame(apply(tmp3[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp3$gene), FUN=max, na.rm=TRUE)))

tmp4$gene <- sapply(strsplit(tmp4$gene_id,"_",fixed = T),"[",2)
tmp4 <- as.data.frame(apply(tmp4[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp4$gene), FUN=max, na.rm=TRUE)))

tmp5$gene <- sapply(strsplit(tmp5$gene_id,"_",fixed = T),"[",2)
tmp5 <- as.data.frame(apply(tmp5[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp5$gene), FUN=max, na.rm=TRUE)))

tmp6 <- as.data.frame(apply(tmp6[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp6$Gene_Name), FUN=max, na.rm=TRUE)))

tmp7 <- as.data.frame(apply(tmp7[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp7$Gene_Name), FUN=max, na.rm=TRUE)))

tmp8$gene <- sapply(strsplit(tmp8$gene_id,"_",fixed = T),"[",2)
tmp8 <- as.data.frame(apply(tmp8[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp8$gene), FUN=max, na.rm=TRUE)))

tmp9$gene <- sapply(strsplit(tmp9$gene_id,"_",fixed = T),"[",2)
tmp9 <- as.data.frame(apply(tmp9[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp9$gene), FUN=max, na.rm=TRUE)))

tmp10 <- as.data.frame(apply(tmp10[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp10$Gene_Name), FUN=max, na.rm=TRUE)))

tmp11 <- as.data.frame(apply(tmp11[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp11$Gene_Name), FUN=max, na.rm=TRUE)))

tmp12$gene <- sapply(strsplit(tmp12$gene_id,"_",fixed = T),"[",2)
tmp12 <- as.data.frame(apply(tmp12[,c("expected_count","TPM","FPKM")], 2, function(x) tapply(x, INDEX=factor(tmp12$gene), FUN=max, na.rm=TRUE)))

comgene.ccl <- intersect(rownames(tmp1),rownames(tmp2))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp3))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp4))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp5))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp6))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp7))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp8))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp9))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp10))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp11))
comgene.ccl <- intersect(comgene.ccl,rownames(tmp12))

ccle.wt.ct <- data.frame("COG-W-408" = tmp1[comgene.ccl,"expected_count"],
                      "17.94" = tmp2[comgene.ccl,"expected_count"],
                      "PCB-00007" = tmp3[comgene.ccl,"expected_count"],
                      "CF-00108" = tmp4[comgene.ccl,"expected_count"],
                      "CF-00136" = tmp5[comgene.ccl,"expected_count"],
                      "CF-00333" = tmp6[comgene.ccl,"expected_count"],
                      "IM-WT-10" = tmp7[comgene.ccl,"expected_count"],
                      "IM-WT-1" = tmp8[comgene.ccl,"expected_count"],
                      "IM-WT-6" = tmp9[comgene.ccl,"expected_count"],
                      "KT-13" = tmp10[comgene.ccl,"expected_count"],
                      "KT-18" = tmp11[comgene.ccl,"expected_count"],
                      "Wit49" = tmp12[comgene.ccl,"expected_count"],
                      check.names = F,
                      row.names = comgene.ccl,
                      stringsAsFactors = F)
ccle.wt.ct <- as.data.frame(round(ccle.wt.ct,0))

ccle.wt <- data.frame("COG-W-408" = tmp1[comgene.ccl,"TPM"],
                      "17.94" = tmp2[comgene.ccl,"TPM"],
                      "PCB-00007" = tmp3[comgene.ccl,"TPM"],
                      "CF-00108" = tmp4[comgene.ccl,"TPM"],
                      "CF-00136" = tmp5[comgene.ccl,"TPM"],
                      "CF-00333" = tmp6[comgene.ccl,"TPM"],
                      "IM-WT-10" = tmp7[comgene.ccl,"TPM"],
                      "IM-WT-1" = tmp8[comgene.ccl,"TPM"],
                      "IM-WT-6" = tmp9[comgene.ccl,"TPM"],
                      "KT-13" = tmp10[comgene.ccl,"TPM"],
                      "KT-18" = tmp11[comgene.ccl,"TPM"],
                      "Wit49" = tmp12[comgene.ccl,"TPM"],
                      check.names = F,
                      row.names = comgene.ccl,
                      stringsAsFactors = F)
ccle.wt <- as.data.frame(log2(ccle.wt + 1))
ccle.wt <- ccle.wt[intersect(Ginfo[which(Ginfo$genetype == "protein_coding"),"genename"],rownames(ccle.wt)),]
ccle.wt <- ccle.wt[rowSums(ccle.wt)>0,]
write.table(ccle.wt,file = file.path(res.path,"Wilms_tumor_CCLE_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

ccle.sinfo <- data.frame(Model_Type = c("Cell line","Cell line","PDX explant","PDX explant","Cell line","Cell line","Cell line","Cell line","Primary cell culture","Primary cell culture","Primary cell culture","Primary cell culture","Primary cell culture"),
                         Histology = rep(c("DAWT","FAWT"),c(4,9)),
                         Sex = c("F","F","F","M","N/A","F","M","F","M","F","M","M","F"),
                         Age = c(2,4,3,5,9,1,1,2,1,1,24,38,5),
                         stringsAsFactors = F,
                         check.names = F,
                         row.names = c("Wit49","17.94","KT-13","KT-18","COG-W-408","IM-WT-1","IM-WT-6","IM-WT-10","Wilms-8","CF-00108","CF-00136","CF-00333","PCB-00007"))
ccle.sinfo <- ccle.sinfo[colnames(ccle.wt),]
write.table(ccle.sinfo,file = file.path(res.path,"Wilms_tumor_CCLE_SINFO.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# batch removal
modcombat = model.matrix(~1, data=annCol.ccle)
batchPCA(indata = t(scale(t(ccle.wt))),
         batch = annCol.ccle[,"Model_Type"],
         fig.dir = fig.path,
         PCA.fig.title = "PCA ccle before combat",
         cols = annColors.ccle$Model_Type,
         showID = F,
         cex = 0.9,
         showLegend = T)

ccle.wt.combat = ComBat(dat=as.matrix(ccle.wt), batch=annCol.ccle[,"Model_Type"], mod=modcombat)
write.table(ccle.wt.combat,file = file.path(res.path,"Wilms_tumor_CCLE_log2TPM_hugo_combat.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

batchPCA(indata = t(scale(t(ccle.wt.combat))),
         batch = annCol.ccle[,"Model_Type"],
         fig.dir = fig.path,
         PCA.fig.title = "PCA ccle after combat",
         cols = annColors.ccle$Model_Type,
         showID = F,
         cex = 0.9,
         showLegend = T)

annCol.ccle <- ccle.sinfo
annCol.ccle$Age <- ifelse(annCol.ccle$Age > median(annCol.ccle$Age),">3.5","<=3.5")
annColors.ccle <- list("Model_Type" = c("Cell line" = "black","PDX explant" = "grey40", "Primary cell culture" = "grey80"),
                       "Histology" = c("DAWT" = "#21498D","FAWT" = "#E53435"),
                       "Sex" = c("F" = jco[2],"M" = jco[1],"N/A" = "white"),
                       "Age" = c(">3.5" = orange,"<=3.5" = alpha(orange,0.6)))
annCol.ccle2 <- annCol.ccle[which(annCol.ccle$Model_Type == "Cell line"),]

# supervised clustering using pathways
paths <- c("HALLMARK_P53_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA","REACTOME_HDACS_DEACETYLATE_HISTONES","REACTOME_HOMOLOGY_DIRECTED_REPAIR")
paths <- MSigDB[which(MSigDB$ont %in% paths),]
paths.sig <- list()
for (i in unique(paths$ont)) {
  paths.sig[[i]] <- intersect(toupper(paths[which(paths$ont == i),"gene"]),rownames(ccle.wt))
}
paths.gsva.ccle <- gsva(as.matrix(ccle.wt),
                        paths.sig,
                        method = "ssgsea")
paths.gsva.ccle <- paths.gsva.ccle[c("HALLMARK_P53_PATHWAY",
                                     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                     "REACTOME_HDACS_DEACETYLATE_HISTONES",
                                     "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
                                     "REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA"),]
hcs <- hclust(distanceMatrix(as.matrix(paths.gsva.ccle[,rownames(annCol.ccle)]), "euclidean"), "ward.D"); hcs.ccle.tme <- hcs
group <- cutree(hcs, k =2)
annCol.ccle$TME <- ifelse(group[rownames(annCol.ccle)] == 1,"cTME-C1","cTME-C2")
annCol.ccle$EZH2_targets <- as.numeric(scale(ezh2.path.ccle[1,rownames(annCol.ccle)]))
annColors.ccle[["TME"]] <- c("cTME-C1" = jco[2],"cTME-C2" = jco[1])
annColors.ccle[["EZH2_targets"]] <- viridisLite::inferno(64)
plotdata <- standarize.fun(paths.gsva.ccle[,rownames(annCol.ccle)],halfwidth = 1)
p1 <- pheatmap(as.matrix(plotdata),
              cluster_rows = F,
              cluster_cols = hcs.ccle.tme,
              border_color = "black",
              annotation_col = annCol.ccle[,c("TME","Histology","Age","Sex","Model_Type")],
              annotation_colors = annColors.ccle[c("TME","Histology","Age","Sex","Model_Type")],
              color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
              #color = viridisLite::inferno(64),
              #color = viridisLite::viridis(64),
              # color = greenred(64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 12,
              cellheight = 12,
              gaps_row = 3,
              show_colnames = T,
              show_rownames = T,
              name = "ssgsea\n[z-scored]",
              cutree_cols = 2,
              fontsize_col = 9)
p1

# clustering using pathway genes
outTab <- NULL
for (i in names(paths.sig)) {
  tmp <- paths.sig[[i]]
  for (j in tmp) {
    tmp1 <- as.numeric(ccle.wt[j,RNA_Clust1.ccle])
    tmp2 <- as.numeric(ccle.wt[j,RNA_Clust2.ccle])
    outTab <- rbind.data.frame(outTab,
                               data.frame(path = i,
                                          gene = j,
                                          avgTME1 = mean(tmp1),
                                          avgTME2 = mean(tmp2),
                                          diff = abs(mean(tmp1)-mean(tmp2)),
                                          p = wilcox.test(tmp1,tmp2)$p.value,
                                          dirct = ifelse(mean(tmp1) > mean(tmp2),"TME1","TME2"),
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
  }
}
outTab <- outTab[order(outTab$diff,decreasing = T),]
selgene1 <- outTab[which(outTab$p < 0.05 & outTab$dirct == "TME1" & outTab$path == "HALLMARK_P53_PATHWAY"),][1:5,]
selgene2 <- outTab[which(outTab$p < 0.05 & outTab$dirct == "TME1" & outTab$path == "HALLMARK_INTERFERON_GAMMA_RESPONSE"),][1:5,]
selgene3 <- outTab[which(outTab$p < 0.05 & outTab$dirct == "TME1" & outTab$path == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),][1:5,]
selgene4 <- outTab[which(outTab$p < 0.5 & outTab$dirct == "TME2" & outTab$path == "REACTOME_HDACS_DEACETYLATE_HISTONES"),][1:5,]
selgene5 <- outTab[which(outTab$p < 0.5 & outTab$dirct == "TME2" & outTab$path == "REACTOME_HOMOLOGY_DIRECTED_REPAIR"),][1:5,]
selgene6 <- outTab[which(outTab$p < 0.5 & outTab$dirct == "TME2" & outTab$path == "REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA"),][1:5,]

plotdata <- ccle.wt[c(selgene1$gene,
                      selgene2$gene,
                      selgene3$gene,
                      selgene4$gene,
                      selgene5$gene,
                      selgene6$gene),]
annRow <- data.frame(class = rep(c("HALLMARK_P53_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                   "REACTOME_HDACS_DEACETYLATE_HISTONES","REACTOME_HOMOLOGY_DIRECTED_REPAIR","REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA"),
                                 c(5,5,5,5,5,5)),
                     row.names = rownames(plotdata))
annColors.ccle[["class"]] = c("HALLMARK_P53_PATHWAY" = darkblue,
                         "HALLMARK_INTERFERON_GAMMA_RESPONSE" = alpha(darkblue,0.6),
                         "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" =alpha(darkblue,0.3),
                         "REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA" = darkred,
                         "REACTOME_HOMOLOGY_DIRECTED_REPAIR" = alpha(darkred,0.6),
                         "REACTOME_HDACS_DEACETYLATE_HISTONES" = alpha(darkred,0.3))

plotdata <- standarize.fun(plotdata,halfwidth = 2)
p2 <- pheatmap(as.matrix(plotdata),
              cluster_rows = F,
              cluster_cols = hcs.ccle.tme,
              border_color = "black",
              annotation_row = annRow,
              annotation_colors = annColors.ccle["class"],
              color = colorpanel(64,low=blue,mid = "black",high=gold),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 12,
              cellheight = 12,
              gaps_row = c(5,10,15,20,25),
              show_colnames = T,
              show_rownames = T,
              name = "Expr.\n[z-scored]",
              cutree_cols = 2,
              fontsize_col = 9)
pdf(file.path(fig.path, "supervised clustering using msigdb pathways on wilms tumors in ccle cohort before combat.pdf"), height=12,width = 10)
draw(p1 %v% p2, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
invisible(dev.off())


# check EZH2
RNA_Clust1.ccle <- rownames(annCol.ccle[which(annCol.ccle$TME == "cTME-C1"),])
RNA_Clust2.ccle <- rownames(annCol.ccle[which(annCol.ccle$TME == "cTME-C2"),])

ezh2.ccle.c1 <- as.numeric(ccle.wt["EZH2",RNA_Clust1.ccle])
ezh2.ccle.c2 <- as.numeric(ccle.wt["EZH2",RNA_Clust2.ccle])
wilcox.test(ezh2.ccle.c1,ezh2.ccle.c2) #0.5697
tmp <- data.frame(exp = c(ezh2.ccle.c1,ezh2.ccle.c2),
                  mut = rep(c("cTME-C1","cTME-C2"),c(4,8)))
tmp$mut <- factor(tmp$mut,levels = c("cTME-C1","cTME-C2"))
pdf(file.path(fig.path,"boxplot for EZH2 expression in ccle WT cohort regarding TME cluster.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "log2(TPM) of EZH2",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,6.5,"P = 0.57",cex = 1)
invisible(dev.off())

# EZH2 pathway
tmp <- MSigDB[grep("KAMMINGA_EZH2_TARGETS",MSigDB$ont),]
ezh2.path.ccle <- gsva(expr = as.matrix(ccle.wt),
                       gset.idx.list = list("KAMMINGA_EZH2_ccleS" = tmp$gene),
                       method = "ssgsea")
ezh2.path.ccle.c1 <- as.numeric(ezh2.path.ccle[1,RNA_Clust1.ccle])
ezh2.path.ccle.c2 <- as.numeric(ezh2.path.ccle[1,RNA_Clust2.ccle])
wilcox.test(ezh2.path.ccle.c1,ezh2.path.ccle.c2) # 0.008
tmp <- data.frame(exp = c(ezh2.path.ccle.c1,ezh2.path.ccle.c2),
                  mut = rep(c("cTME-C1","cTME-C2"),c(8,4)))
tmp$mut <- factor(tmp$mut,levels = c("cTME-C1","cTME-C2"))
pdf(file.path(fig.path,"boxplot for EZH2 target in ccle WT cohort regarding TME cluster.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "EZH2 targets\nssGSEA score",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,1.1,"P = 0.008",cex = 1)
invisible(dev.off())

# supervised clustering using RS signature
HRS.score.ccle <- gsva(as.matrix(ccle.wt),
                       HRS.signature,
                       method = "ssgsea")
hcs <- hclust(distanceMatrix(as.matrix(HRS.score.ccle[,rownames(annCol.ccle)]), "pearson"), "ward.D"); hcs.ccle.hrs <- hcs
HRS.group.ccle <- cutree(hcs, k =2)
annCol.ccle$HRS <- ifelse(HRS.group.ccle[rownames(annCol.ccle)] == 1,"RS-Low","RS-High")
annColors.ccle[["HRS"]] <- c("RS-Low" = "#4558A3","RS-High" = "#C86C0D")
plotdata <- standarize.fun(HRS.score.ccle[,rownames(annCol.ccle)],halfwidth = 1)
p <- pheatmap(as.matrix(plotdata),
              cluster_rows = F,
              cluster_cols = hcs,
              border_color = "black",
              annotation_col = annCol.ccle[,c("HRS","TME","Histology","Age","Sex","Model_Type")],
              annotation_colors = annColors.ccle[c("HRS","TME","Histology","Age","Sex","Model_Type")],
              color = viridisLite::viridis(64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 15,
              cellheight = 10,
              show_colnames = T,
              show_rownames = T,
              cutree_cols = 2,
              fontsize_col = 9)
pdf(file.path(fig.path, "HRS signature heatmap of wilms tumors in ccle cohort before combat.pdf"), height=8,width = 15)
draw(p, annotation_legend_side = "bottom", heatmap_legend_side = "left")
invisible(dev.off())

table(annCol.ccle$TME,annCol.ccle$HRS)
fisher.test(table(annCol.ccle$TME,annCol.ccle$HRS)) # 0.01

# drug sensitivity
trainExpr <- gdsc.expr
trainPtype <- as.data.frame(gdsc.auc.knn)
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

testExpr <- ccle.wt
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) {
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- as.vector(trainPtype[,d])
  
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = T,
                                  selection = 1))
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
gdsc.pred.auc.ccle <- as.data.frame(t(outTab))
gdsc.pred.auc.ccle$HRS <- annCol.ccle[rownames(gdsc.pred.auc.ccle),"HRS"]
gdsc.pred.auc.ccle$TME <- annCol.ccle[rownames(gdsc.pred.auc.ccle),"TME"]
write.table(gdsc.pred.auc.ccle,file = file.path(res.path,"predicted ln_ic50 using GDSC dataset for ccle cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

tmp1 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$HRS == "RS-High"),2]
tmp2 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$HRS == "RS-Low"),2]
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),2]
tmp2 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C1"),2]
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"GSK343"]
tmp2 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"GSK343"]
wilcox.test(tmp1,tmp2)
boxplot(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"Vorinostat"]
wilcox.test(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"Entinostat"]
tmp2 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"Entinostat"]
wilcox.test(tmp1,tmp2)

tmp1 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"PCI-34051"]
tmp2 <- gdsc.pred.auc.ccle[which(gdsc.pred.auc.ccle$TME == "cTME-C2"),"PCI-34051"]
wilcox.test(tmp1,tmp2)

#------------------------------------#
# drug sensitivity using CTRP cohort #
library(impute)
library(SimDesign)
library(pRRophetic)
ctrp.ccl.anno <- read.delim(file.path(data.path,"CTRP_ccl_anno.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 1
ctrp.cpd.anno <- read.delim(file.path(data.path,"CTRP_cpd_anno.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 2

ctrp.auc <- read.delim(file.path(data.path,"CTRP_AUC.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
prism.auc <- read.delim(file.path(data.path,"PRISM_AUC.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 数据来自https://depmap.org/portal/download/ Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4
prism.ccl.anno <- prism.auc[,1:5]
prism.auc <- prism.auc[,-c(1:5)]

ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.auc)]
prism.auc <- prism.auc[,apply(prism.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(prism.auc)]

rmccl <- paste0("CCL",na.omit(ctrp.ccl.anno[which(ctrp.ccl.anno$ccle_primary_site == "haematopoietic_and_lymphoid_tissue"),"master_ccl_id"]))
rownames(ctrp.auc) <- paste0("CCL",rownames(ctrp.auc))
ctrp.auc <- ctrp.auc[setdiff(rownames(ctrp.auc),rmccl),]

ctrp.auc.knn <- impute.knn(as.matrix(ctrp.auc))$data
prism.auc.knn <- impute.knn(as.matrix(prism.auc))$data

ctrp.auc.knn <- ctrp.auc.knn/ceiling(max(ctrp.auc.knn))
prism.auc.knn <- prism.auc.knn/ceiling(max(prism.auc.knn))

ccl.expr <- read.table(file.path(data.path,"CCLE_RNAseq_rsem_genes_tpm_20180929.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 
ccl.expr <- ccl.expr[,-1]; rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)
Ginfo2 <- Ginfo
rownames(Ginfo2) <- sapply(strsplit(rownames(Ginfo2),".",fixed = T),"[",1)
comgene <- intersect(rownames(ccl.expr),rownames(Ginfo2))
ccl.expr <- ccl.expr[comgene,]
ccl.expr$gene <- Ginfo2[comgene,"genename"]; ccl.expr <- ccl.expr[!duplicated(ccl.expr$gene),]; rownames(ccl.expr) <- ccl.expr$gene; ccl.expr <- ccl.expr[,-ncol(ccl.expr)]

## drug sensitivity
comgene <- intersect(rownames(gaby.wt),rownames(target.wt))
comgene <- intersect(comgene,rownames(ccle.wt))
comgene <- intersect(comgene,rownames(ccl.expr))

ctrp.drug <- data.frame(
  drug = c("ETP-46464", # ATR
           "MK-1775", # WEE1
           "BRD-K42260513","BRD-K51831558","BRD1835", # EZH2
           "BRD-A94377914","tubastatin A","panobinostat","Merck60","BRD-K11533227","belinostat", # HDAC
           "BRD-K24690302","BRD-K29313308","pandacostat","Repligen 136","BRD-K51490254",
           "tacedinaline","isonicotinohydroxamic acid","apicidin","BRD-K66532283","ISOX",
           "entinostat","BRD-K80183349","vorinostat","BRD-K85133207","BRD-K88742110"),
  target = c("ATR",
             "WEE1",
             rep("EZH2",3),
             rep("HDAC",21)),
  stringsAsFactors = F,
  check.names = F
)


## CTRP
trainExpr <- log2(ccl.expr[comgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) 
trainPtype <- as.data.frame(ctrp.auc.knn)
ccl.name <- ccl.miss <- c()
for (i in rownames(trainPtype)) {
  if(!is.element(gsub("CCL","",i),ctrp.ccl.anno$master_ccl_id)) {
    cat(i,"\n")
    ccl.miss <- c(ccl.miss, i)
    ccl.name <- c(ccl.name, i)
  } else {
    ccl.name <- c(ccl.name,  ctrp.ccl.anno[which(ctrp.ccl.anno$master_ccl_id == gsub("CCL","",i)),"ccl_name"])
  }
}

cpd.name <- cpd.miss <- c()
for (i in colnames(trainPtype)) {
  if(!is.element(i,ctrp.cpd.anno$master_cpd_id)) {
    cat(i,"\n")
    cpd.miss <- c(cpd.miss, i)
    cpd.name <- c(cpd.name, i)
  } else {
    cpd.name <- c(cpd.name,  ctrp.cpd.anno[which(ctrp.cpd.anno$master_cpd_id == i),"cpd_name"]) 
  }
}

rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),]
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)]
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,intersect(ctrp.drug$drug,colnames(trainPtype))]

### Gaby
testExpr <- gaby.wt
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) {
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.001) 
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.001
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc.gaby <- as.data.frame(t(outTab))
ctrp.pred.auc.gaby$HRS <- annCol.wt.gaby.genomic.alter[rownames(ctrp.pred.auc.gaby),"HRS"]
ctrp.pred.auc.gaby$TME <- annCol.wt.gaby.genomic.alter[rownames(ctrp.pred.auc.gaby),"RNA_Clust"]
ctrp.pred.auc.gaby$TP53 <- annCol.wt.gaby.genomic.alter[rownames(ctrp.pred.auc.gaby),"TP53 mutations"]
write.table(ctrp.pred.auc.gaby,file = file.path(res.path,"predicted auc using CTRP dataset in gaby cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

outTab.TME <- outTab.HRS <- NULL
for (i in intersect(ctrp.drug$drug,colnames(ctrp.pred.auc.gaby))) {
  tmp <- ctrp.pred.auc.gaby[,c(i,"TME","HRS")]
  tmp1 <- tmp[which(tmp$TME == "TP53_Wild"),i]
  tmp2 <- tmp[which(tmp$TME == "TP53_Mutated"),i]
  wt <- wilcox.test(tmp1,tmp2)
  
  outTab.TME <- rbind.data.frame(outTab.TME,
                             data.frame(drug = i,
                                        class = ctrp.drug[which(ctrp.drug$drug == i),"target"],
                                        avgTME1 = mean(tmp1),
                                        avgTME2 = mean(tmp2),
                                        sensitivity = ifelse(mean(tmp1) > mean(tmp2),"TME2","TME1"),
                                        p.TME = wt$p.value),
                             stringsAsFactors = F)
  
  tmp1 <- tmp[which(tmp$HRS == "RS-High"),i]
  tmp2 <- tmp[which(tmp$HRS == "RS-Low"),i]
  wt <- wilcox.test(tmp1,tmp2)
  outTab.HRS <- rbind.data.frame(outTab.HRS,
                             data.frame(drug = i,
                                        class = ctrp.drug[which(ctrp.drug$drug == i),"target"],
                                        avgRS.High = mean(tmp1),
                                        avgRS.Low = mean(tmp2),
                                        sensitivity = ifelse(mean(tmp1) > mean(tmp2),"RS.Low","RS.High"),
                                        p.HRS = wt$p.value),
                             stringsAsFactors = F)
}
outTab.TME[which(outTab.TME$p.TME >= 0.05), "sensitivity"] <- ""
outTab.HRS[which(outTab.HRS$p.HRS >= 0.05), "sensitivity"] <- ""
write.table(outTab.TME,file.path(res.path,"CTRP selected drug sensitivity regarding TME in gaby cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab.HRS,file.path(res.path,"CTRP selected drug sensitivity regarding HRS in gaby cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

## TARGET
testExpr <- target.wt
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.001)
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.001 
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc.target <- as.data.frame(t(outTab))
ctrp.pred.auc.target$HRS <- annCol.target114[rownames(ctrp.pred.auc.target),"HRS"]
ctrp.pred.auc.target$TME <- annCol.target114[rownames(ctrp.pred.auc.target),"tune_Clust"]
ctrp.pred.auc.target$TP53 <- annCol.target114[rownames(ctrp.pred.auc.target),"TP53_MUT"]
write.table(ctrp.pred.auc.target,file = file.path(res.path,"predicted auc using CTRP dataset in target cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

outTab.TME <- outTab.HRS <- NULL
for (i in intersect(ctrp.drug$drug,colnames(ctrp.pred.auc.target))) {
  tmp <- ctrp.pred.auc.target[,c(i,"TME","HRS")]
  tmp1 <- tmp[which(tmp$TME == "C1"),i]
  tmp2 <- tmp[which(tmp$TME == "C2"),i]
  wt <- wilcox.test(tmp1,tmp2)
  
  outTab.TME <- rbind.data.frame(outTab.TME,
                                 data.frame(drug = i,
                                            class = ctrp.drug[which(ctrp.drug$drug == i),"target"],
                                            avgTME1 = mean(tmp1),
                                            avgTME2 = mean(tmp2),
                                            sensitivity = ifelse(mean(tmp1) > mean(tmp2),"TME2","TME1"),
                                            p.TME = wt$p.value),
                                 stringsAsFactors = F)
  
  tmp1 <- tmp[which(tmp$HRS == "RS-High"),i]
  tmp2 <- tmp[which(tmp$HRS == "RS-Low"),i]
  wt <- wilcox.test(tmp1,tmp2)
  outTab.HRS <- rbind.data.frame(outTab.HRS,
                                 data.frame(drug = i,
                                            class = ctrp.drug[which(ctrp.drug$drug == i),"target"],
                                            avgRS.High = mean(tmp1),
                                            avgRS.Low = mean(tmp2),
                                            sensitivity = ifelse(mean(tmp1) > mean(tmp2),"RS.Low","RS.High"),
                                            p.HRS = wt$p.value),
                                 stringsAsFactors = F)
}
outTab.TME[which(outTab.TME$p.TME >= 0.05), "sensitivity"] <- ""
outTab.HRS[which(outTab.HRS$p.HRS >= 0.05), "sensitivity"] <- ""
write.table(outTab.TME,file.path(res.path,"CTRP selected drug sensitivity regarding TME in target cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab.HRS,file.path(res.path,"CTRP selected drug sensitivity regarding HRS in target cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

## CCLE
testExpr <- ccle.wt
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.001) 
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.001 
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc.ccle <- as.data.frame(t(outTab))
ctrp.pred.auc.ccle$HRS <- annCol.ccle[rownames(ctrp.pred.auc.ccle),"HRS"]
ctrp.pred.auc.ccle$TME <- annCol.ccle[rownames(ctrp.pred.auc.ccle),"TME"]
write.table(ctrp.pred.auc.ccle,file = file.path(res.path,"predicted auc using CTRP dataset in ccle cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

outTab.TME <- outTab.HRS <- NULL
for (i in intersect(ctrp.drug$drug,colnames(ctrp.pred.auc.ccle))) {
  tmp <- ctrp.pred.auc.ccle[,c(i,"TME","HRS")]
  tmp1 <- tmp[which(tmp$TME == "cTME-C1"),i]
  tmp2 <- tmp[which(tmp$TME == "cTME-C2"),i]
  wt <- wilcox.test(tmp1,tmp2,alternative = "greater")
  
  outTab.TME <- rbind.data.frame(outTab.TME,
                                 data.frame(drug = i,
                                            class = ctrp.drug[which(ctrp.drug$drug == i),"target"],
                                            avgTME1 = mean(tmp1),
                                            avgTME2 = mean(tmp2),
                                            sensitivity = ifelse(mean(tmp1) > mean(tmp2),"TME2","TME1"),
                                            p.TME = wt$p.value),
                                 stringsAsFactors = F)
  
  tmp1 <- tmp[which(tmp$HRS == "RS-High"),i]
  tmp2 <- tmp[which(tmp$HRS == "RS-Low"),i]
  wt <- wilcox.test(tmp1,tmp2,alternative = "less")
  outTab.HRS <- rbind.data.frame(outTab.HRS,
                                 data.frame(drug = i,
                                            class = ctrp.drug[which(ctrp.drug$drug == i),"target"],
                                            avgRS.High = mean(tmp1),
                                            avgRS.Low = mean(tmp2),
                                            sensitivity = ifelse(mean(tmp1) > mean(tmp2),"RS.Low","RS.High"),
                                            p.HRS = wt$p.value),
                                 stringsAsFactors = F)
}
outTab.TME[which(outTab.TME$p.TME >= 0.05), "sensitivity"] <- ""
outTab.HRS[which(outTab.HRS$p.HRS >= 0.05), "sensitivity"] <- ""
write.table(outTab.TME,file.path(res.path,"CTRP selected drug sensitivity regarding TME in ccle cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab.HRS,file.path(res.path,"CTRP selected drug sensitivity regarding HRS in ccle cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

#----------------------------------------------------------#
# copy number analysis to compare the chromosome landscape #
## Gaby cohort
Segment <- read.table(file.path(data.path,"aWT-CNVSum-GM-01-Seg.txt"),sep="\t",header=T,stringsAsFactors = F)
colnames(Segment)<- c("Sample","Chrom","Start","Stop","Mark","Seg.CN")
Segment$Sample <- gsub("PED","PEDrna",Segment$Sample)
Segment1 <- Segment[which(Segment$Sample %in% RNA_Clust1),]
Segment2 <- Segment[which(Segment$Sample %in% RNA_Clust2),]
unique(Segment1$Sample)
unique(Segment2$Sample)

write.table(Segment1,file.path(res.path,"Gaby_wilms_TME1_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(Segment2,file.path(res.path,"Gaby_wilms_TME2_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)

marker <- Segment1[,1:4]
head(marker)

# create marker file for GISTIC2.0
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chrom"],2))
  c <- c(c,marker[i,"Start"],marker[i,"Stop"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"Gaby_wilms_TME1_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- Segment2[,1:4]
head(marker)

# create marker file for GISTIC2.0
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chrom"],2))
  c <- c(c,marker[i,"Start"],marker[i,"Stop"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"Gaby_wilms_TME2_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

## TARGET cohort
seg <- read.table(file.path(data.path,"wt_target_2018_pub/data_cna_hg19.seg"),sep="\t",header=T,stringsAsFactors = F)
colnames(seg)<- c("ID","Chr","Start","End","Probes","Log2Ratio")
head(seg)

seg$ID <- paste0(seg$ID,"A")
com_sam <- intersect(unique(seg$ID),rownames(matchID))
seg <- seg[which(seg$ID %in% com_sam),]
seg$SRA_RUN <- matchID[match(seg$ID,rownames(matchID)),"SRA_RUN"]

tmp1 <- seg[which(seg$SRA_RUN %in% RNA_Clust1.target),] # 67
tmp2 <- seg[which(seg$SRA_RUN %in% RNA_Clust2.target),] # 27

tmp1$ID <- tmp1$SRA_RUN
tmp2$ID <- tmp2$SRA_RUN

Segment <- rbind.data.frame(tmp1,tmp2)
Segment <- Segment[,1:6]
colnames(Segment)<- c("Sample","Chrom","Start","Stop","Mark","Seg.CN")
Segment1 <- Segment[which(Segment$Sample %in% RNA_Clust1.target),]
Segment2 <- Segment[which(Segment$Sample %in% RNA_Clust2.target),]
unique(Segment1$Sample)
unique(Segment2$Sample)
write.table(Segment1,file.path(res.path,"TARGET_wilms_TME1_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(Segment2,file.path(res.path,"TARGET_wilms_TME2_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)

marker <- Segment1[,1:4]
head(marker)

# create marker file for GISTIC2.0
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chrom"],2))
  c <- c(c,marker[i,"Start"],marker[i,"Stop"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"TARGET_wilms_TME1_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- Segment2[,1:4]
head(marker)

# create marker file for GISTIC2.0
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chrom"],2))
  c <- c(c,marker[i,"Start"],marker[i,"Stop"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"TARGET_wilms_TME2_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

#-------------------------------#
# generate chromosome landscape #
chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}

# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)

scores <- read.table(file.path(data.path,"GISTIC/Gaby_Wilms_TME1/Gaby_wilms_TME1_0.2_0.75.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)
# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-5,max(scores$frequency)+5)
title="SFCE-WT cohort: TME-C1 (n = 5)"

pdf(file.path(fig.path,"gistic frequency TME in gaby cohort.pdf"), width = 10,height = 8)
par(mfrow = c(2,1))
plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 1, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "Frequency", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))

scores <- read.table(file.path(data.path,"GISTIC/Gaby_Wilms_TME2/Gaby_wilms_TME2_0.2_0.75.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)
# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-5,max(scores$frequency)+5)
title="SFCE-WT cohort: TME-C2 (n = 5)"

#pdf(file.path(fig.path,"gistic frequency TME-C2 in gaby cohort.pdf"), width = 10,height = 4)
plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 1, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "Frequency", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
invisible(dev.off())

scores <- read.table(file.path(data.path,"GISTIC/TARGET_Wilms_TME1/TARGET_wilms_TME1_0.2_0.75.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)
# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-5,max(scores$frequency)+5)
title="TARGET-WT cohort: tTME-C1 (n = 67)"

pdf(file.path(fig.path,"gistic frequency TME in target cohort.pdf"), width = 10,height = 8)
par(mfrow = c(2,1))
plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 1, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "Frequency", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))

scores <- read.table(file.path(data.path,"GISTIC/TARGET_Wilms_TME2/TARGET_wilms_TME2_0.2_0.75.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)
# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-5,max(scores$frequency)+5)
title="TARGET-WT cohort: tTME-C2 (n = 27)"

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 1, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "Frequency", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
invisible(dev.off())

# check arm-level CNA difference
gaby.cna.arm <- read.table(file.path(data.path,"GISTIC/Gaby_Wilms_all/Gaby_wilms_0.2_0.75.broad_values_by_arm.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gaby.cna.arm <- as.data.frame(t(gaby.cna.arm))
# gaby.cna.arm[gaby.cna.arm < 0] <- -1
# gaby.cna.arm[gaby.cna.arm > 0] <- 1
gaby.cna.arm[gaby.cna.arm != 0] <- "1"
gaby.cna.arm$TP53 <- mut.matrix["TP53",rownames(gaby.cna.arm)]
gaby.cna.arm$TP53 <- ifelse(gaby.cna.arm$TP53 == "","Wild","Mutant")

p <- c()
for (i in setdiff(colnames(gaby.cna.arm),"TP53")) {
  tmp <- gaby.cna.arm[,c(i,"TP53")]
  p <- c(p,fisher.test(table(tmp[,i],tmp[,2]),alternative = "less")$p.value)
}
names(p) <- setdiff(colnames(gaby.cna.arm),"TP53")
p
# 1p         1q         7p         7q        12p        14q        16q        17p        17q        18q        19q        22q 
# 1.00000000 1.00000000 0.18181818 1.00000000 1.00000000 0.06060606 1.00000000 0.01515152 0.08008658 0.45454545 0.06060606 0.24242424

# check common arm-level CNA regarding TME cluster
gaby.cna.arm <- read.table(file.path(data.path,"GISTIC/Gaby_Wilms_all/Gaby_wilms_0.2_0.75.broad_values_by_arm.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gaby.cna.arm <- as.data.frame(t(gaby.cna.arm))
gaby.cna.arm[gaby.cna.arm > 0] <- "GAIN"
gaby.cna.arm[gaby.cna.arm < 0] <- "LOSS"
gaby.cna.arm[gaby.cna.arm == 0] <- "NORMAL"
gaby.cna.arm <- gaby.cna.arm[rownames(annCol.wt.gaby.genomic.alter),]
gaby.cna.arm$TME <- ifelse(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Wild","TME-C1","TME-C2")

target.cna.arm <- read.table(file.path(data.path,"GISTIC/TARGET_Wilms_all/TARGET_wilms_all_0.2_0.75.broad_values_by_arm.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(target.cna.arm) <- paste0(colnames(target.cna.arm),"A")
com_sam <- intersect(colnames(target.cna.arm),rownames(matchID))
target.cna.arm <- target.cna.arm[,which(colnames(target.cna.arm) %in% com_sam)]
colnames(target.cna.arm) <- matchID[match(colnames(target.cna.arm),rownames(matchID)),"SRA_RUN"]

tmp1 <- target.cna.arm[,which(colnames(target.cna.arm) %in% RNA_Clust1.target)] # 67
tmp2 <- target.cna.arm[,which(colnames(target.cna.arm) %in% RNA_Clust2.target)] # 27
target.cna.arm <- cbind.data.frame(tmp1,tmp2)
target.cna.arm <- as.data.frame(t(target.cna.arm))
target.cna.arm[target.cna.arm > 0] <- "GAIN"
target.cna.arm[target.cna.arm < 0] <- "LOSS"
target.cna.arm[target.cna.arm == 0] <- "NORMAL"
target.cna.arm$TME <- rep(c("TME-C1","TME-C2"),c(67,27))

fisher.test(table(gaby.cna.arm$`1q`,gaby.cna.arm$TME))
fisher.test(table(gaby.cna.arm$`14q`,gaby.cna.arm$TME)) # 0.04762
# TME-C1 TME-C2
# LOSS        0      4
# NORMAL      5      1
fisher.test(table(gaby.cna.arm$`16q`,gaby.cna.arm$TME))
fisher.test(table(gaby.cna.arm$`17q`,gaby.cna.arm$TME))
fisher.test(table(gaby.cna.arm$`17p`,gaby.cna.arm$TME))
fisher.test(table(gaby.cna.arm$`22q`,gaby.cna.arm$TME))

table(target.cna.arm$`14q`,target.cna.arm$TME)
# TME-C1 TME-C2
# GAIN        1      0
# LOSS        7      7
# NORMAL     59     20
fisher.test(matrix(c(7,7,60,20),byrow = T,ncol = 2),alternative = "less") # 0.06

table(target.cna.arm$`17q`,target.cna.arm$TME)
fisher.test(matrix(c(8,4,59,23),byrow = T,ncol = 2),alternative = "less")

# check focal-level CNA regarding TME cluster
gaby.cna.focal <- read.table(file.path(data.path,"GISTIC/Gaby_Wilms_all/Gaby_wilms_0.2_0.75.all_lesions.conf_75.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gaby.cna.focal <- gaby.cna.focal[2:3,c(1,9:20)]
rownames(gaby.cna.focal) <- gaby.cna.focal$Descriptor; gaby.cna.focal <- gaby.cna.focal[,-1]
gaby.cna.focal <- as.data.frame(t(gaby.cna.focal))
gaby.cna.focal <- gaby.cna.focal[rownames(annCol.wt.gaby.genomic.alter),]
gaby.cna.focal$TME <- ifelse(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Wild","TME-C1","TME-C2")

fisher.test(table(gaby.cna.focal$`11q14.3`,gaby.cna.focal$TME))
fisher.test(table(gaby.cna.focal$`16q22.1`,gaby.cna.focal$TME))

#------------------#
# regulon analysis #
write.table(target.wt[,com_sam],"D:/Project/Wilms/Regulon/Wilms_tumor_TARGET_log2TPM_hugo96.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
load(file.path(res.path,"chromatinRemodeling_regact_gaby.RData"))
load(file.path(res.path,"chromatinRemodeling_regact_target.RData"))
load(file.path(res.path,"chromatinRemodeling_regact_ccle.RData"))

regulon.gaby <- chromatinRemodeling_regact_gaby$differential
regulon.target <- chromatinRemodeling_regact_target$differential
regulon.ccle <- chromatinRemodeling_regact_ccle$differential

comregulon <- intersect(colnames(regulon.gaby),colnames(regulon.target))
comregulon <- intersect(comregulon,colnames(regulon.ccle))

regulon.gaby.sel <- regulon.gaby[,comregulon]
regulon.target.sel <- regulon.target[,comregulon]
regulon.ccle.sel <- regulon.ccle[,comregulon]

pdf(file.path(fig.path, "common chromotin remodeling regulon of wilms tumors in gaby cohort.pdf"), height=6,width = 6)
pheatmap(as.matrix(t(regulon.gaby.sel)),
         cluster_rows = F,
         cluster_cols = as.hclust(cluster2.wt.gaby.mRNA.ans$dendro),
         border_color = "black",
         annotation_col = annCol.wt.gaby.genomic.alter[,c("RNA_Clust","Histology"),drop = F],
         annotation_colors = annColors.wt.gaby.genomic.alter[c("RNA_Clust","Histology")],
         color = NMF:::ccRamp(x = heatmap.L.BlYlRd,n = 64),
         treeheight_row = 15,
         treeheight_col = 15,
         cellwidth = 15,
         cellheight = 10,
         show_colnames = T,
         show_rownames = T,
         cutree_cols = 2,
         #gaps_row = c(14,22),
         fontsize_col = 8)
invisible(dev.off())

outTab <- NULL
for (i in comregulon) {
  tmp1 <- regulon.gaby.sel[RNA_Clust1,i]
  tmp2 <- regulon.gaby.sel[RNA_Clust2,i]
  
  wt <- wilcox.test(tmp1,tmp2,alternative = "less")
  outTab <- rbind.data.frame(outTab,
                             data.frame(reglon = i,
                                        avgTME1 = mean(tmp1),
                                        avgTME2 = mean(tmp2),
                                        p = wt$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
write.table(outTab,file.path(res.path,"wilcox test of common chromotin remodeling regulons in TME of gaby wilms cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)


p <- pheatmap(as.matrix(t(regulon.target.sel[com_sam,])),
              cluster_rows = F,
              cluster_cols = dendsort(hcs.target114),
              border_color = NA,
              annotation_col = annCol.target114[com_sam,c("tune_Clust","HISTOLOGY"),drop = F],
              annotation_colors = annColors.target114[c("tune_Clust","HISTOLOGY")],
              #color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
              color = NMF:::ccRamp(x = heatmap.L.BlYlRd,n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 5,
              cellheight = 13,
              #gaps_row = c(14,22),
              cutree_cols = 2,
              show_colnames = F,
              show_rownames = T,
              fontsize_col = 8)
pdf(file.path(fig.path, "common chromotin remodeling regulon wilms tumors in target cohort.pdf"), height=6,width = 12)
draw(p, annotation_legend_side = "bottom")
invisible(dev.off())

outTab <- NULL
for (i in comregulon) {
  tmp1 <- regulon.target.sel[RNA_Clust1.target,i]
  tmp2 <- regulon.target.sel[RNA_Clust2.target,i]
  
  wt <- wilcox.test(tmp1,tmp2,alternative = "less")
  outTab <- rbind.data.frame(outTab,
                             data.frame(reglon = i,
                                        avgTME1 = mean(tmp1),
                                        avgTME2 = mean(tmp2),
                                        p = wt$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
write.table(outTab,file.path(res.path,"wilcox test of common chromotin remodeling regulons in TME of target wilms cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

hcs <- hclust(distanceMatrix(as.matrix(t(regulon.ccle.sel)), "euclidean"), "average"); hcs.ccle.regulon <- hcs
group <- cutree(hcs,2)
annCol.ccle$chromatin_remodeling <- group
fisher.test(table(annCol.ccle$chromatin_remodeling,annCol.ccle$TME))

p <- pheatmap(as.matrix(t(regulon.ccle.sel)),
              cluster_rows = F,
              cluster_cols = hcs.ccle.tme,
              border_color = "black",
              annotation_col = annCol.ccle[,c("TME","Histology")],
              annotation_colors = annColors.ccle[c("TME","Histology")],
              color = NMF:::ccRamp(x = heatmap.L.BlYlRd,n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 15,
              cellheight = 10,
              show_colnames = T,
              show_rownames = T,
              cutree_cols = 2,
              fontsize_col = 9)
pdf(file.path(fig.path, "common chromotin remodeling regulon wilms tumors in ccle cohort.pdf"), height=6,width = 6)
draw(p, annotation_legend_side = "bottom", heatmap_legend_side = "left")
invisible(dev.off())

outTab <- NULL
for (i in comregulon) {
  tmp1 <- regulon.ccle.sel[RNA_Clust1.ccle,i]
  tmp2 <- regulon.ccle.sel[RNA_Clust2.ccle,i]
  
  wt <- wilcox.test(tmp1,tmp2,alternative = "less")
  outTab <- rbind.data.frame(outTab,
                             data.frame(reglon = i,
                                        avgTME1 = mean(tmp1),
                                        avgTME2 = mean(tmp2),
                                        p = wt$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
write.table(outTab,file.path(res.path,"wilcox test of common chromotin remodeling regulons in TME of ccle wilms cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

# box plot for regulons between TME clusters
dd <- as.data.frame(t(regulon.gaby.sel))
dd$gene = rownames(dd)
d2 <- gather(dd, sample, PCT, 1:(ncol(dd)-1))
tmp <- data.frame(sample=rownames(regulon.gaby.sel),TME=annCol.wt.gaby.genomic.alter[rownames(regulon.gaby.sel),"RNA_Clust"],stringsAsFactors = F)
tmp$TME <- ifelse(tmp$TME == "TP53_Wild","TME-C1","TME-C2")
d2 <- merge(d2,tmp,by="sample",all.x=T)
d2$gene <- factor(d2$gene,colnames(regulon.gaby.sel))
pvalues <- sapply(d2$gene, function(x) {
  res <- wilcox.test(PCT ~ TME, data = subset(d2, gene == x))$p.value
})
pv <- data.frame(gene = d2$gene, pvalue = pvalues)

pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

p1 <- ggplot(d2, aes(gene, PCT, fill=TME)) + 
  geom_boxplot(aes(color = TME),outlier.shape = NA) + scale_fill_manual(values = jco[2:1]) + scale_color_manual(values = jco[2:1]) +
  geom_text(aes(gene, y=max(d2$PCT) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("Regulon activity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10,colour = "black"),
        axis.text.y = element_text(colour = "black",size = 10),
        panel.grid = element_blank(),
        legend.position = "top",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dat <- ggplot_build(p1)$data[[1]]
p1 <- p1 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p1
ggsave(file.path(fig.path,"box plot for regulons in wilms tumor gaby cohort.pdf"),width = 8,height = 4)

dd <- as.data.frame(t(regulon.target.sel))
dd$gene = rownames(dd)
d2 <- gather(dd, sample, PCT, 1:(ncol(dd)-1))
tmp <- data.frame(sample=rownames(regulon.target.sel),TME=annCol.target114[rownames(regulon.target.sel),"tune_Clust"],stringsAsFactors = F)
tmp$TME <- ifelse(tmp$TME == "C1","tTME-C1","tTME-C2")
d2 <- merge(d2,tmp,by="sample",all.x=T)
d2$gene <- factor(d2$gene,colnames(regulon.gaby.sel))
pvalues <- sapply(d2$gene, function(x) {
  res <- wilcox.test(PCT ~ TME, data = subset(d2, gene == x))$p.value
})
pv <- data.frame(gene = d2$gene, pvalue = pvalues)

pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

p1 <- ggplot(d2, aes(gene, PCT, fill=TME)) + 
  geom_boxplot(aes(color = TME),outlier.shape = NA) + scale_fill_manual(values = jco[2:1]) + scale_color_manual(values = jco[2:1]) +
  geom_text(aes(gene, y=max(d2$PCT) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("Regulon activity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10,colour = "black"),
        axis.text.y = element_text(colour = "black",size = 10),
        panel.grid = element_blank(),
        legend.position = "top",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dat <- ggplot_build(p1)$data[[1]]
p1 <- p1 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p1
ggsave(file.path(fig.path,"box plot for regulons in wilms tumor target cohort.pdf"),width = 8,height = 4)

tmp <- as.data.frame(t(paths.gsva.ccle))
tmp$TME <- annCol.ccle[rownames(tmp),"TME"]
wilcox.test(tmp$HALLMARK_INTERFERON_GAMMA_RESPONSE~tmp$TME) # 0.004
wilcox.test(tmp$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION~tmp$TME) # 0.004
wilcox.test(tmp$HALLMARK_P53_PATHWAY~tmp$TME) # 0.004
wilcox.test(tmp$REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA~tmp$TME) # 0.004
wilcox.test(tmp$REACTOME_HDACS_DEACETYLATE_HISTONES~tmp$TME) # 0.07
wilcox.test(tmp$REACTOME_HOMOLOGY_DIRECTED_REPAIR~tmp$TME) # 0.048

onco.signature <- read.table(file.path(comAnn.path,"Oncogenetic_signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(onco.signature$Pathway)
onco.sig.ccr <- list()
for (i in cell.type) {
  onco.sig.ccr[[i]] <- intersect(toupper(onco.signature[which(onco.signature$Pathway == i),"Symbol"]),rownames(ccle.wt))
}
onco.score.ccle.wt <- gsva(as.matrix(ccle.wt),onco.sig.ccr,method="gsva",parallel.sz = 1)
dd <- as.data.frame(t(onco.score.ccle.wt))
dd$TME <- annCol.ccle[rownames(dd),"TME"]
dd$sample <- rownames(dd)
d2 <- gather(dd, cell, expr, 1:11)
d2$cell <- factor(d2$cell,levels = rownames(onco.score.ccle.wt))

pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(expr ~ TME, data = subset(d2, cell == x))
  return(res$p.value) #
})
pv <- data.frame(gene = d2$cell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

ggplot(d2, aes(cell, expr, fill=TME)) + scale_fill_manual(values = c(jco[2],jco[1])) +
  geom_boxplot() + 
  geom_text(aes(gene, y=max(d2$expr) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("Enrichment score") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(fig.path,"boxplot for oncogenetic pathways in two clusters in ccle wt cohort.pdf"),width = 8,height = 4)

tmp <- as.data.frame(t(onco.score.gaby.wt))
tmp$TME <- annCol.wt.gaby.genomic.alter[rownames(tmp),"RNA_Clust"]
tmp$HRS <- annCol.wt.gaby.genomic.alter[rownames(tmp),"HRS"]
tmp1 <- tmp[which(tmp$HRS == "RS-High"),"MYC"]
tmp2 <- tmp[which(tmp$HRS == "RS-Low"),"MYC"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(4,6)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for MYC in gaby WT cohort regarding HRS subtype.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "MYC oncogenic pathway\nGSVA score",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(annColors.target114$HRS[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(annColors.target114$HRS[2:1]))
text(1.5,-0.4,"P = 0.038",cex = 1)
invisible(dev.off())

tmp <- as.data.frame(t(onco.score.gaby.wt))
tmp$TME <- annCol.wt.gaby.genomic.alter[rownames(tmp),"RNA_Clust"]
tmp$HRS <- annCol.wt.gaby.genomic.alter[rownames(tmp),"HRS"]
tmp1 <- tmp[which(tmp$TME == "TP53_Wild"),"MYC"]
tmp2 <- tmp[which(tmp$TME == "TP53_Mutated"),"MYC"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for MYC in gaby WT cohort regarding TME subtype.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "MYC oncogenic pathway\nGSVA score",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2:1]))
text(1.5,-0.35,"P = 0.095",cex = 1)
invisible(dev.off())

tmp <- as.data.frame(t(onco.score.target.wt))
tmp$TME <- annCol.target114[rownames(tmp),"tune_Clust"]
tmp$HRS <- annCol.target114[rownames(tmp),"HRS"]
tmp1 <- tmp[which(tmp$HRS == "RS-High"),"MYC"]
tmp2 <- tmp[which(tmp$HRS == "RS-Low"),"MYC"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(65,31)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for MYC in target WT cohort regarding HRS subtype.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "MYC oncogenic pathway\nGSVA score",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(annColors.target114$HRS[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(annColors.target114$HRS[2:1]))
text(1.5,-0.4,"P = 0.178",cex = 1)
invisible(dev.off())

tmp <- as.data.frame(t(onco.score.target.wt))
tmp$TME <- annCol.target114[rownames(tmp),"tune_Clust"]
tmp$HRS <- annCol.target114[rownames(tmp),"HRS"]
tmp1 <- tmp[which(tmp$TME == "C1"),"MYC"]
tmp2 <- tmp[which(tmp$TME == "C2"),"MYC"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("tTME-C1","tTME-C2"),c(69,27)))
tmp$mut <- factor(tmp$mut,levels = c("tTME-C1","tTME-C2"))
pdf(file.path(fig.path,"boxplot for MYC in target WT cohort regarding TME subtype.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "MYC oncogenic pathway\nGSVA score",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2:1]))
text(1.5,-0.35,"P = 0.038",cex = 1)
invisible(dev.off())

#-------------------------#
# merge immune annotation #
immune.bk <- read.delim(file.path(comAnn.path,"Immune.gene.background.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
deres <- read.table(file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.TP53_Mutated_vs_TP53_Wild.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
colnames(deres)[2] <- "HUGO Name"
deres <- merge(deres, immune.bk, by = "HUGO Name", all.x = TRUE)
write.table(deres,file.path(res.path,"Gaby_WT_mRNA_edgeR_test_result.TP53_Mutated_vs_TP53_Wild_withImmuneAnnotation.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

deres <- read.table(file.path(res.path,"TARGET_WT_mRNA_edgeR_test_result.C2_vs_Others.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
colnames(deres)[2] <- "HUGO Name"
deres <- merge(deres, immune.bk, by = "HUGO Name", all.x = TRUE)
write.table(deres,file.path(res.path,"TARGET_WT_mRNA_edgeR_test_result.C2_vs_Others_withImmuneAnnotation.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

deres <- read.table(file.path(res.path,"TARGET_DAWT_mRNA_edgeR_test_result.C2_vs_Others.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
colnames(deres)[2] <- "HUGO Name"
deres <- merge(deres, immune.bk, by = "HUGO Name", all.x = TRUE)
write.table(deres,file.path(res.path,"TARGET_DAWT_mRNA_edgeR_test_result.C2_vs_Others_withImmuneAnnotation.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

#-----------------------------------------------------#
# REDO ANALYSIS FOR USING DAWT IN TARGETE COHORT ONLY #
#-----------------------------------------------------#

#-----------------------------------------------------------#
# separate km curve in target cohort according to histology #
tmp <- data.frame(OS.time = annCol.target114$OS.time/365,
                  OS = ifelse(annCol.target114$OS == "Dead",1,0),
                  histology = annCol.target114$HISTOLOGY,
                  TP53.mut = annCol.target114$TP53_MUT,
                  TP53.cna = annCol.target114$TP53_LOSS,
                  TP53.alter = annCol.target114$TP53_ALTER,
                  tune_Clust = annCol.target114$tune_Clust,
                  Age = annCol.target114$Age/365,
                  Stage = annCol.target114$Stage,
                  stringsAsFactors = F)
tmp$Stage <- factor(tmp$Stage,levels = c("I_II",">III"))
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1],seagreen),
                pval = T,
                data=tmp,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Years)",
                ylab = "Overall Survival (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"KM for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort 2 clusters.pdf"),width = 5,height = 5)
p
invisible(dev.off())

# DAWT
tmp1 <- tmp[which(tmp$histology == "DAWT"),]
table(tmp1$tune_Clust, tmp1$TP53.alter)
fisher.test(table(tmp1$tune_Clust, tmp1$TP53.alter))
table(tmp1$tune_Clust, tmp1$TP53.mut)
fisher.test(table(tmp1$tune_Clust, tmp1$TP53.mut))
table(tmp1$tune_Clust, tmp1$TP53.cna)
fisher.test(table(tmp1$tune_Clust, tmp1$TP53.cna))

fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp1, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp1, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1],seagreen),
                pval = T,
                data=tmp1,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Years)",
                ylab = "Overall Survival (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"KM of DAWT for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort 2 clusters.pdf"),width = 5,height = 5)
p
invisible(dev.off())

# FHWT
tmp2 <- tmp[which(tmp$histology == "FHWT"),]
table(tmp2$tune_Clust, tmp2$TP53.alter)
fisher.test(table(tmp2$tune_Clust, tmp2$TP53.alter))
table(tmp2$tune_Clust, tmp2$TP53.mut)
fisher.test(table(tmp2$tune_Clust, tmp2$TP53.mut))
table(tmp2$tune_Clust, tmp2$TP53.cna)
fisher.test(table(tmp2$tune_Clust, tmp2$TP53.cna))

fitd <- survdiff(Surv(OS.time, OS) ~ tune_Clust, data=tmp2, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ tune_Clust, data=tmp2, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1],seagreen),
                pval = T,
                data=tmp2,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Years)",
                ylab = "Overall Survival (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"KM of FHWT for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort 2 clusters.pdf"),width = 5,height = 5)
p
invisible(dev.off())

dawt.target <- rownames(annCol.target114[which(annCol.target114$HISTOLOGY == "DAWT"),])
fhwt.target <- rownames(annCol.target114[which(annCol.target114$HISTOLOGY == "FHWT"),])

com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- target.pure$es.pure[,intersect(com_sam,dawt.target)]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D"); hcs.target114.dawt <- hcs
group <- cutree(hcs,2)
group <- ifelse(group == 1,2,1)
annCol.target.dawt <- annCol.target114[intersect(com_sam, dawt.target),]
annCol.target.fhwt <- annCol.target114[intersect(com_sam, fhwt.target),]

annCol.target.dawt$Agem <- ifelse(annCol.target.dawt$Age > 1771, ">4.9","<=4.9")
annCol.target.dawt$TME <- paste0("C",group)
annColors.target.dawt <- annColors.target114
annColors.target.dawt$TME <- c("C1" = jco[2],"C2" = jco[1])
annColors.target.dawt$Agem <- c(">4.9" = "#fa7921","<=4.9" = "#FA792199")
p <- pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
              cluster_rows = F,
              cluster_cols = dendsort(hcs.target114.dawt),
              border_color = NA,
              annotation_col = annCol.target.dawt[colnames(plotdata),c("TME","Stage","Gender","Agem","OS","TP53_LOSS","TP53_MUT")],
              annotation_colors = annColors.target.dawt[c("TME","Stage","Gender","Agem","OS","TP53_LOSS","TP53_MUT")],
              color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 10,
              cellheight = 13,
              gaps_row = c(14,22),
              cutree_cols = 2,
              show_colnames = F,
              show_rownames = T,
              fontsize_col = 8)
pdf(file.path(fig.path, "CCR immune signature heatmap of 114 wilms tumors in target cohort DAWT (purified by tp) 2 clusters.pdf"), height=10,width = 15)
draw(p, annotation_legend_side = "bottom")
invisible(dev.off())

tmp <- data.frame(OS.time = annCol.target.dawt$OS.time/365,
                  OS = ifelse(annCol.target.dawt$OS == "Dead",1,0),
                  TME = annCol.target.dawt$TME,
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ TME, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ TME, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1]),
                pval = T,
                data=tmp,
                size=1,
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Years)",
                ylab = "Overall Survival (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"KM of new cluster of DAWT for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort 2 clusters.pdf"),width = 3.5,height = 5)
p
invisible(dev.off())

TLS.target.dawt <- apply(target.wt[TLS,colnames(plotdata)], 2, geometric.mean.rm0)
annCol.target.dawt$TLS <- annTrackScale(TLS.target.dawt, halfwidth = 2)
indata <- target.wt[TLS,colnames(plotdata)]
plotdata <- standarize.fun(indata,halfwidth = 2)
hm1 <- pheatmap(as.matrix(plotdata),
              cluster_rows = F,
              cluster_cols = dendsort(hcs.target114.dawt),
              border_color = NA,
              annotation_col = annCol.target.dawt[colnames(plotdata),c("TME","TLS")],
              annotation_colors = annColors.target.dawt[c("TME","TLS")],
              color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 10,
              cellheight = 13,
              show_colnames = F,
              show_rownames = T,
              cutree_cols = 2,
              fontsize_col = 8)
indata <- target.wt[c("TMEM173","MB21D1"),colnames(plotdata)]
plotdata <- standarize.fun(indata,halfwidth = 2)
hm2 <- pheatmap(as.matrix(plotdata),
                cluster_rows = F,
                cluster_cols = dendsort(hcs.target114.dawt),
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                treeheight_row = 15,
                treeheight_col = 15,
                cellwidth = 10,
                cellheight = 13,
                show_colnames = F,
                show_rownames = T,
                cutree_cols = 2,
                fontsize_col = 8)
plotdata <- standarize.fun(target.cgassting,halfwidth = 2)
hm3 <- pheatmap(as.matrix(plotdata),
                cluster_rows = F,
                cluster_cols = dendsort(hcs.target114.dawt),
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                treeheight_row = 15,
                treeheight_col = 15,
                cellwidth = 10,
                cellheight = 13,
                show_colnames = F,
                show_rownames = T,
                cutree_cols = 2,
                fontsize_col = 8)
pdf(file.path(fig.path, "TLS heatmap of wilms dawt in target cohort.pdf"), height=6,width = 15)
draw(hm1 %v% hm2 %v% hm3, annotation_legend_side = "bottom")
invisible(dev.off())

# C2 vs Others
tmp <- annCol.target.dawt$TME; names(tmp) <- rownames(annCol.target.dawt)
tmp <- na.omit(tmp)
tmp <- ifelse(tmp == "C2","C2","Others")
Groupinfo <- annCol.target.dawt[names(tmp),]
colnames(Groupinfo)[51] <- "group"
Groupinfo$group <- ifelse(Groupinfo$group == "C2","C2","Others")
PASSFlag <- rep(TRUE,length(Mids)); names(PASSFlag) <- Mids
complist <- createList.RNAclust.targetC2(group=tmp)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[Mids, names(tmp)], tailrows, Groupinfo=Groupinfo, features=Mids, featType="TARGET_DAWT_mRNA", complist,PASSFlag=PASSFlag, overwt=TRUE)

deres <- read.table(file.path(res.path,"TARGET_DAWT_mRNA_edgeR_test_result.C2_vs_Others.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- deres$logFC
names(geneList) <- deres$id
geneList <- sort(geneList,decreasing = T)
gsea_target_dawt_c2vsc1 <- GSEA(geneList = geneList,
                                TERM2GENE=MSigDB.hallmark,
                                pvalueCutoff = 0.25,
                                nPerm = 10000,
                                seed = T,
                                verbose=F)
res <- data.frame(gsea_target_dawt_c2vsc1)
write.table(as.data.frame(gsea_target_dawt_c2vsc1),file.path(res.path,"GSEA hallmark results for TARGET DAWT mRNA C2 vs C1 DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

gsea_target_dawt_c2vsc1_hm.rt <- GSEA(geneList = geneList,
                                      TERM2GENE=MSigDB.hm.rt,
                                      pvalueCutoff = 2,
                                      nPerm = 10000,
                                      seed = T,
                                      verbose=F)
res <- data.frame(gsea_target_dawt_c2vsc1_hm.rt)
write.table(as.data.frame(gsea_target_dawt_c2vsc1_hm.rt),file.path(res.path,"GSEA hallmark and reactome results for TARGET DAWT mRNA C2 vs C1 DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

gseaplot2(x = gsea_target_dawt_c2vsc1_hm.rt,
          geneSetID = c("HALLMARK_P53_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
          pvalue_table = F,
          color = c(alpha(darkblue,0.3),alpha(darkblue,0.6),darkblue))
ggsave(filename = file.path(fig.path,"GSEA downregulated pathway in target DAWT cohort.pdf"), width = 5.5,height = 4)

gseaplot2(x = gsea_target_dawt_c2vsc1_hm.rt,
          geneSetID = c("REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_CHROMATIN_MODIFYING_ENZYMES","REACTOME_HDACS_DEACETYLATE_HISTONES"),
          pvalue_table = F,
          color = c(alpha(darkred,0.3),alpha(darkred,0.6),darkred))
ggsave(filename = file.path(fig.path,"GSEA upregulated pathway in target DAWT cohort.pdf"), width = 5.5,height = 4)

gseaplot2(x = gsea_target_dawt_c2vsc1,
          geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          pvalue_table = F,
          color = seagreen)
ggsave(filename = file.path(fig.path,"GSEA inteferon gamma in target DAWT.pdf"), width = 5.5,height = 4)

gsea_target_dawt_c2vsc1_all <- GSEA(geneList = geneList,
                                    TERM2GENE=MSigDB,
                                    pvalueCutoff = 0.25,
                                    nPerm = 10000,
                                    seed = T,
                                    verbose=F)
res <- data.frame(gsea_target_dawt_c2vsc1_all)
write.table(as.data.frame(gsea_target_dawt_c2vsc1_all),file.path(res.path,"GSEA all results for TARGET DAWT mRNA TP53 Mutated C2 vs Wild C1 DEGs.txt"),row.names = T,col.names = NA,sep = "\t",quote = F)

RNA_Clust1.target.dawt <- rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C1"),])
RNA_Clust2.target.dawt <- rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C2"),])

ezh2.target.dawt.c1 <- as.numeric(target.wt["EZH2",RNA_Clust1.target.dawt])
ezh2.target.dawt.c2 <- as.numeric(target.wt["EZH2",RNA_Clust2.target.dawt])
wilcox.test(ezh2.target.dawt.c1,ezh2.target.dawt.c2) #0.0156
tmp <- data.frame(exp = c(ezh2.target.dawt.c1,ezh2.target.dawt.c2),
                  mut = rep(c("C1","C2"),c(24,12)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for EZH2 expression in target DAWT cohort regarding RNA cluster.pdf"),width = 2.5,height = 2.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "RNA Cluster",
        ylab = "log2(TPM) in target",
        #ylim = c(1,5),
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,7.5,"P = 0.016",cex = 1)
invisible(dev.off())

# EZH2 pathway
tmp <- MSigDB[grep("KAMMINGA_EZH2_TARGETS",MSigDB$ont),]
ezh2.path.target.dawt <- gsva(expr = as.matrix(target.wt[,rownames(annCol.target.dawt)]),
                              gset.idx.list = list("KAMMINGA_EZH2_TARGETS" = tmp$gene),
                              method = "ssgsea")
ezh2.path.target.dawt.c1 <- as.numeric(ezh2.path.target.dawt[1,RNA_Clust1.target.dawt])
ezh2.path.target.dawt.c2 <- as.numeric(ezh2.path.target.dawt[1,RNA_Clust2.target.dawt])
wilcox.test(ezh2.path.target.dawt.c1,ezh2.path.target.dawt.c2, alternative = "less")
tmp <- data.frame(exp = c(ezh2.path.target.dawt.c1,ezh2.path.target.dawt.c2),
                  mut = rep(c("C1","C2"),c(24,12)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for EZH2 target in target DAWT cohort regarding RNA cluster.pdf"),width = 2.5,height = 2.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "RNA Cluster",
        ylab = "ssGSEA score",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,2.2,"P = 0.039",cex = 1)
invisible(dev.off())

tmp1 <- annCol.cbio[RNA_Clust1.target.dawt,"Fraction Genome Altered"]
tmp2 <- annCol.cbio[RNA_Clust2.target.dawt,"Fraction Genome Altered"]
wilcox.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("Desert","Infiltrated"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("Infiltrated","Desert"))
pdf(file.path(fig.path,"boxplot for Fraction Genome Altered in cbioportal 101 DAWT regarding CCR immune signature group.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Fraction Genome Altered",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.6,"P = 0.458",cex = 1)
invisible(dev.off())

indata <- target.pure$es.pure[,rownames(annCol.target.dawt)]
tmp <- as.data.frame(t(indata))
tmp$TME <- annCol.target.dawt$TME
d2 <- gather(tmp, gene, expr, `B.cells.naive`:`Fibroblasts`, factor_key=TRUE)
pvalues <- sapply(d2$gene, function(x) {
  res <- wilcox.test(expr ~ TME, data = subset(d2, gene == x))$p.value
})
pv <- data.frame(gene = d2$gene, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
wilcox.test(target.pure$es.pure["T.cells.CD4.memory.activated",RNA_Clust1.target.dawt],
            target.pure$es.pure["T.cells.CD4.memory.activated",RNA_Clust2.target.dawt])

p <- ggplot(d2, aes(gene, expr, fill=TME)) + 
  geom_boxplot(aes(col = TME),outlier.shape = NA,alpha = 1) + 
  geom_text(aes(gene, y=max(expr)), 
            label=pv$sigcode,
            data=d2, 
            inherit.aes=F) + 
  scale_fill_manual(values = jco[2:1]) + 
  scale_color_manual(values = jco[2:1]) + 
  xlab(NULL) + ylab("Tumor microenvironment\n[GSVA]") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90,size = 12, color = "black"),
        axis.title.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
dat <- ggplot_build(p)$data[[1]]
p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p
ggsave(filename = file.path(fig.path,"boxplot for tme in target DAWT.pdf"), width = 6,height = 5)

indata <- target.pure$es.pure[,c(RNA_Clust1.target,RNA_Clust2.target)]
tmp <- as.data.frame(t(indata))
tmp$TME <- annCol.target114[c(RNA_Clust1.target,RNA_Clust2.target),"tune_Clust"]
d2 <- gather(tmp, gene, expr, `B.cells.naive`:`Fibroblasts`, factor_key=TRUE)
pvalues <- sapply(d2$gene, function(x) {
  res <- wilcox.test(expr ~ TME, data = subset(d2, gene == x))$p.value
})
pv <- data.frame(gene = d2$gene, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

p <- ggplot(d2, aes(gene, expr, fill=TME)) + 
  geom_boxplot(aes(col = TME),outlier.shape = NA,alpha = 1) + 
  geom_text(aes(gene, y=max(expr)), 
            label=pv$sigcode,
            data=d2, 
            inherit.aes=F) + 
  scale_fill_manual(values = jco[2:1]) + 
  scale_color_manual(values = jco[2:1]) + 
  xlab(NULL) + ylab("Tumor microenvironment\n[GSVA]") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90,size = 12, color = "black"),
        axis.title.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
dat <- ggplot_build(p)$data[[1]]
p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p
ggsave(filename = file.path(fig.path,"boxplot for tme in target cohort.pdf"), width = 6,height = 5)

indata <- imm.score.gaby.wt3.pure[,rownames(annCol.wt.gaby)]
tmp <- as.data.frame(t(indata))
tmp$TME <- ifelse(annCol.wt.gaby$RNA_Clust == "TP53_Mutated","C2","C1")
d2 <- gather(tmp, gene, expr, `B.cells.naive`:`Fibroblasts`, factor_key=TRUE)
pvalues <- sapply(d2$gene, function(x) {
  res <- wilcox.test(expr ~ TME, data = subset(d2, gene == x))$p.value
})
pv <- data.frame(gene = d2$gene, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

p <- ggplot(d2, aes(gene, expr, fill=TME)) + 
  geom_boxplot(aes(col = TME),outlier.shape = NA,alpha = 1) + 
  geom_text(aes(gene, y=max(expr)), 
            label=pv$sigcode,
            data=d2, 
            inherit.aes=F) + 
  scale_fill_manual(values = jco[2:1]) + 
  scale_color_manual(values = jco[2:1]) + 
  xlab(NULL) + ylab("Tumor microenvironment\n[GSVA]") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90,size = 12, color = "black"),
        axis.title.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
dat <- ggplot_build(p)$data[[1]]
p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p
ggsave(filename = file.path(fig.path,"boxplot for tme in gabt wt.pdf"), width = 6,height = 5)

dd <- as.data.frame(t(regulon.target.sel[rownames(annCol.target.dawt),]))
dd$gene = rownames(dd)
d2 <- gather(dd, sample, PCT, 1:(ncol(dd)-1))
tmp <- data.frame(sample=rownames(regulon.target.sel[rownames(annCol.target.dawt),]),TME=annCol.target.dawt$TME,stringsAsFactors = F)
tmp$TME <- ifelse(tmp$TME == "C1","tTME-C1","tTME-C2")
d2 <- merge(d2,tmp,by="sample",all.x=T)
d2$gene <- factor(d2$gene,colnames(regulon.gaby.sel))
pvalues <- sapply(d2$gene, function(x) {
  res <- wilcox.test(PCT ~ TME, data = subset(d2, gene == x))$p.value
})
pv <- data.frame(gene = d2$gene, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))

p1 <- ggplot(d2, aes(gene, PCT, fill=TME)) + 
  geom_boxplot(aes(color = TME),outlier.shape = NA) + scale_fill_manual(values = jco[2:1]) + scale_color_manual(values = jco[2:1]) +
  geom_text(aes(gene, y=max(d2$PCT) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("Regulon activity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10,colour = "black"),
        axis.text.y = element_text(colour = "black",size = 10),
        panel.grid = element_blank(),
        legend.position = "top",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dat <- ggplot_build(p1)$data[[1]]
p1 <- p1 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p1
ggsave(file.path(fig.path,"box plot for regulons in wilms DAWT target.pdf"),width = 8,height = 4)

#-----------------------#
# same analysis in FHWT #
com_sam <- intersect(target.pure$keep.sam,rownames(annCol.target114))
indata <- target.pure$es.pure[,intersect(com_sam,fhwt.target)]
plotdata <- standarize.fun(indata,halfwidth = 2)
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D"); hcs.target114.fhwt <- hcs
group <- cutree(hcs,2)
group <- ifelse(group == 1,2,1)
annCol.target.fhwt <- annCol.target114[intersect(com_sam, fhwt.target),]
annCol.target.fhwt$Agem <- ifelse(annCol.target.fhwt$Age > 1364.5, ">3.7","<=3.7")
annCol.target.fhwt$TME <- paste0("C",group)
annColors.target.fhwt <- annColors.target114
annColors.target.fhwt$TME <- c("C1" = jco[2],"C2" = jco[1])
annColors.target.fhwt$Agem <- c(">3.7" = "#fa7921","<=3.7" = "#FA792199")
p <- pheatmap(as.matrix(plotdata[immune.sig.ccr.order,]),
              cluster_rows = F,
              cluster_cols = dendsort(hcs.target114.fhwt),
              border_color = NA,
              annotation_col = annCol.target.fhwt[colnames(plotdata),c("TME","Stage","Gender","Agem","OS","TP53_LOSS","TP53_MUT")],
              annotation_colors = annColors.target.fhwt[c("TME","Stage","Gender","Agem","OS","TP53_LOSS","TP53_MUT")],
              color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 6,
              cellheight = 13,
              gaps_row = c(14,22),
              cutree_cols = 2,
              show_colnames = F,
              show_rownames = T,
              fontsize_col = 10)
pdf(file.path(fig.path, "CCR immune signature heatmap of 114 wilms tumors in target cohort fhwt (purified by tp) 2 clusters.pdf"), height=10,width = 15)
draw(p, annotation_legend_side = "bottom")
invisible(dev.off())

tmp <- data.frame(OS.time = annCol.target.fhwt$OS.time/365,
                  OS = ifelse(annCol.target.fhwt$OS == "Dead",1,0),
                  TME = annCol.target.fhwt$TME,
                  stringsAsFactors = F)
tmp <- as.data.frame(na.omit(tmp))
fitd <- survdiff(Surv(OS.time, OS) ~ TME, data=tmp, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ TME, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = c(jco[2],jco[1]),
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Years)",
                ylab = "Overall Survival (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"KM of new cluster of FHWT for group of CCR immune signature by target 114 (purified by tp) data in TARGET cohort 2 clusters.pdf"),width = 4,height = 4)
p
invisible(dev.off())

#------------------------------#
# ERV (use supervised in DAWT) #
indata <- log2(erv.normcount[,rownames(annCol.target.dawt)] + 1)
hcs <- hclust(distanceMatrix(as.matrix(indata), "pearson"), "ward.D")

erv.group.target.dawt <- cutree(hcs, 2)
erv.group.target.dawt <- ifelse(erv.group.target.dawt == 1,2,1)
erv.group.target.dawt <- paste0("ERV-C",erv.group.target.dawt); names(erv.group.target.dawt) <- colnames(indata)
annCol.target.dawt$ERV <- erv.group.target.dawt[rownames(annCol.target.dawt)]
annColors.target.dawt[["ERV"]] <- c("ERV-C1" = orange,"ERV-C2" = green)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = erv.count[,names(erv.group.target.dawt)],
                                      colData = data.frame(samID = names(erv.group.target.dawt),
                                                           erv = as.character(erv.group.target.dawt),
                                                           row.names = names(erv.group.target.dawt)),
                                      design = as.formula("~ erv"))

dds$erv <- relevel(dds$erv,ref = "ERV-C1")
dds <- DESeq(dds)
res <- DESeq2::results(dds, contrast = c("erv","ERV-C1","ERV-C2"))
resData <- as.data.frame(res[order(res$padj),])
resData$id <- rownames(resData)
resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
colnames(resData) <- c("id","baseMean","log2fc","lfcSE","stat","pvalue","padj")
resData$fc <- 2^resData$log2fc
resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
resData <- cbind.data.frame(resData,erv.anno[rownames(resData),])
write.table(resData,file.path(res.path,"deseq2 results between target DAWT ERV 2 groups using all features.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

resData <- read.table(file.path(res.path,"deseq2 results between target DAWT ERV 2 groups using all features.txt"),sep = "\t",row.names = 1,header = T,check.names = F)
plotdata <- standarize.fun(log2(erv.normcount[,rownames(annCol.target.dawt)] + 1)[rownames(resData[which(resData$padj < 0.05 & abs(resData$log2fc) > log2(1.5)),]),], halfwidth = 2)
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")
p <- pheatmap(as.matrix(plotdata),
              cluster_rows = hcg,
              cluster_cols = hcs,
              border_color = NA,
              annotation_col =  annCol.target.dawt[colnames(plotdata),c("ERV","TME","OS","TP53_MUT", "TP53_LOSS")],
              annotation_colors = annColors.target.dawt[c("ERV","TME","OS","TP53_MUT", "TP53_LOSS")],
              color = NMF:::ccRamp(x = heatmap.GrWtRd,n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 6,
              cellheight = 0.8,
              show_colnames = F,
              show_rownames = F,
              cutree_cols = 2,
              fontsize_col = 8)
pdf(file.path(fig.path, "unsupervised clustering of ERVs in target wilms DAWT.pdf"), height=12,width = 10)
draw(p, annotation_legend_side = "left")
invisible(dev.off())

resData <- read.table(file.path(res.path,"deseq2 results between gaby ERV 2 groups using all features.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp1 <- resData[which(resData$padj < 0.05 & abs(resData$log2fc) > log2(1.5)),]

indata <- log2(erv.normcount[rownames(tmp1),rownames(annCol.target.dawt)] + 1)
hcs <- hclust(distanceMatrix(as.matrix(indata), "pearson"), "ward.D")
hcs <- click_rotate(as.dendrogram(hcs))
hcs <- as.hclust(hcs)
erv.group.target.dawt2 <- cutree(hcs, 4)
table(erv.group.target.dawt2)
erv.group.target.dawt2 <- paste0("ERV-C",erv.group.target.dawt2); names(erv.group.target.dawt2) <- colnames(indata)
erv.group.target.dawt2 <- ifelse(erv.group.target.dawt2 == "ERV-C2","ERV-C1","ERV-C2")
annCol.target.dawt$ERV2 <- erv.group.target.dawt2[rownames(annCol.target.dawt)]
annColors.target.dawt[["ERV2"]] <- c("ERV-C1" = orange,"ERV-C2" = green)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = erv.count[,names(erv.group.target.dawt2)],
                                      colData = data.frame(samID = names(erv.group.target.dawt2),
                                                           erv = as.character(erv.group.target.dawt2),
                                                           row.names = names(erv.group.target.dawt2)),
                                      design = as.formula("~ erv"))

dds$erv <- relevel(dds$erv,ref = "ERV-C1")
dds <- DESeq(dds)
res <- DESeq2::results(dds, contrast = c("erv","ERV-C1","ERV-C2"))
resData <- as.data.frame(res[order(res$padj),])
resData$id <- rownames(resData)
resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
colnames(resData) <- c("id","baseMean","log2fc","lfcSE","stat","pvalue","padj")
resData$fc <- 2^resData$log2fc
resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
resData <- cbind.data.frame(resData,erv.anno[rownames(resData),])
write.table(resData,file.path(res.path,"deseq2 results between supervised target DAWT ERV 2 groups using all features.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
resData <- read.table(file.path(res.path,"deseq2 results between supervised target DAWT ERV 2 groups using all features.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
plotdata <- standarize.fun(log2(erv.normcount[,rownames(annCol.target.dawt)] + 1)[rownames(resData[which(resData$padj < 0.05 & abs(resData$log2fc) > log2(2)),]),], halfwidth = 2)

hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "euclidean"), "ward.D")
p <- pheatmap(as.matrix(plotdata),
              cluster_rows = hcg,
              cluster_cols = hcs,
              border_color = NA,
              annotation_col =  annCol.target.dawt[colnames(plotdata),c("ERV2","TME","HISTOLOGY","OS","TP53_MUT", "TP53_LOSS")],
              annotation_colors = annColors.target.dawt[c("ERV2","TME","HISTOLOGY","OS","TP53_MUT", "TP53_LOSS")],
              color = NMF:::ccRamp(x = heatmap.GrWtRd,n = 64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 3.5,
              cellheight = 2.5,
              show_colnames = F,
              show_rownames = F,
              fontsize_col = 8)
pdf(file.path(fig.path, "supervised clustering of ERVs in target wilms DAWT.pdf"), height=12,width = 10)
draw(p, annotation_legend_side = "left")
invisible(dev.off())
fisher.test(table(annCol.target.dawt$TME,annCol.target.dawt$ERV2)) # p = 0.00037

resData <- read.table(file.path(res.path,"deseq2 results between supervised target DAWT ERV 2 groups using all features.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dim(resData[which(resData$padj < 0.05 & resData$log2fc > 1.5),])
dim(resData[which(resData$padj < 0.05 & resData$log2fc < -1.5),])

plot_mode <- "advanced"
logFCcut <- log2(1.5) 
adjPcut <- 0.05 
pvalCut <- 0.05

x <- resData
xmin <- (range(x$log2fc)[1]- (range(x$log2fc)[1]+ 10))
xmax <- (range(x$log2fc)[1]+ (10-range(x$log2fc)[1]))
ymin <- 0
ymax <- 5

n1 <- length(x[, 1])
size <- rep(2, n1)
cols <- rep("grey30", n1)
names(cols)<- rownames(x)

cols[x$padj < adjPcut & x$log2fc >logFCcut]<- "#FB9A99"
cols[x$padj < adjPcut & x$log2fc > 1]<- "#ED4F4F"
cols[x$padj < adjPcut & x$log2fc < -logFCcut]<- "#B2DF8A"
cols[x$padj < adjPcut & x$log2fc < -1]<- "#329E3F"
color_transparent <- adjustcolor(cols, alpha.f = 0.8)
x$color_transparent <- color_transparent
size[x$padj < adjPcut & x$log2fc > logFCcut]<- 4
size[x$padj < adjPcut & x$log2fc > 1]<- 6
size[x$padj < adjPcut & x$log2fc < -logFCcut]<- 4
size[x$padj < adjPcut & x$log2fc < -1]<- 6

# Construct the plot object
p1 <- ggplot(data=x, aes(log2fc, -log10(padj))) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  
  labs(x="log2FoldChange", y="-log10FDR", title="") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(-3,-2, -1, 0, 1, 2,3), 
    labels = c(-3,-2, -1, 0, 1, 2,3),
    limits = c(-3, 3) 
  ) +
  geom_vline(xintercept = c(-1,-logFCcut, logFCcut,1), color="grey40", 
             linetype="longdash", lwd = 0.5) + 
  geom_hline(yintercept = -log10(pvalCut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12
  ) +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size = 12, color = "black"))

p1
ggsave(file.path(fig.path,"volcano plot of differential erv in target DAWT.pdf"), width = 4,height = 4)

#----------------------------------#
# high replication stress analysis #
library(dendextend)
HRS.score.target.dawt <- gsva(as.matrix(target.wt[,rownames(annCol.target.dawt)]),
                              HRS.signature,
                              method = "ssgsea")

plotdata <- standarize.fun(HRS.score.target.dawt,halfwidth = 1)
hcs <- hclust(distanceMatrix(as.matrix(plotdata), "euclidean"), "ward.D")
hcs <- click_rotate(as.dendrogram(hcs))
hcs <- as.hclust(hcs)
HRS.group.target.dawt <- cutree(hcs, k =2)
annCol.target.dawt[colnames(HRS.score.target.dawt),"HRS"] <- ifelse(HRS.group.target.dawt[colnames(HRS.score.target.dawt)] == 1,"RS-High","RS-Low")
annColors.target.dawt[["HRS"]] <- c("RS-Low" = "#4558A3","RS-High" = "#C86C0D")

p <- pheatmap(as.matrix(plotdata),
              cluster_rows = F,
              cluster_cols = hcs,
              border_color = NA,
              annotation_col = annCol.target.dawt[colnames(plotdata),c("HRS","TME", "HISTOLOGY", "OS","TP53_MUT","TP53_LOSS")],
              annotation_colors = annColors.target.dawt[c("HRS","TME", "HISTOLOGY", "OS","TP53_MUT","TP53_LOSS")],
              color = viridisLite::viridis(64),
              treeheight_row = 15,
              treeheight_col = 15,
              cellwidth = 5,
              cellheight = 13,
              cutree_cols = 2,
              show_colnames = F,
              show_rownames = T,
              fontsize_col = 8)
pdf(file.path(fig.path, "HRS signature heatmap of wilms tumors in target DAWT cohort.pdf"), height=10,width = 20)
draw(p, annotation_legend_side = "bottom", heatmap_legend_side = "left")
invisible(dev.off())

fisher.test(table(annCol.target.dawt$HRS,annCol.target.dawt$TP53_MUT)) # 0.6457
fisher.test(table(annCol.target.dawt$HRS,annCol.target.dawt$TME)) # 0.1464

#---------------------------#
# drug sensitivity analysis #
#---------------------------#

testExpr <- target.wt[,rownames(annCol.target.dawt)]
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
for (i in 1:ncol(trainPtype)) {
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- as.vector(trainPtype[,d])
  
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = T,
                                  selection = 1))
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
gdsc.pred.auc.target.dawt <- as.data.frame(t(outTab))
gdsc.pred.auc.target.dawt$HRS <- annCol.target.dawt[rownames(gdsc.pred.auc.target.dawt),"HRS"]
gdsc.pred.auc.target.dawt$TME <- annCol.target.dawt[rownames(gdsc.pred.auc.target.dawt),"TME"]
gdsc.pred.auc.target.dawt$TP53 <- annCol.target.dawt[rownames(gdsc.pred.auc.target.dawt),"TP53_MUT"]
write.table(gdsc.pred.auc.target.dawt,file = file.path(res.path,"predicted ln_ic50 using GDSC dataset for target DAWT cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

tmp1 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$HRS == "RS-High"),2]
tmp2 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$HRS == "RS-Low"),2]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(31,5)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD6739 ATR inhibitor in target DAWT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "AZD6738\nEstimated ln(IC50)",
        main = "ATR inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2.0,"P < 0.001",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$HRS == "RS-High"),5]
tmp2 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$HRS == "RS-Low"),5]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("RS-High","RS-Low"),c(31,5)))
tmp$mut <- factor(tmp$mut,levels = c("RS-High","RS-Low"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD1775 WEE1 inhibitor in target DAWT cohort regarding replication stress group.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Replication stress",
        ylab = "AZD1775\nEstimated ln(IC50)",
        main = "WEE1 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(c("#C86C0D","#4558A3"),0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c("#C86C0D","#4558A3"))
text(1.5,2,"P < 0.001",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C1"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C2"),"Vorinostat"]
wilcox.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("tTME-C1","tTME-C2"),c(24,12)))
tmp$mut <- factor(tmp$mut,levels = c("tTME-C1","tTME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Vorinostat HDAC inhibitor in target DAWT cohort regarding TME.pdf"),width = 4,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "Vorinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.4,"P = 0.024",cex = 1)
invisible(dev.off())

#-----------------------#
# redo clinical samples #
tmp1 <- wt3[which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD4:CD8 ratio","OS","OS.time","PFS","PFS.time")]
colnames(tmp1)[1] <- c("ratio")
tmp1$group <- ifelse(tmp1$ratio > median(tmp1$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = F,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os in 12 anaplatic wt out of 55 wt cohort regarding cd4_8 ratio.pdf"), width = 4, height = 4)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = F,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs in 12 anaplatic wt out of 55 wt cohort regarding cd4_8 ratio.pdf"), width = 4, height = 4)
print(p)
dev.off()


tmp2 <- wt3[-which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD4:CD8 ratio","OS","OS.time","PFS","PFS.time")]
colnames(tmp2)[1] <- c("ratio")
tmp2$group <- ifelse(tmp2$ratio > median(tmp2$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = F,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os in 43 non-anaplatic wt out of 55 wt cohort regarding cd4_8 ratio.pdf"), width = 4, height = 4)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = F,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs in 43 non-anaplatic wt out of 55 wt cohort regarding cd4_8 ratio.pdf"), width = 4, height = 4)
print(p)
dev.off()

# multicox
tmp <- wt3[,c("CD4:CD8 ratio","OS","OS.time","PFS","PFS.time","Overall Stage (review/local)","Risk (review/local)","Histology type (review/local)")]
colnames(tmp)[1] <- c("ratio")
tmp$ratio <- tmp$ratio/10
tmp$group <- factor(ifelse(tmp$ratio > median(tmp$ratio),"High","Low"),levels = c("Low","High"))
tmp$stage <- ifelse(tmp$`Overall Stage (review/local)` %in% c("I","II"),"I_II","III_IV")
tmp$risk <- factor(ifelse(tmp$`Risk (review/local)` %in% c("High risk","High Risk"),"HRisk","LRisk"),levels = c("LRisk","HRisk"))
tmp$histology <- ifelse(tmp$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)"),"AWT","NAWT")
tmp$histology <- factor(tmp$histology, levels = c("NAWT","AWT"))
cox1 <- summary(coxph(Surv(PFS.time,PFS) ~ ratio + risk + stage + histology, data = tmp))
cox2 <- summary(coxph(Surv(OS.time,OS) ~ ratio + risk + stage + histology, data = tmp))

summary(coxph(Surv(PFS.time,PFS) ~ stage, data = tmp))
summary(coxph(Surv(OS.time,OS) ~ risk, data = tmp))

# forest plot
mulcox <- data.frame(variable = rownames(cox1$conf.int),
                     HR = cox1$conf.int[,1],
                     lower.95CI = cox1$conf.int[,3],
                     upper.95CI = cox1$conf.int[,4],
                     p = cox1$coefficients[,5],
                     stringsAsFactors = F)

hrtable <- rbind(c("RFS",NA,NA,NA,NA),
                 mulcox)
tabletext <- cbind(c("Variable",hrtable$variable),
                   c("HR (95% CI)",paste0(round(as.numeric(hrtable$HR),2)," (",round(as.numeric(hrtable$lower.95CI),2),"-",round(as.numeric(hrtable$upper.95CI),2),")")),
                   c("pvalue",round(as.numeric(hrtable$p),3)))
tabletext[2,2:3] <- NA

pdf(file.path(fig.path,"forestplot of risk table of 55 WTs regarding RFS.pdf"), width = 8, height = 6)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))),
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=3,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-2,-1,0,1,2,3,4),
           lwd.xaxis=2,
           xlab="log(Hazard Ratio)",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "7" = gpar(lwd=1, col="black", lty=1)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1.2,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())

mulcox <- data.frame(variable = rownames(cox2$conf.int),
                     HR = cox2$conf.int[,1],
                     lower.95CI = cox2$conf.int[,3],
                     upper.95CI = cox2$conf.int[,4],
                     p = cox2$coefficients[,5],
                     stringsAsFactors = F)

hrtable <- rbind(c("OS",NA,NA,NA,NA),
                 mulcox)
tabletext <- cbind(c("Variable",hrtable$variable),
                   c("HR (95% CI)",paste0(round(as.numeric(hrtable$HR),2)," (",round(as.numeric(hrtable$lower.95CI),2),"-",round(as.numeric(hrtable$upper.95CI),2),")")),
                   c("pvalue",round(as.numeric(hrtable$p),3)))
tabletext[2,2:3] <- NA

pdf(file.path(fig.path,"forestplot of risk table of 55 WTs regarding OS.pdf"), width = 8, height = 6)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),#log2(HR)
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))), 
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=3,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-2,-1,0,1,2,3,4,5,6),
           lwd.xaxis=2,
           xlab="log(Hazard Ratio)",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "7" = gpar(lwd=1, col="black", lty=1)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1.2,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())

#-----------------#
# redo cell lines #
paths <- c("HALLMARK_P53_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "REACTOME_CHROMATIN_MODIFYING_ENZYMES","REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_HDACS_DEACETYLATE_HISTONES")
paths <- MSigDB[which(MSigDB$ont %in% paths),]
paths.sig <- list()
for (i in unique(paths$ont)) {
  paths.sig[[i]] <- intersect(toupper(paths[which(paths$ont == i),"gene"]),rownames(ccle.wt))
}
paths.gsva.ccle <- gsva(as.matrix(ccle.wt),
                        paths.sig,
                        method = "ssgsea")
paths.gsva.ccle <- paths.gsva.ccle[c("HALLMARK_P53_PATHWAY",
                                     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                     "REACTOME_CHROMATIN_MODIFYING_ENZYMES",
                                     "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
                                     "REACTOME_HDACS_DEACETYLATE_HISTONES"),]
hcs <- hclust(distanceMatrix(as.matrix(paths.gsva.ccle[,rownames(annCol.ccle)]), "euclidean"), "ward.D"); hcs.ccle.tme <- hcs
group <- cutree(hcs, k =2)
annCol.ccle$TME <- ifelse(group[rownames(annCol.ccle)] == 1,"cTME-C1","cTME-C2")
annCol.ccle$EZH2_targets <- as.numeric(scale(ezh2.path.ccle[1,rownames(annCol.ccle)]))
annColors.ccle[["TME"]] <- c("cTME-C1" = jco[2],"cTME-C2" = jco[1])
annColors.ccle[["EZH2_targets"]] <- viridisLite::inferno(64)
plotdata <- standarize.fun(paths.gsva.ccle[,rownames(annCol.ccle)],halfwidth = 1)
p1 <- pheatmap(as.matrix(plotdata),
               cluster_rows = F,
               cluster_cols = hcs.ccle.tme,
               border_color = "black",
               annotation_col = annCol.ccle[,c("TME","Histology","Age","Sex","Model_Type")],
               annotation_colors = annColors.ccle[c("TME","Histology","Age","Sex","Model_Type")],
               color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
               #color = viridisLite::inferno(64),
               #color = viridisLite::viridis(64),
               # color = greenred(64),
               treeheight_row = 15,
               treeheight_col = 15,
               cellwidth = 12,
               cellheight = 12,
               gaps_row = 3,
               show_colnames = T,
               show_rownames = T,
               name = "ssgsea\n[z-scored]",
               cutree_cols = 2,
               fontsize_col = 9)
p1
# pdf(file.path(fig.path, "supervised clustering using msigdb pathways on wilms tumors in ccle cohort before combat.pdf"), height=6,width = 10)
# draw(p, annotation_legend_side = "bottom", heatmap_legend_side = "left")
# invisible(dev.off())

# clustering using pathway genes
outTab <- NULL
for (i in names(paths.sig)) {
    tmp1 <- as.numeric(paths.gsva.ccle[i,RNA_Clust1.ccle])
    tmp2 <- as.numeric(paths.gsva.ccle[i,RNA_Clust2.ccle])
    outTab <- rbind.data.frame(outTab,
                               data.frame(path = i,
                                          gene = j,
                                          avgTME1 = mean(tmp1),
                                          avgTME2 = mean(tmp2),
                                          diff = abs(mean(tmp1)-mean(tmp2)),
                                          p = wilcox.test(tmp1,tmp2)$p.value,
                                          dirct = ifelse(mean(tmp1) > mean(tmp2),"TME1","TME2"),
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
}

outTab <- NULL
for (i in names(paths.sig)) {
  tmp <- paths.sig[[i]]
  for (j in tmp) {
    tmp1 <- as.numeric(ccle.wt[j,RNA_Clust1.ccle])
    tmp2 <- as.numeric(ccle.wt[j,RNA_Clust2.ccle])
    outTab <- rbind.data.frame(outTab,
                               data.frame(path = i,
                                          gene = j,
                                          avgTME1 = mean(tmp1),
                                          avgTME2 = mean(tmp2),
                                          diff = abs(mean(tmp1)-mean(tmp2)),
                                          p = wilcox.test(tmp1,tmp2)$p.value,
                                          dirct = ifelse(mean(tmp1) > mean(tmp2),"TME1","TME2"),
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
  }
}
outTab <- outTab[order(outTab$diff,decreasing = T),]
selgene1 <- outTab[which(outTab$p < 0.05 & outTab$dirct == "TME1" & outTab$path == "HALLMARK_P53_PATHWAY"),][1:5,]
selgene2 <- outTab[which(outTab$p < 0.05 & outTab$dirct == "TME1" & outTab$path == "HALLMARK_INTERFERON_GAMMA_RESPONSE"),][1:5,]
selgene3 <- outTab[which(outTab$p < 0.05 & outTab$dirct == "TME1" & outTab$path == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),][1:5,]
selgene4 <- outTab[which(outTab$p < 0.5 & outTab$dirct == "TME2" & outTab$path == "REACTOME_CHROMATIN_MODIFYING_ENZYMES"),][1:5,]
selgene5 <- outTab[which(outTab$p < 0.5 & outTab$dirct == "TME2" & outTab$path == "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR"),][1:5,]
selgene6 <- outTab[which(outTab$p < 1 & outTab$dirct == "TME2" & outTab$path == "REACTOME_HDACS_DEACETYLATE_HISTONES"),][2:6,]

plotdata <- ccle.wt[c(selgene1$gene,
                      selgene2$gene,
                      selgene3$gene,
                      selgene4$gene,
                      selgene5$gene,
                      selgene6$gene),]
annRow <- data.frame(class = rep(c("HALLMARK_P53_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                   "REACTOME_CHROMATIN_MODIFYING_ENZYMES","REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_HDACS_DEACETYLATE_HISTONES"),
                                 c(5,5,5,5,5,5)),
                     row.names = rownames(plotdata))
annColors.ccle[["class"]] = c("HALLMARK_P53_PATHWAY" = darkblue,
                              "HALLMARK_INTERFERON_GAMMA_RESPONSE" = alpha(darkblue,0.6),
                              "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" =alpha(darkblue,0.3),
                              "REACTOME_HDACS_DEACETYLATE_HISTONES" = darkred,
                              "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR" = alpha(darkred,0.6),
                              "REACTOME_CHROMATIN_MODIFYING_ENZYMES" = alpha(darkred,0.3))

plotdata <- standarize.fun(plotdata,halfwidth = 2)
p2 <- pheatmap(as.matrix(plotdata),
               cluster_rows = F,
               cluster_cols = hcs.ccle.tme,
               border_color = "black",
               annotation_row = annRow,
               annotation_colors = annColors.ccle["class"],
               #color = NMF:::ccRamp(x = heatmap.L.BlYlRd,n = 64),
               # color = viridisLite::inferno(64),
               #color = viridisLite::viridis(64),
               # color = greenred(64),
               color = colorpanel(64,low=blue,mid = "black",high=gold),
               treeheight_row = 15,
               treeheight_col = 15,
               cellwidth = 12,
               cellheight = 12,
               gaps_row = c(5,10,15,20,25),
               show_colnames = T,
               show_rownames = T,
               name = "Expr.\n[z-scored]",
               cutree_cols = 2,
               fontsize_col = 9)
pdf(file.path(fig.path, "supervised clustering using msigdb pathways on wilms tumors in ccle cohort before combat2.pdf"), height=12,width = 10)
draw(p1 %v% p2, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
invisible(dev.off())

#--------------------------------------------#
# comparison of drug sensitivity between tme #
tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),2]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),2]
t.test(tmp1,tmp2, alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD6739 ATR inhibitor in gaby WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "AZD6738\nEstimated ln(IC50)",
        main = "ATR inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,2.8,"P = 0.210",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),5]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),5]
t.test(tmp1,tmp2, alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD1775 WEE1 inhibitor in gaby WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "AZD1775\nEstimated ln(IC50)",
        main = "WEE1 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.18,"P = 0.011",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C1"),2]
tmp2 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C2"),2]
t.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("C1","C2"),c(24,12)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD6739 ATR inhibitor in target DAWT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "AZD6738\nEstimated ln(IC50)",
        main = "ATR inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,2.0,"P = 0.75",cex = 3.5)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C1"),5]
tmp2 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C2"),5]
t.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("C1","C2"),c(24,12)))
tmp$mut <- factor(tmp$mut,levels = c("C1","C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of AZD1775 WEE1 inhibitor in target DAWT cohort regarding TME.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "AZD1775\nEstimated ln(IC50)",
        main = "WEE1 inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.5,"P = 0.046",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Wild"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.gaby[which(gdsc.pred.auc.gaby$TME == "TP53_Mutated"),"Vorinostat"]
t.test(tmp1,tmp2,alternative = "greater")
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("TME-C1","TME-C2"),c(5,5)))
tmp$mut <- factor(tmp$mut,levels = c("TME-C1","TME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Vorinostat HDAC inhibitor in gaby WT cohort regarding TME.pdf"),width = 3.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "Vorinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.5,"P = 0.039",cex = 1)
invisible(dev.off())

tmp1 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C1"),"Vorinostat"]
tmp2 <- gdsc.pred.auc.target.dawt[which(gdsc.pred.auc.target.dawt$TME == "C2"),"Vorinostat"]
t.test(tmp1,tmp2)
tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("tTME-C1","tTME-C2"),c(24,12)))
tmp$mut <- factor(tmp$mut,levels = c("tTME-C1","tTME-C2"))
pdf(file.path(fig.path,"boxplot for estimated lnIC50 of Vorinostat HDAC inhibitor in target DAWT cohort regarding TME.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,4,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "",
        ylab = "Vorinostat\nEstimated ln(IC50)",
        main = "HDAC inhibitor sensitivity",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = jco[2:1])
text(1.5,1.4,"P = 0.0099",cex = 1)
invisible(dev.off())

#-------------------------#
# statistic recalculation #
tables1 <- read.delim(file = file.path(data.path,"table s1 21samples.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

tmp1 <- tables1[which(tables1$`Histology type` == "ANA diffuse" & tables1$`Overall Stage` %in% c("I","II","III")),]
tmp1$TP53 <- ifelse(tmp1$`TP53 gene (NM_001276697) mutations` != "No",1,0)

fitd <- survdiff(Surv(OS.time, OS) ~ TP53, data=tmp1, na.action=na.exclude)
p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) 
fit <- survfit(Surv(OS.time, OS)~ TP53, data=tmp1, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
p <- ggsurvplot(fit, conf.int=F,
                risk.table=T, risk.table.col="strata",
                palette = "jama",
                pval = T,
                data=tmp,
                size=1,
                #legend = "none",
                tables.height = 0.3,
                surv.median.line = "hv",
                xlab = "Time (Days)",
                ylab = "Survival probability (%)",
                risk.table.y.text = F)
p$plot <- p$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),labels = seq(0,100,25))

#------------------------------------#
# cGAS-STING related immune pathways #
cgassting <- list(KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY = MSigDB[which(MSigDB$ont == "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY"),"gene"],
                  REACTOME_INNATE_IMMUNE_SYSTEM = MSigDB[which(MSigDB$ont == "REACTOME_INNATE_IMMUNE_SYSTEM"),"gene"])

gaby.cgassting <- GSVA::gsva(as.matrix(gaby.wt[,rownames(cluster2.wt.gaby.mRNA.ans$annCol)]),
                             gset.idx.list = cgassting,
                             method = "ssgsea")
wilcox.test(gaby.cgassting[1,rownames(annCol.wt.gaby.genomic.alter)]~annCol.wt.gaby.genomic.alter$RNA_Clust)
wilcox.test(gaby.cgassting[2,rownames(annCol.wt.gaby.genomic.alter)]~annCol.wt.gaby.genomic.alter$RNA_Clust)


target.cgassting <- GSVA::gsva(as.matrix(target.wt[,rownames(annCol.target.dawt)]),
                               gset.idx.list = cgassting,
                               method = "ssgsea")
wilcox.test(target.cgassting[1,rownames(annCol.target.dawt)]~annCol.target.dawt$TME)
wilcox.test(target.cgassting[2,rownames(annCol.target.dawt)]~annCol.target.dawt$TME)

#-------------------------------------#
# cibersortx for immune deconvolution #
tmp <- 2^gaby.wt-1
write.table(tmp, file = file.path(res.path,"gaby_wt_cibersort_input.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

tmp <- 2^target.wt-1
write.table(tmp, file = file.path(res.path,"target_wt_cibersort_input.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# ciber.gaby <- read.delim(file = file.path(res.path,"CIBERSORTx_Wilms_Gaby.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# ciber.target <- read.delim(file = file.path(res.path,"CIBERSORTx_Wilms_Target.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber.gaby <- read.delim(file = file.path(res.path,"CIBERSORTx_SFCE_WT_LM22.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber.target <- read.delim(file = file.path(res.path,"CIBERSORTx_TARGET_WT_LM22.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

ciber.gaby <- ciber.gaby[,1:22]
ciber.target <- ciber.target[,1:22]

outTab.ciber <- NULL
for (i in colnames(ciber.gaby)) {
  
  tmp1 <- ciber.gaby[rownames(annCol.wt.gaby.genomic.alter[which(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Mutated"),]),i]
  tmp2 <- ciber.gaby[rownames(annCol.wt.gaby.genomic.alter[which(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Wild"),]),i]
  
  outTab.ciber <- rbind.data.frame(outTab.ciber,
                                   data.frame(cell = i,
                                              avg.cold = mean(tmp1),
                                              avg.hot = mean(tmp2),
                                              p = wilcox.test(tmp1,tmp2)$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
write.table(outTab.ciber, file = file.path(res.path,"comparison of cibersort deconvolusion between wt subtypes in gaby cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

outTab.ciber <- NULL
for (i in colnames(ciber.target)) {
  
  tmp1 <- ciber.target[rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C1"),]),i]
  tmp2 <- ciber.target[rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C2"),]),i]
  
  outTab.ciber <- rbind.data.frame(outTab.ciber,
                                   data.frame(cell = i,
                                              avg.cold = mean(tmp2),
                                              avg.hot = mean(tmp1),
                                              p = wilcox.test(tmp1,tmp2)$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
write.table(outTab.ciber, file = file.path(res.path,"comparison of cibersort deconvolusion between wt subtypes in target dawt cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

# TR4 was obtained from FACS-sorted profiles of epithelial cells (EPCAM+), fibroblasts (CD10+), endothelial cells (CD31+) and immune cells (CD45+), obtained from lung cancer specimens
ciber.gaby <- read.delim(file = file.path(res.path,"CIBERSORTx_SFCE_WT_TR4.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber.target <- read.delim(file = file.path(res.path,"CIBERSORTx_TARGET_WT_TR4.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber.gaby <- ciber.gaby[,1:4]
ciber.target <- ciber.target[,1:4]

outTab.ciber <- NULL
for (i in colnames(ciber.gaby)) {
  
  tmp1 <- ciber.gaby[rownames(annCol.wt.gaby.genomic.alter[which(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Mutated"),]),i]
  tmp2 <- ciber.gaby[rownames(annCol.wt.gaby.genomic.alter[which(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Wild"),]),i]
  
  outTab.ciber <- rbind.data.frame(outTab.ciber,
                                   data.frame(cell = i,
                                              avg.cold = mean(tmp1),
                                              avg.hot = mean(tmp2),
                                              p = wilcox.test(tmp1,tmp2)$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
write.table(outTab.ciber, file = file.path(res.path,"comparison of tr4 deconvolusion between wt subtypes in gaby cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

tmp1 <- ciber.gaby[rownames(annCol.wt.gaby.genomic.alter[which(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Mutated"),]),"CD45"]
tmp2 <- ciber.gaby[rownames(annCol.wt.gaby.genomic.alter[which(annCol.wt.gaby.genomic.alter$RNA_Clust == "TP53_Wild"),]),"CD45"]

tmp <- data.frame(exp = c(tmp2,tmp1),
                  mut = rep(c("iWT","dWT"),c(length(tmp2),length(tmp1))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for CD45 immune proportion in gaby WT regarding tme subtypes.pdf"),width = 4,height = 3)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Immune cells (CD45+)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.0159",cex = 1)
invisible(dev.off())

outTab.ciber <- NULL
for (i in colnames(ciber.target)) {
  
  tmp1 <- ciber.target[rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C1"),]),i]
  tmp2 <- ciber.target[rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C2"),]),i]
  
  outTab.ciber <- rbind.data.frame(outTab.ciber,
                                   data.frame(cell = i,
                                              avg.cold = mean(tmp2),
                                              avg.hot = mean(tmp1),
                                              p = wilcox.test(tmp1,tmp2)$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
write.table(outTab.ciber, file = file.path(res.path,"comparison of tr4 deconvolusion between wt subtypes in target dawt cohort.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

tmp1 <- ciber.target[rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C1"),]),"CD45"]
tmp2 <- ciber.target[rownames(annCol.target.dawt[which(annCol.target.dawt$TME == "C2"),]),"CD45"]

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for CD45 immune proportion in target WT regarding tme subtypes.pdf"),width = 2,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "Immune cells (CD45+)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.2,"P = 0.0007",cex = 1)
invisible(dev.off())

dwls.res <- ciber.gaby[c(RNA_Clust1,RNA_Clust2),c(1,2,3,4)]
colnames(dwls.res) <- c("Fibroblast (CD10+)",  "Endothelial (CD31+)",  "Bulk immune cell (CD45+)",  "Epithelial (EPCAM+)")
dwls.res$ID <- c(RNA_Clust1,RNA_Clust2)
df_long <- dwls.res %>%
  pivot_longer(cols = -ID, names_to = "Feature", values_to = "Value")
df_long$Feature <- factor(df_long$Feature, levels = c("Bulk immune cell (CD45+)","Epithelial (EPCAM+)","Fibroblast (CD10+)","Endothelial (CD31+)"))
df_long$ID <- factor(df_long$ID, levels = c(RNA_Clust1,RNA_Clust2))
ggplot(df_long, aes(x = ID, y = Value, fill = Feature)) +
  geom_col(position = position_fill(), width = 1) +
  labs(x = "", y = "Relative Abundance") +
  scale_fill_manual(values = c("#E64B35","#4DBBD5","#00A087","#3C5488")) +
  theme_classic() +
  coord_fixed(ratio = (0.9 * nrow(df_long)) / ncol(df_long))
ggsave(filename = file.path(fig.path,"stacked plot for 4 cell types in each sample in gaby using cibersortx.pdf"),width = 15,height = 5)

dwls.res <- ciber.target[c(RNA_Clust1.target.dawt,RNA_Clust2.target.dawt),c(1,2,3,4)]
colnames(dwls.res) <- c("Fibroblast (CD10+)",  "Endothelial (CD31+)",  "Bulk immune cell (CD45+)",  "Epithelial (EPCAM+)")
dwls.res$ID <- c(RNA_Clust1.target.dawt,RNA_Clust2.target.dawt)
df_long <- dwls.res %>%
  pivot_longer(cols = -ID, names_to = "Feature", values_to = "Value")
df_long$Feature <- factor(df_long$Feature, levels = c("Bulk immune cell (CD45+)","Epithelial (EPCAM+)","Fibroblast (CD10+)","Endothelial (CD31+)"))
df_long$ID <- factor(df_long$ID, levels = c(RNA_Clust1.target.dawt,RNA_Clust2.target.dawt))
ggplot(df_long, aes(x = ID, y = Value, fill = Feature)) +
  geom_col(position = position_fill(), width = 1) +
  labs(x = "", y = "Relative Abundance") +
  scale_fill_manual(values = c("#E64B35","#4DBBD5","#00A087","#3C5488")) +
  theme_classic() +
  coord_fixed(ratio = (0.9 * nrow(df_long)) / ncol(df_long))
ggsave(filename = file.path(fig.path,"stacked plot for 4 cell types in each sample in target using cibersortx.pdf"),width = 15,height = 5)

# Calculate the median follow-up time
library(readxl)
library(prodlim)
data <- read_xlsx("E:/IGBMC/myproject/Wilms/Manuscript/New version/Submission/Nat Commun/Supplementary Tables-2023-01-13.xlsx", sheet = 1, skip = 1)
quantile(prodlim(Hist(time = OS.time, event = OS)~1, data = data, reverse = T))

# submap to compare expression similarity
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

indata <- gaby.wt
indata <- indata[apply(indata, 1, function(x) sum(x > 1) > 0.9*ncol(indata)),]
gaby.var <- rownames(indata)

indata <- target.wt[,rownames(annCol.target.dawt)]
indata <- indata[apply(indata, 1, function(x) sum(x > 1) > 0.9*ncol(indata)),]
target.var <- rownames(indata)

indata <- ccle.wt
indata <- indata[apply(indata, 1, function(x) sum(x > 1) > 0.9*ncol(indata)),]
ccle.var <- rownames(indata)
GENELIST <- intersect(intersect(intersect(ccle.var,gaby.var),target.var), Ginfo[which(Ginfo$genetype == "protein_coding"),"genename"])

ccle.info <- annCol.ccle
ccle.info <- ccle.info[order(ccle.info$TME,decreasing = F),]
ccle.info$rank <- rep(c(1,2),times=as.character(table(ccle.info$TME)))
sam_info <- ccle.info
in_gct <- ccle.wt[GENELIST,rownames(sam_info)]
gct_file <- file.path(res.path,"ccle.tme.for.SubMap.gct")
cls_file <- file.path(res.path,"ccle.tme.for.SubMap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

gaby.info <- annCol.wt.gaby.genomic.alter
gaby.info <- gaby.info[order(gaby.info$RNA_Clust,decreasing = T),]
gaby.info$rank <- rep(c(1,2),times=as.character(table(gaby.info$RNA_Clust)))
sam_info <- gaby.info
in_gct <- gaby.wt[GENELIST,rownames(sam_info)]
gct_file <- file.path(res.path,"gaby.tme.for.SubMap.gct")
cls_file <- file.path(res.path,"gaby.tme.for.SubMap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

target.info <- annCol.target.dawt
target.info <- target.info[order(target.info$TME,decreasing = T),]
target.info$rank <- rep(c(1,2),times=as.character(table(target.info$TME)))
sam_info <- target.info
in_gct <- target.wt[GENELIST,rownames(sam_info)]
gct_file <- file.path(res.path,"target.tme.for.SubMap.gct")
cls_file <- file.path(res.path,"target.tme.for.SubMap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# $Bonferroni.SA.matrix
# B1         B2
# A1 0.003996004 1.00000000
# A2 1.000000000 0.04795205
# 
# $FDR.SA.matrix
# B1         B2
# A1 0.003996004 1.00000000
# A2 0.999000999 0.02397602
# 
# $nominal.p.matrix.Fisher
# B1         B2
# A1 0.000999001 0.99600400
# A2 0.999000999 0.01198801

tmp <- matrix(c(0.001,0.996, # nominal p value
                0.999,0.012, 

                0.004,1, # Bonferroni adjusted p value
                1,0.048,

                0.004,1, # FDR adjusted p value
                0.999,0.024), 
              nrow = 6,byrow = T,
              dimnames = list(c("SFCE-iWT","SFCE-dWT","SFCE-iWT ","SFCE-dWT ","SFCE-iWT  ","SFCE-dWT  "),c("CCLE-WT (no anaplasia)","CCLE-WT (anaplasia)")))

hm <- pheatmap(tmp, 
               border_color = "black",
               number_format = "%.3f",
               cellwidth = 40, cellheight = 40,
               cluster_rows = F,cluster_cols = F,
               color = rev(NMF:::ccRamp(c("#E6EAF7","#B6D1E8","#498EB9","#204F8D"),64)),
               display_numbers = T,
               number_color = "black",
               fontsize_number = 10,
               name = "Statitic",
               annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni adjusted","Bonferroni adjusted","FDR adjusted","FDR adjusted"),
                                           row.names = rownames(tmp)),
               annotation_colors = list(pvalue=c("Nominal p value"="black","Bonferroni adjusted"="grey40","FDR adjusted" = "grey70")))
pdf(file.path(fig.path,"submap heatmap of similarity between gaby cohort and ccle.pdf"),width = 8,height = 8)
draw(hm, heatmap_legend_side = "left",annotation_legend_side = "right")
invisible(dev.off())

# $Bonferroni.SA.matrix
# B1        B2
# A1 0.03996004 1.0000000
# A2 1.00000000 0.1438561
# 
# $FDR.SA.matrix
# B1         B2
# A1 0.03996004 1.00000000
# A2 1.00000000 0.07192807
# 
# $nominal.p.matrix.Fisher
# B1         B2
# A1 0.00999001 1.00000000
# A2 0.97502498 0.03596404
tmp <- matrix(c(0.001,1, # nominal p value
                0.975,0.036, 
                
                0.040,1, # Bonferroni adjusted p value
                1,0.144,
                
                0.040,1, # FDR adjusted p value
                1,0.072), 
              nrow = 6,byrow = T,
              dimnames = list(c("SFCE-iWT","SFCE-dWT","SFCE-iWT ","SFCE-dWT ","SFCE-iWT  ","SFCE-dWT  "),c("TARGET-iWT","TARGET-dWT")))

hm <- pheatmap(tmp, 
               border_color = "black",
               number_format = "%.3f",
               cellwidth = 40, cellheight = 30,
               cluster_rows = F,cluster_cols = F,
               color = rev(NMF:::ccRamp(c("#E6EAF7","#B6D1E8","#498EB9","#204F8D"),64)),
               display_numbers = T,
               number_color = "black",
               fontsize_number = 10,
               name = "Statitic",
               annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni adjusted","Bonferroni adjusted","FDR adjusted","FDR adjusted"),
                                           row.names = rownames(tmp)),
               annotation_colors = list(pvalue=c("Nominal p value"="black","Bonferroni adjusted"="grey40","FDR adjusted" = "grey70")))
pdf(file.path(fig.path,"submap heatmap of similarity between gaby cohort and target.pdf"),width = 8,height = 8)
draw(hm, heatmap_legend_side = "left",annotation_legend_side = "right")
invisible(dev.off())


# immune cells analysis
wt3 <- read.table(file.path(data.path,"clinical data for 55 wts - 27052023.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp <- wt182[rownames(wt3),"EZH2 S21"]; wt3$EZH2 <- tmp
wt3$EZH2.curated <- as.numeric(wt3$EZH2)
#wt3$EZH2.curated[38] <- 50
#wt3$EZH2.curated[53] <- 90

cor.test(wt3$EZH2.curated,wt3$CD3_10HPF,method = "spearman",alternative = "less")
cor.test(wt3$EZH2.curated,wt3$CD8_10HPF,method = "spearman",alternative = "less")

mycol2 <- brewer.pal(12, "Paired")

tmp <- data.frame(count = c(wt3$CD3_10HPF, wt3$CD8_10HPF),
                  gene = c(wt3$EZH2.curated,wt3$EZH2.curated),
                  lymphocyte = rep(c("CD3","CD8"),each = 55),
                  cohort = "gaby")
tmp$lymphocyte <- factor(tmp$lymphocyte, levels = c("CD3","CD8"))
ggplot(tmp, aes(x=gene, y=count, shape = lymphocyte, color = lymphocyte)) + 
  geom_point(size = 3)+
  geom_smooth(method = "lm", se = FALSE,aes(color = lymphocyte)) +
  # geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0, # or, use fill = NA
  #             aes(color = group, fill = group),linetype = "dashed") +
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0.2,
              aes(color = lymphocyte, fill = lymphocyte),linetype = "dashed") +
  scale_color_manual(values = c(darkblue,darkred)) +
  scale_fill_manual(values = c(darkblue,darkred)) +
  scale_shape_manual(values = c(19,17)) + 
  scale_y_continuous(limits=c(0,700), breaks=c(0,100,200,300,400,500,600,700)) + 
  scale_x_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100)) + 
  theme_bw() + 
  ylab("Lymphocyte per 10 HPF") +
  xlab("EZH2") + theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  annotate(geom="text", x=75, y=500, label="CD3: R = -0.22, p = 0.073",
           color=darkblue, size = 5) +
  annotate(geom="text", x=75, y=470, label="CD8: R = -0.36, p = 0.007",
           color=darkred, size = 5)
ggsave(filename = file.path(fig.path,"correlation scatter plot of lymphocyte count and ezh2 in 55 wts.pdf"), width = 5,height = 5)

ggplot(tmp, aes(x=gene, y=count, shape = lymphocyte, color = lymphocyte)) + 
  geom_point(size = 3)+
  geom_smooth(method = "lm", se = FALSE,aes(color = lymphocyte)) +
  # geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0, # or, use fill = NA
  #             aes(color = group, fill = group),linetype = "dashed") +
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0.2,
              aes(color = lymphocyte, fill = lymphocyte),linetype = "dashed") +
  scale_color_manual(values = c(darkblue,darkred)) +
  scale_fill_manual(values = c(darkblue,darkred)) +
  scale_shape_manual(values = c(19,17)) + 
  scale_y_continuous(limits=c(0,1000), breaks=c(0,100,200,300,400,500,1000)) + 
  scale_x_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100)) + 
  theme_bw() + 
  ylab("Lymphocyte per 10 HPF") +
  xlab("EZH2") + theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  annotate(geom="text", x=75, y=500, label="CD3: R = -0.22, p = 0.073",
           color=darkblue, size = 5) +
  annotate(geom="text", x=75, y=470, label="CD8: R = -0.36, p = 0.007",
           color=darkred, size = 5)
ggsave(filename = file.path(fig.path,"correlation scatter plot of lymphocyte count and ezh2 in 55 wts with long y-axis.pdf"), width = 5,height = 5)

# do the correlation differently
tmp1 <- wt3[which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","EZH2.curated","OS","OS.time","PFS","PFS.time")]
cor.test(tmp1$EZH2.curated,tmp1$CD3_10HPF,method = "pearson",alternative = "less")
cor.test(tmp1$EZH2.curated,tmp1$CD8_10HPF,method = "pearson",alternative = "less")

tmp <- data.frame(count = c(tmp1$CD3_10HPF, tmp1$CD8_10HPF),
                  gene = c(tmp1$EZH2.curated,tmp1$EZH2.curated),
                  lymphocyte = rep(c("CD3","CD8"),each = 12),
                  cohort = "gaby")
tmp$lymphocyte <- factor(tmp$lymphocyte, levels = c("CD3","CD8"))
ggplot(tmp, aes(x=gene, y=count, shape = lymphocyte, color = lymphocyte)) + 
  geom_point(size = 3)+
  geom_smooth(method = "lm", se = FALSE,aes(color = lymphocyte)) +
  # geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0, # or, use fill = NA
  #             aes(color = group, fill = group),linetype = "dashed") +
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0.2,
              aes(color = lymphocyte, fill = lymphocyte),linetype = "dashed") +
  scale_color_manual(values = c(blue,red)) +
  scale_fill_manual(values = c(blue,red)) +
  scale_shape_manual(values = c(19,17)) + 
  scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) + 
  scale_x_continuous(limits=c(0,80), breaks=c(0,20,40,60,80)) + 
  theme_bw() + 
  ylab("Lymphocyte per 10 HPF") +
  xlab("EZH2") + theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  annotate(geom="text", x=55, y=500, label="CD3: R = -0.15, p = 0.336",
           color=blue, size = 5) +
  annotate(geom="text", x=55, y=470, label="CD8: R = -0.17, p = 0.323",
           color=red, size = 5)
ggsave(filename = file.path(fig.path,"correlation scatter plot of lymphocyte count and ezh2 in 12 anaplastic out of 55 wts.pdf"), width = 5,height = 5)

tmp1 <- wt3[-which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","EZH2.curated","OS","OS.time","PFS","PFS.time")]
cor.test(tmp1$EZH2.curated,tmp1$CD3_10HPF,method = "spearman",alternative = "less")
cor.test(tmp1$EZH2.curated,tmp1$CD8_10HPF,method = "spearman",alternative = "less")

tmp <- data.frame(count = c(tmp1$CD3_10HPF, tmp1$CD8_10HPF),
                  gene = c(tmp1$EZH2.curated,tmp1$EZH2.curated),
                  lymphocyte = rep(c("CD3","CD8"),each = 43),
                  cohort = "gaby")
tmp$lymphocyte <- factor(tmp$lymphocyte, levels = c("CD3","CD8"))
ggplot(tmp, aes(x=gene, y=count, shape = lymphocyte, color = lymphocyte)) + 
  geom_point(size = 3)+
  geom_smooth(method = "lm", se = FALSE,aes(color = lymphocyte)) +
  # geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0, # or, use fill = NA
  #             aes(color = group, fill = group),linetype = "dashed") +
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0.2,
              aes(color = lymphocyte, fill = lymphocyte),linetype = "dashed") +
  scale_color_manual(values = c(blue,red)) +
  scale_fill_manual(values = c(blue,red)) +
  scale_shape_manual(values = c(19,17)) + 
  scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) + 
  scale_x_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100)) + 
  theme_bw() + 
  ylab("Lymphocyte per 10 HPF") +
  xlab("EZH2") + theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.title.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  annotate(geom="text", x=75, y=500, label="CD3: R = -0.25, p = 0.074",
           color=blue, size = 5) +
  annotate(geom="text", x=75, y=470, label="CD8: R = -0.42, p = 0.005",
           color=red, size = 5)
ggsave(filename = file.path(fig.path,"correlation scatter plot of lymphocyte count and ezh2 in 43 non-anaplastic out of 55 wts.pdf"), width = 5,height = 5)

tmp <- data.frame(exp = c(wt3$`CD3_10HPF`,wt3$`CD8_10HPF`),
                  mut = rep(c("CD3","CD8"),c(55,55)),
                  sam = rep(rownames(wt3),2))

pdf(file.path(fig.path,"boxplot for lymphcytes of new 55 wts.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(3,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data=tmp,
        xlab = "",
        ylab = "Number of cells/10 HPF",
        names = c("CD3","CD8"),
        outline = F,
        ylim = c(0,1800),
        col = c(alpha(red,0.8),alpha(blue,0.8)))
stripchart(exp~mut,data= tmp, vertical = TRUE, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col =  c(alpha(red,0.8),alpha(blue,0.8)))
invisible(dev.off())

par(bty="o", mgp = c(1.5,.33,0), mar=c(3,3,1,0.1), las=1, tcl=-.25,las = 1)
barplot(wt3$`CD3_10HPF`,col = alpha(red,0.6), border = NA)
barplot(wt3$`CD8_10HPF`,add = T, col = alpha(blue,0.6), border = NA)
#tmp2 <- tmp[which(tmp$mut %in% c("CD4","CD8")),]
tmp2 <- tmp
samorder <- rownames(wt3)
tmp2$sam <- factor(tmp2$sam,levels = samorder)
ggplot (tmp2, aes(x=sam, y=exp, fill=mut)) + 
  geom_bar (stat="identity", position = position_dodge(width = 0.8)) +   
  scale_fill_manual(values = c(alpha(red,0.8),alpha(blue,0.8),alpha(jco[2],0.8))) +
  theme_classic() + ylab("Number of cells/10 HPF") + 
  theme(axis.line.y = element_line(size = 0.8),
        axis.ticks.y = element_line(size = 0.2),
        axis.text.x = element_text(size = 8, color = "black", angle = 90),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(
    breaks = seq(0,1800,100),
    labels = seq(0,1800,100)
  )
ggsave(file.path(fig.path,"distribution of lymphcytes in new 55 wts.pdf"),width = 6,height = 6)

# differencec of EZH2 and TP53 mutation
data <- read_xlsx("E:/IGBMC/myproject/Wilms/Manuscript/New version/Submission/Nat Commun/Supplementary Tables-R1.xlsx", sheet = 1, skip = 1)

tmp <- data.frame(ezh2 = data$EZH2,
                  tp53 = ifelse(data$`TP53 gene (NM_001276697) mutations` == "No","no","yes"),
                  row.names = data$`n siop`)
mean(na.omit(tmp[rownames(tp53.yes),"ezh2"]))
mean(na.omit(tmp[rownames(tp53.no),"ezh2"]))

tmp <- na.omit(tmp)
wilcox.test(tmp$ezh2~tmp$tp53)

# check histology for 53 cases
wt3.curated <- wt3[!(is.na(wt3$CD3_10HPF) & is.na(wt3$CD8_10HPF)),]
table(wt3.curated$`Histology type (review/local)`)

# use median cutoff
## anaplastic tumors
tmp1 <- wt3[which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
coxph(Surv(OS.time, OS) ~ CD3_10HPF, data = tmp1)
coxph(Surv(PFS.time, PFS) ~ CD3_10HPF, data = tmp1)
coxph(Surv(OS.time, OS) ~ CD8_10HPF, data = tmp1)
coxph(Surv(PFS.time, PFS) ~ CD8_10HPF, data = tmp1)
colnames(tmp1)[1] <- c("ratio")
tmp1 <- as.data.frame(na.omit(tmp1))
tmp1$group <- ifelse(tmp1$ratio > median(tmp1$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using median cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using median cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

tmp1 <- wt3[which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
colnames(tmp1)[2] <- c("ratio")
tmp1 <- as.data.frame(na.omit(tmp1))
tmp1$group <- ifelse(tmp1$ratio > median(tmp1$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using median cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using median cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

## non-anaplastic tumors
tmp2 <- wt3[-which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
coxph(Surv(OS.time, OS) ~ CD3_10HPF, data = tmp2)
coxph(Surv(PFS.time, PFS) ~ CD3_10HPF, data = tmp2)
coxph(Surv(OS.time, OS) ~ CD8_10HPF, data = tmp2)
coxph(Surv(PFS.time, PFS) ~ CD8_10HPF, data = tmp2)

tmp2 <- tmp2[!is.na(tmp2$CD3_10HPF),]
colnames(tmp2)[1] <- c("ratio")
tmp2$group <- ifelse(tmp2$ratio > median(tmp2$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using median cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using median cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

tmp2 <- wt3[-which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
tmp2 <- tmp2[!is.na(tmp2$CD8_10HPF),]
colnames(tmp2)[2] <- c("ratio")
tmp2$group <- ifelse(tmp2$ratio > median(tmp2$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using median cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using median cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

# all 55 cases
tmp3 <- wt3[,c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
coxph(Surv(OS.time, OS) ~ CD3_10HPF, data = tmp3)
coxph(Surv(PFS.time, PFS) ~ CD3_10HPF, data = tmp3)
coxph(Surv(OS.time, OS) ~ CD8_10HPF, data = tmp3)
coxph(Surv(PFS.time, PFS) ~ CD8_10HPF, data = tmp3)

tmp3 <- tmp3[!is.na(tmp3$CD3_10HPF),]
colnames(tmp3)[1] <- c("ratio")
tmp3$group <- ifelse(tmp3$ratio > median(tmp3$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using median cutoff in updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using median cutoff in updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

tmp3 <- wt3[,c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
tmp3 <- tmp3[!is.na(tmp3$CD8_10HPF),]
colnames(tmp3)[2] <- c("ratio")
tmp3$group <- ifelse(tmp3$ratio > median(tmp3$ratio),"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using median cutoff in updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using median cutoff in updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

# use best cutoff
## anaplastic tumors
tmp1 <- wt3[which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
tmp1$OS.time <- tmp1$OS.time/12
tmp1$PFS.time <- tmp1$PFS.time/12
colnames(tmp1)[1] <- c("ratio")
tmp1 <- as.data.frame(na.omit(tmp1))
bestcutoff <- surv_cutpoint(data = tmp1, variables = "ratio",time = "OS.time",event = "OS",minprop = 0.2)
tmp1$group <- ifelse(tmp1$ratio > bestcutoff$cutpoint[1,1],"High","Low")
tmp1$TP53 <- NA
tmp1[intersect(rownames(tmp1),rownames(tp53.yes)),"TP53"] <- "yes"
tmp1[intersect(rownames(tmp1),rownames(tp53.no)),"TP53"] <- "no"
fisher.test(table(tmp1$TP53,tmp1$group),alternative = "greater")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Years)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using best cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 3.5, height = 5)
print(p)
dev.off()

bestcutoff <- surv_cutpoint(data = tmp1, variables = "ratio",time = "PFS.time",event = "PFS",minprop = 0.2)
tmp1$group <- ifelse(tmp1$ratio > bestcutoff$cutpoint[1,1],"High","Low")

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Years)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using best cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 3.5, height = 5)
print(p)
dev.off()

tmp1 <- wt3[which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
tmp1$OS.time <- tmp1$OS.time/12
tmp1$PFS.time <- tmp1$PFS.time/12
colnames(tmp1)[2] <- c("ratio")
tmp1 <- as.data.frame(na.omit(tmp1))
bestcutoff <- surv_cutpoint(data = tmp1, variables = "ratio",time = "OS.time",event = "OS",minprop = 0.2)
tmp1$group <- ifelse(tmp1$ratio > bestcutoff$cutpoint[1,1],"High","Low")
tmp1$TP53 <- NA
tmp1[intersect(rownames(tmp1),rownames(tp53.yes)),"TP53"] <- "yes"
tmp1[intersect(rownames(tmp1),rownames(tp53.no)),"TP53"] <- "no"
fisher.test(table(tmp1$TP53,tmp1$group),alternative = "greater")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Years)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using best cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 3.5, height = 5)
print(p)
dev.off()

bestcutoff <- surv_cutpoint(data = tmp1, variables = "ratio",time = "PFS.time",event = "PFS",minprop = 0.2)
tmp1$group <- ifelse(tmp1$ratio > bestcutoff$cutpoint[1,1],"High","Low")
fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp1,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp1,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Years)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using best cutoff in 12 anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 3.5, height = 5)
print(p)
dev.off()

## non-anaplastic tumors
tmp2 <- wt3[-which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
tmp2$OS.time <- tmp2$OS.time/12
tmp2$PFS.time <- tmp2$PFS.time/12
coxph(Surv(OS.time, OS) ~ CD3_10HPF, data = tmp2)
coxph(Surv(PFS.time, PFS) ~ CD3_10HPF, data = tmp2)
coxph(Surv(OS.time, OS) ~ CD8_10HPF, data = tmp2)
coxph(Surv(PFS.time, PFS) ~ CD8_10HPF, data = tmp2)

tmp2 <- tmp2[!is.na(tmp2$CD3_10HPF),]
colnames(tmp2)[1] <- c("ratio")
bestcutoff <- surv_cutpoint(data = tmp2, variables = "ratio",time = "OS.time",event = "OS",minprop = 0.2)
tmp2$group <- ifelse(tmp2$ratio > bestcutoff$cutpoint[1,1],"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Years)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using best cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

bestcutoff <- surv_cutpoint(data = tmp2, variables = "ratio",time = "PFS.time",event = "PFS",minprop = 0.2)
tmp2$group <- ifelse(tmp2$ratio > bestcutoff$cutpoint[1,1],"High","Low")

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Years)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using best cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

tmp2 <- wt3[-which(wt3$`Histology type (review/local)` %in% c("ANA diffuse (312)","Focal anaplasia (311)")),c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time")]
tmp2 <- tmp2[!is.na(tmp2$CD8_10HPF),]
tmp2$OS.time <- tmp2$OS.time/12
tmp2$PFS.time <- tmp2$PFS.time/12
colnames(tmp2)[2] <- c("ratio")
bestcutoff <- surv_cutpoint(data = tmp2, variables = "ratio",time = "OS.time",event = "OS",minprop = 0.2)
tmp2$group <- ifelse(tmp2$ratio > bestcutoff$cutpoint[1,1],"High","Low")

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using best cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

bestcutoff <- surv_cutpoint(data = tmp2, variables = "ratio",time = "PFS.time",event = "PFS",minprop = 0.2)
tmp2$group <- ifelse(tmp2$ratio > bestcutoff$cutpoint[1,1],"High","Low")

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp2,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp2,
                size              = 1,
                break.time.by     = 1,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Years)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using best cutoff in 43 non-anaplatic wt out of updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

# all 55 cases
tmp3 <- wt3[,c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time","Histology type (review/local)")]
coxph(Surv(OS.time, OS) ~ CD3_10HPF, data = tmp3)
coxph(Surv(PFS.time, PFS) ~ CD3_10HPF, data = tmp3)
coxph(Surv(OS.time, OS) ~ CD8_10HPF, data = tmp3)
coxph(Surv(PFS.time, PFS) ~ CD8_10HPF, data = tmp3)

tmp3 <- tmp3[!is.na(tmp3$CD3_10HPF),]
colnames(tmp3)[1] <- c("ratio")
bestcutoff <- surv_cutpoint(data = tmp3, variables = "ratio",time = "OS.time",event = "OS",minprop = 0.2)
tmp3$group <- ifelse(tmp3$ratio > bestcutoff$cutpoint[1,1],"High","Low")
table(tmp3$`Histology type (review/local)`,tmp3$group)
#                       High Low
# ANA diffuse (312)        6   4
# Blastemal (212)          5   0
# Epithelial (211)         1   2
# Focal anaplasia (311)    1   0
# Mixed (214)             16   2
# Regressive (216)         7   0
# Stromal (213)            4   2
fisher.test(matrix(c(6+1,4+0,5+1+16+7+4,2+2+2),byrow = T,nrow = 2))

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using best cutoff in updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

bestcutoff <- surv_cutpoint(data = tmp3, variables = "ratio",time = "PFS.time",event = "PFS",minprop = 0.2)
tmp3$group <- ifelse(tmp3$ratio > bestcutoff$cutpoint[1,1],"High","Low")
table(tmp3$`Histology type (review/local)`,tmp3$group)
#                       High Low
# ANA diffuse (312)        6   4
# Blastemal (212)          5   0
# Epithelial (211)         1   2
# Focal anaplasia (311)    1   0
# Mixed (214)             16   2
# Regressive (216)         7   0
# Stromal (213)            3   3
fisher.test(matrix(c(6+1,4+0,5+1+16+7+3,2+2+3),byrow = T,nrow = 2))

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using best cutoff in updated 55 wt cohort regarding cd3 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

tmp3 <- wt3[,c("CD3_10HPF","CD8_10HPF","OS","OS.time","PFS","PFS.time","Histology type (review/local)")]
tmp3 <- tmp3[!is.na(tmp3$CD8_10HPF),]
colnames(tmp3)[2] <- c("ratio")
bestcutoff <- surv_cutpoint(data = tmp3, variables = "ratio",time = "OS.time",event = "OS",minprop = 0.2)
tmp3$group <- ifelse(tmp3$ratio > bestcutoff$cutpoint[1,1],"High","Low")
table(tmp3$`Histology type (review/local)`,tmp3$group)
#                       High Low
# ANA diffuse (312)        3   7
# Blastemal (212)          3   2
# Epithelial (211)         3   0
# Focal anaplasia (311)    1   0
# Mixed (214)             12   7
# Regressive (216)         5   2
# Stromal (213)            5   2
fisher.test(matrix(c(3+1,7+0,3+3+12+5+5,2+7+2+2),byrow = T,nrow = 2))

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                ylim              = c(0,1),
                surv.median.line  = "hv",
                xlab              = "Time (Months)",
                ylab              = "Overall survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of os using best cutoff in updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

bestcutoff <- surv_cutpoint(data = tmp3, variables = "ratio",time = "PFS.time",event = "PFS",minprop = 0.2)
tmp3$group <- ifelse(tmp3$ratio > bestcutoff$cutpoint[1,1],"High","Low")
table(tmp3$`Histology type (review/local)`,tmp3$group)
#                       High Low
# ANA diffuse (312)        2   8
# Blastemal (212)          2   3
# Epithelial (211)         3   0
# Focal anaplasia (311)    1   0
# Mixed (214)             12   7
# Regressive (216)         5   2
# Stromal (213)            3   4
fisher.test(matrix(c(2+1,8+0,2+3+12+5+3,3+7+2+4),byrow = T,nrow = 2))

fitd <- survdiff(Surv(PFS.time, PFS) ~ group,
                 data      = tmp3,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(PFS.time, PFS)~ group,
               data      = tmp3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = T,
                risk.table.col    = "strata",
                palette           = c(darkred,darkblue),
                data              = tmp3,
                size              = 1,
                break.time.by     = 12,
                legend.title      = "",
                surv.median.line  = "hv",
                ylim              = c(0,1),
                xlab              = "Time (Months)",
                ylab              = "Recurrence-free survival",
                risk.table.y.text = FALSE)
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 4, y = 0.1,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file.path(fig.path,"km of pfs using best cutoff in updated 55 wt cohort regarding cd8 count.pdf"), width = 4, height = 5)
print(p)
dev.off()

#--------------------#
# all 182 wt samples #
wt182 <- read.delim(file = file.path(data.path,"clinical data for all 182 wts - 27052023.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 169 T has a negative recurrence time
# Column ag and ai to see prognostic relevance, clinical features, histology, stage, age.
tmp <- wt182[,c("Overall Stage (review/local)","Age at diagnostic (years)","Sexe","H3K36ME3","H3K4ME","rechute oui=1 non=0","TTP months","0= alive; 1= death","os (YEARS)")]
colnames(tmp) <- c("stage","age","sex","H3K36ME3","H3K4ME","PFS","PFS.time","OS","OS.time")

#-----------------------------------#
# compare PDL1, PD1, CTLA4 and LAG3 #

## SFCE
# PD1
tmp1 <- as.numeric(gaby.wt["PDCD1",RNA_Clust1])
tmp2 <- as.numeric(gaby.wt["PDCD1",RNA_Clust2])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for PD1 in Gaby WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "PDCD1 (PD1)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 1",cex = 1)
invisible(dev.off())

# PDL1
tmp1 <- as.numeric(gaby.wt["CD274",RNA_Clust1])
tmp2 <- as.numeric(gaby.wt["CD274",RNA_Clust2])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for PDL1 in Gaby WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "CD274 (PDL1)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.151",cex = 1)
invisible(dev.off())

# CTLA4
tmp1 <- as.numeric(gaby.wt["CTLA4",RNA_Clust1])
tmp2 <- as.numeric(gaby.wt["CTLA4",RNA_Clust2])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for CTLA4 in Gaby WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "CTLA4 (CD152)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.091",cex = 1)
invisible(dev.off())

# LAG3
tmp1 <- as.numeric(gaby.wt["LAG3",RNA_Clust1])
tmp2 <- as.numeric(gaby.wt["LAG3",RNA_Clust2])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for LAG3 in Gaby WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "LAG3 (CD223)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.222",cex = 1)
invisible(dev.off())

## TARGET
# PD1
tmp1 <- as.numeric(target.wt["PDCD1",RNA_Clust1.target.dawt])
tmp2 <- as.numeric(target.wt["PDCD1",RNA_Clust2.target.dawt])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for PD1 in target WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "PDCD1 (PD1)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.679",cex = 1)
invisible(dev.off())

# PDL1
tmp1 <- as.numeric(target.wt["CD274",RNA_Clust1.target.dawt])
tmp2 <- as.numeric(target.wt["CD274",RNA_Clust2.target.dawt])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for PDL1 in target WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "CD274 (PDL1)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.615",cex = 1)
invisible(dev.off())

# CTLA4
tmp1 <- as.numeric(target.wt["CTLA4",RNA_Clust1.target.dawt])
tmp2 <- as.numeric(target.wt["CTLA4",RNA_Clust2.target.dawt])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for CTLA4 in target WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "CTLA4 (CD152)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,0.3,"P = 0.090",cex = 1)
invisible(dev.off())

# LAG3
tmp1 <- as.numeric(target.wt["LAG3",RNA_Clust1.target.dawt])
tmp2 <- as.numeric(target.wt["LAG3",RNA_Clust2.target.dawt])
wilcox.test(tmp1,tmp2)

tmp <- data.frame(exp = c(tmp1,tmp2),
                  mut = rep(c("iWT","dWT"),c(length(tmp1),length(tmp2))))
tmp$mut <- factor(tmp$mut,levels = c("iWT","dWT"))
pdf(file.path(fig.path,"boxplot for LAG3 in target WT regarding tme subtypes.pdf"),width = 2.5,height = 3.5)
par(bty="o", mgp = c(1.5,.33,0), mar=c(4,3,1,0.1), las=1, tcl=-.25,las = 1)
boxplot(exp~mut,data = tmp,
        xlab = "Immunity status",
        ylab = "LAG3 (CD223)",
        outline = F,
        col = ggplot2::alpha(jco[2:1],0.7))
stripchart(exp~mut, vertical = TRUE, data = tmp, cex = 1,
           method = "jitter", add = TRUE, pch = 20, col = c(jco[2],jco[1]))
text(1.5,3.5,"P = 0.631",cex = 1)
invisible(dev.off())

#-------------------------------------------------------#
# differential expression for transcriptomic similarity #
target.wt.ct <- read.delim(file = file.path(data.path,"Wilm_TARGET_Project_rawCounts_V15.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
target.wt.ct <- target.wt.ct[rownames(Ginfo),]
target.wt.ct$gene <- Ginfo$genename
target.wt.ct <- target.wt.ct[!duplicated(target.wt.ct$gene),]
rownames(target.wt.ct) <- target.wt.ct$gene
target.wt.ct <- target.wt.ct[,setdiff(colnames(target.wt.ct),"gene")]

gaby.wt.ct <- target.wt.ct[,rownames(annCol.wt.gaby.genomic.alter)]
target.wt.ct <- target.wt.ct[,rownames(annCol.target.dawt)]

gaby.wt.ct <- gaby.wt.ct[rowSums(gaby.wt.ct) > 0,]
target.wt.ct <- target.wt.ct[rowSums(target.wt.ct) > 0,]
ccle.wt.ct <- ccle.wt.ct[rowSums(ccle.wt.ct) > 0,]

comgene <- intersect(Ginfo[which(Ginfo$genetype == "protein_coding"),"genename"], rownames(gaby.wt.ct))
comgene <- intersect(comgene, rownames(ccle.wt.ct))
comgene <- intersect(comgene, rownames(target.wt.ct))

# GABY
tmp <- annCol.wt.gaby.genomic.alter
tmp$TME <- ifelse(tmp$RNA_Clust == "TP53_Mutated","dWT","iWT")
group <- factor(tmp$TME, levels = c("dWT","iWT"))
design <- model.matrix(~group)
rownames(design) <- rownames(tmp)
y <- DGEList(counts=gaby.wt.ct[comgene,], group=group)
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
ordered_tags <- topTags(lrt, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
deg.gaby <- allDiff

# TARGECT
tmp <- annCol.target.dawt
tmp$TME <- ifelse(tmp$TME == "C1","iWT","dWT")
group <- factor(tmp$TME, levels = c("dWT","iWT"))
design <- model.matrix(~group)
rownames(design) <- rownames(tmp)
y <- DGEList(counts=target.wt.ct[comgene,], group=group)
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
ordered_tags <- topTags(lrt, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
deg.target <- allDiff

# CCLE
tmp <- annCol.ccle
tmp$TME <- ifelse(tmp$TME == "cTME-C1","iWT","dWT")
group <- factor(tmp$TME, levels = c("dWT","iWT"))
design <- model.matrix(~group)
rownames(design) <- rownames(tmp)
y <- DGEList(counts=ccle.wt.ct[comgene,], group=group)
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
ordered_tags <- topTags(lrt, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
deg.ccle <- allDiff
write.table(deg.gaby, file = file.path(res.path,"deg.gaby.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(deg.target, file = file.path(res.path,"deg.target.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(deg.ccle, file = file.path(res.path,"deg.ccle.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

up.gaby <- rownames(deg.gaby[which(deg.gaby$logFC > 1 & deg.gaby$PValue < 0.05 & deg.gaby$FDR < 0.25),])
dn.gaby <- rownames(deg.gaby[which(deg.gaby$logFC < -1 & deg.gaby$PValue < 0.05 & deg.gaby$FDR < 0.25),])

up.target <- rownames(deg.target[which(deg.target$logFC > 1 & deg.target$PValue < 0.05 & deg.target$FDR < 0.25),])
dn.target <- rownames(deg.target[which(deg.target$logFC < -1 & deg.target$PValue < 0.05 & deg.target$FDR < 0.25),])

up.ccle <- rownames(deg.ccle[which(deg.ccle$logFC > 1 & deg.ccle$PValue < 0.05 & deg.ccle$FDR < 0.25),])
dn.ccle <- rownames(deg.ccle[which(deg.ccle$logFC < -1 & deg.ccle$PValue < 0.05 & deg.ccle$FDR < 0.25),])

# save image
save.image(file = "WILMS.RData")
