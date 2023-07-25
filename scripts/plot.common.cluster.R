plot.common.cluster = function(indata=NULL, tumorname=NULL, N.cluster=NULL, N.bootstrap=NULL, N.gene.per.bootstrap=NULL, N.sample.per.bootstrap=NULL, map.res.path=NULL, fig.path=NULL, featType="LNC_Filter", annCol=NULL, annColors=list(), seed=63947813, dist0="pearson", link0="ward.D", dist="euclidean", link="ward.D", clstCol=NULL, namecvt=NULL, height=7, method="HC", p=2, dendsort = F)
{
## data format: row:gene; column:sample
# method: c("HC", "KM") # HC: Hirarchical Clustering; KM: K-Means

if (!is.element(el=method, set=c("HC", "KM"))) { stop("method must be HC or KM!") }

if (method=="HC") { methodlabel <- "" }
if (method=="KM") { methodlabel <- "_KM" }

indata=as.data.frame(na.omit(indata))
gene.list=rownames(indata)
sample.list=colnames(indata)
indata=as.matrix(indata)
set.seed(seed)
if(!file.exists(map.res.path)) {dir.create(map.res.path)}
for (k in 1:N.bootstrap) {
    if (k%%1==0) { cat(paste(as.character(k), ",", sep="")) }
        
    outFile=file.path(map.res.path, paste(featType, methodlabel, "_ClusterNum", as.character(N.cluster), "_GeneNum", as.character(N.gene.per.bootstrap), "_Matrix_Permu_", as.character(k), ".txt", sep=""))
    
    if( file.exists(outFile) ) { next }
    
    gene.subset = sample(gene.list, size = N.gene.per.bootstrap, replace = FALSE)
    sample.subset = sample(sample.list, size = N.sample.per.bootstrap, replace = FALSE)
    perm.data=indata[gene.subset,sample.subset]    
    
    if (method=="HC") {
      hc <- hclust(distanceMatrix(perm.data, dist0, p=p), link0)
      hc.ans=cutree(hc, k=N.cluster)
      map.matrix = common.cluster.map (hc.ans, sample.list)  
    }
    
    if (method=="KM") {
      kmeans.ans = kmeans(t(perm.data), N.cluster)$cluster
      map.matrix = common.cluster.map (kmeans.ans, sample.list)
    }
    
    write.table(map.matrix, file=outFile, row.names=F, col.names=F, sep="\t", quote=F)    
}

pattern <- paste(featType, methodlabel, "_ClusterNum", as.character(N.cluster), "_GeneNum", as.character(N.gene.per.bootstrap), "_Matrix_Permu_", ".*.txt$", sep="")


matrix.Files=dir(map.res.path, pattern=pattern)
matrix.Files.N = length(matrix.Files)
matrix.sum = as.data.frame(matrix(0, nrow=length(sample.list), ncol=length(sample.list)))
cat("\nloading matrix file:", fill=T)
for (k in 1:matrix.Files.N) {
    display.progress(k, matrix.Files.N)
    datak = read.table(file.path(map.res.path, matrix.Files[k]),header=F, row.names=NULL, sep="\t", quote="")
    matrix.sum = matrix.sum + datak
}
matrix.sum = matrix.sum*(ncol(indata)/N.sample.per.bootstrap)/matrix.Files.N

row.names(matrix.sum) = names(matrix.sum) <- sample.list
outFigFile <- paste(tumorname,"_ConsensusMap_", featType, methodlabel, "_Cluster", as.character(N.cluster), "_Permu", as.character(N.bootstrap), "_geneNum", as.character(N.gene.per.bootstrap), ".pdf", sep="")  
hc.res <- hclust(distanceMatrix(as.matrix(matrix.sum), dist), link)
if(isTRUE(dendsort)) {
  dendro <- dendsort(as.dendrogram(hc.res))
} else {dendro <- as.dendrogram(hc.res)}

group = cutree(hc.res, k=N.cluster)
tmp <- names(group)
group <- paste("cluster", group, sep="")
if (!is.null(namecvt)) {
  group <- namecvt[group]
}
names(group) <- tmp

#col.lab.color = rainbow(N.cluster)[group[sample.list]]
if(is.null(annCol)) {
  annCol <- as.data.frame(group[colnames(indata)])
  colnames(annCol) <- "Clustering"
  rownames(annCol) <- colnames(indata)
}else{
  annCol$Clustering <- group[rownames(annCol)]
}
if (is.null(clstCol)) {
  clstCol <- rainbow(N.cluster)
}
annColors$Clustering <- clstCol
names(annColors$Clustering) <- unique(annCol$Clustering)

pdf(file.path(fig.path, outFigFile), height=height)
if(sum(hc.res$labels!=colnames(matrix.sum))>0) {stop("colnames mismatch for aheatmap!")}
pheatmap(as.matrix(matrix.sum),
         cluster_rows = as.hclust(dendro),
         cluster_cols = as.hclust(dendro),
         border_color = NA,
         annotation_col = annCol[colnames(indata),,drop = F],
         annotation_colors = annColors,
         color = NMF:::ccRamp(c("#2165AA","#F2F1F3","#B91F33"),64),
         show_colnames = T,
         show_rownames = F)
invisible(dev.off())

return(list(sum=matrix.sum, hc=hc.res, dendro=dendro, files=matrix.Files, clusterN=N.cluster, permN=N.bootstrap, geneN=N.gene.per.bootstrap, sampleN=N.sample.per.bootstrap, group=group, annCol=annCol[colnames(matrix.sum),], annColors=annColors, method=method))

}