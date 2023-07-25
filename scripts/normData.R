normData <- function(countsTable, tailrows, res.path, tumorname) {

  conds <- factor( rep("samples", times=ncol(countsTable)) )
  cds <- newCountDataSet( countsTable[setdiff(rownames(countsTable), tailrows),], conds )
  cds <- estimateSizeFactors( cds )
  countsNorm <- counts( cds, normalized=TRUE )
  countsNorm <- round(countsNorm, digits=1)
  write.table(countsNorm, file=file.path(res.path, paste(tumorname, ".normalized.count.all.features.txt", sep="")), row.names=T, col.names=NA, sep="\t", quote=F)
  
  options(warn=2)
  tmp <- try (cdsBlind <- estimateDispersions( cds, method="blind"), silent=TRUE)
  if (attributes(tmp)$class=="try-error") {
    cdsBlind <- estimateDispersions( cds, method="blind", fitType="local")
    print ("fitType=local for estimateDispersions")
  }
  options(warn=0)
  vsd <- varianceStabilizingTransformation( cdsBlind )
  
  return(list(countsNorm=countsNorm, vsd=vsd))
  
}