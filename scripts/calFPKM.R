calFPKM <- function(countsTable, Vids, tailrows, Ginfo, outfile) {
  
  TotalRdNum <- colSums(countsTable[setdiff(rownames(countsTable), Vids),])
  indata <- countsTable[setdiff(rownames(countsTable), c(Vids, tailrows)),]
  
  totalmtx <- (matrix(rep(1,times=nrow(indata)), ncol=1) %*% matrix(TotalRdNum, nrow=1)) / 1000000
  lenmtx <- (matrix(Ginfo[rownames(indata),"unqlen"], ncol=1) %*% matrix(rep(1,times=ncol(indata)), nrow=1)) / 1000
  FPKM <- round (indata / totalmtx / lenmtx, digits=2)
  write.table(FPKM, file=outfile, row.names=T, col.names=NA, sep="\t", quote=F)
  
  return (FPKM)
  
}
