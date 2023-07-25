LowExpFilter <- function(FPKM, countsNorm, lowcut=10, lowpct=NULL, Mids, Lids, Ginfo, res.path, tumorname, lownum=NULL, quantile=NULL, overwrite=TRUE) {

  # overwrite: overwrite files?
  
  nullsum <- sum( c(is.null(quantile), is.null(lowpct), is.null(lownum)) )
  if (nullsum!=2) { stop("quantile, lowpct, lownum can only be set once!") }
  
  if (!is.null(lownum)) {
    PASSFlag <- rowSums(countsNorm>=lowcut) >= lownum
    names(PASSFlag) <- rownames(countsNorm)
  }
  if (!is.null(lowpct)) {
    PASSFlag <- rowSums(FPKM<lowcut) < ceiling(lowpct *ncol(FPKM))
    names(PASSFlag) <- rownames(FPKM)
  }
  if (!is.null(quantile)) {
    rowmax <- apply(FPKM, 1, max, na.rm=T)
    rowmax <- rowmax[rowmax>0]
    
    Lrowmax <- rowmax[intersect(names(rowmax), Lids)]
    Mrowmax <- rowmax[intersect(names(rowmax), Mids)]
    
    Lqtl <- quantile(Lrowmax, probs=quantile, na.rm=TRUE)
    Mqtl <- quantile(Mrowmax, probs=quantile, na.rm=TRUE)
    
    LFPKM <- FPKM[intersect(rownames(FPKM), Lids), ]
    MFPKM <- FPKM[intersect(rownames(FPKM), Mids), ]
    
    LPASSFlag <- rowSums(LFPKM>=Lqtl) >= 3
    MPASSFlag <- rowSums(MFPKM>=Mqtl) >= 3
    
    PASSFlag <- c(LPASSFlag, MPASSFlag)
    if (length(PASSFlag) != nrow(FPKM)) {stop("row numbers of PASSFlag and FPKM differ!")}
    PASSFlag <- PASSFlag[rownames(FPKM)]
  }
  
  if (overwrite) {
    ### write FPKM data ###
    outdata <- FPKM[intersect(Mids, names(PASSFlag[PASSFlag==TRUE])),]
    outdata <- outdata[!duplicated(Ginfo[rownames(outdata), "genename"]),]
    outdata <- outdata[!is.na(Ginfo[rownames(outdata), "genename"]), ]
    rownames(outdata) <- Ginfo[rownames(outdata), "genename"]
    write.table(outdata, file=file.path(res.path, paste(tumorname, ".FPKM.mRNA.LowExpfiltered.txt", sep="")), row.names=T, col.names=NA, sep="\t", quote=F)
    
    outdata <- FPKM[intersect(Lids, names(PASSFlag[PASSFlag==TRUE])),]
    outdata <- outdata[!duplicated(Ginfo[rownames(outdata), "genename"]),]
    outdata <- outdata[!is.na(Ginfo[rownames(outdata), "genename"]), ]
    rownames(outdata) <- Ginfo[rownames(outdata), "genename"]
    write.table(outdata, file=file.path(res.path, paste(tumorname, ".FPKM.lncRNA.LowExpfiltered.txt", sep="")), row.names=T, col.names=NA, sep="\t", quote=F)
    
    ### write count data ###
    outdata <- countsNorm[intersect(Mids, names(PASSFlag[PASSFlag==TRUE])),]
    outdata <- outdata[!duplicated(Ginfo[rownames(outdata), "genename"]),]
    outdata <- outdata[!is.na(Ginfo[rownames(outdata), "genename"]), ]  
    rownames(outdata) <- Ginfo[rownames(outdata), "genename"]
    write.table(outdata, file=file.path(res.path, paste(tumorname, ".Count.mRNA.LowExpfiltered.txt", sep="")), row.names=T, col.names=NA, sep="\t", quote=F)
    
    outdata <- countsNorm[intersect(Lids, names(PASSFlag[PASSFlag==TRUE])),]
    outdata <- outdata[!duplicated(Ginfo[rownames(outdata), "genename"]),]
    outdata <- outdata[!is.na(Ginfo[rownames(outdata), "genename"]), ]  
    rownames(outdata) <- Ginfo[rownames(outdata), "genename"]
    write.table(outdata, file=file.path(res.path, paste(tumorname, ".Count.lncRNA.LowExpfiltered.txt", sep="")), row.names=T, col.names=NA, sep="\t", quote=F)
  }
  return (PASSFlag)
}