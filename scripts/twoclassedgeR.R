twoclassedgeR <- function(res.path=NULL, Ginfo=NULL, countsTable=NULL, tailrows=NULL, Groupinfo=NULL, features=NULL, featType=NULL, complist=NULL, PASSFlag=NULL, overwt=FALSE) {
  
  #Groupinfo could contain "batch", which will be considered by edgeR design matrix
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  
  options(warn=1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="")
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(featType, "_edgeR_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) {
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    saminfo <- Groupinfo[samples,]
    colnames(saminfo) <- gsub(pattern="group|Sample_group|Sample_Group|Sample_Group2", replacement="group", colnames(saminfo))
    
    group=factor(saminfo$group,levels = c(ctrlname,treatname))    

    if(is.element("batch",colnames(Groupinfo))) {
      colnames(saminfo) <- gsub(pattern="batch|Batch|BatchEffect|batheffect", replacement="batch", colnames(saminfo))
      cat(paste0("You considered ",length(levels(factor(saminfo$batch)))," batches: ",paste0(levels(factor(saminfo$batch)),collapse = " ")))
      batch=factor(saminfo$batch)
      design <- model.matrix(~batch+group)
    } else {design <- model.matrix(~group)}
    rownames(design) <- samples
    
    y <- DGEList(counts=countsTable[setdiff(rownames(countsTable), tailrows),samples],group=saminfo$group)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit)
    ordered_tags <- topTags(lrt, n=100000)
    allDiff=ordered_tags$table
    allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
    diff=allDiff
    
    features <- intersect(features,intersect(rownames(diff),names(PASSFlag[PASSFlag==TRUE])))
    
    diff=diff[features,]
    diff$id <- Ginfo[rownames(diff),"genename"]
    #res <- diff[!duplicated(diff$id),]
    res <- diff[,c("id","logFC","logCPM","LR","PValue","FDR")]
    res <- res[order(res$FDR),]
    write.table(res, file=outfile, row.names=T, col.names=NA, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
  
}