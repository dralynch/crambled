
library(Rsamtools)
library(FDb.UCSC.snp137common.hg19)

# tested with Rsamtools_1.20.5, FDb.UCSC.snp137common.hg19_1.0.0

CrambledScan<-function(normalbam,tumourbam,title,window=51,redline=F,...){
  
  test<-require(Rsamtools)
  if(!test){message("This function requires Rsamtools")}
  if(test){
    
    ## Get the list of loci to be examined
    mywhich<-CrambledGetSNPs137()
    
    ## Use Rsamtools to pull out the depths and allele fractions at those locations
    fls <- PileupFiles(c(normalbam, tumourbam))
    param <- ApplyPileupsParam(which=mywhich, what="seq")
    res <- applyPileups(fls, CrambledScanInfo, param=param)
    
    ## process the list to allow smoothing
    res<-unlist(res)
    mydepthlist<-res[seq(1,length(res),2)]
    myaflist<-res[seq(2,length(res),2)]
    myaflist[mydepthlist==0]<-0.5    
    nalist<-is.na(mydepthlist)|is.na(myaflist)
    mydepthlist<-mydepthlist[!nalist]
    myaflist<-myaflist[!nalist]
    
    ## now apply a running median
    rmdepth<-runmed(mydepthlist,window)
    rmaf<-runmed(myaflist,window)
    
    CrambledPlot(rmdepth,rmaf,title,redline=redline,...)
  }
}



CrambledScanCellline<-function(celllinebam,title,window=51,...){
  
  test<-require(Rsamtools)
  if(!test){message("This function requires Rsamtools")}
  if(test){
    ## Get the list of loci to be examined
    mywhich<-CrambledGetSNPs137()
    
    ## Use Rsamtools to pull out the depths and allele fractions at those locations
    fls <- PileupFiles(celllinebam)
    param <- ApplyPileupsParam(which=mywhich, what="seq")
    res <- applyPileups(fls, CrambledScanInfoCellline, param=param)
    
    ## process the list to allow smoothing
    res<-unlist(res)
    mydepthlist<-res[seq(1,length(res),2)]
    myaflist<-res[seq(2,length(res),2)]
    
    ## Need to separate out the heterozygous sites from the others
    rmdepth1<-runmed(mydepthlist[which(myaflist<0.1)],21)
    rmdepth2<-runmed(mydepthlist[which(myaflist>=0.1)],21)
    rmaf1<-runmed(myaflist[which(myaflist<0.1)],21)
    rmaf2<-runmed(myaflist[which(myaflist>=0.1)],21)
    
    
    ## In the cell line, LOH events are not distinguishable from 
    ## the far more numerous germline homozygous loci
    ## we down-sample to avoid the homozygous loci dominating
    mypoints<-sample(length(rmdepth1),length(rmdepth2))
    
    CrambledPlot(c(rmdepth1[mypoints],rmdepth2),c(rmaf1[mypoints],rmaf2),title,redline=redline,...)
  }
}


CrambledPlot<-function(depthvec,afvec,title,xlimmin=0,xlimmax=100,redline=F){
  ##############################
  ## Produces a plot that can ##
  ## be loaded into the Shiny ##
  ## app                      ##
  ##############################
  png(width=600,height=400,file=paste(title,"-shiny.png",sep="",collapse=""))
  smoothScatter(depthvec,afvec,transformation = function(x){x},main=title,xlab="Depth",ylab="B allele frequency",ylim=c(0,0.5),xlim=c(xlimmin,xlimmax))
  if(redline){
    afclass<-as.numeric(cut(afvec,seq(-0.005,0.505,0.01)))
    xbounds<-rep(NA,51)
    for(i in 1:51){
      xbounds[i]<-min(c(depthvec[afclass==i],xlimmax))
    }
    lines(smooth(xbounds),seq(0,0.5,0.01),lwd=2,col="red")
  }      
  dev.off()
}

CrambledScanInfo <- function(x)
{
  ##############################
  ## Checks to see if a locus ##
  ## is heterozygous in the   ##
  ## germline sample, and     ##
  ## reports depth and minor  ##
  ## allele fraction in the   ##
  ## tumour if so             ##
  ##############################
  myAF<-NA
  myDepth<-NA
  if(all(dim(x[["seq"]])==c(5,2,1))){
    basecounts<-x[["seq"]][c("A", "C", "G", "T"),,1,drop=F]
    z<-ifelse(sum(basecounts[basecounts[,1,1]>4,1,1])>14,paste(which(basecounts[,1,1]>4),sep="",collapse=""),"")
    if(z %in% c("12","13","14","23","24","34")){
      y<-basecounts[as.numeric(unlist(strsplit(z,""))),2,]
      myDepth<-sum(y)
      myAF<-min(y)/myDepth
    }
  }
  list(myDepth=myDepth,AF=myAF)
}

CrambledScanInfoCellline<-function(x)
{
  ##############################
  ## For each locus, returns  ##
  ## the depth and allele     ##
  ## fraction based on the    ##
  ## two most frequent bases  ##
  ##############################
  myAF<-NA
  myDepth<-NA
  if(all(dim(x[["seq"]])==c(5,1,1))){
    basecounts<-x[["seq"]][c("A", "C", "G", "T"),,1,drop=F]
    y<-sort(basecounts)[3:4]
    myDepth<-sum(y)
    myAF<-min(y)/myDepth
    
  }
  list(myDepth=myDepth,AF=myAF)
}

CrambledGetSNPs137<-function(){
  ###############################
  ## Returns a GRanges object  ##
  ## giving the locations of   ##
  ## the BAM file to examine   ##
  ###############################
  
  ## In this case we choose all 
  ## features with a decent chance
  ## of being heterozygous in the 
  ## germline sample
  test<-require(FDb.UCSC.snp137common.hg19)
  if(!test){message("This function requires FDb.UCSC.snp137common.hg19")}
  if(test){
        
    snp137common <- features(FDb.UCSC.snp137common.hg19)
    usesnps<-which(as.numeric(snp137common@elementMetadata[,1])>0.4)
    
    mysnps<-snp137common[usesnps,]
    myseqs<-as.character(snp137common@seqnames[usesnps])
    myranges<-(snp137common@ranges[usesnps])
    
    mywhich <- GRanges(myseqs, myranges)
    return(mywhich)
  }
}