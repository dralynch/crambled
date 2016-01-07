\name{CrambledGRangesHG19}
\alias{CrambledGRangesHG19}
\docType{data}
\title{
GRanges object defining likely heterozygous sites in the hg19 genome build
}
\description{
177299 rows, each defining a location among the 22 autosomal chromosomes in the hg19 genome build that has a high probability of being heterozygous. For use with the \code{CrambledSuggestSNPsHG19} function.
}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
  ..@ seqnames       :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 54 levels "chr1","chr10",..: 1 15 19 20 23 24 32 34 35 2 ...
  .. .. ..@ lengths        : int [1:22] 13844 14558 12586 12548 11966 11141 10377 9903 7784 9184 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:177299] 768253 808631 840753 850371 852875 853954 854978 856476 862866 864938 ...
  .. .. ..@ width          : int [1:177299] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ NAMES          : chr [1:177299] "rs2977608" "rs11240779" "rs4970382" "rs6679046" ...
  .. .. ..@ elementType    : chr "integer"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3
  .. .. ..@ lengths        : int 177299
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 177299
  .. .. ..@ listData       : Named list()
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomeInfoDb"] with 4 slots
  .. .. ..@ seqnames   : chr [1:54] "chr1" "chr10" "chr11" "chr12" ...
  .. .. ..@ seqlengths : int [1:54] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ is_circular: logi [1:54] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:54] NA NA NA NA ...
  ..@ metadata       : list()
}
\references{
Lynch A. Crambled: A Shiny application to enable intuitive resolution of conflicting cellularity estimates F1000Research 2015, 4:1407 (doi: 10.12688/f1000research.7453.1)

https://github.com/dralynch/crambled
}
\examples{
MyGRangesObject<-CrambledSuggestSNPsHG19()
#appDir <- system.file("Data", package = "Crambled")
#myGrangesObject<-readRDS(file.path(appDir,"CrambledGRangesHG19.rds"))
}
\keyword{datasets}
