\name{CrambledScanInfo}
\alias{CrambledScanInfo}
\title{
Function used internally by \code{applyPileups} in \code{Rsamtools} when called by \code{CrambledScan}.
}
\description{
For a genomic location, this function determines whether the normal sample is heterozygous and, if so, returns the depth of sequencing and allele frequency observed in the tumour.
}
\usage{
CrambledScanInfo(x)
}
\arguments{
  \item{x}{
Only used within \code{CrambledScan}}
}
\details{
This is separated out into its own function to facilitate changing to a different function if so desired. It has no value outside of the call to applyPileups.
}
\references{
Lynch A. Crambled: A Shiny application to enable intuitive resolution of conflicting cellularity estimates F1000Research 2015, 4:1407 (doi: 10.12688/f1000research.7453.1)

https://github.com/dralynch/crambled
}
\author{
%%  ~~who you} are~~
Andy Lynch
}
\seealso{
\code{\link{CrambledScan}}
}
\examples{

x<-NULL
tmp<-array(c(12,15,0,1,0,22,10,3,0,0),dim=c(5,2,1))
## Depth = 23+10 = 32
## AF = 10/32 = 0.3125

rownames(tmp)<-c("A","C","G","T","N")
x[["seq"]]<-tmp

CrambledScanInfo(x)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
