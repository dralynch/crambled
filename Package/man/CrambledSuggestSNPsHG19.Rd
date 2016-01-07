\name{CrambledSuggestSNPsHG19}
\alias{CrambledSuggestSNPsHG19}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reads in a file containing approximately 175,000 sites that are commonly heterozygous in the hg19 genome.
}
\description{
Returns a GRanges object for use by \code{CrambledScan} etc.}
\usage{
CrambledSuggestSNPsHG19()
}
\value{
Returns a Granges object
}
\references{
Lynch A. Crambled: A Shiny application to enable intuitive resolution of conflicting cellularity estimates F1000Research 2015, 4:1407 (doi: 10.12688/f1000research.7453.1)

https://github.com/dralynch/crambled
}
\author{
Andy Lynch
}
\seealso{
\code{\link{CrambledScan}}
}
\examples{
MyGRangesObject<-CrambledSuggestSNPsHG19()
}

