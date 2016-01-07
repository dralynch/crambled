\name{CrambledScanCellline}
\alias{CrambledScanCellline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High level function to produce plots for the Crambled Shiny App if looking at cell line data.
}
\description{
This performs the same job as \code{\link{CrambledScan}}, but in the circumstance that cell line data are being considered and thus no matched normal sample is available.
}
\usage{
CrambledScanCellline(celllinebam, title, window = 51, AFthreshold = 0.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celllinebam}{
A bam file generated from sequencing the cell line.
}
  \item{title}{
A title to be used for captioning plots and naming files (e.g. tumour ID).
}
  \item{window}{
The size of the window on which a running median is to be calculated.
}
  \item{AFthreshold}{
The threshold above which we can assume both alleles are present. Since we anticipate the cell line to be pure, and we don't have a normal sample, we will not be able to distinguish somatic LOH from germline homozygosity. We separate the data at this threshold a) so that germline homozygous sites don't contaminate the running median calculations and b) so that the LOH states are not over-represented in the plot. }
  \item{\dots}{
Other arguments to pass to the plotting functions.
}
}
\value{
A plot to load into the shiny app (via \code{CrambledPlot})
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

## No working example as no bam files are provided

## CrambledScan("celllinebam", "title")

}