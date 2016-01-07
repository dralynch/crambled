\name{CrambledScan}
\alias{CrambledScan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High level function to produce plots for the Crambled Shiny App.
}
\description{
This function obtains a list of sites likely to be heterozygous in the germline (normal) sample being considered, reduces this to a list that do seem to be heterozygous and then interrogates the tumour sample at these locations to obtain details of sequencing depth and allele frequency. The values are then smoothed using a runnig median before the data are plotted.
}
\usage{
CrambledScan(normalbam, tumourbam, title, window = 51, redline = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{normalbam}{
A bam file that should be representative of the germline genotypes associated with the tumour.
}
  \item{tumourbam}{
A bam file generated from sequencing the tumour.
}
  \item{title}{
A title to be used for captioning plots and naming files (e.g. tumour ID).
}
  \item{window}{
The size of the window on which a running median is to be calculated.
}
  \item{redline}{
Whether a red boundary line should be added by the plotting function}
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
\code{\link{CrambledScanCellline}}
}
\examples{

## No working example as no bam files are provided

## CrambledScan("normalbam", "tumourbam", "title")

}