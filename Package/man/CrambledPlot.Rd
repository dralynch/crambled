\name{CrambledPlot}
\alias{CrambledPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to produce a plot for the Crambled Shiny App
}
\description{
Given a vector of sequence coverage and allele frequencies for many SNPs, produces a plot to feed into the crambled Shiny App.
}
\usage{
CrambledPlot(depthvec, afvec, title, xlimmin = 0, xlimmax = 100, redline = F, xpower = 1)
}
\arguments{
  \item{depthvec}{
A vector of depths of coverage associated with a list of (presumed) heterozygous SNPs
}
  \item{afvec}{
A vector of minor allele frequencies associated with a list of (presumed) heterozygous SNPs (note ordered to match the depths of coverage)
}
  \item{title}{
A title to be added to the plot and used in the file name
}
  \item{xlimmin}{
The lower bound for plotting on the x-axis. If changed from the default then the settings of the crambled shiny app will need to be changed also.
}
  \item{xlimmax}{
The upper bound for plotting on the x-axis. If changed from the default then the settings of the crambled shiny app will need to be changed also.
}
  \item{redline}{
If \code{TRUE}, then a smoothed line is added indicating the limits of where data points are observed in the scatter plot. Often this edge is difficult to detect in the plot due to the concentration of observations elsewhere, but can be informative.}
  \item{xpower}{
The power in the transformation of the observation density in the smoothScatter function.
}
}
\value{
The output is a PNG file.
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


## Simulate some plausible looking data
SimCell<-0.62
SimSingleDepth<-30

SimDepths<-NULL
SimAFs<-NULL

##0
DepthA<-(1-SimCell)*SimSingleDepth
DepthB<-(1-SimCell)*SimSingleDepth

SimA<-rpois(500,DepthA)
SimB<-rpois(500,DepthB)

SimDepths<-c(SimDepths,SimA+SimB)
SimAFs<-c(SimAFs,pmin(SimA,SimB)/(SimA+SimB))

##A
DepthA<-(1)*SimSingleDepth
DepthB<-(1-SimCell)*SimSingleDepth

SimA<-rpois(500,DepthA)
SimB<-rpois(500,DepthB)

SimDepths<-c(SimDepths,SimA+SimB)
SimAFs<-c(SimAFs,pmin(SimA,SimB)/(SimA+SimB))

##AB
DepthA<-(1)*SimSingleDepth
DepthB<-(1)*SimSingleDepth

SimA<-rpois(2000,DepthA)
SimB<-rpois(2000,DepthB)

SimDepths<-c(SimDepths,SimA+SimB)
SimAFs<-c(SimAFs,pmin(SimA,SimB)/(SimA+SimB))

##AA
DepthA<-(1+SimCell)*SimSingleDepth
DepthB<-(1-SimCell)*SimSingleDepth

SimA<-rpois(2000,DepthA)
SimB<-rpois(2000,DepthB)

SimDepths<-c(SimDepths,SimA+SimB)
SimAFs<-c(SimAFs,pmin(SimA,SimB)/(SimA+SimB))

##AAB
DepthA<-(1+SimCell)*SimSingleDepth
DepthB<-(1)*SimSingleDepth

SimA<-rpois(3000,DepthA)
SimB<-rpois(3000,DepthB)

SimDepths<-c(SimDepths,SimA+SimB)
SimAFs<-c(SimAFs,pmin(SimA,SimB)/(SimA+SimB))

##AAA
DepthA<-(1+2*SimCell)*SimSingleDepth
DepthB<-(1-SimCell)*SimSingleDepth

SimA<-rpois(1000,DepthA)
SimB<-rpois(1000,DepthB)

SimDepths<-c(SimDepths,SimA+SimB)
SimAFs<-c(SimAFs,pmin(SimA,SimB)/(SimA+SimB))

##AABB
DepthA<-(1+SimCell)*SimSingleDepth
DepthB<-(1+SimCell)*SimSingleDepth

SimA<-rpois(1000,DepthA)
SimB<-rpois(1000,DepthB)

SimDepths<-c(SimDepths,SimA+SimB)
SimAFs<-c(SimAFs,pmin(SimA,SimB)/(SimA+SimB))

SimDepths<-runmed(SimDepths,25)
SimAFs<-runmed(SimAFs,25)


## Draw the figure

CrambledPlot(SimDepths,SimAFs,"ExampleCrambledPlot",redline=TRUE)

}