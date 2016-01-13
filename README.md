# crambled

The crambled tool (so named since it is part of the crambled-seg pipeline) allows one to interactively explore the cellularity, depth of coverage and clonality of tumour samples that have gone for whole-genome sequencing. In particular, it allows the user to assess the competing solutions for such data that can arise from using different packages.

The original scripts and Shiny App folders are preserved for consistency with version 1 of the manuscript "Crambled: A Shiny application to enable intuitive resolution of conflicting cellularity estimates" by Lynch AG, however since the 7th of January 2016 the primary version of Crambled should be viewed as being the R package.

The Crambled R package contains the Shiny application, functions to generate plots from the user's data, and the example plots that feature in the manuscript. It also features a vignette giving more detailed explanation of how the functions work.

## Installation

Install the latest version from github:

```R
install.packages("devtools")
devtools::install_github("dralynch/crambled/Package")
```

## Usage

The two main functions for generating plots from BAM files are 

```R
CrambledScan() #and
CrambledScanCellline()
```

while if depths and allele fractions have already been catalogued using e.g. GATK then they can be converted to a plot using 

```R
CrambledPlot()
```

The Crambled Shiny application can be launched using

```R
RunCrambledShinyApp()
```

For more details see the vignette:

```R
browseVignettes(package = "Crambled")
```


