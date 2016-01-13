# crambled

The crambled tool (so named since it is part of the crambled-seg pipeline) allows one to interactively explore the cellularity, depth of coverage and clonality of tumour samples that have gone for whole-genome sequencing. In particular, it allows the user to assess the competing solutions for such data that can arise from using different packages.

Crambled is a shiny app for use in R, and so if the crambled_app folder is contained within the R working directory, it is run by typing 

library(shiny)
runApp("crambled_app/")

Crambled is described in the manuscript "Crambled: A Shiny application to enable intuitive resolution of conflicting cellularity estimates" by Lynch AG. Also accompanying the application are two example plots that feature in that manuscript and functions to generate plots from the user's own data.


As of 7th January 2016, for the convenience of the user there is an R package incorporating the Shiny application and associated scripts. 

## Installation

Install the latest version from github:

```R
install.packages("devtools")
devtools::install_github("dralynch/crambled/Package")
```