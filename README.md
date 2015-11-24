# crambled

The crambled tool (so named since it is part of the crambled-seg pipeline) allows one to interactively explore the cellularity, depth of coverage and clonality of tumour samples that have gone for whole-genome sequencing. In particular, it allows the user to assess the competing solutions for such data that can arise from using different packages.

Crambled is a shiny app for use in R, and so if the crambled_app folder is contained within the R working directory, it is run by typing 

library(shiny)
runApp("crambled_app/")

Crambled is described in the manuscript "Crambled: A Shiny application to enable intuitive resolution of conflicting cellularity estimates" by Lynch AG. Also accompanying the application are two example plots that feature in that manuscript and functions to generate plots from the user's own data.
