library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Cellularity Estimation"),
  sidebarLayout(position = "left",
                
                sidebarPanel(fileInput('file1', 'Choose png file \n (defaults to cell line example if no file selected)',
                                       accept=c('.png')),
                             sliderInput("mycell",
                                          "Cellularity:",
                                          min = 0.01,
                                          max = 0.99,
                                          value = .6),
                             sliderInput("mydepth",
                                         "Depth of coverage for one copy:",
                                         min = 1,
                                         max = 70,
                                         step=.2,
                                         value = 30),
                checkboxInput("MD", "MetaData"),
                conditionalPanel(
                  condition = "input.MD == true",
                  numericInput("mywidth", "Image width",600, min = 0),
                  numericInput("myheight", "Image height",400, min = 0),
                  numericInput("myxlim1", "X lower limit",0, min = 0),
                  numericInput("myxlim2", "X upper limit",100, min = 0)                  
                )),
                
         mainPanel(
           imageOutput("myImage"),
           plotOutput("predplot")
         )
  )
))