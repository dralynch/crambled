RunCrambledShinyApp <- function() 
{
# taken from http://www.r-bloggers.com/supplementing-your-r-package-with-a-shiny-app-2/
# by Dean Attali  
  
  #test <- require(shiny)
  #if (!test) {
  #  message("This function requires Shiny")
  #}
  appDir <- system.file("ShinyApp", "crambled_app", package = "Crambled")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `Crambled`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}