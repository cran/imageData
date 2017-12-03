cat("#### Test corrPlot\n")
test_that("corrPlot", {
  skip_if_not_installed("imageData")
  skip_on_cran()
  library(imageData)
  library(ggplot2)
  library(GGally)
  
  data(exampleData)
  responses <- c("Area","Area.SV","Area.TV", "Image.Biomass", "Max.Height","Centre.Mass",
                 "Density", "Compactness.TV", "Compactness.SV")
 testthat::expect_silent(corrPlot(responses, longi.dat, pairs.sets=list(c(1:4),c(5:7))))
})