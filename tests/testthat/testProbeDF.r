
cat("#### Test probeDF with Judith Atieno 0278\n")
test_that("chickpea_imageData", {
  skip_if_not_installed("imageData")
  skip_on_cran()
  library(imageData)
  library(ggplot2)
  
  data(dat1)
  
  #'## Values of df for which to obtain plots
  df <- c(4,7)
  
  
  #'## Obtain plots
  t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", response="ShootArea1000",
               which.plots=c("bothseparately"), 
               ggplotFuncs=list(scale_x_continuous(breaks=seq(30, 50, by=2))))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_equal(ncol(t), 15)
  
  t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", response="ShootArea1000",
               which.traits=c("AGR"))
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_equal(ncol(t), 13)
  
  t <- probeDF(dat1, df=df, xname="TimeAfterPlanting", response="ShootArea1000",
               which.traits=c("response"), which.plots="compare", get.rates=FALSE)
  testthat::expect_equal(nrow(t), 8480)
  testthat::expect_equal(ncol(t), 8)
  
})

