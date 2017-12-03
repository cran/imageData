#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test using Rice baby example\n")
test_that("Rice2015_imageData", {
  skip_if_not_installed("imageData")
  skip_on_cran()
  library(dae)
  library(ggplot2)
  library(imageData)
  library(knitr)
  opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
                 "fig.width" = 9, "fig.height" = 9)
  
  
  ## A dummy example to illustrate the use of imageData
  #'# Step 1: Import the raw data
  data(RiceRaw.dat)
  raw.dat <- RiceRaw.dat[1:280, ]
  raw.dat$Smarthouse <- 1
  
  #'# Step 2: Select imaging variables and add covariates and factors (produces longi.dat)
  longi.prime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1)
  
  longi.dat <- designFactors(longi.prime.dat, insertName = "xDays",
                             nzones = 1, nlanesperzone = 1, nmainplotsperlane = 10, 
                             designfactorMethod="StandardOrder")
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 52)
  
  
  #'## Plot the imaging times for 20 carts from a single Lane
  testthat::expect_silent(imagetimesPlot(longi.dat, intervals = "Days", 
                                         timePositions = "Hour",
                                         ggplotFuncs=list(ggplot2::geom_line(aes(group=Snapshot.ID.Tag, 
                                                                                 colour=Lane)))))
  
  #'## Particular edits to longi.dat
  longi.dat <- within(longi.dat, 
                      { 
                        Days.after.Salting <- as.numfac(Days) - 29
                      })
  
  #'# Step 3: Form derived traits that result in a value for each observation
  #'### Set responses
  responses.image <- c("Area")
  responses.smooth <- paste(responses.image, "smooth", sep=".")
  
  #'## Form growth rates for each observation of a subset of responses by differencing
  longi.dat <- splitContGRdiff(longi.dat, responses.image, 
                               INDICES="Snapshot.ID.Tag",
                               which.rates = c("AGR","RGR"))
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 56)
  
  #'## Form Area.WUE 
  longi.dat <- within(longi.dat, 
                      { 
                        Area.WUE <- WUI(Area.AGR*Days.diffs, Water.Loss)
                      })
  
  #'## Add cumulative responses 
  longi.dat <- within(longi.dat, 
                      { 
                        Water.Loss.Cum <- unlist(by(Water.Loss, Snapshot.ID.Tag, 
                                                    cumulate, exclude.1st=TRUE))
                        WUE.cum <- Area / Water.Loss.Cum 
                      })
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 59)
  
  #'# Step 4: Fit splines to smooth the longitudinal trends in the primary traits and calculate their growth rates
  #'
  #'## Smooth responses
  #+
  for (response in c(responses.image, "Water.Loss"))
    longi.dat <- splitSplines(longi.dat, response, x="xDays", INDICES = "Snapshot.ID.Tag", 
                              df = 4)
  longi.dat <- with(longi.dat, longi.dat[order(Snapshot.ID.Tag, xDays), ])
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 61)
  
  #'## Loop over smoothed responses, forming growth rates by differences
  #+
  responses.GR <- paste(responses.smooth, "AGR", sep=".")
  longi.dat <- splitContGRdiff(longi.dat, responses.smooth, 
                               INDICES="Snapshot.ID.Tag",
                               which.rates = c("AGR","RGR"))
  testthat::expect_equal(nrow(longi.dat), 280)
  testthat::expect_equal(ncol(longi.dat), 63)
  
  #'## Finalize longi.dat
  longi.dat <- with(longi.dat, longi.dat[order(Snapshot.ID.Tag, xDays), ])
  
  #'# Step 5: Do exploratory plots on unsmoothed and smoothed longitudinal data
  responses.longi <- c("Area","Area.AGR","Area.RGR", "Area.WUE")
  responses.smooth.plot <- c("Area.smooth","Area.smooth.AGR","Area.smooth.RGR")
  titles <- c("Total area (1000 pixels)", 
              "Total area AGR (1000 pixels per day)", "Total area RGR (per day)",
              "Total area WUE (1000 pixels per mL)")
  titles.smooth<-titles
  nresp <- length(responses.longi)
  limits <- list(c(0,1000), c(-50,125), c(-0.05,0.40), c(0,30))
  
  #' ### Plot unsmoothed profiles for all longitudinal  responses 
  #+ "01-ProfilesAll"
  klimit <- 0
  for (k in 1:nresp)
  { 
    klimit <- klimit + 1
    plt <- longiPlot(data = longi.dat, response = responses.longi[k], 
                     y.title = titles[k], x="xDays+35.42857143", printPlot=FALSE)
    plt <- plt + geom_vline(xintercept=29, linetype="longdash", size=1) +
      scale_x_continuous(breaks=seq(28, 42, by=2)) + 
      scale_y_continuous(limits=limits[[klimit]])
    if (k < 4)
      testthat::expect_silent(print(plt))
    else
      testthat::expect_warning(print(plt))
  }
  
  
  #' ### Plot smoothed profiles for all longitudinal responses - GRs by difference
  #+ "01-SmoothedProfilesAll"
  nresp.smooth <- length(responses.smooth.plot)
  limits <- list(c(0,1000), c(0,100), c(0.0,0.40))
  for (k in 1:nresp.smooth)
  { 
    plt <- longiPlot(data = longi.dat, response = responses.smooth.plot[k], 
                     y.title = titles.smooth[k], x="xDays+35.42857143", printPlot=FALSE)
    plt <- plt + geom_vline(xintercept=29, linetype="longdash", size=1) +
      scale_x_continuous(breaks=seq(28, 42, by=2)) + 
      scale_y_continuous(limits=limits[[k]])
    testthat::expect_silent(print(plt))
  }
  
  
  #'### AGR anomalies - plot without anomalous plants followed by plot of anomalous plants
  #+ "01-0254-AGRanomalies"
  anom.ID <- vector(mode = "character", length = 0L)
  response <- "Area.smooth.AGR"
  cols.output <- c("Snapshot.ID.Tag", "Smarthouse", "Lane", "Position", 
                   "Treatment.1", "Genotype.ID", "Days")
  anomalous <- anomPlot(longi.dat, response=response, lower=2.5, start.time=40, 
                        x = "xDays+35.42857143", vertical.line=29, breaks=seq(28, 42, by=2), 
                        whichPrint=c("innerPlot"), y.title=response)
  subs <- subset(anomalous$data, Area.smooth.AGR.anom & Days==42)
  if (nrow(subs) == 0)
  { cat("\n#### No anomalous data here\n\n")
  } else
  { 
    subs <- subs[order(subs["Smarthouse"],subs["Treatment.1"], subs[response]),]
    print(subs[c(cols.output, response)])
    anom.ID <- unique(c(anom.ID, subs$Snapshot.ID.Tag))
    outerPlot <- anomalous$outerPlot  + geom_text(data=subs,
                                                  aes_string(x = "xDays+35.42857143", 
                                                             y = response, 
                                                             label="Snapshot.ID.Tag"), 
                                                  size=3, hjust=0.7, vjust=0.5)
    print(outerPlot)
  }
  testthat::expect_equal(nrow(subs), 0)  
  
  #'# Step 6: Form single-value plant responses in Snapshot.ID.Tag order.
  #'
  #'## 6a) Set up a data frame with factors only
  #+
  cart.dat <- longi.dat[longi.dat$Days == 31, 
                        c("Smarthouse","Lane","Position","Snapshot.ID.Tag",
                          "xPosn","xMainPosn",
                          "Zones","xZones","SHZones","ZLane","ZMainplots", "Subplots",
                          "Genotype.ID","Treatment.1")]
  cart.dat <- cart.dat[do.call(order, cart.dat), ]
  
  #'## 6b) Get responses based on first and last date.
  #'
  #'### Observation for first and last date
  cart.dat <- cbind(cart.dat, getDates(responses.image, data = longi.dat, 
                                       which.times = c(31), suffix = "first"))
  cart.dat <- cbind(cart.dat, getDates(responses.image, data = longi.dat, 
                                       which.times = c(42), suffix = "last"))
  cart.dat <- cbind(cart.dat, getDates(c("WUE.cum"), 
                                       data = longi.dat, 
                                       which.times = c(42), suffix = "last"))
  responses.smooth <- paste(responses.image, "smooth", sep=".")
  cart.dat <- cbind(cart.dat, getDates(responses.smooth, data = longi.dat, 
                                       which.times = c(31), suffix = "first"))
  cart.dat <- cbind(cart.dat, getDates(responses.smooth, data = longi.dat, 
                                       which.times = c(42), suffix = "last"))
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 19)
  
  #'### Growth rates over whole period.
  #+
  tottime <- 42 - 31
  cart.dat <- within(cart.dat, 
                     { 
                       Area.AGR <- (Area.last - Area.first)/tottime
                       Area.RGR <- log(Area.last / Area.first)/tottime
                     })
  
  #'### Calculate water index over whole period
  cart.dat <- merge(cart.dat, 
                    intervalWUI("Area", water.use = "Water.Loss", 
                                start.times = c(31), 
                                end.times = c(42), 
                                suffix = NULL, 
                                data = longi.dat, include.total.water = TRUE),
                    by = c("Snapshot.ID.Tag"))
  names(cart.dat)[match(c("Area.WUI","Water.Loss.Total"),names(cart.dat))] <- c("Area.Overall.WUE", 
                                                                                "Water.Loss.Overall")
  cart.dat$Water.Loss.rate.Overall <- cart.dat$Water.Loss.Overall / (42 - 31)
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 25)
  
  #'## 6c) Add growth rates and water indices for intervals
  #'### Set up intervals
  #+
  start.days <- list(31,35,31,38)
  end.days <- list(35,38,38,42)
  suffices <- list("31to35","35to38","31to38","38to42")
  
  #'### Rates for specific intervals from the smoothed data by differencing
  #+
  for (r in responses.smooth)
  { for (k in 1:length(suffices))
  { 
    cart.dat <- merge(cart.dat, 
                      intervalGRdiff(r, 
                                     which.rates = c("AGR","RGR"), 
                                     start.times = start.days[k][[1]], 
                                     end.times = end.days[k][[1]], 
                                     suffix.interval = suffices[k][[1]], 
                                     data = longi.dat),
                      by = "Snapshot.ID.Tag")
  }
  }
  
  #'### Water indices for specific intervals from the unsmoothed and smoothed data
  #+
  for (k in 1:length(suffices))
  { 
    cart.dat <- merge(cart.dat, 
                      intervalWUI("Area", water.use = "Water.Loss", 
                                  start.times = start.days[k][[1]], 
                                  end.times = end.days[k][[1]], 
                                  suffix = suffices[k][[1]], 
                                  data = longi.dat, include.total.water = TRUE),
                      by = "Snapshot.ID.Tag")
    names(cart.dat)[match(paste("Area.WUI", suffices[k][[1]], sep="."), 
                          names(cart.dat))] <- paste("Area.WUE", suffices[k][[1]], sep=".")
    cart.dat[paste("Water.Loss.rate", suffices[k][[1]], sep=".")] <- 
      cart.dat[[paste("Water.Loss.Total", suffices[k][[1]], sep=".")]] / 
      ( end.days[k][[1]] - start.days[k][[1]])
  }
  
  cart.dat <- with(cart.dat, cart.dat[order(Snapshot.ID.Tag), ])
  testthat::expect_equal(nrow(cart.dat), 20)
  testthat::expect_equal(ncol(cart.dat), 49)
  
  #'# Step 7: Form continuous and interval SIITs
  #'
  #'## 7a) Calculate continuous
  #+
  cols.retained <-  c("Snapshot.ID.Tag","Smarthouse","Lane","Position",
                      "Days","Snapshot.Time.Stamp", "Hour", "xDays",
                      "Zones","xZones","SHZones","ZLane","ZMainplots",
                      "xMainPosn", "Genotype.ID")
  responses.GR <- c("Area.smooth.AGR","Area.smooth.AGR","Area.smooth.RGR")
  suffices.results <- c("diff", "SIIT", "SIIT")
  responses.SIIT <- unlist(Map(paste, responses.GR, suffices.results,sep="."))
  
  longi.SIIT.dat <- 
    twoLevelOpcreate(responses.GR, longi.dat, suffices.treatment=c("C","S"),
                     operations = c("-", "/", "/"), suffices.results = suffices.results, 
                     columns.retained = cols.retained, 
                     by = c("Smarthouse","Zones","ZMainplots","Days"))
  longi.SIIT.dat <- with(longi.SIIT.dat, 
                         longi.SIIT.dat[order(Smarthouse,Zones,ZMainplots,Days),])
  testthat::expect_equal(nrow(longi.SIIT.dat), 140)
  testthat::expect_equal(ncol(longi.SIIT.dat), 22)
  
  #' ### Plot SIIT profiles 
  #' 
  #+ "03-SIITProfiles"
  k <- 2
  nresp <- length(responses.SIIT)
  limits <- with(longi.SIIT.dat, list(c(min(Area.smooth.AGR.diff, na.rm=TRUE),
                                        max(Area.smooth.AGR.diff, na.rm=TRUE)),
                                      c(0,3),
                                      c(0,1.5)))
  #Plots
  for (k in 1:nresp)
  { 
    plt <- longiPlot(data = longi.SIIT.dat, x="xDays+35.42857143", 
                     response = responses.SIIT[k], 
                     y.title=responses.SIIT[k], 
                     facet.x="Smarthouse", facet.y=".", printPlot=FALSE, )
    plt <- plt + geom_vline(xintercept=29, linetype="longdash", size=1) +
      scale_x_continuous(breaks=seq(28, 42, by=2)) + 
      scale_y_continuous(limits=limits[[k]])
    testthat::expect_silent(print(plt))
  }
  
  #'## 7b) Calculate interval SIITs and check for large values for SIIT for Days 31to35
  #+ "01-SIITIntClean"
  suffices <- list("31to35","35to38","31to38","38to42")
  response <- "Area.smooth.RGR.31to35"
  SIIT <- paste(response, "SIIT", sep=".")
  responses.SIITinterval <- as.vector(outer("Area.smooth.RGR", suffices, paste, sep="."))
  
  cart.SIIT.dat <- twoLevelOpcreate(responses.SIITinterval, cart.dat,
                                    suffices.treatment=c("C","S"), 
                                    suffices.results="SIIT", 
                                    columns.suffixed="Snapshot.ID.Tag")
  testthat::expect_equal(nrow(cart.SIIT.dat), 10)
  testthat::expect_equal(ncol(cart.SIIT.dat), 23)
  tmp<-na.omit(cart.SIIT.dat)
  print(summary(tmp[SIIT]))
  big.SIIT <- with(tmp, tmp[tmp[SIIT] > 1.15, c("Snapshot.ID.Tag.C","Genotype.ID",
                                                paste(response,"C",sep="."), 
                                                paste(response,"S",sep="."), SIIT)])
  big.SIIT <- big.SIIT[order(big.SIIT[SIIT]),]
  testthat::expect_equal(nrow(big.SIIT), 0)
  testthat::expect_equal(ncol(big.SIIT), 5)
  print(big.SIIT)
  plt <- ggplot(tmp, aes_string(SIIT)) +
    geom_histogram(aes(y = ..density..), binwidth=0.05) +
    geom_vline(xintercept=1.15, linetype="longdash", size=1) +
    theme_bw() + facet_grid(Smarthouse ~.)
  print(plt)
  plt <- ggplot(tmp, aes_string(x="Smarthouse", y=SIIT)) +
    geom_boxplot() + theme_bw()
  testthat::expect_silent(print(plt))

})
