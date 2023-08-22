\newcommand{\packageVer}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Version")}}
\newcommand{\packageDate}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Date")}}
\name{imageData-package}
\alias{imageData-package}
\alias{imageData}
\docType{package}
\title{
\packageTitle{imageData}
}
\description{\packageDescription{imageData}

\bold{Version:} {\packageVer{imageData}}

\bold{Date:} \packageDate{imageData}
}


\section{Index}{
For an overview of the use of these functions and an example see below.
\tabular{ll}{
(i) Data\cr\cr

\code{\link{RiceRaw.dat}} 
  \tab Data for an experiment to investigate a rice \cr
       germplasm panel.\cr\cr

(ii) Data frame manipulation\cr\cr

\code{\link{designFactors}} 
  \tab Adds the factors and covariates for a blocked,\cr
  \tab split-plot design.\cr
\code{\link{getDates}} 
  \tab Forms a subset of 'responses' in 'data' that\cr
  \tab contains their values for the nominated times.\cr
\code{\link{importExcel}} 
  \tab Imports an Excel imaging file and allows some\cr
  \tab renaming of variables.\cr
\code{\link{longitudinalPrime}} 
  \tab Selects a set variables to be retained in a\cr
  \tab data frame of longitudinal data.\cr
\code{\link{twoLevelOpcreate}} 
  \tab Creates a data.frame formed by applying, for\cr
  \tab each response, abinary operation to the values of \cr
  \tab two different treatments.\cr
\cr
(iii) Plots\cr\cr

\code{\link{anomPlot}} 
  \tab Identifies anomalous individuals and produces\cr
  \tab longitudinal plots without them and with just them.\cr
\code{\link{corrPlot}} 
  \tab Calculates and plots correlation matrices for a\cr
  \tab set of responses.\cr
\code{\link{imagetimesPlot}} 
  \tab Plots the time within an interval versus the interval.\cr 
  \tab For example, the hour of the day carts are imaged\cr
  \tab against the days after planting (or some other\cr
  \tab number of days after an event).\cr
\code{\link{longiPlot}} 
  \tab Plots longitudinal data from a Lemna Tec \cr
  \tab Scananalyzer.\cr
\code{\link{probeDF}} 
  \tab Compares, for a set of specified values of df,\cr
  \tab a response and the smooths of it, possibly along\cr
  \tab  with growth rates calculated from the smooths.\cr
\cr
(iv) Calculations value-by-value\cr\cr

\code{\link{GrowthRates}} 
  \tab Calculates growth rates (AGR, PGR, RGRdiff) \cr
  \tab between pairs of values in a vector.\cr
\code{\link{WUI}} 
  \tab Calculates the Water Use Index (WUI).\cr
\code{\link{anom}} 
  \tab Tests if any values in a vector are anomalous\cr
  \tab in being outside specified limits.\cr
\code{\link{calcTimes}} 
  \tab Calculates for a set of times, the time intervals \cr
  \tab after an origin time and the position of each with \cr
  \tab in that time.\cr
\code{\link{calcLagged}} 
  \tab Replaces the values in a vector with the result\cr
  \tab of applying an operation to it and a lagged value.\cr
\code{\link{cumulate}} 
  \tab Calculates the cumulative sum, ignoring the\cr
  \tab first element if exclude.1st is TRUE.\cr
\cr
(v) Calculations over multiple values\cr\cr

\code{\link{fitSpline}} 
  \tab Produce the fits from a natural cubic smoothing\cr
  \tab spline applied to a response in a 'data.frame'.\cr
\code{\link{intervalGRaverage}} 
  \tab Calculates the growth rates for a specified\cr
  \tab time interval by taking weighted averages of\cr
  \tab growth rates for times within the interval.\cr
\code{\link{intervalGRdiff}} 
  \tab Calculates the growth rates for a specified\cr
  \tab time interval.\cr
\code{\link{intervalValueCalculate}}
  \tab Calculates a single value that is a function of\cr
  \tab an individual's values for a response over a\cr
  \tab specified time interval.\cr
\code{\link{intervalWUI}} 
  \tab Calculates water use indices (WUI) over a\cr
  \tab specified time interval to a data.frame.\cr
\cr
(vi) Caclulations in each split of a 'data.frame' \cr
\cr
\code{\link{splitContGRdiff}} 
  \tab Adds the growth rates calculated continuously\cr
  \tab over time for subsets of a response to a\cr
  \tab 'data.frame'.\cr
\code{\link{splitSplines}} 
  \tab Adds the fits after fitting a natural cubic\cr
  \tab smoothing spline to subsets of a response to a\cr
  \tab 'data.frame'.\cr
\code{\link{splitValueCalculate}}   
  \tab Calculates a single value that is a function of\cr
  \tab an individual's values for a response.\cr
\cr
(vii) Principal variates analysis (PV A)\cr\cr

\code{\link{intervalPVA}} 
  \tab Selects a subset of variables observed within a\cr
  \tab specified time interval using PVA.\cr
\code{\link{PVA}} 
  \tab Selects a subset of variables using PVA.\cr
\code{\link{rcontrib}} 
  \tab Computes a measure of how correlated each\cr
  \tab variable in a set is with the other variable,\cr
  \tab conditional on a nominated subset of them.\cr
}}
\author{
\packageAuthor{imageData}

Maintainer: \packageMaintainer{imageData}
}
\keyword{ hplot }
\keyword{ manip }
\keyword{ package }
\seealso{\code{\link[dae:dae-package]{dae}}
}
\section{Overview}{
This package can be used to carry out a full seven-step process to produce phenotypic traits from measurements made in a high-throughput phenotyping facility, such as one based on a Lemna-Tec Scananalyzer 3D system and described by Al-Tamimi et al. (2016). Otherwise, individual functions can be used to carry out parts of the process.

The basic data consists of imaging data obtained from a set of pots or carts over time. The carts are arranged in a grid of Lanes \eqn{\times}{x} Positions. There should be a unique identifier for each cart, which by default is \code{Snapshot.ID.Tag}, and variable giving the Days after Planting for each measurement, by default {\code{Time.after.Planting..d.}}. In some cases, it is expected that  there will be a column labelled \code{Snapshot.Time.Stamp}, which reflects the time of the imaging from which a particular data value was obtained.

The full seven-step process is as follows:

\enumerate{
 \item Use \code{\link{importExcel}} to import the raw data from the Excel file. This step should also involve any editing of the data needed to take account of mishaps during the data collection and the need to remove faulty data (produces \code{raw.dat}). Generally, data can be removed by replacing only values for responses with  missing values (\code{NA}) for carts whose data is to be removed, leaving the identifying information intact.
 \item Use \code{\link{longitudinalPrime}} to select a subset of the imaging variables produced by the Lemna Tec Scanalyzer and, if the design is a blocked, split-plot design, use \code{designFactors} to add covariates and factors that might be used in the analysis (produces the data frame \code{longi.prime.dat}).
 \item Add derived traits that result in a value for each observation: use \code{\link{splitContGRdiff}} to obtain continuous growth rates i.e. a growth rate for each time of observation, except the first; \code{\link{WUI}} to produce continuous Water Use Efficiency Indices (WUE) and \code{cumulate} to produce cumulative responses. (Produces the data frame \code{longi.dat}.)
 \item Use \code{\link{splitSplines}} to fit splines to smooth the longitudinal trends in the primary traits and calculate continuous growth rates from the smoothed response (added to the data frame \code{longi.dat}). There are two options for calculating continuous smoothed growth rates: (i) by differencing --- use  \code{\link{splitContGRdiff}}; (ii) from the first derivatives of the splines --- in \code{\link{splitSplines}} include \code{1} in the \code{deriv} argument, include \code{"AGR"} in \code{suffices.deriv} and set the \code{RGR} to say \code{"RGR"}. Optionally, use \code{\link{probeDF}} to compare the smooths for a number of values of \code{df} and, if necessary, re-run \code{\link{splitSplines}} with a revised value of \code{df}.
 \item Perform an exploratory examination of the unsmoothed data by using \code{longiPlot} to produce longitudinal plots of unsmoothed imaging traits and continuous growth rates. Also, use \code{\link{longiPlot}} to plot the smoothed imaging traits and continuous growth rates and \code{\link{anomPlot}} to check for anomalies in the data.
 \item Produce cart data: traits for which there is a single value for each \code{Snapshot.ID.Tag} or cart. (produces the data frame \code{cart.dat})
 \enumerate{
   \item Set up a cart data.frame with the factors and covariates for a single observation from all carts. This can be done by subsetting \code{longi.dat} so that there is one entry for each cart.
   \item Use \code{\link{getDates}} to add traits at specific times to the cart \code{data.frame}, often the first and last day of imaging for each \code{Snapshot.ID.Tag}. The times need to be selected so that there is one and only one observation for each cart. Also form traits, such as growth rates over the whole imaging period, based on these values
   \item Based on the longitudinal plots, decide on the intervals for  which growth rates and WUEs are to be calculated. The growth rates for intervals are calculated from the continuous growth rates, using \code{\link{intervalGRdiff}}, if the continuous growth rates were calculated by differencing, or \code{\link{intervalGRaverage}}, if they were calculated from first derivatives. To calculate WUEs for intervals, use \code{\link{intervalWUI}}, The interval growth rates and WUEs are added to the cart \code{data.frame}.
 }
   \item (Optional) There is also the possibility that, for experiments investigating salinity, the Shoot Ion Independent Tolerance (SIIT) index can be calculated using \code{\link{twoLevelOpcreate}}. 
}
}

\examples{\dontrun{
### This example can be run because the data.frame RiceRaw.dat is available with the package
#'# Step 1: Import the raw data
data(RiceRaw.dat)

#'# Step 2: Select imaging variables and add covariates and factors (produces longi.dat)
longi.dat <- longitudinalPrime(data=RiceRaw.dat, smarthouse.lev=c("NE","NW"))

longi.dat <- designFactors(longi.dat, insertName = "xDays",
                           designfactorMethod="StandardOrder")

#'## Particular edits to longi.dat
longi.dat <- within(longi.dat, 
                    { 
                      Days.after.Salting <- as.numfac(Days) - 29
                    })
longi.dat <- with(longi.dat, longi.dat[order(Snapshot.ID.Tag,Days), ])

#'# Step 3: Form derived traits that result in a value for each observation
#'### Set responses
responses.image <- c("Area")
responses.smooth <- paste(responses.image, "smooth", sep=".")

#'## Form growth rates for each observation of a subset of responses by differencing
longi.dat <- splitContGRdiff(longi.dat, responses.image, 
                             INDICES="Snapshot.ID.Tag",
                             which.rates = c("AGR","RGR"))

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

#'# Step 4: Fit splines to smooth the longitudinal trends in the primary traits and
#'# calculate their growth rates
#'
#'## Smooth responses
#+
for (response in c(responses.image, "Water.Loss"))
  longi.dat <- splitSplines(longi.dat, response, x="xDays", INDICES = "Snapshot.ID.Tag", 
                            df = 4, na.rm=TRUE)
longi.dat <- with(longi.dat, longi.dat[order(Snapshot.ID.Tag, xDays), ])

#'## Loop over smoothed responses, forming growth rates by differences
#+
responses.GR <- paste(responses.smooth, "AGR", sep=".")
longi.dat <- splitContGRdiff(longi.dat, responses.smooth, 
                             INDICES="Snapshot.ID.Tag",
                             which.rates = c("AGR","RGR"))

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
  longiPlot(data = longi.dat, response = responses.longi[k], 
            y.title = titles[k], x="xDays+35.42857143", 
            ggplotFuncs = list(geom_vline(xintercept=29, linetype="longdash", size=1), 
                               scale_x_continuous(breaks=seq(28, 42, by=2)),
                               scale_y_continuous(limits=limits[[klimit]])))
}


#' ### Plot smoothed profiles for all longitudinal responses - GRs by difference
#+ "01-SmoothedProfilesAll"
nresp.smooth <- length(responses.smooth.plot)
limits <- list(c(0,1000), c(0,100), c(0.0,0.40))
for (k in 1:nresp.smooth)
{ 
  longiPlot(data = longi.dat, response = responses.smooth.plot[k], 
            y.title = titles.smooth[k], x="xDays+35.42857143", 
            ggplotFuncs = list(geom_vline(xintercept=29, linetype="longdash", size=1), 
                               scale_x_continuous(breaks=seq(28, 42, by=2)),
                               scale_y_continuous(limits=limits[[klimit]])))
  print(plt)
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
names(cart.dat)[match(c("Area.WUI","Water.Loss.Total"),names(cart.dat))] <- 
        c("Area.Overall.WUE", "Water.Loss.Overall")
cart.dat$Water.Loss.rate.Overall <- cart.dat$Water.Loss.Overall / (42 - 31)

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
  longiPlot(data = longi.SIIT.dat, x="xDays+35.42857143", 
            response = responses.SIIT[k], 
            y.title=responses.SIIT[k], 
            facet.x="Smarthouse", facet.y=".", 
            ggplotFuncs = list(geom_vline(xintercept=29, linetype="longdash", size=1), 
                               scale_x_continuous(breaks=seq(28, 42, by=2)),
                               scale_y_continuous(limits=limits[[klimit]])))
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
tmp<-na.omit(cart.SIIT.dat)
print(summary(tmp[SIIT]))
big.SIIT <- with(tmp, tmp[tmp[SIIT] > 1.15, c("Snapshot.ID.Tag.C","Genotype.ID",
                                              paste(response,"C",sep="."), 
                                              paste(response,"S",sep="."), SIIT)])
big.SIIT <- big.SIIT[order(big.SIIT[SIIT]),]
print(big.SIIT)
plt <- ggplot(tmp, aes_string(SIIT)) +
           geom_histogram(aes(y = ..density..), binwidth=0.05) +
           geom_vline(xintercept=1.15, linetype="longdash", size=1) +
           theme_bw() + facet_grid(Smarthouse ~.)
print(plt)
plt <- ggplot(tmp, aes_string(x="Smarthouse", y=SIIT)) +
           geom_boxplot() + theme_bw()
print(plt)
remove(tmp)
}}
\references{
Al-Tamimi, N, Brien, C.J., Oakey, H., Berger, B., Saade, S., Ho, Y. S., Schmockel, S. M., Tester, M. and Negrao, S. (2016) New salinity tolerance loci revealed in rice using high-throughput non-invasive phenotyping. \emph{Nature Communications}, \bold{7}, 13342.
}
