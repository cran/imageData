\name{news}
\title{News for Package imageData}
\section{imageData Version 0.1-41 (2017-12-03)}{
  \itemize{
    \item{Rebuild under R-3.5.0.}
  }}
\section{imageData Version 0.1-40 (2017-08-31)}{
  \itemize{
    \item{Modify splitContGRdiff to ensure that initial columns are INDICEs followed by times.factor.}
  }}
\section{imageData Version 0.1-39 (2017-07-31)}{
  \itemize{
    \item{Change default for \code{smarthouse.lev} in longitudinalPrime to \code{NULL}, in which case unique values in the Smarthouse column will be used.}
  }}
\section{imageData Version 0.1-37 (2017-05-22)}{
  \itemize{
    \item{Recompiled the package because of problems loading it.}
  }}
\section{imageData Version 0.1-36 (2017-04-22)}{
  \itemize{
    \item{Added argument \code{correctBoundaries} to probeDF to match fitSpline and splitSplines.}
  }}
\section{imageData Version 0.1-35 (2017-04-20)}{
  \itemize{
    \item{Added argument \code{sep} to importExcel to allow for separators other than `,' in a \code{csv} file.}
  }}
\section{imageData Version 0.1-34 (2017-04-19)}{
  \itemize{
    \item{Added argument \code{correctBoundaries} to fitSpline and splitSplines to implement the Huang (2001) method for correcting bias in the boundary estimates of a fitted spline.}
  }}
\section{imageData Version 0.1-33 (2017-03-20)}{
  \itemize{
    \item{Fixed bug in splitContGRdiff when there are multiple missing values for some \code{INDICES} or for the \code{time.factor}. This causes an identifiability problem for the merge function.}
  }}
\section{imageData Version 0.1-32 (2017-02-23)}{
  \itemize{
    \item{Extensively revised the handling of missing values (\code{NA}) in fitSpline so that it and calling functions, such as splitSplines and probeDF, will give a warning, rather than an error, when all data are missing or there are less than 4 distinct x values. All fitted values and derivates will be set to \code{NA}. Introduced the arguments \code{na.x.action} and \code{na.y.action} that allow control of what happens when there are missing values. In particular options have been introduced to allow for predictions and derivatives to be obtained when there are missing response values.}
  }}
\section{imageData Version 0.1-31 (2016-12-21)}{
  \itemize{
    \item{Revised some documentation}
  }}
\section{imageData Version 0.1-30 (2016-12-09)}{
  \itemize{
    \item{Rebuild with ggplot2 2.2.0.}
  }}
\section{imageData Version 0.1-29 (2016-11-16)}{
  \itemize{
    \item{For longitudinalPrime, added arguments \code{cartId}, \code{timeAfterStart} and \code{imageTimes}.}
    \item{For probeDF, revised \code{which.plots} argument and added \code{which.traits} argument. The result is that there is more flexibility in which plots are produced for which traits.}
        \item{Fixed a bug in probeDF that resulted in an error when other than a single factor was entered for \code{facet.x} or \code{facet.y}.}
  }}
\section{imageData Version 0.1-28 (2016-10-12)}{
  \itemize{
    \item{For getDates, added arguments include.individuals and individuals to allow a column of identifiers for individuals to be included in the returned data.frame.}
  }}
\section{imageData Version 0.1-27 (2016-10-07)}{
  \itemize{
    \item{Centre xZones by removing the mean of xZones in designFactors.}
  }}
\section{imageData Version 0.1-26 (2016-10-06)}{
  \itemize{
    \item{Rebuild under R 3.4.0.}
  }}
\section{imageData Version 0.1-25 (2016-08-23)}{
  \itemize{
    \item{For importExcel, changed argument \code{timeAfterPlanting} to \code{timeAfterStart}  and argument \code{planting.time} to \code{startTime}.}
    \item{For calcTimes, changed argument \code{originTime} to \code{startTime}.}
  }}
\section{imageData Version 0.1-24 (2016-08-19)}{
  \itemize{
    \item{Added a new function calcTimes to calculate intervals and time positions within intervals from imaging times. It allows a choice of unit for the interval.}
    \item{Removed the calculation of time values from imagtimesPlot and added arguments \code{intervals} and \code{timePositions} to allow the use to specify the names of the columns to be plotted.}
  }}
\section{imageData Version 0.1-23 (2016-08-17)}{
  \itemize{
    \item{Add \code{...} argument to importExcel to allow passing of arguments to imagetimesPlot.}
    \item{Add \code{ggplotFuncs} argument to imagetimesPlot and dealt with a \code{factor} being supplied to its \code{intervals} argument.}
  }}
\section{imageData Version 0.1-22 (2016-07-01)}{
  \itemize{
    \item{Add arguments to importExcel that enables the user to supply names for cartId, imageTime and timeAfterPlanting.}
    \item{Add idvars argument to longitudinalPrime to allow some flexibility in the naming and inclusion of genotype and treatment names.}
    \item{Fix a bug in probeDF where supplying fact as "." caused a failure.}
    \item{Allow no facets in longiPlot.}
  }}
\section{imageData Version 0.1-21 (2016-05-16)}{
  \itemize{
    \item{Further revised package documentation.}
    \item{Included RiceRaw.dat, a data set from the Plant Accelerator, and revised the example in imageData-pkg to be based on it.}
    \item{Included a small data set and modified examples to use it so that most examples in the manual can be run.}
    \item{Fixed bugs in PVA, intervalPVA, rcontrib and corrPlot.}
  }}
\section{imageData Version 0.1-20 (2016-05-11)}{
  \itemize{
    \item{Revised fitSpline to convert \code{Inf} values into \code{NA}.}
    \item{Fix bug in designFactors when Lane and Position are not \code{factors}.}
    \item{Revised package documentation.}
  }}
\section{imageData Version 0.1-19 (2016-03-22)}{
  \itemize{
    \item{Revised documentation for intervalGRaverage and intervalGRdiff.}
  }}
\section{imageData Version 0.1-18 (2016-03-16)}{
  \itemize{
    \item{Modified splitContGRdiff so that it does not assume \code{data} is ordered by \code{INDICES} then \code{times.factor}.}
    \item{Fix bug in fitSpline when \code{suffices.deriv} is \code{NULL}.}
  }}
\section{imageData Version 0.1-17 (2016-03-09)}{
  \itemize{
    \item{Fix bug in dealing with missing values in splitSplines.}
    \item{Revise anomPlot, intervalGRaverage, intervalGRdiff, intervalPVA, 
          intervalValueCalculate and intervalGRWUI documentation.}
  }}
\section{imageData Version 0.1-16 (2016-03-02)}{
  \itemize{
    \item{Add message identifying absent variables when calling longitudinalPrime.}
    \item{Fix bug in passing \code{times.factor} in intervalGRaverage, intervalGRdiff and intervalWUI.}
  }}
\section{imageData Version 0.1-15 (2016-02-23)}{
  \itemize{
    \item{Fix bug in importing \code{csv} files using importExcel.}
  }}
\section{imageData Version 0.1-14 (2016-02-15)}{
  \itemize{
    \item{Added the ability to import \code{csv} files using importExcel and to specify the 
          format of times.}
  }}
\section{imageData Version 0.1-13 (2015-12-18)}{
  \itemize{
    \item{Added ggplotFuncs and labeller arguments to longiPlot.}
    \item(Fixed bug in designFactors)
    \item{Forced \code{Time after Planting [d]} to be of class \code{difftime}.}
  }}
\section{imageData Version 0.1-12 (2015-11-11)}{
  \itemize{
    \item{Removed the creation of design factors and covariates from longitudinalPrime and add function designFactors to do this separately, if needed.}
    \item{replaced `Water' in `Water.Total' in intervalWUI with the value of \code{water.name}.}
  }}
\section{imageData Version 0.1-11 (2015-11-02)}{
  \itemize{
    \item{Added \code{na.rm} argument to intervalWUI.}
    \item{replaced `Water' in `Water.Total' in intervalWUI with the value of \code{water.name}.}
  }}
\section{imageData Version 0.1-10 (2015-10-19)}{
  \itemize{
    \item{Fixed bug in calcLagged in which \code{operation} was being ignored.}
    \item{changed `Water.Use' to Water.`Loss'.}
  }}
\section{imageData Version 0.1-9 (2015-10-16)}{
  \itemize{
    \item Revised longitudinalPrime to include new arguments designfactorMethod and csvFilename to provide several methods for specifying the design factors Smarthouse, Zones, ZLane, Mainplots and Subplots. Added arguments to allow the specification of the numbers of levels of these factors. Added argument calcWaterUse  to allow Water.Use to have been provided via the raw data file.
    \item Minor bug fix to imagetimesPlot.
    \item Added \code{na.rm} argument to probeDF.
  }}
\section{imageData Version 0.1-8 (2015-09-30)}{
  \itemize{
    \item Revised probeDF to allow more flexibility in what is calculated and plotted.
  }}
\section{imageData Version 0.1-7 (2015-09-24)}{
  \itemize{
    \item Fixed a bug in probeDF.
  }}
\section{imageData Version 0.1-6 (2015-09-15)}{
  \itemize{
    \item Add groupsFactor argument to anomPlot.
    \item Generalized calcLagged to allow for negative lag values.
    \item Fixed a bug in twoLevelOpcreate when there are repeats in the vector of responses and add columns.suffixed argument.
  }}
\section{imageData Version 0.1-5 (2015-08-16)}{
  \itemize{
    \item Add probeDF that allows the comparison for the smooths resulting from different values of df used for the smoothing spline.
    \item Fix bug in calculating sum of Water Use and add water.use argument to allow for different names for the column.
    \item Make NULL default for lower and upper in anomPlot so can have intervals that are open at one end
  }}
\section{Changes in imageData version 0.1-4 (2015-07-25)}{
  \itemize{
    \item Added function anomPlot to identify anomalous individuals and to produce plots with and without them.
    \item Renamed function intervalStat to splitValueCalculate and added function intervalValueCalculate.
    \item Renamed functions splineIndividual and splineResponse to fitSpline and 
          splitSplines, respectively.
    \item Added function splitContGRdiff to compute continuous growth rates by differencing 
           for a response split into subsets.
    \item Added function twoLevelOpcreate to compute Osmotic Stress Tolerance (OST) indices.
    \item Changed output argument in longiPlot to printPlot and added colours.column and colours to allow plotting with multiple colours.
    \item Revised the documentation in the manual.
    \item NEWS added.
  }
}
\section{imageData Version 0.1-3}{
  \itemize{
    \item Early developmental version of the package.
  }}
\section{imageData Version 0.1-2}{
  \itemize{
    \item Early developmental version of the package.
  }}
\section{imageData Version 0.1-1}{
  \itemize{
    \item First developmental version of the package.
  }
}
