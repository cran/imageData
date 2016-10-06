\name{news}
\title{News for Package imageData}
\section{imageData Version 0.1-26 (2016-10-06)}{
  \itemize{
    \item{Rebuild under R 3.4.0.}
  }}
\section{imageData Version 0.1-25 (2016-08-23)}{
  \itemize{
    \item{For importExcel, changed argument \code{timeAfterPlanting} to \code{timeAfterStart}  and argument \code{plenting.time} to \code{startTime}.}
    \item{For calcTimes, changed argument \code{originTime} to \code{startTime}.}
  }}
\section{imageData Version 0.1-24 (2016-08-19)}{
  \itemize{
    \item{Added a new function calcTimes to calculate intervals and time positions within intervals from imaging times. It allows a choice of unit for the interval.}
    \item{Removed the calcuation of time values from imagtimesPlot and added arguments \code{intervals} and \code{timePositions} to allow the use to specify the names of the columns to be plotted.}
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
    \item{Removed the creation of design factors and covariates from longitudinalPrime and ade function designFactors to dot his separately, if needed.}
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
    \item Add probeDF that allows the comparison fo the smooths resulting from different values of df used for the smoothing spline.
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
