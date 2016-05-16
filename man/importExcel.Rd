\name{importExcel}
\alias{importExcel}
\title{Imports an Excel imaging file and allows some renaming of variables}
\description{Uses \code{XLConnect} to import a sheet of imaging data produced by the 
  Lemna Tec Scanalyzer. Basically, the data consists of imaging data obtained from a 
  set of pots or carts over time. There should be a unique identifier for each cart, 
  which by default is \code{Snapshot.ID.Tag} and a column labelled \code{Snapshot.Time.Stamp}, 
  which reflects the time of the imaging pertaining to the data values values in the 
  row of a sheet. Also, \code{Time.After.Planting..d.} can be caclulated, or recaclulated if 
  already prsent, from \code{Snapshot.Time.Stamp} by subtracting a supplied time of planting. 
  
  If the data is prefixed by 'RGB_SV1', 'RGB_SV2' or 'RGB_TV', 
  the 'RGB_' can be removed and the 'SV1', 'SV2' or 'TV' becomes a suffix. 
  }
\usage{
importExcel(file, sheet="raw data", prefix2suffix = TRUE, planting.time = NULL,
            timeformat = "\%Y-\%m-\%d \%H:\%M", imagetimesPlot = TRUE)}
\arguments{
 \item{file}{A \code{\link{character}} giving the path and name of the file containing 
             the data.}
 \item{sheet}{A \code{\link{character}} giving the name of the sheet containing 
             the data, that must include columns named \code{Snapshot.ID.Tag}, 
             which uniquely indexes the carts in the experiment, and 
             \code{Snapshot.Time.Stamp}, which reflects the time of the imaging from 
             which a particular data value was obtained. It is also assumed that a 
             column named \code{Time.after.Planting..d.} is in the sheet or that it will 
             be calculated from \code{Snapshot.Time.Stamp} using the value of 
             \code{planting.time} supplied in the function call.}
 \item{prefix2suffix}{A \code{logical} specifying whether the variables names prefixed 
             by "RGB_SV1", "RGB_SV2" or "RGB_TV" have these prefixes replaced by the 
             suffices ".SV1", ".SV2" or ".TV", respectively.}
 \item{planting.time}{A \code{character} giving the time of planting, in the \code{POSIXct} 
       format timeformat, to be subtracted from \code{Snapshot.Time.Stamp} 
       in recalculating \code{Time.After.Planting..d.}. If \code{planting.time} is 
       \code{NULL} then \code{Time.After.Planting..d.} is not recalculated.} 
 \item{timeformat}{A \code{character} giving the \code{POSIXct} format of characters 
       containing times, in particular \code{Snapshot.Time.Stamp} and \code{planting.time}.} 
 \item{imagetimesPlot}{A \code{logical} indicating whether a plot of the imaging times against 
             the recalculated \code{Time.After.Planting..d.}. It aids in checking 
             \code{Time.After.Planting..d.} and what occurred in imaging the plants.}
}
\value{A \code{\link{data.frame}} containing the data.}
\note{\code{XLConnect} uses Java and can take a lot of memory. You may have to increase 
      the amount of RAM allocated to the Java heap. For example, 
      \code{options(java.parameters = "-Xmx10g")} allocates 10 gigabytes of RAM, which 
      has been needed on occasion to read Excel imaging files. Java itself must also be 
      configured to allocate at least this amount of memory.
      
      Also, note that XLConnect cannot import cells with some Excel functions. For example,
      it cannot import cells that use {ROUNDOWN} in a formula. The problem can be avoided by 
      converting the formulae into values.}
\author{Chris Brien}
\seealso{\code{\link{as.POSIXct}}}
\examples{\dontrun{
raw.0169.dat <- importExcel(file = "0169 analysis_20140603.xlsx", 
                            planting.time = "2013-05-23 8:00 AM")
}}
\keyword{data}
\keyword{manip}