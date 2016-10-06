\name{importExcel}
\alias{importExcel}
\title{Imports an Excel imaging file and allows some renaming of variables}
\description{Uses \code{XLConnect} to import a sheet of imaging data produced by the 
  Lemna Tec Scanalyzer. Basically, the data consists of imaging data obtained from a 
  set of pots or carts over time. There should be a column, which by default is called 
  \code{Snapshot.ID.Tag}, containing a unique identifier for each cart and a column, 
  which by default is labelled \code{Snapshot.Time.Stamp}, containing 
  the time of imaging for each observation in a row of the sheet. Also, if 
  \code{startTime} is not \code{NULL}, \code{calcTimes} is called to 
  calculate, or recaclulate if already present, \code{timeAfterStart} from 
  \code{imageTimes} by subtracting a supplied \code{startTime}. 
  
  If the data is prefixed by 'RGB_SV1', 'RGB_SV2' or 'RGB_TV', 
  the 'RGB_' can be removed and the 'SV1', 'SV2' or 'TV' becomes a suffix. 
  }
\usage{
importExcel(file, sheet="raw data", 
            cartId = "Snapshot.ID.Tag", 
            imageTimes = "Snapshot.Time.Stamp", 
            timeAfterStart = "Time.after.Planting..d.", 
            prefix2suffix = TRUE, 
            startTime = NULL,
            timeFormat = "\%Y-\%m-\%d \%H:\%M", 
            imagetimesPlot = TRUE, ...)}
\arguments{
 \item{file}{A \code{\link{character}} giving the path and name of the file containing 
             the data.}
 \item{sheet}{A \code{\link{character}} giving the name of the sheet containing 
             the data, that must include columns whose names are as specified by \code{cartId}, 
             which uniquely indexes the carts in the experiment, and 
             \code{imageTimes}, which reflects the time of the imaging from 
             which a particular data value was obtained. It is also assumed that a 
             column whhose name is specified by \code{timeAfterStart} is in the sheet or that 
             it will be calculated from \code{imagTimes} using the value of 
             \code{startTime} supplied in the function call.}
 \item{cartId}{A \code{\link{character}} giving the name of the column that contains
               the unique Id for each cart. Note that in importing data into R, spaces 
               and nonalphanumeric characters in names are converted to full stops.}
 \item{imageTimes}{A \code{\link{character}} giving the name of the column that contains
               the time that each cart was imaged. Note that in importing data into R, spaces 
               and nonalphanumeric characters in names are converted to full stops.}
 \item{timeAfterStart}{A \code{\link{character}} giving the name of the column that 
               contains or is to contain the difference between \code{imageTimes} and 
               \code{startTime}. The function \code{calcTimes} is called to caclulate the
               differences. For example, it might contain the number of days after 
               planting. Note that in importing data into R, spaces 
               and nonalphanumeric characters in names are converted to full stops.}
 \item{prefix2suffix}{A \code{logical} specifying whether the variables names prefixed 
             by "RGB_SV1", "RGB_SV2" or "RGB_TV" have these prefixes replaced by the 
             suffices ".SV1", ".SV2" or ".TV", respectively.}
 \item{startTime}{A \code{character} giving the time of planting, in the \code{POSIXct} 
       format timeFormat, to be subtracted from \code{imageTimes} 
       in recalculating \code{timeAfterStart}. If \code{startTime} is 
       \code{NULL} then \code{timeAfterStart} is not recalculated.} 
 \item{timeFormat}{A \code{character} giving the \code{POSIXct} format of characters 
       containing times, in particular \code{imageTimes} and \code{startTime}.} 
 \item{imagetimesPlot}{A \code{logical} indicating whether a plot of the imaging times against 
             the recalculated \code{Time.After.Planting..d.}. It aids in checking 
             \code{Time.After.Planting..d.} and what occurred in imaging the plants.}
 \item{... }{allows for arguments to be pasesd to \code{\link{imagetimesPlot}}. However, if 
             \code{intervals} is passed an error will occur; use timeAfterStart instead.}
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
\seealso{\code{\link{as.POSIXct}}, \code{\link{calcTimes}}, \code{\link{imagetimesPlot}}}
\examples{\dontrun{
raw.0169.dat <- importExcel(file = "0169 analysis_20140603.xlsx", 
                            startTime = "2013-05-23 8:00 AM")
}}
\keyword{data}
\keyword{manip}