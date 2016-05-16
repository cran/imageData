\name{imagetimesPlot}
\alias{imagetimesPlot}
\title{Plots the hour of the day carts are imaged against the days after planting 
(or some other number of days after an event)}
\description{Uses \code{ggplot} to produce a plot of the hour of the day carts 
             are imaged against the days after planting (or some other number of 
             days after an event). A line is produced for each \code{Lane} and a 
             separate plot for each \code{Smarthouse}. It aids in checking 
             \code{times} and what occurred in imaging the plants.}
\usage{
imagetimesPlot(data, times="Time.after.Planting..d.")}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data to be plotted and 
             including the columns \code{Snapshot.ID.Tag}, which uniquely indexes 
             the carts in the experiment, and \code{Snapshot.Time.Stamp}, which 
             reflects the time of the imaging from which a particular data value 
             was obtained.}
 \item{times}{A \code{\link{character}} giving the \code{\link{vector}} of calculated times 
              after planting (or after some time origin) to be plotted on the x-axis.}
}
\value{An object of class "\code{ggplot}", which can be plotted using \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{ggplot}}.}
\examples{
data(exampleData)
imagetimesPlot(data = raw.dat)
}
\keyword{hplot}

