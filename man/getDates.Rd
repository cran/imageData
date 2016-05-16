\name{getDates}
\alias{getDates}
\title{Forms a subset of \code{responses} in \code{data} that contains their values for the nominated times}
\description{Forms a subset of \code{responses} in \code{data} that contains their values for the 
             nominated times.}
\usage{
getDates(responses, times.factor = "Days", data, which.times, 
         suffix = NULL, include.times.factor = FALSE)}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} whose values are to be subsetted.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels will be used to identify the subset and 
             should be numeric values stored as characters.}
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             growth rates are to be calculated.}
 \item{which.times}{A \code{\link{vector}} giving the times that are to be selected.}
 \item{suffix}{A \code{\link{character}} giving the suffix to be appended to 
               \code{responses} to form the names of the columns containing the 
               subset.}
 \item{include.times.factor}{A \code{\link{logical}} indicating whether or not 
             to include the \code{times.factor} in the result.}
}
\value{A \code{\link{data.frame}} containing the subset of \code{responses}, the 
       names of the columns being the concatenation of  \code{responses} and 
       \code{suffix} separated by a full stop.}
\author{Chris Brien}
\examples{
data(exampleData)
AreaLast <- getDates("Area.smooth", data = longi.dat, 
                     which.times = c(42), suffix = "last")
}
\keyword{data}
\keyword{manip}