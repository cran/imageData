\name{splitSplines}
\alias{splitSplines}
\title{Adds the fits after fitting a natural cubic smoothing spline to subsets 
of a response to a \code{data.frame}}
\description{Uses \code{smooth.spline} to fit a spline to a subset of the values 
             of \code{response} and stores the fitted values in \code{data}. 
             The subsets are those values with the same levels combinations 
             of the factors listed in \code{INDICES} and the degrees of 
             smoothing is controlled by \code{df}. The derivatives
             of the fitted spline can also be obtained, as can the Relative 
             Growth Rates (RGR).}
\usage{
splitSplines(data, response, x, INDICES, df = NULL, deriv = NULL, 
             suffices.deriv=NULL, RGR=NULL, na.rm = FALSE, sep=".")}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column to be smoothed.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} that is to be smoothed.}
 \item{x}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the predictor variable.}
 \item{INDICES}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of \code{response} 
             that are to be smoothed separately. If the columns corresponding to 
             \code{INDICES} are not \code{\link{factor}}(s) then they will be 
             coerced to \code{\link{factor}}(s). The subsets are formed 
             using \code{\link{split}}.}
 \item{df}{A \code{numeric} specifying the desired equivalent number of degrees 
      of freedom of the smooth (trace of the smoother matrix). Lower values 
      result in more smoothing. If \code{df = NULL}, ordinary leave-one-out 
      cross-validation is used to determine the amount of smooth.}
 \item{deriv}{A \code{numeric} specifying one or more orders of derivatives 
      that are required.}
 \item{suffices.deriv}{A \code{\link{character}} giving the characters to be 
                       appended to the names of the derivatives.}
 \item{RGR}{A \code{\link{character}} giving the character to be appended 
            to the name of the RGR. If \code{RGR} is not \code{NULL}, 
            \code{deriv} must include 1 so that the the first derivative is 
            available for caluclating the RGR. If \code{RGR} is \code{NULL}, 
            the RGR is not calculated.}
 \item{na.rm}{A \code{logical} indicating whether or not rows of \code{data} 
              with \code{NA}s are to be removed prior to smoothing.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{INDICES} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a factor to delimit 
            levels when the levels of \code{INDICES} are combined to identify 
            subsets.}
}
\value{A \code{\link{data.frame}} containing  \code{data} to which has been 
       added a column with the fitted smooth, the name of the column being 
       \code{response} with \code{.smooth} apended to it. If \code{deriv} is 
       not \code{NULL}, columns containing the values of the derivative(s) 
       will be added to \code{data}; the name each of these columns will 
       be the value of \code{response} with \code{.smooth.dvf} appended, 
       where \code{f} is the order of the derivative, or  the value of 
       \code{response} with \code{.smooth.} and the corresponding 
       element of \code{suffices.deriv} appended. If \code{RGR} is not 
       \code{NULL}, the RGR is calculated as the ratio of value of the first 
       derivative of the fitted spline and the fitted value for the spline.}
\author{Chris Brien}
\seealso{\code{\link{fitSpline}}, \code{\link{smooth.spline}}, 
         \code{\link{predict.smooth.spline}}, \code{\link{splitContGRdiff}}, \code{\link{split}}}
\examples{
data(exampleData)
longi.dat <- splitSplines(longi.dat, response="Area", x="xDays", 
                          INDICES = "Snapshot.ID.Tag", 
                          df = 4, deriv=1, suffices.deriv="AGRdv", RGR="RGRdv")
}
\keyword{data}
\keyword{manip}