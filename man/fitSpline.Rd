\name{fitSpline}
\alias{fitSpline}
\title{Produce the fits from a natural cubic smoothing spline applied to a response in a \code{data.frame}}
\description{Uses \code{smooth.spline} to fit a spline to all the values of 
             \code{response} stored in \code{data}. 
             The amount of smoothing is controlled by \code{df}. The derivatives 
             of the fitted spline can also be obtained, as can the Relative 
             Growth Rate (RGR).}
\usage{
fitSpline(data, response, x, df=NULL, deriv=NULL, suffices.deriv=NULL, 
          RGR=NULL, na.rm=FALSE)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column to be smoothed.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} that is to be smoothed.}
 \item{x}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the predictor variable.}
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
} 
\value{A \code{data.frame} containing \code{x} and the fitted smooth. The names 
       of the columns will be the value of \code{x} and the value of \code{response} 
       with \code{.smooth} appended. If \code{deriv} is not \code{NULL}, columns 
       containing the values of the derivative(s) will be added to the 
       \code{data.frame}; the name each of of these columns will be the value of 
       \code{response} with \code{.smooth.dvf} appended, where 
       \code{f} is the order of the derivative, or  the value of \code{response} 
       with \code{.smooth.} and the corresponding element of 
       \code{suffices.deriv} appended. If \code{RGR} is not \code{NULL}, the RGR 
       is calculated as the ratio of value of the first derivative of the fitted 
       spline and the fitted value for the spline.}
\author{Chris Brien}
\seealso{\code{\link{splitSplines}}, \code{\link{smooth.spline}},
         \code{\link{predict.smooth.spline}}, \code{\link{splitContGRdiff}}}
\examples{
data(exampleData)
fit <- fitSpline(longi.dat, response="Area", , x="xDays", df = 4,
                 deriv=c(1,2), suffices.deriv=c("AGRdv","Acc"))
}
\keyword{data}
\keyword{manip}