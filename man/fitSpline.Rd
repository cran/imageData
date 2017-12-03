\name{fitSpline}
\alias{fitSpline}
\title{Produce the fits from a natural cubic smoothing spline applied to a response in a \code{data.frame}}
\description{Uses \code{smooth.spline} to fit a spline to all the values of 
      \code{response} stored in \code{data}. 

      The amount of smoothing can be controlled by \code{df}. 
      If \code{df = NULL}, the amount of 
      smoothing is controlled by the default arguments and those you supply 
      for \code{smooth.spline}. The method of Huang (2001) for correcting the 
      fitted spline for estimation bias at the end-points will be applied if 
      \code{correctBoundaries} is \code{TRUE}.

      The derivatives of the fitted spline can also be obtained, and the 
      Relative Growth Rate (RGR) computed using them, provided 
      \code{correctBoundaries} is \code{FALSE}. Otherwise, growth rates can be 
      obtained by difference using \code{\link{splitContGRdiff}}.

       By default, \code{smooth.spline} will issue an error if there are not 
       at least four distinct x-values. On the other hand, \code{fitSplines} 
       issues a warning and sets all smoothed values and derivatives to 
       \code{NA}. The handling of missing values in the observations is 
       controlled via \code{na.x.action} and \code{na.y.action}.
}
\usage{
fitSpline(data, response, x, df=NULL, correctBoundaries = FALSE, 
          deriv=NULL, suffices.deriv=NULL, RGR=NULL, 
          na.x.action="exclude", na.y.action = "exclude", ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column to be smoothed.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} that is to be smoothed.}
 \item{x}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the predictor variable.}
 \item{df}{A \code{\link{numeric}} specifying the desired equivalent number of degrees 
      of freedom of the smooth (trace of the smoother matrix). Lower values 
      result in more smoothing. If \code{df = NULL}, the amount of smoothing 
      is controlled by the default arguments for and those that you supply to 
      \code{smooth.spline}.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted spline 
            values are to have the method of Huang (2001) applied 
            to them to correct for estimation bias at the end-points. Note that 
            \code{deriv} must be \code{NULL} for \code{correctBoundaries} to be 
            set to \code{TRUE}.}
 \item{deriv}{A \code{numeric} specifying one or more orders of derivatives 
      that are required.}
 \item{suffices.deriv}{A \code{\link{character}} giving the characters to be 
                       appended to the names of the derivatives.}
 \item{RGR}{A \code{\link{character}} giving the character to be appended 
            to the name of the RGR. If \code{RGR} is not \code{NULL}, 
            \code{deriv} must include 1 so that the first derivative is 
            available for calculating the RGR. If \code{RGR} is \code{NULL}, 
            the RGR is not calculated.}
 \item{na.x.action}{A \code{\link{character}} string that specifies the action to 
            be taken when values of \code{x} are \code{NA}. The possible 
            values are \code{fail}, \code{exclude} or \code{omit}. 
            For \code{exclude} and \code{omit}, predictions and derivatives 
            will only be obtained for nonmissing values of \code{x}. 
            The difference between these two codes is that for \code{exclude} 
            the returned \code{data.frame} will have as many rows as 
            \code{data}, the missing values have been incorporated. }
 \item{na.y.action}{A \code{\link{character}} string that specifies the action to 
            be taken when values of \code{y}, or the \code{response}, are 
            \code{NA}.  The possible values are \code{fail}, \code{exclude}, 
            \code{omit}, \code{allx}, \code{trimx}, \code{ltrimx} or 
            \code{rtrimx}. For  all options, except \code{fail}, missing 
            values in \code{y} will be removed before smoothing. 
            For \code{exclude} and \code{omit}, predictions 
            and derivatives will be obtained only for nonmissing values of 
            \code{x} that do not have missing \code{y} values. Again, the 
            difference between these two is that, only for \code{exclude} 
            will the missing values will be incorporated into the 
            returned \code{data.frame}. For \code{allx}, predictions and 
            derivatives will be obtained for all nonmissing \code{x}. 
            For \code{trimx}, they will be obtained for all nonmissing 
            \code{x} between the first and last nonmissing \code{y} values 
            that have been ordered for \code{x}; for \code{ltrimx} and 
            \code{utrimx} either the lower or upper missing \code{y} 
            values, respectively, are trimmed.}
 \item{...}{allows for arguments to be passed to \code{smooth.spline}.}
} 
\value{A \code{data.frame} containing \code{x} and the fitted smooth. The names 
       of the columns will be the value of \code{x} and the value of \code{response} 
       with \code{.smooth} appended. The number of rows in the \code{data.frame} 
       will be equal to the number of pairs that have neither a missing \code{x} or 
       \code{response} and it will have the same order of code{x} as \code{data}. 
       If \code{deriv} is not \code{NULL}, columns 
       containing the values of the derivative(s) will be added to the 
       \code{data.frame}; the name each of these columns will be the value of 
       \code{response} with \code{.smooth.dvf} appended, where 
       \code{f} is the order of the derivative, or  the value of \code{response} 
       with \code{.smooth.} and the corresponding element of 
       \code{suffices.deriv} appended. If \code{RGR} is not \code{NULL}, the RGR 
       is calculated as the ratio of value of the first derivative of the fitted 
       spline and the fitted value for the spline.}
\author{Chris Brien}
\seealso{\code{\link{splitSplines}}, \code{\link{smooth.spline}},
         \code{\link{predict.smooth.spline}}, \code{\link{splitContGRdiff}}}
\references{
Huang, C. (2001). Boundary corrected cubic smoothing splines. \emph{Journal of Statistical Computation and Simulation}, \bold{70}, 107-121.}
\examples{
data(exampleData)
fit <- fitSpline(longi.dat, response="Area", , x="xDays", df = 4,
                 deriv=c(1,2), suffices.deriv=c("AGRdv","Acc"))
}
\keyword{data}
\keyword{manip}