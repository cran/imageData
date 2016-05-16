\name{probeDF}
\alias{probeDF}
\title{Compares, for a set of specified values of df, a response and the smooths of it, possibly along with growth rates calculated from the smooths}
\description{Takes a \code{response} and, for each \code{individual}, uses 
             \code{smooth.spline} to smooth its values for each of the degrees 
             of freedom values in \code{df}. Provided \code{get.rates} is \code{TRUE},
             both the Absolute Growth Rates (AGR) and the Relative Growth Rates (RGR) 
             are calculated for each smooth, either using differences or first 
             derivatives. A combination of the unsmoothed and smoothed 
             values, as well as the AGR and RGR, can be plotted for each value in 
             \code{df}. Note that the arguments that modify the plots apply to all 
             plots that are produced.}
\usage{
probeDF(data, response = "Area", xname="xDays", individuals="Snapshot.ID.Tag", 
        na.rm = FALSE, df, get.rates = TRUE, rates.method="differences", 
        times.factor = "Days", x = NULL, 
        facet.x = "Treatment.1", facet.y =   "Smarthouse", 
        which.plots = c("smoothed", "AGR", "RGR"), ggplotFuncs = NULL, 
        ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data.}
 \item{response}{A \code{character} specifying the response variable to be 
                 supplied to \code{\link{smooth.spline}} and that 
                 is to be plotted on the y-axis.}
 \item{xname}{A \code{\link{character}} giving the name of the
             \code{\link{numeric}} that contains the values of the predictor 
              variable to be supplied to \code{\link{smooth.spline}}.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of the \code{data} 
             for which each subset corresponds to the \code{response} values for 
             an individual.}
 \item{na.rm}{A \code{logical} indicating whether or not rows of \code{data} 
              with \code{NA}s are to be removed prior to smoothing.}
 \item{df}{A \code{\link{numeric}} specifying the set of degrees of freedom to 
              be probed.}
 \item{get.rates}{A \code{\link{logical}} specifying whether or not the growth 
              rates (AGR and RGR) are to be computed and stored.}
 \item{rates.method}{A \code{\link{character}} specifying the method to use in 
               calculating the growth rates. The two possibilities are 
               \code{"differences"} and \code{"derivates"}.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels will be used in calculating growth rates and 
             should be numeric values stored as characters.}
 \item{x}{A \code{\link{character}} giving the variable to be plotted on the 
          x-axis. If \code{x} is \code{NULL} then \code{xname} is used.}
 \item{facet.x}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate columns of plots. 
                Use \code{"."} if a split into columns is not wanted.}
 \item{facet.y}{A \code{\link{data.frame}} giving the variable to be used to 
                form subsets to be plotted in separate rows of plots.
                Use \code{"."} if a split into columns is not wanted.}
 \item{which.plots}{A \code{\link{character}} giving the plots that are to be 
             produced. If \code{none}, no plots are produced. If  \code{all}, 
             plots of the unsmoothed response and, for each value of \code{df}, 
             plots of the smoothed response, the AGR and  RGR are produced.  If 
             \code{responseComparison}, a combined plot of  the unsmoothed 
             response  and the smoothed response is produced for for each value 
             of \code{df}.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element. Note that 
                these functions are applied to all three plots produced.}
 \item{...}{allows passing of arguments to \code{\link{longiPlot}}.}
}
\value{A \code{\link{data.frame}} containing \code{individuals}, 
       \code{times.factor}, \code{facet.x}, \code{facet.y}, \code{xname}, 
       \code{response}, and, for each \code{df}, the smoothed 
       response, the AGR and the RGR. It is returned invisibly. The names of 
       the new data are constructed by joining elements separated by full 
       stops (\code{.}). In all cases, the last element is the value of 
       \code{df}. For the smoothed response, the other elements are 
       \code{response} and \code{"smooth"}; for AGR and RGR, the other elements 
       are the name of the smoothed response and either \code{"AGR"} or 
       \code{"RGR"}.}
\author{Chris Brien}
\seealso{\code{\link{splitSplines}}, \code{\link{splitContGRdiff}}, \code{\link{smooth.spline}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
vline <- list(ggplot2::geom_vline(xintercept=20, linetype="longdash", size=1),
              ggplot2::scale_x_continuous(breaks=seq(12, 36, by=2)))
probeDF(data = longi.dat, response = "Area", df = c(4,7), x="xDays+24.16666667", 
        ggplotFuncs=vline)
}
\keyword{hplot}
\keyword{manip}