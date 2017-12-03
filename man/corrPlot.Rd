\name{corrPlot}
\alias{corrPlot}
\title{Calculates and plots correlation matrices for a set of responses}
\description{Having calculated the correlations a heat map indicating the magnitude of the 
       correlations is produced using \code{ggplot}. In this heat map, the darker the red in 
       a cell then the closer the correlation is to -1, while the deeper the blue in the cell,
       then the closer the correlation is to 1. Also produced is a matrix plot of all pairwise 
       combinations of the variables. The matrix plot contains a scatter diagram for each pair, 
       as well as the value of the correlation coefficient. The argument \code{pairs.sets} 
       can be used to restrict the pairs in the matrix plot to those combinations within each 
       set. }
\usage{
corrPlot(responses, data, show.sig = FALSE, title = NULL, 
         pairs.sets = NULL, labelSize = 4, ...)}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} containing the variables to be correlated.}
 \item{data}{A \code{\link{data.frame}} containing the columns of variables to be correlated.}
 \item{show.sig}{A \code{\link{logical}} indicating whether or not to give asterisks indicating 
                 significance on the plot.}
 \item{title}{Title for the plot.}
 \item{pairs.sets}{A \code{\link{list}} each of whose components is a \code{\link{numeric}} giving 
                the position of the variable names in \code{responses} that are to be included 
                in the set. All pairs of variables in this pairs.set will be included in a 
                matrix plot.}
 \item{labelSize}{A \code{\link{numeric}} giving the size of the labels in the matrix plot.}
 \item{...}{allows passing of arguments to other functions}
}
\value{\code{NULL}.}
\author{Chris Brien}
\seealso{\code{\link{ggplot}}.}
\examples{\dontrun{
data(exampleData)
responses <- c("Area","Area.SV","Area.TV", "Image.Biomass", "Max.Height","Centre.Mass",
               "Density", "Compactness.TV", "Compactness.SV")
corrPlot(responses, longi.dat, pairs.sets=list(c(1:4),c(5:7)))

}}
\keyword{hplot}

