\name{longitudinalPrime}
\alias{longitudinalPrime}
\title{Selects a set variables to be retained in a data frame of longitudinal data}
\description{Forms the prime traits by selecting a subset of the traits in a data.frame of imaging data produced by the 
  Lemna Tec Scanalyzer. Add factors and explanatory variates that might be of use, including the maximum 
  of the heigh SV1 and SV2 andf the mean of the two SV measurements for several traits.}
\usage{
longitudinalPrime(data, idcolumns = c("Genotype.ID","Treatment.1"), 
                  smarthouse.lev = c("SW"), 
                  calcWaterLoss = TRUE, pixelsPERcm = 18)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns specified in \code{idcolumns} 
             along with  the following columns: 
          
            Snapshot.ID.Tag, Time.after.Planting..d., Smarthouse, Lane, Position, 
            
            Snapshot.Time.Stamp, Weight.Before, 
            
            Weight.After, Water.Amount, Projected.Shoot.Area..pixels., Area.SV1, 
            
            Boundary.Points.To.Area.Ratio.SV1, Caliper.Length.SV1, 
            
            Center.Of.Mass.Y.SV1, Compactness.SV1, Convex.Hull.Area.SV1, 
            
            Max.Dist.Above.Horizon.Line.SV1, Area.SV2, 
            
            Boundary.Points.To.Area.Ratio.SV2, Caliper.Length.SV2, 
            
            Center.Of.Mass.Y.SV2, Compactness.SV2, Convex.Hull.Area.SV2, 
            
            Max.Dist.Above.Horizon.Line.SV2, Area.TV, 
            
            Boundary.Points.To.Area.Ratio.TV,  Caliper.Length.TV, 
            
            Compactness.TV, Convex.Hull.Area.TV.}
 \item{idcolumns}{A \code{\link{character}} vector giving the names of the columns that identify the plants in a cart e.g. Genotype.ID, Treatment.1, Treatment.2.}
 \item{smarthouse.lev}{A \code{\link{character}} vector giving the levels to use 
            for the Smarthouse factor.}
 \item{calcWaterLoss}{A \code{\link{logical}} indicating whether to calculate the Water.Loss. If it is \code{FALSE}, it is assumed that it is already in \code{data}.}
 \item{pixelsPERcm }{A \code{\link{numeric}} giving the number of pixels per cm 
           for the images.}
}
\value{A \code{\link{data.frame}} containing the columns:
       
       \enumerate{
         \item Smarthouse: factor with levels for the Smarthouse
         \item Lane: factor for lane number in a smarthouse
         \item Position: factor for east/west position in a lane
         \item Days: factor for the number of Days After Planting (DAP)
         \item Snapshot.ID.Tag: unique code for each cart 
         \item Reps: factor indexing the replicates for each genotype (calculated)
         \item Snapshot.Time.Stamp: time at which an image was taken in POSIXct format
         \item xPosn: numeric for the Positions within a Lane (calculated)
         \item Hour: hour of the day, to 2 decimal places, at which the image was taken (calculated)
         \item xDays: numeric for the DAP that is centred by subtracting the mean of the unique days (calculated)
         \item idcolumns: the columns listed in \code{idcolumns} that have beeen converted to factors
         \item Weight.Before: weight of the pot before watering
         \item Weight.After: weight of the pot after watering
         \item Water.Amount: the weight of the water added (= Water.After - Water.Before)  (calculated)
         \item Water.Loss: the difference between Weight.Before for the current imaging and the Weight.After for the previous imaging (calculated unless \code{calcWaterLoss} is \code{FALSE})
         \item Area: the Projected.Shoot.Area..pixels. divided by 1000 i.e. the sum of the 3 Shoot Area Views divides by 1000 (calculated)
         \item Area.SV: the mean of Area.SV1 and Area.SV2 divided by a 1000 (calculated)
         \item Area.TV: the Projected.Shoot.Area from Top View divided by 1000 (calculated)
         \item Area.SV1: the Projected.Shoot.Area from Side View 1
         \item Area.SV2: the Projected.Shoot.Area from Side View 2
         \item Image.Biomass: Area.SV x sqrt(Area.TV) (calculated)
         \item Max.Height: maximum of Max.Dist.Above.Horizon.Line.SV1 and Max.Dist.Above.Horizon.Line.SV2 (calculated)
         \item Max.Dist.Above.Horizon.Line.SV1: the Max.Dist.Above.Horizon.Line.SV1 converted to cm (calculated)
         \item Max.Dist.Above.Horizon.Line.SV2: the Max.Dist.Above.Horizon.Line.SV2 converted to cm (calculated)
         \item Density: Area divided by Max.Height (calculated)
         \item Volume: Convex.Hull.TV x Max.Height (calculated)
         \item Centre.Mass: mean of Centre of Mass SV1 and SV2 expressed in cm (calculated)
         \item Center.Of.Mass.Y.SV1: Centre of Mass from Side View 1
         \item Center.Of.Mass.Y.SV2: Centre of Mass from Side View 2
         \item Convex.Hull.SV: mean of Convex.Hull.Area.SV1 divided by 1000 (calculated)
         \item Convex.Hull.TV: Convex.Hull.Area.TV divided by 1000 (calculated)
         \item Convex.Hull.Area.TV: area of Top View Convex Hull
         \item Convex.Hull.Area.SV1: area of Side View 1 Convex Hull
         \item Convex.Hull.Area.SV2: area of Side View 2 Convex Hull
         \item Boundary.To.Area.Ratio.SV1 
         \item Boundary.To.Area.Ratio.SV2 
         \item Boundary.To.Area.Ratio.TV 
         \item Compactness.SV1: Area.SV1 / Convex.Hull.Area.SV1 (calculated)
         \item Compactness.SV2 Area.SV2 / Convex.Hull.Area.SV2 (calculated)
         \item Compactness.SV: Area.SV / Convex.Hull.SV (calculated)
         \item Compactness.TV: Area.TV / Convex.Hull.TV (calculated)
         \item Caliper.Length.SV1 
         \item Caliper.Length.SV2 
         \item Caliper.Length.TV}}
\details{
The columns are copied from data, except for those columns in the list under \bold{Value} that have `(calculated)' appended.
}
\author{Chris Brien}
\examples{
data(exampleData)
longi.prime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1)}
\keyword{data}
\keyword{manip}