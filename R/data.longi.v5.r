globalVariables(c("Snapshot.ID.Tag", "Snapshot.Time.Stamp", "Time.after.Planting..d.", 
                  "Projected.Shoot.Area..pixels.", 
                  "Smarthouse", "Days", "xDays", "xPosn", "xMainPosn", "Area", 
                  "Genotype.ID", "Treatment.1", "Zones", "Lane", "ZLane",
                  "SHZones", "Mainplots", "ZMainplots", "xZones", 
                  "Hour", "Area.SV", "Area.SV1", "Area.SV2", "Area.TV", "Image.Biomass", "Centre.Mass", 
                  "Convex.Hull.SV", "Convex.Hull.TV", "Compactness.SV", "Max.Height", "Density", "Volume",
                  "Center.Of.Mass.Y.SV1", "Center.Of.Mass.Y.SV2", "Convex.Hull.Area.SV1",
                  "Convex.Hull.Area.SV1", "Convex.Hull.Area.TV", "Max.Dist.Above.Horizon.Line.SV1", 
                  "Max.Dist.Above.Horizon.Line.SV2", 
                  "Weight.Before", "Weight.After", "Water.Amount", "r", "Cumulative.Propn"),
                "imageData", add = TRUE)

"check.arg.values" <- function(arg.val, options)
  #Function to check that arg.val is one of the allowed values
  #and to return the position of the argument in the set of values
  #that is stored in options
{ kopt <- pmatch(arg.val, options)
  if (is.na(kopt))
    stop("Value for argument, ",arg.val,", is not an allowed option")
  return(kopt)
}

"ginv" <- function(x, tol = .Machine$double.eps ^ 0.5)
{ 
  # computes Moore-Penrose inverse of a matrix
  if (!is.matrix(x) | length(dim(x)) != 2 )
  {
    if (length(x) == 1)
      if (abs(x) < tol)
        geninv.x <- Inf
      else
        geninv.x <- 1/x
    else
      stop("x must be a matrix")
  }
  else
  {
    svd.x <- svd(x)
    nonzero.x <- (svd.x$d > svd.x$d[1] * tol)
    rank.x <- sum(nonzero.x)
    geninv.x <- matrix(0, dim(x)[1], dim(x)[2])
    if (rank.x)
    { i <- matrix((1:length(nonzero.x))[nonzero.x], rank.x, 2)
    geninv.x[i] <- 1/svd.x$d[nonzero.x]
    if (all(nonzero.x))
      geninv.x <- svd.x$v %*% geninv.x %*% t(svd.x$u)
    else 
      geninv.x <- svd.x$v[, nonzero.x] %*% geninv.x[nonzero.x, nonzero.x] %*% 
      t(svd.x$u[, nonzero.x])
    }
  }
  geninv.x
}

#Function to move RGB_SV1, RB_SV2, RGB_TV prefix to a suffix without the RBG_ 
"pre2suffix" <- function(name)
{ prefix <- (strsplit(name, ".", fixed=TRUE))[[1]][1]
  if (prefix == "RGB_SV1" | prefix == "RGB_SV2" | prefix == "RGB_TV")
  { suffix <- (strsplit(prefix, "_", fixed=TRUE))[[1]][2]
    fst <- nchar(prefix)+2
    name <- paste(substring(name, first=fst), suffix, sep=".")
  }
  return(name)
}

"importExcel" <- function(file, sheet="raw data", cartId = "Snapshot.ID.Tag", 
                          imageTimes = "Snapshot.Time.Stamp", 
                          timeAfterStart = "Time.after.Planting..d.", 
                          prefix2suffix = TRUE, 
                          startTime = NULL, timeFormat = "%Y-%m-%d %H:%M", 
                          imagetimesPlot = TRUE, ...)
{ 
  #Check arguments
  impArgs <- match.call()
  if ("intervals" %in% names(impArgs))
    stop(paste("importExcel assumes that intervals are specified by timeAfterStart; \n", 
               "to have different intervals, call imagetimesPlot separately"))
  if ("timeAfterPlanting"%in% names(impArgs))
    stop("timeAfterPlanting has been deprecated; use timeAfterStart")
  if ("planting.time"%in% names(impArgs))
    stop("planting.time has been deprecated; use startTime")
  
  #Input the raw imaging data
  if (grepl("csv", file))
  { raw.dat <- read.csv(file, as.is=TRUE)
    raw.dat[imageTimes] <- as.POSIXct(raw.dat[[imageTimes]], format = timeFormat)
  }
  else if(grepl("xlsx", file))
    raw.dat <- readWorksheetFromFile(file, sheet=sheet)
  else
    stop("File name does not include csv or xlsx")
  ncinput <- ncol(raw.dat)
  
  #Rename columns to move prefix for camera View to become a suffix without the RGB_
  if (prefix2suffix)
  { vars <- names(raw.dat)
    newvars <- sapply(vars[1:length(vars)], pre2suffix)
    names(newvars) <- NULL
    names(raw.dat)[match(vars, names(raw.dat))] <- newvars
  }
  
  #Change day calculation to take away a time origin and truncate to the nearest whole day
  if (!is.null(startTime))
    raw.dat <- calcTimes(raw.dat, imageTimes = imageTimes, timeFormat = timeFormat,
                         intervals = timeAfterStart , startTime = startTime,
                         intervalUnit = "days", timePositions = "Hour")

  #Plot the imaging times if required
  if (imagetimesPlot)
    imagetimesPlot(raw.dat, intervals=timeAfterStart, timePositions = "Hour", 
                   ...)

  #Check unique for Snapshot.ID.Tag, Time.after.Planting..d.
  combs <- as.vector(table(raw.dat[[cartId]], raw.dat[[timeAfterStart]]))
  if (any(combs != 1))
    warning(paste("There is not a single observation for",  
                  length(combs[combs != 1]), 
                  "combination(s) of",cartId,"and", timeAfterStart))
  
  #Sort data into cartId, Time.after.Planting..d. order and store
  # - may need to be reordered for analysis purposes
  raw.dat <- raw.dat[order(raw.dat[cartId], raw.dat[timeAfterStart]), ]
  return(raw.dat)
}


#Function to reduce imaging responses to those to be retained, forming longi.prime.dat
"longitudinalPrime" <- function(data, idcolumns = c("Genotype.ID","Treatment.1"), 
                                smarthouse.lev = c("SW"), 
                                calcWaterLoss = TRUE, pixelsPERcm = 18)
{ #Extract variables from data to form data frame of longitudinal data
  posndatevars <- c("Snapshot.ID.Tag","Time.after.Planting..d.",
                    "Smarthouse","Lane","Position",
                    "Snapshot.Time.Stamp")
  imagevars <- c("Weight.Before","Weight.After","Water.Amount",
                 "Projected.Shoot.Area..pixels.",
                 "Area.SV1", "Boundary.Points.To.Area.Ratio.SV1", "Caliper.Length.SV1",
                 "Center.Of.Mass.Y.SV1", "Compactness.SV1", "Convex.Hull.Area.SV1", 
                 "Max.Dist.Above.Horizon.Line.SV1", 
                 "Area.SV2", "Boundary.Points.To.Area.Ratio.SV2", "Caliper.Length.SV2", 
                 "Center.Of.Mass.Y.SV2", "Compactness.SV2", "Convex.Hull.Area.SV2", 
                 "Max.Dist.Above.Horizon.Line.SV2", 
                 "Area.TV", "Boundary.Points.To.Area.Ratio.TV",  "Caliper.Length.TV", 
                 "Compactness.TV", "Convex.Hull.Area.TV")
  vars <- c(posndatevars, idcolumns, imagevars)
  
  #Check that vars are in data
  if (any(is.na(match(vars, names(data)))))
  { miss <- vars[is.na(match(vars, names(data)))]
    stop(paste("The following variables are not present in data:  ",
               paste(miss, collapse = ", "), sep = ""))
  }
  
  longi.prime.dat <- data[, vars]

  #Add factors and variates needed in the analysis
  longi.prime.dat <- longi.prime.dat[do.call(order, longi.prime.dat), ]
  longi.prime.dat <- within(longi.prime.dat, 
                            { Smarthouse <- factor(Smarthouse, levels=smarthouse.lev)
                              Days <- as.numeric(Time.after.Planting..d.)
                              xDays <- Days - mean(unique(Days))
                              xPosn <- Position - mean(unique(Position))
                              Area <- Projected.Shoot.Area..pixels./1000
                            })

  longi.prime.dat <- within(longi.prime.dat, {Position <- factor(Position, levels=sort(unique(Position)))}) 
  facs <- c("Lane", idcolumns, "Days")
  longi.prime.dat[facs] <- as.data.frame(lapply(longi.prime.dat[facs], FUN = factor))
  

#Now derive a Reps factor 
#+
if (all(c("Genotype.ID","Treatment.1") %in% posndatevars))
{
  longi.prime.dat <- within(longi.prime.dat, 
                            { Reps <- 1
                            trts <- fac.combine(list(Genotype.ID,Treatment.1))
                            })
  for (t in levels(longi.prime.dat$trts))
  { which.indiv <- with(longi.prime.dat, 
                        sort(unique(longi.prime.dat[trts==t, "Snapshot.ID.Tag"])))
  for (k in 1:length(which.indiv))
    longi.prime.dat[longi.prime.dat$trts == t & 
                      longi.prime.dat$Snapshot.ID.Tag == which.indiv[k], "Reps"] <- k
  }
  longi.prime.dat$Reps <- factor(longi.prime.dat$Reps)
} else 
  longi.prime.dat$Reps <- NA

#Form responses that can be calculated by row-wise  operations: 
longi.prime.dat <- calcTimes(longi.prime.dat, imageTimes = "Snapshot.Time.Stamp",
                              timePositions = "Hour")
longi.prime.dat <- within(longi.prime.dat, 
                          { 
                            Area.SV <- (Area.SV1 + Area.SV2)/1000/2
                            Area.TV <- Area.TV/1000
                            Image.Biomass <- Area.SV* sqrt(Area.TV)
                            Centre.Mass <- ((2157 - Center.Of.Mass.Y.SV1) + 
                                              (2157 - Center.Of.Mass.Y.SV2))/pixelsPERcm/2
                            Convex.Hull.SV <- (Convex.Hull.Area.SV1 + Convex.Hull.Area.SV1)/1000/2
                            Convex.Hull.TV <- Convex.Hull.Area.TV/1000
                            Compactness.SV <- Area.SV / Convex.Hull.SV
                            Max.Height <- pmax(Max.Dist.Above.Horizon.Line.SV1, 
                                               Max.Dist.Above.Horizon.Line.SV2)/pixelsPERcm
                            Density <- Area/Max.Height
                            Density <- ifelse(is.infinite(Density), NA, Density)
                            Volume <- Convex.Hull.TV*Max.Height
                          })

out.posndatevars <- c("Smarthouse","Lane","Position","Days",
                      "Snapshot.ID.Tag", "xPosn", "Reps",
                      "Snapshot.Time.Stamp", "Hour", "xDays")
out.imagevars <- c("Weight.Before","Weight.After","Water.Amount", "Water.Loss", 
                   "Area", "Area.SV", "Area.TV", 
                   "Area.SV1", "Area.SV2", "Image.Biomass", 
                   "Max.Height", "Max.Dist.Above.Horizon.Line.SV1", 
                   "Max.Dist.Above.Horizon.Line.SV2", 
                   "Density", "Volume", "Centre.Mass", 
                   "Center.Of.Mass.Y.SV1", "Center.Of.Mass.Y.SV2", 
                   "Convex.Hull.SV", "Convex.Hull.TV", "Convex.Hull.Area.TV",
                   "Convex.Hull.Area.SV1", "Convex.Hull.Area.SV2", 
                   "Boundary.Points.To.Area.Ratio.SV1", 
                   "Boundary.Points.To.Area.Ratio.SV2", 
                   "Boundary.Points.To.Area.Ratio.TV",  
                   "Compactness.SV1", "Compactness.SV2", 
                   "Compactness.SV", "Compactness.TV", 
                   "Caliper.Length.SV1", "Caliper.Length.SV2", 
                   "Caliper.Length.TV")
out.vars <- c(out.posndatevars, idcolumns, out.imagevars)

#'## Calculate Water Use
#+
if (calcWaterLoss)
  longi.prime.dat <- within(longi.prime.dat, 
                            { Water.Loss <-   unlist(by(Weight.After, list(Snapshot.ID.Tag), 
                                                        FUN=calcLagged)) - Weight.Before
                            })
else
  out.vars <- out.vars[-(match("Water.Loss",out.vars))]  


#Re-order rows and response columns
longi.prime.dat <- with(longi.prime.dat, longi.prime.dat[order(Snapshot.ID.Tag, Time.after.Planting..d.), ])
longi.prime.dat <- longi.prime.dat[out.vars]
return(longi.prime.dat)
}

#Function to add design factors for blocked, split plot design
"designFactors" <- function(data, insertName = NULL, designfactorMethod = "LanePosition", 
                            nzones = 6, nlanesperzone = 4, nmainplotsperlane = 11, nsubplotspermain = 2)
{ options <- c("LanePosition","StandardOrder")
  desfactor <- options[unlist(lapply(designfactorMethod, check.arg.values, options=options))]

  #Extract variables from data
  vars <- names(data)
  required <- c("Smarthouse", "xPosn", "Snapshot.ID.Tag", "xDays")
  if (desfactor == "LanePosition")
    required <- c("Lane", "Position", required)
  if (any(is.na(match(required, vars))))
    stop("Some required columns are not in data")

  n <- nrow(data)
  smarthouse.lev <- levels(data$Smarthouse)
  nshouse <- length(smarthouse.lev)
  ncarts = length(unique(data$Snapshot.ID.Tag))
  nexpcarts <- nshouse * nzones * nlanesperzone * nmainplotsperlane * nsubplotspermain

  if (desfactor == "StandardOrder" )
  { if (ncarts != nexpcarts)
    stop(paste("The number of unique Snapshot.ID.Tag values must be equal to the product of the numbers of\n",
               "      smarthouses, zones, lanes per zone, mainplots per zone and subplots per mainplot"))
  } else
  { if (ncarts != nexpcarts)
    warning(paste("The number of unique Snapshot.ID.Tag values is not equal to the product of the numbers of\n",
                  "      smarthouses, zones, lanes per zone, mainplots per zone and subplots per mainplot"))
  }
  if (n %% ncarts != 0)
    warning("There is not the same number imagings for each cart")

  #Add factors and variates needed in the analysis
  data <- data[do.call(order, data), ]

  #Generate design factors
  if (desfactor == "LanePosition")
  { if (!is.factor(data$Lane))
    {
      levs <- unique(data$Lane)
      levs <- levs[order(levs)]
      data <- cbind(data, 
                    with(data, fac.divide(factor(Lane, levels = levs), 
                                          list(Zones=nzones, ZLane = nlanesperzone))))
    } else
    data <- cbind(data, 
                  with(data, fac.divide(Lane, 
                                        list(Zones=nzones, ZLane = nlanesperzone))))
    if (!is.factor(data$Position))
    {
      levs <- unique(data$Position)
      levs <- levs[order(levs)]
      data <- cbind(data, 
                    with(data, 
                         fac.divide(factor(Position, levels = levs), 
                                    list(Mainplots=nmainplotsperlane, 
                                         Subplots = nsubplotspermain))))
    } else
    data <- cbind(data, 
                  with(data, 
                       fac.divide(Position, 
                                  list(Mainplots=nmainplotsperlane, 
                                       Subplots = nsubplotspermain))))
  } else
    if (desfactor == "StandardOrder")
    { id <- unique(data$Snapshot.ID.Tag)
      data <- merge(data, 
                    data.frame(fac.gen(list(Smarthouse=smarthouse.lev, 
                                            Zones=nzones, ZLane = nlanesperzone, 
                                            Mainplots=nmainplotsperlane, 
                                            Subplots = nsubplotspermain))[,-1],
                               Snapshot.ID.Tag = id), 
                    all.x=TRUE, sort=FALSE)
    } 
  if (nshouse == 1)
  {
    xMain <- with(data, aggregate(xPosn, by=list(Zones, ZLane, Mainplots), mean))
    names(xMain) <- c("Zones", "ZLane", "Mainplots", "xMainPosn") 
    data <- merge(data, xMain, all.x =TRUE, by = c("Zones", "ZLane", "Mainplots"), sort=FALSE)
    
  } else
  {
    xMain <- with(data, aggregate(xPosn, by=list(Smarthouse, Zones, ZLane, Mainplots), mean))
    names(xMain) <- c("Smarthouse", "Zones", "ZLane", "Mainplots", "xMainPosn") 
    data <- merge(data, xMain, all.x =TRUE, sort=FALSE)
  }
    data <- with(data, data[order(Snapshot.ID.Tag, xDays), ])
    data <- within(data, { SHZones <- fac.combine(list(Smarthouse,Zones))
                           ZMainplots <- fac.combine(list(ZLane,Mainplots))
                           xZones <- as.numeric(Zones)
                         })


  facs <- c("Zones","xZones","SHZones","ZLane","ZMainplots",
            "Subplots", "xMainPosn")
  out.vars <- c(vars,facs)
  if (!is.null(insertName))
  { k <- match(insertName, vars)
    if (!is.na(k))
      out.vars <- c(vars[1:k],facs,vars[(k+1):length(vars)])
  }

  #Re-order rows and response columns
  data <- with(data, data[order(Snapshot.ID.Tag, xDays), ])
  data <- data[out.vars]
  return(data)
}

"splitContGRdiff" <- function(data, responses, INDICES, 
                              which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
                              times.factor = "Days")
{ options <- c("AGR","PGR","RGR")

  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) & length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  vars <- c(INDICES, times.factor, responses)
  times.diffs <- paste(times.factor, "diffs", sep=".")
  if (times.diffs %in% names(data))
    vars <- c(vars, times.diffs)
  if (any(is.na(match(vars, names(data)))))
    stop("One or more of response, INDICES and times.factor are not in data")

  #Get columns needed for calculation and order for INDICES, then times.factor
  tmp <- data[vars]
  tmp <- tmp[do.call(order, tmp), ]
  
  #Form time differences in a way that every first day is the same Day
  # - setting first time point to missing results in the growth rates also being NA
  if (!(times.diffs %in% names(tmp)))
  { tmp[times.diffs] <- as.numfac(tmp[[times.factor]])
    tmp[times.diffs] <- calcLagged(tmp[[times.diffs]], operation ="-")
    tmp <- by(tmp, INDICES = as.list(tmp[INDICES], simplify=FALSE), 
               function(tmp, times.diffs)
               { tmp[[times.diffs]][1] <- NA
                 return(tmp)                   
               }, 
               times.diffs = times.diffs)  
    tmp <- do.call(rbind, tmp)
  }

  #Form AGR (because Day.diffs is NA for first day, so will the growth rates)
  if ("AGR" %in% opt)
  { if (is.null(suffices.rates))
      responses.GR <- paste(responses, "AGR", sep=".")
    else
      responses.GR <- paste(responses, suffices.rates[match("AGR",opt)], sep=".")
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                               FUN = AGRdiff, 
                                               time.diffs = tmp[[times.diffs]]))
  }

  #Form PGR (because Day.diffs is NA for first day, so will the growth rates)
  if ("PGR" %in% opt)
  { if (is.null(suffices.rates))
    responses.GR <- paste(responses, "PGR", sep=".")
    else
      responses.GR <- paste(responses, suffices.rates[match("PGR",opt)], sep=".")
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                               FUN = PGR, 
                                               time.diffs = tmp[[times.diffs]]))
  }
  
  #Form RGR (because Day.diffs is NA for first day, so will the growth rates)
  if ("RGR" %in% opt)
  { if (is.null(suffices.rates))
    responses.GR <- paste(responses, "RGR", sep=".")
    else
      responses.GR <- paste(responses, suffices.rates[match("RGR",opt)], sep=".")
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                               FUN = RGRdiff, 
                                               time.diffs = tmp[[times.diffs]]))
  }
  data <- merge(data, tmp, sort = FALSE, all.x = TRUE)
  return(data)
}


#Function to fit a spline using smooth.spline
"fitSpline" <- function(data, response, x, df=NULL, deriv=NULL, suffices.deriv=NULL, 
                        RGR=NULL, na.rm=FALSE)
{ #check input arguments
  if (!is.null(deriv) & !is.null(suffices.deriv))
    if (length(deriv) != length(suffices.deriv))
      stop("The number of names supplied must equal the number of derivatives specified")
  if (!is.null(RGR)) 
    if (!(1 %in% deriv))
      stop("To form the RGR, 1 must be included in deriv so that the first derivative is obtained")
  else
    kagr <- match(1, deriv)
  #Convert any infinite values to missign
  if (any(is.infinite(data[[response]])))
  {
    data[[response]][is.infinite(data[[response]])] <- NA
    warning("Some infinite values have been converted to missing")
  }
  #Remove NAs if na.rm is TRUE
  if (na.rm & (sum(is.na(data[[response]])) < (length(data[[response]]-1))))
    data <- data[!is.na(data[[response]]), ]
  if (nrow(data)==0)
    stop("Zero-length response supplied")
  if (any(is.na(data[response])))
  { #Have some NAs and so all fitted values are set to NA
    warning("Some response values are NA - all fitted values have been set to NA")
    fit <- list(data[[x]], rep(NA, nrow(data)))
    names(fit) <- c(x, paste(response,"smooth",sep="."))
    if (!is.null(deriv))
    { for (d in deriv)
      if (is.null(suffices.deriv))
        fit[[paste(response,".smooth.dv",d,sep="")]] <- rep(NA, nrow(data))
      else
      { k <- match(d, deriv)
        fit[[paste(response,"smooth", suffices.deriv[k], sep=".")]] <- rep(NA, nrow(data))
      }
      #Add RGR if required
      if (!is.null(RGR))
        fit[[paste(response,"smooth",RGR,sep=".")]] <- rep(NA, nrow(data))
    }
  } else
  { #smooth and obtain any derivatives required
    if (is.null(df))
    { fit.spline <- with(data, 
                         smooth.spline(data[c(x, response)], all.knots=TRUE))
    } else
    { fit.spline <- with(data, 
                         smooth.spline(data[c(x, response)], all.knots=TRUE, df=df))
    }
    fit <- list(fit.spline$x, fit.spline$y)
    names(fit) <- c(x, paste(response,"smooth",sep="."))
    if (!is.null(deriv))
    { for (d in deriv)
      if (is.null(suffices.deriv))
        fit[[paste(response,".smooth.dv",d,sep="")]] <- predict(fit.spline, deriv=d)$y
      else
      { k <- match(d, deriv)
        fit[[paste(response,"smooth", suffices.deriv[k], sep=".")]] <- predict(fit.spline, deriv=d)$y
      }
      #Add RGR if required
      if (!is.null(RGR))
      { rsmooth <- paste(response,"smooth",sep=".")
        if (is.null(suffices.deriv))
          fit[[paste(rsmooth,RGR,sep=".")]] <- fit[[paste(rsmooth,".dv",d,sep="")]]/fit[[rsmooth]]
        else
        { k <- match(1, deriv)
          fit[[paste(rsmooth,RGR,sep=".")]] <- fit[[paste(rsmooth, suffices.deriv[k], sep=".")]]/fit[[rsmooth]]
        }
      }
    }
  }
  fit <- as.data.frame(fit)
  return(fit)    
}

#Fit splines to smooth the longitudinal trends in the primary responses
#Specify responses to be smoothed and then loop over them
"splitSplines" <- function(data, response, x, INDICES, df = NULL, deriv = NULL, suffices.deriv=NULL, 
                           RGR=NULL, na.rm = FALSE, sep=".")
{ #Split data frame by each combination of the INDICES factors
  tmp <- split(data, as.list(data[INDICES]), sep=sep)
  #Fit splines for each combination of the INDICES factors
  tmp <- lapply(tmp, fitSpline, response=response, x = x, df=df, 
                deriv=deriv, suffices.deriv=suffices.deriv,  RGR=RGR, na.rm=na.rm)
  tmp <- do.call(rbind, tmp)
  ncols <- ncol(tmp)
  indices <- rownames(tmp)
  indices <- strsplit(indices, split=sep, fixed=TRUE)
  for (fac in 1:length(INDICES))
  { tmp[[INDICES[fac]]] <- unlist(lapply(indices, 
                                         function(x, fac)
                                         { x[fac]}, 
                                         fac))
  if (is.factor(data[[INDICES[fac]]]))
    tmp[[INDICES[fac]]] <- factor(tmp[[INDICES[fac]]])
  else
    if (is.numeric(data[[INDICES[fac]]]))
      tmp[[INDICES[fac]]] <- as.numeric(tmp[[INDICES[fac]]])
  }
  tmp <- tmp[, c((ncols+1):length(tmp),1:ncols)]
  tmp <- na.omit(tmp)
  data <- merge(data, tmp, all.x = TRUE, sort=FALSE)
  return(data)
}

#Functions to do calculations between successive dates 
# - does not assume same number time points for all individuals
#"Replace"  <- function(x, y) {z <- y}
"calcLagged" <- function(x, operation = NULL, lag=1)
  #This function replaces the observations with values calculated  
  # (i) for positive lag, itself and the value lag observations before it, 
  # (ii) for negative lag, itself and the value lag observations after it.
  #operation specifies calculation to be made on the pair of  values 
  #It returns as many values as are in data, the 1st lag values being NA
{ n <- length(x)
  nl <- n-abs(lag)
  if (is.null(operation))
  { if (lag > 0)
      x[(lag+1):n] <- x[1:nl]
     else
     { if (lag < 0)
         x[1:nl] <- x[(abs(lag)+1):n]
     }
  }
  else
  { FUN <- get(operation)
    FUN <- match.fun(FUN)
    if (lag > 0)
      x[(lag+1):n] <- FUN(x[(lag+1):n], x[1:nl])
    else
    { if (lag < 0)
        x[1:nl] <- FUN(x[1:nl], x[(abs(lag)+1):n])
      else
        x[1:n] <- FUN(x[1:n], x[1:n])
    }
  }
  if (lag > 0)
    x[1:lag] <- NA
  else
    x[(nl+1):n] <- NA
  return(x)
}

#Function to test if a any values in a set of values are anomalous
# in being outside specified limits 
"anom" <- function(x, lower=NULL, upper=NULL, na.rm = TRUE)
{ if (is.null(lower))
  { if (is.null(upper))
      stop("Must supply at least a lower or an  upper limit below or above which values are anomalous")
    else
      anom <- any(x > upper, na.rm=na.rm)
  } else
  { if (is.null(upper))
      anom <- any(x < lower, na.rm=na.rm)
    else
      anom <- any(x < lower, na.rm=na.rm) || any(x > upper, na.rm=na.rm)
  }
  return(anom)
}

#Function to calculate the cumulative sum, ignoring the first element if exclude.1st is TRUE
"cumulate" <- function(x, exclude.1st=FALSE)
{ sum <- x
  if (exclude.1st)
    sum[-1] <- cumsum(x[-1])
  else
    sum <- cumsum(x)
  return(sum)
}

#Functions to calculate growth rates between successive imagings
"AGRdiff" <- function(x, time.diffs, lag=1)
{ x.diffs <- calcLagged(x, operation = "-", lag = lag)
  x.diffs <- x.diffs / time.diffs
}
"PGR" <- function(x, time.diffs, lag=1)
{ x.rates <- calcLagged(x, operation = "/", lag = lag)
  x.rates <- x.rates ^ (1/time.diffs)
}

"RGRdiff" <- function(x, time.diffs, lag=1)
{ x.rates <- log(PGR(x, time.diffs, lag = lag))
}

"WUI" <- function(response, water)
{ response.WUI <- ifelse(water != 0, 
                         response / water, 
                         NA)
  return(response.WUI)
}


#Function that produces a longitudinal plot
"longiPlot" <- function(data, x = "xDays+44.5", response = "Area", individuals="Snapshot.ID.Tag", 
                        x.title = "Days", y.title = "Area (1000 pixels)", title = NULL, 
                        facet.x = "Treatment.1", facet.y =   "Smarthouse", labeller = NULL,  
                        colour = "black", colour.column=NULL, colour.values=NULL, alpha = 0.1, 
                        ggplotFuncs = NULL, printPlot = TRUE)
{ 
  data <- data[!is.na(data[response]),]
  longi.plot <- ggplot(data=data, aes_string(x = x, y = response)) +
                theme_bw() +
                theme(panel.grid.major = element_line(colour = "grey60", size = 0.5), 
                      panel.grid.minor = element_line(colour = "grey80", size = 0.5)) +
                xlab(x.title) + ylab(y.title) + ggtitle(title)
  
  #Do facet if have any
  if (facet.x != "." | facet.y != ".")
  {
    facet.string <- paste(facet.y,facet.x,sep="~")
    if (is.null(labeller))
      longi.plot <- longi.plot + facet_grid(facet.string)
    else
      longi.plot <- longi.plot + facet_grid(facet.string, labeller = labeller)
    longi.plot <- longi.plot + theme(strip.text = element_text(size=12, face="bold"),
                                     axis.title = element_text(face="bold"), legend.position="none")
  }
  if (is.null(colour.column))
    longi.plot <- longi.plot + geom_line(aes_string(group=individuals),  
                                         colour=colour, alpha=alpha)
  else
    longi.plot <- longi.plot + geom_line(aes_string(group=individuals, colour=colour.column), 
                                         alpha=alpha)
  if (!(is.null(colour.values)))
    longi.plot <- longi.plot + scale_colour_manual(values = colour.values)
  
  if (!is.null(ggplotFuncs))
    for (f in ggplotFuncs)
      longi.plot <- longi.plot + f
  
  if (printPlot)
    print(longi.plot)
  invisible(longi.plot)
}


#Function that calculates intervals and imageTimes from imageTimes
"calcTimes" <- function(data, imageTimes = NULL, 
                        timeFormat = "%Y-%m-%d %H:%M",
                        intervals = "Time.after.Planting..d.", startTime = NULL, 
                        intervalUnit = "days", timePositions = NULL)
{
  if (!is.null(imageTimes))
  {
    if (!(imageTimes %in% names(data)))
      stop("A column for imageTimes is not present in data")
    if (any(class(data[[imageTimes]])[1] %in% c("character", "factor")))
      data[imageTimes] <- as.POSIXct(data[[imageTimes]], format = timeFormat)
    units <- c("secs", "mins", "hours", "days")
    unit <- units[check.arg.values(intervalUnit, options=units)]
    if (unit == "secs")
    {
      d <- getOption("digits.secs")
      if (d == 0)
        warning(paste("Fractions of sections will not be stored or extracted unless: \n",
                      "(i) option(digits.secs) has been set to the number of decimal places required \n",
                      "and (ii) %OS is used for seconds in timeFormat",
                      sep=""))
    }
    if (!is.null(startTime))
    { 
      startTime <- as.POSIXct(startTime, format = timeFormat)
      data[[intervals]] <- difftime(data[[imageTimes]], startTime, units=intervalUnit)
      data[[intervals]] <- as.numeric(trunc(data[[intervals]], units=intervalUnit))
    }
    if (!is.null(timePositions))
    {
      data[[timePositions]] <- trunc(data[[imageTimes]], units=unit)
      if (unit == "secs")
      {
        data[[timePositions]] <- as.numeric(format(data[[imageTimes]], "%OS"))
        data[[timePositions]] <- data[[timePositions]] - floor(data[[timePositions]])
      }
      else
        data[[timePositions]] <- as.numeric(difftime(data[[imageTimes]], 
                                                     data[[timePositions]], 
                                                     units=units[(match(unit, units) - 1)]))
    }
  }
  return(data)
}


#Function that produces a plot of the imaging times
"imagetimesPlot" <- function(data, intervals = "Time.after.Planting..d.", 
                             timePositions = "Hour", 
                             groupVariable = "Snapshot.ID.Tag", colourVariable = "Lane", 
                             ggplotFuncs = NULL)
{ 
  #Check whether have enough information to do the calculations
  if (!all(c(intervals, timePositions, groupVariable, colourVariable) %in% names(data)))
    stop(paste("At least one of the columns for intervals, timePositions", 
               "groupVariable or colourVariable is not present in data", sep=""))
  if (!(is.numeric(data[[intervals]])))
    data[intervals] <- dae::as.numfac(data[[intervals]])
  if (!(is.numeric(data[[colourVariable]])))
    data[colourVariable] <- dae::as.numfac(data[[colourVariable]])
  
  #Do plot
  start <- min(data[intervals], na.rm=TRUE)
  end <- max(data[intervals], na.rm=TRUE)
  time.plot <- ggplot(data, aes_string(x=intervals, y=timePositions)) +
    geom_line(aes_string(group=groupVariable, colour=colourVariable), alpha=0.05) + 
    scale_colour_gradient(low="grey60", high="grey20") + 
    geom_point(aes_string(group=groupVariable), size=0.5) +
    facet_grid(Smarthouse ~ .) + theme_bw() +
    scale_x_continuous(breaks=seq(start, end, by=2)) +
    ylab("Hour of day")
  
  if (!is.null(ggplotFuncs))
    for (f in ggplotFuncs)
      time.plot <- time.plot + f
  
  print(time.plot)
  invisible(time.plot)
}

"anomPlot" <- function(data, x="xDays+24.16666667", response="Area.smooth.RGR", 
                       individuals="Snapshot.ID.Tag", 
                       breaks=seq(12, 36, by=2), vertical.line=NULL, 
                       groupsFactor=NULL, lower=NULL, upper=NULL, 
                       start.time=NULL, end.time=NULL, times.factor = "Days", 
                       suffix.interval=NULL, 
                       columns.retained=c("Snapshot.ID.Tag", "Smarthouse", "Lane", "Position", 
                                          "Treatment.1", "Genotype.ID"),
                       whichPrint=c("anomalous","innerPlot","outerPlot"), na.rm=TRUE, ...)
{ 
  if (!all(individuals %in% columns.retained))
    stop("The individuals column(s) is (are) not in the columns.retained")
  if (is.null(lower) & is.null(upper))
    stop("Must set at least one of lower and upper")
  options <- c("anomalous","innerPlot","outerPlot")
  opt <- options[unlist(lapply(whichPrint, check.arg.values, options=options))]
  
  #Determine anomalous individuals
  if (is.null(groupsFactor))
  { 
    anomalous.individuals <- intervalValueCalculate(response=response, FUN = "anom", data=data,
                                                    individuals=individuals, 
                                                    lower=lower, upper=upper, 
                                                    start.time=start.time, end.time=end.time, 
                                                    times.factor=times.factor, 
                                                    suffix.interval=suffix.interval, 
                                                    na.rm=na.rm)
    data <- merge(data, anomalous.individuals[,1:2], by=individuals, sort=FALSE)
  }
  else
  { 
    tmp <- split(data, data[[groupsFactor]])
    ngrps <- length(tmp)
    nstart <- length(start.time)
    nend <- length(end.time)
    nlow <- length(lower)
    nup <- length(upper)
    if (nstart == 0)
    { 
      kstart <- NULL
      if (nend != 1 & nend != ngrps)
        stop("Number of end.time values must be equal 1 or the number of levels in groupsFactor")
      kend <- end.time[1]
    } else
      if (nend ==0)
      { 
        kend <- NULL
        if (nstart != 1 & nstart != ngrps)
          stop("Number of start.time values must be equal 1 or the number of levels in groupsFactor")
        kstart <- start.time[1]
      } else
      {     
        if (nstart != nend | (nstart != ngrps & nstart != 1))
          stop("Number of start.time and end.time values must be equal and equal to 1 \n",
               "or the number of levels in groupsFactor")
        kstart <- start.time[1]
        kend <- end.time[1]
      }
    if (!(nlow == 0 | nlow  == ngrps |  nlow == 1))
      stop("Number of lower values must equal to 1 or the number of levels in groupsFactor")
    if (!(nup == 0 | nup  == ngrps | nup == 1))
      stop("Number of upper values must equal to 1 or the number of levels in groupsFactor")
    klow <- lower[1]
    kup <- upper[1]
    for (k in 1:ngrps)
    { 
      if (nstart > 1)
        kstart <- start.time[k]
      if (nend > 1)
        kend <- end.time[k]
      if (nlow > 1)
        klow <- lower[k]
      if (nup > 1)
        kup <- upper[k]
      anomalous.individuals <- intervalValueCalculate(response=response, FUN = "anom", data=tmp[[k]],
                                                      individuals=individuals, 
                                                      lower=klow, upper=kup, 
                                                      start.time=kstart, end.time=kend, 
                                                      times.factor=times.factor, 
                                                      suffix.interval=suffix.interval, 
                                                      na.rm=na.rm)
      tmp[[k]] <- merge(tmp[[k]], anomalous.individuals[,1:2], by=individuals, sort=FALSE)
    }
    data <- do.call(rbind, tmp)
  }
  response.anom <- names(anomalous.individuals)[2]
  
  #Plot without anomalous individuals
  if (sum(!data[[response.anom]] > 0))
  { 
    innerPlot <- longiPlot(data = subset(data, !data[[response.anom]]), 
                           x=x, response = response, 
                           printPlot=FALSE, ...)
    innerPlot <- innerPlot + scale_x_continuous(breaks=breaks)
    if (!is.null(vertical.line))
      innerPlot <- innerPlot + geom_vline(xintercept=vertical.line, linetype="longdash", size=1)
    if ("innerPlot" %in% opt)
      print(innerPlot)
  } else
    innerPlot <- NULL
  
  #Print out anomalous individuals
  if ("anomalous" %in% opt)
  { 
    anom.dat <- data[c(columns.retained, response.anom)] 
    anom.dat <- split(anom.dat, anom.dat[[individuals]])
    anom.dat <- lapply(anom.dat, 
                       function(dat)
                         dat <- dat[1,])
    anom.dat <- do.call(rbind, anom.dat)
    anom.dat <- anom.dat[anom.dat[[response.anom]],]
    anom.dat <- anom.dat[order(anom.dat[[individuals]]), columns.retained]
    print(anom.dat)
  }  
  
  #Plot anomalous individuals, adding Snapshot.ID.Tag
  if (sum(data[[response.anom]] > 0))
  { 
    outerPlot <- longiPlot(data = subset(data, data[[response.anom]]), 
                           x=x, response = response, alpha=0.5, colour="purple", 
                           printPlot=FALSE, ...)
    outerPlot <- outerPlot + scale_x_continuous(breaks=breaks)
    
    if (!is.null(vertical.line))
      outerPlot <- outerPlot + geom_vline(xintercept=vertical.line, linetype="longdash", size=1)
    
    if ("outerPlot" %in% opt)
      print(outerPlot)
  } else
    outerPlot <- NULL
  
  invisible(list(data = data, innerPlot = innerPlot, outerPlot = outerPlot))
}

"probeDF" <- function(data, response = "Area", xname="xDays", individuals="Snapshot.ID.Tag", 
                      na.rm = FALSE, df, get.rates = TRUE, rates.method="differences", 
                      times.factor = "Days", x = NULL, 
                      facet.x = "Treatment.1", facet.y = "Smarthouse", 
                      which.plots = c("smoothed", "AGR", "RGR"), ggplotFuncs = NULL, ...)
{ 
  options <- c("differences","derivative")
  opt <- options[check.arg.values(rates.method, options=options)]
  options <- c("none", "responseComparison", "unsmoothed", "smoothed", "AGR", "RGR", "all")
  plots <- options[unlist(lapply(which.plots, check.arg.values, options=options))]
  if ("all" %in% plots)
    plots <- c("unsmoothed", "smoothed", "AGR", "RGR")
  if (any(c("AGR","RGR") %in% plots) & !get.rates)
    stop("get.rates is FALSE but growth-rate plots have been requested")
  
  #Form data.frame with just columns needed 
  v <- c(individuals, times.factor, xname, response)
  if (facet.x != ".")
    v <- c(v, facet.x)
  if (facet.y != ".")
    v <- c(v, facet.y)
  if (is.null(x))
    x <- xname
  #  else
  #    v <- c(v,x)
  if (!all(v %in% names(data)))
    stop(paste("Do not have the following required columns in data: ", 
               paste(v[!(v %in% names(data))],sep=", "), "\n", sep=""))
  tmp <- data[v]
  
  #Smooth response and form growth rates
  response.smooth <- paste(response, "smooth", sep=".")
  responses <- response.smooth 
  for (degfree in df)
  { 
    if (opt == "differences")
    { 
      tmp <- splitSplines(tmp, response, x=xname, INDICES = individuals, df = degfree, 
                          na.rm = na.rm)
      if (get.rates)
      { 
        responses <- c(responses, paste(response.smooth, c("AGR","RGR"), sep="."))
        tmp <- splitContGRdiff(tmp, response.smooth, INDICES=individuals,
                               which.rates = c("AGR","RGR"), times.factor=times.factor)
      } 
    } else #derivatives
    { 
      if (get.rates)
      { 
        responses <- c(responses, paste(response.smooth, c("AGR","RGR"), sep="."))
        tmp <- splitSplines(tmp, response, x=xname, INDICES = individuals, deriv=1, 
                            suffices.deriv="AGR", RGR="RGR", df = degfree, na.rm = na.rm)
      }
      else
        tmp <- splitSplines(tmp, response, x=xname, INDICES = individuals, df = degfree,
                            na.rm = na.rm)
    }
    responses.df <- paste(responses, as.character(degfree), sep=".")
    names(tmp)[match(responses, names(tmp))] <- responses.df
  }
  
  #Plot some combination of unsmoothed and smoothed response, AGR and RGR
  if (!("none" %in% plots))
  { 
    responses.tmp <- names(tmp)
    if (any(c("responseComparison", "unsmoothed") %in% plots))
    { 
      pltu <- longiPlot(data = tmp, x=x, response = response, 
                        facet.x=facet.x, facet.y=facet.y, 
                        title="Unsmoothed data", y.title = response, 
                        printPlot=FALSE, ...)
      if (!is.null(ggplotFuncs))
        for (f in ggplotFuncs)
          pltu <- pltu + f
        if (!("responseComparison" %in% plots))
          print(pltu)
    }
    if ("responseComparison" %in% plots)
      for (r in paste(response.smooth, df, sep="."))
      { 
        plt <- longiPlot(data = tmp, x=x, response = r, 
                         facet.x=facet.x, facet.y=facet.y, 
                         title="Smoothed response", y.title = r, 
                         printPlot=FALSE, ...)
        if (!is.null(ggplotFuncs))
          for (f in ggplotFuncs)
            plt <- plt + f
          grid.newpage()
          pushViewport(viewport(layout = grid.layout(1, 2)))
          print(pltu, vp=viewport(layout.pos.row=1, layout.pos.col=1))
          print(plt, vp=viewport(layout.pos.row=1, layout.pos.col=2))
      }
    else
      if ("smoothed" %in% plots)
        for (r in paste(response.smooth, df, sep="."))
        { 
          plt <- longiPlot(data = tmp, x=x, response = r, 
                           facet.x=facet.x, facet.y=facet.y, 
                           title="Smoothed response", y.title = r, 
                           printPlot=FALSE, ...)
          if (!is.null(ggplotFuncs))
            for (f in ggplotFuncs)
              plt <- plt + f
            print(plt)
        }
    if ("AGR" %in% plots)
    { 
      responses.plt <- responses.tmp[grep("AGR", responses.tmp,fixed=T)]
      for (r in responses.plt)
      { 
        plt <- longiPlot(data = tmp, x=x, response = r, 
                         facet.x=facet.x, facet.y=facet.y, 
                         title=paste("Growth rates by ", opt, sep=""), 
                         y.title = r, printPlot=FALSE, ...)
        if (!is.null(ggplotFuncs))
          for (f in ggplotFuncs)
            plt <- plt + f
          print(plt)
      }
    }
    if ("RGR" %in% plots)
    { 
      responses.plt <- responses.tmp[grep("RGR", responses.tmp,fixed=T)]
      for (r in responses.plt)
      {
        plt <- longiPlot(data = tmp, x=x, response = r, 
                         facet.x=facet.x, facet.y=facet.y, 
                         title=paste("Growth rates by ", opt, sep=""), 
                         y.title = r, printPlot=FALSE, ...)
        if (!is.null(ggplotFuncs))
          for (f in ggplotFuncs)
            plt <- plt + f
          print(plt)
      }
    }
  }
  invisible(tmp)
}
