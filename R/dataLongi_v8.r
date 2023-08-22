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

#Function to move camera prefix to a suffix without the characters up to the first _ (eg RBG_) 
#If no _ then moves whole prefix.
#Assumes that prefix ends at first full.stop
"pre2suffix" <- function(name, labsCamerasViews, keepCameraType)
{ 
  prefix <- (strsplit(name, ".", fixed=TRUE))[[1]][1]
  if (any(labsCamerasViews %in% prefix))
  { 
    if (grepl("_", prefix) && !keepCameraType)
    {
      suffix <- (strsplit(prefix, "_"))[[1]]
      if (length(suffix) == 2)
        suffix <- suffix[2]
      else
        suffix <- paste(suffix[2:length(suffix)], collapse = "_")
    }
    else
      suffix <- prefix
    fst <- nchar(prefix)+2
    name <- paste(substring(name, first=fst), suffix, sep=".")
  }
  return(name)
}

"importExcel" <- function(file, sheet="raw data", sep = ",", cartId = "Snapshot.ID.Tag", 
                          imageTimes = "Snapshot.Time.Stamp", 
                          timeAfterStart = "Time.after.Planting..d.", 
                          cameraType = "RGB", keepCameraType = FALSE, 
                          labsCamerasViews = NULL, prefix2suffix = TRUE, 
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
  { 
    raw.dat <- read.csv(file, sep = sep, as.is=TRUE)
    raw.dat[imageTimes] <- as.POSIXct(raw.dat[[imageTimes]], format = timeFormat)
  }
  else if(grepl("xlsx", file))
  {
    raw.dat <- as.data.frame(read_excel(file, sheet=sheet))
    colnames(raw.dat) <- make.names(colnames(raw.dat))
    #raw.dat <- readWorksheetFromFile(file, sheet=sheet)
  }
  else
    stop("File name does not include csv or xlsx")
  ncinput <- ncol(raw.dat)
  if (!(cartId %in% names(raw.dat)))
    stop("cartId not in imported data")
  if (!(imageTimes%in% names(raw.dat)))
    stop("imageTimes not in imported data")
  

  #Rename image columns 
  #Change cameras and views as specified by labsCamerasViews
  vars <- names(raw.dat)
  if (!is.null(names(labsCamerasViews)))
  {
    old.names <- names(labsCamerasViews)
    for (old in old.names)
      vars <- gsub(old, labsCamerasViews[old], vars, fixed = TRUE)
    names(raw.dat) <- vars
  } else #find
  {
    labsCamerasViews <- vars[grep(cameraType, vars, fixed = TRUE)]
    if (length(labsCamerasViews) == 0)
      warning(paste("No imaging variables for a camera of type ", cameraType, " found", sep=""))
    labsCamerasViews <- strsplit(labsCamerasViews, ".", fixed=TRUE)
    labsCamerasViews <- unique(unlist(lapply(labsCamerasViews, 
                                        function(name) {return(name[[1]][1])})))
    names(labsCamerasViews) <- labsCamerasViews
  }
  
  #Move prefix for camera View to become a suffix without the RGB_
  if (prefix2suffix)
  { 
    newvars <- sapply(vars[1:length(vars)], pre2suffix, 
                      labsCamerasViews = labsCamerasViews, keepCameraType = keepCameraType)
    names(newvars) <- NULL
    names(raw.dat)[match(vars, names(raw.dat))] <- newvars
  } else
  {
    if (!keepCameraType)
    {
      vars <- names(raw.dat)
      vars <- gsub(paste(cameraType, "_",sep = ""), "", vars, fixed = TRUE) 
      names(raw.dat) <- vars
    }
  }
  
  #Change day calculation to take away a time origin and truncate to the nearest whole day
  if (!is.null(startTime))
    raw.dat <- calcTimes(raw.dat, imageTimes = imageTimes, timeFormat = timeFormat,
                         intervals = timeAfterStart , startTime = startTime,
                         intervalUnit = "days", timePositions = "Hour")

  #Plot the imaging times if required
  if (imagetimesPlot)
    imagetimesPlot(raw.dat, intervals=timeAfterStart, timePositions = "Hour", 
                   groupVariable = cartId, ...)

  #Check unique for Snapshot.ID.Tag, Time.after.Planting..d.
  combs <- as.vector(table(raw.dat[[cartId]], raw.dat[[timeAfterStart]]))
  if (any(combs != 1))
    warning(paste("There is not a single observation for",  
                  length(combs[combs != 1]), 
                  "combination(s) of",cartId,"and", timeAfterStart))
  
  #Sort data into cartId, Time.after.Planting..d. order and store
  # - may need to be reordered for analysis purposes
  raw.dat <- raw.dat[order(raw.dat[[cartId]], raw.dat[[timeAfterStart]]), ]
  return(raw.dat)
}


#Function to reduce imaging responses to those to be retained, forming longi.prime.dat
"longitudinalPrime" <- function(data, cartId = "Snapshot.ID.Tag", 
                                imageTimes = "Snapshot.Time.Stamp", 
                                timeAfterStart = "Time.after.Planting..d.", 
                                idcolumns = c("Genotype.ID","Treatment.1"),
                                traits = list(all = c("Area", "Boundary.Points.To.Area.Ratio", 
                                                      "Caliper.Length", "Compactness", 
                                                      "Convex.Hull.Area"), 
                                              side = c("Center.Of.Mass.Y", 
                                                       "Max.Dist.Above.Horizon.Line")),
                                labsCamerasViews = list(all = c("SV1", "SV2", "TV"),
                                                      side = c("SV1", "SV2")), 
                                smarthouse.lev = NULL, 
                                calcWaterLoss = TRUE, pixelsPERcm)
{ 
  #Extract variables from data to form data frame of longitudinal data
  posndatevars <- c(cartId,timeAfterStart,
                    "Smarthouse","Lane","Position",imageTimes)
  imagevars <- NULL
  if (is.list(traits) && !is.null(traits))
  {
    if (is.null(labsCamerasViews))
      imagevars <- unlist(traits)
    else
    {
      if (!is.list(labsCamerasViews) || length(traits) != length(labsCamerasViews))
        stop(paste0("When traits is a list then labsCamerasViews must also be a list ",
                    "with the same number of components as traits"))
      if (length(labsCamerasViews) == 1)
        imagevars <- as.vector(outer(traits[[1]], labsCamerasViews[[1]], paste, sep = "."))
      else
        imagevars <- unlist(mapply(FUN =  function(traits, names) {
          if (is.null(names))
             t <- traits
          else
            t <- as.vector(t(outer(traits, names, paste, sep = ".")))
          invisible(t)
        }, 
        traits, labsCamerasViews))
      names(imagevars) <- NULL
    }
  } else
  {
    if (is.character(traits))
    {
      if (is.null(labsCamerasViews))
        imagevars <- unlist(traits)
      else
        imagevars <- as.vector(outer(traits, labsCamerasViews, paste, sep = "."))
    } else
      stop("traits is neither a list nor a character")
  }
  if (calcWaterLoss)
    vars <- c(posndatevars, idcolumns, "Weight.Before","Weight.After","Water.Amount",
              "Projected.Shoot.Area..pixels.", imagevars)
  else
    vars <- c(posndatevars, idcolumns, "Projected.Shoot.Area..pixels.", imagevars)

  #Check that vars are in data
  if (!all(vars %in% names(data)))
    stop(paste("The following variables are not present in data:  ",
               paste(vars[!(vars %in% names(data))], collapse = ", "), sep = ""))

  longi.prime.dat <- data[, vars]

  #Add factors and variates needed in the analysis
  longi.prime.dat <- longi.prime.dat[do.call(order, longi.prime.dat), ]
  if (is.null(smarthouse.lev))
    smarthouse.lev <- unique(longi.prime.dat$Smarthouse)
  longi.prime.dat$Days <- as.numeric(longi.prime.dat[[timeAfterStart]])
  longi.prime.dat <- within(longi.prime.dat, 
                            { 
                              Smarthouse <- factor(Smarthouse, levels=smarthouse.lev)
                              xDays <- Days - mean(unique(Days))
                              xPosn <- Position - mean(unique(Position))
                              Position <- factor(Position, levels=sort(unique(Position)))
                              Area <- Projected.Shoot.Area..pixels./1000
                            })

  facs <- c("Lane", idcolumns, "Days")
  longi.prime.dat[facs] <- as.data.frame(lapply(longi.prime.dat[facs], FUN = factor))
  

#Now derive a Reps factor 
#+
if (all(idcolumns %in% vars))
{
  longi.prime.dat <- within(longi.prime.dat, 
                            { 
                              Reps <- 1
                              trts <- fac.combine(as.list(longi.prime.dat[idcolumns]))
                            })
  for (t in levels(longi.prime.dat$trts))
  { 
    which.indiv <- with(longi.prime.dat, 
                        sort(unique(longi.prime.dat[trts==t, cartId])))
  for (k in 1:length(which.indiv))
    longi.prime.dat[longi.prime.dat$trts == t & 
                      longi.prime.dat$Snapshot.ID.Tag == which.indiv[k], "Reps"] <- k
  }
  longi.prime.dat$Reps <- factor(longi.prime.dat$Reps)
} else 
  longi.prime.dat$Reps <- NA

  #Form responses that can be calculated by row-wise  operations: 
  longi.prime.dat <- calcTimes(longi.prime.dat, imageTimes = imageTimes,
                               timePositions = "Hour")
  kpx.vars <- imagevars[c(grep("Area.", imagevars, fixed = TRUE), 
                          grep("Convex.Hull.Circumference", imagevars, fixed = TRUE))]
  kpx.vars <- kpx.vars[!grepl("Ratio", kpx.vars, fixed = TRUE)]
  longi.prime.dat[kpx.vars] <- longi.prime.dat[kpx.vars]/1000
  
  #'## Calculate Water Use
  #+
  if (calcWaterLoss)
    longi.prime.dat <- within(longi.prime.dat, 
                              { 
                                Water.Loss <- unlist(by(Weight.After, list(Snapshot.ID.Tag), 
                                                        FUN=calcLagged)) - Weight.Before
                              })

  
  out.posndatevars <- c("Smarthouse","Lane","Position","Days",
                        cartId, imageTimes,"xPosn", "Reps", "Hour", "xDays")
  imagevars <- c("Area", imagevars)
  if (calcWaterLoss)
    imagevars <- c("Weight.Before","Weight.After","Water.Amount", "Water.Loss", imagevars)
  out.vars <- c(out.posndatevars, idcolumns, imagevars)
  
  #Re-order rows and response columns
  longi.prime.dat <- longi.prime.dat[order(longi.prime.dat[[cartId]], longi.prime.dat$Days), ]
  longi.prime.dat <- longi.prime.dat[out.vars]
  return(longi.prime.dat)
}

#Function to add design factors for blocked, split plot design
"designFactors" <- function(data, insertName = NULL, designfactorMethod = "LanePosition", 
                            nzones = 6, nlanesperzone = 4, nmainplotsperlane = 11, nsubplotspermain = 2)
{ 
  options <- c("LanePosition","StandardOrder")
  desfactor <- options[check.arg.values(designfactorMethod, options=options)]

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
                           xZones <- xZones - mean(unique(xZones))
                           
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
{ 
  options <- c("AGR","PGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) & length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  vars <- c(INDICES, times.factor, responses)
  times.diffs <- paste(times.factor, "diffs", sep=".")
  if (times.diffs %in% names(data))
    vars <- c(vars, times.diffs)
  if (any(is.na(match(vars, names(data)))))
    stop("One or more of response, INDICES and times.factor are not in data")
  if (any(is.na(data[[times.factor]])))
    warning(paste("Some values of ",times.factor,
                  " are missing, which can result in merge producing a large data.frame", 
                  sep = ""))
  if (any(unlist(lapply(as.list(data[INDICES]), 
                        function(f)
                          any(is.na(f))))))
    warning(paste("Some values of the factors in INDICES are missing, ",
                  "which can result in merge producing a large data.frame", sep = ""))
  
  #Get columns needed for calculation and order for INDICES, then times.factor
  tmp <- data[vars]
  tmp <- tmp[do.call(order, tmp), ]
  
  #Form time differences in a way that every first day is the same Day
  # - setting first time point to missing results in the growth rates also being NA
  if (!(times.diffs %in% names(tmp)))
  { 
    tmp[times.diffs] <- as.numfac(tmp[[times.factor]])
    tmp[times.diffs] <- calcLagged(tmp[[times.diffs]], operation ="-")
    tmp <- split(tmp, f = as.list(tmp[INDICES], simplify=FALSE))
    tmp <- lapply(tmp, 
                  function(tmp, times.diffs)
                  { 
                    if (nrow(tmp) > 0)
                      tmp[[times.diffs]][1] <- NA
                    return(tmp)                   
                  }, 
                  times.diffs = times.diffs)
    tmp <- do.call(rbind, tmp)
    rownames(tmp) <- NULL
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
  #Remove NAs in INDICES and time.factor in tmp
  if (any(is.na(tmp[[times.factor]])))
    tmp <- tmp[!is.na(tmp[[times.factor]]), ]
  if (any(unlist(lapply(as.list(tmp[INDICES]), 
                        function(f)
                          any(is.na(f))))))
    for (f in INDICES)
      tmp <- tmp[!is.na(tmp[f]),]
  tmp <- tmp[,-match(responses,names(tmp))]
  if (times.diffs %in% names(data))
    tmp <- tmp[,-match(times.diffs,names(tmp))]
  data <- merge(data, tmp, by = c(INDICES, times.factor), sort = FALSE, all.x = TRUE)
  data  <- data[do.call(order, data),]
  return(data)
}

#Function to fit splines, including possible boundary correction from Huang (2001)
ncsSpline <- function(vars, correctBoundaries = FALSE, 
                      df = NULL, cv = FALSE,  ...)
{
  if (ncol(vars) != 2)
    stop("Must supply a two-column matrix or data.frame")
  if (!correctBoundaries)
  {
    if (is.null(df))
    { 
      fity <- smooth.spline(vars, all.knots=TRUE, ...)
    } else
    { 
      fity <- smooth.spline(vars, all.knots=TRUE, df=df, ...)
    }
    fit.spline <- list(x = fity$x, 
                       y = fity$y, 
                       lev = fity$lev,
                       lambda = fity$lambda,
                       df = fity$df,
                       uncorrected.fit = fity)  
  } else
  {
    nval <- nrow(vars)
    W <- matrix(NA, nrow = nval, ncol = 4)
    W2 <- matrix(NA, nrow = nval, ncol = 4)
    W3 <- matrix(NA, nrow = nval, ncol = 4)
    x <- vars[,1]
    y <- vars[,2]
    
    # construct the four polynomials
    W[,1] <- x^2/2-x^4/4+x^5/10
    W[,2] <- x^3/6-x^4/6+x^5/20
    W[,3] <- x^4/4-x^5/10
    W[,4] <- -x^4/12+x^5/20
    if (is.null(df))
    {
      lam <- vector(mode = "numeric", length = nval)
      rss <- vector(mode = "numeric", length = nval)
      tr <- vector(mode = "numeric", length = nval)
      gcv <- vector(mode = "numeric", length = nval)
      # use GCV to search for best lambda
      for(i in 1:nval) 
      {
        u <- -7+7.0*(i-1)/(nval-1)
        lam[i] <- 10^u
        slam <- lam[i]
        # get regular fit for y and four polynomials
        fity <- smooth.spline(x,y,all.knots=TRUE,spar=slam, ...)
        W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,spar=slam, ...)$y
        W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,spar=slam, ...)$y
        W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,spar=slam, ...)$y
        W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,spar=slam, ...)$y
        #Apply boundary correction to the fitted spline
        W2 <- W-W2
        h <- solve(t(W2)%*%W2,t(W2)%*%(y-fity$y))
        # get the newfit, rss and calculate the trace of new S
        newfit <- fity$y+W2%*%h
        rss[i] <- sum((newfit-y)^2)
        W3[,1] <- smooth.spline(x,W2[,1],all.knots=TRUE,spar=slam, ...)$y
        W3[,2] <- smooth.spline(x,W2[,2],all.knots=TRUE,spar=slam, ...)$y
        W3[,3] <- smooth.spline(x,W2[,3],all.knots=TRUE,spar=slam, ...)$y
        W3[,4] <- smooth.spline(x,W2[,4],all.knots=TRUE,spar=slam, ...)$y
        W3 <- W2-W3
        K <- solve(t(W2)%*%W2,t(W3))
        tr[i] <- sum(fity$lev)+sum(diag(K%*%W2))
        gcv[i] <- nval*rss[i]/(nval-tr[i])^2
      }
      # get the optimal lambda and apply to y and the four polynomials
      lam.opt<- lam[order(gcv)[1]]
      slam <- lam.opt
      fity <- smooth.spline(x,y,all.knots=TRUE,spar=slam)
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,spar=slam, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,spar=slam, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,spar=slam, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,spar=slam, ...)$y
    } else #df is specified
    {
      # get regular fit for y and four polynomials for specified df
      fity <- smooth.spline(x,y,all.knots=TRUE,df=df, ...)
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,df=df, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,df=df, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,df=df, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,df=df, ...)$y
      # calculate the trace of new S
      W3[,1] <- smooth.spline(x,W2[,1],all.knots=TRUE,df=df, ...)$y
      W3[,2] <- smooth.spline(x,W2[,2],all.knots=TRUE,df=df, ...)$y
      W3[,3] <- smooth.spline(x,W2[,3],all.knots=TRUE,df=df, ...)$y
      W3[,4] <- smooth.spline(x,W2[,4],all.knots=TRUE,df=df, ...)$y
      W3 <- W2-W3
      K <- solve(t(W2)%*%W2,t(W3))
    }
    #Apply boundary correction to the fitted spline
    W2 <- W-W2
    h <- solve(t(W2)%*%W2,t(W2)%*%(y-fity$y))
    newy <- as.vector(fity$y+W2%*%h)
    # get the newfit leverage values of new S
    lev <- as.vector(fity$lev+diag(W2%*%K))
    tr <- sum(lev)
    rss <- sum((y - newy)^2)
    lambda <- nval*rss/(nval-tr)^2
    # Set up list to return
    fit.spline <- list(x = fity$x, 
                       y = newy, 
                       lev = lev,
                       lambda = lambda,
                       df = df,
                       uncorrected.fit = fity)  
  }
  class(fit.spline) <- "ncsSpline"
  return(fit.spline)
}

predict.ncsSpline <- function(object, x, correctBoundaries = FALSE, 
                             df = NULL, cv = FALSE,  ...)
{
  if (!inherits(object, "ncsSpline"))
    stop("Must supply a an object of class ncsSpline")
  fit <- predict(object$uncorrected.fit, x = x)
  
  #Correct boundaries
  if (correctBoundaries)
  {
    nval <- length(x)
    W <- matrix(NA, nrow = nval, ncol = 4)
    W2 <- matrix(NA, nrow = nval, ncol = 4)
    vars <- data.frame(x = object$uncorrected.fit$x, 
                       yin = object$uncorrected.fit$yin)
    vars <- merge(as.data.frame(fit), vars, all.x = TRUE)
    vars$yin[is.na(vars$yin)] <- vars$y[is.na(vars$yin)]
    vars <- vars[order(vars$x), ]
    x <- vars$x

    # construct the four polynomials
    W[,1] <- x^2/2-x^4/4+x^5/10
    W[,2] <- x^3/6-x^4/6+x^5/20
    W[,3] <- x^4/4-x^5/10
    W[,4] <- -x^4/12+x^5/20
    if (is.null(df))
    {
      slam <- object$lambda
      # get regular fit for y and four polynomials
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,spar=slam, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,spar=slam, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,spar=slam, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,spar=slam, ...)$y
    } else #df is specified
    {
      # get regular fit for y and four polynomials for specified df
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,df=df, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,df=df, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,df=df, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,df=df, ...)$y
    }
    #Apply boundary correction to the fitted spline
    W2 <- W-W2
    h <- solve(t(W2)%*%W2,t(W2)%*%(vars$yin-vars$y))
    fit <- list(x = vars$x, 
                y = as.vector(vars$y+W2%*%h))  
  }
  return(fit)
}

#Function to fit a spline using smooth.spline
"fitSpline" <- function(data, response, x, df=NULL, smoothing.scale = "identity", 
                        correctBoundaries = FALSE, 
                        deriv=NULL, suffices.deriv=NULL, RGR=NULL, AGR = NULL, 
                        na.x.action="exclude", na.y.action = "exclude", ...)
{ 
  #check input arguments
  impArgs <- match.call()
  if ("na.rm"%in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  smscales <- c("identity", "logarithmic")
  smscale <- smscales[check.arg.values(smoothing.scale, options=smscales)]
  
  na.x <- na.y <- c("exclude", "omit", "fail")
  na.y <- c(na.y, "allx", "trimx", "ltrimx", "utrimx")
  na.act.x <- na.x[check.arg.values(na.x.action, options=na.x)]
  na.act.y <- na.y[check.arg.values(na.y.action, options=na.y)]  
  
  if (correctBoundaries & !is.null(deriv))
    stop("Unable to estimate derivatives when correctBoundaries = TRUE")
  
  if (!is.null(deriv) & !is.null(suffices.deriv))
    if (length(deriv) != length(suffices.deriv))
      stop("The number of names supplied must equal the number of derivatives specified")
  if (!is.null(RGR)) 
  {
    if (!(1 %in% deriv))
      stop("To form the RGR, 1 must be included in deriv so that the first derivative is obtained")
    else
      kagr <- match(1, deriv)
  }
  
  tmp <- data
  #Transform data, if required
  if (smscale == "logarithmic")
    tmp[[response]] <- log(tmp[[response]])

  #Convert any infinite values to missing
  if (any(is.infinite(tmp[[response]])))
  {
    tmp[[response]][is.infinite(tmp[[response]])] <- NA
    warning("Some infinite values have been converted to missing")
  }
  
  #Set up fit data.frame
  fit.names <- c(x, paste(response,"smooth",sep="."))
  if (!is.null(deriv))
  {
    if (is.null(suffices.deriv))
      fit.names <- c(fit.names, paste(response,".smooth.dv",deriv,sep=""))
    else
      fit.names <- c(fit.names, paste(response,"smooth", suffices.deriv, sep="."))
  }
  #Add RGR if required
  if (!is.null(RGR))
    fit.names <- c(fit.names, paste(response,"smooth",RGR,sep="."))
  fit <- vector(mode = "list", length = length(fit.names))
  names(fit) <- fit.names
  
  #Process missing values
  #tmp will have no missing values and is what is used in the fitting
  #data retains missing values and is used to obtain the returned data.frame
  #x.pred has the x-values for which predictions are required 
  #  - it should include all the x values that are returned by smooth.spline;
  #    it will not include any x values that are missing, but may include
  #    x values for which there are missing y values, depending on the settings 
  #    of na.y.action, and these x values wil not have been supplied to smooth.spline.
  if (nrow(tmp) == 0)
  {
    warning("A response with no data values supplied")
  } else
  {
    #remove any observations with missing x values
    if (na.act.x %in% c("omit", "exclude"))
    {
      tmp <- tmp[!is.na(tmp[[x]]), ]
      if (na.act.x == "omit")
        data <- tmp
    }
    x.pred <- tmp[[x]]
    #Are there any missing response values now
    if (na.act.y != "fail")
    {
      if (na.act.y %in% c("omit", "exclude"))
      {
        x.pred <- tmp[!is.na(tmp[[response]]), ][[x]]
      } else
      {
        if (grepl("trimx", na.act.y, fixed = TRUE))
        {
          tmp <- tmp[order(tmp[[x]]), ]
          y.nonmiss <- c(1:nrow(tmp))[!is.na(tmp[[response]])]
          x.pred <- tmp[[x]]
          if (length(y.nonmiss) > 0)
          {
            if (na.act.y %in% c("trimx", "ltrimx"))
            {
              x.pred <- x.pred[y.nonmiss[1]:length(tmp[[x]])]
              y.nonmiss <- y.nonmiss - y.nonmiss[1] + 1
            }
            if (na.act.y %in% c("trimx", "utrimx"))
              x.pred <- x.pred[1:y.nonmiss[length(y.nonmiss)]]
            y.nonmiss <- diff(y.nonmiss)
            if (any(y.nonmiss >= 3))
              warning("There are runs of 3 or more contiguous y values that are missing")
          } else
            x.pred <- NULL
        }
      }
      tmp <- tmp[!is.na(tmp[[response]]), ]
      if (na.act.y == "omit")
        data <- tmp
    }
  }
  
  #What are the distinct x values
  distinct.xvals <- sort(unique(tmp[[x]]))
  tol <- 1e-06 * IQR(distinct.xvals)
  distinct.xvals <- remove.repeats(distinct.xvals, tolerance = tol)
  if (length(distinct.xvals) < 4) 
  { #< 4 distinct values and so all fitted values are set to NA
    warning(paste("Need at least 4 distinct x values to fit a spline",
                  "- all fitted values set to NA", sep = " "))
    fit <- as.data.frame(matrix(NA, nrow=nrow(data), ncol = length(fit.names)))
    colnames(fit) <- fit.names
    fit[x] <- data[[x]]
  } else
  { 
    #smooth and obtain predictions corresponding to x.pred
    fitcorrectBoundaries <- correctBoundaries
    if (length(distinct.xvals) <= 5)
    {
      warning(paste("Need more than 5 distinct x values to correct the end-points of a spline",
                    "- no corrections made", sep = " "))
      fitcorrectBoundaries <- FALSE
      fit.spline <- with(tmp, 
                         ncsSpline(tmp[c(x, response)], 
                                   correctBoundaries = fitcorrectBoundaries, 
                                   df = df, ...))
    } else
    {
      fit.spline <- with(tmp, 
                         ncsSpline(tmp[c(x, response)], 
                                   correctBoundaries = correctBoundaries, 
                                   df = df, ...))
    }
    x.pred <- remove.repeats(sort(x.pred), tolerance = tol)
    fit <- NULL
    if (length(x.pred) == length(fit.spline$x))
    {
      if (all(abs(x.pred - fit.spline$x) < tol))
        fit <- list(fit.spline$x, fit.spline$y)
    }
    #Need to refit for current x.pred
    if (is.null(fit))
      fit <- predict.ncsSpline(fit.spline, x = x.pred, 
                               correctBoundaries = fitcorrectBoundaries)
    rsmooth <- paste(response,"smooth",sep=".")
    names(fit) <- c(x, rsmooth)
    #backtransform if transformed
    if (smscale == "logarithmic")
      fit[[rsmooth]] <- exp(fit[[rsmooth]])
    
    #get derivatives if required
    if (!correctBoundaries & !is.null(deriv))
    {
      for (d in deriv)
      {
        if (is.null(suffices.deriv))
          rsmooth.dv <- paste(response,".smooth.dv",d,sep="")
        else
        { 
          k <- match(d, deriv)
          rsmooth.dv <- paste(response,"smooth", suffices.deriv[k], sep=".")
        }
        fit[[rsmooth.dv]] <- predict(fit.spline$uncorrected.fit, x = x.pred, deriv=d)$y
      }
      #Add RGR if required
      if (!is.null(RGR) && smscale == "identity")
      { 
        if (is.null(suffices.deriv))
          rsmooth.dv <- paste(response,".smooth.dv",1,sep="")
        else
        { 
          k <- match(1, deriv)
          rsmooth.dv <- paste(response,"smooth", suffices.deriv[k], sep=".")
        }
        if (!(rsmooth.dv %in% names(fit)))
          stop("First derivative not available to calculate RGR")
        fit[[paste(rsmooth,RGR,sep=".")]] <- fit[[rsmooth.dv]]/fit[[rsmooth]]
      }
      #Add AGR if required
      if (!is.null(AGR) && smscale == "logarithmic")
      { 
        if (is.null(suffices.deriv))
          rsmooth.dv <- paste(response,".smooth.dv",1,sep="")
        else
        { 
          k <- match(1, deriv)
          rsmooth.dv <- paste(response,"smooth", suffices.deriv[k], sep=".")
        }
        if (!(rsmooth.dv %in% names(fit)))
          stop("First derivative not available to calculate AGR")
        fit[[paste(rsmooth,AGR,sep=".")]] <- fit[[rsmooth.dv]]*fit[[rsmooth]]
      }
    }
    fit <- as.data.frame(fit)
    #Merge data and fit, preserving x-order in data
    x.ord <- order(data[[x]])
    fit <- merge(data[c(x,response)],fit, all.x = TRUE, sort = FALSE)
    fit <- fit[,-match(response, names(fit))]
    fit <- fit[order(fit[[x]]),]
    fit[x.ord,] <- fit
  }
  rownames(fit) <- NULL
  return(fit)    
}

#Fit splines to smooth the longitudinal trends in the primary responses
#Specify responses to be smoothed and then loop over them
"splitSplines" <- function(data, response, x, INDICES, df = NULL, 
                           smoothing.scale = "identity", 
                           correctBoundaries = FALSE, 
                           deriv = NULL, suffices.deriv=NULL, RGR=NULL, AGR = NULL, sep=".", 
                           na.x.action="exclude", na.y.action = "exclude", ...)
{ 
  impArgs <- match.call()
  if ("na.rm"%in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  smscales <- c("identity", "logarithmic")
  smscale <- smscales[check.arg.values(smoothing.scale, options=smscales)]

  #Split data frame by each combination of the INDICES factors
  old.names <- names(data)
  tmp <- split(data, as.list(data[INDICES]), sep=sep)
  #Fit splines for each combination of the INDICES factors
  tmp <- lapply(tmp, fitSpline, response=response, x = x, df=df, 
                smoothing.scale = smscale, 
                correctBoundaries = correctBoundaries, 
                deriv=deriv, suffices.deriv=suffices.deriv,  RGR=RGR, AGR=AGR, 
                na.x.action=na.x.action, na.y.action=na.y.action, ...)
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
  #Remove any pre-existing smoothed cols in data
  tmp.smooth <- names(tmp)[-match(c(INDICES,x), names(tmp))]
  tmp.smooth <- na.omit(match(tmp.smooth, names(data)))
  if (length(tmp.smooth) > 0)
    data <- data[ ,-tmp.smooth]
  tmp <- tmp[!is.na(tmp[[x]]), ]
  data <- merge(data, tmp, all.x = TRUE, sort=FALSE)
  #Rearrange columns so original column are followed by new columns
  new.names <- names(data)
  new.names <- new.names[-match(old.names, new.names)]
  data <- data[c(old.names, new.names)]
  rownames(data) <- NULL
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

#Function to test if any values in a set of values are anomalous
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


"fac.getinFormula" <- function(formula = NULL, data = NULL, ...)
#Get a list of the factors in a formula
{ 
  if (is.character(formula))
  { 
    formula <- as.formula(paste("~ ",formula, sep=""))
  }
  else
  { 
    if (!is.null(terms))  
      formula <- as.formula(formula)
  }
  if (is.null(terms))
    facs <- NULL
  else
    facs <- rownames(attr(terms(with(data, formula)), which="factor"))
  return(facs)
}

plotDeviations <- function(dat, x, obs, smoothed, devnplots = "none", x.title = NULL, 
                           facet.x = "Treatment.1", facet.y = "Smarthouse", 
                           labeller = NULL, df)  
{
  ggfacet <- list()
  #Set up facet if have any
  if (facet.x != "." | facet.y != ".")
  {
    facet.string <- paste(facet.y,facet.x,sep="~")
    if (is.null(labeller))
      ggfacet <- list(facet_grid(facet.string))
    else
      ggfacet <- list(facet_grid(facet.string, labeller = labeller))
  }
  
  if ("absolute" %in% devnplots)
  {
    dat$deviations <- dat[[obs]] - dat[[smoothed]]
    print(ggplot(data = dat, aes_string(x=x, y = "deviations")) +
            theme_bw() + 
            geom_boxplot() + 
            geom_hline(yintercept=0, colour="blue") +
            ylab(paste("Absolute", obs, "deviations for df", df, sep = " ")) + 
            xlab(x.title) + ggfacet)
  }
  if ("relative" %in% devnplots)
  {
    dat$deviations <- (dat[[obs]] - dat[[smoothed]])/dat[[smoothed]]
    print(ggplot(data = dat, aes_string(x=x, y = "deviations")) +
            theme_bw() + 
            geom_boxplot() + 
            geom_hline(yintercept=0, colour="blue") +
            ylab(paste("Relative", obs, "deviations for df", df, sep = " ")) + 
            xlab(x.title) + ggfacet)
  }
  invisible()
}

"probeDF" <- function(data, response = "Area", xname="xDays", 
                      individuals="Snapshot.ID.Tag", 
                      na.x.action="exclude", na.y.action = "exclude", 
                      df, smoothing.scale = "identity", correctBoundaries = FALSE, 
                      get.rates = TRUE, rates.method="differences", 
                      times.factor = "Days", x = NULL, x.title = NULL, 
                      facet.x = "Treatment.1", facet.y = "Smarthouse", 
                      labeller = NULL, colour = "black", colour.column=NULL, 
                      colour.values=NULL, alpha = 0.1, 
                      which.traits = c("response", "AGR", "RGR"), 
                      which.plots = "smoothedonly",
                      deviations.boxplots = "none", 
                      ggplotFuncs = NULL, ...)
{ 
  #check input arguments
  impArgs <- match.call()
  if ("na.rm"%in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  options <- c("differences","derivative")
  opt <- options[check.arg.values(rates.method, options=options)]
  options <- c("none", "smoothedonly", "bothseparately", "compare", "deviationsboxplot")
  plots <- options[check.arg.values(which.plots, options=options)]
  if (any(c("bothseparately", "compare") %in% plots))
    plotunsmooth <- TRUE
  else
    plotunsmooth <- FALSE
  options <- c("none", "absolute", "relative")
  devnplots <- options[unlist(lapply(deviations.boxplots, check.arg.values, 
                                     options=options))]
  if ("none" %in% devnplots & length(devnplots) > 1)
    devnplots <- "none"
  if (is.null(x.title))
    x.title <- times.factor
  options <- c("response", "AGR", "RGR", "all")
  traits <- options[unlist(lapply(which.traits, check.arg.values, options=options))]
  if ("all" %in% traits)
    traits <- c("response", "AGR", "RGR")
  
  if (any(c("AGR","RGR") %in% traits) & !get.rates)
    stop("get.rates is FALSE but growth-rate plots have been requested")
  else
    if (!any(c("AGR","RGR") %in% traits) & get.rates)
      get.rates <- FALSE
  
  #Form data.frame with just columns needed 
  v <- c(individuals, times.factor, xname, response)
  if (facet.x != ".")
    v <- c(v, fac.getinFormula(facet.x))
  if (facet.y != ".")
    v <- c(v, fac.getinFormula(facet.y))
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
  responses.smooth <- response.smooth 
  if ((plotunsmooth | !("none" %in% devnplots)) & get.rates)
  {
    tmp <- splitContGRdiff(tmp, response, INDICES=individuals,
                           which.rates = c("AGR","RGR"), times.factor=times.factor)
  }
  for (degfree in df)
  { 
    if (opt == "differences")
    { 
      tmp <- splitSplines(tmp, response, x=xname, INDICES = individuals, 
                          df = degfree, smoothing.scale = smoothing.scale, 
                          correctBoundaries = correctBoundaries, 
                          na.x.action = na.x.action, na.y.action = na.y.action)
      if (get.rates)
      { 
        responses.smooth <- c(responses.smooth, 
                              paste(response.smooth, c("AGR","RGR"), sep="."))
        tmp <- splitContGRdiff(tmp, response.smooth, INDICES=individuals,
                               which.rates = c("AGR","RGR"), times.factor=times.factor)
      } 
    } else #derivatives
    { 
      if (get.rates)
      { 
        responses.smooth <- c(responses.smooth, 
                              paste(response.smooth, c("AGR","RGR"), sep="."))
        tmp <- splitSplines(tmp, response, x=xname, INDICES = individuals, deriv=1, 
                            suffices.deriv="AGR", RGR="RGR", df = degfree, 
                            smoothing.scale = smoothing.scale, 
                            na.x.action = na.x.action, na.y.action = na.y.action)
      }
      else
        tmp <- splitSplines(tmp, response, x=xname, INDICES = individuals, 
                            df = degfree, smoothing.scale = smoothing.scale, 
                            correctBoundaries = correctBoundaries, 
                            na.x.action = na.x.action, na.y.action = na.y.action)
    }
    responses.df <- paste(responses.smooth, as.character(degfree), sep=".")
    names(tmp)[match(responses.smooth, names(tmp))] <- responses.df
  }
  
  #Plot some combination of unsmoothed and smoothed response, AGR and RGR
  if (!("none" %in% plots) | !("none" %in% devnplots))
  { 
    responses.tmp <- names(tmp)
    #Plot response
    if (plotunsmooth & "response" %in% traits)
    { 
      pltu <- longiPlot(data = tmp, x=x, response = response, individuals = individuals, 
                        facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                        colour = colour, colour.column = colour.column, 
                        colour.values = colour.values, alpha = alpha, 
                        title="Unsmoothed response", x.title = x.title, y.title = response, 
                        printPlot=FALSE, ...)
      if (!is.null(ggplotFuncs))
        for (f in ggplotFuncs)
          pltu <- pltu + f
        if (!("compare" %in% plots))
          print(pltu)
    }
    if ("compare" %in% plots)
      for (degfree in df)
      { 
        r <- paste(response.smooth, degfree, sep=".")
        plt <- longiPlot(data = tmp, x=x, response = r, individuals = individuals, 
                         facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                         colour = colour, colour.column = colour.column, 
                         colour.values = colour.values, alpha = alpha, 
                         title="Smoothed response", x.title = x.title, y.title = r, 
                         printPlot=FALSE, ...)
        if (!is.null(ggplotFuncs))
          for (f in ggplotFuncs)
            plt <- plt + f
          grid.newpage()
          pushViewport(viewport(layout = grid.layout(1, 2)))
          print(pltu, vp=viewport(layout.pos.row=1, layout.pos.col=1))
          print(plt, vp=viewport(layout.pos.row=1, layout.pos.col=2))
          plotDeviations(dat = tmp, x = times.factor, obs = response, smoothed = r, 
                         devnplots = devnplots, x.title = x.title, 
                         facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                         df = degfree)
      } else
        if ("compare" %in% plots)
          for (degfree in df)
          { 
            r <- paste(response.smooth, degfree, sep=".")
            plt <- longiPlot(data = tmp, x=x, response = r, individuals = individuals, 
                             facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                             colour = colour, colour.column = colour.column, 
                             colour.values = colour.values, alpha = alpha, 
                             title="Smoothed response", x.title = x.title, y.title = r, 
                             printPlot=FALSE, ...)
            if (!is.null(ggplotFuncs))
            {
              for (f in ggplotFuncs)
                plt <- plt + f
            }
            print(plt)
            plotDeviations(dat = tmp, x = times.factor, obs = response, smoothed = r, 
                           devnplots = devnplots, x.title = x.title, 
                           facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                           df = degfree)
          }
    
    #Plot AGRs
    if ("AGR" %in% traits)
    { 
      if (plotunsmooth & "AGR" %in% traits)
      { 
        pltu <- longiPlot(data = tmp, x=x, response = paste(response,"AGR",sep="."), 
                          individuals = individuals, 
                          facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                          colour = colour, colour.column = colour.column, 
                          colour.values = colour.values, alpha = alpha, 
                          title="Unsmoothed AGR by difference", 
                          x.title = x.title, y.title = paste(response,"AGR",sep="."), 
                          printPlot=FALSE, ...)
        if (!is.null(ggplotFuncs))
          for (f in ggplotFuncs)
            pltu <- pltu + f
          if (!("compare" %in% plots))
            print(pltu)
      }
      if ("compare" %in% plots)
        for (degfree in df)
        { 
          r <- paste(response.smooth, "AGR", degfree, sep=".")
          plt <- longiPlot(data = tmp, x=x, response = r, individuals = individuals, 
                           facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                           colour = colour, colour.column = colour.column, 
                           colour.values = colour.values, alpha = alpha, 
                           title="Smoothed AGR", x.title = x.title, y.title = r, 
                           printPlot=FALSE, ...)
          if (!is.null(ggplotFuncs))
            for (f in ggplotFuncs)
              plt <- plt + f
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(1, 2)))
            print(pltu, vp=viewport(layout.pos.row=1, layout.pos.col=1))
            print(plt, vp=viewport(layout.pos.row=1, layout.pos.col=2))
            plotDeviations(dat = tmp, x = times.factor, 
                           obs = paste(response,"AGR",sep="."), smoothed = r, 
                           devnplots = devnplots, x.title = x.title, 
                           facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                           df = degfree)
        } else
          for (degfree in df)
          { 
            r <- paste(response.smooth, "AGR", degfree, sep=".")
            plt <- longiPlot(data = tmp, x=x, response = r, individuals = individuals, 
                             facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                             colour = colour, colour.column = colour.column, 
                             colour.values = colour.values, alpha = alpha, 
                             title=paste("Smoothed AGR by ", opt, sep=""), 
                             x.title = x.title, y.title = r, printPlot=FALSE, ...)
            if (!is.null(ggplotFuncs))
              for (f in ggplotFuncs)
                plt <- plt + f
              print(plt)
              plotDeviations(dat = tmp, x = times.factor, 
                             obs = paste(response,"AGR",sep="."), smoothed = r, 
                             devnplots = devnplots, x.title = x.title, 
                             facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                             df = degfree)
          }
    }
      
    #Plot RGR
    if ("RGR" %in% traits)
    { 
      if (plotunsmooth & "RGR" %in% traits)
      { 
        pltu <- longiPlot(data = tmp, x=x, response = paste(response,"RGR",sep="."), 
                          individuals = individuals, 
                          facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                          colour = colour, colour.column = colour.column, 
                          colour.values = colour.values, alpha = alpha, 
                          title="Unsmoothed RGR by difference", 
                          x.title = x.title, y.title = paste(response,"RGR",sep="."), 
                          printPlot=FALSE, ...)
        if (!is.null(ggplotFuncs))
          for (f in ggplotFuncs)
            pltu <- pltu + f
          if (!("compare" %in% plots))
            print(pltu)
      }
      if ("compare" %in% plots)
        for (degfree in df)
        { 
          r <- paste(response.smooth, "RGR", degfree, sep=".")
          plt <- longiPlot(data = tmp, x=x, response = r, individuals = individuals, 
                           facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                           colour = colour, colour.column = colour.column, 
                           colour.values = colour.values, alpha = alpha, 
                           title="Smoothed RGR", x.title = x.title, y.title = r, 
                           printPlot=FALSE, ...)
          if (!is.null(ggplotFuncs))
            for (f in ggplotFuncs)
              plt <- plt + f
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(1, 2)))
            print(pltu, vp=viewport(layout.pos.row=1, layout.pos.col=1))
            print(plt, vp=viewport(layout.pos.row=1, layout.pos.col=2))
            plotDeviations(dat = tmp, x = times.factor, 
                           obs = paste(response,"RGR",sep="."), smoothed = r, 
                           devnplots = devnplots, x.title = x.title, 
                           facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                           df = degfree)
        } else
          for (degfree in df)
          { 
            r <- paste(response.smooth, "RGR", degfree, sep=".")
            plt <- longiPlot(data = tmp, x=x, response = r, individuals = individuals, 
                             facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                             colour = colour, colour.column = colour.column, 
                             colour.values = colour.values, alpha = alpha, 
                             title=paste("Smoothed RGR by ", opt, sep=""), 
                             x.title = x.title, y.title = r, printPlot=FALSE, ...)
            if (!is.null(ggplotFuncs))
              for (f in ggplotFuncs)
                plt <- plt + f
              print(plt)
              plotDeviations(dat = tmp, x = times.factor, 
                             obs = paste(response,"RGR",sep="."), smoothed = r, 
                             devnplots = devnplots, x.title = x.title, 
                             facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                             df = degfree)
          }
    }
  }
  invisible(tmp)
}
