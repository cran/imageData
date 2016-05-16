#Functions to do a Principal Variables Analysis (PVA)
"S22update" <- function(k, i, S22, nvar)
#Function that adds a variable and calucllates the partial correlation matrix
#Called by PVA and PVA.manual
{ #Update S22
  si <- S22[i, i]
  if (k < (nvar - 1))
  { S21 <- S22[i, -i]
    S22 <- S22[,-i] 
    S22 <- S22[-i, ]
  }
  else
  { if (k == (nvar-1))
  { noti <- c(1,2)[!(i == c(1,2))]
    S21 <- S22[i, noti]
    S22 <- S22[noti,noti]
  }
  }
  S22 <- S22 - (S21 %*% t(S21))/si
  return(S22)
}

"PVA" <- function(responses, data, nvarselect = NULL, p.variance = 1, include = NULL, 
                  plot = TRUE, ...)
  #Automatically selects variables using Principla Variables Analysis (PVA)
  #nvarselect is the number of to variables to select, counting those in include.
  #It is possible to use several criteria to control the selection:
  # 1) select all varaibles is increasing order of amount of information provided
  # 2) select nvarselect variables (this will default to the number of variables in R)
  # 3) select just enough variables, up to a maximum of nvarselect variables, to explain at least 
  #    p.variance*100 per cent of the total varaince  
{ #Get correlation matrix
  R <- rcorr(as.matrix(data[responses]))$r
  
  #Initialize
  nvar <- nrow(R)
  if (is.null(nvarselect))
    nvarselect <- nvar
  nvarselect <- min(nvarselect, nvar)
  varnotselect <- rownames(R)
  varselect <- vector(mode = "character", length = nvarselect)
  h.ordered <- vector(mode = "numeric", length = nvarselect)
  p.var <- data.frame(Variable = 1:nvarselect,
                      Cumulative.Propn = rep(0, (nvarselect))) 
  
  #Deal with variables that must be included
  ivarselect <- 0
  pinclude <- 0
  S22 <- R
  if (!(is.null(include)))
  { ivarselect <- length(include)
    nvarselect <- nvarselect - ivarselect
    varselect[1:ivarselect] <- include
    varnotselect <- varnotselect[!(varnotselect %in% include)]
    
    if (ivarselect == 1)
    { 
      S22 <- R[varnotselect, varnotselect] - (matrix(R[varnotselect, include], ncol=ivarselect) %*% 
                                                ginv(matrix(R[include, include], ncol=ivarselect)) %*% 
                                                matrix(R[include, varnotselect], nrow=ivarselect))
      h.ordered[1] <- sum(R[, include]*R[, include])
      p.var$Cumulative.Propn[1:ivarselect] <- 1 - (sum(diag(S22))/nvar)
    }
    
    else
    { for (k in 1:ivarselect)      
    { i <- match(include[k], colnames(S22))
      h.ordered[k] <-  sum((S22*S22)[,i])
      S22 <- S22update(k, i, S22, nvar)
      p.var$Cumulative.Propn[k] <- 1 - (sum(diag(S22))/nvar)
    }
    }
    pinclude <- p.var$Cumulative.Propn[ivarselect]
  }
  
  # If still have not explained enough variance or variables, select some more
  k <- 0
  if (pinclude  <= p.variance & nvarselect > 0)
  { #Adjust for variables that must be included
    kvarselect <- nvarselect
    if (nvarselect+ivarselect >= nvar)
      kvarselect <- nvarselect-1
    for (k in 1:kvarselect)
    { #Select variable with maximum h
      h <- colSums(S22*S22)
      h.ordered[k+ivarselect] <- max(h)
      mh <- match(h.ordered[k+ivarselect], h)
      varselect[k+ivarselect] <- varnotselect[mh]
      varnotselect <- varnotselect[-mh]
      #Update S22
      S22 <- S22update(k+ivarselect, mh, S22, nvar)
      
      #Caclulate proportion of variance
      p.var$Cumulative.Propn[k+ivarselect] <- 1 - (sum(diag(S22))/nvar)
      if (p.var$Cumulative.Propn[k+ivarselect] >= p.variance) 
        break()
    }
    
    #Add last variable if required
    if (nvarselect+ivarselect == nvar & p.var$Cumulative.Propn[k+ivarselect] <= p.variance)
    { varselect[nvar] <- varnotselect[1]
      varnotselect <- varnotselect[-1]
      h.ordered[nvar] <- h[varselect[nvar]]
      p.var$Cumulative.Propn[nvar] <- 1
    } else #crop results if required
    { if (k+ivarselect < length(varselect))
    { varselect <- varselect[1:(k+ivarselect)]
      h.ordered <- h.ordered[1:(k+ivarselect)]
      p.var <- p.var[1:(k+ivarselect),]
    }
    }
  }
  
  #Finalize p.var data.frame
  nvarselect <- nrow(p.var)
  p.var <- within(p.var, 
{ Added.Propn <- Cumulative.Propn
  h.partial <- h.ordered
  Selected <- factor(varselect, levels=varselect)
})
p.var$Added.Propn[2:nvarselect] <- p.var$Added.Propn[2:nvarselect] - p.var$Added.Propn[1:(nvarselect-1)]
p.var <- p.var[c("Variable","Selected", "h.partial", "Added.Propn", "Cumulative.Propn")]

#Plot the variance explained
if (plot)
{ plt <- ggplot(p.var, aes(x=Selected, y=Cumulative.Propn)) + 
    geom_bar(stat="identity") + theme_bw()  +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=20)) +
    scale_y_continuous(breaks=seq(0,1,0.1)) + ylab("Cumulative variance proportion") +
    xlab("Variable selected")
  print(plt)
}

return(p.var)
}

"rcontrib" <- function(responses, data, include = NULL)
  #Allows the manual selection of a set of variables
  #include specifies the set of variables
  #h is returned and from this you can decide which variable(s) to add to the include list
{ R <- Hmisc::rcorr(as.matrix(data[responses]))$r
  
  if (is.null(include))
    h <- colSums(R*R)
  else
  { #Initialize
    nvar <- nrow(R)
    ivarselect <- length(include)
    h.ordered <- vector(mode = "numeric", length = ivarselect)
    p.var <- data.frame(Variable = 1:ivarselect,
                        Cumulative.Propn = rep(0, (ivarselect))) 
    if (ivarselect == 1)
    { varnotselect <- rownames(R)
      varnotselect <- varnotselect[!(varnotselect %in% include)]
      S22 <- R[varnotselect, varnotselect] - (matrix(R[varnotselect, include], ncol=ivarselect) %*% 
                                                ginv(matrix(R[include, include], ncol=ivarselect)) %*% 
                                                matrix(R[include, varnotselect], nrow=ivarselect))
      h.ordered[1] <- sum(R[, include]*R[, include])
      p.var$Cumulative.Propn[1:ivarselect] <- 1 - (sum(diag(S22))/nvar)
    }
    else
    { S22 <- R
      for (k in 1:ivarselect)      
      { i <- match(include[k], colnames(S22))
        h.ordered[k] <-  sum((S22*S22)[,i])
        S22 <- S22update(k, i, S22, nvar)
        p.var$Cumulative.Propn[k] <- 1 - (sum(diag(S22))/nvar)
      }
    }
    h <- colSums(S22*S22)
  }
  
  return(h)
}

"intervalPVA" <- function(responses, data, times.factor = "Days", start.time, end.time, 
                           nvarselect = NULL, p.variance = 1, include = NULL, 
                           plot = TRUE, ...)
  #Call PVA to perform a PVA on all data in a time interval
{ d <- subset(data, data[[times.factor]] %in% as.character(start.time[1]:end.time[1]), 
              select=responses)
  p.var <- PVA(responses, data=d, 
               nvarselect = nvarselect, p.variance = p.variance, include = include, 
               plot=plot, ...)
  return(p.var)
}

#Functions to calculate and plot correlation matrices for a set of responses,
"corrPlot" <- function(responses, data, show.sig = FALSE, title = NULL, 
                       pairs.sets = NULL, labelSize = 4, ...)
{ red <- RColorBrewer::brewer.pal(3, "Set1")[1]
  blu <- RColorBrewer::brewer.pal(3, "Set1")[2]
  corr <- Hmisc::rcorr(as.matrix(data[responses]))
  p <- within(reshape::melt.array(corr$P), 
              { X1 <- factor(X1, levels=dimnames(corr$r)[[1]])
                X2 <- factor(X2, levels=rev(levels(X1)))
              })
names(p)[match("value", names(p))] <- "p"
corr <- within(reshape::melt.array(corr$r), 
               { X1 <- factor(X1, levels=dimnames(corr$r)[[1]])
                 X2 <- factor(X2, levels=rev(levels(X1)))
                })
names(corr)[match("value", names(corr))] <- "r"
corr <- merge(corr, p, by=c("X1", "X2"))
plt <- ggplot(corr, aes(X1, X2)) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradient2(low=red, high=blu, limits=c(-1, 1)) +
  labs(x=NULL, y=NULL, ggtitle=title) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=20),
        axis.text.y=element_text(size=20),
        plot.title=element_text(face="bold")) +
  guides(fill=guide_legend(title="r"))
if (show.sig)
{ corr <- within(corr,
                 sig <- ifelse(p > 0.1, "",
                               ifelse(p > 0.05, paste0(round(r, 2), "(.)"),
                                      ifelse(p > 0.01, paste0(round(r, 2), "(*)"),
                                             ifelse(p > 0.001, paste0(round(r, 2), "(**)"),
                                                    paste0(round(r, 2), "(***)"))))))
  plt <- plt + geom_text(data=corr, aes(label=sig), size=3)
} else
  plt <- plt + geom_text(data=corr, aes(label=round(r, 2)), size=4)
print(plt)
if (is.null(pairs.sets))
  print(GGally::ggpairs(data=data[responses], axisLabels = "internal"))
else
  for(k in 1:length(pairs.sets))
    print(GGally::ggpairs(data=data[responses[pairs.sets[[k]]]], axisLabels = "internal")) 
#Cannot work out how to use ggally_diagAxis to set the label parameters
#                  params=list(labelSize=labelSize, gridlLabelSize=2)))
invisible()
}
