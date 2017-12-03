"check.arg.values" <- function(arg.val, options)
  #Function to check that arg.val is one of the allowed values
  #and to return the position of the argument in the set of values
  #that is stored in options
{ kopt <- pmatch(arg.val, options)
if (is.na(kopt))
  stop("Value for argument, ",arg.val,", is not an allowed option")
return(kopt)
}

remove.repeats <- function(x, column = 1, tolerance = 1E-06)
# function to remove repeated values that differ by no more than tolerance  
#If x is two-dimensional then column specifies which column of x is to be tested for repeats
{ 
  ndim <- length(dim(x))
  if (ndim <= 1)
  {
    n <- length(x)
    if (n > 1)
    { 
      repeats <-   c(FALSE, abs(x[2:n] - x[1:(n-1)]) < tolerance)
      x <- x[!repeats]
    }
  } else
  {
    if (ndim == 2)
    {
      repeats <-   c(FALSE, abs(x[2:n, column] - x[1:(n-1), column]) < tolerance)
      x <- x[!repeats,]
      
    } else
    {
      stop("remove.repeats cannot handle objects of more than two dimensions")
    }
  }
  return(x)
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

