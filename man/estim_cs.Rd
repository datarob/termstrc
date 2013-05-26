\name{estim_cs}
\Rdversion{1.1}
\alias{estim_cs}

\title{Cubic Splines Term Structure Estimation
}
\description{
Function for estimating the term structure of coupon bonds based on
cubic splines.
}
\usage{
estim_cs(bonddata, group, matrange="all", rse=TRUE)
}

\arguments{
 \item{bonddata}{a data set of bonds in list format.}
 \item{group}{vector defining the group of bonds used for the estimation,
 
e.g., \command{c("GERMANY","AUSTRIA")}.}
  \item{matrange}{use \command{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \command{c(lower,upper)}. }

  \item{rse}{ \command{TRUE} (default) calculates robust standard erros for the confidence intervalls of the discount curve}
}

\seealso{ \code{\link{estim_cs.couponbonds}}}

