\name{estim_nss.dyncouponbonds}
\Rdversion{1.1}
\alias{estim_nss.dyncouponbonds}

\title{S3 Estim_nss method}
\description{The method performs an iterative term structure estimation procedure
  on a dynamic bond data set of the class
  \command{"dyncouponbonds"}. Available methods are Nelson/Siegel,
  Diebold/Li and (adjusted) Svensson. 

}
\usage{
\method{estim_nss}{dyncouponbonds}(dataset, group, matrange = "all", method = "ns",
 lambda = 0.0609 * 12, tauconstr = NULL, optimtype = "firstglobal", 
constrOptimOptions = list(control = list(maxit = 2000), 
outer.iterations = 200, outer.eps = 1e-04), ...)
}

\arguments{
  \item{dataset}{dynamic bond data set of the class \command{"dyncouponbonds"}.}
  \item{group}{vector defining the group of bonds used for the estimation,

e.g., \command{c("GERMANY","AUSTRIA")}.
}
  \item{matrange}{use \command{"all"} for no restrictions, or restrict
    the maturity range (in years) used for the estimation with \command{c(lower,upper)}.
}
  \item{method}{ \command{"ns"} for Nelson/Siegel (default), \command{"dl"} for
    Diebold/Li, \command{"sv"} for Svensson or \command{"asv"} for adjusted Svensson.
}
  \item{lambda}{parameter on a yearly time scale with fixed value for \command{"dl"} spot rate
    function (default: 0.0609*12)
}
  \item{tauconstr}{
}
  \item{optimtype}{use \command{"firstglobal"} for an inital search for
  globally optimal start parameters or \command{"allglobal"} for a
  search at every iteration. 
}
  \item{constrOptimOptions}{list with solver control parameters
    (default: control=list(), outer.interations=30,
    outer.eps.=1e-04). For further documentation please refer to
    \command{\link{optim}}}
  \item{\dots}{further arguments}
   
}

\details{
The method iteratively applies the method \command{"estim_nss.couponbonds".}
}
\value{ The method returns an object of the class
  \command{"dyntermstrc_nss"}. The object is a list with sublists of the
  class \command{"termstrc_nss"}. 
}

\seealso{ \code{\link{estim_nss.couponbonds}}
}

\examples{## Run: demos(nss_dynamic)} 