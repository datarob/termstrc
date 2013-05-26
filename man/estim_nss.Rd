\name{estim_nss}
\Rdversion{1.1}
\alias{estim_nss}

\title{Parametric Term Structure Estimation
}
\description{Function for estimating the term structure of coupon bonds
  and yields, with the spot rate function of Diebold/Li, Nelson/Siegel
  or Svensson.
}
\usage{
estim_nss(dataset, ...)
}

\arguments{
  \item{dataset}{object of the class \command{"zeroyields"}, \command{"couponbonds"} or \command{"dyncouponbonds"}
}
  \item{\dots}{further arguments
}
}

\seealso{\code{\link{estim_nss.zeroyields}}, \code{\link{estim_nss.couponbonds}}, \code{\link{estim_nss.dyncouponbonds}}
}

