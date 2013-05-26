\name{rm_bond.dyncouponbonds}
\alias{rm_bond.dyncouponbonds}

\title{S3 Remove Bond Method
}
\description{
 Specified bonds and their associated data are removed from a dynamic bond data set of the class \command{"dyncouponbonds"}.
}
\usage{
\method{rm_bond}{dyncouponbonds}(bonddata, group, ISIN)
}

\arguments{
  \item{bonddata}{bond data set. }
   \item{group}{the group where the bonds to be removed belong to.}
  \item{ISIN}{the ISIN numbers of the bonds to remove. } 
}

\value{The function returns the new bond data set.}


\examples{
data(datadyncouponbonds)
newdynbonddata <- rm_bond(datadyncouponbonds,"GERMANY","DE0001135150")
}