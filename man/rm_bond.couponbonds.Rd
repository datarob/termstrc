\name{rm_bond.couponbonds}
\alias{rm_bond.couponbonds}

\title{S3 Remove Bond Method
}
\description{
 Specified bonds and their associated data are removed from a static bond data set of the class \command{"couponbonds"}.
}
\usage{
\method{rm_bond}{couponbonds}(bonddata, group, ISIN)
}


\arguments{
  \item{bonddata}{bond data set. }
   \item{group}{the group where the bonds to be removed belong to.}
  \item{ISIN}{the ISIN numbers of the bonds to remove. } 
}


\value{The function returns the new bond data set.}

\examples{
data(govbonds) 
newgovbonds <- rm_bond(govbonds,"GERMANY","DE0001135150")
}