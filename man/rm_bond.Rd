\name{rm_bond}
\alias{rm_bond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bond Removal Function}
\description{
Specified bonds and their associated data are removed from a static or dynamic bond data set
}
\usage{
rm_bond(bonddata, group, ISIN)
}

\arguments{
  \item{bonddata}{bond data set of the class \command{"couponbond"} or \command{"dyncouponbond"}}
  \item{group}{the group where the bonds to be removed belong to.}
  \item{ISIN}{the ISIN numbers of the bonds to remove.
}
}
