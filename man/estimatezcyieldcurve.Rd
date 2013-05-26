\name{estimatezcyieldcurve}
\Rdversion{1.1}
\alias{estimatezcyieldcurve}

\title{Estimate Zero-coupon Yield Curves
}
\description{
Estimate Zero-coupon Yield curves assuming a certain spot rate function
}
\usage{
estimatezcyieldcurve(method, startparam, lambda, objfct, grad_objfct,
 constraints, constrOptimOptions, m, cf, weights, p)
}

\arguments{
  \item{method}{form of the spot rate function}
  \item{startparam}{start parameter vector}
  \item{lambda}{parameter for Diebold/Li}
  \item{objfct}{objective function, e.g., sum of the weighted squared price errors}
  \item{grad_objfct}{gradient}
  \item{constraints}{constraints for the solver}
  \item{constrOptimOptions}{solver options}
  \item{m}{maturities}
  \item{cf}{cash flows}
  \item{weights}{weights}
  \item{p}{prices}
}
\details{Used as internal helper function 
}



