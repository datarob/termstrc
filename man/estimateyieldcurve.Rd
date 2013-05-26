\name{estimateyieldcurve}
\alias{estimateyieldcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate Zero-coupon Yield Curves
}
\description{
Estimate Zero-coupon Yield curves assuming a certain spot rate function
}
\usage{
estimateyieldcurve(method, y, m, beta, lambda,
 objfct, grad_objfct, constraints, constrOptimOptions)
}

\arguments{
\item{method}{form of the spot rate function}
  \item{y}{yields
}
  \item{m}{maturities
}
  \item{beta}{parameter vector
}
  \item{lambda}{parameter for Diebold/Li}
  \item{objfct}{objective function}
  \item{grad_objfct}{grad_objfct
}
  \item{constraints}{parameter constraints
}
  \item{constrOptimOptions}{solver options
}
}
\details{internal helper function}
