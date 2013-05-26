\name{grad_ns}
\alias{grad_ns}

\title{Gradient of the Nelson/Siegel Loss Function for Yields
}
\description{
Calculates the gradient of the Nelson/Siegel Loss Function for Yields
}


\usage{
grad_ns(beta, m, y)
}


\arguments{
  \item{beta}{Parameter of the Nelson/Siegel spot rate function (for
  details see \code{\link{spr_ns}}). 
}
  \item{m}{maturity vector
}
  \item{y}{yield vector
}
}


\seealso{ \code{\link{objfct_ns}},\code{\link{spr_ns}}
}