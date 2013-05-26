\name{objfct_ns}
\Rdversion{1.1}
\alias{objfct_ns}

\title{Nelson/Siegel Loss Function for Yields}
\description{
Calculates the sum of the squared spot rate error.
}
\usage{
objfct_ns(beta, m, y)
}

\arguments{
  \item{beta}{Paramter vector of the Nelson/Siegel spot rate function
  (for details see: \code{\link{spr_ns}}),
}
  \item{m}{maturity vector
}
  \item{y}{observed yield vector
}
}


\seealso{ \code{\link{spotrates}}, \code{\link{spr_ns}}
}
