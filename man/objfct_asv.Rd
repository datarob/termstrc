\name{objfct_asv}
\alias{objfct_asv}

\title{Adjusted Svensson Loss Function for Yields
}
\description{Calculates the sum of the squared spot rate error.
}
\usage{
objfct_asv(beta, m, y)
}

\arguments{
  \item{beta}{Paramter vector of the adjusted Svensson spot rate function
  (for details see: \code{\link{spr_asv}}),
}
  \item{m}{maturity vector
}
  \item{y}{observed yield vector
}
}


\seealso{ \code{\link{spotrates}}, \code{\link{spr_asv}}
}