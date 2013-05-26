\name{objfct_sv}
\Rdversion{1.1}
\alias{objfct_sv}

\title{Svensson Loss Function for Yields
}
\description{Calculates the sum of the squared spot rate error.
}
\usage{
objfct_sv(beta, m, y)
}

\arguments{
  \item{beta}{Paramter vector of the Svensson spot rate function
  (for details see: \code{\link{spr_sv}}),
}
  \item{m}{maturity vector
}
  \item{y}{observed yield vector
}
}


\seealso{ \code{\link{spotrates}}, \code{\link{spr_sv}}
}