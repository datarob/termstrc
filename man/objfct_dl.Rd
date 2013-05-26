\name{objfct_dl}
\alias{objfct_dl}

\title{Diebold/Li Loss Function for Yields}
\description{
Calculates the sum of the squared spot rate error.
}

\usage{
objfct_dl(beta, lambda, m, y)
}

\arguments{
\item{beta}{Paramter vector of the Diebold/Li spot rate function
  (for details see: \code{\link{spr_dl}}),
}
  \item{lambda}{Fixed spot rate function parameter
}
  \item{m}{maturity vector
}
  \item{y}{observed yield vector
}
}