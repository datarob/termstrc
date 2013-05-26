\name{objfct_dl_bonds}
\alias{objfct_dl_bonds}

\title{Diebold/Li Loss Function for Bonds}
\description{Calculates the sum of the weighted squared price error.}
\usage{
objfct_dl_bonds(beta, lambda, m, cf, w, p)
}

\arguments{
  \item{beta}{Paramter vector of the Diebold/Li spot rate function
}
  \item{lambda}{Lambda of the Diebold/Li spot rate function
}
  \item{m}{maturity matrx
}
  \item{cf}{cashflow matrix
}
  \item{w}{weights vector
}
  \item{p}{price vector
}
}

