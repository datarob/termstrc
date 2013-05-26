\name{grad_ns_bonds}
\alias{grad_ns_bonds}

\title{Nelson/Siegel Gradient Function
}
\description{
Calculates the gradient of the objective function. The objective function minimizes the sum of the weighted squared price errors. The spot rate function is based on Nelson/Siegel. 
}
\usage{
grad_ns_bonds(beta, m, cf, w, p)
}

\arguments{
  \item{beta}{Spot rate parameter vector
}

  \item{m}{maturity matrix
}
  \item{cf}{cashflow matrix
}
  \item{w}{weights vector
} 
  \item{p}{price vector
}
}
 
\value{ returns the gradient vector}