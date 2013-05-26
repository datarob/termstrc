\name{grad_asv_bonds}
\alias{grad_asv_bonds}


\title{adjusted Svensson Gradient Function
}
\description{
Calculates the gradient of the objective function. The objective function minimizes the sum of the weighted squared price errors. The spot rate function is based on an adjusted version of Svensson. 
}
\usage{
grad_asv_bonds(beta,m, cf, w, p)
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
