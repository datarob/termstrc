\name{grad_asv_bonds_grid}
\alias{grad_asv_bonds_grid}

\title{adjusted Svensson Gradient Function for the Grid Search
}
\description{
Calculates the gradient of the objective function for the grid search. The objective function minimizes the sum of the weighted squared price errors. The spot rate function is based on an adjusted version of Svensson. 
}
\usage{
grad_asv_bonds_grid(beta, tau, m, cf, w, p)
}

\arguments{
  \item{beta}{Spot rate parameter vector
}
  \item{tau}{fixed parameters
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
