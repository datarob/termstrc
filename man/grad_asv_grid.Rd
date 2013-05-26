\name{grad_asv_grid}
\alias{grad_asv_grid}


\title{Adjusted Svensson Gradient Function for the Grid Search
}
\description{
Calculates the gradient of the objective function for the grid search. The objective function minimizes the sum of the squared yield errors. The spot rate function is based on the adjusted version of Svensson. 
}


\usage{
grad_asv_grid(beta, tau, m, y)
}

\arguments{
  \item{beta}{Spot rate function parameter vector
}
  \item{tau}{fixed parameters
}
  \item{m}{maturity vector
}
  \item{y}{yield vector 
}
}

