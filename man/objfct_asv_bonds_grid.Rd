\name{objfct_asv_bonds_grid}
\alias{objfct_asv_bonds_grid}

\title{Adjusted Svensson Grid Loss Function for Bonds}
\description{Calculates the sum of the weighted squared price error.}
\usage{
objfct_asv_bonds_grid(beta, tau, m, cf, w, p)
}

\arguments{
  \item{beta}{Beta parameters of adjusted Svensson spot rate function 
}
  \item{tau}{Tau parameters of adjusted Svensson spot rate function
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

