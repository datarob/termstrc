\name{grad_dl_bonds}
\alias{grad_dl_bonds}

\title{Diebold/Li Gradient function}
\description{Calculates the gradient of the objective function. The objective function minimizes the sum of the weighted squared price errors. The spot rate function is based on Diebold/Li. }

\usage{
grad_dl_bonds(beta, lambda, m, cf, w, p)
}

\arguments{
  \item{beta}{Spot rate parameter vector}
  \item{lambda}{fixed spot rate parameter
}
  \item{m}{maturity matrix}
  \item{cf}{cashflow matrix
}
  \item{w}{weights vector
}
  \item{p}{price vector 
}
}

\value{returns the gradient vector }


