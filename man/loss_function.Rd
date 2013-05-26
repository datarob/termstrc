\name{loss_function}
\Rdversion{1.1}
\alias{loss_function}
\title{Loss Function used for the Term Structure Estimation }
\description{
  The loss function defines the objective function used for the optimisation. In case of term structure estimation
 the objective function is the sum of the weighted squared price errors.
}
\usage{
loss_function(p, phat, omega)
}

\arguments{
  \item{p}{vector of observed prices.}
  \item{phat}{vector of estimated prices.}
  \item{omega}{weights vector, e.g., duration based weights.}
}
