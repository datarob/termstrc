\name{rmse}
\Rdversion{1.1}
\alias{rmse}
\title{Root Mean Squared Error}
\description{Calulates the root mean squared error (RMSE).}

\usage{
rmse(actual, estimated)
}

\arguments{
  \item{actual}{vector, consisting of the observed values.}
  \item{estimated}{vector, consisting of the estimated values.}
}
\details{
  Calculation of the RMSE according to the formula:
\deqn{ \mbox{RMSE}=\sqrt{\frac{1}{m}\bm{\epsilon}^2\bm{\iota}}, }
whereas \eqn{\bm{\epsilon}} is the vector of the yield or price errors of the bonds and \eqn{\bm{\iota}} is a column vector filled with ones. \eqn{m} is the number of bonds, for which \eqn{\bm{\epsilon}} has been calculated. 
}

\seealso{\command{\link{aabse}}}