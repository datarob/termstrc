\name{aabse}
\Rdversion{1.1}
\alias{aabse}

\title{Average Absolute Mean Error}

\description{Calculation of the average absolute mean error (AABSE). The AABSE is also called mean absolute error (MAE).}
\usage{
aabse(actual, estimated)
}

\arguments{
  \item{actual}{vector, consisting of the observed values.}
  \item{estimated}{vector, consisting of the estimated values.}
}


\details{
  Calculation of the AABSE according to the formula:
\deqn{\mbox{AABSE}=\frac{1}{m}|\bm{\epsilon}|\bm{\iota},}



whereas \eqn{\bm{\epsilon}} is the vector of the yield or price errors of the bonds and \eqn{\bm{\iota}} is a column vector filled with ones. \eqn{m} is the number of bonds, for which \eqn{\bm{\epsilon}} has been calculated. 

}

\seealso{\code{\link{rmse}}}

