\name{forwardrates}
\Rdversion{1.1}
\alias{forwardrates}


\title{Forward Rate Calculation }

\description{Calculates forward rates according to the Diebold/Li, Nelson/Siegel, Svensson approach.}

\usage{
forwardrates(method, beta, m, lambda)
}



\arguments{
 \item{method}{forward rate function type: \command{"dl"} for Diebold/Li, \command{"ns"} for Nelson/Siegel, \command{"sv"} for Svensson,\command{"asv"} for adjusted Svensson.}
  \item{beta}{ parameter vector \eqn{\bm{\beta}}.}
  \item{m}{maturity or a vector of maturities.}
  \item{lambda}{\eqn{=1/\tau_1}, a scalar; only required for Diebold/Li forward rate function
}
}
\value{The function returns a vector with the calculated forward rates.}
  
\seealso{\command{\link{fwr_dl}}, \command{\link{fwr_ns}}, \command{\link{fwr_sv}}}
\examples{
forwardrates(method="ns",beta=c(0.03,0.02,0.01,5),m=1:30)
}

