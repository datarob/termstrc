\name{spotrates}
\Rdversion{1.1}
\alias{spotrates}

\title{Function for the Calculation of the Spot Rates }

\description{ The function calculates the spot rates for the chosen spot rate function (Diebold/Li, Nelson/Siegel, Svensson) , a provided maturity and parameter vector.

}


\usage{
spotrates(method, beta, m, lambda)
}



\arguments{
 \item{method}{spot rate function type: \command{"dl"} for Diebold/Li, \command{"ns"} for Nelson/Siegel, \command{"sv"} for Svensson, \command{"asv"} for adjusted Svensson.}
  \item{beta}{ parameter vector \eqn{\bm{\beta}}.}
  \item{m}{maturity or a vector of maturities.}
  \item{lambda}{\eqn{=1/\tau_1}, a scalar; only required for Diebold/Li spot rate function
}
}

\value{Returns a vector with the calculated spot rates.}

\seealso{ \code{\link{spr_dl}},  \code{\link{spr_ns}}, \code{\link{spr_sv}}
}

\examples{
spotrates(method="ns",beta=c(0.03,0.02,0.01,5),m=1:30)
  }

