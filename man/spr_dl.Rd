\name{spr_dl}
\Rdversion{1.1}
\alias{spr_dl}

\title{Spot Rate Function according to the Diebold and Li Version of the Nelson/Siegel Spot Rate Function
}
\description{
This function calculates the spot rates for certain maturity dates and a parameter vector according to Diebold/Li (2006).
}
\usage{
spr_dl(beta, m, lambda)
}


\arguments{
  \item{beta}{a vector of parameters  \eqn{{\bm{\beta}} = \left(\beta_0,\beta_1,\beta_2\right)}.}
  \item{m}{one maturity (or a vector of maturities).}
  \item{lambda}{\eqn{=\frac{1}{\tau_1}}, a scalar}
}


\details{The spot rate according to Diebold/Li for a maturity \eqn{m} is defined as:
\deqn{
  s(m,\bm{\beta},\lambda) = \beta_0 + \beta_1\frac{1-\exp(-m\lambda)}{m\lambda} + \beta_2\left(\frac{1-\exp(-m\lambda)}{m\lambda} - \exp(-m\lambda)\right).}}



\value{Returns a vector consisting of the calculated spot rates.
}
\references{ F.X. Diebold and C. Li: Forecasting the Term Structure of Government Bond Yields. \emph{Journal 
of Econometrics}, \bold{130}:337--364.}

\seealso{
code{\link{spr_ns}}
}


\examples{
spr_dl(c(0.1,0.03,0.01),1:30,0.0609)
}


