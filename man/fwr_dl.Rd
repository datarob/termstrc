\name{fwr_dl}
\Rdversion{1.1}
\alias{fwr_dl}

\title{Forward Rate Calculation according to Diebold/Li.}
\description{Calculate forward rates according to Diebold/Li(2006).
}
\usage{
fwr_dl(beta, m, lambda)
}

\arguments{
  \item{beta}{parameter vector  \eqn{{\bm{\beta}} = \left(\beta_0,\beta_1,\beta_2\right)}. }
  \item{m}{maturity or maturity vector.}
  \item{lambda}{\eqn{=\frac{1}{\tau_1}}, a scalar}
}

\details{The forward rate for a maturity \eqn{m} is calculated according to the following formula:

\deqn{f(m,\bm{\beta},\lambda) = \beta_0+\beta_1\exp\left(-m\lambda\right)+\beta_2\left[\left(m\lambda\right)\exp\left(-m\lambda\right)\right].}
}
\value{The function returns the calculated forward rate (vector).}


\references{ F.X. Diebold and C. Li: Forecasting the Term Structure of Government Bond Yields. \emph{Journal 
of Econometrics}, \bold{130}:337--364.}


\seealso{\code{\link{fwr_sv}},\code{\link{fwr_ns}}, \code{\link{forwardrates}}}
\examples{fwr_dl(beta=c(0.03,0.02,0.01),1:30,lambda=1/5)}

