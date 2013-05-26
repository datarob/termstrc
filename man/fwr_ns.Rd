\name{fwr_ns}
\Rdversion{1.1}
\alias{fwr_ns}

\title{Forward Rate Calculation according to Nelson/Siegel.}
\description{Calculate forward rates according to Nelson/Siegel(1987).}
\usage{
fwr_ns(beta, m)
}



\arguments{
  \item{beta}{parameter vector  \eqn{{\bm{\beta}} = \left(\beta_0,\beta_1,\beta_2,\tau_1\right)}. }
  \item{m}{maturity or maturity vector.}
}

\details{The forward rate for a maturity \eqn{m} is calculated using the following relation:

\deqn{f(m,\bm{\beta}) = \beta_0+\beta_1\exp\left(-\frac{m}{\tau_1}\right)+\beta_2\left[\left(\frac{m}{\tau_1}\right)\exp\left(-\frac{m}{\tau_1}\right)\right].}
}
\value{The function returns the calculated forward rate (vector).}


\references{  
Charles R. Nelson and Andrew F. Siegel (1987): 
Parsimonious Modeling of Yield Curves.
\emph{The Journal of Business}, \bold{60(4)}:473--489.}


\seealso{\code{\link{fwr_sv}},\code{\link{fwr_dl}}, \code{\link{forwardrates}}}
\examples{
fwr_ns(beta=c(0.03,0.02,0.01,5),1:30)
}




