\name{fwr_sv}
\Rdversion{1.1}
\alias{fwr_sv}

\title{Forward Rate Calculation according to Svensson (1994).}
\description{
Calculate forward rates according to Svensson (1994). 
}
\usage{
fwr_sv(beta, m)
}


\arguments{
  \item{beta}{parameter vector \eqn{{\bm{\beta}} = \left(\beta_0,\beta_1,\beta_2,\tau_1,\beta_3,\tau_2\right)}.}
  \item{m}{maturity or vector of maturities. }
}
\details{The forward rate for a maturity \eqn{m} is calculated according to the following formula:

\deqn{f(m,\bm{\beta}) = \beta_0+\beta_1\exp\left(-\frac{m}{\tau_1}\right)+\beta_2\left[\left(\frac{m}{\tau_1}\right)\exp\left(-\frac{m}{\tau_1}\right)\right]
+\beta_3\left[\left(\frac{m}{\tau_2}\right)\exp\left(-\frac{m}{\tau_2}\right)\right].}

}
\value{Returns the a vector with the calculated forward rate (vector).}


\references{Lars E.O. Svensson (1994): 
Estimating and Interpreting Forward Interest Rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}.   }

\seealso{\command{\link{fwr_ns}},\command{\link{fwr_dl}} \command{\link{forwardrates}}}
\examples{
fwr_sv(c(0.03,0.02,0.01,5,0.01,10),1:30)
}
