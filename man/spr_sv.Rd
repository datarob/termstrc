\name{spr_sv}
\Rdversion{1.1}
\alias{spr_sv}

\title{Spot Rate Function according to Svensson }
\description{This function calculates the spot rates for certain maturity dates and a parameter vector
  according to Svensson (1994).}
\usage{
spr_sv(beta, m)
}


\arguments{
  \item{beta}{ a vector of parameters \eqn{\bm{\beta} = (\beta_0, \beta_1,\beta_2,\tau_1, \beta_3, \tau_2)}. }
  \item{m}{ one maturity (or a vector of maturities). }
}


\details{ The spot rate according to Svensson for a maturity \eqn{m} is calculated using the following function:


\deqn{s(m,\bm{\beta}) = \beta_0 + \beta_1\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} + \beta_2\left(\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} - \exp(-\frac{m}{\tau_1})\right) + \beta_3\left(\frac{1-\exp(-\frac{m}{\tau_2})}{\frac{m}{\tau_2}} - \exp(-\frac{m}{\tau_2})\right)}}

\value{
Returns a vector consisting of the calculated spot rates. 
}


\references{
Lars E.O. Svensson (1994):
Estimating and Interpreting Forward Interest Rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}. }

\examples{
 spr_sv(c(0.07,0.3,0.05,2,0.08,7),1:30)
}



