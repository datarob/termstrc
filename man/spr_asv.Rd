\name{spr_asv}
\alias{spr_asv}

\title{Adjusted Svensson Spot rate function
}
\description{This function calculates the spot rates for certain maturity dates and a parameter vector
  according to an adjusted version of Svensson (1994).
}
\usage{
spr_asv(beta, m)
}



\arguments{
  \item{beta}{ a vector of parameters \eqn{\bm{\beta} = (\beta_0, \beta_1,\beta_2,\tau_1, \beta_3, \tau_2)}. }
  \item{m}{ one maturity (or a vector of maturities). }
}


\details{ The adjusted Svensson spot rate function is defined as:


\deqn{s(m,\bm{\beta}) = \beta_0 + \beta_1\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} + \beta_2\left(\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} - \exp(-\frac{m}{\tau_1})\right) + \beta_3\left(\frac{1-\exp(-\frac{m}{\tau_2})}{\frac{m}{\tau_2}} - \exp(-\frac{2m}{\tau_2})\right)}}

\value{
Returns a vector consisting of the calculated spot rates. 
}


\references{
Lars E.O. Svensson (1994):
Estimating and Interpreting Forward Interest Rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}. }

\examples{
 spr_asv(c(0.07,0.03,0.05,2,0.08,7),1:30)
}

