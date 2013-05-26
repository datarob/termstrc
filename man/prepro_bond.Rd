\name{prepro_bond}
\Rdversion{1.1}
\alias{prepro_bond}

\title{Bonddata preprocess function}
\description{
Preprocessing a static coupon bond data set, i.e., calculation of cashflows, maturities matrices, price, accrued interest vectors, yield-to-maturity and duration matrices. 
}
\usage{
prepro_bond(group, bonddata, matrange = "all")
}
\arguments{
  \item{group}{character, specifies group of a bond data set.}
  \item{bonddata}{Static bond data set.}
  \item{matrange}{bond data set is filtered according to chosen maturity spectrum \code{c(min,max)}.}
}

\value{
 \item{n_group}{group length}
 \item{sgroup}{sequence of the group length}
  \item{cf}{list with cashflows matrices}
 \item{cf_p}{list with cashflows matrices including the current dirty prices}
  \item{m}{list with maturites matrices}
  \item{m_p}{list with cashflows matrices including the maturities of the current dirty prices}
  \item{p}{list with the dirty price vectors}
  \item{ac}{list with the accrued interest vectors}
  \item{y}{list with the yield-to-maturity matrices}
  \item{duration}{list with the duration, duration based weights matrices}
  \item{timestamp}{date of the data}
}

\examples{
data(govbonds)
bdata <- prepro_bond("GERMANY",govbonds,c(0,10))
## print maturites matrix
bdata$m
}
