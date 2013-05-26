\name{bond_prices}
\Rdversion{1.1}
\alias{bond_prices}

\title{Bond Price Calculation}

\description{Function for the calculation of bond prices according to the chosen approach (Diebold/Li, Nelson/Siegel, Svensson) based on the cashflows and maturities matrix of the bonds.
}





\usage{
bond_prices(method = "ns", beta, m, cf, lambda)
}

\arguments{
  \item{method}{ defines the desired method: \code{"ns"} for the Nelson/Siegel,\code{"dl"} for Diebold/Li, \code{"sv"} for the Svensson approach.}
  \item{beta}{parameter vector, is linked to the chosen approach.}
  \item{m}{maturities matrix, consists of the maturity dates which are appended to the cashflows of the bonds.}
  \item{cf}{ cashflows matrix.}
  \item{lambda}{additional parameter for the \code{"dl"} spot rate function }
  

}

\value{ 
  Returns a list with:
  \item{spot_rates}{spot rates matrix}
  \item{discount_factors}{discount factors matrix}
  \item{bond_prices}{bond prices vector}
}



\seealso{ \code{\link{spotrates}}}

\examples{
data(govbonds)
cf <- create_cashflows_matrix(govbonds[[1]])
m <- create_maturities_matrix(govbonds[[1]])
beta <- c(0.0511,-0.0124,-0.0303,2.5429)
bond_prices(method="ns",beta,m,cf)$bond_prices
}





