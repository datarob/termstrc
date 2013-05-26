\name{bond_yields}
\Rdversion{1.1}
\alias{bond_yields}



\title{Bond Yield Calculation}

\description{Function for the calculation of bond yields.}

\usage{
bond_yields(cashflows, m, searchint = c(-1, 1), tol = 1e-10)
}

\arguments{
  \item{cashflows}{ matrix with the cashflows of the bonds, including the current dirty price.}
  \item{m}{maturity matrix of the bonds}
  \item{searchint}{search interval for root finding.}
  \item{tol}{desired accuracy for function \command{uniroot}.}
}


\value{ The function returns a matrix with the yields of the bonds and the associated maturities. 
}

\seealso{ \command{\link{uniroot}}}

\examples{
data(govbonds)
cf_p <- create_cashflows_matrix(govbonds[[1]],include_price=TRUE)
m_p <- create_maturities_matrix(govbonds[[1]],include_price=TRUE)
bond_yields(cf_p,m_p)
}



