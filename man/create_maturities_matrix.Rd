\name{create_maturities_matrix}
\Rdversion{1.1}
\alias{create_maturities_matrix}


\title{Maturity Matrix Creation }
\description{
  Creates a matrix of maturities for a specified group of bonds for a static bond data set. The number of rows is the number of cashflows for the bond with the longest maturity.

}

\usage{
create_maturities_matrix(group, include_price = FALSE)
}

\arguments{
  \item{group}{static bond data set for a certain group of bonds.}
  \item{include_price}{if \command{TRUE} the dirty price is included (default: \command{FALSE}).}
}


\value{The maturity matrix is returned. 
}


\seealso{\command{\link{create_cashflows_matrix}}}

\examples{
data(govbonds)
(m <- create_maturities_matrix(govbonds[[1]]))

## maturities matrix with included maturity of the current 
## dirty price, i.e., zero. 
(m_p <- create_maturities_matrix(govbonds[[1]],include_price=TRUE))

}
