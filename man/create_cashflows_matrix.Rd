\name{create_cashflows_matrix}
\Rdversion{1.1}
\alias{create_cashflows_matrix}


\title{Cashflows Matrix Creation}

\description{
  Creates a matrix of cashflows for a specified group of bonds for a static bond data set. The number of rows is the number of cashflows for the bond with the longest maturity.
}

\usage{
create_cashflows_matrix(group, include_price = FALSE)
}



\arguments{

  \item{group}{static bond data set for a certain group of bonds.}
  \item{include_price}{if \command{TRUE} the dirty price is included (default: \command{FALSE}).}

}





\value{Returns a matrix which consits of the calculated cashflows.}


\seealso{\command{\link{create_maturities_matrix}} }

\examples{
data(govbonds)
(cf <- create_cashflows_matrix(govbonds[[1]]))

## cf matrix with included current dirty price 
(cf_p <- create_cashflows_matrix(govbonds[[1]],include_price=TRUE))
}


