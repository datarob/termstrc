\name{zeroyields}
\alias{zeroyields}

\title{Zeroyields Data Set Generation}
\description{
The function generates a zeroyield data set out of yield, date and
maturities data. 
}
\usage{
zeroyields(maturities, yields, dates)
}

\arguments{
  \item{maturities}{maturities vector of the yields
}
  \item{yields}{yields matrix
}
  \item{dates}{vector of the observations dates in the format "%d.%m.%Y".
}
}

\value{returns a list, which belongs to the class \command{"zeroyields"}. For the class plot, print and summary methods are offered.}
\seealso{ \code{\link{print.zeroyields}},\code{\link{summary.zeroyields}}, \code{\link{plot.zeroyields}}}

