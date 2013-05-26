\name{cSums}
\alias{cSums}

\title{Form Column Sums}
\description{Calculates column sums for numer matrices}
\usage{
cSums(x, na.rm = FALSE, dims = 1L)
}

\arguments{
  \item{x}{ matrix
}
  \item{na.rm}{logical. Should missing values (including NaN) be omitted from the calculations?
}
  \item{dims}{integer. Which dimensions are regarded as rows or columns to sum over
}
}


\note{The function is an optimized version of \code{colSums} and only used for internaly calculations.}

\seealso{ \code{\link{colSums}}}
