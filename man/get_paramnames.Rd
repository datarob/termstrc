\name{get_paramnames}
\alias{get_paramnames}

\title{Parameter Names
}
\description{Parameter Names for term structure estimation methods
}
\usage{
get_paramnames(method)
}
\arguments{
  \item{method}{form of the spot rate function,i.e., one of the
    following \command{"ns","sv","asv","dl"}
}
}

\value{Returns a character string with the names of the elements of the
  parameter vector. 
}

\examples{get_paramnames("ns")
}
