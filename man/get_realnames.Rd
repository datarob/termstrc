\name{get_realnames}
\alias{get_realnames}

\title{Name Conversion
}
\description{
Converts Term Structure Method Into Real Name
}
\usage{
get_realnames(method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{form of the spot rate function, i.e., \command{"ns","sv","asv","dl"} }
}

\value{Returns a character string with the real name}

\examples{
get_realnames("asv")
}

