\name{param}
\Rdversion{1.1}
\alias{param}

\title{Term Structure Parameter Extraction}
\description{The function extracts the estimated term structure parameters}
\usage{
param(x, ...)
}
\arguments{
  \item{x}{object of the class \command{"dyntermstrc_yields"} or \command{"dyntermstrc_nss"}}
  \item{\dots}{further arguments}
}

\value{Returns a list of the class \command{"dyntermstrc_param"}}
\details{For the class \command{"dyntermstrc_param"} print, summary and
  plot methods are offered.}

\seealso{ \code{\link{param.dyntermstrc_nss}}, \code{\link{param.dyntermstrc_yields}},
\code{\link{summary.dyntermstrc_param}}, \code{\link{plot.dyntermstrc_param}} }

