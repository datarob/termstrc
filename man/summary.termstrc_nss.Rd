\name{summary.termstrc_nss}
\Rdversion{1.1}
\alias{summary.termstrc_nss}

\title{ S3 Summary Method }
\description{
  S3 summary method for objects of the class \command{"termstrc_nss"}.
}
\usage{
\method{summary}{termstrc_nss}(object,...)
}

\arguments{
  \item{object}{object of the class \command{"termstrc_nss"}.}
  \item{\dots}{other arguments.}
}
\details{
  The summary method for an object of the class \command{"termstrc_nss"} prints the solution of
 the goodness of fit statistics (RMSE,AABSE) of the optimisation. Moreover a convergence information of the used
 optimiser (\command{optim}) is printed. 
}

\seealso{ \command{\link{nlminb}}, \command{\link{plot.termstrc_nss}}, \command{\link{print.termstrc_nss}}, \command{\link{rmse}}, \command{\link{aabse}}}
