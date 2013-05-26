\name{summary.termstrc_cs}
\Rdversion{1.1}
\alias{summary.termstrc_cs}

\title{S3 Summary Method for Termstrc_cs }
\description{
  S3 summary method for objects of the class \command{"termstrc_cs"}.
}
\usage{
\method{summary}{termstrc_cs}(object,...)
}

\arguments{
  \item{object}{object of the class \command{"termstrc_cs"}.}
  \item{\dots}{other arguments.}
}
\details{
 The summary method for an object of the class \command{"termstrc_cs"} calculates goodness of fit statistics (RMSE, AABSE) of the price and yield errors. Additionally, summary statistics of the regression analysis of the parameters are printed. 
}

\seealso{\command{\link{plot.termstrc_cs}}, \command{\link{print.termstrc_cs}}, \command{\link{rmse}}, \command{\link{aabse}}, \command{\link{summary.lm}}}
