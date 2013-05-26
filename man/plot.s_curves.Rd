\name{plot.s_curves}
\Rdversion{1.1}
\alias{plot.s_curves}

\title{S3 Plot Method}
\description{
 S3 plot method for an object of the class \command{"s_curves"}. 
}
\usage{
\method{plot}{s_curves}(x, xlim = c(range(mapply(function(i) range(x[[i]][, 1]), seq(x)))),
	ylim = c(range(mapply(function(i) range(x[[i]][, 2]), 
	seq(x)))) * 10000, expoints = NULL, 
	xlab = "Maturity (years)", ylab = "Spread (basis points)",
	lwd = 2, lty = 1, main = "Spread curves", ...)
}

\arguments{
  \item{x}{object of the class \command{"s_curves"}.}
  \item{ylim}{the y limits of the plot, for details see \command{\link{plot.default}}.}
  \item{xlim}{the x limits of the plot, for details see \command{\link{plot.default}}.}
  \item{lty}{the line type, for details see \command{\link{par}}.}
  \item{lwd}{the line width, for details see \command{\link{par}}.}
  \item{expoints}{extrapolation points (default: \command{NULL}).}
  \item{ylab}{a label for the y axis, for details see \command{\link{plot.default}}.}
  \item{xlab}{a label for the x axis, for details see \command{\link{plot.default}}.}
  \item{main}{a main title for the plot, for details see \command{\link{title}}.}
  \item{\dots}{other graphical parameters, see \command{\link{par}}.}
}
\details{
The spread curves (the difference of zero-coupon yield curves) are plotted, if at least two groups of bonds
were specified.
}

\seealso{\command{\link{plot.df_curves}}, \command{\link{plot.fwr_curves}}, \command{\link{plot.spot_curves}}}

