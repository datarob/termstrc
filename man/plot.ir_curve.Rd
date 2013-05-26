\name{plot.ir_curve}
\Rdversion{1.1}
\alias{plot.ir_curve}


\title{ S3 Plot Method }
\description{
S3 plot method for an object of the class \command{"ir_curve".}
}
\usage{
\method{plot}{ir_curve}(x, ylim = c(), xlim = c(), lwd = 2, type = "l",
	xlab = "Maturity (years)", ylab = "Zero-coupon yields (in percent)",
	col = "steelblue", lty = 1, ...)
}

\arguments{
  \item{x}{object of the class \command{"ir_curve"}.}
  \item{ylim}{the y limits of the plot, for details see \command{\link{plot.default}}.}
  \item{xlim}{the x limits of the plot, for details see \command{\link{plot.default}}.}
  \item{lwd}{the line width, for details see \command{\link{par}}.}
  \item{type}{1-character string giving the type of plot desired, for details see \command{\link{plot.default}}.}
  \item{xlab}{a label for the x axis, for details see \command{\link{plot.default}}.}
  \item{ylab}{a label for the y axis, for details see \command{\link{plot.default}}.}
  \item{col}{ the colors for lines and points.}
  \item{lty}{the line type, for details see \command{\link{par}}.}
  \item{\dots}{ other graphical parameters, see \command{\link{par}}.}
}

