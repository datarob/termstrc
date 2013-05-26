\name{plot.fwr_curves}
\Rdversion{1.1}
\alias{plot.fwr_curves}

\title{S3 Plot Method}
\description{
 S3 plot method for an object of the class \command{"fwr_curves"}.
}
\usage{
\method{plot}{fwr_curves}(x, multiple = FALSE, 
	ylim = c(range(mapply(function(i) range(x[[i]][, 2]),
	seq(x)))) * 100, xlim = c(), type = "l", lty = 1,
	lwd = 2, expoints = NULL, ylab = "Forward rate (percent)", 
	xlab = "Maturity (years)", main = "Forward rate curves", ...)
}

\arguments{
   \item{x}{object of the class \command{"fwr_curves"}. }
  \item{multiple}{if \command{TRUE} all forward rate curves are plotted together (default: \command{FALSE}). }
  \item{ylim}{the y limits of the plot, for details see \command{\link{plot.default}}. }
  \item{xlim}{the x limits of the plot, for details see \command{\link{plot.default}}.  }
  \item{type}{1-character string giving the type of plot desired, for details see \command{\link{plot.default}}.  }
  \item{lty}{the line type, for details see \command{\link{par}}. }
  \item{lwd}{the line width, for details see \command{\link{par}}.}
  \item{expoints}{extrapolation points (default: \command{NULL}).}
  \item{ylab}{a label for the y axis, for details see \command{\link{plot.default}}.}
  \item{xlab}{a label for the x axis, for details see \command{\link{plot.default}}.}
  \item{main}{a main title for the plot, for details see \command{\link{title}}.}
  \item{\dots}{ other graphical parameters, see \command{\link{par}}.}
}
  

\seealso{\command{\link{plot.df_curves}},  \command{\link{plot.s_curves}},  \command{\link{plot.spot_curves}}}

