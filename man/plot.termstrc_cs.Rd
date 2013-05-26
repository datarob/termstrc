\name{plot.termstrc_cs}
\Rdversion{1.1}
\alias{plot.termstrc_cs}

\title{S3 Plot Method for Cubic Splines }
\description{
 S3 plot method for an object of the class \command{"termstrc_cs"}.
}

\usage{
\method{plot}{termstrc_cs}(x, matrange = c(min(mapply(function(i) min(x$y[[i]][, 1]), 
	seq(x$n_group))), max(mapply(function(i) max(x$y[[i]][, 1]), 
	seq(x$n_group)))), multiple = FALSE, 
	ctype = "spot", lwd=2, lty=1, type = "l", 
	errors = "none", inset = c(0.1, 0.3), ask=TRUE, ...)
}

\arguments{
 \item{x}{object of the class \command{"termstrc_cs"}.}
  \item{matrange}{ maturity range for the plot, e.g. c(2,10).}
  \item{multiple}{if \command{TRUE} all curves are plotted together (default: \command{FALSE}).}
  \item{ctype}{ parameter setting for the desired curve type, \command{"spot"} (\command{"forward"}, \command{"discount"}, \command{"spread"}) for the spot rate (forward rate, discount factor, spread) curves. Use \command{"none"} if no curve plot is desired.  }
 \item{errors}{Specify the type of the error plot. If \command{"price"} (\command{"yield"}) the price (yield) errors will be plot. Use \command{"none"} if no error plot is desired. }
\item{lwd}{the line width, for details see \command{\link{par}}.}
\item{lty}{the line type, for details see \command{\link{par}}.}
  \item{type}{1-character string giving the type of plot desired, for details see \command{\link{plot.default}}.}
  \item{inset}{ inset distance(s) from the margins as a fraction of the plot region, for details see \command{\link{legend}}.}
  \item{ask}{if \command{TRUE} (and the \R session is interactive) the user is asked for input, before a new figure is drawn, see \command{\link{par}} for details.}
  \item{\dots}{other graphical parameters, see \command{\link{par}}. }
}
\details{
Depending on the choice of the curve type (\command{"spot"}, \command{"forward"}, \command{"discount"}, \command{"spread"}) the corresponding
curves will be plotted. Either separately or together (\command{mulitple = TRUE}). If the curves are plotted separately also the knot points
used for the estimation of the cubic splines and the yield-to-maturities will be plotted. In addition, with a zero-coupon yield curve plot the 95 \% confidence interval of the curve will be plotted. To ease the analysis of the goodness of the estimation, serveral error plots for the yield and price error are offered. }  


\seealso{ \command{\link{plot.df_curves}},  \command{\link{plot.error}},  \command{\link{plot.fwr_curves}},  \command{\link{plot.ir_curve}},  \command{\link{plot.s_curves}},  \command{\link{plot.spot_curves}},  \command{\link{plot.termstrc_cs}} }
