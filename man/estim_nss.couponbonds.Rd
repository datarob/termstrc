\name{estim_nss.couponbonds}
\Rdversion{1.1}
\alias{estim_nss.couponbonds}

\title{S3 Estim_nss Method
}
\description{
Zero-coupon yield curve estimation with the parametric Nelson/Siegel
(1987), Svensson (1994) and Diebold/Li (2006) method. The method
requires an object of the class \command{"couponbonds"}. 
}
\usage{
\method{estim_nss}{couponbonds}(dataset, group, matrange = "all", method = "ns",
 startparam = NULL, lambda = 0.0609 * 12, tauconstr = NULL, 
constrOptimOptions = list(control = list(maxit = 2000),
 outer.iterations = 200, outer.eps = 1e-04),...)
}


\arguments{
  \item{dataset}{a static coupon bond data set of the class \command{"couponbonds"}
}
  \item{group}{vector defining the group of bonds used for the estimation,

               e.g., \command{c("GERMANY","AUSTRIA")}. The spot rate
		 curve of the first group element will be used as the reference curve for the
		 spread curve calculation.
}
  \item{matrange}{use \command{"all"} for no restrictions, or restrict
    the maturity range (in years) used for the estimation with \command{c(lower,upper)}.
}
  \item{method}{\command{"ns"} for Nelson/Siegel (default), \command{"dl"} for
    Diebold/Li, \command{"sv"} for Svensson or \command{"asv"} for adjusted Svensson.
}
  \item{startparam}{matrix of start parameters (number of columns is the
    number of parameters). If no start parameters are given, globally
    optimal parameters are searched automatically (default: \command{NULL})
}
  \item{lambda}{parameter on a yearly time scale with fixed value for \command{"dl"} spot rate
    function (default: 0.0609*12)
}
  \item{tauconstr}{
}
  \item{constrOptimOptions}{list with solver control parameters
    (default: control=list(), outer.interations=30,
    outer.eps.=1e-04). For further documentation please refer to \command{\link{optim}}}
  \item{\dots}{further arguments}
}

\value{
  The function \command{nelson_estim} returns an object of the class \command{"nelson"}. The object
contains the following items (mainly lists):
  \item{group}{group of bonds (e.g. countries) used for the estimation.}
  \item{matrange}{\command{"none"} or a vector with the maturity range.}
  \item{method}{estimation method (\command{"Nelson/Siegel"} or
    \command{"Svensson"}).}
  \item{startparam}{calculated starparameters.}
  \item{n_group}{length of object \command{group}, i.e. the number of
    countries.}
  \item{lambda}{lambda parameter of \command{"dl"} spot rate function.}
  \item{spsearch}{detailed data from the start parameter search algorithm}
 \item{spot}{zero-coupon yield curves as object of the class \command{"spot_curves"}.}
  \item{spread}{spread curves as object of the class \command{"s_curves"}.}
  \item{forward}{forward curves as object of the class \command{"fwr_curves"}.} 
 \item{discount}{discount curves as object of the class \command{"df_curves"}.} 
 \item{expoints}{extrapolation points for Nelson/Siegel method.}
 \item{cf}{cashflow matrices.}
 \item{m}{maturity matrices.}
 \item{duration}{duration matrix, including the modified duration and
  duration based weights.}
\item{p}{dirty prices.}
  \item{phat}{estimated bond prices.}
  \item{perrors}{pricing errors and maturities, object of the class
    \command{"error"}.}
  \item{ac}{accrued interest}
  \item{y}{bond yields.}
  \item{yhat}{one list for each group with the theoretical bond yields calculated with the estimated bond prices \command{phat}.}
  \item{yerrors}{yield errors and maturities as object of the class \command{"error"}.}
  \item{opt_result}{optimization results from \command{optim}, e.g. optimal parameters, convergence info.}

  
}


\references{ 
Charles R. Nelson and Andrew F. Siegel (1987): 
Parsimonious Modeling of Yield Curves.
\emph{The Journal of Business}, \bold{60(4)}:473--489.

Lars E.O. Svensson (1994): 
Estimating and Interpreting Forward Interest Rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}.


F.X. Diebold and C. Li: Forecasting the Term Structure of Government Bond Yields. \emph{Journal 
of Econometrics}, \bold{130}:337--364.

 }




\note{An error message concerning the function \command{uniroot()} is in
  general caused by wrongly specified start parameters or by data issues. 
 
For objects of the class \command{"spot_curves"}, \command{"s_curves"}, \command{"df_curves"}, \command{"fwr_curves"}, \command{"error"} appropriate plot methods are offered.  For objects of the class \command{"termstrc_nss"} print, summary and plot methods are available. Another term structure estimation method is provided by the function \command{\link{estim_cs}}.} 

\seealso{\command{\link{print.termstrc_nss}}, \command{\link{summary.termstrc_nss}},  \command{\link{plot.termstrc_nss}},  \command{\link{estim_cs}},  \command{\link{plot.spot_curves}},  \command{\link{plot.s_curves}},  \command{\link{plot.df_curves}},  \command{\link{plot.fwr_curves}},  

\command{\link{plot.error}},  \command{\link{uniroot}}.}


\examples{
## Run: demo(nss_static)}