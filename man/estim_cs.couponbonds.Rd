\name{estim_cs.couponbonds}
\Rdversion{1.1}
\alias{estim_cs.couponbonds}

\title{S3 Estim_cs Method}
\description{S3 estim.cs method for an object of the class \command{"couponbonds"}. The method estimates the discount curve with the cubic splines approach by McCulloch (1975).}


\usage{
\method{estim_cs}{couponbonds}(bonddata, group, matrange = "all", rse = TRUE)
}

\arguments{
 \item{bonddata}{a data set of bonds in list format.}
 \item{group}{vector defining the group of bonds used for the estimation, 

e.g., \command{c("GERMANY","AUSTRIA")}.}
  \item{matrange}{use \command{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \command{c(lower,upper)}. }

  \item{rse}{ \command{TRUE} (default) calculates robust standard erros for the confidence intervalls of the discount curve}
}


\details{\itemize{
 \item{\command{group}}{The first element of the vector will be used as the
  			    	  	reference country for the spread curve calculation.
					 \command{group} can be either a vector of bond groups or a scalar.}
\item{\command{bonddata}}{The package is designed to work with a certain list data structure. For more information use the function \command{str()} to explore the structure of the example data sets.}
}}

\value{
  The function returns an object of the class \command{"termstrc_cs"}. The object
contains the following items (mainly lists):
  \item{group}{group of bonds (e.g. countries) used for the estimation.}
  \item{matrange}{\command{"none"} or a vector with the maturity range.}
  \item{n_group}{length of object \command{group}, i.e. the number of countries.}
  \item{knotpoints}{selected knot points for the cubic splines estimation.}
 \item{spot}{zero-coupon yield curves as object of the class \command{"spot_curves"}.}
  \item{spread}{spread curves as object of the class \command{"s_curves"}.}
  \item{forward}{forward curves as object of the class \command{"fwr_curves"}.} 
 \item{discount}{discount curves as object of the class \command{"df_curves"}.} 
 \item{cf}{cashflow matrices.}
  \item{m}{maturity matrices.}
  \item{p}{dirty prices.}
  \item{phat}{estimated bond prices.}
  \item{perrors}{pricing errors and maturities as object of the class \command{"error"}.}
  \item{y}{bond yields.}
  \item{yhat}{one list for each group with the theoretical bond yields calculated with the estimated bond prices \command{phat}.}
  \item{yerrors}{yield errors and maturities as object of the class \command{"error"}.}
  \item{alpha}{OLS coefficients of cubic splines estimation.}
  \item{regout}{OLS estimation results as object of the class \command{"lm"}.}
  \item{rse}{robust standard errors for confidence interval calculation}
}
\references{
J.Huston McCulloch (1971):
Measuring the Term Structure of Interest Rates.
\emph{The Journal of Business, \bold{44} 19--31.} 

J. Huston McCulloch (1975):
The Tax-Adjusted Yield Curve.
\emph{The Journal of Finance, \bold{30} 811--830.} 
}



\note{ For objects of the class \command{"spot_curves"}, 
\command{"s_curves"}, \command{"df_curves"}, \command{"fwr_curves"}, \command{"error"} appropriate plot methods are offered. For objects of the list item \command{regout} standard \command{lm} methods apply. For objects of the class \command{"termstrc_cs"} print, summary and plot methods are available. Another term structure estimation method is provided by the method 

\command{\link{estim_nss.couponbonds}}.} 

\seealso{\command{\link{print.termstrc_cs}}, \command{\link{summary.termstrc_cs}},  \command{\link{plot.termstrc_cs}},

  \command{\link{estim_nss.couponbonds}},  \command{\link{plot.spot_curves}},  \command{\link{plot.s_curves}},  \command{\link{plot.df_curves}}, 

\command{\link{plot.fwr_curves}},  \command{\link{plot.error}},  \command{\link{summary.lm}},  \command{\link{print.lm}},  \command{\link{plot.lm}}.}

\examples{
# load data set
data(govbonds)

# define countries, for which the estimation 
# of the zero-coupon yield curves will be carried out
group <- c("GERMANY", "AUSTRIA")

# set maturtiy range
matrange <- c(0, 19)  

# perform estimation
x <- estim_cs(govbonds, group, matrange)

# print the obtained parameters of the estimation
print(x)

# goodness of fit measures
summary(x)

# plot the zero-coupon yield curve for each country
plot(x,errors="none")

# plot all zero-coupon yield curves together
plot(x,multiple=TRUE,errors="none")

# spread curve splot
plot(x,ctype="spread",errors="none")

# price error plot for all countries
plot(x,ctype="none")
}
