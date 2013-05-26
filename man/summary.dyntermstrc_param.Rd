\name{summary.dyntermstrc_param}
\Rdversion{1.1}
\alias{summary.dyntermstrc_param}

\title{S3 Summary Method}

\description{S3 summary method for objects of the class \command{"dyntermstrc_param"}.}
\usage{
\method{summary}{dyntermstrc_param}(object, type = "none", lags = 1, selectlags = "Fixed", ...)
}

\arguments{
  \item{object}{object of the class \command{"dyntermstrc_param"}.}
  \item{type}{use \command{"trend"} and a trend is considered for the
    unit root test (default: \command{"none"}).}
  \item{lags}{ number of lags for unit root test function
    \command{ur.df} from package \command{urca} (default:\command{1})}
  \item{selectlags}{leg selection flag for function \command{ur.df} from
    package \command{urca} (default: \command{"Fixed"}).}
  \item{\dots}{further arguments}}
\details{
The function calculates from an object of the class
    \command{"dyntermstrc_param"} the augmented Dickey Fuller test for
    the levels and first differences for each of the estimated
    term structure parameters. Additionally, the correlation of the
    parameter levels and differences are calculated.  }
\value{an object with the class \command{"summary.dyntermstrc_param"}}



