\name{plot.dyntermstrc_param}
\Rdversion{1.1}
\alias{plot.dyntermstrc_param}

\title{S3 Plot Method}
\description{Plot method for objects of the class
  \command{"dyntermstrc_param"}. The method is able to plot the
   time series of the parameter levels and first differences and the
   empirical (partial) autocorrelation function. 
}
\usage{
\method{plot}{dyntermstrc_param}(x, type = "param", ...)
}

\arguments{
  \item{x}{object of the class \command{"dyntermstrc_param"}}
  \item{type}{\command{"param"} (default) for the parameters,
  \command{"diffparam"} for the parameter differences and
  \command{"acf"} for the plot of the (partial) autocorrelation function
  of the parameters.}
  \item{\dots}{further arguments}
}

