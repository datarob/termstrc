\name{plot.spsearch}
\Rdversion{1.1}
\alias{plot.spsearch}

\title{S3 Plot Method
}
\description{S3 plot method for objects of the class
  \command{"spsearch"}. The methods plot details on the objective
  function of the start parameter search. 
}
\usage{
\method{plot}{spsearch}(x, main = "Start parameter search", rgl=TRUE, ...)
}

\arguments{
  \item{x}{object of the class \command{"spsearch"}.}
  \item{main}{title}.
  \item{rgl}{if \command{TRUE} (default) the rgl device will be used for the plot.  }
  \item{\dots}{further arguments.}
}
