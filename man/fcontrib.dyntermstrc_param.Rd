\name{fcontrib.dyntermstrc_param}
\alias{fcontrib.dyntermstrc_param}
\title{S3 fcontrib Method
}
\description{S3 fcontrib method for objects of the class \command{"dyntermstrc_param"}. The function plots the factor contribution of the parameters of the
different spot rate functions at a certain point in time. 
}
\usage{
\method{fcontrib}{dyntermstrc_param}(x, method = "ns", lambda = 0.0609 * 12, index = 1, 
m = 1:10, ylim = NULL, ...)
}

\arguments{
  \item{x}{ object of the class \command{dyntermstrc_param}

}
  \item{method}{Spot rate function, one of the following \command{"ns","sv","dl","sv"}
}
  \item{lambda}{additional parameter for \command{"dl"} spot rate
  function. 
}
  \item{index}{ specific point in time
}
  \item{m}{maturity spectrum for the plot, e.g., \command{"c(min,max)"}
}
  \item{ylim}{range of the y axis. 
}
  \item{\dots}{further arguments
}
}

