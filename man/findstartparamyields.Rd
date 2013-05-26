\name{findstartparamyields}
\alias{findstartparamyields}

\title{Find Globally Optimal Start Parameters}
\description{Start parameter search routine for term structure estimation based on a yield data set. 
The algorithm searches for the parameters over a grid spanned over tau1 (tau2). 
}


\usage{
findstartparamyields(y,m, method, tauconstr, control = list(), outer.iterations = 30, outer.eps = 1e-04)
}
\arguments{
  \item{y}{yields
}
  \item{m}{maturities
}
  \item{method}{type of spot rate function
}
  \item{tauconstr}{tau parameter constraints
  }
    \item{control}{solver control parameters, for details see \code{\link{optim}} 
}
  \item{outer.iterations}{see \code{\link{constrOptim}}
}
  \item{outer.eps}{see \code{\link{constrOptim}}
}
}
\details{Used as internal helper function
}
\value{Returns an object of the class \command{"spsearch"}, which
  includes the startparameters and details concerning the optimization. 
}