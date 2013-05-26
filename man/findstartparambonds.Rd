\name{findstartparambonds}
\Rdversion{1.1}
\alias{findstartparambonds}

\title{Find Globally Optimal Startparameters
}
\description{Start parameter search routine for term structure estimation based on a coupon bond data set. 
The algorithm searches for the parameters over a grid spanned over tau1 (tau2). 
}
\usage{
findstartparambonds(p, m, cf, weights, method, tauconstr,
 control = list(), outer.iterations = 30, outer.eps = 1e-04)
}

\arguments{
  \item{p}{price vector
}
  \item{m}{maturites matrix
}
  \item{cf}{cashflows matrix
}
  \item{weights}{duration based weights
}
  \item{method}{form of the spot rate function
}
  \item{tauconstr}{
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
