\name{grad_asv}
\alias{grad_asv}

\title{Gradient of the adjusted Svensson Loss Function for Yields
}
\description{Calculates the gradient of the adjusted Svensson Loss Function for Yields
}
\usage{
grad_asv(beta, m, y)
}

\arguments{
  \item{beta}{Parameter of the adjusted Svensson spot rate function (for
  details see \code{\link{spr_asv}}). 
}
  \item{m}{maturity vector
}
  \item{y}{yield vector
}
}


\seealso{ \code{\link{objfct_asv}},\code{\link{spr_asv}}
}