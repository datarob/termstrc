\name{summary.dyntermstrc_nss}
\Rdversion{1.1}
\alias{summary.dyntermstrc_nss}

\title{S3 Summary Method}
\description{Summary method for objects of the class
  \command{"dyntermstrc_nss"}. The average RMSE and AABSE for the prices
  and yields is calculated. Additionally convergence information and the
  message from the used solver is printed. 
}
\usage{
\method{summary}{dyntermstrc_nss}(object, ...)
}

\arguments{
  \item{object}{object of the class \command{"dyntermstrc_nss"}.}
  \item{\dots}{further arguments}
}

\value{The method returns an object of the class \command{"summary.dyntermstrc_nss"}.}
