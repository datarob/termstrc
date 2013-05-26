\name{plot.error}
\Rdversion{1.1}
\alias{plot.error}

\title{S3 Plot Method }
\description{
S3 plot method for an object of the class \command{error}.}
\usage{
\method{plot}{error}(x, type = "b", main = "", mar = c(7, 6, 6, 2) + 0.1, 
	oma = c(4, 2, 2, 2) + 0.1, ylab = "Error", ...)
}

\arguments{
  \item{x}{object of the class \command{error}. }
  \item{type}{1-character string giving the type of plot desired, for details see \command{\link{plot.default}}.}
  \item{main}{a main title for the plot, for details see \command{\link{title}}.}
  \item{mar}{ A numerical vector of the form 'c(bottom, left, top, right)'
          which gives the number of lines of margin to be specified on
          the four sides of the plot, for details see \command{\link{par}}.}
  \item{oma}{A vector of the form 'c(bottom, left, top, right)' giving
          the size of the outer margins in lines of text. }
  \item{ylab}{a label for the y axis, for details see \command{\link{plot.default}}.}
  \item{\dots}{ other graphical parameters, see \command{\link{par}}.}
}
\details{
Absolute yield and price errors as a resuslt of the term structure estimation can be plotted.
The scaling of the x axis depends on the maturity of the bonds, each bond is labeled with its
ISIN number. The error plots seems especially useful in identifying misspriced bonds. For removing them,
the function \command{\link{rm_bond}} may be applied.  
}

\seealso{\command{\link{rm_bond}}}
