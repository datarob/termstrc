\name{duration}
\Rdversion{1.1}
\alias{duration}

\title{Duration, modified Duration and Duration based Weights}
\description{
  The function calculates the Macauly duration, modified duration and duration based weights.}
\usage{
duration(cf_p, m_p, y)
}

\arguments{
  \item{cf_p}{cashflows matrix including the current dirty prices of the bonds.}
  \item{m_p}{maturity matrix, the first row is filled with zeros.}
  \item{y}{yields of the bonds.}
}
\details{The duration vector is calculated using the following formula:

\deqn{\bm{d}= \frac{\bm{\iota}'(\bm{C}\cdot\bm{M}\cdot\bm{D})}{\bm{\iota}'(\bm{C}\cdot\bm{D})},}
 whereas \eqn{\bm{C}} is the cashflow matrix and \eqn{\bm{M}} is the maturity matrix. \eqn{\bm{\iota}} is a column vector filled with ones.
\eqn{(\cdot)} denotes a elementwise matrix mulitplication and \eqn{" ' "}  the transpose of a vector (matrix). 

The weight \eqn{\omega_j} for one bond \eqn{j} is defined as
\deqn{  \omega_j=\frac{\frac{1}{d_j}}{\sum_{i=1}^m\frac{1}{d_i}},}
where \eqn{d_j} is the duration of the j-th bond. 
}

\value{
  The function returns a matrix with three columns, i.e,. duration,
  modified duration and duration based weights. 
}


\examples{ 
data(govbonds)
cf_p <- create_cashflows_matrix(govbonds[[1]],include_price=TRUE)
m_p <- create_maturities_matrix(govbonds[[1]],include_price=TRUE)
y <- bond_yields(cf_p,m_p)
duration(cf_p,m_p,y[,2])

}