\name{postpro_bond}
\Rdversion{1.1}
\alias{postpro_bond}

\title{Post Processing of Term Structure Estimation Results}
\description{The function calculates based on the term structure estimation results the errors for prices and yields and differnt curves (spot, forward, discount curve).}
\usage{
postpro_bond(opt_result, m, cf, sgroup, n_group, y, p,
 ac, m_p, method, lambda)
}

\arguments{
  \item{opt_result}{parameter vector}
  \item{m}{maturities matrices}
  \item{cf}{cahsflows matrices}
  \item{sgroup}{sequence of the group length}
  \item{n_group}{lenght of the group}
  \item{y}{yield-to-maturity matrices}
  \item{p}{dirty price vectors}
  \item{ac}{accrued interest vectors}
  \item{m_p}{maturity matrices including the maturities for the current dirty prices}
  \item{method}{form of the spot rate function}
  \item{lambda}{additional paramter for the Diebold/Li spot rate function}
}

\note{Used as internal helper function}

