\name{govbonds}
\alias{govbonds}
\docType{data}
\title{European Government Bonds}
\description{European government bonds.}

\usage{data(govbonds)}

\details{
The data set \command{bonds} consists of German, Austrian and French government bonds.
}

\note{
If you use your own data set, make sure that the structure is identical to the provided data sets.
Use the function \command{str()} to explore the data set.}

\examples{
data(govbonds)
str(govbonds) 



# The following code may be used to generate an empty data set,
# which can then be filled with bond data:

ISIN <- vector()
MATURITYDATE <- vector()
STARTDATE <- vector()
COUPONRATE <- vector()
PRICE <- vector()
ACCRUED <- vector()

CFISIN <- vector()
CF <- vector()
DATE <- vector()

CASHFLOWS <- list(CFISIN,CF,DATE)
names(CASHFLOWS) <- c("ISIN","CF","DATE")

TODAY <- vector()

mycountry1 <- list(ISIN,MATURITYDATE,STARTDATE,
                   COUPONRATE,PRICE,ACCRUED,CASHFLOWS,TODAY)
mycountry2 <- list(ISIN,MATURITYDATE,STARTDATE,
                   COUPONRATE,PRICE,ACCRUED,CASHFLOWS,TODAY)

names(mycountry1) <- c("ISIN","MATURITYDATE","STARTDATE","COUPONRATE",
                       "PRICE","ACCRUED","CASHFLOWS","TODAY")
names(mycountry2) <- c("ISIN","MATURITYDATE","STARTDATE","COUPONRATE",
                       "PRICE","ACCRUED","CASHFLOWS","TODAY")

mybonds <- list(mycountry1,mycountry2)

names(mybonds) <- c("mycountry1","mycountry2")

}

\keyword{datasets}
