termstrc: Zero-coupon Yield Curve Estimation
=============================================

The package offers a wide range of functions for term structure estimation based on static and dynamic coupon bond and yield data sets. The implementation focuses on the cubic splines approach of McCulloch (1971, 1975) and the Nelson and Siegel (1987) method with extensions by Svensson (1994), Diebold and Li (2006) and De Pooter (2007). We propose a weighted constrained optimization procedure with analytical gradients and a globally optimal start parameter search algorithm. Extensive summary statistics and plots are provided to compare the results of the different estimation methods. Several demos are available using data from European government bonds and yields.

You can install the stable version from [CRAN](http://cran.r-project.org/package=termstrc) with:

```r
install.packages("termstrc")
```

## Development

To install the development version of *termstrc*, it's easiest to use the [devtools](http://cran.r-project.org/package=devtools) package:

```r
library(devtools)
install_github("termstrc", "datarob")
```

You can also download the [zip ball](https://github.com/datarob/termstrc/zipball/master) or [tar ball](https://github.com/datarob/termstrc/tarball/master), decompress and run `R CMD INSTALL` on it. (On Windows first install [Rtools](http://cran.rstudio.com/bin/windows/Rtools/)).

