`zil` The Zero Inflated Laplace Distribution
===

The `zil` package extends R with the Zero Inflated Laplace Distribution and provides functions for calculation of the density, distribution and quantile functions including random number generation. The Zero Inflated Laplace is a flexible distribution used for e.g. meassurment error modeling where where a fraction (1/(1+p)) of errors can be zero, while also allowing for potential bias (parameter 'mu') and large spread in errors (parameter 'b').


## Installation

Currently there isn't a release on [CRAN](http://cran.r-project.org/),
though there may one day be one. You can still  
download the [zip](https://github.com/MarcoDVisser/zil/zipball/master) 
or [tar ball](https://github.com/MarcoDVisser/zil/tarball/master).
Then decompress and run `R CMD INSTALL` on it, 
or use the **devtools** package to install the development version.

```r
## Make sure your current packages are up to date
update.packages()
## devtools is required
library(devtools)
install_github("zil", "MarcoDVisser")
```
## Examples

```r
# using the density function and random numbers
x<-rzil(1000,mu=5)
hist(x,freq=FALSE)
lines(-5:15,dzil(-5:15,mu=5),lwd=3,col='green')
```
![](http://i.imgur.com/QxeArJm.png)

Here parameters include, the  mean 'mu' sets the central tendency of the Laplace, with "b"  defining the spread around 'mu'. The zero inflation parameter  "p", defines the point mass density at 0. The density at zero relates to p as (1/(p+1)) while the probabilty density not  equal to zero is given by 1-(1/(p+1)). The function "rzil" is based  on the inverse cummulative function, and translates uniform random numbers to their zil equivalents. 

