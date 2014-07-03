##' Density, distribution function, quantile function and random
##' generation for the Zero Inflated Laplace distribution with mean equal
##' to 'mu', shape 'b' and zero inflation 'p'. 
##'
##' The mean 'mu' sets the central tendency of the Laplace, with "b"
##' defining the spread around 'mu'. The zero inflation parameter
##' "p", defines the point mass density above 0. The density under
##' zero relates to p as (1/(p+1)) while the probabilty density not
##' equal to zero is given by 1-(1/(p+1)). The function "rzil" is based
##' on the inverse cummulative function, and translates uniform
##' random numbers to their zil equivalents. 
##' 
##' @param mu the mean of the non-zero density
##' @param b the shape parameter defining spread around mu
##' (must be >0).
##' @param p the zero inflation parameter (must be >0) 
##' @author Marco D. Visser
##' @rdname dzil
##' @examples
##' curve(dnorm(x,mu=5),-2,10)
##' x<-rzil(100,mu=5)
##' hist(x)
##' @title The Zero Inflated Laplace 
##' @seealso \code{\link{dnorm}}, \code{\link{dunif}},
##' or \code{\link{runif}}.
##' @concept probabilty density
##' @export
dzil <- function(x,mu=0,b=1,p=1){

  if(sum(p<0)+sum(b<0)){stop("b or p < 0")}

  {(p*0^abs(x))+{exp(-abs(x-mu)/b)/(2*b)}}/(1+p)
 
}

##' Distribution function \code{pzil}
##' @rdname dzil
##' @export
pzil <- function(x,mu=0,b=1,p=1){

  if(sum(p<0)+sum(b<0)){stop("b or p < 0")}
  
  plp <- ifelse(x<mu,.5*exp((x-mu)/b),
                1-(0.5*exp(-((x-mu)/b)))) *(1-(p+1)^(-1))

  ifelse(x>=0,plp+(1/(p+1)),plp)
}

##' Quantile function \code{qzil} 
##' @rdname dzil
##' @export
qzil <- function(q,mu=0,b=1,p=1){
  
  if(sum(p<0)+sum(b<0)){stop("b or p < 0")}
  if(sum(q<0|q>1)){stop("q not a valid quantile")}

  qden <- c(pzil(0-.Machine$double.eps,mu=mu,b=b,p=p),
            pzil(0,mu=mu,b=b,p=p))

  x <- q
  x[q>qden[1]&q<qden[2]] <- 0
  x[q<=qden[1]] <- qlp(q[q<=qden[1]]/(1/(p+1)),mu=mu,b=b)
  x[q>=qden[2]] <- qlp((q[q>=qden[2]]-diff(qden))/(1/(p+1)),mu=mu,b=b)

  return(x)
}

##' Random number generator \code{rzil}
##'
##' @rdname dzil
##' @export
rzil <- function(n,mu=0,b=1,p=1){

  qzil(runif(n,0+.Machine$double.eps,
             1-.Machine$double.eps),mu=mu,b=b,p=p)
  
}

## Helper function in qzil
qlp <-function(q,mu,b) {
  mu-b*sign(q-0.5)*log(1-2*abs(q-0.5))
}



##' Density, distribution function, quantile function and random
##' generation for the Inflated Laplace distribution with mean equal
##' to 'mu', shape 'b' and inflation 'p' at point 'mu'.
##'
##' The mean 'mu' sets the central tendency of the inflated Laplace, with "b"
##' defining the spread around 'mu'. The inflation parameter
##' "p", defines the point mass density inflation above mu.
##' The function "ril" is based
##' on the inverse cummulative function, and translates uniform
##' random numbers to their zil equivalents. 
##' 
##' @param mu the mean of the non-zero density
##' @param b the shape parameter defining spread around mu
##' (must be >0).
##' @param p the inflation parameter (must be >0) 
##' @author Marco D. Visser
##' @rdname dil
##' @examples
##' curve(dnorm(x,mu=5),-2,10)
##' x<-ril(100,mu=5)
##' hist(x)
##' @title The Zero Inflated Laplace 
##' @seealso \code{\link{dnorm}}, \code{\link{dunif}},
##' or \code{\link{runif}}.
##' @concept probabilty density
##' @export
dil <- function(x,mu=0,b=1,p=1){

  if(sum(p<0)+sum(b<0)){stop("b or p < 0")}

  {(0.5*p*0^abs(x-mu))+{exp(-abs(x-mu)/b)/(2*b)}}/(1+p)
 
}

##' Distribution function \code{pil}
##' @rdname dil
##' @export
pil <- function(x,mu=0,b=1,p=1){

  if(sum(p<0)+sum(b<0)){stop("b or p < 0")}
  
  plp <- ifelse(x<mu,.5*exp((x-mu)/b),
                1-(0.5*exp(-((x-mu)/b)))) *(1-(p+1)^(-1))

  ifelse(x>=mu,plp+(1/(p+1)),plp)
}

##' Quantile function \code{qzil} 
##' @rdname dzil
##' @export
qil <- function(q,mu=0,b=1,p=1){
  
  if(sum(p<0)+sum(b<0)){stop("b or p < 0")}
  if(sum(q<0|q>1)){stop("q not a valid quantile")}

  qden <- c(pil(0-.Machine$double.eps,mu=mu,b=b,p=p),
            pil(0,mu=mu,b=b,p=p))

  x <- q
  x[q>qden[1]&q<qden[2]] <- 0
  x[q<=qden[1]] <- qlp(q[q<=qden[1]]/(1/(p+1)),mu=mu,b=b)
  x[q>=qden[2]] <- qlp((q[q>=qden[2]]-diff(qden))/(1/(p+1)),mu=mu,b=b)

  return(x)
}

##' Random number generator \code{ril}
##'
##' @rdname dil
##' @export
ril <- function(n,mu=0,b=1,p=1){

  qzil(runif(n,0+.Machine$double.eps,
             1-.Machine$double.eps),mu=mu,b=b,p=p)
  
}

## Helper function in qil
qlp <-function(q,mu,b) {
  mu-b*sign(q-0.5)*log(1-2*abs(q-0.5))
}


