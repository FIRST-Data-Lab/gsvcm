#' Generating samples in Simulation.
#'
#' \code{Datagenerator} is used to generate samples on horseshoe domain for Scenario 1 (Gaussian), Scenario 2 (Poisson), and Scenario 3 (negative binomial).
#'
#' @import mgcv
#' @importFrom MASS rnegbin
#'
#' @param family  The family object, specifying the distribution and link to use. Choose "gaussian" for Gaussian distribution, "poisson" for poisson distribution, and "nb_bps" for negative binomial distribution.
#' \cr
#' @param grid The distance between grid points -- default is set to 0.005.
#' \cr
#' @return A data matrix with a response ('y'), true coefficient functions ('m1' and 'm2'), covariates ('x1' and 'x2'), and locations ('u' and 'v').
#' \cr
#' @details This R package is the implementation program for manuscript entitled "Generalized Spatially Varying Coefficinet Models" by Myungjin Kim and Li Wang.
#' @export
#'

Datagenerator=function(family, grid = 0.005){
  # Construct the HorseShoe boundary
  fsb=list(fs.boundary())

  # simulate some fitting data, inside boundary
  uu=seq(-1,3.5, grid)
  vv=seq(-1,1, grid)
  n1=length(uu)
  n2=length(vv)
  u=rep(uu,n2)
  v=rep(vv,rep(n1,n2))
  m1=fs.test(u,v,b=1)
  N=length(m1)
  m2=m1
  m2[!is.na(m1)]=mapply(beta1_geo, round(u[!is.na(m1)],2), round(v[!is.na(m1)],2))


  set.seed(1234)
  x1=runif(N,0,1)
  x2=runif(N,0,1)


  # Scenario 1. Gaussian
  if(family$family=="gaussian"){
    sigma=0.5
    eps=rnorm(N,mean=0,sd=sigma)
    y=m1*x1+m2*x2+eps

  }

  # Scenarios 2 and 3. Poisson and negative binomial
  if(family$family=="nb_bps"){
    m1=m1
    m2=m2
    log.mu=m1*x1+m2*x2
    mu=exp(log.mu)
    y=mu
    for(i in 1:N){
      y[i]=rnegbin(1,mu=mu[i],theta=6)
    }
  }
  if(family$family=="poisson"){
    m1=m1
    m2=m2
    log.mu=m1*x1+m2*x2
    mu=exp(log.mu)
    y=mu
    for(i in 1:N){
      y[i]=rpois(1,mu[i])
    }
  }
  pop=cbind(y,m1,m2,x1,x2,u,v)
  return(pop)
}


# modified beta 1 coefficient function:
beta1_geo=function(x,y,r=0.5,r0=0.1){
  q=pi/2*r
  if(x>=0 && y>=r0 && y<=(2*r-r0)){
    z = 7/45*(-5*q-2*x+(y-r)^2)
  } else if(x>=0 && -r0>=y && -(2*r-r0)<=y){
    z = 7/45*(5*q+2*x+(y-r)^2)
  } else if(x<0){
    z = 7/45*(5*atan(y/x)*r+(sqrt(x^2+y^2)-r)^2)
  } else {
    z = NA
  }
  return(z+1)
}
