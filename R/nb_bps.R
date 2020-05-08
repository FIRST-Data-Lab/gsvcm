#' Negative Binomial Family
#' @export

nb_bps=function(link='log',theta){
  if(link=='log') {
    linkfun=function(mu) log(mu)
    linkinv=function(eta) exp(eta)
    mu.eta=function(eta) exp(eta)
    gprime=function(mu) 1/mu
  }
  variance=function(mu,theta) mu+mu^2/theta
  initialize=function(y){
    if (any(y < 0))
      stop("negative values not allowed for the 'Negative binomial' family")
    y+0.1
  }
  structure(list(family = "nb_bps", link=link,
                 linkfun = linkfun, linkinv = linkinv, variance = variance,
                 mu.eta = mu.eta, gprime = gprime,
                 initialize = initialize), class = "family")
}



