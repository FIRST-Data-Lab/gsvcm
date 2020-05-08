#' Generalized Quasi-Likelihood Ratio (GQLR) test for generalized spatially varying coefficient regression
#'
#' \code{test.gsvcm} implements Generalized Quasi-Likelihood Ratio (GQLR) test using wild bootstrap.
#'
#' @importFrom BPST basis
#'
#' @param y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param S The cooridinates of dimension \code{n} by two. Each row is the coordinates of an observation.
#' \cr
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 2.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param lambda The vector of the candidates of penalty parameter -- default is grid points of 10 to the power of a sequence from -6 to 6 by 0.5.
#' \cr
#' @param test_iter The function to be tested -- default is the first function, which corresponds to the first covariate.
#' \cr
#' @param family The family object, specifying the distribution and link to use.
#' \cr
#' @param off offset -- default is 0.
#' \cr
#' @param r.theta smoothness penalty parameter candidates -- default is c(2,8).
#' \cr
#' @param nB The number of boostrap replication -- default is 100.
#' \cr
#' @param initial The seed used for bootstrap sample in the GQLR test -- default is 123.
#' \cr
#' @param eps.sigma Error tolerance for the Pearson estimate of the scale parameter, which is as close as possible to 1, when estimating an additional parameter \code{theta} for negative binomial scenario -- default is 0.01.
#' \cr
#' @return The function returns a list with the following items:
#' \item{b.GQLR}{The GQLR test statistics based on sample.}
#' \item{obs.GQLR}{The GQLR test statistics based on each boostrap iteration.}
#' \item{pvalue}{The p-value for the GQLR test.}
#'
#' @details This R package is the implementation program for manuscript entitled "Generalized Spatially Varying Coefficinet Models" by Myungjin Kim and Li Wang.
#' \cr
#' @examples
#' # See an example of fit.gsvcm.
#'
#' @export
#'
test.gsvcm =
function(y, X, S, V, Tr, d, r, lambda = 10^seq(-6, 6, by = 0.5),
         test_iter = 1, family, off = 0, r.theta = c(2, 8), nB = 100, initial = 123, eps.sigma = 0.01)
{

  if(!is.vector(y)){
    warning("The response variable, y, should be a vector.")
    y = as.vector(y)
  }


  if(!is.matrix(X)){
    warning("The explanatory variable, X, should be a matrix.")
    X = as.matrix(X)
  }
  if(!is.matrix(S)){
    warning("The coordinates, S, should be a matrix.")
    S = as.matrix(S)
  }
  if(ncol(X) == 1){
    stop("The number of covariates should be larger than at least two for a test.")
  }

  # Environmental Variables
  variance=family$variance
  linkinv=family$linkinv
  linkfun=family$linkfun
  mu.eta=family$mu.eta
  initialize=family$initialize

  cat("\n Model fitting.....\n")
  mfit=fit.gsvcm(y, X, S, V, Tr, d, r, lambda, family, off, r.theta, eps.sigma)
  y_bpst = predict.gsvcm(mfit,X,S)
  B=mfit$B
  K=mfit$K
  Q2=mfit$Q2
  Ind=mfit$Ind.inside

  yi=mfit$y
  Xi=mfit$X
  Si=mfit$S
  n=length(yi)
  np=ncol(Xi)

  result=mfit
  BQ2=B%*%Q2
  P=t(Q2)%*%K%*%Q2


  eta1 = as.matrix(kr(Xi,BQ2,byrow=TRUE)%*%as.vector(result$theta_hat))

  y_hat=linkinv(eta1)
  if(family$family != "nb_bps") tt=(yi-y_hat)^2/variance(y_hat) else tt=(yi-y_hat)^2/variance(y_hat,result$theta)
  sigma_2=1/(length(yi)-result$df)*sum(tt)


  ######## Wild Boostrap Test ########
  # Reduced model test using y
  ind.l = test_iter
  can.nl = seq(1:dim(Xi)[2])
  ind.nl = setdiff(can.nl,ind.l)
  X.l=as.matrix(Xi[,ind.l])
  X.nl=as.matrix(Xi[,ind.nl])

  J=dim(BQ2)[2]
  n.l=ncol(X.l)
  n.nl=ncol(X.nl)

  result2=gsvcm.est.reduce(yi,X.l,X.nl,BQ2,P,result$lambdac,family,off,result$theta)

  eta2=X.l*result2$theta_hat[1:n.l]+kr(X.nl,BQ2,byrow=TRUE)%*%result2$theta_hat[(n.l+1):(n.l+n.nl*J)]
  eta2=as.matrix(eta2)

  y_hat2=linkinv(eta2)
  obs.GQLR = gsvcm.GQLR(yi = yi, y_hat1 = y_hat, y_hat2 = y_hat2, family = family, theta = result$theta)

  b.GQLR=c()
  no.neg=c()
  start.time=Sys.time()

  set.seed(initial)
  cat('\014')
  for(biter in 1:nB){
    cat(paste0(round(biter / nB * 100), '% completed'))

    res.full=yi-y_hat

    vi=c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
    mammen.prob=c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5)))

    eps.star.mammen=sample(vi, size=length(yi), replace=TRUE, prob=mammen.prob)*res.full
    y.star.mammen=y_hat2+eps.star.mammen
    y.star=y.star.mammen

    if(family$family == "gaussian"){
        y.star = y.star
      } else {
        y.star[y.star < 0] = 0
    }


    # boostrap full model
    result3 = fit.gsvcm(y.star,Xi,Si,V,Tr,d,r,lambda,family,off,c(result$theta, result$theta))
    eta3 = as.matrix(kr(Xi,BQ2,byrow=TRUE)%*%as.vector(result3$theta_hat))
    y_hat3=linkinv(eta3)

    # boostrap reduced model
    result4=gsvcm.est.reduce(y.star, X.l, X.nl, BQ2, P, result3$lambdac, family,off=0, result$theta)
    eta4 = X.l*result4$theta_hat[1:n.l]+kr(X.nl,BQ2,byrow=TRUE)%*%result4$theta_hat[(n.l+1):(n.l+n.nl*J)]
    eta4 = as.matrix(eta4)
    y_hat4=linkinv(eta4)

    b.obs.GQLR = gsvcm.GQLR(yi = y.star, y_hat1 = y_hat3, y_hat2 = y_hat4, family = family, theta = result$theta)
    b.GQLR=c(b.GQLR, b.obs.GQLR)

    if (biter == nB) cat(': Done \n') else cat('\014')
  }
  end.time=Sys.time()
  end.time-start.time
  est.pvalue=mean(b.GQLR>obs.GQLR)
  list(b.GQLR = b.GQLR, obs.GQLR = obs.GQLR, pvalue = est.pvalue)
}


#' Estimation for reduced model (GSVCMs with partially linear part)
#'
#' This is an internal function of package \code{gsvcm}.
#'
gsvcm.est.reduce=function(y,X.l,X.nl,BQ2,P,lambda,family,off=0,theta=0){
  #Environmental Variables
  variance=family$variance
  linkinv=family$linkinv
  linkfun=family$linkfun
  mu.eta=family$mu.eta
  gprime=family$gprime
  initialize=family$initialize

  J=dim(BQ2)[2]
  nobs=length(y)
  n.l=ncol(X.l)
  n.nl=ncol(X.nl)

  if(n.l==0 | is.null(X.l)){
    n.l=0
    X.l=NULL
  }
  if(n.nl==0 | is.null(X.nl)){
    n.nl=0
    X.nl=NULL
    XB=NULL
  }else{
    XB=kr(X.nl,BQ2,byrow=TRUE)
  }

  Z = as.matrix(cbind(X.l, XB))

  D=matrix(0,(n.l+n.nl*J),(n.l+n.nl*J))
  D[(n.l+1):(n.l+n.nl*J),(n.l+1):(n.l+n.nl*J)]=as.matrix(kronecker(diag(rep(1,n.nl)),P))


  weights=rep(1,nobs)
  nl=length(lambda)

  alpha_all=matrix(rep(0,(n.l+n.nl*J)*nl),ncol=nl)
  gcv_all=rep(0,nl)

  for(il in 1:nl){
    mustart = y + 0.1
    etastart=linkfun(mustart)
    if(family$family != "nb_bps") var=variance(mustart) else var=variance(mustart,theta)
    mevg=mu.eta(etastart)

    Y_iter=(etastart-off)+(y-mustart)/mevg
    W_iter=(mevg^2)/var
    W_iter=as.vector(W_iter)
    temp1=as.matrix(W_iter*Z)
    temp2=as.matrix(W_iter*Y_iter)
    temp3= t(Z) %*% temp1 + lambda[il]*D
    alpha_old=solve(temp3,t(Z) %*% temp2)

    step=0
    delta1=1
    delta2=1
    while(delta1>1e-5 & sum(is.infinite(delta2))==0 & step<=10 ){
      step=step+1
      eta=Z%*%alpha_old+off

      mu=linkinv(eta)
      if(family$family != "nb_bps") var=variance(mu) else var=variance(mu,theta)
      mevg=mu.eta(eta)
      Y_iter=(eta-off)+(y-mu)/mevg
      W_iter=as.vector((mevg^2)/var)
      temp1=as.matrix(W_iter*Z)
      temp2=as.matrix(W_iter*Y_iter)
      temp3=t(Z) %*% temp1+lambda[il]*D
      alpha_new=solve(temp3,t(Z) %*% temp2)

      eta_new=Z%*%alpha_new
      delta1=sqrt(mean((alpha_new-alpha_old)^2))
      delta2=exp(eta_new)
      if(sum(is.infinite(delta2))==0){
        alpha_old=alpha_new
      }
    }
    alpha_all[,il]=as.matrix(alpha_old)

    #calculate gcv
    DD=solve(temp3)
    Slambda=Z%*%tcrossprod(DD,temp1)
    df=sum(diag(Slambda))
    yhat=Z%*%alpha_old
    gcv=nobs*sum(W_iter*(Y_iter-yhat)^2)/(nobs-df)^2;

    gcv_all[il]=gcv
  }

  j=which.min(gcv_all)
  gcv=gcv_all[j]
  lambdac=lambda[j]
  alpha_hat=alpha_all[,j]

  theta_hat=alpha_hat
  list(theta_hat=theta_hat,
       lambdac=lambdac,gcv=gcv,df=df)
}


#' Quasi-likelihood test statistics based on a variety of \code{family}, evaluated by full and reduced estimates.
#'
#' This is an internal function of package \code{gsvcm}.

gsvcm.GQLR =
  function(yi, y_hat1, y_hat2, family, theta)
  {
    if(family$family == "gaussian"){
      obs.GQLR =
        sum(quasi.gaussian.theta(yi, y_hat1) - quasi.gaussian.theta(yi, y_hat2), na.rm = TRUE)
    } else if(family$family == "poisson"){
      obs.GQLR =
        sum(quasi.poi.theta(yi, y_hat1) - quasi.poi.theta(yi, y_hat2), na.rm = TRUE)
    } else if (family$family == "Gamma"){
      obs.GQLR =
        sum(quasi.Gamma(yi, y_hat1) - quasi.Gamma(yi, y_hat2), na.rm = TRUE)
    } else if (family$family == "inverse.gaussian"){
      obs.GQLR = sum(quasi.inv.gauss(yi, y_hat1) - quasi.inv.gauss(yi, y_hat2), na.rm = TRUE)
    } else if (family$family == "binomial"){
      obs.GQLR =
        sum(quasi.binom (yi, y_hat1) - quasi.binom(yi, y_hat2), na.rm = TRUE)
    } else if(family$family == "nb_bps"){
      obs.GQLR =
        sum(quasi.nb.theta(yi, y_hat1, theta) - quasi.nb.theta(yi, y_hat2, theta), na.rm = TRUE)
    } else {
      stop("Unavailable family is used.")
    }
    return(obs.GQLR)
  }


#' Quasi-likelihoods
#'
#' This is an internal function of package \code{gsvcm}.
#'
quasi.gaussian.theta = function(y, mu){
  return(-(y - mu)^2 / 2)
}

quasi.poi.theta = function(y, mu){
  return(y * log(mu) - mu)
}

quasi.Gamma = function(y, mu){
  return(-y / mu - log(mu))
}

quasi.inv.gauss = function(y, mu){
  return(-y / (2 * mu^2) + 1 / mu)
}

quasi.binom = function(y, mu){
  return(y * log(mu / (1 - mu)) + log(1 - mu))
  # y/m where y follows a binomial distribution and m is the number of trials.
}

quasi.nb.theta = function(y, mu, theta){
  return(y * log(mu / (mu + theta)) + theta * log(theta / (mu + theta)))
}

