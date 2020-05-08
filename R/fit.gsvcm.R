#' Fitting generalized regression models
#'
#' \code{fit.gsvcm} fits the generalized spatially varying coeffcient models.
#'
#' @importFrom MGLM kr
#' @import BPST
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
#' @param family The family object, specifying the distribution and link to use.
#' \cr
#' @param off offset -- default is 0.
#' \cr
#' @param r.theta The endpoints of an interval to search for an additional parameter \code{theta} for negative binomial scenario -- default is c(2,8).
#' \cr
#' @param eps.sigma Error tolerance for the Pearson estimate of the scale parameter, which is as close as possible to 1, when estimating an additional parameter \code{theta} for negative binomial scenario -- default is 0.01.
#' \cr
#' @param method GSVCM or GSVCMQR. GSVCM is based on Algorithm 1 in Subsection 3.1 and GSVCMQR is based on Algorithm 2 in Subsection 3.2 -- default is GSVCM.
#' \cr
#' @param Cp TRUE or FALSE. There are two modified measures based on the QRGSVCM method for smoothness parameters in the manuscript. TRUE is for Cp measure and FALSE is for GCV measure.
#' \cr
#' @return The function returns an object with S3 class "gsvcm" with the following items:
#' \item{beta}{The estimated coefficient functions.}
#' \item{theta_hat}{The estimated spline coefficient functions.}
#' \item{lambdac}{Selected tuning (penalty) parameter for bivariate penalized spline based on GCV.}
#' \item{gcv}{Generalized cross-validation (GCV).}
#' \item{df}{Effective degree of freedom.}
#' \item{theta}{The estimated additional parameter for a negative binomial random component.}
#' \item{V}{The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.}
#' \item{Tr}{The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.}
#' \item{d}{The degree of piecewise polynomials}
#' \item{r}{The smoothness parameter}
#' \item{B}{The spline basis function of dimension \code{n} by \code{nT}*\code{{(d+1)(d+2)/2}}, where \code{n} is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. The length of points means the length of ordering indices of observation points. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.}
#' \item{Q2}{The Q2 matrix after QR decomposition of the smoothness matrix \code{H}.}
#' \item{K}{The thin-plate energy function.}
#' \item{ind.inside}{A vector contains the indexes of all the points which are inside the triangulation.}
#' \item{tria.all}{The area of each triangle within the given triangulation.}
#' \item{X}{The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.}
#' \item{y}{The response of dimension \code{n} by one, where \code{n} is the number of observations.}
#' \item{S}{The cooridinates of dimension \code{n} by two. Each row is the coordinates of an observation.}
#' \item{family}{The family object}
#' \cr
#' @details This R package is the implementation program for manuscript entitled "Generalized Spatially Varying Coefficinet Models" by Myungjin Kim and Li Wang.
#' \cr
#' @examples
# # Population:
#' family=poisson()
#' ngrid = 0.05
#'
#' # Data generation:
#' all_pop = as.matrix(Datagenerator(family, ngrid))
#' pop.r=all_pop[!is.na(all_pop[,'m1']),]
#' N=nrow(pop.r)
#'
#' # Triangulations and setup:
#' Tr = Tr0_horse; V = V0_horse; n = 1000; d = 2; r = 1
#'
#' # set up for smoothing parameters in the penalty term:
#' lambda_start=0.0001; lambda_end=10; nlambda=10
#' lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))
#'
#' # Generate Sample:
#' ind.s=sample(N,n,replace=FALSE)
#' data=as.matrix(pop.r[ind.s,])
#' y=data[,1]; beta0=data[,c(2:3)]; X=data[,c(4:5)]; S=data[,c(6:7)]
#'
#' # Fit the model:
#' mfit0 = fit.gsvcm(y, X, S, V, Tr, d, r, lambda, family)
#'
#' # Prediction:
#' y_hat = predict(mfit0, X, S)
#'
#' # 10 Crossvalidation:
#' MSPE = cv.gsvcm(y, X, S, V, Tr, d = d, r = r, lambda, family, off = 0,  r.theta = c(4, 8), nfold = 10)
#'
#' # GQLR Test:
#' test.gsvcm(y, X, S, V, Tr, d, r, lambda, test_iter = 1, family = family, nB = 20)
#' @export
#'

fit.gsvcm =
function(y, X, S, V, Tr, d = 2, r = 1, lambda = 10^seq(-6, 6, by = 0.5),
         family, off = 0, r.theta = c(2, 8), eps.sigma = 0.01, method = "GSVCM", Cp = TRUE){
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
  this.call = match.call()

  Basis.full = basis(V, Tr, d, r, S)
  K = Basis.full$K
  Q2 = Basis.full$Q2
  B = Basis.full$B
  ind.inside = Basis.full$Ind.inside
  tria.all = Basis.full$tria.all

  y = y[ind.inside]
  X = as.matrix(X[ind.inside, ])
  S = as.matrix(S[ind.inside, ])

  BQ2 = B %*% Q2
  P = t(Q2) %*% K %*% Q2

  if(method == "GSVCM"){
    mfit = gsvcm.est(y, X, BQ2, P, lambda, family, off, r.theta, eps.sigma)
  } else if (method =="GSVCMQR"){
    mfit = gsvcm.est.qr(y, X, BQ2, P, lambda, family, off, r.theta, eps.sigma, Cp)
  }

  mfit$V = V; mfit$Tr = Tr; mfit$d = d; mfit$r = r; mfit$B = B;
  mfit$Q2 = Q2; mfit$K = K; mfit$ind.inside = ind.inside;
  mfit$tria.all = tria.all; mfit$X = X; mfit$y = y; mfit$S = S;
  mfit$family = family;  mfit$call = this.call;
  class(mfit) = "gsvcm"
  return(mfit)
}
