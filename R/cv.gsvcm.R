#' k-fold cross-validation MSPE for generalized spatially varying coefficient regression
#'
#' \code{cv.gsvcm} implements k-fold cross-validation MSPE for generalized spartially varying coefficient regression, and returns the mean squared prediction error (MSPE).
#'
#' @importFrom BPST basis
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
#' @param r.theta  The endpoints of an interval to search for an additional parameter \code{theta} for negative binomial scenario -- default is c(2,8).
#' \cr
#' @param nfold The number of folds -- default is 10. Although \code{nfold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable for \code{nfolds} is 3.
#' \cr
#' @param initial The seed used for cross-validation sample -- default is 123.
#' \cr
#' @param eps.sigma Error tolerance for the Pearson estimate of the scale parameter, which is as close as possible to 1, when estimating an additional parameter \code{theta} for negative binomial scenario -- default is 0.01.
#' \cr
#' @param method GSVCM or GSVCMQR. GSVCM is based on Algorithm 1 in Subsection 3.1 and GSVCMQR is based on Algorithm 2 in Subsection 3.2 -- default is GSVCM.
#' \cr
#' @param Cp TRUE or FALSE. There are two modified measures based on the QRGSVCM method for smoothness parameters in the manuscript. TRUE is for Cp measure and FALSE is for GCV measure.
#' \cr
#' @return The k-fold cross-validation (CV) mean squared prediction error (MSPE).
#' \cr
#' @details This R package is the implementation program for manuscript entitled "Generalized Spatially Varying Coefficinet Models" by Myungjin Kim and Li Wang.
#' @examples
#' # See an example of fit.gsvcm.
#' @export
#'

cv.gsvcm =
function(y, X, S, V, Tr, d = 2, r = 1, lambda = 10^seq(-6, 6, by = 0.5),
                     family, off = 0, r.theta = c(2, 8), nfold = 10, initial = 123, eps.sigma = 0.01, method = "GSVCM", Cp =TRUE)
{
  linkinv = family$linkinv;

  if(nfold < 3){
    warning("The number of folds in CV is too small. Instead, the default 10-fold CV is used.")
    nfold = 10
  }

  if(!is.matrix(X)){
    warning("The explanatory variable, X, should be a matrix.")
    X = as.matrix(X)
  }
  if(!is.matrix(S)){
    warning("The coordinates, S, should be a matrix.")
    S = as.matrix(S)
  }

  Ball = basis(V, Tr, d, r, S)
  K = Ball$K
  Q2 = Ball$Q2
  B = Ball$B
  ind.inside = Ball$Ind.inside
  tria.all = Ball$tria.all

  BQ2 = B %*% Q2
  P = t(Q2) %*% K %*% Q2
  y = y[ind.inside]
  X = X[ind.inside, ]
  S = S[ind.inside, ]
  n = length(y)

  sfold = round(n / nfold)
  set.seed(initial)
  Test = sample(1:n)
  cv.error = c()

  for(ii in 1:nfold){
    if(ii < nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):(ii * sfold)])
    }
    if(ii == nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):n])
    }
    Train.set = setdiff(1:n, Test.set)

    # Consider univariate case.
    if(is.vector(X) == 1){
      X.test = as.matrix(X[Test.set])
      X.train = as.matrix(X[Train.set])
    } else {
      X.test = X[Test.set, ]
      X.train = X[Train.set, ]
    }

    B.test = B[Test.set, ]
    B.train = B[Train.set, ]
    BQ2.test = BQ2[Test.set, ]
    BQ2.train = BQ2[Train.set, ]
    y.test = y[Test.set]
    y.train = y[Train.set]

    if(method == "GSVCM"){
      mfit.ii = gsvcm.est(y.train, X.train, BQ2.train, P, lambda, family, off, r.theta, eps.sigma)
    } else if (method =="GSVCMQR"){
      mfit.ii = gsvcm.est.qr(y.train, X.train, BQ2.train, P, lambda, family, off, r.theta, eps.sigma, Cp)
    }


    W.test = as.matrix(kr(X.test, BQ2.test, byrow = TRUE))
    eta = W.test %*% as.vector(mfit.ii$theta_hat)
    ypred.ii = linkinv(eta)


    pred.error = mean((y.test - ypred.ii)^2)
    cv.error = c(cv.error, pred.error)
  }
  return(cv.error)
}
