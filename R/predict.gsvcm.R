#' Making predictions from a fitted generalized spatially varying coefficient model
#'
#' This function is used to make predictions of the generalized spatially varying coefficient models.
#'
#' @importFrom BPST basis
#'
#' @param mfit Fitted ``gsvcm" object.
#' \cr
#' @param Xpred The design matrix for prediction.
#' \cr
#' @param Spred The cooridinates for prediction.
#' \cr
#' @return A vector of predicted response is returned.
#' \cr
#' @details This R package is the implementation program for manuscript entitled "Generalized Spatially Varying Coefficinet Models" by Myungjin Kim and Li Wang.
#' \cr
#' @examples
#' # See an example of fit.gsvcm.
#' @export
#'

predict.gsvcm = function(mfit, Xpred, Spred){
  if(!is.matrix(Xpred)){
    warning("The explanatory variable, Xpred, should be a matrix.")
    Xpred = as.matrix(Xpred)
  }
  if(!is.matrix(Spred)){
    warning("The coordinates, Spred, should be a matrix.")
    Spred = as.matrix(Spred)
  }

  family = mfit$family; linkinv = family$linkinv;

  if(identical(Spred, mfit$S)){
      W = as.matrix(kr(mfit$X, mfit$B %*% mfit$Q2, byrow = TRUE))
      eta = W %*% as.vector(mfit$theta_hat)
      ypred = linkinv(eta)
    } else {
      V = mfit$V; Tr = mfit$Tr; d = mfit$d; r = mfit$r; Q2 = mfit$Q2

      Basis.full = basis(V, Tr, d, r, Spred, FALSE, FALSE)
      ind.inside.pred = Basis.full$Ind.inside
      Bpred = Basis.full$B

      Xpred = as.matrix(Xpred[ind.inside.pred, ])

      W = as.matrix(kr(Xpred, Bpred %*% Q2, byrow = TRUE))
      eta = W %*% as.vector(mfit$theta_hat)
      ypred = linkinv(eta)
    }
  return(ypred)
}
