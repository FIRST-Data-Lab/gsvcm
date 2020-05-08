#' Produces coefficient function plots for a fitted "gsvcm" object.
#'
#' This function produces the plots of the estimated coefficient function for an generalized spatially varying coefficient regression based on sample locations.
#'
#' @importFrom Triangulation TriPlot
#' @importFrom graphics points
#' @importFrom plot3D image2D
#'
#' @param mfit Fitted ``gsvcm" object.
#' \cr
#' @param gridnumber the number of grid points for plots.
#' \cr
#' @return None
#'
#' @export
#'

plot.gsvcm = function(mfit, gridnumber = 50){

  triplot = TriPlot(mfit$V, mfit$Tr)
  readline(prompt = "Press [enter] to continue")
  np = ncol(mfit$X)

  S1b = cbind(min(mfit0$V[,1]), max(mfit0$V[,1]))
  S2b = cbind(min(mfit0$V[,2]), max(mfit0$V[,2]))

  # Generate dist.
  dist = max(S1b[2]-S1b[1], S2b[2]-S2b[1])/gridnumber
  uu = seq(S1b[1], S1b[2], dist);
  vv = seq(S2b[1], S2b[2], dist);
  n1 = length(uu); n2 = length(vv)
  u = rep(uu,n2); v = rep(vv,rep(n1,n2));
  uvpop = cbind(u,v)

  xran = (range(S1b)[2]-range(S1b)[1])/10
  yran = (range(S2b)[2]-range(S2b)[1])/10

  # Generate population basis for plot
  B_pred = basis(mfit$V,mfit$Tr,mfit$d,mfit$r,uvpop)
  BQ2_pred = as.matrix(B_pred$Bi%*%B_pred$Q2)
  Ind_all=B_pred$Ind

  for(ii in 1:np){
    beta.r1.all=matrix(NA,n1*n2,1)
    beta.r1.all[Ind_all,]=BQ2_pred%*%mfit$theta_hat[,ii]
    beta.r1.mtx=matrix(beta.r1.all,n1,n2)
    image2D(uu,vv,z=matrix(beta.r1.mtx,n1,n2), xlim=c(min(S1b)-xran,max(S1b)+xran),
            ylim=c(min(S2b)-yran,max(S2b)+yran))
    readline(prompt = "Press [enter] to continue")
  }

}

