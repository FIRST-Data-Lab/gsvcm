#' Estimation for GSVCMs and QRGSVCMs, which is used in R funtion "fit.gsvcm".

gsvcm_est=
function(y, X, BQ2, P, lambda, family, off = 0, theta = 0)
{
  # Environmental Variables
  variance = family$variance
  linkinv = family$linkinv
  linkfun = family$linkfun
  mu.eta = family$mu.eta
  initialize = family$initialize

  nobs = length(y)
  np = ncol(X)
  J = dim(BQ2)[2]
  Z = as.matrix(kr(X, BQ2, byrow=TRUE))
  D = as.matrix(kronecker(diag(rep(1, np)), P))
  weights = rep(1, nobs)

  nl = length(lambda)
  alpha_all = matrix(rep(0, (np * J) * nl), ncol = nl)
  gcv_all = rep(0, nl)
  df_all=c()


  for(il in 1:nl){
    mustart = y + 0.1
    etastart = linkfun(mustart)
    if(theta != 0) var = variance(mustart, theta) else var = variance(mustart)
    mevg = mu.eta(etastart)

    Y_iter = (etastart - off) + (y - mustart) / mevg
    W_iter = (mevg^2) / var
    W_iter = as.vector(W_iter)
    temp1 = as.matrix(W_iter * Z)
    temp2 = as.vector(W_iter * Y_iter)

    temp3 = t(Z) %*% temp1 + lambda[il] * D
    alpha_old = solve(temp3, t(Z) %*% temp2)
    step = 0
    delta1 = 1
    delta2 = 1

    while(delta1 > 1e-5 & sum(is.infinite(delta2)) == 0 & step <= 10){
      step = step + 1
      eta = Z %*% alpha_old + off
      mu = linkinv(eta)
      if(theta != 0) var = variance(mu, theta) else var = variance(mu)

      mevg = mu.eta(eta)
      Y_iter = (eta - off) + (y - mu) / mevg
      W_iter = as.vector((mevg^2) / var)
      temp1 = as.matrix(W_iter * Z)
      temp2 = as.vector(W_iter * Y_iter)


      temp3 = t(Z) %*% temp1 + lambda[il] * D

      alpha_new = solve(temp3, t(Z) %*% temp2)
      eta_new = Z %*% alpha_new
      delta1 = sqrt(mean((alpha_new - alpha_old)^2))
      delta2 = exp(eta_new)

      if(sum(is.infinite(delta2)) == 0){
        alpha_old = alpha_new
      }
    }
    alpha_all[,il] = as.matrix(alpha_old)

    DD = solve(temp3)
    Slambda = Z %*% tcrossprod(DD,temp1)
    df = sum(diag(Slambda))
    df_all=c(df_all,df)
    yhat = Z %*% alpha_old
    gcv = nobs * sum(W_iter * (Y_iter-yhat)^2) / (nobs - df)^2
    gcv_all[il] = gcv
  }
  j = which.min(gcv_all)
  gcv = gcv_all[j]
  df = df_all[j]
  lambdac = lambda[j]
  alpha_hat = alpha_all[,j]

  theta_hat = matrix(alpha_hat, J, np)
  beta = as.matrix(BQ2 %*% theta_hat)

  list(beta = beta, theta_hat = theta_hat, lambdac = lambdac, gcv = gcv, df = df)
}

gsvcm_est_nb =
  function(y, X, BQ2, P, lambda, family, off = 0, r.theta, eps.sigma)
  {
    # Environmental Variables
    variance = family$variance
    linkinv = family$linkinv
    linkfun = family$linkfun
    mu.eta = family$mu.eta
    initialize = family$initialize

    lower = r.theta[1]
    upper = r.theta[2]
    result1 = gsvcm_est(y, X, BQ2, P, lambda, family, off, upper)
    y_hat = linkinv(kr(X, BQ2, byrow = TRUE) %*% as.vector(result1$theta_hat))
    tt = (y - y_hat)^2 / variance(y_hat, upper)
    U_sigma = 1 / (length(y) - result1$df) * sum(tt)

    result2 = gsvcm_est(y, X, BQ2, P, lambda, family, off, lower)
    y_hat = linkinv(kr(X, BQ2, byrow = TRUE) %*% as.vector(result2$theta_hat))
    tt = (y - y_hat)^2 / variance(y_hat, lower)
    L_sigma = 1 / (length(y) - result2$df) * sum(tt)


    M_sigma = 2
    step = 1
    if ( (L_sigma - 1) > 0 & (U_sigma - 1) > 0 ) {
      M_sigma = 1
      result = result2
      middle = lower
    }

    if ( (L_sigma - 1) < 0 & (U_sigma - 1) < 0 ) {
      M_sigma = 1
      result = result1
      middle = upper
    }


    while(abs(M_sigma - 1) > eps.sigma & step <= 10){
      middle = (upper + lower) / 2
      result = gsvcm_est(y, X, BQ2, P, lambda, family, off, middle)
      y_hat = linkinv(kr(X, BQ2, byrow = TRUE) %*% as.vector(result$theta_hat))
      tt = (y - y_hat)^2 / variance(y_hat, middle)
      M_sigma = 1 / (length(y) - result$df) * sum(tt)
      if((M_sigma - 1) * (U_sigma - 1) < 0){
        lower = middle
        L_sigma = M_sigma
      } else {
        upper = middle
        U_sigma = M_sigma
      }
      print(c(middle, M_sigma))
      step = step + 1
    }
    list(beta = result$beta, theta_hat = result$theta_hat,
         lambdac = result$lambdac, gcv = result$gcv, df = result$df, theta = middle)
  }

gsvcm.est =
  function(y, X, BQ2, P, lambda, family, off = 0, r.theta = c(2, 8), eps.sigma = 0.01)
  {
    if(family$family != "nb_bps"){
      return(gsvcm_est(y, X, BQ2, P, lambda, family, off, theta = 0))
    }
    if(family$family == "nb_bps"){
      return(gsvcm_est_nb(y, X, BQ2, P, lambda, family, off, r.theta, eps.sigma))
    }
  }


gsvcm_est_qr=
  function(y, X, BQ2, P, lambda, family, off = 0, theta = 0, Cp = TRUE)
  {
    # Environmental Variables
    variance = family$variance
    linkinv = family$linkinv
    linkfun = family$linkfun
    mu.eta = family$mu.eta
    initialize = family$initialize

    nobs = length(y)
    np = ncol(X)
    J = dim(BQ2)[2]
    Z = as.matrix(kr(X, BQ2, byrow=TRUE))

    D = as.matrix(kronecker(diag(rep(1, np)), P))
    weights = rep(1, nobs)

    nl = length(lambda)
    alpha_all = matrix(rep(0, (np * J) * nl), ncol = nl)
    gcv_all = rep(0, nl)
    df_all = rep(0, nl)

    M_num = 10
    delta1 = 1
    delta2 = 1
    q = 0;

    set.seed(123)
    ind = sample(rep(1:M_num,each = length(y)/M_num))
    sub_y = lapply(split(1:length(y),ind), function(i) y[i])
    sub_X = lapply(split(1:nrow(Z),ind), function(i) Z[i,])

    while(delta1 > 1e-5 & sum(is.infinite(delta2)) == 0 & q <= 20){
      r = 0
      f = NULL
      R = NULL

      for(m in 1:M_num){
        # part (a)
        f0 = f
        R0 = R

        # part (b)
        sub_y_iter = sub_y[[m]]
        sub_X_iter = sub_X[[m]]

        # part (c)
        if (q == 0){
          xi = 0.1
          eta_old = linkfun(sub_y_iter+xi) #xi: small number
          alpha_old = 0 # initial for alpha_old
        } else {
          eta_old = as.vector(sub_X_iter%*%alpha_old)
        }

        # part (d)
        mu_i = linkinv(eta_old)
        mevg_i = mu.eta(eta_old)

        z_i = (mevg_i)^(-1)*(sub_y_iter-mu_i)+eta_old

        # NB
        if(theta != 0) var_i = variance(mu_i, theta) else var_i = variance(mu_i)
        w_i = (var_i)^(-1/2)*mevg_i

        # part (e)
        r = r + norm(as.matrix(w_i*z_i), type='2')^2

        # part (f)
        target_X_iter = rbind(R0, as.matrix(w_i*sub_X_iter))
        QR_WX_iter = qr(target_X_iter)
        Q = qr.Q(QR_WX_iter)
        R = qr.R(QR_WX_iter)
        f = t(Q)%*%as.vector(c(f0, w_i*z_i))

        Q = NULL
      }

      # Step 3:
      r.vec=r-sum((as.matrix(f))^2)

      # Step 4:

      for(il in 1:nl){
        temp = crossprod(R)+lambda[il] * D
        alpha_new = solve(temp, t(R) %*% f)

        alpha_all[,il] = as.matrix(alpha_new)

        #calculate gcv / Cp

        F_lambda = solve(temp)%*%crossprod(R)
        df = sum(diag(F_lambda))
        df_all[il]=df

        if (Cp == TRUE){
          # phi estimate
          phi =(sum((f-R%*%alpha_new)^2)+ r.vec)/(nobs-df)
          # Cp
          gcv = sum((f-R%*%alpha_new)^2)+ r.vec + 2*phi*df
        } else if (Cp == FALSE){
          gcv = nobs*(sum((f-R%*%alpha_new)^2)+r.vec) / (nobs - df)^2
        }
        gcv_all[il] = gcv
      }

      j = which.min(gcv_all)
      gcv = gcv_all[j]
      df = df_all[j]
      lambdac = lambda[j]
      alpha_hat = alpha_all[,j]
      alpha_new = as.matrix(alpha_hat)
      eta_new = R %*% alpha_new
      delta1 = sqrt(mean((alpha_new - alpha_old)^2))
      delta2 = exp(eta_new)
      if(sum(is.infinite(delta2)) == 0){
        alpha_old = alpha_new
      }

      q = q + 1
    }


    alpha_hat = alpha_new
    theta_hat = matrix(alpha_hat, J, np)
    beta = as.matrix(BQ2 %*% theta_hat)

    list(beta = beta, theta_hat = theta_hat, lambdac = lambdac, gcv = gcv, df = df)
  }


gsvcm_est_nb_qr =
  function(y, X, BQ2, P, lambda, family, off = 0, r.theta, eps.sigma, Cp = TRUE)
  {
    # Environmental Variables
    variance = family$variance
    linkinv = family$linkinv
    linkfun = family$linkfun
    mu.eta = family$mu.eta
    initialize = family$initialize

    lower = r.theta[1]
    upper = r.theta[2]
    result1 = gsvcm_est_qr(y, X, BQ2, P, lambda, family, off, upper, Cp)

    y_hat = linkinv(kr(X, BQ2, byrow = TRUE) %*% as.vector(result1$theta_hat))
    tt = (y - y_hat)^2 / variance(y_hat, upper)
    U_sigma = 1 / (length(y) - result1$df) * sum(tt)

    result2 = gsvcm_est_qr(y, X, BQ2, P, lambda, family, off, lower, Cp)
    y_hat = linkinv(kr(X, BQ2, byrow = TRUE) %*% as.vector(result2$theta_hat))
    tt = (y - y_hat)^2 / variance(y_hat, lower)
    L_sigma = 1 / (length(y) - result2$df) * sum(tt)


    M_sigma = 2
    step = 1
    if ( (L_sigma - 1) > 0 & (U_sigma - 1) > 0 ) {
      M_sigma = 1
      result = result2
      middle = lower
    }

    if ( (L_sigma - 1) < 0 & (U_sigma - 1) < 0 ) {
      M_sigma = 1
      result = result1
      middle = upper
    }


    while(abs(M_sigma - 1) > eps.sigma & step <= 10){
      middle = (upper + lower) / 2
      result = gsvcm_est_qr(y, X, BQ2, P, lambda, family, off, middle, Cp)
      y_hat = linkinv(kr(X, BQ2, byrow = TRUE) %*% as.vector(result$theta_hat))
      tt = (y - y_hat)^2 / variance(y_hat, middle)
      M_sigma = 1 / (length(y) - result$df) * sum(tt)
      if((M_sigma - 1) * (U_sigma - 1) < 0){
        lower = middle
        L_sigma = M_sigma
      } else {
        upper = middle
        U_sigma = M_sigma
      }
      print(c(middle, M_sigma))
      step = step + 1
    }
    list(beta = result$beta, theta_hat = result$theta_hat,
         lambdac = result$lambdac, gcv = result$gcv, df = result$df, theta = middle)
  }

gsvcm.est.qr =
  function(y, X, BQ2, P, lambda, family, off = 0, r.theta = c(2, 8), eps.sigma = 0.01, Cp = TRUE)
  {
    if(family$family != "nb_bps"){
      return(gsvcm_est_qr(y, X, BQ2, P, lambda, family, off, theta = 0, Cp))
    }
    if(family$family == "nb_bps"){
      return(gsvcm_est_nb_qr(y, X, BQ2, P, lambda, family, off, r.theta, eps.sigma, Cp))
    }
  }
