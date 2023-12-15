# ==========================================================================================================================
# A and B functions for affine yield R(t,T) = -1/(T-t)*[ln A(t,T) - rt*B(t,T)]
# ==========================================================================================================================

# Functions to compute CIR2-implied yields

A_tT = function(mat, kappa, lambda, theta, sigma) {
  
  h = sqrt((kappa-lambda)^2 + 2*sigma^2)
  A = ((2*h*exp((kappa-lambda+h)*mat/2))/(2*h+(kappa-lambda+h)*(exp(mat*h)-1)))^(2*kappa*theta/sigma^2)
  return(A)
  
}

vA_tT = Vectorize(A_tT)

B_tT = function(mat, kappa, lambda, theta, sigma) {
  
  h = sqrt((kappa-lambda)^2 + 2*sigma^2)
  B = (2*(exp(mat*h)-1))/(2*h+(kappa-lambda+h)*(exp(mat*h)-1))
  return(B)
  
}

vB_tT = Vectorize(B_tT)

# ==========================================================================================================================
# KF single run for the general form of filter (but var-covar state-dependent)
# ==========================================================================================================================

# This is Kalman recursion modified in two ways:
# 1) negative filtered states are replaced with 0
# 2) var-covar of the shock in the transition equation is state-dependent
# transition equation is y(t) = a + Phi*y(t-1) + v(t)
# obseravtion equation is R(t) = A + B*y(t) + e(t)
# shocks v(t) are uncorrelated with observation errors e(t)
# the var-covar of v(t) is Q(t), diagonal, that depends on the filtered y(t-1) here
# v(t)[i,i] = w0[i]+w1[i]*y(t-1)[i], where w0 and w1 are implied by the CIR2 model

# A single run of the Kalman filter: filtered states and negative log-likelihood
KF_run = function(Rt, A, B, a, Phi, U, y0, Q0, w0, w1, Nobs, Nfcs, dt) {
  
  # define a matrix to collect filtered states
  y = matrix(c(y0, rep(0, Nobs*Nfcs)), nrow = Nfcs)
  
  # define an array to collect var-covar of filtered states
  Sigma = array(NA, c(Nfcs, Nfcs, Nobs+1))
  Sigma[,,1] = Q0  
  
  # define an array to collect state-dependent Q
  Q = Sigma
  
  # keep track of Log-Likelihood
  LL = 0
  
  # define a Kalman recursion
  for (t in 2:(Nobs+1)) {
    
    # project the state
    y[,t] = a+Phi%*%y[,t-1]
    
    # calculate var-covar of state innovations that depend on y[,t-1]
    Q[,,t] = diag(w0+w1*y[,t-1])
    
    # project var-covar of the state
    #Sigma[,,t] = Phi%*%Sigma[,,t-1]%*%t(Phi)+Q[,,t]
    Sigma[,,t] = tcrossprod(tcrossprod(Phi, Sigma[,,t-1]), Phi)+Q[,,t]
    
    # calculate innovation in the observation relative to expectation
    u = Rt[t-1,]-(A+B%*%y[,t])
    
    #track the var-covar of this innovation
    #Omega = B%*%Sigma[,,t]%*%t(B)+U
    Omega = tcrossprod(tcrossprod(B, Sigma[,,t]), B) + U
    invOmega = solve(Omega)
    
    # calculate Kalman gain
    #Kt = Sigma[,,t]%*%t(B)%*%invOmega
    Kt = tcrossprod(Sigma[,,t],B)%*%invOmega
    
    # update the state
    y[,t] = y[,t]+Kt%*%u
    
    # replace negative predcited state with zeroes
    y[y[,t]<0, t] = 0
    
    # update the var-covar of the state
    Sigma[,,t] = (diag(Nfcs)-Kt%*%B)%*%Sigma[,,t]

    # update negative log-likelihood (skip  M/2*log(2*pi) part)
    LL = LL + 1/2*log(det(Omega)) + 1/2*crossprod(u,invOmega)%*%u
    
  }
  
  return(list(y = y, Sigma = Sigma, LL = LL))
  
  
}

# ==========================================================================================================================
# A single run of the Kalman filter given CIR2 model parameters
# ==========================================================================================================================

# Coordinates of parameter vector:
# 1: kappa_1
# 2: theta_1
# 3: lambda_1
# 4: sigma_1
# 5: kappa_2
# 6: theta_2
# 7: lambda_2
# 8: sigma_2
# 9: alphat

CIR2_KF_run = function(pars, data, delta, tenors, Nobs, Nfcs, Nylds, me, smoother) {
  
  # define SS matrices and other SS model parameters
  
  a_in = matrix(c(pars[1]*pars[2]/(pars[1]-pars[3])*(1-exp(-(pars[1]-pars[3])*delta)),
                  pars[5]*pars[6]/(pars[5]-pars[7])*(1-exp(-(pars[5]-pars[7])*delta))),
                nrow = Nfcs, ncol = 1)
  
  Phi_in = matrix(c(exp(-(pars[1]-pars[3])*delta), 0, 0, exp(-(pars[5]-pars[7])*delta)),
                  nrow = Nfcs)
  
  w0_in = c(pars[4]^2*pars[1]*pars[2]/(2*(pars[1]-pars[3])^2)*(1-exp(-(pars[1]-pars[3])*delta))^2,
            pars[8]^2*pars[5]*pars[6]/(2*(pars[5]-pars[7])^2)*(1-exp(-(pars[5]-pars[7])*delta))^2)
  w1_in = c(pars[4]^2/(pars[1]-pars[3])*(exp(-(pars[1]-pars[3])*delta)-exp(-2*(pars[1]-pars[3])*delta)),
            pars[8]^2/(pars[5]-pars[7])*(exp(-(pars[5]-pars[7])*delta)-exp(-2*(pars[5]-pars[7])*delta)))
  
  y0_in = matrix(c(pars[1]*pars[2]/(pars[1]-pars[3]),
                   pars[5]*pars[6]/(pars[5]-pars[7])), 
                 nrow = Nfcs, ncol = 1)
  Q0_in = matrix(c(pars[4]^2*pars[1]*pars[2]/(2*(pars[1]-pars[3])^2),0,
                   0,pars[8]^2*pars[5]*pars[6]/(2*(pars[5]-pars[7])^2)),
                 nrow = Nfcs)
  
  A_in = matrix(-log(vA_tT(mat = tenors, kappa = pars[1], lambda = pars[3], theta = pars[2],
                           sigma = pars[4]))/tenors-
                  log(vA_tT(mat = tenors, kappa = pars[5], lambda = pars[7], theta = pars[6],
                            sigma = pars[8]))/tenors+
                  pars[9], nrow = Nylds)
  B_in = matrix(c(vB_tT(mat = tenors, kappa = pars[1], lambda = pars[3], theta = pars[2],
                        sigma = pars[4])/tenors,
                  vB_tT(mat = tenors, kappa = pars[5], lambda = pars[7], theta = pars[6],
                        sigma = pars[8])/tenors),
                nrow = Nylds)
  
  # run the filter for given parameters
  kfrun = KF_run(Rt = data, A = A_in, B = B_in, a = a_in, Phi = Phi_in, y0 = y0_in, Nobs = Nobs, Nfcs = Nfcs,
                 Q0 = Q0_in, w0 = w0_in, w1 = w1_in, dt = delta, U = diag(me^2, Nylds))
  
  if (smoother == F) {
    
    return(kfrun)
    
  } else {
    
    ksrun = KF_run(Rt = data[Nobs:1,],
                   A = A_in, B = B_in, a = a_in, Phi = Phi_in,
                   y0 = matrix(kfrun$y[, Nobs+1], nrow = Nfcs, ncol = 1), Nobs = Nobs, Nfcs = Nfcs, 
                   Q0 = Q0_in, w0 = w0_in, w1 = w1_in, dt = delta, U = diag(me^2, Nylds))
    ksrun$y = ksrun$y[, (Nobs+1):1]
    ksrun$Sigma = ksrun$Sigma[,, (Nobs+1):1]
    
    return(ksrun)
    
  }
  
}

# ==========================================================================================================================
# A single run of the Kalman filter given CIR2 model parameters (!!!) CPP version of the Kalman loop
# ==========================================================================================================================

# Coordinates of parameter vector:
# 1: kappa_1
# 2: theta_1
# 3: lambda_1
# 4: sigma_1
# 5: kappa_2
# 6: theta_2
# 7: lambda_2
# 8: sigma_2
# 9: alphat

CIR2_KF_run_cpp = function(pars, data, delta, tenors, Nobs, Nfcs, Nylds, me, smoother) {
  
  # define SS matrices and other SS model parameters
  
  a_in = matrix(c(pars[1]*pars[2]/(pars[1]-pars[3])*(1-exp(-(pars[1]-pars[3])*delta)),
                  pars[5]*pars[6]/(pars[5]-pars[7])*(1-exp(-(pars[5]-pars[7])*delta))),
                nrow = Nfcs, ncol = 1)
  
  Phi_in = matrix(c(exp(-(pars[1]-pars[3])*delta), 0, 0, exp(-(pars[5]-pars[7])*delta)),
                  nrow = Nfcs)
  
  w0_in = c(pars[4]^2*pars[1]*pars[2]/(2*(pars[1]-pars[3])^2)*(1-exp(-(pars[1]-pars[3])*delta))^2,
            pars[8]^2*pars[5]*pars[6]/(2*(pars[5]-pars[7])^2)*(1-exp(-(pars[5]-pars[7])*delta))^2)
  w1_in = c(pars[4]^2/(pars[1]-pars[3])*(exp(-(pars[1]-pars[3])*delta)-exp(-2*(pars[1]-pars[3])*delta)),
            pars[8]^2/(pars[5]-pars[7])*(exp(-(pars[5]-pars[7])*delta)-exp(-2*(pars[5]-pars[7])*delta)))
  
  y0_in = matrix(c(pars[1]*pars[2]/(pars[1]-pars[3]),
                   pars[5]*pars[6]/(pars[5]-pars[7])), 
                 nrow = Nfcs, ncol = 1)
  Q0_in = matrix(c(pars[4]^2*pars[1]*pars[2]/(2*(pars[1]-pars[3])^2),0,
                   0,pars[8]^2*pars[5]*pars[6]/(2*(pars[5]-pars[7])^2)),
                 nrow = Nfcs)
  
  A_in = matrix(-log(vA_tT(mat = tenors, kappa = pars[1], lambda = pars[3], theta = pars[2],
                           sigma = pars[4]))/tenors-
                  log(vA_tT(mat = tenors, kappa = pars[5], lambda = pars[7], theta = pars[6],
                            sigma = pars[8]))/tenors+
                  pars[9], nrow = Nylds)
  B_in = matrix(c(vB_tT(mat = tenors, kappa = pars[1], lambda = pars[3], theta = pars[2],
                        sigma = pars[4])/tenors,
                  vB_tT(mat = tenors, kappa = pars[5], lambda = pars[7], theta = pars[6],
                        sigma = pars[8])/tenors),
                nrow = Nylds)
  
  # run the filter for given parameters
  kfrun = KF_run_cpp(Rt = data, A = A_in, B = B_in, a = a_in, Phi = Phi_in, y0 = y0_in, Nobs = Nobs, Nfcs = Nfcs,
                     Q0 = Q0_in, w0 = w0_in, w1 = w1_in, dt = delta, U = diag(me^2, Nylds))
  
  if (smoother == F) {
    
    return(kfrun)
    
  } else {
    
    ksrun = KF_run_cpp(Rt = data[Nobs:1,],
                       A = A_in, B = B_in, a = a_in, Phi = Phi_in,
                       y0 = matrix(kfrun$y[, Nobs+1], nrow = Nfcs, ncol = 1), Nobs = Nobs, Nfcs = Nfcs, 
                       Q0 = Q0_in, w0 = w0_in, w1 = w1_in, dt = delta, U = diag(me^2, Nylds))
    ksrun$y = ksrun$y[, (Nobs+1):1]
    ksrun$Sigma = ksrun$Sigma[,, (Nobs+1):1]
    
    return(ksrun)
    
  }
  
}

# ==========================================================================================================================
# function to compute a dirty price of a coupon bond as a function of YTM
# ==========================================================================================================================

## coup = coup in decimals
## coup_dates = vector of payment dates (last one = maturity)
## coup_freq = number of coupons per year (must be consistent with coupon dates = check yourself)

cb_price = function(coup, coup_dates, coup_freq, ytm, notional) {
  
  sum((1+ytm/coup_freq)^(-coup_freq*coup_dates)*coup/coup_freq*notional) +
    (1+ytm/coup_freq)^(-coup_freq*max(coup_dates))*notional
  
}

# ==========================================================================================================================
# function to compute yield to maturity as a function of a dirty price
# ==========================================================================================================================

ytm = function(coup, coup_dates, coup_freq, price, notional) {
  
  A = try(uniroot(function(z) {
    
    cb_price(coup = coup, coup_dates = coup_dates, coup_freq = coup_freq, notional = notional, ytm = z) - price

  }, interval = c(-2,3))$root, silent = T)
  
  if (is(A, 'try-error')) {A = -9999}
  
  return(A)
  
}

# ==========================================================================================================================
# function to compute yield to maturity as a function of a dirty price
# ==========================================================================================================================

ytm_subcpp = function(coup, coup_dates, coup_freq, price, notional) {
  
  A = try(uniroot(function(z) {
    
    cb_price_full_cpp(coup = coup, coup_dates = coup_dates, coup_freq = coup_freq, notional = notional, ytm = z) - price
    
  }, interval = c(-1,2))$root, silent = T)
  
  if (is(A, 'try-error')) {A = -9999}
  
  return(A)
  
}

# ==========================================================================================================================
# Prices of 1) ZC riskless bond in a CIR2 model, 2) ZC risky bond in a Duffee (99) model
# ==========================================================================================================================

PZC_CIR2_rf = function(mat, s1, s2, alpha, kappa1, lambda1, theta1, sigma1, kappa2, lambda2, theta2, sigma2) {
  
  exp(-alpha*mat + 
        log(A_tT(mat = mat, kappa = kappa1, theta = theta1, lambda = lambda1, sigma = sigma1)) +
        log(A_tT(mat = mat, kappa = kappa2, theta = theta2, lambda = lambda2, sigma = sigma2))  -
        s1*B_tT(mat = mat, kappa = kappa1, theta = theta1, lambda = lambda1, sigma = sigma1) -
        s2*B_tT(mat = mat, kappa = kappa2, theta = theta2, lambda = lambda2, sigma = sigma2))
  
}

PZC_D99_def = function(mat, s1, s2, sh, alpha, kappa1, lambda1, theta1, sigma1, kappa2, lambda2, theta2, sigma2,
                         alpha_def, s1bar, s2bar, beta1, beta2, kappah, lambdah, thetah, sigmah, recov) {
  
  V0 = exp(-(alpha+alpha_def-beta1*s1bar - beta2*s2bar)*mat + 
           log(A_tT(mat = mat, kappa = kappa1, theta = theta1*(1+beta1), lambda = lambda1, sigma = sigma1*sqrt(1+beta1))) +
           log(A_tT(mat = mat, kappa = kappa2, theta = theta2*(1+beta2), lambda = lambda2, sigma = sigma2*sqrt(1+beta2))) +
           log(A_tT(mat = mat, kappa = kappah, theta = thetah, lambda = lambdah, sigma = sigmah))-
           s1*B_tT(mat = mat, kappa = kappa1, theta = theta1*(1+beta1), lambda = lambda1, sigma = sigma1*sqrt(1+beta1)) -
           s2*B_tT(mat = mat, kappa = kappa2, theta = theta2*(1+beta2), lambda = lambda2, sigma = sigma2*sqrt(1+beta2)) -
           sh*B_tT(mat = mat, kappa = kappah, theta = thetah, lambda = lambdah, sigma = sigmah))
  
  P =   exp(-alpha*mat + 
              log(A_tT(mat = mat, kappa = kappa1, theta = theta1, lambda = lambda1, sigma = sigma1)) +
              log(A_tT(mat = mat, kappa = kappa2, theta = theta2, lambda = lambda2, sigma = sigma2))  -
              s1*B_tT(mat = mat, kappa = kappa1, theta = theta1, lambda = lambda1, sigma = sigma1) -
              s2*B_tT(mat = mat, kappa = kappa2, theta = theta2, lambda = lambda2, sigma = sigma2))
  
  return(recov*P + (1-recov)*V0)
  
  
}

vPZC_D99_def = Vectorize(PZC_D99_def)

# ==========================================================================================================================
# YTM of a coupon risky bond in a Duffee (99) model as a function of the default state + coupon dates + CIR2 states + 
# the rest in pars_s (parameters of the CIR2), pars_h (parameters of default int), or pars_bonds (bond characteristics)
# ==========================================================================================================================

# parameters of the YTM function (what one needs to compute YTM beside the default intensity state)
# [A] pars_bonds: a list where each element (=bond) contains coupon, coupon frequency, coupon dates, notional and recovery
# 1. coup
# 2. coup_freq
# 3. c_dates
# 4. notional
# 5. recov
# 6. accrued interest
# [B] pars_s:
# 1. alpha
# 2. kappa1
# 3. lambda1
# 4. theta1
# 5. sigma1
# 6. kappa2
# 7. lambda2
# 8. theta2
# 9. sigma2
# 10. s1bar
# 11. s2bar
# [C] pars_h
# 1. alpha_def
# 2. beta1
# 3. beta2
# 4. kappah
# 5. lambdah
# 6. thetah
# 7. sigmah

# Mind that s_1, s_2, and c_dates within pars_bonds are full-sample vectors/lists and only t_ind observation is taken
# this is needed to avoid updating the parameter vector within an EKF later; just supply a different time index

#YTM_D99_def = function(s_def, s_1, s_2, pars_s, pars_h, pars_bonds, t_ind) {
#  
#  sapply(pars_bonds, function(x) {
#    
#    V = vPZC_D99_def(mat = x[[3]][[t_ind]], s1 = s_1[t_ind], s2 = s_2[t_ind], sh = s_def,
#                     alpha = pars_s[1], kappa1 = pars_s[2], lambda1 = pars_s[3], theta1 = pars_s[4], sigma1 = pars_s[5],
#                     kappa2 = pars_s[6], lambda2 = pars_s[7], theta2 = pars_s[8], sigma2 = pars_s[9],
#                     s1bar = pars_s[10], s2bar = pars_s[11],
#                     alpha_def = pars_h[1], beta1 = pars_h[2], beta2 = pars_h[3], kappah = pars_h[4],
#                     lambdah = pars_h[5], thetah = pars_h[6], sigmah = pars_h[7],
#                     recov = x[[5]])
#    
#    P_coup = sum(V)*x[[1]]/x[[2]]+V[length(V)]
#    ytm_coup = ytm(coup = x[[1]], coup_dates = x[[3]][[t_ind]], coup_freq = x[[2]], price = P_coup,
#                   notional = x[[4]])
#    return(ytm_coup)
#    
#  })
#  
#}

YTM_D99_def_single = function(s_def, s_1, s_2, pars_s, pars_h, pars_bonds, t_ind) {
  
    V = vPZC_D99_def(mat = pars_bonds[[3]][[t_ind]], s1 = s_1[t_ind], s2 = s_2[t_ind], sh = s_def,
                     kappa1 = pars_s[2], lambda1 = pars_s[3], theta1 = pars_s[4], sigma1 = pars_s[5],
                     kappa2 = pars_s[6], lambda2 = pars_s[7], theta2 = pars_s[8], sigma2 = pars_s[9],
                     alpha = pars_s[1], s1bar = pars_s[10], s2bar = pars_s[11],
                     alpha_def = pars_h[1], beta1 = pars_h[2], beta2 = pars_h[3],
                     kappah = pars_h[4], lambdah = pars_h[5], thetah = pars_h[6], sigmah = pars_h[7],
                     recov = pars_bonds[[5]])
    
    P_coup = sum(V*pars_bonds[[1]][[t_ind]])
    ytm_coup = ytm(coup = pars_bonds[[1]], coup_freq = pars_bonds[[2]], coup_dates = pars_bonds[[3]][[t_ind]], 
                   notional = pars_bonds[[4]], price = P_coup)
    return(ytm_coup)
  
}

YTM_D99_def = function(pb, s_def, s_1, s_2, pars_s, pars_h, t_ind) {
  
  mapply(YTM_D99_def_single, pars_bonds = pb,
         MoreArgs = list(s_def=s_def, s_1=s_1, s_2=s_2, pars_s=pars_s, pars_h=pars_h, t_ind=t_ind))
  
}

Price_D99_def_single = function(s_def, s_1, s_2, pars_s, pars_h, pars_bonds, t_ind) {
  
  V = vPZC_D99_def(mat = pars_bonds[[3]][[t_ind]], s1 = s_1[t_ind], s2 = s_2[t_ind], sh = s_def,
                   kappa1 = pars_s[2], lambda1 = pars_s[3], theta1 = pars_s[4], sigma1 = pars_s[5],
                   kappa2 = pars_s[6], lambda2 = pars_s[7], theta2 = pars_s[8], sigma2 = pars_s[9],
                   alpha = pars_s[1], s1bar = pars_s[10], s2bar = pars_s[11],
                   alpha_def = pars_h[1], beta1 = pars_h[2], beta2 = pars_h[3],
                   kappah = pars_h[4], lambdah = pars_h[5], thetah = pars_h[6], sigmah = pars_h[7],
                   recov = pars_bonds[[5]])
  
  P_coup = sum(V*pars_bonds[[1]][[t_ind]]) - pars_bonds[[6]][t_ind]
#  ytm_coup = ytm(coup = pars_bonds[[1]], coup_freq = pars_bonds[[2]], coup_dates = pars_bonds[[3]][[t_ind]], 
#                 notional = pars_bonds[[4]], price = P_coup)
  return(P_coup)
  
}

Price_D99_def = function(pb, s_def, s_1, s_2, pars_s, pars_h, t_ind) {
  
  mapply(Price_D99_def_single, pars_bonds = pb,
         MoreArgs = list(s_def=s_def, s_1=s_1, s_2=s_2, pars_s=pars_s, pars_h=pars_h, t_ind=t_ind))
  
}



# ==========================================================================================================================
# Single EKF run for the general form of filter (but var-covar state-dependent)
# ==========================================================================================================================

# This is an Extended Kalman Filter recursion for a 1-factor CIR with irregularly spaced obseravtions
# as in Duffee (1999)
# 
# transition equation is y(t) = a + Phi*y(t-1) + v(t) [i.e. same as in the CIR2 model], y is def. intensity
# obseravtion equation is R(t) = H(y(t)) + e(t), i.e., there's a non-linear dependence of YTM on default
# FUN = H(y) is linearized around a filtered state y(t|t-1)
# U is the measurement error matrix

# a, Phi, w0, and w1 need to be time-dependent (if there're missing obs.) hence are arrays with the 3rd dimension = time

# A single run of the Kalman filter: filtered states and negative log-likelihood
EKF_run = function(Rt, a, Phi, y0, w0, w1, FUN, U, dt, Nobs, Nfcs, Nbonds) {

  # define a matrix to collect filtered states
  y = matrix(c(y0, rep(0, Nobs*Nfcs)), nrow = Nfcs)
  
  # define an array to collect var-covar of filtered states
  Sigma = array(NA, c(Nfcs, Nfcs, Nobs+1))
  
  # adjust starting values
  linh = jacobian(f = FUN, var = y[,1], params = list(ind = 1))
  y[,1] = 1/Phi[,,1]*(t(linh)%*%(Rt[1,] - FUN(y[,1], ind = 1) + linh*y[,1])/(t(linh)%*%linh) - a[,,1])
  Sigma[,,1] = U*(t(linh)%*%linh)*(1/Phi[,,1])^2 + as.matrix(w0[,,1]+w1[,,1]*y[,1])*(1/Phi[,,1])^2
  
  # define an array to collect state-dependent Q
  Q = Sigma
  
  # keep track of Log-Likelihood
  LL = 0
  
  # define a Kalman recursion
  for (t in 2:(Nobs+1)) {
    
    # project the state
    y[,t] = a[,,t-1]+Phi[,,t-1]%*%y[,t-1]
    
    # calculate var-covar of state innovations that depends on y[,t-1]
    Q[,,t] = as.matrix(w0[,,t-1]+w1[,,t-1]*y[,t-1])
    
    # project var-covar of the state
    Sigma[,,t] = Phi[,,t-1]%*%Sigma[,,t-1]%*%t(Phi[,,t-1])+Q[,,t]
    
    # calculate innovation in the observation relative to expectation
    u = Rt[t-1,]-FUN(y[,t], ind = t-1)
    
    # linearize the measurement
    B = jacobian(f = FUN, var = y[,t], params = list(ind = t-1))
    
    #track the var-covar of this innovation
    Omega = t(B)%*%B*Sigma[,,t]+U
    
    # calculate Kalman gain
    Kt = Sigma[,,t]*solve(Omega)%*%t(B)
    
    # update the state
    y[,t] = y[,t]+Kt%*%u
    
    # replace negative predcited state with zeroes
    y[y[,t]<0, t] = 0
    
    # update the var-covar of the state
    Sigma[,,t] = (diag(Nfcs)-Kt%*%B)%*%Sigma[,,t]
    
    # update negative log-likelihood (skip  M/2*log(2*pi) part)
    LL = LL + 1/2*log(det(Omega)) + 1/2*t(u)%*%u*solve(Omega)
    
  }
  
  return(list(y = y, Sigma = Sigma, LL = LL))
  
  
}

# ==========================================================================================================================
# A single run of the EKF given Duffee 99 model parameters
# ==========================================================================================================================

D99_EKF_run = function(Price_func, pars_s, pars_h, pars_b, data, days_ind, delta, Nobs, Nfcs, Nbonds, me) {
  
  # define SS matrices and other SS model parameters
  a_in = array(pars_h[4]*pars_h[6]/(pars_h[4]-pars_h[5])*(1-exp(-(pars_h[4]-pars_h[5])*delta*days_ind)),
               dim = c(1,1, Nobs))
  
  Phi_in = array(exp(-(pars_h[4]-pars_h[5])*delta*days_ind), dim = c(1,1, Nobs))

  w0_in = array(pars_h[7]^2*pars_h[4]*pars_h[6]/(2*(pars_h[4]-pars_h[5])^2)*
                 (1-exp(-(pars_h[4]-pars_h[5])*delta*days_ind))^2,
                dim = c(1,1, Nobs)) 
  
  w1_in = array(pars_h[7]^2/(pars_h[4]-pars_h[5])*(exp(-(pars_h[4]-pars_h[5])*delta*days_ind)-
                                                     exp(-2*(pars_h[4]-pars_h[5])*delta*days_ind)),
                dim = c(1,1, Nobs)) 
  
  y0_in = matrix(pars_h[4]*pars_h[6]/(pars_h[4]-pars_h[5]), nrow = 1, ncol = 1)
  
  # run the filter for given parameters
  kfrun = EKF_run(Rt = data, a = a_in, Phi = Phi_in, y0 = y0_in, w0 = w0_in, w1 = w1_in, FUN = Price_func,
                  U = as.matrix(me^2), dt = delta, Nobs = Nobs, Nfcs = Nfcs, Nbonds = Nbonds)
  
  return(kfrun)
  
}