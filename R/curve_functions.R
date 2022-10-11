# load("/Users/gregorymatthews/Dropbox/shapeanalysisgit/data/ptsTrainList.RData")
# ref_file <- read.csv("/Users/gregorymatthews/Dropbox/shapeanalysisgit/refFile.csv")
# source("https://raw.githubusercontent.com/cran/fdasrvf/master/R/curve_functions.R")
# source("https://raw.githubusercontent.com/cran/fdasrvf/master/R/utility_functions.R")

#This version of curve_functions.R removes scaling from the distance measures.  

resamplecurve <- function (x, N = 100, mode = "O") 
{
  n = nrow(x)
  T1 = ncol(x)
  xn = matrix(0, n, N)
  delta = rep(0, T1)
  for (r in 2:T1) {
    delta[r] = pvecnorm(x[, r] - x[, r - 1], 2)
  }
  cumdel = cumsum(delta)/sum(delta)
  newdel = seq(0, 1, length.out = N)
  for (r in 1:n) {
    xn[r, ] = spline(cumdel, x[r, ], xout = newdel)$y
  }
  
  #Removed project curve.  
  #if (mode == "C") {
    #q = curve_to_q(xn)
    #qn = project_curve(q)
    #xn = q_to_curve(qn)
  #}
  return(xn)
}


curve_to_q <- function (beta, scale = F) 
{
  n = nrow(beta)
  T1 = ncol(beta)
  v = apply(beta, 1, gradient, 1/(T1 - 1))
  v = t(v)
  q = matrix(0, n, T1)
  for (i in 1:T1) {
    L = sqrt(pvecnorm(v[, i], 2))
    if (L > 1e-04) {
      q[, i] = v[, i]/L
    }
    else {
      q[, i] = v[, i] * 1e-04
    }
  }
  
  if (scale == TRUE){
    q = q/sqrt(innerprod_q2(q, q))
    }
  
  return(q)
}

find_basis_normal <- function(q){
  n = nrow(q)
  T1 = ncol(q)
  
  f1 = matrix(0,n,T1)
  f2 = matrix(0,n,T1)
  for (i in 1:T1){
    f1[,i] = q[1,i]*q[,i]/pvecnorm(q[,i])+c(pvecnorm(q[,i]),0)
    f2[,i] = q[2,i]*q[,i]/pvecnorm(q[,i])+c(0,pvecnorm(q[,i]))
  }
  h3 = f1
  h4 = f2
  integrandb3 = rep(0,T1)
  integrandb4 = rep(0,T1)
  for (i in 1:T1){
    integrandb3[i] = t(q[,i])%*%h3[,i]
    integrandb4[i] = t(q[,i])%*%h4[,i]
  }
  b3 = h3 - q*trapz(seq(0,1,length.out=T1),integrandb3)
  b4 = h4 - q*trapz(seq(0,1,length.out=T1),integrandb4)
  
  basis = list(b3,b4)
  
  return(basis)
}

project_curve <- function(q, scale = F){
  T1 = ncol(q)
  n = nrow(q)
  if (n==2){
    dt = 0.35
  }
  
  if (n==3) {
    dt = 0.2
  }
  
  epsilon =- 1e-6
  
  e = diag(1,n)
  iter = 1
  res = rep(1,n)
  J = matrix(0,n,n)
  s = seq(0,1,length.out=T1)
  qnorm = rep(0,T1)
  G = rep(0,n)
  C = rep(0,301)
  
  qnew = q
  if (scale == TRUE){
  qnew = qnew / sqrt(innerprod_q2(qnew,qnew))
  }
  
  while (pvecnorm(res) > epsilon){
    if (iter > 300){
      break
    }
    
    # Compute Jacobian
    for (i in 1:n){
      for (j in 1:n){
        J[i,j]  = 3 * trapz(s, qnew[i,]*qnew[j,])
      }
    }
    J = J + diag(1,n)
    
    for (i in 1:T1){
      qnorm[i] = pvecnorm(qnew[,i])
    }
    
    # Compute the residue
    for (i in 1:n){
      G[i] = trapz(s,qnew[i,]*qnorm)
    }
    
    res = -1*G
    
    if (pvecnorm(res)<epsilon)
      break
    
    x = solve(J,res)
    C[iter] = pvecnorm(res)
    
    delG = find_basis_normal(qnew)
    tmp = 0
    for (i in 1:n){
      tmp = tmp + x[i]*delG[[i]]*dt
    }
    qnew = qnew + tmp
    
    iter = iter + 1
  }
  
  if (scale == TRUE){
    qnew = qnew / sqrt(innerprod_q2(qnew,qnew))
  }
  
  
  return(qnew)
}

#library(fdasrvf)
# beta1 <- t(ptsTrainList[["LM1"]][[1]])
# beta2 <- t(ptsTrainList[["LM1"]][[2]])
# 
# 
# beta1 <- partial_shape_closed_obs
# beta2 <- complete_shape_obs

calculatecentroid <- function(beta){
  n = nrow(beta)
  T1 = ncol(beta)
  
  betadot = apply(beta,1,gradient,1.0/(T1-1))
  betadot = t(betadot)
  
  normbetadot = apply(betadot,2,pvecnorm,2)
  integrand = matrix(0, n, T1)
  for (i in 1:T1){
    integrand[,i] = beta[,i] * normbetadot[i]
  }
  
  scale = trapz(seq(0,1,length.out=T1), normbetadot)
  centroid = trapz(seq(0,1,length.out=T1), integrand, 2)/scale
  
  return(centroid)
}

# #A Test
# plot(ptsTrainList[[1]][[1]])
# points(ptsTrainList[[1]][[2]],col = "red", type = "l")
# points(ptsTrainList[[1]][[100]],col = "red", type = "l")
# 
# beta1 <- t(ptsTrainList[[1]][[1]])
# beta2 <- t(ptsTrainList[[1]][[2]])
# 
# 
# plot(scale(ptsTrainList[["LM1"]][["DSCN3732"]], center = TRUE, scale = FALSE), xlim = c(-500,500), asp = 1, col = "red", type = "l")
# points(scale(ptsTrainList[["LM1"]][["DSCN2616"]], center = TRUE, scale = FALSE), col = "blue", type = "l")
# points(t(results_list[["DSCN2871"]]$imputed_partial_shape[[1]][[1]]), col = "orange", type = "l")
# 
# part <- results_list[["DSCN2871"]]$imputed_partial_shape[[1]][[1]]
# whole <- t(ptsTrainList[["LM1"]][["DSCN3732"]])
# whole <- resamplecurve(whole, N = 199, mode = "C")
# calc_shape_dist(part,whole,mode = "C")
# 
# points(t(whole), col = "red", type = "l")
# 
# part <- results_list[["DSCN2871"]]$imputed_partial_shape[[1]][[1]]
# whole <- t(ptsTrainList[["LM1"]][["DSCN2616"]])
# whole <- resamplecurve(whole, N = 199, mode = "C")
# calc_shape_dist(part,whole,mode = "C")
# 
# points(t(whole), col = "blue", type = "l")
# 
# beta1 <- results_list[["DSCN2871"]]$imputed_partial_shape[[1]][[1]]
# beta2 <- resamplecurve(t(ptsTrainList[["LM1"]][["DSCN2616"]]),N=199)
# beta3 <- resamplecurve(t(ptsTrainList[["LM1"]][["DSCN3732"]]),N=199)
# 
# centroid2 = calculatecentroid(beta3)
# dim(centroid2) = c(length(centroid2),1)
# beta3 = beta3 - repmat(centroid2, 1, T1)
# 
# plot(t(beta1), col = "orange", xlim = c(-200,200), ylim = c(-200,200), asp = 1)
# points(t(beta2),col = "blue")
# points(t(beta3),col = "red")
# points(t(beta6),col = "green")




inverse_exp_coord <- function(beta1, beta2, mode = "O", rotated = T, scale = FALSE){
  T1 = ncol(beta1)
  centroid1 = calculatecentroid(beta1)
  dim(centroid1) = c(length(centroid1),1)
  beta1 = beta1 - repmat(centroid1, 1, T1)
  centroid2 = calculatecentroid(beta2)
  dim(centroid2) = c(length(centroid2),1)
  beta2 = beta2 - repmat(centroid2, 1, T1)
  
  q1 = curve_to_q(beta1)
  
  if (mode=="C"){
    isclosed = TRUE
  }
  
  # Iteratively optimize over SO(n) x Gamma using old DP
  out = reparam_curve(beta1, beta2, rotated=rotated, isclosed=isclosed, mode=mode)
  
  if (mode=="C"){
    beta2 = shift_f(beta2, out$tau)
  }
  
  
  beta2 = out$R %*% beta2
  gamI = invertGamma(out$gam)
  
  
  
  #This needs to only be run for open curves?.
  if (mode == "O") {
  beta2 = group_action_by_gamma_coord(beta2, gamI)
  }
  
  if (rotated){
    out = find_rotation_seed_coord(beta1, beta2, mode = mode)
    q2n = curve_to_q(out$beta2new)
  } else {
    q2n = curve_to_q(beta2)
  }

  
  #I have removed this step.  
  #It's causing errors and I don't think it's really necessary.
  # if (mode=="C"){
  #  q2n = project_curve(q2n)
  # }
  
  # Compute geodesic distance
  if (scale == FALSE){
    dist  = innerprod_q2(q1-q2n, q1-q2n)
  }
  
  if (scale == TRUE){
    q1dotq2 = innerprod_q2(q1, q2n)
    if (q1dotq2>1){
    q1dotq2 = 1.
  }
    dist = acos(q1dotq2)
  }
  
   
  
  
  
  # u = q2n - q1dotq2 * q1
  # normu = sqrt(innerprod_q2(u,u))
  # 
  # if (normu > 1e-4){
  #   v = u*acos(q1dotq2)/normu
  # } else {
  #   v = matrix(0, nrow(beta1), T1)
  # }
  
  return(dist=dist)
}


#Always need to speciify the mode. 
calc_shape_dist <- function (beta1, beta2, mode = "O", scale = FALSE) 
{
  out = inverse_exp_coord(beta1, beta2, mode, scale = scale)
  return(out)
}


group_action_by_gamma <- function(q, gamma){
  n = nrow(q)
  T1 = ncol(q)
  gammadot = gradient(gamma, 1.0/T1)
  qn = matrix(0, n, T1)
  timet = seq(0, 1, length.out = T1)
  
  for (j in 1:n){
    qn[j,] = spline(timet, q[j,], xout=gamma)$y * sqrt(gammadot)
  }
  
  #This line is removed to prevent scaling.  
  #qn = qn/sqrt(innerprod_q2(qn,qn))
  
  return(qn)
}

group_action_by_gamma_coord <- function(f, gamma){
  n = nrow(f)
  T1 = ncol(f)
  fn = matrix(0, n, T1)
  timet = seq(0, 1, length.out = T1)
  
  for (j in 1:n){
    fn[j,] = spline(timet, f[j,], xout=gamma)$y
  }
  
  return(fn)
}


find_rotation_seed_coord <- function(beta1, beta2, mode="O"){
  n = nrow(beta1)
  T1 = ncol(beta1)
  q1 = curve_to_q(beta1)
  Ltwo = rep(0,T1)
  Rlist = array(0,c(n,n,T1))
  #I don't need these lines if mode is "O"
  for (ctr in 1:T1){
    beta2n = shift_f(beta2, ctr)
    out = find_best_rotation(beta1, beta2n)
    q2new = curve_to_q(out$q2new)
    Ltwo[ctr] = innerprod_q2(q1-q2new,q1-q2new)
    Rlist[,,ctr] = out$R
  }
  
  tau = which.min(Ltwo)
  O_hat = Rlist[,,tau]
  if (mode=="C")
    beta2new = shift_f(beta2,tau)
  else
    beta2new = beta2
  
  beta2new = O_hat %*% beta2new
  
  return(list(beta2new=beta2new,O_hat=O_hat,tau=tau))
}


find_rotation_and_seed_q <- function(q1,q2){
  n = nrow(q1)
  T1 = ncol(q1)
  Ltwo = rep(0,T1)
  Rlist = array(0,c(n,n,T1))
  for (ctr in 1:T1){
    q2n = shift_f(q2,ctr)
    out = find_best_rotation(q1,q2n)
    Ltwo[ctr] = innerprod_q2(q1-out$q2new,q1-out$q2new)
    Rlist[,,ctr] = out$R
  }
  
  tau = which.min(Ltwo)
  O_hat = Rlist[,,tau]
  q2new = shift_f(q2,tau)
  q2new = O_hat %*% q2new
  
  return(list(q2new=q2new,O_hat=O_hat,tau=tau))
}


# library(fdasrvf)
# beta1 <- t(ptsTrainList[["LM1"]][[1]])
# beta2 <- 3*t(ptsTrainList[["LM1"]][[2]])
# 
# calc_shape_dist(beta1,beta2, mode = "C")


shift_f <- function(f, tau){
  n = nrow(f)
  T1 = ncol(f)
  fn = matrix(0, n, T1)
  fn[,1:(T1-1)] = circshift(f[,1:(T1-1)], c(0,tau))
  fn[,T1] = fn[,1]
  
  return(fn)
}


innerprod_q2 <- function(q1, q2){
  T1 = ncol(q1)
  val = sum(sum(q1*q2))/T1
  
  return(val)
}




project_tangent <- function(w, q, basis){
  w = w - innerprod_q2(w,q)*q
  bo = gram_schmidt(basis)
  
  wproj = w - innerprod_q2(w, bo[[1]])*bo[[1]] - innerprod_q2(w,bo[[2]])*bo[[2]]
  
  return(wproj)
}