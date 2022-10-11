repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  if (is.null(mx)){
    mx = 1
    nx = length(X)
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }else {
    nx = dim(X)[2]
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  
  return(mat)
}

ndims <- function(x){
  return(length(dim(x)))
}


pvecnorm <-function(v,p=2){
  sum(abs(v)^p)^(1/p)
}

pvecnorm2 <-function(dt,x){
  sqrt(sum(abs(x)*abs(x))*dt)
}

get_cumdel <- function(X){
  n <- dim(X)[1]
  T <- dim(X)[2]
  del <- c()
  del[1] <- 0
  for (r in 2:T){
    del[r] <- norm(as.matrix(X[,r]) - as.matrix(X[,r-1]), type = "F")
  }
  cumdel <- cumsum(del)/sum(del)
  cumdel
}


ShiftF <- function(vec, rotations){
  vec_index <- 1:ncol(vec)
  vec_index <- vec_index - rotations
  vec <- cbind(vec[,sign(vec_index) == 1], vec[,sign(vec_index) != 1])
  vec
}

find_best_rotation <- function(q1, q2){
  n <- dim(q1)[1]
  T <- dim(q1)[2]
  
  A <- q1%*%t(q2)
  svdA <- svd(A)
  #If A is pos def made the diagonl elements all 1?
  if (det(A) > 0){
    svdA$d <- rep(1,n)
  }
  else {
    #don't know what is happening here?  
    #Make the last element negative?  
    svdA$d <- rep(1,n)
    svdA$d[length(svdA$d)] <- -svdA$d[length(svdA$d)]
  }
  R <- svdA$u%*%diag(svdA$d)%*%svdA$v
  q2new <- R%*%q2
  return(list(R = R, q2new = q2new))
}


library(caTools)
#Finds the inner product of functions.  
InnerProd_Q <- function(q1, q2){
  n <-  dim(q1)[1]
  T <-  dim(q1)[2]
  val <- trapz(seq(0,1,length = T),apply(q1*q2,2,sum))
  return(val)
}


trapz <- function(x,y,dims=1){
  if ((dims-1)>0){
    perm = c(dims:max(ndims(y),dims), 1:(dims-1))
  } else {
    perm = c(dims:max(ndims(y),dims))
  }
  
  if (ndims(y) == 0){
    m = 1
  } else {
    if (length(x) != dim(y)[dims])
      stop('Dimension Mismatch')
    y = aperm(y, perm)
    m = nrow(y)
  }
  
  if (m==1){
    M = length(y)
    out = sum(diff(x)*(y[-M]+y[-1])/2)
  } else {
    slice1 = y[as.vector(outer(1:(m-1), dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+')) ]
    dim(slice1) = c(m-1, length(slice1)/(m-1))
    slice2 = y[as.vector(outer(2:m, dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+'))]
    dim(slice2) = c(m-1, length(slice2)/(m-1))
    out = t(diff(x)) %*% (slice1+slice2)/2.
    siz = dim(y)
    siz[1] = 1
    out = array(out, siz)
    perm2 = rep(0, length(perm))
    perm2[perm] = 1:length(perm)
    out = aperm(out, perm2)
    ind = which(dim(out) != 1)
    out = array(out, dim(out)[ind])
  }
  
  return(out)
}


circshift <- function(a, sz) {
  if (is.null(a)) return(a)
  
  if (is.vector(a) && length(sz) == 1) {
    n <- length(a)
    s <- sz %% n
    a <- a[(1:n-s-1) %% n + 1]
    
  } else if (is.matrix(a) && length(sz) == 2) {
    n <- nrow(a); m <- ncol(a)
    s1 <- sz[1] %% n
    s2 <- sz[2] %% m
    a <- a[(1:n-s1-1) %% n + 1, (1:m-s2-1) %% m + 1]
  } else
    stop("Length of 'sz' must be equal to the no. of dimensions of 'a'.")
  
  return(a)
}
