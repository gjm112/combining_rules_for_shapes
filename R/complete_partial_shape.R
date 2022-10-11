# plot(t(partial_shape_closed_obs), xlim = c(-500,300))
# points(t(partial_shape_closed_mis))
# points(t(complete_shape_obs),col = "red")
# points(t(complete_shape_mis),col = "blue")

#source("/Users/gregorymatthews/Dropbox/shapeanalysisgit/R/calc_shape_dist_partial.R")
##################################################################
#This function takes in ONE complete donors (i.e. completed shapes) and ONE partial shape.  
#It then return the completed shape
##################################################################
#complete_shape <- t(ptsTrainList[[1]][[2]])
#partial_shape <- t(ptsTrainList[[1]][[1]][11:42,])

#write.csv(complete_shape, file = "/Users/gregorymatthews/Desktop/complete_shape_2983.csv")
#write.csv(partial_shape, file = "/Users/gregorymatthews/Desktop/partial_shape_closed.csv")


 #complete_partial_shape(complete_shape,partial_shape)

# complete_shape <- complete_shape_list[["DSCN3144"]] 
# complete_shape <- resamplecurve(complete_shape, 1000, mode = "C")
# partial_shape <- resamplecurve(partial_shape, 1000)
# plot(t(complete_shape))
# plot(t(partial_shape))

complete_partial_shape <- function(complete_shape, partial_shape, plot = FALSE, scale = FALSE){
  
  #Dimension
  d <- dim(complete_shape)[1]
  #Number of points for complete_shape and partial_shape
  N_complete <- dim(complete_shape)[2]
  N_partial <- dim(partial_shape)[2]
  
  t <- seq(0,1,length = 100)
  x0 <- matrix(NA,ncol = 100,nrow = 2)
  for (j in 1:d){
    x0[j,] <- (1-t)*partial_shape[j,N_partial] + t*partial_shape[j,1]
  }
  
  partial_shape_closed <- cbind(partial_shape,x0[,2:100])
  
  N_complete_new <- 500
  t <- seq(0,1,length = N_complete_new)
  
  olddel <- get_cumdel(partial_shape_closed)
  
  N <- 100
  #Note: resamplecurve is using splines 
  #Does this sampl ing need to 
  #partial_shape_closed <- resamplecurve(partial_shape_closed,N_complete_new)
  
  #newpt <- which(t < olddel[N_partial])
  #newpt <- newpt[length(newpt)]
  
  #Resample the complete and missing parts
  #partial_shape_closed_obs <- resamplecurve(partial_shape_closed[,1:newpt],N)
  #partial_shape_closed_mis <- resamplecurve(partial_shape_closed[,newpt:dim(partial_shape_closed)[2]],N)
  
  partial_shape_closed_obs <- resamplecurve(partial_shape_closed[,1:(dim(partial_shape_closed)[2] - (dim(x0)[2] - 1))],N)
  
  partial_shape_closed_mis <- resamplecurve(partial_shape_closed[,(dim(partial_shape_closed)[2] - (dim(x0)[2] - 1)):dim(partial_shape_closed)[2]],N)
  
  #Find the centroid of the observed part
  cent1 <- apply(partial_shape_closed_obs,1, mean)
  
  #Centering
  partial_shape_closed_obs <- partial_shape_closed_obs - cent1
  partial_shape_closed_mis <- partial_shape_closed_mis - cent1
  
  
  if (scale == TRUE){
  #scale factor
  sc1 <- norm(partial_shape_closed_obs, type = "F")
  
  #Scaling the shape
  partial_shape_closed_obs <- partial_shape_closed_obs/sc1
  partial_shape_closed_mis <- partial_shape_closed_mis/sc1
  }
  
  minE <- Inf
  #I think we are looking across all strting points around the curve?
  for (j in 0:(N_complete-1)){
    #What does shiftF do??
    #Why N_complete - 1 and not just N_complete??????
    mu <- ShiftF(complete_shape[,1:(N_complete-1)],j) 
    mu <- cbind(mu,mu[,1])
    
    olddel1 <- get_cumdel(mu)
    
    N <- 100
    library(fdasrvf)
    mu <- resamplecurve(mu,N_complete_new)
    
    newpt1 <- which(t < olddel1[N_partial]) #olddel1 changed to olddel per suggestion of sebastian. Simulations were run with this line.
    #newpt1 <- which(t < olddel[N_partial])
    newpt1 <- newpt1[length(newpt1)]
    
    mu1 <- resamplecurve(mu[,1:newpt1],N) 
    mu2 <- resamplecurve(mu[,newpt1:dim(mu)[2]],N) 
    
    cent2 <- apply(mu1,1,mean)
    mu1 <- mu1 - cent2
    mu2 <- mu2 - cent2
    
    if (scale == TRUE){
    sc2=norm(mu1, type = "F")
    mu1=mu1/sc2
    mu2=mu2/sc2
    }
    
    #Finding the best rotation
    out <- find_best_rotation(partial_shape_closed_obs,mu1)
    R <- out$R
    q2new <- out$q2new
    
    mu1n <- R%*%mu1
    
    Ec <- InnerProd_Q(partial_shape_closed_obs-mu1n,partial_shape_closed_obs-mu1n)
    # if (Ec < minE){
    #   Rbest <- R
    #   complete_shape_obs <- Rbest%*%mu1
    #   complete_shape_mis <- Rbest%*%mu2
    #   minE <- Ec
    # }
    
    if (Ec < minE){
      Rbest <- R
      complete_shape_obs <- mu1
      complete_shape_mis <- mu2
      minE <- Ec
    }
    
    
  }
  
  
  
  
  donor <- cbind(complete_shape_obs,complete_shape_mis[,2:dim(complete_shape_mis)[2]])
  
  partial_shape_closed_new <- cbind(partial_shape_closed_obs,partial_shape_closed_mis[,2:dim(partial_shape_closed_mis)[2]])
  
  
  
  
  n <- 40
  tn <- seq(0,1, length = N)
  
  #Why minus 1 in cos?  
  #Defining basic functions
  b <- matrix(NA, nrow = 2*n, ncol = length(tn))
  for (j in 1:n){
    b[j,] <- sin(2*pi*j*tn)/(sqrt(2)*pi*j)
    b[j+n,] <- (cos(2*pi*j*tn)-1)/(sqrt(2)*pi*j)
  }
  
  
  n <- dim(b)[1]
  iter <- 1
  eps <- 15 #This controls the step size in the gradient descent algorith.
  
  while(iter < 500){
    
    v <- partial_shape_closed_mis - complete_shape_mis
    
    #Computing each basis component and then adding them.  
    gradE <- matrix(0,nrow = 2,ncol = N)
    for (j in 1:n){
      for (k in 1:d){
        val <- trapz(tn,v[k,]*b[j,])*b[j,]
        gradE[k,] <- gradE[k,]+val
      }}
    
    ngE <- trapz(tn,apply(gradE*gradE,2,sum))
    
    partial_shape_closed_mis <- partial_shape_closed_mis-eps*gradE
    
    partial_shape_closed_new <- cbind(partial_shape_closed_obs,partial_shape_closed_mis[,2:dim(partial_shape_closed_mis)[2]])
    
    ngE
    
    
    iter=iter+1
    
    if (plot == TRUE & iter%%25 == 0){
    #plot(t(partial_shape_closed_new), type = "l", xlim = c(-400,400))
    plot(t(partial_shape_closed_new), type = "l")
    points(t(donor), type = "l" , col= "red")
    }
    
    
  }
  
  out <- list(partial_shape_imputed = partial_shape_closed_new, partial_obs = partial_shape_closed_obs,  donor = donor)
  return(out)
  
  
}
