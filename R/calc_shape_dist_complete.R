#complete_shape and partial shape are both completed shapes.  
#partial shape should be in the imputed shape.  
#This function is no symmetric.  
calc_shape_dist_complete <- function(complete_shape,partial_shape,scale = FALSE){

#Dimension
d <- dim(complete_shape)[1]
#Number of points for complete_shape and partial_shape
N_complete <- dim(complete_shape)[2]
N_partial <- dim(partial_shape)[2]

N <- 100
partial_shape <- resamplecurve(partial_shape,N)

#Find the centroid of the observed part
cent1 <- apply(partial_shape,1, mean)

#Centering
partial_shape <- partial_shape - cent1


if (scale == TRUE){
  #scale factor
  sc1 <- norm(partial_shape, type = "F")
  
  #Scaling the shape
  partial_shape <- partial_shape/sc1
  
}


minE <- Inf
#I think we are looking across all strting points around the curve?
for (j in 0:(N_complete-1)){
  #What does shiftF do??
  #Why N_complete - 1 and not just N_complete??????
  mu <- ShiftF(complete_shape[,1:(N_complete-1)],j) 
  #mu <- cbind(mu,mu[,1])
  
  mu <- resamplecurve(mu,N) 
  
  cent2 <- apply(mu,1,mean)
  mu <- mu - cent2
  
  
  if (scale == TRUE){
    sc2=norm(mu, type = "F")
    mu=mu/sc2
  }
  
  #Finding the best rotation
  out <- find_best_rotation(partial_shape,mu)
  R <- out$R
  q2new <- out$q2new
  
  mu_n <- R%*%mu
  
  Ec <- InnerProd_Q(partial_shape-mu_n,partial_shape-mu_n)
  if (Ec < minE){
    Rbest <- R
    complete_shape_new <- Rbest%*%mu
    minE <- Ec
  }
  
  
}






#mu <- cbind(complete_shape_obs,complete_shape_mis[,2:dim(complete_shape_mis)[2]])

#ccompnew <- cbind(partial_shape_obs,partial_shape_mis[,2:dim(partial_shape_mis)[2]])

dist <- calc_shape_dist(complete_shape_new,partial_shape, mode = "C")  
return(dist)

}