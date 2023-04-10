library(fdasrvf)
library(parallel)
source("/Users/gregorymatthews/Dropbox/gladysvale/R/utility.R")
source("/Users/gregorymatthews/Dropbox/gladysvale/R/curve_functions.R")

setwd("/Users/gregorymatthews/Dropbox/shape_completion_Matthews_et_al/")
source("./R/utility.R")
source("./R/curve_functions.R")
source("./R/calc_shape_dist_partial.R")
source("./R/calc_shape_dist_complete.R")
source("./R/complete_partial_shape.R")
source("./R/impute_partial_shape.R")

#Find a set of teeth to work with.  
load("/Users/gregorymatthews/Dropbox/gladysvale/RData/teeth_BW_train_20210622.RData")
ref <- read.csv("/Users/gregorymatthews/Dropbox/gladysvale/reference_file_20210622.csv")

#Pull out Alcelaphini from LM1
ind <- unique(ref$image[ref$tribe == "Alcelaphini" & ref$type == "LM1"])
ind <- intersect(ind, names(teeth_BW_train))

#Population
#each shape is represented by 60 points as a starting point
teeth <- list()
pop <- array(NA, c(2,60,length(ind)))
for (i in ind){
  teeth[[i]] <- resamplecurve(t(teeth_BW_train[[i]]),60)
}
for (i in 1:length(teeth)){
  pop[,,i] <-teeth[[i]] 
}


#convert pop to matrix
teeth_out <- do.call(rbind,teeth)
write.csv(teeth_out, file = "/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/teeth_out.csv", row.names = FALSE)

#run our MATLAB script
system("/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r \"run('/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/compute_mean_for_combining_rules.m'); exit\"")

pop_mean_q <- read.csv("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/mean_q.csv", header = FALSE)
pop_mean_beta <- read.csv("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/mean_beta.csv", header = FALSE)

save(pop_mean_beta, file = "/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/pop_mean_beta.RData")
save(pop_mean_q, file = "/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/pop_mean_q.RData")

#Now project the population down onto the tangent space of the mean 
#run our MATLAB script
system("/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r \"run('/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/project_to_tangent_space.m'); exit\"")


VV <- read.csv("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/teeth_pop_coords_in_tan_space_of_overall_mean.csv", header = FALSE)
VV <- as.matrix(VV)
numPCs <- 25
K <- cov(VV)
svd_out <- svd(K)
PC_feat <- (as.matrix(VV) %*% svd_out$u[,1:numPcs])

#I'm going to start with 10 pcs. 
#svd_out$d - eigen values.  These are the variances.  
#rnorm() requires sd so we need the sqrt.  

set.seed(1234)
nsim <- 5
x_sim <- matrix(NA, nrow = nsim, ncol = numPCs)
for (i in 1:nsim){
  x_sim[i,] <- rnorm(numPCs, 0, sqrt(svd_out$d[1:numPCs]))
}

simulated_shapes_in_tan_space <- x_sim %*% t(svd_out$u[,1:numPCs])

n_col <- ncol(simulated_shapes_in_tan_space)

simulated_shapes_in_tan_space_list <- list()
for (i in 1:nrow(simulated_shapes_in_tan_space)){
  simulated_shapes_in_tan_space_list[[i]] <- rbind(simulated_shapes_in_tan_space[i,1:(n_col/2)],
                                                   simulated_shapes_in_tan_space[i,((n_col/2) + 1):n_col])
}

#stacked and ready to use in matlab
simulated_shapes_in_tan_space_out <- do.call(rbind,simulated_shapes_in_tan_space_list)

#Now send these back to shape space
write.csv(simulated_shapes_in_tan_space_out,file = "/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/simulated_shapes_in_tan_space.csv", row.names = FALSE)

#Now send back to shape space. 
system("/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r \"run('/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/send_sim_shapes_back_to_shape_space.m'); exit\"")

#Read these simulated shapes back in to R. 
sim_shapes_data_shapespace_beta <- read.csv("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/sim_shapes_data_shapespace_beta.csv", header = FALSE)
sim_shapes_data_shapespace_beta <- as.matrix(sim_shapes_data_shapespace_beta)

#Simulated shapes!
i <- 4
plot(sim_shapes_data_shapespace_beta[i,1:100],sim_shapes_data_shapespace_beta[i,101:200])

