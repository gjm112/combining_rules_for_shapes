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
numPcs <- 10
K <- cov(VV)
svd_out <- svd(K)
PC_feat <- (as.matrix(VV) %*% svd_out$u[,1:numPcs])

test <- prcomp(VV, center = FALSE, scale = FALSE)



