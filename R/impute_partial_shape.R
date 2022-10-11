
######################################
#This function takes in a list of potential donors and one partial shape 
#It then returns M completed shapes
######################################

impute_partial_shape <- function(complete_shape_list, partial_shape, k = 10, M = 5, scale = FALSE){
  
  # dist_vec <- c(NA, length(complete_shape_list))
  # for (i in 1:length(complete_shape_list)){print(i)
  #   dist_vec[i] <- calc_shape_dist_partial(complete_shape_list[[i]],partial_shape)
  # }
  #set.seed(1234, kind = "L'Ecuyer-CMRG")
  #dist_vec <- unlist(mclapply(complete_shape_list,calc_shape_dist_partial,partial_shape = partial_shape, scale = scale, mc.cores = 20))
  dist_vec <- unlist(lapply(complete_shape_list,calc_shape_dist_partial,partial_shape = partial_shape, scale = scale))

  #plot(t(complete_shape_list[["DSCN4687"]]), ylim = c(0,2000))
  #points(t(partial_shape), col = "red")
  
  # library(foreach)
  # library(doParallel)
  # 
  # #setup parallel backend to use many processors
  # cores=detectCores()
  # cl <- makeCluster(cores[1]-1) #not to overload your computer
  # registerDoParallel(cl)
  # dist_vec <- foreach(i=1:length(complete_shape_list), .combine=c) %dopar% {
  #   source("./R/calc_shape_dist_partial.R")
  #   calc_shape_dist_partial(complete_shape_list[[i]],partial_shape)
  # }
  # #stop cluster
  # stopCluster(cl)
  
  potential_donor_indices <- c(1:length(dist_vec))[order(dist_vec)][1:k]
  #subset(ref_file, ref %in% names(dist_vec)[potential_donor_indices])
  #Randomly pick M of these with replacement
  donor_indices <- sample(potential_donor_indices, M, replace = TRUE)
  
  
  #I only need to impute the unique indices and then i'll repeat them when the data gets output
  imputed_partial_list <- list()
  for (q in unique(donor_indices)){print(q)
    imp <- complete_partial_shape(complete_shape_list[[q]],partial_shape, scale = scale)
    imputed_partial_list[[q]] <- imp$partial_shape_imputed
  }
  
  imputed_partial_list_out <- list()
  for (j in 1:M){
    imputed_partial_list_out[[j]] <- imputed_partial_list[[donor_indices[j]]]
  }
  
  out <- list(imputed = imputed_partial_list_out, partial_obs = imp$partial_obs, dist_vec = dist_vec)
  return(out)
  
}