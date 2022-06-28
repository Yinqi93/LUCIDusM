# use bic to conduct model selection
tune_lucid <- function(K_list, G, Z, Y, ...) {
  # change K_list into a matrix
  K_matrix <- as.matrix(expand.grid(K_list))
  if(min(K_matrix) < 2) {
    stop("minimum K should be 2")
  }
  bic <- rep(0, nrow(K_matrix))
  model_list <- vector(mode = "list", length = nrow(K_matrix))
  for(i in 1:nrow(K_matrix)) {
    model_list[[i]] <- EM_lucid(G = G, Z = Z, Y = Y,
                                K = K_matrix[i, ],
                                ...)
    bic[i] <- cal_bic(model_list[[i]])
  }
  model_opt_index <- which(bic == min(bic))
  K_matrix <- cbind(K_matrix, bic)
  return(list(tune_K = K_matrix,
              model_list = model_list,
              model_opt = model_list[[model_opt_index]]))
}