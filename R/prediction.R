# functions to predict 
predict_lucid <- function(model,
                          G = NULL,
                          Z = NULL,
                          Y = NULL) {
  K <- model$K
  nOmics <- length(K)
  N <- model$N
  # initialize container for predicted value
  pred_X <- vector(mode = "list", length = nOmics)
  pred_z <- vector(mode = "list", length = nOmics)
  # prediction of X and Y based on fitted data
  if(all(is.null(G), is.null(Z), is.null(Y))) {
    r <- model$z
    # 1 - prediction for X
    for (i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , j], margin = i)
      }))
      pred_X[[i]] <- map(r_margin)
      pred_z[[i]] <- r_margin
    }
    # 2 - prediction for Y
    pred_Y <- predict(model$res_Delta$fit)
  } else { # prediction of X and Y based on new data
    # predict X and Y based on G
    
    # predict X and Y based on Z
    
    # predict X and Y based on X and G
    
    # predict X and Y based 
    
  }
  
  return(list(pred_X = pred_X,
              pred_z = pred_z,
              pred_Y = pred_Y))
}