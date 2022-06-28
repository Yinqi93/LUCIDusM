
#' EM algorithm for LUCID with multi-omics data
#'
#' @param G a matrix
#' @param Z a list, each element is a matrix
#' @param Y a vector
#' @param K a vector, each element is a integer, representing number of latent 
#' clusters corresponds to each element in Z
#' @param modelNames 
#' @param useY 
#' @param init_par 
#'
#' @return
#'
#' @examples
EM_lucid <- function(G, Z, Y, K, 
                     modelNames = rep("VVV", length(K)), 
                     useY = TRUE, 
                     init_par = NULL,
                     tol = 1e-3, 
                     family = c("gaussian", "binomial"),
                     max_itr = 1e3, 
                     seed = 123) {
  # check data format
  if(!is.matrix(G)) {
    G <- as.matrix(G)
  }
  for(i in 1:length(Z)) {
    if(!is.matrix(Z[[i]])) {
      Z[[i]] <- as.matrix(Z[[i]])
    }
  }
  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  
  
  # basic setup
  N <- nrow(G)
  nOmics <- length(Z)
  nG <- ncol(G)
  nZ <- as.integer(sapply(Z, ncol))
  family <- match.arg(family)
  
  # initialize model parameters
  set.seed(seed)
  Mu_Sigma <- initialize_Mu_Sigma(K = K, Z = Z, modelNames = modelNames)
  Mu <- Mu_Sigma$Mu
  Sigma <- Mu_Sigma$Sigma
  Beta <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    invisible(capture.output(temp_fit <- nnet::multinom(Mu_Sigma$z[[i]] ~ G)))
    Beta[[i]] <- coef(temp_fit)
  }
  # Beta <- initialize_Beta(K = K, nG = nG)
  Delta <- initialize_Delta(K = K, nCoY = 0, family = family, 
                            z = Mu_Sigma$z, Y = Y)
  loglik <- -Inf
  
  # start EM algorithm
  flag_converge <- FALSE
  itr <- 0
  while(!flag_converge & itr < max_itr) {
    # E-step
    Estep_array <- Estep(G = G, Z = Z, Y = Y, 
                         Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Delta,
                         family = family, useY = useY)
    Estep_r <- Estep_to_r(Estep_array = Estep_array,
                          K = K,
                          N = N)
    
    # M-step
    res_Beta <- Mstep_GtoX(G = G, r = Estep_r, K = K, N = N)
    res_Mu_Sigma <- Mstep_XtoZ(Z = Z, r = Estep_r, K = K, 
                               modelNames = modelNames, N = N)
    if(useY) {
      res_Delta <- Mstep_XtoY(Y = Y, r = Estep_r, K = K, N = N,
                              family = family)
    }
    
    # update parameters
    Beta <- res_Beta$Beta
    Mu <- res_Mu_Sigma$Mu
    Sigma <- res_Mu_Sigma$Sigma
    if(useY) {
      Delta <- res_Delta$Delta
    }
    
    
    
    # check convergence
    loglik_update <- cal_loglik(Estep_array = Estep_array,
                                Estep_r = Estep_r)
    if(abs(loglik - loglik_update) < tol) {
      flag_converge <- TRUE
      cat("converge!\n")
    } else {
      itr <- itr + 1
      loglik <- loglik_update
      cat(paste0("iteration ", itr, ": log-likelihood = ", loglik_update, "\n"))
    }
  }
  
  
  if(!useY) {
    res_Delta <- Mstep_XtoY(Y = Y, r = Estep_r, K = K, N = N,
                            family = family)
    Delta <- res_Delta$Delta
  }
  
  return(list(res_Beta = res_Beta,
              res_Mu_Sigma = res_Mu_Sigma,
              res_Delta = res_Delta,
              loglik = loglik_update,
              z = Estep_r,
              K = K,
              N = N))
}