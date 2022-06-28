Mstep_GtoX <- function(G, r, K, N) {
  nOmics <- length(K)
  # store multinomial logistic regression model with corresponding coefficients
  fit <- vector(mode = "list", length = nOmics)
  Beta <- vector(mode = "list", length = nOmics)
  
  # if 2 omics data
  if(nOmics == 2) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , j], margin = i)
      }))
      invisible(capture.output(temp_fit <- nnet::multinom(r_margin ~ G)))
      fit[[i]] <- temp_fit
      Beta[[i]] <- coef(temp_fit)
    }  
  }
  
  return(list(fit = fit,
              Beta = Beta))
}

Mstep_XtoZ <- function(Z, r, K, modelNames, N) {
  nOmics <- length(K)
  # store GMM model with corresponding model
  fit <- vector(mode = "list", length = nOmics)
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  
  # if 2 omics data
  if(nOmics == 2) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , j], margin = i)
      }))
      r_margin <- round(r_margin, digits = 8)
      temp_fit <- mstep(data = Z[[i]], 
                        G = K[i], 
                        z = r_margin, 
                        modelName = modelNames[i])
      fit[[i]] <- temp_fit
      Mu[[i]] <- temp_fit$parameters$mean
      Sigma[[i]] <- temp_fit$parameters$variance$sigma
    }
  }
  
  
  return(list(fit = fit,
              Mu = Mu,
              Sigma = Sigma))
}

Mstep_XtoY <- function(Y, r, K, N, family) {
  # if 2 omics data
  if(length(K) == 2) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(rowSums(r[, , i]), colSums(r[, , i]))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1)]
    
    if(family == "gaussian") {
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
    }
    
    if(family == "binomial") {
      fit <- glm(Y ~ r_fit, family = "binomial")
      mu <- as.numeric(coef(fit))
      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      # fit <- NULL
      # p <- matrix(rep(0, prod(K)),
      #             nrow = K[1],
      #             ncol = K[2])
      # for(i in 1:K[1]) {
      #   for(j in 1:K[2]) {
      #     p[i, j] <- sum(r[i, j, ] * Y) / sum(r[i, j, ])
      #   }
      # }
      fit <- fit
      mu <- p
      sd <- NULL
    }
    
    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }
  
  Delta <- list(mu = mu,
                sd = sd,
                K = K)
  return(list(fit = fit,
              Delta = Delta))
}