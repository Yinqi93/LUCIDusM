# log-sum-exp trick
LogSumExp <- function(vec) {
  max_vec <- max(vec)
  trick <- max_vec + log(sum(exp(vec - max_vec)))
  return(trick)
}

# K should be >= 2
check_K <- function(K) {
  for(x in K) {
    if(x != as.integer(x)) {
      stop("K should be a vector of integer")
    }
  }
  if(min(K) < 2) {
    stop("each element in K should be greater or equal than 2")
  }
}


# initialize Beta
initialize_Beta <- function(K, nG) {
  nOmics <- length(K)
  Beta <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # row represent cluster, column represents variable
    Beta[[i]] <- matrix(runif((nG + 1) * (K[i] - 1), min = -1, max = 1),
                        nrow = K[i] - 1)
  }
  return(Beta)
}

# initialize Mu
initialize_Mu <- function(K, nZ) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # column represents cluster, row represents variable
    Mu[[i]] <- matrix(runif(K[i] * nZ[i], min = -1, max = 1),
                      nrow = nZ[i])
  }
  return(Mu)
}


# initialize Sigma
initialize_Sigma <- function(K, nZ) {
  nOmics <- length(K)
  Sigma <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # each element is an array of nZ[i] x nZ[i] x K[i]
    Sigma[[i]] <- array(0, dim = c(nZ[i], nZ[i], K[i]))
    for(j in 1:K[i]) {
      Sigma[[i]][, , j] <- diag(nZ[i])
    }
  }
  return(Sigma)
}



# initialize Mu and Sigma
initialize_Mu_Sigma <- function(K, Z, modelNames) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  z <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    temp_fit <- Mclust(data = Z[[i]],
                       G = K[i],
                       modelNames = modelNames[i])
    Mu[[i]] <- temp_fit$parameters$mean
    Sigma[[i]] <- temp_fit$parameters$variance$sigma
    z[[i]] <- temp_fit$z
  }
  return(list(Mu = Mu,
              Sigma = Sigma,
              z = z))
}

# initialize Delta
# for normal outcome, Delta is a list, with beta + sigma
# for binary outcome, Delta is a vector
initialize_Delta <- function(K, nCoY = 0, family = c("gaussian", "binomial"),
                             z, Y) {
  family <- match.arg(family)
  if(family == "gaussian") {
    # if 2 omics data
    r_matrix <- cbind(z[[1]], z[[2]])
    r_fit <- r_matrix[, -c(1, K[1] + 1)]
    fit <- lm(Y ~ r_fit)
    mu <- as.numeric(coef(fit))
    sd <- sd(resid(fit))
    x <- list(mu = mu,
              sd = sd,
              K = K)
  }
  if(family == "binomial") {
    # if 2 omics data
    r_matrix <- cbind(z[[1]], z[[2]])
    r_fit <- r_matrix[, -c(1, K[1] + 1)]
    fit <- glm(Y ~ r_fit, family = "binomial")
    b <- as.numeric(coef(fit))
    # x <- list(mu = b,
    #           K = K)
    b_array <- vec_to_array(K = K, mu = b)
    p <- 1 / (1 + exp(-b_array))
    x <- list(mu = p,
              K = K)
  }
  return(x)
}


# indicator function
indicator <- function(x) {
  m <- 0
  if(x > 1) {
    m <- 1
  }
  return(m)
}

# transform mu to an array, each element corresponds to a mean for a combination
# of clusters
vec_to_array <- function(K, mu) {
  res <- array(data = rep(0, prod(K)),
               dim = K)
  # if nK = 2, transform mu to a matrix
  if(length(K) == 2) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        res[i, j] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1]
      }
    }
  }
  # if nK = 3, transform mu to a 3d array
  if(length(K) == 3) {
    
  }
  return(res)
}
