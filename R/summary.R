# rearrange cluster order
#
# for continuous outcome - use the cluster combination corresponding to smallest
# mean as the reference cluster
get_ref_cluster <- function(Delta) {
  K <- Delta$K
  mu <- Delta$mu
  mu_matrix <- vec_to_array(K = K, mu = mu)
  ref_index <- which(mu_matrix == min(mu_matrix))
  ref <- arrayInd(ref_index, .dim = K)
  return(ref)
}


# re-arrange parameters for Delta
reorder_Delta <- function(ref, Delta) {
  K <- Delta$K
  mu_matrix <- vec_to_array(K = K, mu = Delta$mu)
  mu <- mu_matrix[ref]

  # if 2 omics data
  if(length(K) == 2) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1] - mu_matrix[ref[1], 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i] - mu_matrix[1, ref[2]]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order)
  }

  Delta$mu <- mu
  return(list(Delta = Delta,
              K_order = K_order))
}



reorder_Mu_Sigma <- function(Mu_Sigma, K_order) {
  for(i in 1:length(K_order)) {
    temp_Mu <- Mu_Sigma$Mu[[i]]
    temp_Sigma <- Mu_Sigma$Sigma[[i]]
    # reorder Mu
    Mu_Sigma$Mu[[i]] <- temp_Mu[, K_order[[i]]]
    Mu_Sigma$Sigma[[i]] <- temp_Sigma[, , K_order[[i]]]
  }
  return(Mu_Sigma)
}



reorder_Beta <- function(Beta, K_order) {
  for(i in 1:length(K_order)) {
    temp_Beta <- Beta[[i]]
    temp_Beta <- rbind(rep(0, ncol(temp_Beta)),
                       temp_Beta)
    temp_Beta_reorder <- temp_Beta[K_order[[i]], ]
    ref <- temp_Beta_reorder[1, ]
    for(j in 1:nrow(temp_Beta_reorder)) {
      temp_Beta_reorder[j, ] <- temp_Beta_reorder[j, ] - ref
    }
    Beta[[i]] <- temp_Beta_reorder[-1, ]
  }

  return(Beta)
}



reorder_z <- function(z, K_order) {
  if(length(K_order) == 2) {
    z <- z[K_order[[1]], K_order[[2]], ]
  }
  return(z)
}


#' function to reorder all model parameters
#'
#' @param model A model returned by EM_lucid
#'
#' @return A LUCID model reordered by effect size of outcome
#' @export
#'
reorder_lucid <- function(model) {
  ref <- get_ref_cluster(Delta = model$res_Delta$Delta)
  r_Delta <- reorder_Delta(ref = ref,
                           Delta = model$res_Delta$Delta)
  r_Mu_Sigma <- reorder_Mu_Sigma(model$res_Mu_Sigma,
                                 K_order = r_Delta$K_order)
  r_Beta <- reorder_Beta(Beta = model$res_Beta$Beta,
                         K_order = r_Delta$K_order)
  model$res_Delta$Delta <- r_Delta
  model$res_Mu_Sigma$Mu <- r_Mu_Sigma$Mu
  model$res_Mu_Sigma$Sigma <- r_Mu_Sigma$Sigma
  model$res_Beta$Beta <- r_Beta
  model$z <- reorder_z(model$z, K_order = r_Delta$K_order)
  return(model)
}


# function to calculate BIC
cal_bic <- function(model) {
  nOmics <- length(model$K)
  Beta <- model$res_Beta$Beta
  Mu <- model$res_Mu_Sigma$Mu
  Sigma <- model$res_Mu_Sigma$Sigma
  Delta <- model$res_Delta$Delta

  # calculate number of parameters
  n_Beta <- sum(sapply(1:nOmics, function(i) {
    length(Beta[[i]])
  }))
  n_Mu <- sum(sapply(1:nOmics, function(i) {
    length(Mu[[i]])
  }))
  n_Sigma <- sum(sapply(1:nOmics, function(i) {
    length(Sigma[[i]])
  }))
  n_Delta <- length(Delta$mu) + 1
  n_par <- n_Beta + n_Mu + n_Sigma + n_Delta

  # calculate BIC
  bic <- -2 * model$loglik + n_par * log(model$N)
  return(bic)
}
