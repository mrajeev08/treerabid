# Functions to estimate detection probabilities ----
# Per Cori et al
# adapted from vimes here
# for estimating detection probabilities

#' Estimate detection probabilities given observed distances between case pairs
#'
#' Adapted from [Cori et al. 2019]() and the R package `vimes`
#'
#' @param dist_pairs observed distances between linked case pairs (i.e. from
#'  contact tracing data + transmission tree reconstruction)
#' @param p vector of probabilities corresponding to expected probabilities
#'  of discrete sequential distances (this is still not ironed out, but only
#'  needed for empirical distributions)
#' @param pdf type of distance, options are "temporal", "spatial", "empirical", or "genetic"
#' @param fixed_pars parameters to pass to the d functions (i.e. dempiric, dgamma, dgenetic, dspatial)
#' @param alpha the value at which to cutoff the potential number of cases between two linked cases
#'  see Cori et al. for more details
#'  si_base <- distcrete::distcrete("lnorm", 1L,
#'                                  meanlog = treerabid::params_treerabid$SI_meanlog,
#'                                  sdlog = treerabid::params_treerabid$SI_sdlog)$d(1:1000)
#'
convolve_empirical <- function(p, alpha = 0.001, max_time, min_pi) {

    out <- list()
    p <- fill_with(p, 0, max_time)
    out[[1]] <- p
    max_kappa <- get_kappa(alpha, pi = min_pi)

    for (k in 2:max_kappa) {
      ps <- stats::convolve(out[[k-1]],
                            rev(p),
                            type="open")
      if(length(ps) < max_time) {
        ps <- fill_with(ps, 0, max_time)
      }

      if(length(ps) > max_time) {
        ps <- ps[1:max_time]
      }

      if(sum(ps) > alpha) {
        out[[k]] <- ps
      } else {
        break
      }
    }

    out <- do.call(cbind, out)
    out[out < 0] <- 0
    colnames(out) <- seq_len(ncol(out))

  return(out)
}

get_empirical_probs <- function(p_mat, t_diff) {

  out_probs <- p_mat[t_diff, ]

}

# sim kappa between case dates  (instead of max_gens, set an alpha level!)
est_generations <- function(out_probs,
                            kappa_weights = TRUE) {

   k_max <- ncol(out_probs)
   scaled_probs <- out_probs/rowSums(out_probs)
   gens <- apply(scaled_probs, 1, function(x) sample(k_max, 1, prob = x))

   if(kappa_weights) {
    gens <- tabulate(gens, nbins = k_max)/length(gens)
  }

  return(gens)
}

est_pi <- function(t_diff, nsims = 1000, out_probs,
                   candidate_pis) {

  sims <- lapply(seq_len(nsims),
                 function(x) {
                   est_generations(out_probs,
                                   kappa_weights = TRUE)
                             })

  candidate_weights <- lapply(candidate_pis,
                              function(x) {
                                dgeom(seq_len(ncol(out_probs)) - 1, x)
                              })
  candidate_weights <- do.call(cbind, candidate_weights)
  out <- unlist(lapply(sims, function(x) {
    candidate_pis[which.min(colSums((x - candidate_weights)^2))]
  }))

  return(out)
}

# # get the prob mat
# max_base <- qlnorm(0.999, meanlog = treerabid::params_treerabid$SI_meanlog,
#                    sdlog = treerabid::params_treerabid$SI_sdlog)
# p <- distcrete::distcrete("lnorm", 1L,
#                            meanlog = treerabid::params_treerabid$SI_meanlog,
#                            sdlog = treerabid::params_treerabid$SI_sdlog)$d(1:max_base)
# out_mat <- convolve_empirical(p, alpha = 0.001, max_time = 15000, min_pi = 0.01)
# system.time({
#  tt <- rbindlist(lapply(seq(0.01, 0.99, by = 0.05), function(z) {
#   rbindlist(lapply(seq_len(10), function(x) {
#     t_diff <- sim_times_pi(si_fun, nobs = 572, params = treerabid::params_treerabid, alpha = 0.001,
#                            pi = z)
#     out_probs <- get_empirical_probs(out_mat, t_diff)
#     ests <- est_pi(t_diff, nsims = 10, out_probs,
#                    candidate_pis = seq(0.01, 0.99, by = 0.25))
#     data.table(true = z, estimated = ests, sim = x)}))
#   }))
# })


# Helpers ----
fill_with <- function(x, filling, L = length(x)) {
  if (L <= length(x)) {
    return(x)
  }
  out <- rep(filling, L)
  out[seq_along(x)] <- x
  return(out)
}

# the probability of 1:kappa_max -1 intermediate cases having been unobserved
# and the kth case being observed
# given a reporting probability (that applies equally to cases)
fit_detection_from_kappa <- function(alpha = 0.01,
                                     weights_obs,
                                     pi) {

  max_kappa <- get_kappa(alpha, pi)

  weights <- dgeom(seq_len(max_kappa) - 1, pi)

  # fill in either weights or weigths_obs with zeroes if needed
  if(length(weights) < length(weights_obs)) {
    weights[(length(weights) + 1):length(weights_obs)] <- 0
  }

  if(length(weights) > length(weights_obs)) {
    weights_obs[(length(weights_obs) + 1):length(weights)] <- 0
  }

  # get the least squares
  return(sum((weights_obs - weights)^2))
}

# sim kappa between case dates  (instead of max_gens, set an alpha level!)
sim_generations <- function(t_diff, si_fun, params, max_kappa = 100,
                            kappa_weights = TRUE) {

  out <- matrix(si_fun(length(t_diff) * max_kappa, params), nrow = length(t_diff))

  # starting one sorted from min to max (closest poss match in high rep scenario)
  t_diff <- sort(t_diff)
  out[, 1] <- sort(out[, 1])

  # get the difference
  out_sum <- t(apply(out, 1, cumsum)) - t_diff # get the diff
  out_sum[out_sum > 0] <- -Inf
  gens <- apply(out_sum, 1, function(x) which.max(x)) # select the one before tdiff exceeded
  gens[is.na(gens)] <- 1

  if(kappa_weights) {
    gens <- tabulate(gens, nbins = max_kappa)/length(gens)
  }

  return(gens)
}

## simulate fake data
si_fun <- function(N, params) {

  rlnorm(N, meanlog = params$SI_meanlog, sdlog = params$SI_sdlog)

}

# simulate times given generation function & pi
sim_times_pi <- function(si_fun, nobs, params, alpha = 0.001, pi) {

  max_kappa <- get_kappa(alpha, pi)
  out <- matrix(si_fun(nobs * max_kappa, params), ncol = max_kappa)
  out <- t(apply(out, 1, cumsum))
  weights <- dgeom(seq_len(max_kappa) - 1, pi)
  kappas <- sample(seq_len(max_kappa), nobs, prob = weights, replace = TRUE)
  t_diff <- out[cbind(seq_len(nobs), kappas)]

  return(t_diff)
}

# fit N simulations to kappa and get minimum
fit_sims_pi <- function(t_diff, nsims = 1000,
                        candidate_pis, si_fun, params, alpha = 0.001) {

  max_max_kappa <- get_kappa(alpha, pi = min(candidate_pis))

  candidate_weights <- lapply(candidate_pis,
                              function(x) {
                                dgeom(seq_len(ncol(out_probs)) - 1, x)
                              })
  candidate_weights <- do.call(cbind, candidate_weights)

  out <-
    parallel::mclapply(
      seq_len(nsims),
      function(x) {
        weights_sim <- sim_generations(t_diff, si_fun, params,
                                       max_kappa = max_max_kappa,
                                       kappa_weights = TRUE)
        candidate_pis[which.min(colSums((weights_sim - candidate_weights)^2))]
      }
    )

    return(out)

}

get_kappa <- function(alpha, pi, min_kappa = 2) {
  pmax(qgeom(1 - alpha, pi) + 1L, min_kappa)
}

# tt <- rbindlist(lapply(seq(0.01, 0.99, by = 0.05), function(z) {
#   rbindlist(lapply(seq_len(10), function(x) {
#     t_diff <- sim_times_pi(si_fun, nobs = 572, params = treerabid::params_treerabid, alpha = 0.01,
#                            pi = z)
#     ests <- fit_sims_pi(t_diff, nsims = 10, candidate_pis = seq(0.01, 0.99, by = 0.01),
#                         si_fun, params = treerabid::params_treerabid, alpha = 0.01)
#     data.table(true = z, estimated = ests, sim = x)}))
#   }))



