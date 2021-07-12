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
#'
#' @return an optim object
#' @export
#'
estimate_detection <- function(dist_pairs, p = NULL,
                               pdf = c("temporal", "spatial", "empirical", "genetic"),
                               fixed_pars, alpha = 0.001) {

  # optim
  out <- optim(par = runif(1, 0.1, 1),
               est_detection,
               x = dist_pairs,
               p = p,
               pdf = pdf,
               fixed_pars = fixed_pars,
               alpha = 0.001,
               method = "Brent", upper = 1, lower = 0.1)

  return(out)

}

#' Function to get likelihood of a detection probability given observed distances
#' between case pairs and an underlying distribution
#'
#' @param x distances between case pairs
#' @inheritParams estimate_detection
#' @param par starting value for detection probability
#' @param adj constraining distances to be non-zero
#'
#' @return
#' @export
#'
est_detection <- function(x, p = NULL,
                          alpha,
                          pdf = c("temporal", "spatial", "empirical", "genetic"),
                          fixed_pars = list(),
                          par = c(pi = 0.1),
                          adj = 1e-6) {

  pi <- par
  if(pdf == "temporal") {
    x[x == 0] <- x[x == 0] + 1e-6
    loglik <- -sum(dtemporal(x, shape = fixed_pars$SI_shape, scale = fixed_pars$SI_scale,
                             pi = pi, alpha = alpha, log = TRUE), na.rm = TRUE)

  }

  if(pdf == "spatial") {
    x[x == 0] <- x[x == 0] + 1e-6
    loglik <- -sum(dspatial(x, sd = fixed_pars$sd, pi = pi, alpha = alpha, log = TRUE), na.rm = TRUE)

  }

  if(pdf == "genetic") {

    loglik <- -sum(dgenetic(x, gamma_shape = fixed_pars$SI_shape,
                            gamma_scale =  fixed_pars$SI_scale,
                            poisson_rate = fixed_pars$snp_rate_daily,
                            pi = pi, alpha = alpha, log = TRUE), na.rm = TRUE)

  }

  if(pdf == "empirical") {
    if(is.null(p)) stop("empirical estimates require p (i.e. vector of probabilities)")
    y <- dempiric(p, pi, alpha = 0.001)
    maxl <- max(c(x, length(y)))
    px <- tabulate(x, maxl)/sum(maxl)
    py <- y/sum(y)
    loglik <- kb_stat(x, y)
  }

  return(loglik)

}

dtemporal <- function(x, shape, rate = 1, scale = 1/rate, pi, alpha = 0.001, log = TRUE) {

  max_kappa <- qgeom(1 - alpha, pi) + 1L
  weights <- dgeom(seq_len(max_kappa) - 1, pi)
  weights <- weights / sum(weights)
  distributions <- convolve_gamma(shape, scale = scale,
                                  kappa = max_kappa, keep_all = TRUE, log = log)(x)

  out <- distributions %*% weights
  return(as.vector(out))
}

dgenetic <- function(x, gamma_shape, gamma_rate = 1,
                     gamma_scale = 1 / gamma_rate,
                     poisson_rate,
                     pi,
                     alpha = 0.001, log = TRUE) {

  max_kappa <- qgeom(1 - alpha, pi) + 1L
  weights <- dgeom(seq_len(max_kappa) - 1, pi)
  weights <- weights / sum(weights)
  distributions <- convolve_gamma_poisson(gamma_shape,
                                          gamma_scale = gamma_scale,
                                          poisson_rate = poisson_rate,
                                          kappa = max_kappa,
                                          keep_all = TRUE,  log = log)(x)

  out <- distributions %*% weights
  return(as.vector(out))
}

dspatial <- function(x, sd, pi, alpha = 0.001, log = TRUE) {

  max_kappa <- qgeom(1 - alpha, pi) + 1L
  weights <- dgeom(seq_len(max_kappa) - 1, pi)
  weights <- weights / sum(weights)
  distributions <- convolve_spatial(sd = sd,
                                    kappa = max_kappa,
                                    keep_all = TRUE,
                                    log = log)(x)

  out <- distributions %*% weights
  return(as.vector(out))
}

convolve_gamma <- function(shape, rate = 1, scale = 1 / rate,
                           kappa, keep_all = FALSE,
                           log = TRUE) {

  if (keep_all) {
    f <- function(x) {
      out <- sapply(x, function(e)
        stats::dgamma(e, shape = (1:kappa) * shape, scale = scale, log = log)
      )

      if (is.matrix(out)) {
        return(t(out))
      } else {
        return(matrix(out))
      }
    }
  }  else {
    f <- function(x) stats::dgamma(x, shape = kappa * shape, scale = scale, log = log)
  }
  return(f)
}


convolve_gamma_poisson <- function(gamma_shape, gamma_rate = 1,
                                   gamma_scale = 1 / gamma_rate,
                                   poisson_rate, kappa, keep_all = FALSE,
                                   log = TRUE) {

  ## using prob = 1-p intead of p so that our definition correponds to that of dnbinom
  prob <- 1 - (gamma_scale * (poisson_rate / (gamma_scale * poisson_rate + 1)))

  if (keep_all) {
    f <- function(x) {
      out <- sapply(x, function(e)
        stats::dnbinom(e, size = (1:kappa) * gamma_shape, prob = prob, log = log)
      )

      if (is.matrix(out)) {
        return(t(out))
      } else {
        return(matrix(out))
      }
    }
  } else {
    f <- function(x)
      stats::dnbinom(x, size = kappa * gamma_shape, prob = prob, log = log)
  }
  return(f)
}

convolve_spatial <- function(sd, kappa, keep_all = FALSE, log = TRUE) {

  if (keep_all) {
    f <- function(x) {
      out <- sapply(x, function(e)
        VGAM::drayleigh(e, scale=sd*sqrt(1:kappa), log = log))

      if (is.matrix(out)) {
        return(t(out))
      } else {
        return(matrix(out))
      }

    }
  } else {
    f <- function(x) VGAM::drayleigh(x,scale=sd*sqrt(kappa), log = log)
  }
  return(f)
}

dempiric <- function(p, pi, alpha = 0.001) {

  max_kappa <- qgeom(1 - alpha, pi) + 1L
  weights <- dgeom(seq_len(max_kappa) - 1, pi)
  weights <- weights / sum(weights)
  distributions <- convolve_empirical(p, max_kappa, TRUE)

  out <- distributions %*% weights
  return(as.vector(out))
}


convolve_empirical <- function(x, kappa, keep_all = FALSE) {

  if (kappa == 1) {
    if (keep_all) {
      x <- matrix(x)
      colnames(x) <- "1"
    }
    return(x)
  }

  if (keep_all) {
    out <- list()
    out[[1]] <- x

    for (k in 2:kappa) {
      out[[k]] <- stats::convolve(out[[k-1]],
                                  rev(x),
                                  type="open")
    }

    L <- length(out[[kappa]])
    out <- lapply(out, fill_with, 0, L)
    out <- as.matrix(data.frame(out))
    colnames(out) <- seq_len(ncol(out))
  } else {
    out <- x

    for (k in 2:kappa) {
      out <- stats::convolve(out,
                             rev(x),
                             type="open")
    }
  }
  return(out)
}


# Helpers ----

# inputs should be scaled to 1
kb_stat <- function(px, py) {

  sum(exp(px) * (log(exp(px) / exp(py))))

}


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
  sims <- parallel::mclapply(seq_len(nsims),
                             function(x) {
                               sim_generations(t_diff, si_fun, params,
                                               max_kappa = max_max_kappa,
                                               kappa_weights = TRUE)
                              })
  unlist(parallel::mclapply(sims, function(x) {
    ss <- unlist(lapply(candidate_pis, function(z) {
        max_kappa <- get_kappa(alpha, pi = z)
        fit_detection_from_kappa(alpha, weights_obs = x, pi = z)
      }))
    candidate_pis[which.min(ss)]
  }))
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



