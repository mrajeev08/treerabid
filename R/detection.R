# Functions to estimate detection probabilities

# 1. run trees 100 times (using known links but reconstructing unknown)
# 2. For each tree, run optim 100 times

# Per Cori et al
# adapted from vimes here
# for estimating detection probabilities

estimate_detection <- function(dist_pairs, p = NULL,
                               pdf = c("temporal", "spatial", "empirical"), 
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
  
  out <- data.table(est = out$par, 
                    convergence = out$convergence)
  return(out)
  
}

est_detection <- function(x, p = NULL,
                          alpha, 
                          pdf = c("temporal", "spatial", "empirical"), 
                          fixed_pars = list(), 
                          par = c(pi = 0.1), 
                          adj = 1e-6) {
  
  pi <- par
  x[x == 0] <- x[x == 0] + 1e-6
  if(pdf == "temporal") {
    
    loglik <- -sum(dtemporal(x, shape = fixed_pars$SI_shape, scale = fixed_pars$SI_scale, 
                             pi = pi, alpha = alpha, log = TRUE), na.rm = TRUE)
    
  }
  
  if(pdf == "spatial") {
    
    loglik <- -sum(dspatial(x, sd = fixed_pars$sd, pi = pi, log = TRUE), na.rm = TRUE)
    
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
fill_with <- function(x, filling, L = length(x)) {
  if (L <= length(x)) {
    return(x)
  }
  out <- rep(filling, L)
  out[seq_along(x)] <- x
  return(out)
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


dempiric <- function(p, pi, alpha = 0.001) {
  
  max_kappa <- qgeom(1 - alpha, pi) + 1L
  weights <- dgeom(seq_len(max_kappa) - 1, pi)
  weights <- weights / sum(weights)
  distributions <- convolve_empirical(p, max_kappa, TRUE)
  
  out <- distributions %*% weights
  return(as.vector(out))
}

# inputs should be scaled to 1
kb_stat <- function(px, py) {
  
  sum(exp(px) * (log(exp(px) / exp(py))))
  
}
