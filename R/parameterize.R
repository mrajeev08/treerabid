# Distribution functions -----

#' Gamma distribution for probabilities and cutoff for serial interval
#'
#' These are examples for how to write the si_fun and dist_fun functions.
#' They must take three parameters: ttree, params, and cutoff.
#'
#' If the cutoff is NULL, it will take the ttree and add a new column
#' with the probability (using the column t_diff which is generated inside
#' build_tree). If the cutoff is not NULL then it will return the value at which
#' you should prune the trees--this should either return a single value of length 1
#' or a vector of values of length nrow(ttree). See dist_gamma_mixed for an example.
#'
#' @export
#'
si_gamma1 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, t_prob := dgamma(t_diff, shape = params$SI_shape, scale = params$SI_scale)]
  } else {
    # return the cutoff value given a prob
    qgamma(cutoff, shape = params$SI_shape, scale = params$SI_scale)
  }
}

#' Convolved gamma serial interval
#'
#' See ?si_gamma1 for more details.
#'
#' @export
si_gamma2 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, t_prob := dgamma(t_diff, shape = params$SI2_shape, scale = params$SI2_scale)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qgamma(cutoff, shape = params$SI2_shape, scale = params$SI2_scale)
  }
}

#' Lognormal serial interval
#'
#' See ?si_gamma1 for more details.
#'
#' @export
si_lnorm1 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, t_prob := dlnorm(t_diff, meanlog = params$SI_meanlog, sdlog = params$SI_sdlog)]
  } else {
    # return the cutoff value given a prob
    qlnorm(cutoff, meanlog = params$SI_meanlog, sdlog = params$SI_sdlog)
  }
}

#' Convolved serial interval
#'
#' See ?si_gamma1 for more details.
#'
#' @export
si_lnorm2 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, t_prob := dlnorm(t_diff, meanlog = params$SI2_meanlog, sdlog = params$SI2_sdlog)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qlnorm(cutoff, meanlog = params$SI2_meanlog, sdlog = params$SI2_sdlog)
  }
}

#' Mixture gamma dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_gamma_mixed <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(owned & dist_diff > 100,
                               dgamma(dist_diff, shape= params$DK_shape, scale = params$DK_scale),
                               !owned & dist_diff > 100,
                               dgamma(dist_diff, shape= params$DK2_shape, scale = params$DK2_scale),
                               dist_diff <= 100,
                               pgamma(100, shape = params$DK_shape, scale = params$DK_scale)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    ifelse(ttree$owned, qgamma(cutoff, shape = params$DK_shape, scale = params$DK_scale),
           qgamma(cutoff, shape = params$DK2_shape, scale = params$DK2_scale))
  }
}

#' Mixture weibull dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_weibull_mixed <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(owned & dist_diff > 100,
                               dweibull(dist_diff, shape= params$DK_shape_weibull,
                                        scale = params$DK_scale_weibull),
                               !owned & dist_diff > 100,
                               dweibull(dist_diff, shape= params$DK2_shape_weibull,
                                        scale = params$DK2_scale_weibull),
                               dist_diff <= 100,
                               pweibull(100, shape = params$DK_shape_weibull,
                                        scale = params$DK_scale_weibull)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    ifelse(ttree$owned, qweibull(cutoff, shape = params$DK_shape_weibull,
                                 scale = params$DK_scale_weibull),
           qweibull(cutoff, shape = params$DK2_shape_weibull,
                    scale = params$DK2_scale_weibull))
  }
}
#' Mixture lognormal dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_lnorm_mixed <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(owned & dist_diff > 100,
                               dlnorm(dist_diff, meanlog =  params$DK_meanlog, sdlog = params$DK_sdlog),
                               !owned & dist_diff > 100,
                               dlnorm(dist_diff, meanlog = params$DK2_meanlog, sdlog = params$DK2_sdlog),
                               dist_diff <= 100,
                               plnorm(100, meanlog =  params$DK_meanlog, sdlog = params$DK_sdlog)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    ifelse(ttree$owned, qlnorm(cutoff, meanlog = params$DK_meanlog, sdlog = params$DK_sdlog),
           qlnorm(cutoff, meanlog = params$DK2_meanlog, sdlog = params$DK2_sdlog))
  }
}

#' Gamma dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_gamma1<- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(dist_diff > 100,
                               dgamma(dist_diff, shape= params$DK_shape, scale = params$DK_scale),
                               dist_diff <= 100,
                               pgamma(100, shape = params$DK_shape, scale = params$DK_scale)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qgamma(cutoff, shape = params$DK_shape, scale = params$DK_scale)
  }
}

#' Weibull dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_weibull1<- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(dist_diff > 100,
                               dweibull(dist_diff, shape= params$DK_shape_weibull,
                                        scale = params$DK_scale_weibull),
                               dist_diff <= 100,
                               pweibull(100, shape = params$DK_shape_weibull,
                                        scale = params$DK_scale_weibull)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qweibull(cutoff, shape = params$DK_shape_weibull, scale = params$DK_scale_weibull)
  }
}

#' Lognormal dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_lnorm1 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(dist_diff > 100,
                               dlnorm(dist_diff, meanlog =  params$DK_meanlog, sdlog = params$DK_sdlog),
                               dist_diff <= 100,
                               plnorm(100, meanlog =  params$DK_meanlog, sdlog = params$DK_sdlog)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qlnorm(cutoff, meanlog = params$DK_meanlog, sdlog = params$DK_sdlog)
  }
}

#' Convolved Gamma dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_gamma2 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(dist_diff > 100,
                               dgamma(dist_diff, shape= params$DK2_shape, scale = params$DK2_scale),
                               dist_diff <= 100,
                               pgamma(100, shape = params$DK2_shape, scale = params$DK2_scale)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qgamma(cutoff, shape = params$DK2_shape, scale = params$DK2_scale)
  }
}

#' Convolved Weibull dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_weibull2 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(dist_diff > 100,
                               dweibull(dist_diff, shape= params$DK2_shape_weibull,
                                        scale = params$DK2_scale_weibull),
                               dist_diff <= 100,
                               pweibull(100, shape = params$DK2_shape_weibull,
                                        scale = params$DK2_scale_weibull)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qweibull(cutoff, shape = params$DK2_shape_weibull, scale = params$DK2_scale_weibull)
  }
}

#' Convolved lognormal dispersal kernel
#'
#' See ?si_gamma1 for more details.
#'
#' @export
dist_lnorm2 <- function(ttree, params, cutoff = NULL) {
  if(is.null(cutoff)) {
    ttree[, dist_prob := fcase(dist_diff > 100,
                               dlnorm(dist_diff, meanlog =  params$DK2_meanlog, sdlog = params$DK2_sdlog),
                               dist_diff <= 100,
                               plnorm(100, meanlog =  params$DK2_meanlog, sdlog = params$DK2_sdlog)/100)]
  } else {
    # return the cutoff value given a prob (either length 1 or length of the ttree)
    qlnorm(cutoff, meanlog = params$DK2_meanlog, sdlog = params$DK2_sdlog)
  }
}
