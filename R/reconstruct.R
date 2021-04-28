#' Reconstruct transmission trees from a case line list.
#'
#' @param id_case id of case
#' @param id_biter id of biter (i.e. known biter from contact tracing data, if
#'  unknown should be 0, if no contact tracing data, can pass NULL)
#' @param y_coord y coordinate of case (should be in UTM: to do = use haversine distance for long/lat)
#' @param x_coord x coordinate of case (should be in UTM: to do = use haversine distance for long/lat)
#' @param owned whether animal is owned or not (per Katie, to account for uncertainty in case locations)
#' @param date_symptoms case date (i.e. date symptoms started)
#' @param days_uncertain uncertainty in days around case date
#' @param use_known_source whether to assign known progenitors from contact tracing data
#' @param prune whether to prune links (i.e. based on Cori et al./Mancy et al.)
#'  at a certain cutoff probability. This also results in assignment of incursions
#' @param si_fun the function to get the probability of a given temporal
#'  difference between cases in days (i.e. the serial interval)
#' @param dist_fun the function to get the probability of a given spatial
#' difference between cases in meters (i.e. the dispersal kernel)
#' @param cutoff the probability level at which to prune links (i.e. cases can not be
#'  linked if the distance or time difference is greater than the this %ile of the distribution)
#' @param params list of parameters to pass to si_fun and dist_fun
#'
#' @return a data.table with the reconstructed tree
#' @export
#'
build_tree <- function(id_case,
                         id_biter,
                         y_coord,
                         x_coord,
                         owned,
                         date_symptoms,
                         days_uncertain,
                         use_known_source = FALSE,
                         prune = TRUE,
                         si_fun,
                         dist_fun,
                         cutoff = 0.95,
                         params) {

  if(cutoff >= 1 | cutoff <= 0 | length(cutoff) > 1) {
    stop("Cutoff value should be a single probability between (0, 1)")
  }

  # build line list with uncertainty
  t <- date_symptoms +  add_uncertainty(days_uncertain)
  t_dt <- data.table(id_case, id_biter, x_coord, y_coord, owned, t, join = t)
  t_dt$type <- "reconstructed" # for tracking whether using known sources or not

  setkey(t_dt, join) # this speeds up join

  # Get the si cutoff if pruning (either length 1 or length nrow(t_dt))
  t_dt_on <- t_dt[, .(t_progen = t,
                      x_coord_progen = x_coord,
                      y_coord_progen = y_coord,
                      id_progen = id_case,
                      join_on = join)]

  setkey(t_dt_on, join_on)

  if(use_known_source) {

    # Build the known tree
    known_tree <- build_known_tree(t_dt,
                                   t_dt_on,
                                   si_fun,
                                   dist_fun,
                                   params)

    # Filter out of progenitor assigment (but not out of the candidate progens!)
    t_dt <- t_dt[!(id_case %in% known_tree$id_case)]

    # Filter out any invalid pairs (i.e. progen can't be an option for secondary case)
    # This shouldn't happen if uncertainty is propagated forwards?
    invalid <- paste(known_tree$id_case, known_tree$id_biter)

  }

  # do a inner join to get possible progenitors
  if(prune) {
    # get the si cutoff
    si_cutoff <- si_fun(ttree = t_dt, cutoff = cutoff, params = params)
    t_dt_on[, max := join_on + si_cutoff]

    # pruning at first stage based on time cutoffs to speed up
    ttree <- t_dt[t_dt_on,
                  on = .(join > join_on, join <= max),
                  allow.cartesian = TRUE, nomatch = NULL][, -c("join", "join.1")]
  } else {
    ttree <- t_dt[t_dt_on,
                  on = .(join > join_on),
                  allow.cartesian = TRUE, nomatch = NULL][, -c("join", "max")]
  }

  if(use_known_source) {
    # a candidate progenitor cannot be a known secondary case from contact tracing
    ttree <- ttree[!(paste(id_progen, id_case) %in% invalid)]
  }

  # For each case + progenitor combination, get the time and distance lag
  ttree[, c("dist_diff", "t_diff") := .(sqrt((x_coord - x_coord_progen)^2 + (y_coord - y_coord_progen)^2),
                                      as.numeric(t - t_progen))]

  # prune again for distance if true
  if(prune) {
    # Get the dist cutoff if pruning (either length 1 or length nrow(t_dt))
    dist_cutoff <- dist_fun(ttree = ttree, cutoff = cutoff, params = params)
    ttree <- ttree[dist_diff <= dist_cutoff]

    # and assign incursions as cases where no progenitor was identified
    # based on the distance OR time cutoffs
    incursions <- t_dt[!(t_dt$id_case %in% ttree$id_case)][, -"join"]
  } else {
    incursions <- NULL # otherwise don't return incursions because you're not assigning them
  }

  # This is actually the slow part so limiting # of possibilities speeds things up a lot
  si_fun(ttree = ttree, params = params)
  dist_fun(ttree = ttree, params = params)

  # source probability
  ttree[, source_prob := dist_prob * t_prob][, prob_scale := source_prob/sum(source_prob), by = id_case]
  ttree <- ttree[, selected := assign_progen(prob_scale), by = id_case][selected == 1] # filter to the chosen progenitor
  ttree[, prob_ll := log(source_prob)]

  ttree <- rbind(incursions, ttree, fill = TRUE)

  if(use_known_source) {  ttree <- rbind(ttree, known_tree, fill = TRUE) }
  ttree[, incursion := is.na(id_progen)]

  return(ttree)

}

#' Reconstruct known tree (i.e. from tracing data)
#'
#' @param t_dt made inside build_tree
#' @param t_dt_on made inside build_tree
#' @inheritParams build_tree
#' @param min_time if uncertainty results in negative time difference between
#'  known case pairs, set it to this value (better way to deal with propagating
#'  uncertainty given known case pairs?)
#'
#' @return a data.table with the known contact tracing tree
#' @keywords internal
#'
build_known_tree <- function(t_dt,
                             t_dt_on,
                             si_fun,
                             dist_fun,
                             params,
                             min_time = 1e-6) {

  # (as these get linked regardless and we need to maintain their dates relative to one another)
  t_dt_on[, join_ct := id_progen]

  # Filter ones with non-zero biter ids making sure the ids exist in the progenitor list
  known_biters <- t_dt[id_biter != 0 & id_biter %in% t_dt_on$id_progen]
  known_biters[, join_ct := id_biter]

  # Inner join with the biters
  known_tree <- t_dt_on[known_biters, on = "join_ct"][, -"join_ct"]
  # Also get rid of join_ct from t_dt_on (modify in place)
  t_dt_on[, join_ct := NULL]

  # Get their probs for distance and time
  known_tree[, c("dist_diff", "t_diff") := .(sqrt((x_coord - x_coord_progen)^2 + (y_coord - y_coord_progen)^2),
                                             as.numeric(t - t_progen))]

  # fix any that have negative values due to uncertainty (better way?)
  known_tree$t_diff[known_tree$t_diff <= 0] <- min_time

  # probabilities
  si_fun(ttree = known_tree, params = params)
  dist_fun(ttree = known_tree, params = params)
  known_tree[, source_prob := dist_prob * t_prob]

  # Deal with multiple id's here, selecting the one with the highest source probability
  known_tree <- known_tree[known_tree[, .I[which.max(source_prob)], by = "id_case"]$V1][, -c("join", "join_on")]

  # Set their type
  known_tree$type <- "traced"

  return(known_tree)
}

# Helper functions ----------

#' Assign progenitor
#'
#' @keywords internal
assign_progen <- function(prob_scaled) {

  cum_scaled <- cumsum(prob_scaled)
  rand_var <- runif(1, 0, 1)
  ind <- cum_scaled[cum_scaled > rand_var][1]
  ifelse(cum_scaled == ind, 1L, 0L)

}

#' Add uncertainty
#'
#' @keywords internal
add_uncertainty <- function(uncertainty){

  days_offset <- unlist(lapply(uncertainty, sample, size = 1)) * sample(c(-1, 1), length(uncertainty), replace = TRUE)

}

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

#' Helper for functions to export to foreach loop in boot_trees
#'
#' @param filename path to the R script with functions. You still need to source
#'  this function before calling boot_trees.
#'
#' @return a character vector of names functions in R script
#' @export
#'
list_funs <- function(filename) {
  temp.env <- new.env()
  sys.source(filename, envir = temp.env)
  functions <- ls(envir = temp.env)
  return(functions)
}

#' Wrapper function for bootstrapped trees: reproducibly & in parallel
#'
#' @inheritParams build_tree
#' @param N number of trees to build
#' @param seed seed to pass to doRNG to make trees reproducible
#' @param exp_funs functions to export to foreach loop (i.e. for customized
#'  si_fun/dist_fun which relies on an external scripts), may be needed
#'  for certain types of cluster configs
#' @param exp_pkgs packages to export to foreach loop, defaults to
#'  data.table and treerabid, if other dependencies for si_fun or
#'  dist_fun, then pass here. May be needed for
#'  certain types of cluster configs,
#'
#' @return a data.table with bootstrapped trees
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#' @export
#'
boot_trees <- function(id_case,
                       id_biter,
                       y_coord,
                       x_coord,
                       owned,
                       date_symptoms, # needs to be in a date class
                       days_uncertain,
                       use_known_source = FALSE,
                       prune = TRUE,
                       si_fun,
                       dist_fun,
                       cutoff = 0.95,
                       params,
                       N = 1,
                       seed = 1245,
                       exp_funs = NULL,
                       exp_pkgs = c("data.table", "treerabid")) {

  if(any(is.na(x_coord) | is.na(y_coord) | is.na(date_symptoms))) {
    stop("Missing data in times or locations!")
  }

  # check that date_symptoms is a date class
  if(class(date_symptoms) != "Date") {
    stop("date_symptoms is not of class Date!")
  }

  foreach(i = seq_len(N), .combine = 'rbind', .options.RNG = seed,
          .export = exp_funs,
          .packages = exp_pkgs) %dorng% {
    ttree <-
      build_tree(id_case = id_case, id_biter = id_biter, y_coord = y_coord,
                 x_coord = x_coord,
                 owned = owned, date_symptoms = date_symptoms,
                 days_uncertain = days_uncertain,
                 use_known_source = use_known_source,
                 prune = prune,
                 si_fun,
                 dist_fun,
                 cutoff = cutoff,
                 params = params)
      ttree$sim <- i
      ttree

  }

}
