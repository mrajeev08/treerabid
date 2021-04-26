# Reconstruct known tree (i.e. from tracing data) ----
traced_tree <- function(id_case,
                        id_biter, 
                        y_coord,
                        x_coord,
                        owned, 
                        date_symptoms,
                        si_pdist = c("gamma", "lnorm"), 
                        dist_pdist = c("gamma", "lnorm"), 
                        params = list()) {
  
  # Get the distributions
  dist_fun <- get(paste0("dist_", dist_pdist))
  si_fun <- get(paste0("si_", si_pdist))
  
  # build line list without uncertainty 
  # (as these get linked regardless and we need to maintain their dates relative to one another)
  t_dt <- data.table(id_case, id_biter, x_coord, y_coord, owned, t = date_symptoms, join = id_biter)
  
  t_dt_on <- t_dt[, .(t_progen = t,
                      x_coord_progen = x_coord,
                      y_coord_progen = y_coord,
                      id_progen = id_case, 
                      join = id_case)]
  
  # Filter ones with non-zero biter ids making sure the ids exist in the progenitor list
  known_biters <- t_dt[id_biter != 0 & id_biter %in% t_dt_on$id_progen]
  
  # Inner join with the biters
  known_tree <- t_dt_on[known_biters, on = "join"]
  
  # Get their probs for distance and time
  known_tree[, c("dist_diff", "t_diff") := .(sqrt((x_coord - x_coord_progen)^2 + (y_coord - y_coord_progen)^2), 
                                             as.numeric(t - t_progen))]
  known_tree[, c("dist_prob", "t_prob") := 
               .(dist_fun(dist_diff, pars = params, owned), 
                 si_fun(t_diff, pars = params))]
  known_tree[, source_prob := dist_prob * t_prob]
  
  # Deal with multiple id's here, selecting the one with the highest source probability
  known_tree <- known_tree[known_tree[, .I[which.max(source_prob)], by = "id_case"]$V1]
  
  # Set their type
  known_tree$type <- "traced"
  
  return(known_tree)
}


# Reconstruct transmission trees ----
faster_trees <- function(id_case,
                         id_biter, 
                         y_coord,
                         x_coord,
                         owned, 
                         date_symptoms,
                         days_uncertain,
                         use_known_source = FALSE, 
                         prune = TRUE,
                         si_pdist = c("gamma", "lnorm"), 
                         dist_pdist = c("gamma", "lnorm"), 
                         cutoff = 0.95, 
                         params = list(), 
                         known_tree = NULL) {
  

  # Get the distributions
  dist_fun <- get(paste0("dist_", dist_pdist))
  si_fun <- get(paste0("si_", si_pdist))
  
  # Get the cutoff vals (corresponds to the input function used to generate links, make more flexible?)
  si_cutoff <- get(paste0("si_cutoff_", si_pdist))(cutoff = cutoff, pars = params)
  dist_cutoff <- get(paste0("dist_cutoff_", dist_pdist))(cutoff = cutoff, pars = params)
  
  # build line list with uncertainty
  t <- date_symptoms +  add_uncertainty(days_uncertain)
  t_dt <- data.table(id_case, id_biter, x_coord, y_coord, owned, t, join = t)
  t_dt$type <- "reconstructed" # for tracking whether using known sources or not
  
  setkey(t_dt, join) # this speeds up join
  
  t_dt_on <- t_dt[, .(t_progen = t,
                      x_coord_progen = x_coord,
                      y_coord_progen = y_coord,
                      id_progen = id_case, 
                      join_on = join, 
                      max = join + si_cutoff)]

  setkey(t_dt_on, join_on)
  
  if(use_known_source) {
    # Filter out of progenitor assigment (but not out of the candidate progens!)
    t_dt <- t_dt[!(id_case %in% known_tree$id_case)]
    invalid <- known_tree[, .(id_progen = id_case, id_case = id_biter, exclude = TRUE)]
    
    # generate uncertainty
    
    # update probs based on distribution being used
    known_tree[, c("dist_prob", "t_prob") := 
                 .(dist_fun(dist_diff, pars = params, owned), 
                   si_fun(t_diff, pars = params))]
    known_tree[, source_prob := dist_prob * t_prob]
  }
  
  # do a inner join to get possible progenitors
  if(prune) {
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
    ttree <- invalid[ttree, on = c("id_case", "id_progen")][is.na(exclude)][, -c("exclude")]
  }
  
  # For each case + progenitor combination, get the time and distance lag
  ttree[, c("dist_diff", "t_diff") := .(sqrt((x_coord - x_coord_progen)^2 + (y_coord - y_coord_progen)^2), 
                                      as.numeric(t - t_progen))]
  
  # prune again for distance if true
  # and assign incursions as cases where no progenitor was identified based on the distance OR time cutoff
  if(prune) {
    ttree <- ttree[dist_diff <= dist_cutoff]
    incursions <- t_dt[!(t_dt$id_case %in% ttree$id_case)][, -"join"]
  } else {
    incursions <- NULL # otherwise don't return incursions because you're not assigning them
  }
  
  # This is actually the slow part so limiting # of possibilities speeds things up a lot
  ttree[, c("dist_prob", "t_prob") := 
        .(dist_fun(dist_diff, pars = params, owned), 
          si_fun(t_diff, pars = params))]
  
  # source probability
  ttree[, source_prob := dist_prob * t_prob][, prob_scale := source_prob/sum(source_prob), by = id_case]
  ttree <- ttree[, selected := assign_progen(prob_scale), by = id_case][selected == 1] # filter to the chosen progenitor
  ttree[, prob_ll := log(source_prob)]
  
  ttree <- rbind(incursions, ttree, fill = TRUE)
  
  if(use_known_source) {  ttree <- rbind(ttree, known_tree, fill = TRUE) }
  ttree[, incursion := is.na(id_progen)]
  
  return(ttree)
  
}

# Assign progenitor ----
assign_progen <- function(prob_scaled) {
  
  cum_scaled <- cumsum(prob_scaled)
  rand_var <- runif(1, 0, 1)
  ind <- cum_scaled[cum_scaled > rand_var][1]
  ifelse(cum_scaled == ind, 1L, 0L)
  
}

# Helper functions ----------

# Add uncertainty  -----
add_uncertainty <- function(uncertainty){
  
  days_offset <- unlist(lapply(uncertainty, sample, size = 1)) * sample(c(-1, 1), length(uncertainty), replace = TRUE)
  
}

# Distribution functions (probably way to simplify this) -----
si_gamma1 <- function(x, pars) {
  dgamma(x, shape = pars$SI_shape, scale = pars$SI_scale)
}

si_gamma2 <- function(x, pars) {
  dgamma(x, shape = pars$SI2_shape, scale = pars$SI2_scale)
}

si_lnorm1 <- function(x, pars) {
  dlnorm(x, meanlog = pars$SI_ml, sdlog = pars$SI_sdlog)
}

si_lnorm2 <- function(x, pars) {
  dlnorm(x, meanlog = pars$SI2_ml, sdlog = pars$SI2_sdlog)
}

dist_gamma_mixed <- function(x, pars, owned) {
  
  out <- x
  out[owned] <- dgamma(x[owned], shape= pars$DK_shape, scale = pars$DK_scale)
  out[!owned] <- dgamma(x[!owned], shape= pars$DK2_shape, scale = pars$DK2_scale)
  out[x <= 100] <- pgamma(100, shape = pars$DK_shape, scale = pars$DK_scale)/100
  out
  
}

dist_lnorm_mixed <- function(x, pars, owned) {
  out <- x
  out[owned] <- dlnorm(x[owned], meanlog =  pars$DK_meanlog, sdlog = pars$DK_sdlog)
  out[!owned] <- dlnorm(x[!owned], meanlog =  pars$DK2_meanlog, sdlog = pars$DK2_sdlog)
  out[x <= 100] <- plnorm(100, meanlog =  pars$DK_meanlog, sdlog = pars$DK_sdlog)/100
  out
}

dist_gamma1 <- function(x, pars, owned) {
  
  out <- dgamma(x, shape= pars$DK_shape, scale = pars$DK_scale)
  out[x <= 100] <- pgamma(100, shape = pars$DK_shape, scale = pars$DK_scale)/100
  out
  
}

dist_lnorm1 <- function(x, pars, owned) {
  
  out <- dlnorm(x, meanlog =  pars$DK_meanlog, sdlog = pars$DK_sdlog)
  out[x <= 100] <- plnorm(100, meanlog =  pars$DK_meanlog, sdlog = pars$DK_sdlog)/100
  out
}

dist_gamma2 <- function(x, pars, owned) {
  
  out <- dgamma(x, shape= pars$DK2_shape, scale = pars$DK2_scale)
  out[x <= 100] <- pgamma(100, shape = pars$DK2_shape, scale = pars$DK2_scale)/100
  out
  
}

dist_lnorm2 <- function(x, pars, owned) {
  
  out <- dlnorm(x, meanlog =  pars$DK2_meanlog, sdlog = pars$DK2_sdlog) # use convoled for all
  out[x <= 100] <- plnorm(100, meanlog =  pars$DK2_meanlog, sdlog = pars$DK2_sdlog)/100
  out
  
}

# Get cutoffs -----
si_cutoff_lnorm1 <- function(cutoff, pars) {
  
  cut <- qlnorm(cutoff, meanlog = pars$SI_ml, sdlog = pars$SI_sdlog)

}

si_cutoff_gamma1 <- function(cutoff, pars) {
  
  cut <- qgamma(cutoff, shape = pars$SI_shape, scale = pars$SI_scale)
  
}

si_cutoff_lnorm2 <- function(cutoff, pars) {
  
  cut <- qlnorm(cutoff, meanlog = pars$SI2_ml, sdlog = pars$SI2_sdlog)
  
}

si_cutoff_gamma2 <- function(cutoff, pars) {
  
  cut <- qgamma(cutoff, shape = pars$SI2_shape, scale = pars$SI2_scale)
  
}

dist_cutoff_lnorm1 <- function(cutoff, pars) {
  
  cut <- qlnorm(cutoff, meanlog = pars$DK_meanlog, sdlog = pars$DK_sdlog)
  
}

dist_cutoff_gamma1 <- function(cutoff, pars) {
  
  cut <- qgamma(cutoff,  shape= pars$DK_shape, scale = pars$DK_scale)
  
}


dist_cutoff_lnorm2 <- function(cutoff, pars) {
  
  cut <- qlnorm(cutoff, meanlog = pars$DK2_meanlog, sdlog = pars$DK2_sdlog)
  
}

dist_cutoff_gamma2 <- function(cutoff, pars) {
  
  cut <- qgamma(cutoff,  shape= pars$DK2_shape, scale = pars$DK2_scale)
  
}

# Should these differ by whether owned or not?
dist_cutoff_lnorm_mixed <- dist_cutoff_lnorm2
dist_cutoff_gamma_mixed <- dist_cutoff_gamma2
  
# Helper for functions to export to foreach loop  -----
list_funs <- function(filename) {
  temp.env <- new.env()
  sys.source(filename, envir = temp.env)
  functions <- ls(envir = temp.env)
  return(functions)
}

# Wrapper function for running N times reproducibly & in parallel -----
boot_trees <- function(id_case,
                       id_biter, 
                       y_coord,
                       x_coord,
                       owned, 
                       date_symptoms, # needs to be in a date class
                       days_uncertain,
                       use_known_source = FALSE, 
                       prune = TRUE,
                       si_pdist = c("gamma", "lnorm"), 
                       dist_pdist = c("gamma", "lnorm"), 
                       params = list(), 
                       N = 1, 
                       verbose = TRUE,
                       seed = 1245, 
                       known_tree,
                       exp_funs = list_funs("R/fast_trees.R")) {
  
  if(any(is.na(x_coord) | is.na(y_coord) | is.na(date_symptoms))) {
    stop("Missing data in times or locations!")
  }
  
  # check that date_symptoms is a date class
  if(class(date_symptoms) != "Date") {
    stop("date_symptoms is not of class Date!")
  }
  
  foreach(i = seq_len(N), .combine = 'rbind', .options.RNG = seed, 
          .export = exp_funs, 
          .packages = c("data.table")) %dorng% {
    ttree <- 
      faster_trees(id_case = id_case, id_biter = id_biter, y_coord = y_coord, 
                   x_coord = x_coord,
                   owned = owned, date_symptoms = date_symptoms, 
                   days_uncertain = days_uncertain, use_known_source = use_known_source,
                   prune = prune,
                   si_pdist = si_pdist, dist_pdist = dist_pdist, params = params, 
                   known_tree = known_tree) 
      ttree$sim <- i
      ttree
    
  }
  
}