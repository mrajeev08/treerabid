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
#' @param lineages a data table with two columns, id_case and lineage, designating a lineage
#'  assignment for each case, defaults to NULL which means trees wont be resolved to a phylogeny
#' @param exclude_progen boolean of length id_case or 1, if TRUE then case should be excluded as a potential
#'  progenitor (i.e. if including livestock cases or other species that are dead-end transmissions)
#' @param use_known_source whether to assign known progenitors from contact tracing data
#' @param prune whether to prune links (i.e. based on Cori et al./Mancy et al.)
#'  at a certain cutoff probability. This also results in assignment of incursions
#' @param si_fun the function to get the probability of a given temporal
#'  difference between cases in days (i.e. the serial interval)
#' @param dist_fun the function to get the probability of a given spatial
#'  difference between cases in meters (i.e. the dispersal kernel)
#' @param cutoff the probability level at which to prune links (i.e. cases can not be
#'  linked if the distance or time difference is greater than the this %ile of the distribution);
#'  can either be a single number or a named vector with a cutoff called "dist" and "time"
#' @param params list of parameters to pass to si_fun and dist_fun
#' @param min_time if uncertainty results in negative time difference between
#'  known case pairs, set it to this value (better way to deal with propagating
#'  uncertainty given known case pairs?)
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
                       lineages = NULL,
                       exclude_progen = FALSE,
                       use_known_source = FALSE,
                       known_tree = NULL,
                       prune = TRUE,
                       si_fun,
                       dist_fun,
                       cutoff = 0.95,
                       params,
                       min_time = 1e-6) {

  if(length(cutoff) > 1) {
    if(length(cutoff) > 2 | !(c("time", "dist") %in% names(cutoff))) {
      stop("Should pass two cutoff thresholds in a vector named 'dist' and 'time'")
    }
    time_cut <- cutoff["time"]
    dist_cut <- cutoff["dist"]
  } else {
    time_cut <- cutoff
    dist_cut <- cutoff
  }

  if(any(cutoff > 1) | any(cutoff <= 0)) {
    stop("Cutoff value should be a single probability between (0, 1]")
  }

  # build line list with uncertainty
  t <- add_uncertainty(days_uncertain, date_symptoms, id_biter,
                       id_case, use_known_source)
  case_dt <- data.table(id_case, id_biter, x_coord, y_coord, owned, t, join = t)
  case_dt$type <- "reconstructed" # for tracking whether using known sources or not

  setkey(case_dt, join) # this speeds up join

  # Get the si cutoff if pruning (either length 1 or length nrow(case_dt))
  progen_dt <- case_dt[, .(t_progen = t,
                           x_coord_progen = x_coord,
                           y_coord_progen = y_coord,
                           id_progen = id_case,
                           join_on = join)]

  # exclude progenitors that are invalid (i.e. if any are livestock species, etc)
  progen_dt <- progen_dt[!exclude_progen, ]

  setkey(progen_dt, join_on)

  if(use_known_source) {

    if(is.null(known_tree)) {
      stop("Need to pass a reconstructed tree of known progenitors to `known_tree`,
            see `build_known_tree` function.")
    }

    k_tree <- copy(known_tree) # so working on copy specific to sim

    # get the uncertainty for the simulation (this will be the first match)
    k_tree$t <- case_dt$t[match(k_tree$id_case, case_dt$id_case)]
    k_tree$t_progen <- case_dt$t[match(k_tree$id_progen, case_dt$id_case)]

    # Get their diff for time
    k_tree[, t_diff := as.numeric(t - t_progen)]

    # fix any that have negative values due to uncertainty (better way?)
    k_tree$t_diff[k_tree$t_diff <= 0] <- min_time

    # Deal with multiple id's here (selecting ones that have multiple potential progenitors)
    k_tree <- select_progenitor(tree = k_tree, k_tree = NULL, lineages = NULL,
                                incursions = NULL,
                                si_fun = si_fun, dist_fun = dist_fun,
                                params = params, known = TRUE)

    # Filter out of progenitor assigment (but not out of the candidate progens!)
    case_dt <- case_dt[!(id_case %in% k_tree$id_case)]

  } else {
    k_tree <- NULL
  }

  # do a inner join to get possible progenitors
  if(prune) {
    # get the si cutoff
    si_cutoff <- si_fun(ttree = case_dt, cutoff = time_cut, params = params)
    progen_dt[, max := join_on + si_cutoff]

    # pruning at first stage based on time cutoffs to speed up
    ttree <- case_dt[progen_dt,
                     on = .(join > join_on, join <= max),
                     allow.cartesian = TRUE, nomatch = NULL][, -c("join", "join.1")]
  } else {
    ttree <- case_dt[progen_dt,
                     on = .(join > join_on),
                     allow.cartesian = TRUE, nomatch = NULL][, -"join"]
  }

  if(use_known_source) {
    # a candidate progenitor cannot be a known secondary case from contact tracing
    ttree <- ttree[!(k_tree[, .(id_progen = id_case, id_case = id_progen)]),
                   on = c("id_case", "id_progen")]
  }

  # For each case + progenitor combination, get the time and distance lag
  ttree[, c("dist_diff", "t_diff") := .(sqrt((x_coord - x_coord_progen)^2 + (y_coord - y_coord_progen)^2),
                                        as.numeric(t - t_progen))]

  # prune again for distance if true
  if(prune) {
    # Get the dist cutoff if pruning (either length 1 or length nrow(case_dt))
    dist_cutoff <- dist_fun(ttree = ttree, cutoff = dist_cut, params = params)
    ttree <- ttree[dist_diff <= dist_cutoff]
  }

  # and assign incursions as cases where no progenitor was identified
  # based on the distance OR time cutoffs
  # and also even if not pruning, the first case in the dataset will
  # not be assigned a link so need to include it in the tree here
  incursions <- case_dt[!(case_dt$id_case %in% ttree$id_case)][, -"join"]

  # This is actually the slow part so limiting # of possibilities speeds things up a lot
  # Also joins up with known tree and incursions
  ttree <- select_progenitor(tree = ttree, k_tree = k_tree, lineages = lineages,
                             incursions = incursions,
                             si_fun = si_fun, dist_fun = dist_fun,
                             params = params, known = FALSE)

  return(ttree)

}

#' Reconstruct known tree (i.e. from tracing data)
#'
#' @inheritParams build_tree
#'
#' @return a data.table with the known contact tracing tree
#' @keywords internal
#'
build_known_tree <- function(id_case,
                             id_biter,
                             y_coord,
                             x_coord,
                             owned,
                             date_symptoms) {


  # build line list with uncertainty
  case_dt <- data.table(id_case, id_biter, x_coord, y_coord, owned, t = date_symptoms)

  # Get the si cutoff if pruning (either length 1 or length nrow(case_dt))
  progen_dt <- case_dt[, .(t_progen = t,
                           x_coord_progen = x_coord,
                           y_coord_progen = y_coord,
                           id_progen = id_case,
                           join_ct = id_case)]
  setkey(progen_dt, join_ct)

  # Filter ones with non-zero biter ids
  # making sure the ids exist in the progenitor list
  known_biters <- case_dt[id_biter != 0 & id_biter %in% progen_dt$id_progen]
  known_biters[, join_ct := id_biter]
  setkey(known_biters, join_ct)

  # Inner join with the biters
  known_tree <- progen_dt[known_biters, on = "join_ct"][, -"join_ct"]

  # Get their distance difference and type
  known_tree[, c("dist_diff",
                 "type") := .(sqrt((x_coord - x_coord_progen)^2 + (y_coord - y_coord_progen)^2),
                              "traced")]

  return(known_tree)
}

# Helper functions ----------

# Select progenitor
#' Wrapper to select single progenitor in each case
#'
#' @param tree the data.table with possible case pairs to select from
#' @param k_tree the data.table with the known tree (when use_known_source = TRUE,
#'  NULL otherwise)
#' @param incursion the data.table with incursions (when prune =TRUE, NULL otherwise)
#' @inheritParams build_tree
#'
#' @return a data.table filtered to the selected case-progenitor pair
#' @keywords internal
#'
select_progenitor <- function(tree, lineages, k_tree, incursions,
                              si_fun, dist_fun, params, known = FALSE) {

  # probabilities
  si_fun(ttree = tree, params = params)
  dist_fun(ttree = tree, params = params)
  tree[, source_prob := dist_prob * t_prob][, prob_scale := source_prob/sum(source_prob), by = id_case]
  ttree <- tree[, selected := assign_progen(prob_scale), by = id_case][selected == 1]

  if(!known) {
    # Bind to incursions & known tree (can be NULL)
    ttree <- rbindlist(list(incursions, k_tree, ttree), fill = TRUE)

    # Check which ones have mismatched lineages
    ttree <- ttree[lineages, on = "id_case"]
    tree <- tree[lineages, on = "id_case"]

    known_progens <- k_tree$id_case

    out <- find_lins_to_fix(ttree, known_progens, selector = "prob_scale")

    lins_to_fix <- out$lins_to_fix
    membership_dt <- out$membership_dt

    # set links to fix to NA
    ttree[id_case %in% lins_to_fix]$id_progen <- NA
    nfixes <- length(lins_to_fix)

    if(nfixes > 0) {

      rfixes <- sample(nfixes, nfixes) # fix in random order

      for(i in rfixes) {

        # Join up links with updated membership_dt
        tree <- membership_dt[tree, on = "id_case"]
        setnames(membership_dt, c("membership", "id_case", "lineage_chain"),
                 c("membership_progen", "id_progen", "lineage_progen_chain"))
        tree <- membership_dt[tree, on = "id_progen"]
        tree[, membership_progen := ifelse(is.na(membership_progen), 0,
                                           membership_progen)]
        tree[, lineage_progen_chain := ifelse(is.na(lineage_progen_chain), 0,
                                        lineage_progen_chain)]

        # For those that have no incursions and are also the minimum case replace
        # with the next most likely progen, that is not already in the current chains
        candidate_links <- tree[id_case %in% lins_to_fix[i] & membership != membership_progen]

        # Also filter to those that in chain with same (or totally unsampled lineages)
        candidate_links <- candidate_links[lineage_chain * lineage_progen_chain == 0 | lineage_chain * lineage_progen_chain == lineage_chain^2]

        # Rescale probabilities & select
        fixed_links <- candidate_links[, prob_scale := source_prob/sum(source_prob),
                                       by = id_case][, selected := assign_progen(prob_scale),
                                                     by = id_case][selected == 1]

        # if none then set to NA (prob & links)
        incs <- ifelse(nrow(fixed_links) == 0, lins_to_fix[i], 0)
        set_incs <- ttree[id_case %in% incs][, c("id_biter", "id_case", "x_coord",
                                                 "y_coord", "owned", "t", "type",
                                                 "lineage")]
        set_incs[, id_progen := NA]

        # bind them together
        ttree <-
          rbindlist(list(ttree[!(id_case %in% lins_to_fix[i])],
                         fixed_links, set_incs),
                    fill = TRUE)

        # clean ttree & links_all
        ttree[, c("membership", "membership_progen", "lineage_chain", "lineage_progen_chain") := NULL]
        tree[, c("membership", "membership_progen", "lineage_chain", "lineage_progen_chain") := NULL]

        # update membership_dt
        membership_dt <- get_membership(ttree)

      }
    }

    ttree <- membership_dt[ttree, on = "id_case"]

  }

  ttree[, prob_ll := log(source_prob)]
  ttree[, incursion := is.na(id_progen)]

  return(ttree)
}

#' Assign progenitor
#'
#' @keywords internal
assign_progen <- function(prob_scaled) {

  cumul_scaled <- cumsum(prob_scaled)
  rand_var <- runif(1, 0, 1)
  ind <- which(cumul_scaled > rand_var)[1]
  selected <- rep(0L, length(cumul_scaled))
  selected[ind] <- 1L

  return(selected)
}

#' Add uncertainty
#'
#' @keywords internal
add_uncertainty <- function(uncertainty, date_symptoms, id_biter,
                            id_case, use_known_source, buffer = 7,
                            max_tries = 100){

  if(any(uncertainty > 0)) {

    days_offset <- unlist(lapply(uncertainty, sample, size = 1))
    sign <- sample(c(-1, 1), length(uncertainty), replace = TRUE)
    date_uncertain <- date_symptoms +  days_offset * sign
    niter <- 0

    if(use_known_source) {
      date_min <- date_uncertain[match(id_biter, id_case)] + buffer # this will be relative to the first one
      date_min[is.na(date_min)] <- date_uncertain[is.na(date_min)]

      while(any(date_uncertain < date_min) & niter <= max_tries) {
        niter <- niter + 1
        # constrain known biters
        date_uncertain[date_uncertain < date_min] <- date_min[date_uncertain < date_min]

        # check again
        date_min <- date_uncertain[match(id_biter, id_case)] + buffer
        date_min[is.na(date_min)] <- date_uncertain[is.na(date_min)]
      }

      if(any(date_uncertain < date_min)) {
        warning("Significant date uncertainty means that some case dates may
               not line up with known sequence of events from contact tracing!")
      }
    }

    return(date_uncertain)

  } else {
    return(date_symptoms)
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
#'  certain types of cluster configs
#' @param ncores the number of cores to parallelize over
#'  (defaults to parallel::detectCores() - 1)
#' @param chunk whether to chunk sims or not
#'  (i.e. to increase efficiency of parallelization, may want to run 100/core vs.
#'   1/core, as sims can be mem intensive,
#'   best to try out on your architecture to see which one faster)
#'
#' @return a data.table with bootstrapped trees
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#' @importFrom parallel detectCores
#' @export
#'
boot_trees <- function(id_case,
                       id_biter,
                       y_coord,
                       x_coord,
                       owned,
                       date_symptoms, # needs to be in a date class
                       days_uncertain,
                       lineages = NULL,
                       exclude_progen = FALSE,
                       use_known_source = FALSE,
                       prune = TRUE,
                       si_fun,
                       dist_fun,
                       cutoff = 0.95,
                       params,
                       min_time = 1e-6,
                       N = 1,
                       seed = 1245,
                       exp_funs = NULL,
                       exp_pkgs = c("data.table", "treerabid", "igraph"),
                       ncores = detectCores() - 1,
                       chunk = TRUE) {

  if(any(is.na(x_coord) | is.na(y_coord) | is.na(date_symptoms))) {
    stop("Missing data in times or locations!")
  }

  # check that date_symptoms is a date class
  if(class(date_symptoms) != "Date") {
    stop("date_symptoms is not of class Date!")
  }

  if(use_known_source){
    known_tree <- build_known_tree(id_case,
                                   id_biter,
                                   y_coord,
                                   x_coord,
                                   owned,
                                   date_symptoms)
  } else {
    known_tree <- NULL
  }

  # Check for duplicated id's here
  # (i.e. depending on whether you are using a known source or not)
  msg <- "id_case has duplicated values,
          but you are not using a known set of possible sources for these cases,
          there should only be one record per case id!"
  if(!use_known_source) {
    dups <- tabulate(id_case)
    if(any(dups > 1)) {
      stop(msg)
    }
  } else {
    # any biter ids that are < 1 (i.e. 0, neg numbers, and NAs)
    dups <- tabulate(id_case[id_biter < 1])
    if(any(dups > 1)) {
      stop(msg)
    }

  }

  if(N <= ncores | !chunk) {
    chnks <- N
    grps <- sims <- seq(1, N)
  } else {
    chnks <- floor(N/ncores)
    sims <- seq(1, N)
    grps <- rep(seq(1, chnks), N)[1:length(sims)]
  }

  foreach(i = seq_len(chnks),
          .combine = 'rbind', .options.RNG = seed,
          .export = exp_funs,
          .packages = exp_pkgs) %dorng% {

          nsims <- sims[grps == i]

          rbindlist(
            lapply(
              nsims,
              function (x) {
                ttree <-
                  build_tree(id_case = id_case, id_biter = id_biter, y_coord = y_coord,
                             x_coord = x_coord,
                             owned = owned, date_symptoms = date_symptoms,
                             days_uncertain = days_uncertain,
                             exclude_progen = exclude_progen,
                             use_known_source = use_known_source,
                             known_tree = known_tree,
                             lineages = lineages,
                             prune = prune,
                             si_fun,
                             dist_fun,
                             cutoff = cutoff,
                             params = params,
                             min_time = min_time)
                ttree$sim <- x
                ttree
              }
            )
          )
    }

}
