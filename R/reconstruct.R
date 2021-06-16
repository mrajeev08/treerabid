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
#' @param lineages a vector of integer lineage ids based on a phylogeny (0 = unsequenced cases);
#'  defaults to NULL;
#' @param max_iter integer, when using lineage data, the number of tries to fix links in
#'  trees to be consistent with lineage data
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
#'  linked if the distance or time difference is greater than the this %ile of the distribution)
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
                       max_iter = 100,
                       exclude_progen = FALSE,
                       use_known_source = FALSE,
                       known_tree = NULL,
                       prune = TRUE,
                       si_fun,
                       dist_fun,
                       cutoff = 0.95,
                       params,
                       min_time = 1e-6) {

  if(cutoff >= 1 | cutoff <= 0 | length(cutoff) > 1) {
    stop("Cutoff value should be a single probability between (0, 1)")
  }

  if(!is.null(lineages)) {
    lineages <- data.table(id_case, lineage = lineages)
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
    k_tree <- select_progenitor(tree = k_tree,
                                si_fun = si_fun,
                                dist_fun = dist_fun,
                                params = params)

    # Filter out of progenitor assigment (but not out of the candidate progens!)
    case_dt <- case_dt[!(id_case %in% k_tree$id_case)]

  }

  # do a inner join to get possible progenitors
  if(prune) {
    # get the si cutoff
    si_cutoff <- si_fun(ttree = case_dt, cutoff = cutoff, params = params)
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
    dist_cutoff <- dist_fun(ttree = ttree, cutoff = cutoff, params = params)
    ttree <- ttree[dist_diff <= dist_cutoff]

    # and assign incursions as cases where no progenitor was identified
    # based on the distance OR time cutoffs
    incursions <- case_dt[!(case_dt$id_case %in% ttree$id_case)][, -"join"]
  } else {
    incursions <- NULL # otherwise don't return incursions because you're not assigning them
  }

  # This is actually the slow part so limiting # of possibilities speeds things up a lot
  ttree <- select_progenitor(tree = ttree, incursions = incursions,
                             k_tree = k_tree,
                             si_fun = si_fun, dist_fun = dist_fun,
                             params = params, lineages = lineages, max_iter = max_iter)


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

# Select progenitor incorporating phylogenetic data
#' Wrapper to select single progenitor in each case
#'
#' @param tree the data.table with possible case pairs to select from
#' @inheritParams build_tree
#'
#' @return a data.table filtered to the selected case-progenitor pair
#' @keywords internal
#'
select_progenitor <- function(tree, incursions = NULL, k_tree = NULL,
                              si_fun, dist_fun,
                              params, lineages = NULL, max_iter = 100) {

  # probabilities
  si_fun(ttree = tree, params = params)
  dist_fun(ttree = tree, params = params)
  tree[, source_prob := dist_prob * t_prob][, prob_scale := source_prob/sum(source_prob), by = id_case]
  tree[, prob_ll := log(source_prob)]

  # Select progenitors
  tree[, selected := assign_progen(prob_scale), by = id_case]

  # Bind with incursions and known tree (setting these to selected)
  if(!is.null(incursions)) {
    incursions$selected <- 1
    tree <- rbind(incursions, tree, fill = TRUE)
  }

  if(!is.null(k_tree)) {
    k_tree$selected <- 1
    tree <- rbind(tree, k_tree, fill = TRUE)
  }

  tree[, incursion := is.na(id_progen)]

  if(is.null(lineages)) {
    tr_test <- tree[selected == 1]

  } else {

    tree <- lineages[tree, on = "id_case"]
    tr_test <- tree[selected == 1]

    to_fix <- check_lineages(tr_test, lineages)
    tries <- 1

    while(!is.null(to_fix) & tries < max_iter) {

      # Get the unique broken links
      to_fix <- to_fix[, .(check = .N), by = c("id_case", "i.id_case")]

      # set to NA to remake membership graph (and set other cols to NA as well!)
      tr_test[id_case %in% to_fix$i.id_case]$id_progen <- NA

      # filter to candidates to reselect to
      tr_reselect <- tree[id_case %in% to_fix$i.id_case]

      # Filter out the ones we broke with an antijoin
      tr_reselect <- tr_reselect[!(to_fix[, .(id_progen = id_case, id_case = i.id_case)]),
                                 on = c("id_case", "id_progen")]

      # Remake membership
      mbr_new <- get_membership(tr_test, lineages)

      # Filter to ones that are not in the same chain or in a chain that has
      # a different sampled lineage id
      tr_reselect <- tr_reselect[mbr_new[,
                                         .(id_case,
                                           mbr_case = membership,
                                           lin_case = lineage)], on = "id_case"]
      tr_reselect <- tr_reselect[mbr_new[,
                                         .(id_progen = id_case,
                                           mbr_progen = membership,
                                           lin_progen = lineage)], on = "id_progen"]
      tr_reselect <- tr_reselect[mbr_progen != mbr_case &
                                 lin_case * lin_progen %in% c(lin_case^2, 0)]

      # reselect progenitors
      tr_reselect[, prob_scale := source_prob/sum(source_prob), by = id_case]
      tr_reselect[, selected := assign_progen(prob_scale), by = id_case]

      tr_test <- rbind(tr_test[!(id_case %in% tr_reselect$id_case)],
                       tr_reselect)

      # Repeat
      to_fix <- check_lineages(tr_test, lineages)
      tries <- tries + 1
    }

    if(!is.null(to_fix)) {
      warning("Some inconsistencies in lineage assignments, up max_iter to
              see if these can be resolved (although note it will increase
              computational time!")
    }
  }

  return(tr_test)
}

#' Title
#'
#' @param tr_test
#'
#' @return
#' @export
#'
check_lineages <- function(tr_test, lineages) {

  # build directed graph
  gr <- graph_from_data_frame(d = tr_test[, c("id_progen",
                                              "id_case")][!is.na((id_progen))],
                              vertices = lineages,
                              directed = TRUE)

  # Get the chain membership
  V(gr)$membership <- components(gr)$membership
  membership <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")))
  tr_test <- tr_test[membership, on = "id_case"]

  # Filter to chains that have multiple lineages per chain
  multilins <- tr_test[lineage != 0][, .(check = length(unique(lineage))),
                                     by = "membership"][check > 1]
  reassign <- tr_test[membership %in% multilins$membership]

  # Filter to those sampled
  reassign <- reassign[lineage != 0]

  # join reassign with itself
  reassign <- reassign[reassign, on = .(membership == membership), allow.cartesian = TRUE]

  # filter out same case ids and same lineages
  reassign <- reassign[id_case != i.id_case & lineage != i.lineage]

  if(nrow(reassign) > 0) {

    reassign$row_id <- 1:nrow(reassign)
    pths <- get_edge_dt(gr, lins = reassign)
    pths <- reassign[pths, on = "row_id"]
    browser()
    pths[, .(freq = .N), by = c("from", "to")]
    # Select the edge to break: get a random variate to select by
    pths$selector <- runif(nrow(pths))

    # filter out any pths i.id_case already has a known progenitor!
    browser()

    # Find the minimum selector
    to_fix <- pths[pths[, .I[which.min(selector)],
                        by = c("row_id")]$V1]

  } else {
    to_fix <- NULL
  }

  return(to_fix)
}

#' Title
#'
#' @param tr_test
#'
#' @return
#' @export
#'
get_membership <- function(tr_test, lineages) {

  # build directed graph
  gr <- graph_from_data_frame(d = tr_test[, c("id_case",
                                              "id_progen")][!is.na((id_progen))],
                              vertices = lineages,
                              directed = TRUE)

  # Get the chain membership
  V(gr)$membership <- components(gr)$membership
  membership <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")),
                           lineage = as.numeric(vertex_attr(gr, "lineage")))

  # will only be one non zero lineage per chain as this is after mismatches are broken
  membership[, lineage := sum(lineage), by = "membership"]

  return(membership)
}

#' Title
#'
#' @param gr
#' @param lins
#'
#' @return
#' @export
#'
get_edge_dt <- function(gr, lins) {

  sps <- shortest_paths(gr, from = as.character(lins$id_case),
                        to = as.character(lins$i.id_case),
                        mode = "all",
                        output = "epath")$epath

  rbindlist(
    lapply(seq_len(length(sps)),
           function(x) {
             dt <- as_data_frame(subgraph.edges(gr, sps[[x]],
                                                    delete.vertices = FALSE))
             if(nrow(dt) > 0) {
               dt$row_id <- x
             }
             return(dt)
           }
           ), fill = TRUE
    )

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
#'  certain types of cluster configs,
#' @param ncores number of cores in cluster, for how to split up the parallelization
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
                       lineages = NULL,
                       max_iter = 100,
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
                       exp_pkgs = c("data.table", "treerabid"),
                       ncores = parallel::detectCores() - 1) {

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

  chnks <- seq(1, N, by = floor(N/ncores))
  from <- c(1, chnks[2:ncores] + 1)
  to <- chnks[2:(ncores + 1)]

  foreach(i = seq_len(length(from)),
          .combine = 'rbind', .options.RNG = seed,
          .export = exp_funs,
          .packages = exp_pkgs) %dorng% {

          nsims <- seq(from[i], to[i])

          rbindlist(
            lapply(
              nsims,
              function (x) {
                ttree <-
                  build_tree(id_case = id_case, id_biter = id_biter, y_coord = y_coord,
                             x_coord = x_coord,
                             owned = owned, date_symptoms = date_symptoms,
                             days_uncertain = days_uncertain,
                             lineages = lineages,
                             max_iter = max_iter,
                             exclude_progen = exclude_progen,
                             use_known_source = use_known_source,
                             known_tree = known_tree,
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
