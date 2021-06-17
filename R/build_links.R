#' Summarize links between cases across bootstrapped transmission trees
#'
#' @param ttrees bootstrapped trees from `boot_trees`
#' @param N number of bootstrapped trees
#'
#' @return a data.table with all linked cases across trees by their frequency and
#'  probability
#' @export
#'
build_all_links <- function(ttrees, N) {

  links_all <- ttrees[, .(links = .N,
<<<<<<< HEAD
                          t_diff_median_days = median(t_diff),
                          dist_diff_meters = median(dist_diff)),
=======
                          t_diff_median = median(t_diff),
                          dist_diff = median(dist_diff)),
>>>>>>> 57c51cf975722fe7b76e84441f8d926594f58fa4
                      by = c("id_case", "id_progen")][, prob := links/N]

  return(links_all)

}

#' Build consensus links between cases (i.e. most often selected progenitors for a given case)
#'
#' @param links_all output from `build_all_links`
#' @param case_dates data.table with two columns id_case (id of case) and symptoms_started
#'  (date symptoms started of case without uncertainty)
#' @param lineages a data table with two columns, id_case and lineage, designating a lineage
#'  assignment for each case, defaults to NULL which means trees wont be resolved to a phylogeny
#' @param fix_loops If there are loops in the transmission tree (i.e. indicating uncertainty in
#'  who-infected-whom), the loops can either be broken by reassigning the progenitor
#'  of the case with the earliest case date in the loop (fix_loops = "by_date", the default).
#'  or by breaking and replacing
#'  links in the loop randomly (fix_loops = "random").
#' @param known_progens a numeric vector of case ids for which progenitors are known
#' @param max_tries number of times to iterate through and break loops/inconsistent lineage assignments
#'
#' @return a data.table with the case pairs most often linked in transmission tree and their
#'  probability (i.e. prop of times a given case was selected as a progenitor)
#' @export
#'
build_consensus_links <- function(links_all,
                                  case_dates = NULL,
                                  lineages = NULL,
                                  known_progens = NULL,
                                  fix_loops = c("by_date", "random"),
                                  max_tries = 100) {

  fix_loops <- match.arg(fix_loops)
  # Get the consensus links
  links_consensus <- links_all[links_all[, .I[which.max(links)], by = c("id_case")]$V1] # returns first max

  if(is.null(lineages)) {
    lineages <- data.table(id_case = links_consensus$id_case, lineage = 0)
  }

  links_consensus <- links_consensus[lineages, on = "id_case"]

  # spit out the links, the update consensus links, and the updated membership
  list2env(find_links_to_fix(links_consensus, fix_loops, case_dates,
                             known_progens), envir = environment())
  niter <- 0

  # option to reassign to next likely and check loops
  while(length(to_fix) > 0 & niter < max_tries) {

    niter <- niter + 1

    # Join up links with updated membership
    links_all <- membership[links_all, on = "id_case"]
    setnames(membership, c("membership", "id_case", "lineage"),
             c("membership_progen", "id_progen", "lineage_progen"))
    links_all <- membership[links_all, on = "id_progen"]
    links_all[, membership_progen := ifelse(is.na(membership_progen), 0,
                                            membership_progen)]
    links_all[, lineage_progen := ifelse(is.na(lineage_progen), 0,
                                         lineage_progen)]

    # For those that have no incursions and are also the minimum case replace
    # with the next most likely progen, that is not already in the current chains
    candidate_links <- links_all[id_case %in% to_fix & membership != membership_progen]

    # Also filter to those that in chain with same (or totally unsampled lineages)
    candidate_links <- candidate_links[lineage * lineage_progen %in% c(lineage^2, 0)]

    fixed_links <- candidate_links[candidate_links[, .I[which.max(links)], by = c("id_case")]$V1]

    # if none then set to NA (prob & links)
    incs <- to_fix[!(to_fix %in% fixed_links$id_case)]
    set_incs <- links_consensus[id_case %in% incs]
    set_incs[, c("id_progen", "links", "prob") := .(NA, NA, NA)]

    # bind them together
    links_consensus<-
      rbindlist(list(links_consensus[!(id_case %in% to_fix)],
                     fixed_links, set_incs),
                fill = TRUE)

    # clean links_consensus & links_all
    links_consensus[, c("membership", "membership_progen", "lineage", "lineage_progen") := NULL]
    links_all[, c("membership", "membership_progen", "lineage", "lineage_progen") := NULL]

    # rejoin with lineages
    # to deal with issue that only known progenitors were possible for reassignment)
    links_consensus <- links_consensus[lineages, on = "id_case"]

    # update and recheck for to_fix
    list2env(find_links_to_fix(links_consensus, fix_loops, case_dates,
                               known_progens), envir = environment())

    if(niter == max_tries) {
      # update final membership & check loops if still some unresolved
      gr <- get_gr(links_consensus)
      loops <- check_loops(links_consensus)

      # Get the chain membership
      V(gr)$membership <- components(gr)$membership
      membership <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                               id_case = as.numeric(vertex_attr(gr, "name")))
      links_consensus <- links_consensus[membership, on = "id_case"]
      lins_to_fix <- check_lineages(links_consensus)
      loops <- check_loops(links_consensus)

      links_consensus[, membership := NULL]

    } else {
      lins_to_fix <- loops <- NULL
    }
  }


  # Test to make sure no more to_fix and warn if there are!
  if(length(loops) > 0) {

    warning(
      "There are still loops in the transmission tree, indicating that case date
      uncertainty may be too high to indentify introductions and differentiate chains.
      You can also try increasing the value of max_tries.")

  }

  if(nrow(lins_to_fix > 0)) {
    warning(
      "Couldn't completely resolve tree to phylogeny, try increasing the number
      of bootstrapped trees or max_tries.")
  }

  none_found <- sum(is.na(links_consensus$prob))

  if(none_found > 0) {
    warning(
      paste0(
        none_found,
        " cases were assigned no progenitor because all potential progenitors
        were filtered out by fixing loops or resolving to a phylogeny. If appropriate,
        you might consider increasing the number of bootstrapped trees to make sure
        this is not an artifact of sampling."))
  }

  # Clean up the data.table
  return(links_consensus)
}


#' Internal function for finding loops and mismatched lineages to fix
#'
#' @param links_consensus consensus links build within build_consensus_links function
#' @inheritParams build_consensus_links
#'
#' @return a numeric vector of case ids for which links should be reassigned
#'
#' @importFrom igraph V subgraph.edges E count_multiple girth components
#'  vertex_attr graph_from_data_frame
#' @keywords internal
#'
find_links_to_fix <- function(links_consensus, fix_loops, case_dates,
                              known_progens) {

  # build undirected & find the loops (which_multiple) & any cycles (girth)
  gr <- get_gr(links_consensus)

  loops <- check_loops(links_consensus)

  # Get the chain membership
  V(gr)$membership <- components(gr)$membership
  membership <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")))
  links_consensus <- links_consensus[membership, on = "id_case"]

  if(fix_loops == "by_date") {

    if(is.null(case_dates)) {
      stop("Need to pass case_dates data.table to fix loops by date.")
    }

    # join with case dates
    links_consensus <- links_consensus[case_dates, on = "id_case"]
    setnames(links_consensus, "symptoms_started", "selector")

  } else {

    # get a random variate to select by
    links_consensus$selector <- runif(nrow(links_consensus))

  }

  # Filter to the ones with loops
  # & filter out any that are known from contact tracing
  loops_to_fix <- links_consensus[id_case %in% loops & !(id_case %in% known_progens)]

  # Find the minimum selector (either by case date or randomly select link to break)
  loops_to_fix <- loops_to_fix[loops_to_fix[, .I[which.min(selector)],
                          by = c("membership")]$V1]$id_case

  # Resolve the phylogeny
  if(length(unique(links_consensus[lineage != 0]$lineage)) > 1) {
    lins_to_fix <- find_lineages(gr, links_consensus, known_progens)
  } else {
    lins_to_fix <- NULL
  }

  to_fix <- as.numeric(unique(c(loops_to_fix, lins_to_fix)))

  # Update links consensus and get new membership
  links_consensus[id_case %in% to_fix]$id_progen <- NA
  membership <- get_membership(links_consensus)

  return(list(membership = membership, to_fix = to_fix))
}

#' Helper function to get graph from links
#'
#' @param links consensus links or individual tree
#'
#' @return an igraph object
#'
get_gr <- function(links) {

  # build undirected & find the loops (which_multiple) & any cycles (girth)
  gr <- graph_from_data_frame(d = links[, c("id_case",
                                            "id_progen")][!is.na((id_progen))],
                              vertices = links[, c("id_case", "lineage")],
                              directed = FALSE)

  return(gr)

}


#' Helper function to check for loops
#'
#' @param links either the consensus links or a single tree (with
#'  cols id_case & id_progen)
#'
#' @return a vector of case ids which are part of a loop in the tree
#' @export
#'
check_loops <- function(links) {

  # build undirected & find the loops (which_multiple) & any cycles (girth)
  gr <- graph_from_data_frame(d = links[, c("id_case",
                                            "id_progen")][!is.na((id_progen))],
                              vertices = links[, "id_case"],
                              directed = FALSE)

  loops <- names(V(subgraph.edges(gr, E(gr)[count_multiple(gr) > 1])))
  loops <- as.numeric(c(loops, names(girth(gr)$circle)))

  return(loops)
}

#' Find the links to fix too resolve lineage discrepancies
#'
#' @param gr the graph of the links
#' @param links links with original chain membership and lineages
#' @inheritParams build_consensus_links
#' @return a numeric vector of case ids to fix
#'
find_lineages <- function(gr, links, known_progens) {

  # Filter to chains that have multiple lineages per chain
  multilins <- links[lineage != 0][, .(check = length(unique(lineage))),
                                     by = "membership"][check > 1]
  reassign <- links[membership %in% multilins$membership]

  # Filter to those sampled
  reassign <- reassign[lineage != 0][, c("id_case", "membership", "lineage")]

  # join reassign with itself
  reassign <- reassign[reassign, on = .(membership == membership), allow.cartesian = TRUE]

  # filter out same case ids and same lineages
  reassign <- reassign[id_case != i.id_case & lineage != i.lineage]

  # get unique pairs only
  reassign <- reassign[!duplicated(t(apply(reassign[,
                                                    c("id_case", "i.id_case")],
                                           1, sort))), ]

  if(nrow(reassign) > 0) {

    pths <- get_edge_dt(gr, lins = reassign)

    # when you get the paths you lose the directionality
    # so need to join back up with the actual links
    check <- rbind(pths[, .(id_progen = as.numeric(from),
                            id_case = as.numeric(to), row_id)],
                   pths[, .(id_progen = as.numeric(to),
                            id_case = as.numeric(from), row_id)])
    pths <- links[,
                  c("id_case", "id_progen")][check,
                                             on = c("id_case", "id_progen"),
                                             nomatch = NULL]

    # Filter out any edges that are known (from tracing)
    pths <- pths[!(id_case %in% known_progens)]

    # Filter to the most frequently selected ones for each row id
    pths[, freq := .N, by = c("id_case", "id_progen")]

    pths[, max_freq := max(freq), by = "row_id"]
    pths <- pths[freq == max_freq]

    # Select the edge to break: get a random variate to select by
    pths$selector <- runif(nrow(pths))

    # Find the minimum selector
    to_fix <- pths[pths[, .I[which.min(selector)],
                        by = c("row_id")]$V1]

  } else {
    to_fix <- NULL
  }

  return(unique(to_fix$id_case))
}

#' Get the number of lineages in each chain
#'
#' @param links
#'
#' @return
#' @export
#'
#' @examples
check_lineages <- function(links) {

  # Filter to chains that have multiple lineages per chain
  multilins <- links[lineage != 0][, .(check = length(unique(lineage))),
                                   by = "membership"][check > 1]
  return(multilins)
}

#' Internal function for getting edges between mismatched lineages
#'
#' @param gr
#' @param lins
#'
#' @return
#'
get_edge_dt <- function(gr, lins) {

  sps <-
    suppressWarnings(
      shortest_paths(gr,
                     from = as.character(lins$id_case),
                     to = as.character(lins$i.id_case),  mode = "all",
                     output = "epath")$epath
      )

  sps <- Filter(function(x) length(x) > 0, sps)

  if(length(sps) > 0) {
    out <-
      rbindlist(
        lapply(seq_len(length(sps)),
               function(x) {
                 dt <- as_data_frame(subgraph.edges(gr, sps[[x]],
                                                    delete.vertices = FALSE))
                 dt$row_id <- x
                 return(dt)
               }
        ))
  } else {
    out <- data.table(from = 0, to = 0, row_id = 0)[0]
  }

  return(out)
}

#' Internal function for checking membership of chains AFTER lineage mismatches
#' & loops broken
#'
#' @param links
#'
#' @return
#'
get_membership <- function(links) {

  # build directed graph
  gr <- get_gr(links)

  # Get the chain membership
  V(gr)$membership <- components(gr)$membership
  membership <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")),
                           lineage = as.numeric(vertex_attr(gr, "lineage")))

  # will only be one non zero lineage per chain as this is after mismatches are broken
  membership[, lineage := sum(unique(lineage)), by = "membership"]

  return(membership)
}

#' Build the consensus tree (i.e. the tree with the highest proportion of consensus links)
#'
#' @param links_consensus output from `build_consensus_links`
#' @param ttrees  bootstrapped trees from `boot_trees`
#'
#' @return a data.table with the consensus tree (with score for each link,
#'  1 if the consensus link, 0 if not.
#' @export
#'
build_consensus_tree <- function(links_consensus, ttrees) {

  sim_scores <- ttrees[links_consensus, on = c("id_case", "id_progen")][, score := 1][, .(score = sum(score)), by = "sim"]

  best_pairs <- paste0(links_consensus$id_case, "_", links_consensus$id_progen)
  tree_consensus <- ttrees[sim == sim_scores$sim[which.max(sim_scores$score)]]
  tree_consensus[, score := fifelse(paste0(id_case, "_", id_progen) %in% best_pairs, 1, 0)]

  return(tree_consensus)
}

