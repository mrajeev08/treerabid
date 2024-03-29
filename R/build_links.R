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
                          t_diff_median_days = median(t_diff),
                          dist_diff_meters = median(dist_diff)),
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
#' @param known_progens a numeric vector of case ids for which progenitors are known
#' @param max_cycles how many cycles to look for (if lots of date uncertainty
#'  you may want to increase this but it will likely slow things down a lot!)
#' @param link_all whether you want to try and force links between cases (that is
#'  even if a case is often identified as the beginning of a chain, you want to
#'  see if it can be assigned to another progenitor, i.e. when trying to link
#'  up all cases to a sampled set of lineages)
#'
#' @return a data.table with the case pairs most often linked in transmission tree and their
#'  probability (i.e. prop of times a given case was selected as a progenitor)
#' @export
#'
build_consensus_links <- function(links_all,
                                  case_dates,
                                  lineages = NULL,
                                  known_progens = NULL,
                                  max_cycles = 100,
                                  link_all = FALSE) {

  links_backup <- links_all[is.na(id_progen)]

  if(link_all == TRUE) {
    links_all <- links_all[!is.na(id_progen)]
  }

  # Get the consensus links
  links_consensus <- links_all[links_all[, .I[which.max(links)], by = c("id_case")]$V1] # returns first max
  links_to_add <- case_dates$id_case[!(case_dates$id_case %in% links_consensus$id_case)]
    
   if (length(links_to_add) > 0) {
      links_to_add <- links_backup[id_case %in% links_to_add]
      links_consensus <- rbind(links_consensus, links_to_add)
   }

  if(is.null(lineages)) {
    lineages <- data.table(id_case = links_consensus$id_case, lineage = 0)
  }

  links_consensus <- links_consensus[lineages, on = "id_case"]
  links_all <- links_all[lineages, on = "id_case"]

  # first fix the lineages
  list2env(find_lins_to_fix(links_consensus, known_progens,
                            selector = "prob"), envir = environment())

  # set links to fix to NA
  links_consensus[id_case %in% lins_to_fix]$id_progen <- NA
  nfixes <- length(lins_to_fix)

  if(nfixes > 0) {

    # Fix in order of least -> most likely
    for(i in seq_len(nfixes)) {

      # Join up links with updated membership_dt
      links_all <- membership_dt[links_all, on = "id_case"]
      setnames(membership_dt, c("membership", "id_case", "lineage_chain"),
               c("membership_progen", "id_progen", "lineage_progen_chain"))
      links_all <- membership_dt[links_all, on = "id_progen"]
      links_all[, membership_progen := ifelse(is.na(membership_progen), 0,
                                              membership_progen)]
      links_all[, lineage_progen_chain := ifelse(is.na(lineage_progen_chain), 0,
                                                 lineage_progen_chain)]

      # For those that have no incursions and are also the minimum case replace
      # with the next most likely progen, that is not already in the current chains
      candidate_links <- links_all[id_case %in% lins_to_fix[i] & membership != membership_progen]

      # Also filter to those that in chain with same (or totally unsampled lineages)
      candidate_links <- candidate_links[lineage_chain * lineage_progen_chain == 0 | lineage_chain * lineage_progen_chain == lineage_chain^2]
      fixed_links <- candidate_links[candidate_links[, .I[which.max(links)], by = c("id_case")]$V1]

      # if none then set to NA (prob & links)
      incs <- ifelse(nrow(fixed_links) == 0, lins_to_fix[i], 0)
      set_incs <- links_consensus[id_case %in% incs]
      set_incs[, c("id_progen", "links", "prob") := .(NA, NA, NA)]

      # bind them together
      links_consensus<-
        rbindlist(list(links_consensus[!(id_case %in% lins_to_fix[i])],
                       fixed_links, set_incs),
                  fill = TRUE)

      # clean links_consensus & links_all
      links_consensus[, c("membership", "membership_progen", "lineage_chain", "lineage_progen_chain") := NULL]
      links_all[, c("membership", "membership_progen", "lineage_chain", "lineage_progen_chain") := NULL]

      # update membership_dt
      membership_dt <- get_membership(links_consensus)

    }
  }

  # then fix the loops
  list2env(find_loops_to_fix(links_consensus, case_dates,
                             known_progens, max_cycles = max_cycles), envir = environment())

  # set links to fix to NA
  links_consensus[id_case %in% loops_to_fix]$id_progen <- NA
  nfixes <- length(loops_to_fix)

  if(nfixes > 0) {

    for(i in seq_len(nfixes)) {

      # Join up links with updated membership
      links_all <- membership_dt[links_all, on = "id_case"]
      setnames(membership_dt, c("membership", "id_case", "lineage_chain"),
               c("membership_progen", "id_progen", "lineage_progen_chain"))
      links_all <- membership_dt[links_all, on = "id_progen"]
      links_all[, membership_progen := ifelse(is.na(membership_progen), 0,
                                              membership_progen)]
      links_all[, lineage_progen_chain := ifelse(is.na(lineage_progen_chain), 0,
                                           lineage_progen_chain)]

      # For those that have no incursions and are also the minimum case replace
      # with the next most likely progen, that is not already in the current chains
      candidate_links <- links_all[id_case %in% loops_to_fix[i] & membership != membership_progen]

      # Also filter to those that in chain with same (or totally unsampled lineages)
      candidate_links <- candidate_links[lineage_chain * lineage_progen_chain == 0 | lineage_chain * lineage_progen_chain == lineage_chain^2]
      fixed_links <- candidate_links[candidate_links[, .I[which.max(links)], by = c("id_case")]$V1]

      # if none then set to NA (prob & links)
      incs <- ifelse(nrow(fixed_links) == 0, loops_to_fix[i], 0)
      set_incs <- links_consensus[id_case %in% incs]
      set_incs[, c("id_progen", "links", "prob") := .(NA, NA, NA)]

      # bind them together
      links_consensus<-
        rbindlist(list(links_consensus[!(id_case %in% loops_to_fix[i])],
                       fixed_links, set_incs),
                  fill = TRUE)

      # clean links_consensus & links_all
      links_consensus[, c("membership", "membership_progen", "lineage_chain", "lineage_progen_chain") := NULL]
      links_all[, c("membership", "membership_progen", "lineage_chain", "lineage_progen_chain") := NULL]

      # update membership
      membership_dt <- get_membership(links_consensus)

    }
  }

  # update final membership & check loops if still some unresolved
  links_consensus <- membership_dt[links_consensus, on = "id_case"]
  lins_to_fix <- check_lineages(links_consensus)
  loops <- check_loops(links_consensus)

  # Test to make sure no more to_fix and warn if there are!
  if(length(loops) > 0) {

    warning(
      "There are still loops in the transmission tree, indicating that case date
        uncertainty may be too high to indentify introductions and differentiate chains.")

  }

  if(nrow(lins_to_fix) > 0) {
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

#' Internal function for finding loops to fix
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
find_loops_to_fix <- function(links_consensus, case_dates,
                              known_progens, max_cycles = 100) {

  # build undirected & find the loops (which_multiple) & any cycles (girth)
  loops <- check_loops(links_consensus)
  loops_to_fix_all <- character(0)
  iter <- 0
  membership_dt <- get_membership(links_consensus)

  while(iter < max_cycles & length(loops) > 0) {

    links_consensus <- links_consensus[membership_dt, on = "id_case"]

    if(is.null(case_dates)) {
      stop("Need to pass case_dates data.table to fix loops by date.")
    }

    # join with case dates
    links_consensus <- links_consensus[case_dates, on = "id_case"]

    # Filter to the ones with loops
    # & filter out any that are known from contact tracing
    loops_to_fix <- links_consensus[id_case %in% loops & !(id_case %in% known_progens)]

    # Find the minimum selector case date
    loops_to_fix <- loops_to_fix[loops_to_fix[, .I[which.min(symptoms_started)],
                                              by = c("membership")]$V1]
    loops_to_fix <- loops_to_fix[order(symptoms_started)]$id_case
    loops_to_fix_all <- c(loops_to_fix_all, as.numeric(unique(loops_to_fix)))

    # Update links consensus and get new membership
    links_consensus[id_case %in% loops_to_fix]$id_progen <- NA

    loops <- check_loops(links_consensus)

    membership_dt <- get_membership(links_consensus)

    iter <- iter + 1
  }

  return(list(membership_dt = membership_dt, loops_to_fix = loops_to_fix_all))
}


#' Internal function for finding mismatched lineages to fix
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
find_lins_to_fix <- function(links_consensus, known_progens, selector) {

  # build undirected & find the loops (which_multiple) & any cycles (girth)
  gr <- get_gr(links_consensus)

  # Get the chain membership
  V(gr)$membership <- components(gr)$membership
  membership_dt <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")))
  links_consensus <- links_consensus[membership_dt, on = "id_case"]

  # Resolve the phylogeny
  if(length(unique(links_consensus[lineage != 0]$lineage)) > 1) {
    lins_to_fix <- find_lineages(gr, links_consensus, known_progens,
                                 selector = selector)
  } else {
    lins_to_fix <- NULL
  }

  lins_to_fix <- as.numeric(unique(lins_to_fix))

  # Update links consensus and get new membership_dt
  links_consensus[id_case %in% lins_to_fix]$id_progen <- NA

  membership_dt <- get_membership(links_consensus)

  return(list(membership_dt = membership_dt, lins_to_fix = lins_to_fix))
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
find_lineages <- function(gr, links, known_progens,
                          selector = c("prob", "prob_scale")) {

  # selector
  selector <- match.arg(selector)

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

  if(nrow(reassign) > 0) {

    # get unique pairs only
    reassign <- reassign[!duplicated(t(apply(reassign[,
                                                      c("id_case", "i.id_case")],
                                             1, sort))), ]

    pths <- get_edge_dt(gr, lins = reassign)

    # when you get the paths you lose the directionality
    # so need to join back up with the actual links
    check <- rbind(pths[, .(id_progen = as.numeric(from),
                            id_case = as.numeric(to), row_id)],
                   pths[, .(id_progen = as.numeric(to),
                            id_case = as.numeric(from), row_id)])
    pths <- links[check, on = c("id_case", "id_progen"),
                                             nomatch = NULL]

    # Filter out any edges that are known (from tracing)
    pths <- pths[!(id_case %in% known_progens)]

    # Filter to the most frequently selected ones for each row id
    pths[, freq := .N, by = c("id_case", "id_progen")]

    pths[, max_freq := max(freq), by = "row_id"]
    pths <- pths[freq == max_freq]

    # Select the edge to break: minimum probs
    pths$selector <- pths[, get(selector)]
    to_fix <- pths[pths[, .I[which.min(selector)],
                          by = c("row_id")]$V1]

  } else {
    to_fix <- NULL
  }

  return(unique(to_fix[order(selector)]$id_case))
}

#' Get the number of lineages in each chain
#'
#' @param links the tree or consensus tree
#'
#' @return
#' @export
#'
check_lineages <- function(links) {

  # Filter to chains that have multiple lineages per chain
  multilins <- links[lineage != 0][, .(check = length(unique(lineage))),
                                   by = "membership"][check > 1]
  return(multilins)
}

#' Internal function for getting edges between mismatched lineages
#'
#' @param gr the graph of the tree
#' @param lins the lineages
#'
#' @return
#'
get_edge_dt <- function(gr, lins) {

  sps <-
    suppressWarnings(
      lapply(seq_len(nrow(lins)),
             function(x) {
               shortest_paths(gr,
                       from = as.character(lins$id_case[x]),
                       to = as.character(lins$i.id_case[x]),  mode = "all",
                       output = "epath")$epath[[1]]
               }
      )
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
#' @param links the tree or consensus tree
#'
#' @importFrom igraph V subgraph.edges E count_multiple girth components
#'  vertex_attr graph_from_data_frame
#'
#' @return
#'
get_membership <- function(links) {

  # build directed graph
  gr <- get_gr(links)

  # Get the chain membership
  V(gr)$membership <- components(gr)$membership
  membership_dt <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")),
                           lineage_chain = as.numeric(vertex_attr(gr, "lineage")))

  # will only be one non zero lineage per chain as this is after mismatches are broken
  membership_dt[, lineage_chain := sum(unique(lineage_chain)), by = "membership"]

  return(membership_dt)
}

#' Build the consensus tree (i.e. the tree with the highest proportion of consensus links)
#'
#' @param links_consensus output from `build_consensus_links`
#' @param ttrees  bootstrapped trees from `boot_trees`
#' @param links_all the links summarized over the trees from `build_all_links`
#' @param type how to summarize the trees, either majority rule or Maximum clade
#'  credibility-ish
#' @param output whether to output the consensus tree or the sim number of the
#'  consensus tree
#'
#' @return a data.table with the consensus tree (with score for each link,
#'  1 if the consensus link, 0 if not.
#' @export
#'
build_consensus_tree <- function(links_consensus, ttrees, links_all,
                                 type = c("majority", "mcc"),
                                 output = c("tree", "sim")) {

  type <- match.arg(type)
  output <- match.arg(output)

  tree_consensus <- links_all[ttrees, on = c("id_progen", "id_case")]

  if(type == "majority") {
    sim_scores <- ttrees[links_consensus, on = c("id_case", "id_progen")][, score := 1][, .(score = sum(score)), by = "sim"]
    } else {
    # Join with links all and take the product of those
    sim_scores <- tree_consensus[, .(score = prod(prob, na.rm = TRUE)),
                                 by = "sim"]
  }

  if(output == "tree") {
    tree_consensus <- tree_consensus[sim == sim_scores$sim[which.max(sim_scores$score)]]
    best_pairs <- paste0(links_consensus$id_case, "_", links_consensus$id_progen)
    tree_consensus[, score := fifelse(paste0(id_case, "_", id_progen) %in% best_pairs, 1, 0)]

    return(tree_consensus)

  } else {
    return(sim_scores$sim[which.max(sim_scores$score)])
  }

}

