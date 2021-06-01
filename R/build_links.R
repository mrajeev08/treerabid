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

  links_all <- ttrees[, .(links = .N), by = c("id_case", "id_progen")][, prob := links/N]
  return(links_all)

}

#' Build consensus links between cases (i.e. most often selected progenitors for a given case)
#'
#' @param links_all output from `build_all_links`
#' @param case_dates data.table with two columns id_case (id of case) and symptoms_started
#'  (date symptoms started of case without uncertainty).
#' @param fix_loops If there are loops in the transmission tree (i.e. indicating uncertainty in
#' who-infected-whom), the loops can either be broken by reassigning the progenitor
#' of the case with the earliest case date in the loop (fix_loops = "by_date", the default). or by breaking and replacing
#' links in the loop randomly (fix_loops = "random").
#' @param max_tries number of times to iterate through and break loops
#'
#' @return a data.table with the case pairs most often linked in transmission tree and their
#'  probability (i.e. prop of times a given case was selected as a progenitor)
#' @export
#'
build_consensus_links <- function(links_all, case_dates,
                                  fix_loops = c("by_date", "random"),
                                  max_tries = 100) {

  fix_loops <- match.arg(fix_loops)
  # Get the consensus links
  links_consensus <- links_all[links_all[, .I[which.max(links)], by = c("id_case")]$V1] # returns first max

  # spit out the loops and updated dt
  list2env(find_loops(links_consensus), envir = environment())
  niter <- 0

  # option to reassign to next likely and check loops
  while(length(loops) > 0 & niter <= max_tries) {

    niter <- niter + 1

    links_all <- membership[links_all, on = "id_case"]
    setnames(membership, c("membership", "id_case"), c("membership_progen", "id_progen"))
    links_all <- membership[links_all, on = "id_progen"]
    links_all[, membership_progen := ifelse(is.na(membership_progen), 0,
                                            membership_progen)]

    if(fix_loops == "by_date") {

      # join with case dates
      links_consensus <- links_consensus[case_dates, on = "id_case"]
      setnames(links_consensus, "symptoms_started", "selector")

    } else {

      # get a random variate to select by
      links_consensus$selector <- runif(nrow(links_consensus))

    }

    # Filter to the ones with loops
    to_fix <- links_consensus[id_case %in% loops]

    # Find the minimum selector (either by case date or randomly select link to break)
    to_fix <- to_fix[to_fix[, .I[which.min(selector)],
                            by = c("membership")]$V1]

    # For those that have no incursions and are also the minimum case replace
    # with the next most likely progen, that is not already in the current chains
    candidate_links <- links_all[id_case %in% to_fix$id_case & membership != membership_progen]
    fixed_links <- candidate_links[candidate_links[, .I[which.max(links)], by = c("id_case")]$V1]

    # if none then set to NA (prob & links)
    set_incs <- to_fix[!(id_case %in% fixed_links$id_case)]
    set_incs[, c("id_progen", "links", "prob") := .(NA, NA, NA)]

    # bind them together
    links_consensus<-
      rbindlist(list(links_consensus[!(id_case %in% fixed_links$id_case)],
                     fixed_links, set_incs),
                fill = TRUE)

    # clean links_consensus & links_all
    links_consensus[, c("membership", "selector", "membership_progen") := NULL]
    links_all[, c("membership", "membership_progen") := NULL]

    # update and recheck for loops
    list2env(find_loops(links_consensus), envir = environment())
  }

  # Test to make sure no more loops and warn if there are!
  if(length(loops) > 0) {

    warning("There are still loops in the transmission tree, indicating that case
            date uncertainty may be too high to indentify introductions and differentiate chains.
            You can also try increasing the value of max_tries.")

  }

  # Clean up the data.table
  return(links_consensus[, -"membership"])
}

#' Internal function for finding loops
#'
#' @param links_consensus
#'
#' @return
#' @importFrom igraph V subgraph.edges E count_multiple girth components
#'  vertex_attr graph_from_data_frame
#' @keywords internal
#'
find_loops <- function(links_consensus) {

  # build undirected & find the loops (which_multiple) & any cycles (girth)
  gr <- graph_from_data_frame(d = links_consensus[, c("id_case",
                                                      "id_progen")][!is.na((id_progen))],
                              vertices = links_consensus[, "id_case"],
                              directed = FALSE)
  loops <- names(V(subgraph.edges(gr, E(gr)[count_multiple(gr) > 1])))
  loops <- as.numeric(c(loops, names(girth(gr)$circle)))

  # Get the chain membership
  V(gr)$membership <- components(gr)$membership
  membership <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")))
  links_consensus <- links_consensus[membership, on = "id_case"]

  return(list(links_consensus = links_consensus, membership = membership,
              loops = loops))
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

