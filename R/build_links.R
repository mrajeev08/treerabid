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
#' @param all_links output from `build_all_links`
#' @param case_dates data.table with two columns id_case (id of case) and symptoms_started
#'  (date symptoms started of case without uncertainty).
#'
#' If it is not clear which case seeded a chain, the case with the earliest start
#' date is selected as the incursion.
#'
#' To do: add in parts to break a loop: if you assign to next closest progen,
#' then you get lowest # of intros & higher Re ests, if you break the loops
#' by known case date, then you'll get highest # of intros & lower Re ests.
#'
#' @return a data.table with the case pairs most often linked in transmission tree and their
#'  probability (i.e. prop of times a given case was selected as a progenitor)
#' @export
#'
build_consensus_links <- function(all_links, case_dates) {

  # Get the consensus links
  links_consensus <- links_all[links_all[, .I[which.max(links)], by = c("id_case")]$V1] # returns first max


  # fix_loops = c("chain_break",
  #               "chain_replace",
  #               "all_break",
  #               "all_replace"),
  # min_prob = 0.05
  # Option to break loops
  # get id_progen_id (the progenitor id of the progenitor id)
  # while any id_case == id_progen_id are the loops
  # break them by either assigning second most likely
  # and updating id_progen_id
  # or just breaking by min to known case date (if tied just selects one)

  # only break loops to distinguish incursions

  gr <- graph_from_data_frame(d = links_consensus[, c("id_case",
                                                      "id_progen")][!is.na((id_progen))],
                              vertices = links_consensus[, "id_case"],
                              directed = TRUE)
  V(gr)$membership <- components(gr)$membership
  membership <- data.table(membership = as.numeric(vertex_attr(gr, "membership")),
                           id_case = as.numeric(vertex_attr(gr, "name")))

  # This also accounts for indirect loops which you might get given really large uncertainty!
  links_consensus <- links_consensus[membership, on = "id_case"][case_dates, on = "id_case"]
  links_consensus[, c("test", "min_case") := .(sum(is.na(id_progen)),
                                               min_case = id_case[which.min(symptoms_started)]),
                  by = "membership"]

  # option to reassign to second most likely and check loops
  links_consensus[, id_progen := ifelse(id_case %in% min_case & test == 0,
                                        NA,
                                        id_progen)]
  # set prob to NA!
  check <- links_consensus[, .(check = sum(is.na(id_progen))), by = "membership"]

  if(any(check$check != 1)) {
    warning("Multiple incursions are linked to the same chain or
             some chains are lacking an origin: there may be too much uncertainty
             in your dates to reliably distinguish chains!")
  }

  return(links_consensus[, -c("membership", "test", "min_case")])
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

