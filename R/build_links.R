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
#'
#' @return a data.table with the case pairs most often linked in transmission tree and their
#'  probability (i.e. prop of times a given case was selected as a progenitor)
#' @export
#'
build_consensus_links <- function(all_links) {
  # Get the consensus links
  links_consensus <- links_all[links_all[, .I[which.max(links)], by = c("id_case")]$V1] # returns first max
  return(links_consensus)
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
  best_pairs <- paste0(links_consensus$id_case, "_", links_consensus$id_progen)
  ttrees$score <- ifelse(paste0(ttrees$id_case, "_", ttrees$id_progen) %in% best_pairs, 1, 0)
  sim_scores <- ttrees[, .(score = sum(score)), by = "sim"]
  tree_consensus <- ttrees[sim == sim_scores$sim[which.max(sim_scores$score)]]
  return(ttrees)
}
