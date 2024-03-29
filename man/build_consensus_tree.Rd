% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_links.R
\name{build_consensus_tree}
\alias{build_consensus_tree}
\title{Build the consensus tree (i.e. the tree with the highest proportion of consensus links)}
\usage{
build_consensus_tree(
  links_consensus,
  ttrees,
  links_all,
  type = c("majority", "mcc"),
  output = c("tree", "sim")
)
}
\arguments{
\item{links_consensus}{output from \code{build_consensus_links}}

\item{ttrees}{bootstrapped trees from \code{boot_trees}}

\item{links_all}{the links summarized over the trees from \code{build_all_links}}

\item{type}{how to summarize the trees, either majority rule or Maximum clade
credibility-ish}

\item{output}{whether to output the consensus tree or the sim number of the
consensus tree}
}
\value{
a data.table with the consensus tree (with score for each link,
1 if the consensus link, 0 if not.
}
\description{
Build the consensus tree (i.e. the tree with the highest proportion of consensus links)
}
