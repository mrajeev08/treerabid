% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_links.R
\name{build_consensus_links}
\alias{build_consensus_links}
\title{Build consensus links between cases (i.e. most often selected progenitors for a given case)}
\usage{
build_consensus_links(
  links_all,
  case_dates,
  lineages = NULL,
  known_progens = NULL,
  max_cycles = 100,
  link_all = FALSE
)
}
\arguments{
\item{links_all}{output from \code{build_all_links}}

\item{case_dates}{data.table with two columns id_case (id of case) and symptoms_started
(date symptoms started of case without uncertainty)}

\item{lineages}{a data table with two columns, id_case and lineage, designating a lineage
assignment for each case, defaults to NULL which means trees wont be resolved to a phylogeny}

\item{known_progens}{a numeric vector of case ids for which progenitors are known}

\item{max_cycles}{how many cycles to look for (if lots of date uncertainty
you may want to increase this but it will likely slow things down a lot!)}

\item{link_all}{whether you want to try and force links between cases (that is
even if a case is often identified as the beginning of a chain, you want to
see if it can be assigned to another progenitor, i.e. when trying to link
up all cases to a sampled set of lineages)}
}
\value{
a data.table with the case pairs most often linked in transmission tree and their
probability (i.e. prop of times a given case was selected as a progenitor)
}
\description{
Build consensus links between cases (i.e. most often selected progenitors for a given case)
}
