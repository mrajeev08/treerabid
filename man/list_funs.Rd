% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct.R
\name{list_funs}
\alias{list_funs}
\title{Helper for functions to export to foreach loop in boot_trees}
\usage{
list_funs(filename)
}
\arguments{
\item{filename}{path to the R script with functions. You still need to source
this function before calling boot_trees.}
}
\value{
a character vector of names functions in R script
}
\description{
Helper for functions to export to foreach loop in boot_trees
}
