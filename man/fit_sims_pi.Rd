% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detection.R
\name{fit_sims_pi}
\alias{fit_sims_pi}
\title{Fit detection probability given observed time differences between linked cases
and a serial interval distribution}
\usage{
fit_sims_pi(
  t_diff,
  nsims = 1000,
  candidate_pis,
  si_fun,
  params,
  alpha = 0.001,
  known_kappas = NULL,
  seed = 132,
  sort = TRUE
)
}
\arguments{
\item{t_diff}{the observed time differences between linked cases}

\item{nsims}{the number of estimates to generate for the observed time differences}

\item{candidate_pis}{the candidate values of the detection probability to evaluate}

\item{si_fun}{a function for the serial interval with arguments N (the number
to draw) and params (a list of parameters for the function), see si_fun_lnorm for
an example.}

\item{params}{a list with parameters for the si_fun function to draw serial intervals}

\item{alpha}{probability, the value at which to constrain kappa (i.e. to determine max_kappa
for sim_generations), i.e. the probability of observing this kappa for a given pi is < alpha}

\item{known_kappas}{vector of known kappas (i.e. if some cases are traced, you know
kappa = 1 for these cases)}

\item{sort}{whether to sort the first column and t_diff to weight towards kappa = 1
to deal with higher sensitivity at higher reporting thresholds}
}
\value{
a vector of estimates of the detection probability generated by minimizing
the sum of squares between the observed and the expected (only looks at values passed
into candidate_pis)
}
\description{
Fit detection probability given observed time differences between linked cases
and a serial interval distribution
}
\examples{
# This example shows how to generate simulated data based on a detection estimate
# and a serial interval distribution and see whether the values can be recovered

\dontrun{
system.time({
tt <- rbindlist(lapply(runif(1000), function(z) {
    t_diff <- sim_times_pi(si_fun_lnorm, nobs = 500, params = treerabid::params_treerabid, alpha = 0.01,
                          pi = z)
    ests <- fit_sims_pi(t_diff, nsims = 5, candidate_pis = seq(0.01, 0.99, by = 0.01),
                       si_fun_lnorm, params = treerabid::params_treerabid, alpha = 0.01)
    data.table(true = z, estimated = ests)}))
 })

plot(tt$true, tt$estimated)
abline(a = 0, b = 1, col = "red") # the 1:1 line
}


}
