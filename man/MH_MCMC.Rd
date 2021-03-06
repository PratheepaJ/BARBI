% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MH_MCMC.R
\name{MH_MCMC}
\alias{MH_MCMC}
\title{MH_MCMC}
\usage{
MH_MCMC(itera, k, al_c, be_c, startvalue_lamda_r = 0)
}
\arguments{
\item{itera}{numeric, number of MCMC samples.}

\item{k}{numeric, observed count for a taxa in a sample.}

\item{al_c}{numeric, estimated shape parameter for the contaminant intensity.}

\item{be_c}{numeric, estimated scale parameter for the contaminant intensity.}

\item{startvalue_lamda_r}{numeric, initial value for the true intensity.}
}
\value{
An array of MCMC samples.
}
\description{
Metropolis Hastings MCMC
}
