% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MC_simulate.R
\name{simulate_MC}
\alias{simulate_MC}
\title{Simulate Metacommunity Dynamics}
\usage{
simulate_MC(
  patches,
  species,
  dispersal = 0.01,
  plot = TRUE,
  torus = FALSE,
  kernel_exp = 0.1,
  env1Scale = 500,
  temporal_autocorr = FALSE,
  timesteps = 200,
  burn_in = 100,
  initialization = 0,
  max_r = 5,
  min_env = 0,
  max_env = 1,
  env_niche_breadth = 0.5,
  optima_spacing = "random",
  intra = 1,
  min_inter = 0,
  max_inter = 0.8,
  comp_scaler = 0.05,
  extirp_prob = 0.05,
  extirp_all_pop = F,
  local_start = F,
  local_start_radius = 15,
  landscape,
  disp_mat,
  env.df,
  env_optima,
  int_mat
)
}
\arguments{
\item{patches}{number of patches to simulate}

\item{species}{number of species to simulate}

\item{dispersal}{dispersal probability between 0 and 1}

\item{plot}{option to show plot of landscape}

\item{torus}{whether to model the landscape as a torus}

\item{kernel_exp}{the exponential rate at which dispersal decreases as a function of the distance between patches}

\item{env1Scale}{environmental autocorrelation (e.g., 0.01 : strong autocorr. to larger values, e.g., 100 : no autocorr.)}

\item{timesteps}{number of timesteps to simulate}

\item{burn_in}{length of burn in period}

\item{initialization}{length of initial period before environmental change begins}

\item{max_r}{intrinsic growth rate in optimal environment, can be single value or vector of length species}

\item{min_env}{minium environmental optima}

\item{max_env}{minium environmental optima}

\item{env_niche_breadth}{standard deviation of environmental niche breadth, can be single value or vector of length species}

\item{optima_spacing}{"even" or "random" to specify how optima should be distributed}

\item{intra}{intraspecific competition coefficient, single value or vector of length species}

\item{min_inter}{min interspecific comp. coefficient}

\item{max_inter}{max interspecific comp. coefficient}

\item{comp_scaler}{value to multiply all competition coefficients by}

\item{extirp_prob}{probability of local extirpation for each population in each time step (should be a very small value, e.g. 0 or 0.002)}

\item{extirp_all_pop}{(logical - default False) should extirpation applied to all species in a patch (T) or to species/patch combination (F)}

\item{local_start}{(logical - default False) Are species initially clustered spatially ?}

\item{local_start_radius}{Radius in which initial species positions are initially clustered}

\item{landscape}{optional dataframe with x and y columns for patch coordinates}

\item{disp_mat}{optional matrix with each column specifying the probability that an individual disperses to each other patch (row)}

\item{env.df}{optional dataframe with environmental conditions with columns: env1, patch, time}

\item{env_optima}{optional values of environmental optima, should be a vector of length species}

\item{int_mat}{optional externally generated competition matrix}
}
\value{
list that includes metacommunity dynamics, landscape coordinates, environmental conditions, species environmental traits, dispersal matrix, and the competition matrix
}
\description{
Function is here updated
\itemize{
\item allow to remove temporal environmental heterogeneity (see temporal_autocorr param, default)
\item extirpation can be applied to patch instead of to species/patch (see extirp_all_pop param
}
}
\details{
Carefull, the default behaviour is now without temporal variability (and then, w/o temporal autocorrelation)
}
\examples{
# simulate_MC(patches = 6, species = 10, dispersal = 0.001, min_inter = 1, max_inter = 1, env_niche_breadth = 10)

}
\author{
Patrick L. Thompson, (updated M.D.)
}
