\name{zic.ssvs}

\alias{zic.ssvs}

\title{SSVS for Zero-Inflated Count Models} 

\description{
\code{zic.ssvs} applies SSVS to zero-inflated count models
}

\usage{
zic.ssvs(formula, data, tausq0bar, tausq1bar, omegasq0bar, omegasq1bar,
         ebar, fbar, pbar, qbar, n.burnin, n.mcmc, n.thin)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit specifying the response variable and covariates.}
  \item{data}{A data frame in which to interpret the variables in \code{formula}.}
  \item{tausq0bar}{The ``small'' prior variances of \eqn{\beta}{beta}, a vector of length equal to the number of covariates.}
  \item{tausq1bar}{The ``large'' prior variances of \eqn{\beta}{beta}, a vector of length equal to the number of covariates.}
  \item{omegasq0bar}{The ``small'' prior variances of \eqn{\delta}{delta}, a vector of length equal to the number of covariates.}
  \item{omegasq1bar}{The ``large'' prior variances of \eqn{\delta}{delta}, a vector of length equal to the number of covariates.}
  \item{ebar}{The shape parameter for the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{fbar}{The inverse scale parameter the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{pbar}{The prior inclusion probabilties of \eqn{\beta}{beta}, a vector of length equal to the number of covariates.}
  \item{qbar}{The prior inclusion probabilties of \eqn{\delta}{delta}, a vector of length equal to the number of covariates.}
  \item{n.burnin}{Number of burn-in iterations of the sampler.}
  \item{n.mcmc}{Number of iterations of the sampler.}
  \item{n.thin}{Thinning interval.}
}

\value{
   A list containing the following elements:
   \item{beta}{The posterior draws for \eqn{\beta}{beta}.}
   \item{delta}{The posterior draws for \eqn{\delta}{delta}.}
   \item{sigma2}{The posterior draws for \eqn{\sigma^2}{sigma^2}.}
   \item{gamma}{The posterior draws for \eqn{\gamma}{gamma}.}
   \item{kappa}{The posterior draws for \eqn{\kappa}{kappa}.}
}
    
\details{
  The considered zero-inflated count model is given by
  \deqn{y_i^* \sim \mathrm{Poisson}[\exp(\eta^*_i)],}{%
        y*_i ~ Poisson[exp(eta*_i)],}
  \deqn{\eta^*_i = x_i'\beta + \varepsilon_i,\; \varepsilon_i \sim \mathrm{N}(0,\sigma^2),}{%
        eta*_i = x_i' * beta + epsilon_i, epsilon_i ~ N( 0, sigma^2 ),}
  \deqn{d_i^* = x_i'\delta + \nu_i,\; \nu_i \sim \mathrm{N}(0,1),}{%
        d*_i = x_i' * delta + nu_i, nu_i ~ N( 0, 1 ),}
  \deqn{y_i = 1(d_i^*>0)y_i^*,}{%
        y_i = 1(d*_i>0) y*_i,} 
  where \eqn{y_i} and \eqn{x_i} are observed. The SSVS prior is given by
  \deqn{\beta_j \sim (1-\gamma_j)\mathrm{N}\left(0, \overline{\tau}_{0j}^2\right)+\gamma_j\mathrm{N}\left(0, \overline{\tau}_{1j}^2\right), \quad j=1,\ldots,k,}{%
        beta_j ~ (1-gamma_j) N(0,taubarsq_0j) + gamma_j N(0,taubarsq_1j), j=1,...k,}
  \deqn{\delta_j \sim (1-\kappa_j)\mathrm{N}\left(0, \overline{\omega}_{0j}^2\right)+\kappa_j\mathrm{N}\left(0,\overline{\omega}_{1j}^2\right), \quad j=1,\ldots,k,}{delta_j ~ (1-kappa_j) N(0,omegabarsq_0j) + kappa_j N(0,omegabarsq_1j), j=1,...k,}
  \deqn{\mathrm{P}(\gamma_j=1) = 1-\mathrm{P}(\gamma_j=0) = \overline p_j, \quad j=1,\ldots,k,}{%
        P(gamma_j=1) = 1-P(gamma_j=1) = pbar_j,}
  \deqn{\mathrm{P}(\kappa_j=1) = 1-\mathrm{P}(\kappa_j=0) = \overline q_j, \quad j=1,\ldots,k,}{%
        P(kappa_j=1) = 1-P(kappa_j=0) = qbar_j,}
  \deqn{\sigma^2 \sim \textrm{Inv-Gamma}\left(\bar e,\bar f\right).}{%
        sigma^2 ~ Inv-Gamma(ebar,fbar).}

  The sampling algorithm described in Jochmann (2009) is used.
}

\examples{
# library( zic )
# data( docvisits )
# formula <- docvisits ~ age + agesq + health + handicap + hdegree +
#            married + schooling + hhincome + children + self +
#            civil + bluec + employed + public + addon 

# prior parameters for preliminary run
# bbar <- rep( 10.0, 16 )
# dbar <- rep( 10.0, 16 )
# ebar <- 3.0
# fbar <- 2.0

# set seed and run preliminary MCMC sampler
# set.seed(1)
# prelim <- zic( formula, docvisits, bbar, dbar, ebar, fbar, 10000, 100000, 10 )

# prior parameters SSVS run
# varbeta <- apply( prelim$beta, 2, var )
# vardelta <- apply( prelim$delta, 2, var )
# tausq0bar <- 0.01 * varbeta
# tausq1bar <- 100.0 * varbeta
# omegasq0bar <- 0.01 * vardelta
# omegasq1bar <- 100.0 * vardelta
# ebar <- 3.0
# fbar <- 2.0
# pbar <- rep( 0.5, 16 )
# qbar <- rep( 0.5, 16 )

# set seed and run preliminary MCMC sampler
# set.seed(2)
# results <- zic.ssvs( formula, docvisits, tausq0bar, tausq1bar,
#                      omegasq0bar, omegasq1bar, ebar, fbar,
#                      pbar, qbar, 10000, 100000, 10 )

# print posterior inclusion probabilities means for delta
# apply( results$kappa, 2, mean )
}

\references{
  Jochmann, M. (2009). ``What Belongs Where? Variable Selection for Zero-Inflated Count Models with an Application to the Demand for
  Health Care''. Available at:
  \url{http://personal.strath.ac.uk/markus.jochmann}.
}

