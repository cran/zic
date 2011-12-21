\name{zic}

\alias{zic}

\title{Bayesian Inference for Zero-Inflated Count Models} 

\description{
\code{zic} fits zero-inflated count models via Markov chain Monte Carlo methods.
}

\usage{
zic(formula, data, bbar, dbar, ebar, fbar, n.burnin, n.mcmc, n.thin)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit specifying the response variable and covariates.}
  \item{data}{A data frame in which to interpret the variables in \code{formula}.}
  \item{bbar}{The diagonal elements of the prior variance matrix of \eqn{\beta}{beta}, a vector of length equal to the number of covariates.}
  \item{dbar}{The diagonal elements of the prior variance matrix of \eqn{\delta}{delta}, a vector of length equal to the number of covariates.}
  \item{ebar}{The shape parameter for the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{fbar}{The inverse scale parameter the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{n.burnin}{Number of burn-in iterations of the sampler.}
  \item{n.mcmc}{Number of iterations of the sampler.}
  \item{n.thin}{Thinning interval.}
}

\value{
   A list containing the following elements:
   \item{beta}{The posterior draws for \eqn{\beta}{beta}.}
   \item{delta}{The posterior draws for \eqn{\delta}{delta}.}
   \item{sigma2}{The posterior draws for \eqn{\sigma^2}{sigma^2}.}
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
  where \eqn{y_i} and \eqn{x_i} are observed. The assumed prior distributions are
  \deqn{\beta \sim \mathrm{N}(0,\bar B) \quad \mathrm{with} \quad \bar B = \mathrm{diag}(\bar b_1,\ldots,\bar b_k),}{%
        beta ~ N(0,Bbar) with Bbar = diag(bbar1,...,bbark),}
  \deqn{\delta \sim \mathrm{N}(0,\bar D) \quad \mathrm{with} \quad \bar D = \mathrm{diag}(\bar d_1,\ldots,\bar d_k),}{%
        delta ~ N(0,Dbar) with Dbar = diag(dbar1,...,dbark),}
  \deqn{\sigma^2 \sim \textrm{Inv-Gamma}\left(\bar e,\bar f\right).}{%
        sigma^2 ~ Inv-Gamma(ebar,fbar).}

  The sampling algorithm developed by Jochmann (2009) is used.
}

\examples{
# library( zic )
# data( docvisits )

# prior parameters and formula
# bbar <- rep( 10.0, 16 )
# dbar <- rep( 10.0, 16 )
# ebar <- 3.0
# fbar <- 2.0
# formula <- docvisits ~ age + agesq + health + handicap + hdegree +
#            married + schooling + hhincome + children + self +
#            civil + bluec + employed + public + addon 

# set seed and run MCMC sampler
# set.seed(1)
# results <- zic( formula, docvisits, bbar, dbar, ebar, fbar, 10000, 100000, 10 )

# print posterior means for beta
# apply( results$beta, 2, mean )
}

\references{
  Jochmann, M. (2009). ``What Belongs Where? Variable Selection for Zero-Inflated Count Models with an Application to the Demand for
  Health Care''. Available at: \url{http://personal.strath.ac.uk/markus.jochmann}.
}
