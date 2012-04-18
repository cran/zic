\name{zic.ssvs}

\alias{zic.ssvs}

\title{SVS for Zero-Inflated Count Models} 

\description{
\code{zic.ssvs} applies SVS to zero-inflated count models
}

\usage{
zic.ssvs(formula, data, a0, g0.beta, h0.beta, nu0.beta, e0, f0, 
                        c0, g0.delta, h0.delta, nu0.delta, 
                        n.burnin, n.mcmc, n.thin, tune = 1.0 )
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit specifying the response variable and covariates.}
  \item{data}{A data frame in which to interpret the variables in \code{formula}.}
  \item{a0}{The prior variance of \eqn{\alpha}{alpha}.}
  \item{g0.beta}{The shape parameter for the inverse gamma prior on \eqn{\kappa_j^\beta}{kappa_j^beta}.}
  \item{h0.beta}{The inverse scale parameter for the inverse gamma prior on \eqn{\kappa_j^\beta}{kappa_j^beta}.}
  \item{nu0.beta}{Prior parameter for the spike of the hypervariances for the \eqn{\beta_j}{beta_j}.}
  \item{e0}{The shape parameter for the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{f0}{The inverse scale parameter the inverse gamma prior on \eqn{\sigma^2}{sigma^2}.}
  \item{c0}{The prior variance of \eqn{\gamma}{gamma}.}	  
  \item{g0.delta}{The shape parameter for the inverse gamma prior on \eqn{\kappa_j^\delta}{kappa_j^delta}.}
  \item{h0.delta}{The inverse scale parameter for the inverse gamma prior on \eqn{\kappa_j^\delta}{kappa_j^delta}.}
  \item{nu0.delta}{Prior parameter for the spike of the hypervariances for the \eqn{\delta_j}{delta_j}.}
  \item{n.burnin}{Number of burn-in iterations of the sampler.}
  \item{n.mcmc}{Number of iterations of the sampler.}
  \item{n.thin}{Thinning interval.}
  \item{tune}{Tuning parameter of Metropolis-Hastings step.}
}

\value{
   A list containing the following elements:
   \item{alpha}{Posterior draws of \eqn{\alpha}{alpha} (coda mcmc object).}
   \item{beta}{Posterior draws of \eqn{\beta}{beta} (coda mcmc object).}
   \item{gamma}{Posterior draws of \eqn{\gamma}{gamma} (coda mcmc object).}	
   \item{delta}{Posterior draws of \eqn{\delta}{delta} (coda mcmc object).}
   \item{sigma2}{Posterior draws of \eqn{\sigma^2}{sigma^2} (coda mcmc object).}
   \item{I.beta}{Posterior draws of indicator whether \eqn{\tau_j^\beta}{tau_j^beta} is one (coda mcmc object).}
   \item{I.delta}{Posterior draws of indicator whether \eqn{\tau_j^\delta}{tau_j^delta} is one (coda mcmc object).}
   \item{omega.beta}{Posterior draws of \eqn{\omega^\beta}{omega^beta} (coda mcmc object).}
   \item{omega.delta}{Posterior draws of \eqn{\omega^\delta}{omega^delta} (coda mcmc object).}
   \item{yrep}{Posterior draws of replicated data sets \eqn{y^{rep}}.}
   \item{acc}{Acceptance rate of the Metropolis-Hastings step.}
}

\details{
  The considered zero-inflated count model is given by
  \deqn{y_i^* \sim \mathrm{Poisson}[\exp(\eta^*_i)],}{%
        y*_i ~ Poisson[exp(eta*_i)],}
  \deqn{\eta^*_i = \alpha + x_i'\beta + \varepsilon_i,\; \varepsilon_i \sim \mathrm{N}(0,\sigma^2),}{%
        eta*_i = x_i' * beta + epsilon_i, epsilon_i ~ N( 0, sigma^2 ),}
  \deqn{d_i^* = \gamma + x_i'\delta + \nu_i,\; \nu_i \sim \mathrm{N}(0,1),}{%
        d*_i = x_i' * delta + nu_i, nu_i ~ N(0,1),}
  \deqn{y_i = 1(d_i^*>0)y_i^*,}{%
        y_i = 1(d*_i>0) y*_i,} 
  where \eqn{y_i} and \eqn{x_i} are observed. The assumed prior distributions are
  \deqn{\alpha \sim \mathrm{N}(0,a_0),}{%
         alpha ~ N(0,a0),}
  \deqn{\beta_j\sim \mathrm{N}(0,\tau^\beta_j\kappa^\beta_j),}{%
        beta_j ~ N(0, tau_j^beta * kappa_j^beta),}
  \deqn{\kappa^\beta_j\sim\textrm{Inv-Gamma}(g_0^\beta,h_0^\beta),}{%
        kappa_j^\beta ~ Inv-Gamma(g0^beta,h0^beta),}
  \deqn{\tau_j^\beta \sim (1-\omega^\beta)\delta_{\nu^\beta_0}+\omega^\beta\delta_1,}{%
        tau_j^\beta ~ (1-\omega^beta) delta_(nu0^\beta)+ omega^beta delta_1,}
  \deqn{\omega^\beta\sim\mathrm{Uniform}(0,1),}{%
        omega^beta ~ U(0,1),}
  \deqn{\gamma \sim \mathrm{N}(0,c_0),}{%
         gamma ~ N(0,c0)}
  \deqn{\delta_j\sim \mathrm{N}(0,\tau^\delta_j\kappa^\delta_j),}{%
        delta_j ~ N(0, tau_j^delta * kappa_j^delta),}
  \deqn{\kappa^\delta_j\sim\textrm{Inv-Gamma}(g_0^\delta,h_0^\delta),}{%
        kappa_j^\delta ~ Inv-Gamma(g0^delta,h0^delta),}
  \deqn{\tau_j^\delta \sim (1-\omega^\delta)\delta_{\nu^\delta_0}+\omega^\delta\delta_1,}{%
        tau_j^\delta ~ (1-\omega^delta) delta_(nu_0^\delta)+ omega^delta delta_1,}
  \deqn{\sigma^2 \sim \textrm{Inv-Gamma}\left(e_0,f_0\right).}{%
         sigma^2 ~ Inv-Gamma(e0,f0).}


  The sampling algorithm described in Jochmann (2012) is used.
}

\references{
  Jochmann, M. (2012). ``What Belongs Where? Variable Selection for Zero-Inflated Count Models with an Application to the Demand for
  Health Care''. Available at: \url{http://www.staff.ncl.ac.uk/markus.jochmann}.
}


