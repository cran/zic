#include "zicmodel.h"

ZicModel::ZicModel( const ivec& _y, const mat& _X, const SpikeSlabPrior& _betaprior, const SpikeSlabPrior& _deltaprior, const double _tune )
  : y(_y), X(_X), n(y.n_rows), n0(accu(y==0)), k(X.n_cols-1), eta(zeros<vec>(n)), dstar(zeros<vec>(n)), yrep(zeros<ivec>(n)), 
    beta(X.n_cols-1,_betaprior), delta(X.n_cols-1,_deltaprior), sigma2(1.0), e0(_betaprior.e), f0(_betaprior.f), tune(_tune)
{
}; 
 
void ZicModel::sample( const int nburnin, const int nmcmc, const int nthin, 
                       vec& alphamem, mat& betamem, vec& gammamem, mat& deltamem, vec& sigma2mem, imat& Ibetamem, imat& Ideltamem, vec& omegabetamem, vec& omegadeltamem, imat& yrepmem, vec& accmem )
{
  Rcout << "Gibbs sampler is running.\n";    
  Rcout << "Progress: [                    ]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
  
  int save_idx = 0;
  for( int i=1; i<=nburnin+nmcmc; i++ )
    {
      if( i%((nburnin+nmcmc)/20) == 0 )
	Rcout << "#";

      beta.update( eta, X, sigma2, rnd );
      sigma2_update();
      delta.update( dstar, X, 1.0, rnd );
      latvar_update();     
      
      if( ( i>nburnin ) && ( i%nthin==0 ) )
	{
	  alphamem(save_idx) = beta()(0);
	  betamem.row(save_idx) = beta().subvec(1,k).t(); 
	  gammamem(save_idx) = delta()(0);
	  deltamem.row(save_idx) = delta().subvec(1,k).t(); 
	  sigma2mem(save_idx) = sigma2;
	  Ibetamem.row(save_idx) = conv_to<ivec>::from( ( beta.get_tau() > 0.99 ) ).t(); 
	  Ideltamem.row(save_idx) = conv_to<ivec>::from( ( delta.get_tau() > 0.99 ) ).t(); 
	  omegabetamem(save_idx) = beta.get_omega();
	  omegadeltamem(save_idx) = delta.get_omega();
	  accmem(save_idx) = acc;
	  yrep_update();
	  yrepmem.row(save_idx) = yrep.t();

	  save_idx++;
	}
    }
  
  Rcout << "]\nGibbs sampler is finished.\n";
}

void ZicModel::latvar_update()   // checked 
{
  const vec etaold = eta;
  const vec Xbeta = X * beta();
  const vec Xdelta = X * delta();

  for( int i=0; i<n0; i++ )
    {
      // update eta
      double p = Rcpp::stats::pnorm_0( Xdelta(i), true, false );   
      double eta_prob = eta(i) + tune * rnd.t( 0.0, sigma2, 5.0 );
      double logp_prob = log( 1.0 - p + p * exp(-exp(eta_prob)) ) - 0.5 * (eta_prob-Xbeta(i)) * (eta_prob-Xbeta(i)) / sigma2;
      double logp_old  = log( 1.0 - p + p * exp(-exp(eta(i)  )) ) - 0.5 * (eta(i)  -Xbeta(i)) * (eta(i)  -Xbeta(i)) / sigma2;
      if( exp( logp_prob - logp_old ) > rnd.uniform() )
	eta(i) = eta_prob;
    
      // update dstar
      if( rnd.uniform() > (1.0 - p) / ( 1.0 - p + p * exp(-exp(eta(i))) ) )
	dstar(i) = rnd.normal_lt( 0.0, Xdelta(i), 1.0 );
      else
	dstar(i) = rnd.normal_rt( 0.0, Xdelta(i), 1.0 );
    }

  for( int i=n0; i<n; i++ )
    {
      // update eta
      double eta_prob = eta(i) + tune * rnd.t( 0.0, sigma2, 5.0 );
      double logp_prob = -exp(eta_prob) + eta_prob * y(i) - 0.5 * (eta_prob-Xbeta(i)) * (eta_prob-Xbeta(i)) / sigma2;
      double logp_old  = -exp(eta(i)  ) + eta(i)   * y(i) - 0.5 * (eta(i)  -Xbeta(i)) * (eta(i)  -Xbeta(i)) / sigma2;
      if( exp( logp_prob - logp_old ) > rnd.uniform() )
	eta(i) = eta_prob;
 
      // update dstar
      dstar(i) = rnd.normal_lt( 0.0, Xdelta(i), 1.0 );
    }

  acc = 1.0 - mean( conv_to<vec>::from( eta == etaold ) );
}

void ZicModel::sigma2_update()   // checked
{
  const vec resid = eta - X * beta();
  const double ssr = dot( resid, resid );
  sigma2 = rnd.invGamma( e0 + 0.5 * n, f0 + 0.5 * ssr );
}

/*
void ZicModel::latvar_update2()
{
  const vec etaold( eta );

  const vec Xbeta = X * beta();
  const vec Xdelta = X * delta();

  for( int i=0; i<n0; i++ )
    {
      // update eta

      if( dstar(i) <= 0 )
	{
	  eta(i) = Xbeta(i) + sqrt(sigma2) * rnd.normal();
	}
      else
	{
	  double eta_prob = eta(i) + tune * rnd.t( 0.0, sigma2, 10.0 );
	  double p_prob = exp( logp( eta_prob, y(i), Xbeta(i), sigma2 ) );
	  double p_old = exp( logp( eta(i), y(i), Xbeta(i), sigma2 ) );

	  if( ( p_prob / p_old ) > rnd.uniform() )
	    eta(i) = eta_prob;
	}

      // update dstar

      double p = Rcpp::stats::pnorm_0( Xdelta(i), true, false );   
      if( rnd.uniform() > (1-p)/(1-p+p*exp(-exp(eta(i)))))
	dstar(i) = rnd.normal_lt( 0.0, Xdelta(i), 1.0 );
      else
	dstar(i) = rnd.normal_rt( 0.0, Xdelta(i), 1.0 );
    }
  for( int i=n0; i<n; i++ )
    {
      // update eta

      double eta_prob = eta(i) + tune * rnd.t( 0.0, sigma2, 10.0 );
      double p_prob = exp( logp( eta_prob, y(i), Xbeta(i), sigma2 ) );
      double p_old = exp( logp( eta(i), y(i), Xbeta(i), sigma2 ) );

      if( ( p_prob / p_old ) > rnd.uniform() )
	eta(i) = eta_prob;
 
      // update dstar
      
      dstar(i) = rnd.normal_lt( 0.0, Xdelta(i), 1.0 );
    }

  acc = 1.0 - mean( conv_to<vec>::from( eta == etaold ) );
}
*/

void ZicModel::yrep_update()
{
  vec epsdraw(n);
  rnd.normal( epsdraw );
  vec etastardraw(n);
  etastardraw = X * beta() + sigma2 * epsdraw;
  ivec ystardraw(n);
  for( int i=0; i<n; i++ )
    ystardraw(i) = rnd.poisson( exp( etastardraw(i) ) );
  vec nudraw(n);
  rnd.normal( nudraw );
  vec dstardraw(n);
  dstardraw = X * delta() + nudraw;
  yrep = ( dstardraw > 0 ) % ystardraw;
}

