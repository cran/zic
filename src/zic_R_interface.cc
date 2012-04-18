#include <RcppArmadillo.h>
#include "zicmodel.h"

extern "C" SEXP zic_sample( SEXP _y, SEXP _X, 
			    SEXP _a0, SEXP _b0, SEXP _gbeta0,  SEXP _hbeta0, SEXP _nubeta0, SEXP _e0, SEXP _f0, 
			    SEXP _c0, SEXP _d0, SEXP _gdelta0, SEXP _hdelta0, SEXP _nudelta0,
			    SEXP _svs0, SEXP _nburnin, SEXP _nmcmc, SEXP _nthin, SEXP _tune ) 
{
  try 
     {
       // data matrices and dimensions
       Rcpp::IntegerVector yR( _y );  
       Rcpp::NumericMatrix XR( _X );  
       const int n = XR.nrow();
       const int k = XR.ncol() - 1;
       const ivec y( yR.begin(), yR.size(), false );
       const mat X( XR.begin(), n, k+1, false );   

       // prior
       SpikeSlabPrior betaprior, deltaprior;
       
       betaprior.Va = Rcpp::as<double>( _a0 );
       betaprior.Vb = Rcpp::as<double>( _b0 );
       betaprior.g  = Rcpp::as<double>( _gbeta0 );
       betaprior.h  = Rcpp::as<double>( _hbeta0 );
       betaprior.nu = Rcpp::as<double>( _nubeta0 );
       betaprior.e   = Rcpp::as<double>( _e0 );
       betaprior.f   = Rcpp::as<double>( _f0 );
       betaprior.svs = Rcpp::as<bool>( _svs0 );
       
       deltaprior.Va   = Rcpp::as<double>( _c0 );
       deltaprior.Vb   = Rcpp::as<double>( _d0 );
       deltaprior.g    = Rcpp::as<double>( _gdelta0 );
       deltaprior.h    = Rcpp::as<double>( _hdelta0 );
       deltaprior.nu   = Rcpp::as<double>( _nudelta0 );
       deltaprior.e    = -9.9;
       deltaprior.f    = -9.9;
       deltaprior.svs  = Rcpp::as<bool>( _svs0 );
       
       // mcmc iterations
       const int nburnin = Rcpp::as<int>( _nburnin );
       const int nmcmc = Rcpp::as<int>( _nmcmc );
       const int nthin = Rcpp::as<int>( _nthin );

       // tuning parameter
       const double tune = Rcpp::as<double>( _tune );
              
       // memories
       vec alphamem( nmcmc/nthin ), gammamem( nmcmc/nthin ), sigma2mem( nmcmc/nthin ), omegabetamem( nmcmc/nthin ), omegadeltamem( nmcmc/nthin ), accmem( nmcmc/nthin );
       mat betamem( nmcmc/nthin, k ), deltamem( nmcmc/nthin, k ); 
       imat Ibetamem( nmcmc/nthin, k ), Ideltamem( nmcmc/nthin, k ), yrepmem( nmcmc/nthin, n );
            
       // initialize model and sample
       ZicModel m( y, X, betaprior, deltaprior, tune );
       m.sample( nburnin, nmcmc, nthin, alphamem, betamem, gammamem, deltamem, sigma2mem, Ibetamem, Ideltamem, omegabetamem, omegadeltamem, yrepmem, accmem );
       
       // return memories
       Rcpp::List memlist;

       if( betaprior.svs )
	 {
	   return Rcpp::List::create( Rcpp::Named( "alpha" )       = alphamem,
	                              Rcpp::Named( "beta" )        = betamem,
	                              Rcpp::Named( "gamma" )       = gammamem,
		                      Rcpp::Named( "delta" )       = deltamem,
                                      Rcpp::Named( "sigma2" )      = sigma2mem,
                                      Rcpp::Named( "I.beta" )      = Ibetamem,
		                      Rcpp::Named( "I.delta" )     = Ideltamem,
		                      Rcpp::Named( "omega.beta" )  = omegabetamem,
		                      Rcpp::Named( "omega.delta" ) = omegadeltamem,
		                      Rcpp::Named( "yrep" )        = yrepmem, 
		                      Rcpp::Named( "acc" )         = accmem         );
	 }
       else
         {
	   return Rcpp::List::create( Rcpp::Named( "alpha" )       = alphamem,
	                              Rcpp::Named( "beta" )        = betamem,
	                              Rcpp::Named( "gamma" )       = gammamem,
		                      Rcpp::Named( "delta" )       = deltamem,
                                      Rcpp::Named( "sigma2" )      = sigma2mem,
			              Rcpp::Named( "yrep" )        = yrepmem, 
		                      Rcpp::Named( "acc" )         = accmem         );
	 }
     } 
  catch( std::exception &ex ) 
    {
      forward_exception_to_r( ex );
    } 
  catch( ... ) 
    { 
      ::Rf_error( "C++ exception (unknown reason)" ); 
    }
  return R_NilValue; 
}





