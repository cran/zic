#include "zic.h"
 
using namespace std;
   
extern "C"    
{
  void zic( int* _y, double* _X, const int* _n, const int* _k, const int* _ssvs, const int* _nburnin, const int* _nmcmc, const int* _nthin, double* _p_b, double* _p_d, const double* _p_e, const double* _p_f, double* _p_p, double* _p_q, double* _bbetamem, double* _deltamem, double* _sigma2mem, double* _gammamem, double* _kappamem )
{  
  // dimensions
  int n = *_n, k = *_k, nburnin = *_nburnin, nmcmc = *_nmcmc, nthin = *_nthin;
  
  // ssvs: on or off
  int ssvs = *_ssvs;
  
  // y and X
  iarray y( _y, n );
  rmatrix X( _X, n, k );  
  
  // prior
  rmatrix p_b( _p_b, k, 2 ), p_d( _p_d, k, 2 );
  rvector p_p( _p_p, k ), p_q( _p_q, k );
  Prior prior;
  prior.b.resize(k,2); prior.b = p_b;
  prior.d.resize(k,2); prior.d = p_d;
  prior.e = *_p_e; prior.f = *_p_f;
  prior.p.resize(k); prior.p = p_p;
  prior.q.resize(k); prior.q = p_q;

  // memories
  rmatrix bbetamem( _bbetamem, nmcmc / nthin, k ), deltamem( _deltamem, nmcmc / nthin, k ), gammamem( _gammamem, nmcmc / nthin, k ), kappamem( _kappamem, nmcmc / nthin, k ); 
  rvector sigma2mem( _sigma2mem, nmcmc / nthin );
  
  // intitalize random generator
  Rndgen rnd;	
  
  // sample
  ZipModel model( y, X, prior, rnd );
  model.sample( ssvs, nburnin, nmcmc, nthin, bbetamem, deltamem, sigma2mem, gammamem, kappamem );
}   
      
}    	

void ZipModel::sample( const int ssvs, const int nburnin, const int nmcmc, const int nthin, rmatrix& bbetamem, rmatrix& deltamem, rvector& sigma2mem, rmatrix& gammamem, rmatrix& kappamem )
{
  Rprintf( "MCMC Sampler is running!\n" );
  int save_idx = 1;
  for( int i=1; i<=nburnin+nmcmc; i++ )
    {
      bbeta_sample();
      sigma2_sample();
      delta_sample();
      eta_sample();
      ystar_dstar_sample();
      
      if( ssvs )
	{
	  gamma_sample();
	  kappa_sample();
	}
      
      // save
      if( ( i>nburnin ) && ( (i-nburnin)%nthin==0 ) )
	{
	  bbetamem[save_idx] = bbeta;
	  deltamem[save_idx] = delta;
	  sigma2mem[save_idx] = sigma2;
	  if( ssvs )
	    {
	      gammamem[save_idx] = gamma;
	      kappamem[save_idx] = kappa;
	    }
	  save_idx++;
	}
    }
  Rprintf( "Finished!\n" );
}

void ZipModel::eta_sample( void )
{
  // setting parameters for ARS
  int n_max = 10;                  // maximum number of abscissae
  int n_start = 5;                 // starting number of abscissae
  double emax = 64.0;              // large number for which exp(.) can be calculated
  int lb = 0;                      // true, if lower bound exists
  double xlb = 0.0;                // value of the lower bound
  int ub = 0;                      // true, if upper bound exists
  double xub = 0.0;                // value of the upper bound
  int ifault;                      // return value, 0: successful initialisation, 1: not enough starting points
                                   // 2: ns is less than m, 3: no abscissae to left of mode (if lb = false),
                                   // 4: no abscissae to right of mode (if ub = false), 5: non-log-concavity detected
  int iwv[n_max+7];                // real-valued working space (size: n_max+7)
  int* iwv_p = iwv;                // respective pointer
  double rwv[6*n_max+15];          // real-valued working space (size: 6*n_max+15)
  double* rwv_p = rwv;             // respective pointer

  double draw;
  for( int i=1; i<=n; i++ )
    {
      // parameters of function
      ars::FnPar fnpar;
      fnpar.y = ystar[i];
      fnpar.mean = X[i] * bbeta;
      fnpar.var = sigma2;
	  
      // starting configuration
      double a[5];        // starting abscissae (array of size n_start)
      double f[5];        // log density at starting abscissae
      double fprime[5];   // gradients of log density at starting abcissae
	  
      *a = eta[i] - 1.0;
      fn_eval( a, &fnpar, f, fprime );
      while( fprime[0] < 0.1 )
	{
	  (*a)--;
	  fn_eval( a, &fnpar, f, fprime );
	}
      
      *(a+1) = eta[i] - 0.5;
      fn_eval( a+1, &fnpar, f+1, fprime+1 );
      
      *(a+2) = eta[i];
      fn_eval( a+2, &fnpar, f+2, fprime+2 );
      
      *(a+3) = eta[i] + 0.5;
      fn_eval( a+3, &fnpar, f+3, fprime+3 );
      
      *(a+4) = eta[i] + 1.0;
      fn_eval( a+4, &fnpar, f+4, fprime+4 );
      while( fprime[4] > -0.1 )
	{
	  (*(a+4))++;
	  fn_eval( a+4, &fnpar, f+4, fprime+4 );
	}
      
      // intializing ARS
      ars::initial_( &n_max, &n_start, &emax, a, f, fprime, &lb, &xlb, &ub, &xub, &ifault, iwv_p, rwv_p );
      
      if( ifault == 0 )
	{
	  // sampling
	  typedef double (Rndgen::*memfp)( void );
	  memfp uniform_gen = &Rndgen::uniform;
	  ars::sample_( iwv_p, rwv_p, &fn_eval, &fnpar, &draw, &ifault, uniform_gen, &rnd );
	  eta[i] = draw;
	}
      else
	{
	  cerr << "Error number " << ifault << endl;
	  exit(1);
	}
    }
}

int fn_eval( const double* eta, const ars::FnPar* par, double* hx, double* hpx )
{
  *hx = -exp( *eta ) + *eta * par->y - 0.5 * ( *eta - par->mean ) * ( *eta - par->mean ) / par->var;
  *hpx = -exp( *eta ) + par->y - ( *eta - par->mean ) / par->var;
  return 0;
}

void ZipModel::ystar_dstar_sample()
{
  for( int i=1; i<=n; i++ )
    {
      double mu_dstar = X[i] * delta;
      if( y[i] > 0 )
	{
	  ystar[i] = y[i];
	  dstar[i] = rnd.normal_lt( 0.0, mu_dstar, 1.0 );
	}
      else 
	{
	  double p = pnorm( mu_dstar, 0.0, 1.0, true, false );   
	  if( rnd.uniform() > (1-p) / (1-p+p*exp(-exp(eta[i]))) )
	    {
	      ystar[i] = 0;
	      dstar[i] = rnd.normal_lt( 0.0, mu_dstar, 1.0 );
	    }
	  else
	    {
	      ystar[i] = rnd.poisson( exp(eta[i]) );
	      dstar[i] = rnd.normal_rt( 0.0, mu_dstar, 1.0 );
	    }
	}
    }
}

void ZipModel::bbeta_sample( void )
{
  srsmatrix Var(k);
  rvector mean(k);
  Var = ( XX / sigma2 + Binv ).inv();
  mean = Var * ~X * eta / sigma2;
  bbeta = rnd.mnormal( mean, Var );
}

void ZipModel::delta_sample( void )
{
  srsmatrix Var(k);
  rvector mean(k);
  Var = ( XX + Dinv ).inv();
  mean = Var * ~X * dstar;
  delta = rnd.mnormal( mean, Var );
}

void ZipModel::sigma2_sample( void )
{
  rvector eps(n);
  eps = eta - X * bbeta;
  sigma2 = 1.0 / rnd.gamma( n / 2.0 + prior.e, eps * eps / 2.0 + prior.f );
}

void ZipModel::gamma_sample( void )
{
  for( int i=1; i<=k; i++ )
    {
      double bbeta2 = bbeta[i] * bbeta[i];
      double p_out = ( 1.0 - prior.p[i] ) / sqrt( prior.b(i,1) ) * exp( - 0.5 * bbeta2 / prior.b(i,1) );
      double p_in = prior.p[i] / sqrt( prior.b(i,2) ) * exp( - 0.5 * bbeta2 / prior.b(i,2) );
      if( rnd.uniform() < ( p_in / ( p_in + p_out ) ) )
	{
	  gamma[i] = 1.0;
	  Binv.set( i, i, 1.0 / prior.b(i,2) );
	}
      else
	{
	  gamma[i] = 0.0;
	  Binv.set( i, i, 1.0 / prior.b(i,1) );
	}
    }
}

void ZipModel::kappa_sample( void )
{
  for( int i=1; i<=k; i++ )
    {
      double delta2 = delta[i] * delta[i];
      double p_out = ( 1.0 - prior.q[i] ) / sqrt( prior.d(i,1) ) * exp( - 0.5 * delta2 / prior.d(i,1) );
      double p_in = prior.q[i] / sqrt( prior.d(i,2) ) * exp( - 0.5 * delta2 / prior.d(i,2) );
      if( rnd.uniform() < ( p_in / ( p_in + p_out ) ) )
	{
	  kappa[i] = 1.0;
	  Dinv.set( i, i, 1.0 / prior.d(i,2) );
	}
      else
	{
	  kappa[i] = 0.0;
	  Dinv.set( i, i, 1.0 / prior.d(i,1) );
	}
    }
}
     
ZipModel::ZipModel( const iarray& _y, const rmatrix& _X, Prior& _prior, Rndgen& _rnd ): prior( _prior ), rnd( _rnd )
{
  // dimensions
  n = _X.msize();
  k = _X.nsize();
  
  // resizing and constructing matrices
  y.resize(n); eta.resize(n); ystar.resize(n); dstar.resize(n);
  X.resize(n,k); XX.resize(k); 
  y = _y;
  X = _X;
  XX.syrk( true, 1.0, X, 0.0 );
  
  // resizing parameter vectors and setting starting values
  bbeta.resize(k);
  delta.resize(k);
  sigma2 = 0.1; 
  Binv.resize(k);
  for( int i=1; i<=k; i++ )
    Binv.set( i, i, 1.0 / prior.b(i,2) );
  Dinv.resize(k);
  for( int i=1; i<=k; i++ )
    Dinv.set( i, i, 1.0 / prior.d(i,2) );
  gamma.resize(k); gamma.set( 1.0 );
  kappa.resize(k); kappa.set( 1.0 );
}                         
