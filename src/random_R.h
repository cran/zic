#ifndef _RANDOM_R_HH
#define _RANDOM_R_HH

#include "cvm.h"
#include <Rmath.h>
#include <R_ext/Random.h>

using namespace cvm;

class Rndgen
{
public:
  Rndgen( void ) { GetRNGstate(); };
  ~Rndgen() { PutRNGstate(); };

  /////
  
  double uniform( void ) { return unif_rand(); };
  double normal( void ) { return norm_rand(); };
  double gamma( const double a, const double b ) { return rgamma( a, 1.0 / b ); };
  int poisson( const double a ) { return rpois( a ); };

  /////

  void normal( rvector& x )
  {
    int i;
    for( i=1; i<=x.size(); ++i ) 
      x[i] = normal();
  };

  rvector mnormal( const rvector& mu, const srsmatrix& Sigma ) 
  {
    int k = mu.size();
    srmatrix uchol(k);
    rvector draw(k);
    uchol.cholesky( Sigma ); 
    normal( draw );
    return mu + draw * uchol;
  };

  double normal_lt( const double a, const double var )
  {
    // as in the GSL scientific library; see also Knuth, v2, 3rd ed, pp 123-128

    double s = a / sqrt( var );
  
    if( s < 1.0 )
      {
	double x;
	do
	  {
	    x = normal();
	  } 
	while( x < s );
	return x * sqrt( var );
      }

    else
      {
	double u, v, x;
	do
	  {
	    u = uniform();
	    do
	      {
		v = uniform();
	      }
	    while( v == 0.0 );
	    x = sqrt (s * s - 2 * log (v));
	  }
	while( x * u > s );
	return x * sqrt( var );
      }
  };

  double normal_lt( const double a, const double mean, const double var )
  {
    return mean + normal_lt( a - mean, var );
  };

  double normal_rt( const double a, const double mean, const double var )
  {
    return -normal_lt( -a, -mean, var );
  };
};

#endif

