#ifndef _ZIC_HH
#define _ZIC_HH

#include <R.h>
#include <Rmath.h>
#include "cvm.h"
#include "ars.h"
#include "random_R.h"

int fn_eval( const double* eta, const ars::FnPar* par, double* hx, double* hpx );

struct Prior
{
  rmatrix b, d;
  double e, f;
  rvector p, q;
};

class ZipModel
{
public:

  ZipModel( const iarray& _y, const rmatrix& _X, Prior& _prior, Rndgen& _rnd );
  ~ZipModel( void ) {};
  void sample( const int, const int, const int, const int, rmatrix&, rmatrix&, rvector&, rmatrix&, rmatrix& );
  
private:

  // sampling methods
  void eta_sample( void );
  void ystar_dstar_sample( void );
  void bbeta_sample( void );
  void delta_sample( void );
  void sigma2_sample( void );
  void gamma_sample( void );
  void kappa_sample( void );

  // dimensions
  int n, k;    

  // data and latent variables
  iarray y, ystar;                
  rvector eta, dstar;
  rmatrix X;
  srsmatrix XX;

  // parameters
  rvector bbeta, delta;          
  double sigma2;
  srsmatrix Binv, Dinv;
  rvector gamma, kappa;

  // Prior
  Prior& prior;

  // Random number generator
  Rndgen& rnd;
};

#endif   
