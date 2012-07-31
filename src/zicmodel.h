#ifndef _ZICMODEL_H
#define _ZICMODEL_H

#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include "random_R.h"
#include "spikeslab.h"

using namespace arma;



class ZicModel
{
public:
  ZicModel( const ivec& _y, const mat& _X, const SpikeSlabPrior& _betaprior, const SpikeSlabPrior& _deltaprior, const double _tune );
  void sample( const int nburnin, const int nmcmc, const int nthin, 
               vec& alphamem, mat& betamem, vec& gammamem, mat& deltamem, vec& sigma2mem, imat& Ibetamem, imat& Ideltamem, 
               vec& omegabetamem, vec& omegadeltamem, vec& accmem );
  void sigma2_update();
  void latvar_update();
 
private:
  const ivec y;
  const mat X; 
  const int n, n0, k;
  vec eta, dstar;
 
  SpikeSlabParameter beta, delta;
  double sigma2;

  const double tune;
  double acc;
  const SpikeSlabPrior betaprior, deltaprior;

  Rnd rnd;
};



#endif
