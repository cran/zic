#ifndef _ARS_H_
#define _ARS_H_

#include "random_R.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef /* Subroutine */ int (*S_fp)(...);
typedef double (Rndgen::*memfp)( void );

#define max(a,b) ((a) >= (b) ? (a) : (b))
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define TRUE_ (1)
#define FALSE_ (0)

namespace ars
{
  struct FnPar
  {
    int y;
    double mean, var;
  };

  int initial_(const int* ns, const int* m, double *emax, double *x, double *hx, double *hpx, int *lb, 
	       double *xlb, int *ub, double *xub, int *ifault, int *iwv, double *rwv);

  int sample_( int *iwv, double *rwv, int (*eval)(const double*, const FnPar*, double*, double*), const FnPar*, double *beta, int *ifault, memfp, Rndgen* );

  int spl1_(int *ns, int *n, int *ilow, int * ihigh, int *ipt, double *scum, double *cu, double *x, double *hx, double *hpx, double *z__, double *huz, 
	    double *huzmax, int *lb, double *xlb, double *hulb, int *ub, double *xub, double *huub, int (*eval)(const double*, const FnPar*, double*, double*), const FnPar*, double *beta, int *ifault, double *emax, double *eps, double *alcu, memfp, Rndgen* );

  int splhull_(double *u2, int *ipt, int *ilow, int *lb, double *xlb, double *hulb, double *huzmax, 
	       double *alcu, double *x, double *hx, double *hpx, 
	       double *z__, double *huz, double *scum, double *eps, 
	       double *emax, double *beta, int *i__, int *j);
  
  int intersection_(double *x1, double *y1, double *yp1, double *x2, double *y2, double *yp2, double *z1,
		    double *hz1, double *eps, int *ifault);
  
  int update_(int *n, int *ilow, int *ihigh, 
	      int *ipt, double *scum, double *cu, double *x, 
	      double *hx, double *hpx, double *z__, double *huz, 
	      double *huzmax, double *emax, int *lb, double *xlb, 
	      double *hulb, int *ub, double *xub, double *huub, 
	      int *ifault, double *eps, double *alcu);
  
  double expon_(double *x, double *emax);
}

}

#endif
