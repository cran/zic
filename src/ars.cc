// ars.f -- translated with the help of f2c (version 20061008)
// http://www.netlib.org/f2c/libf2c.zip

#include "ars.h"
#include <stdio.h>

#ifdef __cplusplus				
extern "C" {
#endif

namespace ars{

int initial_( const int *ns, const int *m, double *emax, double *x, double *hx, double *hpx, int *lb, 
	      double *xlb, int *ub, double *xub, int *ifault, int *iwv, double *rwv )
{	
  /* System generated locals */	
  int i__1;	
  double d__1, d__2;	
	
  /* Builtin functions */	
  double log(double);	
	
  /* Local variables */	
  static int i__;	
  static double cu;	
  static int nn, ix, iz;	
  static double eps;	
  static int ihx;	
  static double alcu, hulb, huub;	
  static int iipt, ihpx, ilow, ihuz, ihigh, iscum;	
  extern double expon_(double *, double *);	
  static int horiz;	
  static double huzmax;

  /* This subroutine takes as input the number of starting values m */
  /*  and the starting values x(i),hx(i),hpx(i)  i=1,m */
  /* As output we have pointer ipt along with ilow and ihigh and the lower */
  /* and upper hulls defined  by z,hz,scum,cu,hulb,huub stored in working */
  /* vectors iwv and rwv */
  /* Ifault detects wrong starting points or non-concavity */

  /* DESCRIPTION OF PARAMETERS and place of storage */

  /*     lb   iwv(5) : boolean indicating if there is a lower bound to the */
  /*                    domain */
  /*     ub   iwv(6) : boolean indicating if there is an upper bound */
  /*     xlb  rwv(8) : value of the lower bound */
  /*     xub  rwv(9) : value of the upper bound */
  /*     emax rwv(3) : large value for which it is possible to compute */
  /*                   an exponential, eps=exp(-emax) is taken as a small */
  /*                   value used to test for numerical unstability */
  /*     m    iwv(4) : number of starting points */
  /*     ns   iwv(3) : maximum number of points defining the hulls */
  /*     x    rwv(ix+1)  : vector containing the abscissae of the starting */
  /*                       points */
  /*     hx   rwv(ihx+1) : vector containing the ordinates */
  /*     hpx  rwv(ihpx+1): vector containing the derivatives */
  /*     ifault      : diagnostic */
  /*     iwv,rwv     : int and double working vectors */

  /* Parameter adjustments */	
  --rwv;	
  --iwv;	
  --hpx;		
  --hx;	
  --x;

  /* Function Body */
  d__1 = -(*emax);
  eps = expon_(&d__1, emax);
  *ifault = 0;
  ilow = 1;
  ihigh = 1;
  nn = *ns + 1;
  /* at least one starting point */
  if (*m < 1) {
    *ifault = 1;
  }
  huzmax = hx[1];
  if (! (*ub)) {
    *xub = (float)0.;
  }
  if (! (*lb)) {
    *xlb = (float)0.;
  }
  hulb = (*xlb - x[1]) * hpx[1] + hx[1];
  huub = (*xub - x[1]) * hpx[1] + hx[1];
  /* if bounded on both sides */
  if (*ub && *lb) {
    huzmax = max(huub,hulb);
    horiz = abs(hpx[1]) < eps;
    if (horiz) {
      d__1 = (huub + hulb) * (float).5 - huzmax;
      cu = expon_(&d__1, emax) * (*xub - *xlb);
    } else {
      d__1 = huub - huzmax;
      d__2 = hulb - huub;
      cu = expon_(&d__1, emax) * (1 - expon_(&d__2, emax)) / hpx[1];
    }
  } 
  else if (*ub && ! (*lb)) {
    /* if bounded on the right and unbounded on the left */
    huzmax = huub;
    cu = (float)1. / hpx[1];
  } 
  else if (! (*ub) && *lb) {
    /* if bounded on the left and unbounded on the right */
    huzmax = hulb;
    cu = (float)-1. / hpx[1];
    /* if unbounded at least 2 starting points */
  } 
  else {
    cu = (float)0.;
    if (*m < 2) {
      *ifault = 1;
    }
  }
  if (cu > (float)0.) {
    alcu = log(cu);
  }
  /* set pointers */
  iipt = 6;
  iz = 9;
  ihuz = nn + iz;
  iscum = nn + ihuz;
  ix = nn + iscum;
  ihx = nn + ix;
  ihpx = nn + ihx;
  /* store values in working vectors */
  iwv[1] = ilow;
  iwv[2] = ihigh;
  iwv[3] = *ns;
  iwv[4] = 1;
  if (*lb) {
    iwv[5] = 1;
  } 
  else {
    iwv[5] = 0;
  }
  if (*ub) {
    iwv[6] = 1;
  } 
  else {
    iwv[6] = 0;
  }
  if (*ns < *m) {
    *ifault = 2;
  }
  iwv[iipt + 1] = 0;
  rwv[1] = hulb;
  rwv[2] = huub;
  rwv[3] = *emax;
  rwv[4] = eps;
  rwv[5] = cu;
  rwv[6] = alcu;
  rwv[7] = huzmax;
  rwv[8] = *xlb;
  rwv[9] = *xub;
  rwv[iscum + 1] = (float)1.;
  i__1 = *m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    rwv[ix + i__] = x[i__];
    rwv[ihx + i__] = hx[i__];
    rwv[ihpx + i__] = hpx[i__];
    /* L9: */
  }
  /* create lower and upper hulls */
  i__ = 1;
 L10:
  if (i__ < *m) {
    update_(&iwv[4], &iwv[1], &iwv[2], &iwv[iipt + 1], &rwv[iscum + 1], &
	    rwv[5], &rwv[ix + 1], &rwv[ihx + 1], &rwv[ihpx + 1], &rwv[iz 
								      + 1], &rwv[ihuz + 1], &rwv[7], &rwv[3], lb, &rwv[8], &rwv[1], 
	    ub, &rwv[9], &rwv[2], ifault, &rwv[4], &rwv[6]);
    i__ = iwv[4];
    if (*ifault != 0) {
      return 0;
    }
    goto L10;
  }
  /* test for wrong starting points */
  if (! (*lb) && hpx[iwv[1]] < eps) {
    *ifault = 3;
  }
  if (! (*ub) && hpx[iwv[2]] > -eps) {
    *ifault = 4;
  }
  return 0;
} /* initial_ */


  int sample_( int *iwv, double *rwv, int (*eval)(const double*, const FnPar*, double*, double*), const FnPar* fnpar, double *beta, int *ifault, memfp ugen, Rndgen* rnd )
{
  static int lb, ub;
  static int nn, ns, ix, iz, ihx;
  static int iipt, ihpx, ihuz, iscum;

  /* Parameter adjustments */
  --rwv;
  --iwv;

  /* Function Body */
  iipt = 6;
  iz = 9;
  ns = iwv[3];
  nn = ns + 1;
  ihuz = nn + iz;
  iscum = nn + ihuz;
  ix = nn + iscum;
  ihx = nn + ix;
  ihpx = nn + ihx;
  lb = FALSE_;
  ub = FALSE_;
  if (iwv[5] == 1) {
    lb = TRUE_;
  }
  if (iwv[6] == 1) {
    ub = TRUE_;
  }

  /* call sampling subroutine */

  spl1_(&ns, &iwv[4], &iwv[1], &iwv[2], &iwv[iipt + 1], &rwv[iscum + 1], &rwv[5], &rwv[ix + 1], &rwv[ihx + 1], &rwv[ihpx + 1], &rwv[iz + 1],
	&rwv[ihuz + 1], &rwv[7], &lb, &rwv[8], &rwv[1], &ub, &rwv[9], &	rwv[2], eval, fnpar, beta, ifault, &rwv[3], &rwv[4], &rwv[6], ugen, rnd );

  return 0;
}


int spl1_(int *ns, int *n, int *ilow, int * ihigh, int *ipt, double *scum, double *cu, double *x, double *hx, double *hpx, double *z__, double *huz, 
          double *huzmax, int *lb, double *xlb, double *hulb, int *ub, double *xub, double *huub, 
          int (*eval)(const double*, const FnPar*, double*, double*), const FnPar* fnpar, double *beta, int *ifault, double *emax, double *eps, double *alcu, memfp ugen, Rndgen* rnd )
{
  /* Builtin functions */
  double log(double);
  
  /* Local variables */
  static int i__, j, n1;
  static double u1, u2, fx, alu1, alhl, alhu;
  static int sampld;
  
  /* this subroutine performs the adaptive rejection sampling, it calls */
  /* subroutine splhull to sample from the upper hull ,if the sampling */
  /* involves a function evaluation it calls the updating subroutine */
  /* ifault is a diagnostic of any problem: non concavity, 0 random number */
  /*    or numerical imprecision */

  /* Parameter adjustments */
  --huz;
  --z__;
  --hpx;
  --hx;
  --x;
  --scum;
  --ipt;

  /* Function Body */
  sampld = FALSE_;
 L10:
  if (! sampld) 
    {
      u2 = (rnd->*ugen)();
      /* test for zero random number */
      if (u2 == (float)0.) 
	{
	  *ifault = 6;
	  return 0;
	}
      splhull_(&u2, &ipt[1], ilow, lb, xlb, hulb, huzmax, alcu, &x[1], &hx[1], 
	     &hpx[1], &z__[1], &huz[1], &scum[1], eps, emax, beta, &i__, &j);
      /* sample u1 to compute rejection */
      u1 = (rnd->*ugen)();
      if (u1 == (float)0.) 
	{
	  *ifault = 6;
	}
      alu1 = log(u1);
      /* compute alhu: upper hull at point u1 */
      alhu = hpx[i__] * (*beta - x[i__]) + hx[i__] - *huzmax;
      if (*beta > x[*ilow] && *beta < x[*ihigh]) {
	/* compute alhl: value of the lower hull at point u1 */
      if (*beta > x[i__]) {
	j = i__;
	i__ = ipt[i__];
      }
      alhl = hx[i__] + (*beta - x[i__]) * (hx[i__] - hx[j]) / (x[i__] - x[j]) - *huzmax;
      /* squeezing test */
      if (alhl - alhu > alu1) {
	sampld = TRUE_;
      }
    }
    /* if not sampled evaluate the function ,do the rejection test and update */
    if (! sampld) {
      n1 = *n + 1;
      x[n1] = *beta;
      (*eval)(&x[n1], fnpar, &hx[n1], &hpx[n1]);
      fx = hx[n1] - *huzmax;
      if (alu1 < fx - alhu) {
	sampld = TRUE_;
      }
      /* update while the number of points defining the hulls is lower than ns */
      if (*n < *ns) {
	update_(n, ilow, ihigh, &ipt[1], &scum[1], cu, &x[1], &hx[1], 
		&hpx[1], &z__[1], &huz[1], huzmax, emax, lb, xlb, 
		hulb, ub, xub, huub, ifault, eps, alcu);
      }
      if (*ifault != 0) {
	return 0;
      }
    }
    goto L10;
  }
  return 0;
} /* spl1_ */


/* *********************************************************************** */

/* Subroutine */ int splhull_(double *u2, int *ipt, int *ilow, 
			      int *lb, double *xlb, double *hulb, double *huzmax, 
			      double *alcu, double *x, double *hx, double *hpx, 
			      double *z__, double *huz, double *scum, double *eps, 
			      double *emax, double *beta, int *i__, int *j)
{
  /* System generated locals */
  double d__1, d__2;

  /* Builtin functions */
  double log(double);

  /* Local variables */
  static double eh, sign, logdu, logtg;
  extern double expon_(double *, double *);
  static int horiz;


  /* this subroutine samples beta from the normalised upper hull */


  /* Parameter adjustments */
  --scum;
  --huz;
  --z__;
  --hpx;
  --hx;
  --x;
  --ipt;

  /* Function Body */
  *i__ = *ilow;

  /* find from which exponential piece you sample */
 L20:
  if (*u2 > scum[*i__]) {
    *j = *i__;
    *i__ = ipt[*i__];
    goto L20;
  }
  if (*i__ == *ilow) {
    /* sample below z(ilow),depending on the existence of a lower bound */
    if (*lb) {
      eh = *hulb - *huzmax - *alcu;
      horiz = (d__1 = hpx[*ilow], abs(d__1)) < *eps;
      if (horiz) {
	d__1 = -eh;
	*beta = *xlb + *u2 * expon_(&d__1, emax);
      } else {
	sign = (d__1 = hpx[*i__], abs(d__1)) / hpx[*i__];
	d__2 = (d__1 = hpx[*i__], abs(d__1));
	logtg = log(d__2);
	logdu = log(*u2);
	eh = logdu + logtg - eh;
	if (eh < *emax) {
	  d__1 = sign * expon_(&eh, emax) + (float)1.;
	  *beta = *xlb + log(d__1) / hpx[*i__];
	} else {
	  *beta = *xlb + eh / hpx[*i__];
	}
      }
    } else {
      /*     hpx(i) must be positive , x(ilow) is left of the mode */
      d__1 = hpx[*i__] * *u2;
      *beta = (log(d__1) + *alcu - hx[*i__] + x[*i__] * hpx[*i__] + *
	       huzmax) / hpx[*i__];
    }
  } else {
    /*   sample above(j) */
    eh = huz[*j] - *huzmax - *alcu;
    horiz = (d__1 = hpx[*i__], abs(d__1)) < *eps;
    if (horiz) {
      d__1 = -eh;
      *beta = z__[*j] + (*u2 - scum[*j]) * expon_(&d__1, emax);
    } else {
      sign = (d__1 = hpx[*i__], abs(d__1)) / hpx[*i__];
      d__2 = (d__1 = hpx[*i__], abs(d__1));
      logtg = log(d__2);
      d__1 = *u2 - scum[*j];
      logdu = log(d__1);
      eh = logdu + logtg - eh;
      if (eh < *emax) {
	d__1 = sign * expon_(&eh, emax) + (float)1.;
	*beta = z__[*j] + log(d__1) / hpx[*i__];
      } else {
	*beta = z__[*j] + eh / hpx[*i__];
      }
    }
  }
  return 0;
} /* splhull_ */


/* *********************************************************************** */
/* Subroutine */ int intersection_(double *x1, double *y1, double 
				   *yp1, double *x2, double *y2, double *yp2, double *z1,
				   double *hz1, double *eps, int *ifault)
{
  static double dh, y12, y21;


  /* computes the intersection (z1,hz1) between 2 tangents defined by */
  /*   x1,y1,yp1 and x2,y2,yp2 */


  /* first test for non-concavity */
  y12 = *y1 + *yp1 * (*x2 - *x1);
  y21 = *y2 + *yp2 * (*x1 - *x2);
  if (y21 < *y1 || y12 < *y2) {
    *ifault = 5;
    return 0;
  }
  dh = *yp2 - *yp1;
  /*  IF the lines are nearly parallel, */
  /*  the intersection is taken at the midpoint */
  if (abs(dh) <= *eps) {
    *z1 = (*x1 + *x2) * (float).5;
    *hz1 = (*y1 + *y2) * (float).5;
    /*  Else compute from the left or the right for greater numerical */
    /*       precision */
  } else if (abs(*yp1) < abs(*yp2)) {
    *z1 = *x2 + (*y1 - *y2 + *yp1 * (*x2 - *x1)) / dh;
    *hz1 = *yp1 * (*z1 - *x1) + *y1;
  } else {
    *z1 = *x1 + (*y1 - *y2 + *yp2 * (*x2 - *x1)) / dh;
    *hz1 = *yp2 * (*z1 - *x2) + *y2;
  }
  /*  test for misbehaviour due to numerical imprecision */
  if (*z1 < *x1 || *z1 > *x2) {
    *ifault = 7;
  }
  return 0;
} /* intersection_ */


/* *********************************************************************** */

/* Subroutine */ int update_(int *n, int *ilow, int *ihigh, 
			     int *ipt, double *scum, double *cu, double *x, 
			     double *hx, double *hpx, double *z__, double *huz, 
			     double *huzmax, double *emax, int *lb, double *xlb, 
			     double *hulb, int *ub, double *xub, double *huub, 
			     int *ifault, double *eps, double *alcu)
{
  /* System generated locals */
  double d__1, d__2;

  /* Builtin functions */
  double log(double);

  /* Local variables */
  static int i__, j;
  static double u, dh;
  extern /* Subroutine */ int intersection_(double *, double *, 
					    double *, double *, double *, double *, 
					    double *, double *, double *, int *);
  extern double expon_(double *, double *);
  static int horiz;


  /* this subroutine increments n and updates all the parameters which */
  /* define the lower and the upper hull */


  /* DESCRIPTION OF PARAMETERS and place of storage */

  /*     ilow iwv(1)    : index of the smallest x(i) */
  /*     ihigh iwv(2)   : index of the largest x(i) */
  /*     n    iwv(4)    : number of points defining the hulls */
  /*     ipt  iwv(iipt) : pointer array:  ipt(i) is the index of the x(.) */
  /*                      immediately larger than x(i) */
  /*     hulb rwv(1)    : value of the upper hull at xlb */
  /*     huub rwv(2)    : value of the upper hull at xub */
  /*     cu   rwv(5)    : integral of the exponentiated upper hull divided */
  /*                      by exp(huzmax) */
  /*     alcu rwv(6)    : logarithm of cu */
  /*     huzmax rwv(7)  : maximum of huz(i); i=1,n */
  /*     z    rwv(iz+1) : z(i) is the abscissa of the intersection between */
  /*                      the tangents at x(i) and x(ipt(i)) */
  /*     huz  rwv(ihuz+1): huz(i) is the ordinate of the intersection */
  /*                        defined above */
  /*     scum rwv(iscum): scum(i) is the cumulative probability of the */
  /*                      normalised exponential of the upper hull */
  /*                      calculated at z(i) */
  /*     eps  rwv(4)    : =exp(-emax) a very small number */

  /* Parameter adjustments */
  --huz;
  --z__;
  --hpx;
  --hx;
  --x;
  --scum;
  --ipt;

  /* Function Body */
  ++(*n);
  /* update z,huz and ipt */
  if (x[*n] < x[*ilow]) {
    /* insert x(n) below x(ilow) */
    /*   test for non-concavity */
    if (hpx[*ilow] > hpx[*n]) {
      *ifault = 5;
    }
    ipt[*n] = *ilow;
    intersection_(&x[*n], &hx[*n], &hpx[*n], &x[*ilow], &hx[*ilow], &hpx[*ilow], &z__[*n], &huz[*n], eps, ifault);
    if (*ifault != 0) {
      return 0;
    }
    if (*lb) {
      *hulb = hpx[*n] * (*xlb - x[*n]) + hx[*n];
    }
    *ilow = *n;
  } else {
    i__ = *ilow;
    j = i__;
    /* find where to insert x(n) */
  L10:
    if (x[*n] >= x[i__] && ipt[i__] != 0) {
      j = i__;
      i__ = ipt[i__];
      goto L10;
    }
    if (x[*n] >= x[i__]) {
      /* insert above x(ihigh) */
      /*   test for non-concavity */
      if (hpx[i__] < hpx[*n]) {
	*ifault = 5;
      }
      *ihigh = *n;
      ipt[i__] = *n;
      ipt[*n] = 0;
      intersection_(&x[i__], &hx[i__], &hpx[i__], &x[*n], &hx[*n], &hpx[*n], &z__[i__], &huz[i__], eps, ifault);
      if (*ifault != 0) {
	return 0;
      }
      *huub = hpx[*n] * (*xub - x[*n]) + hx[*n];
      z__[*n] = (float)0.;
      huz[*n] = (float)0.;
    } else {
      /* insert x(n) between x(j) and x(i) */
      /*   test for non-concavity */
      if (hpx[j] < hpx[*n] || hpx[i__] > hpx[*n]) {
	*ifault = 5;
      }
      ipt[j] = *n;
      ipt[*n] = i__;
      /*     insert z(j) between x(j) and x(n) */
      intersection_(&x[j], &hx[j], &hpx[j], &x[*n], &hx[*n], &hpx[*n], &z__[j], &huz[j], eps, ifault);
      if (*ifault != 0) {
	return 0;
      }
      /*     insert z(n) between x(n) and x(i) */
      intersection_(&x[*n], &hx[*n], &hpx[*n], &x[i__], &hx[i__], &hpx[i__], &z__[*n], &huz[*n], eps, ifault);
      if (*ifault != 0) {
	return 0;
      }
    }
  }
  /* update huzmax */
  j = *ilow;
  i__ = ipt[j];
  *huzmax = huz[j];
 L20:
  if (huz[j] < huz[i__] && ipt[i__] != 0) {
    j = i__;
    i__ = ipt[i__];
    /* Computing MAX */
    d__1 = *huzmax, d__2 = huz[j];
    *huzmax = max(d__1,d__2);
    goto L20;
  }
  if (*lb) {
    *huzmax = max(*huzmax,*hulb);
  }
  if (*ub) {
    *huzmax = max(*huzmax,*huub);
  }
  /* update cu */
  /*  scum receives area below exponentiated upper hull left of z(i) */
  i__ = *ilow;
  horiz = (d__1 = hpx[*ilow], abs(d__1)) < *eps;
  if (! (*lb) && ! horiz) {
    d__1 = huz[i__] - *huzmax;
    *cu = expon_(&d__1, emax) / hpx[i__];
  } else if (*lb && horiz) {
    d__1 = *hulb - *huzmax;
    *cu = (z__[*ilow] - *xlb) * expon_(&d__1, emax);
  } else if (*lb && ! horiz) {
    dh = *hulb - huz[i__];
    if (dh > *emax) {
      d__1 = *hulb - *huzmax;
      *cu = -expon_(&d__1, emax) / hpx[i__];
    } else {
      d__1 = huz[i__] - *huzmax;
      *cu = expon_(&d__1, emax) * (1 - expon_(&dh, emax)) / hpx[i__];
    }
  } else {
    *cu = 0.;
  }
  scum[i__] = *cu;
  j = i__;
  i__ = ipt[i__];
 L30:
  if (ipt[i__] != 0) {
    dh = huz[j] - huz[i__];
    horiz = (d__1 = hpx[i__], abs(d__1)) < *eps;
    if (horiz) {
      d__1 = (huz[i__] + huz[j]) * (float).5 - *huzmax;
      *cu += (z__[i__] - z__[j]) * expon_(&d__1, emax);
    } else {
      if (dh < *emax) {
	d__1 = huz[i__] - *huzmax;
	*cu += expon_(&d__1, emax) * (1 - expon_(&dh, emax)) / hpx[i__];
      } else {
	d__1 = huz[j] - *huzmax;
	*cu -= expon_(&d__1, emax) / hpx[i__];
      }
    }
    j = i__;
    i__ = ipt[i__];
    scum[j] = *cu;
    goto L30;
  }
  horiz = (d__1 = hpx[i__], abs(d__1)) < *eps;
  /* if the derivative is very small the tangent is nearly horizontal */
  if (! (*ub || horiz)) {
    d__1 = huz[j] - *huzmax;
    *cu -= expon_(&d__1, emax) / hpx[i__];
  } else if (*ub && horiz) {
    d__1 = (*huub + hx[i__]) * (float).5 - *huzmax;
    *cu += (*xub - x[i__]) * expon_(&d__1, emax);
  } else if (*ub && ! horiz) {
    dh = huz[j] - *huub;
    if (dh > *emax) {
      d__1 = huz[j] - *huzmax;
      *cu -= expon_(&d__1, emax) / hpx[i__];
    } else {
      d__1 = *huub - *huzmax;
      *cu += expon_(&d__1, emax) * (1 - expon_(&dh, emax)) / hpx[i__];
    }
  }
  scum[i__] = *cu;
  if (*cu > 0.) {
    *alcu = log(*cu);
  }
  /* normalize scum to obtain a cumulative probability while excluding */
  /*    unnecessary points */
  i__ = *ilow;
  u = (*cu - scum[i__]) / *cu;
  if (u == (float)1. && hpx[ipt[i__]] > 0.) {
    *ilow = ipt[i__];
    scum[i__] = (float)0.;
  } else {
    scum[i__] = (float)1. - u;
  }
  j = i__;
  i__ = ipt[i__];
 L40:
  if (ipt[i__] != 0) {
    j = i__;
    i__ = ipt[i__];
    u = (*cu - scum[j]) / *cu;
    if (u == (float)1. && hpx[i__] > 0.) {
      *ilow = i__;
    } else {
      scum[j] = (float)1. - u;
    }
    goto L40;
  }
  scum[i__] = (float)1.;
  if (*ub) {
    *huub = hpx[*ihigh] * (*xub - x[*ihigh]) + hx[*ihigh];
  }
  if (*lb) {
    *hulb = hpx[*ilow] * (*xlb - x[*ilow]) + hx[*ilow];
  }
  return 0;
} /* update_ */


/* *********************************************************************** */

double expon_(double *x, double *emax)
{
  /* System generated locals */
  double ret_val;

  /* Builtin functions */
  double exp(double);

  /* performs an exponential without underflow */
  if (*x < -(*emax)) {
    ret_val = (float)0.;
  } else {
    ret_val = exp(*x);
  }
  return ret_val;
} 

}

#ifdef __cplusplus				
}
#endif					       
