#include <stdio.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "cuba.h"
#include "udStruct.h"
#include "rhoFun.h"
#include "sq.h"

static int rho_Integrand(const int *ndim, const double xx[],
			 const int *ncomp, double ff[], void *userdata);

int main(void)
{
  struct intargu ag;
  ag.nvec = 1;
  ag.epsrel = 1e-6;
  ag.epsabs = 1e-6;
  ag.flags = 0 | 8; // 0 | 8
  ag.seed = 0;
  ag.pmineval = 2e4; // 2.5e5
  ag.pmaxeval = 5e4;   // 1e7

  ag.nstart = 1000;
  ag.nincrease = 500;
  ag.nbatch = 1000;
  ag.gridno = 0;

  ag.statefile = NULL;
  ag.spin = NULL;

  static struct userdata ud;
  ud.R = 6.38;
  ud.b = 8;
  ud.d = 0.535;
  ud.n0 = 8.596268e-4;
  ud.Y0 = 5.36;
  ud.flag = '-';
  
  int neval, fail;
  double integral, error, prob;

  Vegas(3, 1, rho_Integrand, &ud, ag.nvec,
	ag.epsrel, ag.epsabs, ag.flags, ag.seed,
	ag.pmineval, ag.pmaxeval, ag.nstart, ag.nincrease,
	ag.nbatch, ag.gridno, ag.statefile, ag.spin,
	&neval, &fail, &integral, &error, &prob);
  printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
  printf("    int_rho:\t%.8f +- %.8f\tp = %.3f\n", integral, error, prob);

  

  return 0;
}

static int rho_Integrand(const int *ndim, const double xx[],
			 const int *ncomp, double ff[], void *userdata) {
  struct userdata *ud = (struct userdata *) userdata;
  static double jacobian;
  
  static double x_p;
  static double y_p;
  static double z_p;
  static double Imin[3];
  static double Imax[3];
  static double extra;
  static double gamma;
  static double b;

  extra = 3;
  gamma = cosh(ud->Y0);

  if (ud->flag == '+') {
    b = ud->b;
  } else {
    b = -ud->b;
  }

  Imin[0] = -(ud->R + b/2.0) - extra;
  Imax[0] = ud->R - b/2.0 + extra;
  Imin[1] = -ud->R - extra;
  Imax[1] = ud->R + extra;
  Imin[2] = -(ud->R + extra) / gamma;
  Imax[2] = (ud->R + extra) / gamma;

  x_p = Imin[0] + (Imax[0] - Imin[0]) * xx[0];
  y_p = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
  z_p = Imin[2] + (Imax[2] - Imin[2]) * xx[2];

  ff[0] = rhoFun(x_p, y_p, z_p, ud->R, ud->b, ud->d, ud->n0, ud->Y0, ud->flag);

  jacobian = 1.0;
  for (int i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }
  
  ff[0] = jacobian * ff[0];

  return 0;
}


