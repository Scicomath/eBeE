#include <stdio.h>
#include "udStruct.h"
#include "eBFun.h"


int main(void)
{
  struct intargu ag;
  double x, y, tau, R, b, Y0, a, Z, eBy;
  ag.nvec = 1;
  ag.epsrel = 1e-6;
  ag.epsabs = 1e-6;
  ag.flags = 0 | 8; // 0 | 8
  ag.seed = 0;
  ag.smineval = 1.5e5;
  ag.smaxeval = 1e6;
  ag.pmineval = 2.5e5;
  ag.pmaxeval = 1e7;

  ag.nstart = 1000;
  ag.nincrease = 500;
  ag.nbatch = 1000;
  ag.gridno = 0;

  ag.statefile = NULL;
  ag.spin = NULL;

  x = 1.0;
  y = 1.0;
  tau = 0.01;
  R = 7.0;
  b = 4;
  Y0 = 5.36;
  a = 0.5;
  Z = 79.0;
  eB(x, y, tau, R, b, Y0, a, Z, &ag, &eBy);

  printf("eBy = %5.3f\n", eBy);
  return 0;
  
}
