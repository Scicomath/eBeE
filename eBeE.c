#include <stdio.h>
#include "udStruct.h"
#include "eBFun.h"


int main(void)
{
  struct intargu ag;
  double x, y, tau, R, b, Y0, a, Z, eBy;
  ag.nvec = 1;
  ag.epsrel = 0.0001;
  ag.epsabs = 0.001;
  ag.flags = 0 | 8; // 0 | 8
  ag.seed = 0;
  ag.mineval = 0;
  ag.maxeval = 100000;

  ag.nstart = 1000;
  ag.nincrease = 500;
  ag.nbatch = 1000;
  ag.gridno = 0;

  ag.statefile = NULL;
  ag.spin = NULL;

  x = 0.0;
  y = 0.0;
  tau = 0.1;
  R = 7;
  b = 4;
  Y0 = 5.4;
  a = 0.5;
  Z = 79.0;
  eB(x, y, tau, R, b, Y0, a, Z, &ag, &eBy);

  return 0;
  
}
