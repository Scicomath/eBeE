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
  ag.smineval = 1.5e4; // 1.5e5
  ag.smaxeval = 4e4;   // 1e6
  ag.pmineval = 2e4; // 2.5e5
  ag.pmaxeval = 5e4;   // 1e7

  ag.nstart = 1000;
  ag.nincrease = 500;
  ag.nbatch = 1000;
  ag.gridno = 0;

  ag.statefile = NULL;
  ag.spin = NULL;

  x = 4.0;
  y = 2.0;
  tau = 0.01;
  R = 7.0;
  b = 10.0;
  Y0 = 5.36;
  a = 0.5;
  Z = 79.0;
  
  int i, j;
  for (i = 0; i <= 200; i++) { // i = 201; i <= 400; i++
    for (j = 0; j <= 400; j++) {
      x = -10. + (20.)*i/400.;
      y = -10. + (20.)*j/400.;
      eB(x, y, tau, R, b, Y0, a, Z, &ag, &eBy);
      printf("%f %f %f\n",x, y, eBy);
    }
  }

  
  /* eB(x, y, tau, R, b, Y0, a, Z, &ag, &eBy); */

  /* printf("eBy = %5.3f\n", eBy); */
  return 0;
  
}
