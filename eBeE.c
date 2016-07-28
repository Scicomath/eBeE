#include <stdio.h>
#include "udStruct.h"
#include "eBFun.h"


int main(void)
{
  struct intargu ag;
  double x, y, t, R, b, Y0, a, Z, eBy, error;
  int verbose;
  ag.nvec = 1;
  ag.epsrel = 1e-6;
  ag.epsabs = 1e-6;
  ag.flags = 0 | 8; // 0 | 8
  ag.seed = 0;
  ag.smineval = 1.5e5; // 1.5e5
  ag.smaxeval = 4e6;   // 1e6
  ag.pmineval = 2e5; // 2.5e5
  ag.pmaxeval = 5e7;   // 1e7

  ag.nstart = 1000;
  ag.nincrease = 500;
  ag.nbatch = 1000;
  ag.gridno = 0;

  ag.statefile = NULL;
  ag.spin = NULL;

  x = 0.0;
  y = 0.0;
  t = 0.1;
  R = 6.38;
  b = 6.0;
  Y0 = 5.36;
  a = 0.5;
  Z = 79.0;
  printf("# parameter:\n");
  //  printf("# t = %g\n", t);
  printf("# R = %g\n", R);
  printf("# b = %g\n", b);
  printf("# Y0 = %g\n", Y0);
  printf("# Z = %g\n", Z);

  
  // 计算原点磁场随时间的变化
  int i, N;
  double tmin, tmax;
  verbose = 0;
  x = 0.0;
  y = 0.0;
  N = 100;
  tmin = 0.0;
  tmax = 3.0;
  printf("# t       \teBy     \tabserr  \trelerr\n");
  for (i = 0; i <= N; i++) {
    t = tmin + (tmax - tmin)*i/N;
    eB(x, y, t, R, b, Y0, a, Z, &ag, &eBy, &error, verbose);
    printf("  %-8g\t%-8g\t%-8g\t%-8g\n", t, eBy, error, error/eBy*100.0);
  }
  

  /* // 计算x-y平面
  int i, j;
  verbose = 0;
  for (i = 0; i <= 400; i++) { 
    for (j = 0; j <= 400; j++) {
      x = -10. + (20.)*i/400.;
      y = -10. + (20.)*j/400.;
      eB(x, y, t, R, b, Y0, a, Z, &ag, &eBy, &error, verbose);
      printf("%f %f %f\n",x, y, eBy);
    }
  }
  */

  /*
  // 计算单点磁场
  verbose = 1;
  x = 0.0;
  y = 0.0;
  t = 0.51;
  eB(x, y, t, R, b, Y0, a, Z, &ag, &eBy, &error, verbose);
  printf("eBy = %g\terror = %g\trelerror = %g%%\n", eBy, error, error/eBy*100.0);
  */
  
  return 0;
}
