#include <stdio.h>
#include <math.h>
#include "udStruct.h"
#include "eBFun.h"
#include "sqrtStoY.h"
#include "mean_eB.h"


int main(void)
{
  struct intargu ag;
  double eBy, mean_eBy, error;
  int verbose;
  ag.nvec = 1;
  ag.epsrel = 1e-3;
  ag.epsabs = 1e-9;
  ag.flags = 0 | 8; // 0 | 8
  ag.seed = 0;
  ag.smineval = 1.5e3; // 1.5e5
  ag.smaxeval = 4e6;   // 1e6
  ag.pmineval = 2e3; // 2.5e5
  ag.pmaxeval = 1e7;   // 1e7

  ag.nstart = 1000;
  ag.nincrease = 500;
  ag.nbatch = 1000;
  ag.gridno = 0;

  ag.statefile = NULL;
  ag.spin = NULL;

  struct userdata ud;
  ud.x = 0.0;
  ud.y = 0.0;
  ud.t = 0.1;
  ud.R = 6.38;
  ud.b = 4.0;
  ud.Y0 = sqrtStoY(200);//质心系能量200GeV
  ud.d = 0.535;
  ud.n0 = 8.596268e-4;
  ud.a = 0.5;
  ud.Z = 79.0;
  printf("# parameter:\n");
  printf("# R = %g\n", ud.R);
  printf("# b = %g\n", ud.b);
  printf("# Y0 = %g\n", ud.Y0);
  printf("# Z = %g\n", ud.Z);

  ud.method = 0; // 0|1|2 : 分别表示Kharzeev,莫玉俊,艾鑫的方法
  /*
  // 计算原点磁场随时间的变化
  int i, N;
  double tmin, tmax;
  verbose = 0;
  ud.x = 0.0;
  ud.y = 0.0;
  N = 200;
  tmin = -3.0;
  tmax = 3.0;
  printf("# t       \teBy     \tabserr  \trelerr\n");
  for (i = 0; i <= N; i++) {
    ud.t = tmin + (tmax - tmin)*i/N;
    if (ud.t >= 0) {
      eB(&ud, &ag, &eBy, &error, verbose);
    } else {
      eBtminus(&ud, &ag, &eBy, &error, verbose);
    }
    printf("  %-8g\t%-8g\t%-8g\t%-8g\n", ud.t, eBy, error, fabs(error/eBy*100.0));
  }
  */

  /*
  // 计算x-y平面
  int i, j;
  verbose = 0;
  for (i = 0; i <= 400; i++) { 
    for (j = 0; j <= 400; j++) {
      ud.x = -10. + (20.)*i/400.;
      ud.y = -10. + (20.)*j/400.;
      if (ud.t >= 0) {
	eB(&ud, &ag, &eBy, &error, verbose);
      } else {
	eBtminus(&ud, &ag, &eBy, &error, verbose);
      }
      printf("  %-8g\t%-8g\t%-8g\t%-8g\t%-8g\n", ud.x, ud.y, eBy, error, fabs(error/eBy*100.0));
    }
  }
  */

  /*
  // 计算单点磁场
  verbose = 1;
  ud.x = -10.0;
  ud.y = -2.95;
  ud.t = 0.1;
  if (ud.t >= 0) {
    eB(&ud, &ag, &eBy, &error, verbose);
  } else {
    eBtminus(&ud, &ag, &eBy, &error, verbose);
  }
  printf("eBy = %g\terror = %g\trelerror = %g%%\n", eBy, error, fabs(error/eBy*100.0));
  */
  /*
  // 计算单点平均磁场
  ud.x = 0.0;
  ud.y = 0.0;
  ud.t = 0.0;
  ud.method = 2;
  mean_eB_homo(&ud, &ag, &mean_eBy);
  printf("<eBy> = %g\n", mean_eBy);
  */
  
 
  // 计算原点平均磁场随时间的变化
  int i, N;
  double tmin, tmax;
  ud.method = 2;
  ud.x = 0.0;
  ud.y = 0.0;
  N = 200;
  tmin = 0.0;
  tmax = 3.0;
  printf("# t       \teBy\n");
  for (i = 0; i <= N; i++) {
    ud.t = tmin + (tmax - tmin)*i/N;
    mean_eB_homo(&ud, &ag, &mean_eBy);
    printf("  %-8g\t%-8g\n", ud.t, mean_eBy);
  }
  
  return 0;
}
