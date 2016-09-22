#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "udStruct.h"
#include "eBFun.h"
#include "sqrtStoY.h"


int main(int argc, char **argv)
{
  struct intargu ag;
  double eBy, error;
  int verbose;
  ag.nvec = 1;
  ag.epsrel = 1e-2;
  ag.epsabs = 1e-8;
  ag.flags = 0 | 8; // 0 | 8
  ag.seed = 0;
  ag.smineval = 1.5e4; // 1.5e5
  ag.smaxeval = 4e6;   // 1e6
  ag.pmineval = 2e4; // 2.5e5
  ag.pmaxeval = 1e7;   // 1e7

  ag.nstart = 1000;
  ag.nincrease = 500;
  ag.nbatch = 1000;
  ag.gridno = 0;

  ag.statefile = NULL;
  ag.spin = NULL;
  
  if (argc < 6) {
    printf("错误: 必须输入大于等于四个参数.\n Usage: eBeE Au/Pb sqrtS b Ai/DK/Mo origin/x-z t\n");
    return 0;
  }

  struct userdata ud;
  ud.x = 0.0;
  ud.y = 0.0;
  ud.t = 0.001;
  ud.R = 6.38; // Au: 6.38 Pb 6.68
  ud.b = 8.0;
  //ud.Y0 = sqrtStoY(200);//质心系能量200GeV // Au: 200 Pb: 2760
  ud.Y0 = sqrtStoY(atof(argv[2]));
  ud.d = 0.535; // Au: 0.535 Pb: 0.546
  ud.n0 = 8.596268e-4; // Au: 8.596268e-4 Pb: 7.51392e-4
  ud.a = 0.5;
  ud.Z = 79.0; // Au: 79.0 Pb: 82
  if (strcmp(argv[1], "Pb") == 0) {
    ud.R = 6.68;
    ud.Y0 = sqrtStoY(2760);
    ud.n0 = 7.51392e-4;
    ud.Z = 82.0;
  }
  ud.b = atof(argv[3]);
  printf("# parameter:\n");
  printf("# R = %g\n", ud.R);
  printf("# b = %g\n", ud.b);
  printf("# Y0 = %g\n", ud.Y0);
  printf("# Z = %g\n", ud.Z);

  ud.method = 2; // 0|1|2 : 分别表示Kharzeev,莫玉俊,艾鑫的方法  
  if (strcmp(argv[4], "Ai") == 0) {
    ud.method = 2;
  } else if (strcmp(argv[4], "DK") == 0) {
    ud.method = 0;
  } else if (strcmp(argv[4], "Mo") == 0) {
    ud.method = 1;
  } else {
    ud.method = 2;
  }


  if (strcmp(argv[5], "origin") == 0) {
  
    // 计算原点磁场随时间的变化
    int i, N;
    double tmin, tmax;
    verbose = 0;
    ud.x = 0.0;
    ud.y = 0.0;
    N = 500;
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
  }
  
  

  if (strcmp(argv[5], "x-z") == 0) {
    // 计算x-z平面
    int i, j, xN, zN;
    double xmin, xmax, zmin, zmax;
    xmin = -20.0;
    xmax = 20.0;
    xN = 200;
    zmin = -0.4;
    zmax = 0.4;
    zN = 200;
    ud.y = 0.0;
    verbose = 0;
    ud.t = atof(argv[6]);
    for (i = 0; i <= xN; i++) {
      for (j = 0; j <= zN; j++) {
	ud.x = xmin + (xmax - xmin)*i/xN;
	ud.z = zmin + (zmax - zmin)*j/zN;
	if (ud.t >= 0) {
	  eB(&ud, &ag, &eBy, &error, verbose);
	} else {
	  eBtminus(&ud, &ag, &eBy, &error, verbose);
	}
	printf("%-8g\t%-8g\t%-8g\t%-8g\t%-8g\n", ud.x, ud.z, eBy, error, fabs(error/eBy*100.0));
      }
    }
  }
  
  
  if (strcmp(argv[5], "originSingle") == 0) {
    verbose = 1;
    ud.x = 0.0;
    ud.y = 0.0;
    ud.t = atof(argv[6]);
    if (ud.t >= 0) {
      eB(&ud, &ag, &eBy, &error, verbose);
    } else {
      eBtminus(&ud, &ag, &eBy, &error, verbose);
    }
    printf("eBy = %g\terror = %g\trelerror = %g%%\n", eBy, error, fabs(error/eBy*100.0));
  }

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
  
  return 0;
}
