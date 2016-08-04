#include <stdio.h>
#include <math.h>
#include "udStruct.h"
#include "eBFun.h"
#include "sq.h"
#include "mean_eB.h"
/* 
   mean_eB_homo: 计算(x,y)处沿z轴的平均磁场， 假设核物质密度均匀分布
 */

int mean_eB_homo(struct userdata *ud, const struct intargu *ag,
	    double *mean_eBy) {
  double zmin, zmax, lenleft, lenright, v, gamma, range, eBy,
    error;
  int N, verbose;
  if (ud->t < 0.0) {
    printf("Error: To calculate average magnetic field, the t must be equal or greater than zero.\n");
    *mean_eBy = 0.0;
    return 0;
  }
  // 点(x,y)到左核， 即z轴正向运动的核的距离
  lenleft = sqrt(Sq(ud->R) - Sq(ud->x + ud->b/2.0) - Sq(ud->y));
  // 点(x,y)到右核， 即z轴负向运动的核的距离
  lenright = sqrt(Sq(ud->R) - Sq(ud->x - ud->b/2.0) - Sq(ud->y));
  v = tanh(ud->Y0);
  gamma = cosh(ud->Y0);
  zmin = fmin(v*ud->t-lenleft/gamma, -v*ud->t-lenleft/gamma);
  zmax = fmax(v*ud->t+lenleft/gamma, -v*ud->t+lenleft/gamma);
  range = zmax - zmin;
  N = 100;
  verbose = 0;
  *mean_eBy = 0.0;
  for (int i = 0; i <= N; i++) {
    ud->z = zmin + range*i/N;
    eB(ud, ag, &eBy, &error, verbose);
    printf("z = %g\teBy = %g\terror = %g\n", ud->z, eBy, error);
    *mean_eBy += eBy;
  }
  *mean_eBy /= (N+1);
  return 0;
}

/*
  mean_eB_inhomo: 计算(x,y)处沿z轴的平均磁场， 假设核物质密度均匀分布
 */

int mean_eB_inhomo(struct userdata *ud, const struct intargu *ag,
		   double *mean_eBy) {
  return 0;
}
