#include <stdio.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "cuba.h"
#include "udStruct.h"
#include "rhoFun.h"
#include "eBFun.h"
#include "sq.h"

const double alpha_EM = 1.0 / 137.0;
const double hbarc = 197.32696;

static int eB_Integrand(const int *ndim, const double xx[],
			 const int *ncomp, double ff[], void *userdata);
static inline double f(double Y, double Y0, double a);


/***************************************************
 * eB_Integrand: 磁场(eB)的被积函数。                 *
 ***************************************************/

static int eB_Integrand(const int *ndim, const double xx[],
			 const int *ncomp, double ff[], void *userdata) {

  struct userdata *ud = (struct userdata *) userdata;
  static double sign;
  static double jacobian;
  static double denominator;
  
  static double x_p; // 积分变量源点x坐标
  static double y_p; // 积分变量源点y坐标
  static double z_p; // 积分变量源点z坐标
  static double Y;   // 积分变量快度(注意不是初始快度Y0)
  static double Imin[4]; // 积分下限
  static double Imax[4]; // 积分上限
  static double extra; // 范围扩大因子

  extra = 3; // 加上extra是因为wood-saxon分布并不是完全在半径为R的球内
  
  // 根据被积区域类型和核标记确定积分上下限
  if (ud->type == 'p') {
    Imin[0] = -(ud->R - ud->b/2.0);
    Imax[0] = ud->R - ud->b/2.0;
    Imin[1] = -sqrt(Sq(ud->R) - Sq(ud->b/2.0));
    Imax[1] = sqrt(Sq(ud->R) - Sq(ud->b/2.0));
    Imin[2] = -(ud->R+extra); 
    Imax[2] = (ud->R+extra);
    Imin[3] = -ud->Y0;
    Imax[3] = ud->Y0;
  } else {
    // 判断核标记
    if (ud->flag == '+') {
      Imin[0] = -(ud->R + ud->b/2.0) - extra;
      Imax[0] = 0.0;
      Imin[1] = -ud->R - extra;
      Imax[1] = ud->R + extra;
      Imin[2] = -(ud->R+extra);
      Imax[2] = (ud->R+extra);
    } else {
      Imin[0] = 0.0;
      Imax[0] = (ud->R + ud->b/2.0) + extra;
      Imin[1] = -ud->R - extra;
      Imax[1] = ud->R + extra;
      Imin[2] = -(ud->R + extra);
      Imax[2] = (ud->R + extra);
    }
  }
  // 变量变换 x -> min + (max - min) * x 将积分区间变为0-1
  x_p = Imin[0] + (Imax[0] - Imin[0]) * xx[0];
  y_p = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
  z_p = Imin[2] + (Imax[2] - Imin[2]) * xx[2];
  if (*ndim == 4)
    Y = Imin[3] + (Imax[3] - Imin[3]) * xx[3];

  

#define eB_y ff[0] // 磁场的y分量
#define eB_x ff[1] // 磁场的x分量
  
  // 判断被积区域类型
  if (ud->type == 'p') { // 对于参与者(p)
    if (ud->flag == '+')
      sign = 1.0;
    else
      sign = -1.0;

    denominator = (pow(Sq(x_p - ud->x) + Sq(y_p - ud->y) + Sq(ud->t * sinh(Y) - ud->z*cosh(Y)) ,1.5));
    //denominator = (pow(Sq(x_p - ud->x) + Sq(y_p - ud->y) + Sq(ud->t * sinh(Y)) ,1.5));
    // 判断是否在被积区域内
    if ( (Sq(x_p + ud->b/2.0) + Sq(y_p) <= Sq(7.0)) &&//(Sq(x_p + ud->b/2.0) + Sq(y_p) <= Sq(ud->R)) &&
	 (Sq(x_p - ud->b/2.0) + Sq(y_p) <= Sq(7.0)) &&//(Sq(x_p - ud->b/2.0) + Sq(y_p) <= Sq(ud->R)) &&
	 (fabs(denominator) > 0.001) ) {
      eB_y = f(Y,ud->Y0,ud->a) * sinh(Y) * rhoFun(x_p, y_p, z_p, ud->R, ud->b, ud->d, ud->n0, ud->flag ) *
	(ud->x - x_p) / denominator;
    }
    else
      eB_y = 0.0;
  } else if (ud->type == 's') {
    // 判断符号正负
    if (ud->flag == '+')
      sign = 1.0;
    else
      sign = -1.0;

    denominator = (pow(Sq(x_p - ud->x) + Sq(y_p - ud->y) + Sq(ud->t * sign * sinh(ud->Y0) - ud->z*cosh(ud->Y0)) ,1.5));
    //denominator = (pow(Sq(x_p - ud->x) + Sq(y_p - ud->y) + Sq(ud->t * sinh(ud->Y0)) ,1.5));
    // 判断是否在被积区域内
    if ( //(Sq(x_p + sign*ud->b/2.0) + Sq(y_p) <= Sq(ud->R)) &&
	(Sq(x_p - sign*ud->b/2.0) + Sq(y_p) >= Sq(7.0)) &&//(Sq(x_p - sign*ud->b/2.0) + Sq(y_p) >= Sq(ud->R)) &&
	 (fabs(denominator) > 0.001) ) 
      eB_y = rhoFun(x_p, y_p, z_p, ud->R, ud->b, ud->d, ud->n0, ud->flag) *
	(ud->x - x_p) / denominator;
    else {
      eB_y = 0.0;
    }
    
  } else {
    printf("[eBpFun.c]error: 被积函数类型(type)只能为，参与者'p'或旁观者's'\n");
  }

  jacobian = 1.0;
  for (int i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }

  eB_y = jacobian * eB_y;

  return 0;
}

/**************************************************
 * f: 参与者碰撞后的快度分布函数                       *
 **************************************************/
static inline double f(double Y, double Y0, double a) {
  return (a*exp(a*Y)) / (2*sinh(a*Y0));
}


/**************************************************
 * eB: 计算(x,y,z)点处的磁场                         *
 **************************************************/
int eB(const double x, const double y, const double z, const double t,
       const double R, const double b, const double Y0, const double d,
       const double n0, const double a, const double Z,
       const struct intargu *ag, double *eBy, double *totalerror, const int verbose) {
  int comp, nregions, neval, fail;
  double integral, error, prob;
  double eBp_plus, eBp_minus, eBs_plus, eBs_minus;
  double constant;
  static struct userdata ud;
  ud.x = x;
  ud.y = y;
  ud.z = z;
  ud.t = t;
  ud.R = R;
  ud.b = b;
  ud.d = d;
  ud.n0 = n0;
  ud.Y0 = Y0;
  ud.a = a;
  ud.Z = Z;

  *totalerror = 0.0;

  ud.type = 'p';
  ud.flag = '+';
  Vegas(4, 1, eB_Integrand, &ud, ag->nvec,
  	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
  	ag->pmineval, ag->pmaxeval, ag->nstart, ag->nincrease,
  	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
  	&neval, &fail, &integral, &error, &prob);
  constant = Sq(hbarc) * Z * alpha_EM;
  eBp_plus = constant * integral;
  error = fabs(constant * error);
  *totalerror += error;
  if (verbose == 1) {
    printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
    printf("    eBp_plus:\t%.8f +- %.8f(%.2f%%)\tp = %.3f\n", eBp_plus, error, error/eBp_plus*100.0, prob);
  }

  ud.type = 'p';
  ud.flag = '-';
  Vegas(4, 1, eB_Integrand, &ud, ag->nvec,
  	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
  	ag->pmineval, ag->pmaxeval, ag->nstart, ag->nincrease,
  	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
  	&neval, &fail, &integral, &error, &prob);
  constant = -constant;
  eBp_minus = constant * integral;
  error = fabs(constant * error);
  *totalerror += error;
  if (verbose == 1) {
    printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
    printf("    eBp_minus:\t%.8f +- %.8f(%.2f%%)\tp = %.3f\n", eBp_minus, error, error/eBp_minus*100.0, prob);
  }

  ud.type = 's';
  ud.flag = '+';
  Vegas(3, 1, eB_Integrand, &ud, ag->nvec,
  	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
  	ag->smineval, ag->smaxeval, ag->nstart, ag->nincrease,
  	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
  	&neval, &fail, &integral, &error, &prob);
  constant = Sq(hbarc) * Z * alpha_EM * sinh(Y0);
  eBs_plus = constant * integral;
  error = fabs(constant * error);
  *totalerror += error;
  if (verbose == 1) {
    printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
    printf("    eBs_plus:\t%.8f +- %.8f(%.2f%%)\tp = %.3f\n", eBs_plus, error, error/eBs_plus*100.0, prob);
  }

  ud.type = 's';
  ud.flag = '-';
  Vegas(3, 1, eB_Integrand, &ud, ag->nvec,
  	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
  	ag->smineval, ag->smaxeval, ag->nstart, ag->nincrease,
  	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
  	&neval, &fail, &integral, &error, &prob);
  constant = -constant;
  eBs_minus = constant * integral;
  error = fabs(constant * error);
  *totalerror += error;
  if (verbose == 1) {
    printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
    printf("   eBs_minus:\t%.8f +- %.8f(%.2f%%)\tp = %.3f\n", eBs_minus, error, error/eBs_minus*100.0, prob);
  }

  *eBy = eBp_plus + eBp_minus + eBs_plus + eBs_minus;
  
  return 0;
}
