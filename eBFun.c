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
  
  static double x_p; // 积分变量源点横坐标
  static double y_p; // 积分变量源点纵坐标
  static double Y;   // 积分变量快度(注意不是初始快度Y0)
  static double Imin[3]; // 积分下限
  static double Imax[3]; // 积分上限
  
  // 根据被积区域类型和核标记确定积分上下限
  if (ud->type == 'p') {
    Imin[0] = -(ud->R - ud->b/2.0);
    Imax[0] = ud->R - ud->b/2.0;
    Imin[1] = -sqrt(Sq(ud->R) - Sq(ud->b/2.0));
    Imax[1] = sqrt(Sq(ud->R) - Sq(ud->b/2.0));
    Imin[2] = -ud->Y0;
    Imax[2] = ud->Y0;
  } else {
    // 判断核标记
    if (ud->flag == '+') {
      Imin[0] = -(ud->R + ud->b/2.0);
      Imax[0] = 0.0;
      Imin[1] = -ud->R;
      Imax[1] = ud->R;
    } else {
      Imin[0] = 0.0;
      Imax[0] = ud->R + ud->b/2.0;
      Imin[1] = -ud->R;
      Imax[1] = ud->R;
    }
  }
  // 变量变换 x -> min + (max - min) * x 将积分区间变为0-1
  x_p = Imin[0] + (Imax[0] - Imin[0]) * xx[0];
  y_p = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
  if (*ndim == 3)
    Y = Imin[2] + (Imax[2] - Imin[2]) * xx[2];

  

#define eB_y ff[0] // 磁场的y分量
#define eB_x ff[1] // 磁场的x分量
  
  // 判断被积区域类型
  if (ud->type == 'p') { // 对于参与者(p)
    denominator = (pow(Sq(x_p - ud->x) + Sq(y_p - ud->y) + Sq(ud->tau * sinh(Y)) ,1.5));
    // 判断是否在被积区域内
    if ( (Sq(x_p + ud->b/2.0) + Sq(y_p) <= Sq(ud->R)) &&
	 (Sq(x_p - ud->b/2.0) + Sq(y_p) <= Sq(ud->R)) ) {// fabs(denominator) > 0.001
      eB_y = f(Y,ud->Y0,ud->a) * sinh(Y) * rhoFun(x_p, y_p, ud->R, ud->b, ud->flag ) *
	(ud->x - x_p) / denominator;
    }
    else
      eB_y = 0.0;
  } else if (ud->type == 's') {
    // 判断符号正负
    if (ud->type == '+')
      sign = 1.0;
    else
      sign = -1.0;

    denominator = (pow(Sq(x_p - ud->x) + Sq(y_p - ud->y) + Sq(ud->tau * sinh(ud->Y0)) ,1.5));
    // printf("denominator = %5.3f \n", denominator);
    // 判断是否在被积区域内
    if ( (Sq(x_p + sign*ud->b/2.0) + Sq(y_p) <= Sq(ud->R)) &&
	 (Sq(x_p - sign*ud->b/2.0) + Sq(y_p) >= Sq(ud->R)) ) // fabs(denominator) > 0.001
      eB_y = rhoFun(x_p, y_p, ud->R, ud->b, ud->flag) *
	(ud->x - x_p) / denominator;
    else
      eB_y = 0.0;
    
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
 * eB: 计算(x,y)点处的磁场                           *
 **************************************************/
int eB(const double x, const double y, const double tau,
       const double R, const double b, const double Y0,
       const double a, const double Z,
       const struct intargu *ag, double *eBy) {
  int comp, nregions, neval, fail;
  double integral, error, prob;
  double eBp_plus, eBp_minus, eBs_plus, eBs_minus;
  double constant;
  static struct userdata ud;
  ud.x = x;
  ud.y = y;
  ud.tau = tau;
  ud.R = R;
  ud.b = b;
  ud.Y0 = Y0;
  ud.a = a;

  ud.type = 'p';
  ud.flag = '+';
  Vegas(3, 1, eB_Integrand, &ud, ag->nvec,
	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
	ag->mineval, ag->maxeval, ag->nstart, ag->nincrease,
	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
	&neval, &fail, &integral, &error, &prob);
  constant = Sq(hbarc) * Z * alpha_EM;
  eBp_plus = constant * integral;
  error = constant * error;
  printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
  printf("    eBp_plus:\t%.8f +- %.8f\tp = %.3f\n", eBp_plus, error, prob);

  ud.type = 'p';
  ud.flag = '-';
  Vegas(3, 1, eB_Integrand, &ud, ag->nvec,
	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
	ag->mineval, ag->maxeval, ag->nstart, ag->nincrease,
	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
	&neval, &fail, &integral, &error, &prob);
  constant = -constant;
  eBp_minus = constant * integral;
  error = constant * error;
  printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
  printf("    eBp_minus:\t%.8f +- %.8f\tp = %.3f\n", eBp_minus, error, prob);

  ud.type = 's';
  ud.flag = '+';
  Vegas(2, 1, eB_Integrand, &ud, ag->nvec,
	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
	ag->mineval, ag->maxeval, ag->nstart, ag->nincrease,
	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
	&neval, &fail, &integral, &error, &prob);
  constant = Sq(hbarc) * Z * alpha_EM * sinh(Y0);
  eBs_plus = constant * integral;
  error = constant * error;
  printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
  printf("    eBs_plus:\t%.8f +- %.8f\tp = %.3f\n", eBs_plus, error, prob);

  ud.type = 's';
  ud.flag = '-';
  Vegas(2, 1, eB_Integrand, &ud, ag->nvec,
	ag->epsrel, ag->epsabs, ag->flags, ag->seed,
	ag->mineval, ag->maxeval, ag->nstart, ag->nincrease,
	ag->nbatch, ag->gridno, ag->statefile, ag->spin,
	&neval, &fail, &integral, &error, &prob);
  constant = -constant;
  eBs_minus = constant * integral;
  error = constant * error;
  printf("Vegas result:\tneval %d\tfail %d\n", neval, fail);
  printf("   eBs_minus:\t%.8f +- %.8f\tp = %.3f\n", eBs_minus, error, prob);

  return 0;
}
