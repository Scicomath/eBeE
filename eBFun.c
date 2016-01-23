#include <gsl/gsl_math.h>
#include "cuba.h"
#include "usStruct.h"
#include "rhoFun.h"

/**************************************************
 * eBp_Integrand: 参与者(p)产生的磁场(eB_p)的被积函数。*
 **************************************************/

static int eBp_Integrand(const int *ndim, const double xx[],
			 const int *ncomp, double ff[], void *userdata) {

  struct userdata *ud = (struct userdata *) userdata;
  static double sign;
  static double jacobian = 1.0;
  // 变量变换 x -> min + (max - min) * x 将积分区间变为0-1
  static double x_p; // 积分变量源点横坐标
  static double y_p; // 积分变量源点纵坐标
  static double Y;   // 积分变量快度(注意不是初始快度Y0)
  x_p = ud->min[0] + (ud->max[0] - ud->min[0]) * xx[0];
  y_p = ud->min[1] + (ud->max[1] - ud->min[1]) * xx[1];
  if (*ndim == 3)
    Y = ud->min[2] + (ud->max[2] - ud->min[2]) * xx[2];

  

#define eB_y ff[0] // 磁场的y分量
#define eB_x ff[1] // 磁场的x分量
  
  // 判断被积区域类型
  if (ud->type == 'p') { // 对于参与者(p)
    // 判断是否在被积区域内
    if ( (Sqrt(x_p + b/2.0) + Sqrt(y_p) <= Sqrt(R)) &&
	 (Sqrt(x_p - b/2.0) + Sqrt(y_p) <= Sqrt(R)) )
      eB_y = f(Y,ud->Y0,ud->a) * sinh(Y) * rhoFun(x_p, y_p, ud->R, ud->b, ud->flag ) *
	(ud->x - x_p) / (pow(Sqrt(x_p - ud->x) + Sqrt(y_p - ud->y) + Sqrt(ud->tau * sinh(Y)) ,1.5));
    else
      eB_y = 0.0;
  } else if (ud->type == 's') {
    // 判断符号正负
    if (us-type == '+')
      sign = 1.0;
    else
      sign = -1.0;

    // 判断是否在被积区域内
    if ( (Sqrt(x_p + sign*b/2.0) + Sqrt(y_p) <= Sqrt(R)) &&
	 (Sqrt(x_p - sign*b/2.0) + Sqrt(y_p) >= Sqrt(R)) )
      eB_y = rhoFun(x_p, y_p, ud->R, ud->b, ud->flag) *
	(ud->x - x_p) / (pow(Sqrt(x_p - ud->x) + Sqrt(y_p - ud->y) + Sqrt(ud->tau * sinh(Y0)) ,1.5));
    else
      eB_y = 0.0;
    
  } else {
    printf("[eBpFun.c]error: 被积函数类型(type)只能为，参与者'p'或旁观者's'\n");
  }

  for (int i = 0; i < *ndim; i++) {
    jacobian = jacobian * (ud->max[i] - ud->min[i]);
  }

  eB_y = jacobian * eB_y;

  return 0;
}

/**************************************************
 * f: 参与者碰撞后的快度分布函数                         *
 **************************************************/
static inline double f(double Y, double Y0, double a) {
  return (a*exp(a*Y)) / (2*sinh(a*Y0));
}
