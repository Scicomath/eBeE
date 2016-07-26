#include <stdio.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "rhoFun.h"
#include "sq.h"


/***************************************************
 * rhoFun: 计算(x_p,y_p,z_p)处的核数密度，核数密度按照均分布 *
 ***************************************************/

double rhoFun(double x_p, double y_p, double z_p, double R, double b, double d, double n0, double Y, char flag) {
  // x_p, y_p, z_p: 源点坐标
  // R: 核的半径
  // b: 碰撞参量
  // d: 核几何参数
  // n0: 核中心数密度
  // Y: 核运动快度
  // flag: 核标记，'+'表示沿z轴正方向的核(即左核)，'-'表示沿z轴负方向的核(即右核)
  // 返回值: 相应核(x_p, y_p, z_p)处的数密度

  double rho;
  double len;
  double gamma;
  gamma = cosh(Y);
  
  // 判断是左核(+)还是右核(-)
  if (flag == '+')
    b = -b; // 左核，b反号，因为核心的x为负
  else if (flag == '-')
    ; // 右核，不变
  else
    printf("[rhoFun]error: 核类型(flag)应该为‘+’或者为'-'\n");

  len = sqrt(Sq(x_p - b/2) + Sq(y_p) + Sq(gamma*z_p));
  rho = gamma * n0 / ( 1 + exp((len-R)/d) ); // 乘以gamma因为密度变大

  return rho;
}
