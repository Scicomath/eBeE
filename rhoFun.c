#include <stdio.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "rhoFun.h"
#include "sq.h"


/**************************************************************************
 * rhoFun_DK: 计算(x_p,y_p)处的核数密度，核数密度按照均匀分布,              *
 *            然后压缩为二维圆盘 (Kharzeev的方法)                          *
 *************************************************************************/

double rhoFun_DK(double x_p, double y_p, double R, double b, char flag) {
  // x_p, y_p: 横纵坐标
  // R: 核的半径
  // b: 碰撞参量
  // flag: 核标记，'+'表示沿z轴正方向的核(即左核)，'-'表示沿z轴负方向的核(即右核)
  // 返回值: 相应核(x_p, y_p)处的数密度

  double rho;
  double len;
  
  // 判断是左核(+)还是右核(-)
  if (flag == '+')
    b = -b; // 左核，b反号，因为核心的x为负
  else if (flag == '-')
    ; // 右核，不变
  else
    printf("[rhoFun]error: 核类型(flag)应该为‘+’或者为'-'\n");

  len = Sq(R) - (Sq(x_p - b/2) + Sq(y_p));
  if (len >= 0)
    rho = 2.0*sqrt(len) / (4.0 * M_PI * Pow3(R)/3.0);
  else // 如果len < 0 则表明在核外
    rho = 0.0;
  
  return rho;
}


/**************************************************************************
 * rhoFun_Mo: 计算(x_p,y_p)处的核数密度，核数密度按照woods-saxon分布        *
 *            然后压缩为二维圆盘 (莫玉俊方法)                              *
 * 注: 此函数本应有一个积分, 为方便其见, 将此积分与x,y积分合在一起          *
 ************************************************************************/

double rhoFun_Mo(double x_p, double y_p, double z_p, double R, double b,
		 double d, double n0, char flag) {
  // x_p, y_p, z_p: 源点坐标
  // R: 核的半径
  // b: 碰撞参量
  // d: 核几何分布参数
  // n0: 核中心数密度
  // flag: 核标记，'+'表示沿z轴正方向的核(即左核)，'-'表示沿z轴负方向的核(即右核)
  // 返回值: 相应核(x_p, y_p, z_p)处的数密度

  double rho;
  double len;
  
  // 判断是左核(+)还是右核(-)
  if (flag == '+')
    b = -b; // 左核，b反号，因为核心的x为负
  else if (flag == '-')
    ; // 右核，不变
  else
    printf("[rhoFun]error: 核类型(flag)应该为‘+’或者为'-'\n");

  len = sqrt(Sq(x_p - b/2) + Sq(y_p) + Sq(z_p));
  rho = n0 / ( 1 + exp((len-R)/d) );
  
  return rho;
}

double rhoFun_Mo_ori(double x_p, double y_p, double R, double b, char flag) {
  // x_p, y_p: 横纵坐标
  // R: 核的半径
  // b: 碰撞参量
  // flag: 核标记，'+'表示沿z轴正方向的核(即左核)，'-'表示沿z轴负方向的核(即右核)
  // 返回值: 相应核(x_p, y_p)处的数密度

  double rho;
  double len;
  
  // 判断是左核(+)还是右核(-)
  if (flag == '+')
    b = -b; // 左核，b反号，因为核心的x为负
  else if (flag == '-')
    ; // 右核，不变
  else
    printf("[rhoFun]error: 核类型(flag)应该为‘+’或者为'-'\n");

  len = Sq(R) - (Sq(x_p - b/2) + Sq(y_p));
  if (len >= 0)
    rho = 2.0*sqrt(len) * (0.17/(157.0*(1+exp((sqrt(len) - R)/0.54))));
  else // 如果len < 0 则表明在核外
    rho = 0.0;
  
  return rho;
}


/**************************************************************************
 * rhoFun_Ai: 计算(x_p,y_p,z_p)处的核数密度，核数密度按照woods-saxon分布    *
 *            然后压缩为三维椭球圆盘 (艾鑫方法)                             *
 **************************************************************************/

double rhoFun_Ai(double x_p, double y_p, double z_p, double R, double b,
		 double d, double n0, double Y, char flag) {
  // x_p, y_p, z_p: 源点坐标
  // R: 核的半径
  // b: 碰撞参量
  // d: 核几何分布参数
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
  rho = gamma * n0 / ( 1 + exp((len-R)/d) ); // 此处乘以gamma是因为密度变大
  
  return rho;
}


