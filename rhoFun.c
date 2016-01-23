#include <math.h>
#include "rhoFun.h"
#include "sqrt.h"

const double pi = 3.141592654;

/***************************************************
 * rhoFun: 计算(x_p,y_p)处的核数密度，核数密度按照均分布 *
 ***************************************************/

double rhoFun(double x_p, double y_p, double R, double b, char flag) {
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

  len = Sqrt(R) - (Sqrt(x_p - b/2) + Sqrt(y_p));
  if (len >= 0)
    rho = 2.0*len / (4.0*pi*Pow3(R)/3.0);
  else // 如果len < 0 则表明在核外
    rho = 0.0;

  return rho;
}
