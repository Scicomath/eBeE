#include <math.h>
#include "sq.h"
#include "sqrtStoY.h"

/******************************
 * 将质心系能量转换为快度      *
 ******************************/

double sqrtStoY(double sqrtS) {
  double E, p, m, Y;
  m = 0.938272;
  E = sqrtS / 2;
  p = sqrt(Sq(E) - Sq(m));
  Y = 0.5 * log((E+p)/(E-p));
  return Y;
}
