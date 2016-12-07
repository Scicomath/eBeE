#include <cassert>
#include <iostream>
#include <string>
#include <cmath>
#include "HIeB.hpp"
#include "sq.h"

const double alpha_EM = 1.0 / 137.0;
const double m = 0.938272; // mass of proton
const double hbarc = 197.32696;
const double mpi = 140.0;
const double hbarcOverMpiSqure = (hbarc*hbarc)/(mpi*mpi);

HIeB::HIeB()
{
  t = 0.1;
  b = 8.0;
  sigma = 5.8/hbarc;
  sigmaChi = 1.5/hbarc;
  SetSqrtS(200.0);
}

void HIeB::SetSqrtS(double sqrtS)
{
  double E, p;
  
  assert (sqrtS >= 0.0);
  mSqrtS = sqrtS;

  // calculate corresponding rapidity
  E = sqrtS / 2.0;
  p = sqrt(E*E - m*m);
  mY0 = 0.5 * log((E+p)/(E-p));
  mY = mY0;

  // calculate corresponding Lorentz contraction factor
  mGamma = cosh(mY0);

  // calculate corresponding velocity
  mV = sqrt(1.0 - 1.0/Sq(mGamma));
}

double HIeB::GetSqrtS() const
{
  return mSqrtS;
}

void HIeB::SetGamma(double gamma)
{
  double p;
  assert (gamma >= 0.0);
  mGamma = gamma;

  // calculate corresponding rapidity
  mY0 = acosh(gamma);
  mY = mY0;

  // calculate corresponding sqrtS
  mV = tanh(mY0); // velocity in nature unit
  p = mGamma * m * mV;
  mSqrtS = 2.0 * sqrt(m*m + p*p);
}

double HIeB::GetGamma() const
{
  return mGamma;
}

void HIeB::SetY0(double Y0)
{
  double p;
  assert (Y0 >= 0.0);
  mY0 = Y0;
  mY = mY0;

  // calculate corresponding Lorentz contraction factor
  mGamma = cosh(mY0);

  // calculate corresponding center-of-mass frame
  mV = tanh(mY0);
  p = mGamma * m * mV;
  mSqrtS = 2.0 * sqrt(m*m + p*p);
}

double HIeB::GetY0() const
{
  return mY0;
}

void HIeB::SetY(double Y)
{
  assert(Y >= 0.0);
  mY = Y;
}

double HIeB::GetY() const
{
  return mY;
}

double HIeB::GetV() const
{
  return mV;
}

void HIeB::SetNucleiType(std::string nuclei)
{
  if (nuclei == "Au") {
    mNucleiType = "Au";
    mR = 6.38;
    md = 0.535;
    mn0 = 8.596268e-4;
    mZ = 79.0;
  } else if (nuclei == "Pb") {
    mNucleiType = "Pb";
    mR = 6.68;
    md = 0.535;
    mn0 = 7.51392e-4;
    mZ = 82.0;
  } else {
    printf("Error: undefind nuclei type!");
  }
}

std::string HIeB::GetNucleiType() const
{
  return mNucleiType;
}

void HIeB::pointEB(double sign)
{
  double denominator;

  denominator = pow( Sq(xp-x) + Sq(yp-y) + Sq(t*sign*sinh(mY)-z*cosh(mY)), 1.5);

  if (fabs(denominator) > 0.001) {
    peBx = sign * alpha_EM * sinh(mY0) *
      (yp-y)/denominator;
    peBy = sign * alpha_EM * sinh(mY0) *
      (x-xp)/denominator;
  } else {
    peBx = 0.0;
    peBy = 0.0;
  }

  peBz = 0.0;

}

void HIeB::pointEB_Wangqun(double sign)
{
  double delta, A, eBphi, eBr, xT, sqrtDelta;
  xT = hypot(x-xp, y-yp); // sqrt((x-xp)^2 + (y-yp)^2)

  delta = Sq(mGamma) * Sq(sign*mV*t-z) + Sq(xT);
  if (delta > 0.001) {
    sqrtDelta = sqrt(delta);
    A = sigma*mV*mGamma*0.5*(mGamma*fabs(sign*mV*t-z) - sqrtDelta);

    eBphi = hbarcOverMpiSqure*alpha_EM*
      ( sign*mV*mGamma*xT )
      /(pow(delta, 1.5))
      *(1 + sigma*mV*mGamma*0.5*sqrtDelta)*exp(A);

    eBr = -hbarcOverMpiSqure*sigmaChi*alpha_EM*0.5*(mV*Sq(mGamma)*xT)/(pow(delta, 1.5)) * (mGamma*fabs(sign*mV*t-z) + A*sqrtDelta)*exp(A);

    peBz = hbarcOverMpiSqure*sigmaChi*alpha_EM*0.5*(sign*mV*mGamma)/(pow(delta, 1.5))*(Sq(mGamma)*Sq(sign*mV*t-z)*(1 + sigma*mV*mGamma*0.5*sqrtDelta) + delta*(1 - sigma*mV*mGamma*0.5*sqrtDelta))*exp(A);

    peBx = eBr * (x - xp)/xT - eBphi * (y - yp)/xT;
    peBy = eBr * (y - yp)/xT + eBphi * (x - xp)/xT;
  } else {
    peBx = 0.0;
    peBy = 0.0;
    peBz = 0.0;
  }
  
}

double HIeB::rhoFun_Ai()
{
  double rho;
  double len;
  double newb = 0.0;

  // 判断是左核(+)还是右核(-)
  if (mFlag == '+')
    newb = -b;
  else if (mFlag == '-')
    newb = b;
  else
    std::cout << "Error: flag must be + or -!\n";

  len = sqrt(Sq(xp - newb/2.0) + Sq(yp) + Sq(mGamma*zp));
  rho = mGamma * mn0 / (1 + exp((len - mR)/md) ); // 此处乘以gamma是因为密度变大

  return rho;
}

int HIeB::eB_Intgrand(const int *ndim, const double xx[], const int *ncomp,
		double ff[], void *userdata)
{
  int i;
  double sign, jacobian, extra, rho;
  double Imin[4];
  double Imax[4];

  extra = 3.0; // 加上extra是因为wood-saxon分布并不是完全在半径为R的球内

  // 根据被积区域类型和核标记确定积分上下限
  if (mRegionType == 'p') {
    Imin[0] = -(mR - b/2.0);
    Imax[0] = mR - b/2.0;
    Imin[1] = -sqrt(Sq(mR) - Sq(b/2.0));
    Imax[1] = sqrt(Sq(mR) - Sq(b/2.0));
    Imin[2] = -(mR + extra) / mGamma;
    Imax[2] = (mR + extra) / mGamma;
    Imin[3] = -mY0;
    Imax[3] = mY0;
  } else {
    // 判断核标记
    if (mFlag == '+') {
      Imin[0] = -(mR + b/2.0) - extra;
      Imax[0] = 0.0;
      Imin[1] = -mR - extra;
      Imax[1] = mR + extra;
      Imin[2] = -(mR+extra) / mGamma;
      Imax[2] = (mR+extra) / mGamma;
    } else {
      Imin[0] = 0.0;
      Imax[0] = (mR + b/2.0) + extra;
      Imin[1] = -mR - extra;
      Imax[1] = mR + extra;
      Imin[2] = -(mR + extra) / mGamma;
      Imax[2] = (mR + extra) / mGamma;
    }
  }

  // 变量变换 x -> min + (max - min) * x 将积分区间变为0-1
  xp = Imin[0] + (Imax[0] - Imin[0]) * xx[0];
  yp = Imin[1] + (Imax[1] - Imin[1]) * xx[1];
  zp = Imin[2] + (Imax[2] - Imin[2]) * xx[2];
  if (*ndim == 4)
    mY = Imin[3] + (Imax[3] - Imin[3]) * xx[3];

  if (mFlag == '+')
    sign = 1.0;
  else
    sign = -1.0;

  rho = rhoFun_Ai();
  
  if (mRegionType == 'p') {
    // 判断是否在被积区域内
    if ( (Sq(xp + b/2.0) + Sq(yp) <= Sq(mR)) &&
	 (Sq(xp - b/2.0) + Sq(yp) <= Sq(mR)) ) {
      pointEB_Wangqun(sign);
      ff[0] = mZ * rho * peBx;
      ff[1] = mZ * rho * peBy;
      ff[2] = mZ * rho * peBz;
    } else {
      ff[0] = 0.0;
      ff[1] = 0.0;
      ff[2] = 0.0;
    }
  } else if (mRegionType == 's') {
    if ( Sq(xp - sign*b/2.0) + Sq(yp) >= Sq(mR) ) {
      pointEB(sign);
      ff[0] = mZ * rho * peBx;
      ff[1] = mZ * rho * peBy;
      ff[2] = mZ * rho * peBz;
    } else {
      ff[0] = 0.0;
      ff[1] = 0.0;
      ff[2] = 0.0;
    }
  }

  jacobian = 1.0;
  for (i = 0; i < *ndim; i++) {
    jacobian = jacobian * (Imax[i] - Imin[i]);
  }

  ff[0] *= jacobian;
  ff[1] *= jacobian;
  ff[2] *= jacobian;

  return 0;
  
}


double HIeB::GetPointEB(int comp)
{
  if (comp == 1) {
    return peBx;
  } else if (comp == 2) {
    return peBy;
  } else if (comp == 3) {
    return peBz;
  } else {
    std::cout << "Error: component must be 1(x) or 2(y) or 3(z)!\n";
  }

  return 0.0;
}
