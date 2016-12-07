#ifndef HIEBHPP

#define HIEBHPP

#include <cmath>
#include <string>

class HIeB
{
public:
  HIeB();
  double eBx, eBy, eBz;
  double x, y, z, t;
  double xp, yp, zp;
  double b; // impact parameter
  double sigma, sigmaChi; // electric conductivity and chiral magnetic conductivity
  void SetSqrtS(double sqrtS);
  double GetSqrtS() const;
  void SetGamma(double gamma);
  double GetGamma() const;
  void SetY0(double Y0);
  double GetY0() const;
  void SetY(double Y);
  double GetY() const;
  double GetV() const;
  void SetNucleiType(std::string nuclei);
  std::string GetNucleiType() const;
  // magnetic field of point charge(vaccum)
  void pointEB(double sign);
  // magnetic field of point charge(Wangqun's method)
  void pointEB_Wangqun(double sign);
  double rhoFun_Ai();
  double GetPointEB(int comp);
  
private:
  double mGamma; // Lorentz contraction factor
  double mSqrtS; // energy in center-of-mass frame
  double mY0; // initial rapidity
  double mY; // rapidity
  double mV; // velocity

  // nucleus parameter
  std::string mNucleiType;
  double mR; // radius of nucleus
  double md; // parameter of nucleus shape
  double mn0; // normalization coefficient
  double mZ; // electric charge of nucleus
  
  double ma = 0.5; // rapidity distribution parameter

  // point charge eB
  double peBx, peBy, peBz;

  char mRegionType, mFlag;

  int eB_Intgrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
  
};


#endif
