#ifndef EBFUN_H
#define EBFUN_H

int eB(const double x, const double y, const double tau,
       const double R, const double b, const double Y0,
       const double a, const double Z,
       const struct intargu *ag, double *eBy, 
       double *totalerror, const int verbose);

#endif
