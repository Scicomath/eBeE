#ifndef RHOFUN_H
#define RHOFUN_H

double rhoFun_DK(double x_p, double y_p, double R, double b, char flag);
double rhoFun_Mo(double x_p, double y_p, double z_p, double R, double b,
		 double d, double n0, char flag);
double rhoFun_Ai(double x_p, double y_p, double z_p, double R, double b,
		 double d, double n0, double Y, char flag);
#endif
