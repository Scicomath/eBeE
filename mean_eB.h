#ifndef MEAN_EB_H
#define MEAN_EB_H

int mean_eB_homo(struct userdata *ud, const struct intargu *ag,
		 double *mean_eBy);
int mean_eB_inhomo(struct userdata *ud, const struct intargu *ag,
		   double *mean_eBy);

#endif
