#ifndef SON_UPD_H
#define SON_UPD_H

#include"son.h"

void ennetodue_SoN(SoN const * const in, int i, int j, double *xi, double *theta);
void duetoenne_SoN(double theta, int i, int j, SoN *out);

void single_heatbath_SoN(SoN *link, SoN const * const staple);
void single_overrelaxation_SoN(SoN *link, SoN const * const staple);
void cool_SoN(SoN *link, SoN const * const staple);

void single_overrelaxation_SoNVecs(SoNVecs *restrict link, SoNVecs const * const staple);


#endif // SON_UPD_H

