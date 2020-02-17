#ifndef SUN_UPD_H
#define SUN_UPD_H

#include"su2.h"
#include"sun.h"

void ennetodue_SuN(SuN const * const in, int i, int j, double *xi, Su2 *u);
void duetoenne_SuN(Su2 const * const in, int i, int j, SuN *out);

void single_heatbath_SuN(SuN *link, SuN const * const staple);
void single_overrelaxation_SuN(SuN *link, SuN const * const staple);
void cool_SuN(SuN *link, SuN const * const staple);

void single_overrelaxation_SuNVecs(SuNVecs *restrict link, SuNVecs const * const staple);


#endif
