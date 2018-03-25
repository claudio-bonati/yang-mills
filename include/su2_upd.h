#ifndef SU2_UPD_H
#define SU2_UPD_H

#include"gparam.h"
#include"su2.h"

void randheat_Su2(double k, double *out);
void single_heatbath_Su2(Su2 *link, Su2 const * const staple, GParam const * const param);
void single_overrelaxation_Su2(Su2 *link, Su2 const * const staple);
void cool_Su2(Su2 *link, Su2 const * const staple);

#endif
