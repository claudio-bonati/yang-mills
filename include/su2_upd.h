#ifndef SU2_UPD_H
#define SU2_UPD_H

#include"su2.h"

void randheat_Su2(double k, double *out);
void single_heatbath_Su2(Su2 *link, Su2 const * const staple);
void single_overrelaxation_Su2(Su2 *link, Su2 const * const staple);
void cool_Su2(Su2 *link, Su2 const * const staple);

void single_overrelaxation_Su2Vecs(Su2Vecs *restrict link, Su2Vecs const * const staple);

#endif
