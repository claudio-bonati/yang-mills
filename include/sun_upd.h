#ifndef SUN_UPD_H
#define SUN_UPD_H

#include"gparam.h"
#include"su2.h"
#include"sun.h"

void single_heatbath_SuN(SuN *link, SuN const * const staple, GParam const * const param);
void single_heatbath_aux_SuN(SuN *link, SuN const * const staple, double beta);
void single_overrelaxation_SuN(SuN *link, SuN const * const staple);
void cool_SuN(SuN *link, SuN const * const staple);

#endif
