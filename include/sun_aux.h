#ifndef SUN_AUX_H
#define SUN_AUX_H

#include"su2.h"
#include"sun.h"

void ennetodue(SuN const * const in, int i, int j, double *xi, Su2 *u);
void duetoenne(Su2 const * const in, int i, int j, SuN *out);


#endif
