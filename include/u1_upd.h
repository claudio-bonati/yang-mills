#ifndef U1_UPD
#define U1_UPD

#include"u1.h"

void single_heatbath_U1(U1 * restrict link, U1 const * const restrict staple);
void single_overrelaxation_U1(U1 * restrict link, U1 const * const restrict staple);
void cool_U1(U1 * restrict link, U1 const * const restrict staple);

void single_overrelaxation_U1Vecs(U1Vecs *restrict link, U1Vecs const * const staple);

#endif // U1_UPD

