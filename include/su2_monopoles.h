#ifndef SU2_MONOPOLES_H
#define SU2_MONOPOLES_H

#include"gauge_conf.h"
#include"su2.h"
#include"sun.h"
#include"macro.h"

void comp_MAG_gauge_transformation_Su2(Su2 X_links[2*STDIM],
                                       double const lambda[2],
                                       double overrelaxparam,
                                       Su2 *G_mag);

void diagonalize_X_Su2_aux(double OverRelaxParam,
                           double X[3],
                           Su2 *G);

void comp_outdiagnorm_of_X_Su2(Su2 X_links[2*STDIM],
                               double const lambda[2],
                               double *non_diag_contr);

void comp_functional_fmag_Su2(Su2 X_links[2*STDIM], 
                              double const lambda[2],
                              double *fmag);

void diag_projection_single_site_Su2(Gauge_Conf *GC,
                                     Su2 *link,
                                     long r,
                                     int dir);

#endif

