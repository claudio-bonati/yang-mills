#ifndef SUN_MONOPOLES_H
#define SUN_MONOPOLES_H

#include"macro.h"

#include<complex.h>
#include<stdio.h>

#include"gparam.h"
#include"geometry.h"
#include"su2.h"

void comp_MAG_gauge_transformation_SuN(SuN X_links[2*STDIM],
                                       double const lambda[NCOLOR],
                                       double OverRelaxParam,
                                       SuN *G_mag);

void comp_outdiagnorm_of_X_SuN(SuN X_links[2*STDIM],
                               double const lambda[NCOLOR],
                               double *outdiagnorm2);

void comp_functional_fmag_SuN(SuN X_links[2*STDIM],
                              double const lambda[NCOLOR],
                              double *fmag);

void diag_projection_single_site_SuN(Gauge_Conf *GC,
                                     SuN *link,
                                     long r,
                                     int dir);
  
#endif
