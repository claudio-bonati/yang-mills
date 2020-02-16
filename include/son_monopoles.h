#ifndef SON_MONOPOLES_H
#define SON_MONOPOLES_H

#include"gauge_conf.h"
#include"macro.h"
#include"son.h"

void comp_MAG_gauge_transformation_SoN(SoN X_links[2*STDIM],
                                       double const lambda[NCOLOR],
                                       double OverRelaxParam,
                                       SoN *G_mag);

void comp_outdiagnorm_of_X_SoN(SoN X_links[2*STDIM],
                               double const lambda[NCOLOR],
                               double *outdiagnorm2);

void comp_functional_fmag_SoN(SoN X_links[2*STDIM],
                              double const lambda[NCOLOR],
                              double *fmag);

void diag_projection_single_site_SoN(Gauge_Conf *GC,
                                     SoN *link,
                                     long r,
                                     int dir);


#endif // SON_MONOPOLES_H

