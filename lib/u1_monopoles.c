#ifndef U1_MONOPOLES_C
#define U1_MONOPOLES_C

#include"../include/macro.h"

#include<complex.h>
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/sun_monopoles.h"
#include"../include/su2_monopoles.h"
#include"../include/u1_monopoles.h"


// This function compute the MAG gauge transformation in the SU(N) case
void comp_MAG_gauge_transformation_U1(U1 X_links[2*STDIM],
                                      double const lambda[1],
                                      double OverRelaxParam,
                                      U1 *G_mag)
   {  
   (void) G_mag;
   (void) OverRelaxParam;
   (void) lambda;
   (void) X_links;
   }


// compte the squared absolute values of out-of-diagonal terms of X(n)
void comp_outdiagnorm_of_X_U1(U1 X_links[2*STDIM],
                              double const lambda[1],
                              double *non_diag_contr)
   {
   (void)non_diag_contr;
   (void)lambda[0];
   (void)X_links[0];
   }


// compute the value of the functional to be maximized in MAG
void comp_functional_fmag_U1(U1 X_links[2*STDIM],
                             double const lambda[1],
                             double *fmag)
   {
   (void) fmag;
   (void) X_links;
   (void) lambda;
   }


// extract the diagonal part of the link and save it in GC->diag_proj
// in the U(1) we just need the argument of the complex number
void diag_projection_single_site_U1(Gauge_Conf *GC,
                                    U1 *link,
                                    long r,
                                    int dir)
   {
   double phi;
   
   phi = atan2(cimag(link->comp), creal(link->comp));
   (GC->diag_proj[r][dir][0]) = phi;
   }




 #endif




