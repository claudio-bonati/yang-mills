#ifndef SUN_MONOPOLES_C
#define SUN_MONOPOLES_C

#include"../include/macro.h"

#include<complex.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/son_monopoles.h"

// This function compute the MAG gauge transformation in the SU(N) case
void comp_MAG_gauge_transformation_SoN(SoN X_links[2*STDIM],
                                       double const lambda[NCOLOR],
                                       double overrelaxparam,
                                       SoN *G_mag)

   {
   (void) X_links;
   (void) lambda;
   (void) overrelaxparam;
   (void) G_mag;           // just to avoid warnings

   fprintf(stderr, "The function comp_MAG_gauge_transformation_SoN has to be written! (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }


// compute the squared absolute values of out-of-diagonal terms of X(n)
void comp_outdiagnorm_of_X_SoN(SoN X_links[2*STDIM],
                               double const lambda[NCOLOR],
                               double *outdiagnorm2)
   {
   (void) X_links;
   (void) lambda;
   (void) outdiagnorm2; // just to avoid warnings

   fprintf(stderr, "The function comp_outdiagnorm_of_X_SoN has to be written! (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }


// compute the value of the functional to be maximized in MAG
void comp_functional_fmag_SoN(SoN X_links[2*STDIM],
                              double const lambda[NCOLOR],
                              double *fmag)
   {
   (void) X_links;
   (void) lambda;
   (void) fmag; // just to avoid warnings

   fprintf(stderr, "The function comp_functional_fmag_SoN has to be written! (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }

// extract the diagonal part of the link and save it in GC->diag_proj
void diag_projection_single_site_SoN(Gauge_Conf *GC,
                                     SoN *link,
                                     long r,
                                     int dir)
   {
   (void) GC;
   (void) link;
   (void) r;
   (void) dir; // just to avoid warnings

   fprintf(stderr, "The function diag_projection_single_site_SoN has to be written! (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }


#endif




