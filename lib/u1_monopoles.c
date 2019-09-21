#ifndef U1_MONOPOLES_C
#define U1_MONOPOLES_C

#include"../include/macro.h"

#include<complex.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/sun_monopoles.h"
#include"../include/su2_monopoles.h"
#include"../include/u1_monopoles.h"


void comp_MAG_gauge_transformation_U1  (U1 X_links[2*STDIM],
                                        double lambda[NCOLOR],
                                        double OverRelaxParam,
                                        U1 *G_mag)
   {  
   (void)G_mag;
   (void)OverRelaxParam;
   (void)lambda[0];
   (void)X_links[0];
   }


void comp_non_diagonal_contribution_U1  (U1 X_links[2*STDIM], 
                                         double *lambda,
                                         double *non_diag_contr)

   {
   (void)non_diag_contr;
   (void)lambda[0];
   (void)X_links[0];
   }

void comp_functional_fmag_U1 (U1 X_links[2*STDIM], 
                              double *lambda,
                              double *fmag)

   {
   (void)fmag;
   (void) X_links[0];
   (void) lambda[0];
   }



 #endif




