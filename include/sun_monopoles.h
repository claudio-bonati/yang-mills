#ifndef SUN_MONOPOLES_H
#define SUN_MONOPOLES_H

#include"macro.h"

#include<complex.h>
#include<openssl/md5.h>
#include<stdio.h>

#include"gparam.h"
#include"geometry.h"
#include"su2.h"


void comp_MAG_gauge_transformation_SuN (SuN X_links[2*STDIM],
                                        double lambda[NCOLOR],
                                        double OverRelaxParam,
                                        SuN *G_mag);

void comp_non_diagonal_contribution_SuN (SuN X_links[2*STDIM], 
                                         double *lambda,
                                         double *non_diag_contr);
  

void comp_functional_fmag_SuN(SuN X_links[2*STDIM],
                              double lambda[NCOLOR],
                              double *fmag);

#endif
