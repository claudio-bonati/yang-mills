#ifndef U1_MONOPOLES_H
#define U1_MONOPOLES_H

#include"u1.h"
#include"u1.h"
#include"macro.h"
#include<complex.h>
#include<openssl/md5.h>
#include<stdio.h>
#include"gparam.h"
#include"geometry.h"



void comp_non_diagonal_contribution_U1 (U1 X_links[2*STDIM], 
                                        double *lambda, 
                                        double *non_diag_contr);

void comp_MAG_gauge_transformation_U1 (U1 helper_X[2*STDIM],
                                       double lambda[NCOLOR],
                                       double OverRelaxParam,
                                       U1 *G_mag);


void comp_functional_fmag_U1 (U1 X_links[2*STDIM], 
                              double *lambda,
                              double *fmag);
 

#endif

