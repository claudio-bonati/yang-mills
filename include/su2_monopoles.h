#ifndef SU2_MONOPOLES_H
#define SU2_MONOPOLES_H

#include"su2.h"
#include"sun.h"
#include"macro.h"
#include<complex.h>
#include<openssl/md5.h>
#include<stdio.h>
#include"gparam.h"
#include"geometry.h"



void comp_non_diagonal_contribution_Su2 (Su2 X_links[2*STDIM], 
                                         double lambda[2], 
                                         double *non_diag_contr);

void max_X_comp_G_Su2_aux (double OverRelaxParam, 
                           double X[3], 
                           Su2 *G);


void comp_MAG_gauge_transformation_Su2 (Su2 X_links[2*STDIM],
                                        double lambda[NCOLOR],
                                        double OverRelaxParam,
                                        Su2 *G_mag);



void comp_functional_fmag_Su2(Su2 X_links[2*STDIM], 
                              double lambda[2],
                              double *fmag);

#endif

