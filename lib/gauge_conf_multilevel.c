#ifndef GAUGE_CONF_MULTI_C
#define GAUGE_CONF_MULTI_C

#include<complex.h>
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/macro.h"
#include"../include/mymalloc.h"
#include"../include/tens_prod.h"

void multihit(Gauge_Conf const * const restrict GC,
              Geometry const * const restrict geo,
              GParam const * const restrict param,
              long r,
              int dir,
              int num_hit,
              GAUGE_GROUP *G)
  {
  if(num_hit>0)
    {
    int i;
    const double inv_num_hit = 1.0/(double) num_hit;
    GAUGE_GROUP staple, partial;

    zero(G);
    equal(&partial, &(GC->lattice[r][dir]));
    calcstaples_w(GC, geo, r, dir, &staple);

    for(i=0; i<num_hit; i++)
       {
       single_heatbath(&partial, &staple, param);
       plus_equal(G, &partial);

       single_overrelaxation(&partial, &staple);
       plus_equal(G, &partial);

       single_overrelaxation(&partial, &staple);
       plus_equal(G, &partial);

       single_overrelaxation(&partial, &staple);
       plus_equal(G, &partial);

       single_overrelaxation(&partial, &staple);
       plus_equal(G, &partial);

       unitarize(&partial);
       }
    times_equal_real(G, 0.2*inv_num_hit);
    }
  else
    {
    equal(G, &(GC->lattice[r][dir]));
    (void) *param;
    }
  }

// compute polyakov loop on a single slice
void compute_local_poly(Gauge_Conf const * const restrict GC,
                        Geometry const * const restrict geo,
                        GParam const * const restrict param,
                        int t_start,
                        int dt,
                        GAUGE_GROUP *loc_poly)
  {
  int i;
  long r;

  GAUGE_GROUP matrix;

  for(r=0; r<param->d_space_vol; r++)
     {
     // polyakov loop correlators     // note: if param->d_dist_poly<2 the multihit is never used
                                      //       see the first line of the multihit function
     one(&(loc_poly[r]));
     for(i=0; i<dt; i++)
        {
        multihit(GC,
                 geo,
                 param,
                 sisp_and_t_to_si(r, t_start+i, param),
                 0,
                 param->d_multihit,
                 &matrix);
        times_equal(&(loc_poly[r]), &matrix);
        }
     }// end of the loop on r
  }


// multilevel with just level 1
void multilevel1(Gauge_Conf * restrict GC,
                 Geometry const * const restrict geo,
                 GParam const * const restrict param,
                 int t_start,
                 int dt)
  {
  const double inv_upd_level1= 1.0/(double) param->d_up_level1;

  // LEVEL0 (maybe not used), do not average
  if(dt > param->d_ml_step1)
    {
    int i;
    long int r;

    // initialyze ml_polycorr_ris_level1 to 1
    for(r=0; r<param->d_ml_size; r++)
       {
       one_TensProd(&(GC->ml_polycorr_ris_level1[r]));
       }

    // call lower levels
    for(i=0; i<(param->d_size[0])/(param->d_ml_step1); i++)
       {
       multilevel1(GC,
                   geo,
                   param,
                   t_start+i*param->d_ml_step1,
                   param->d_ml_step1);
       }
    }

  // LEVEL1
  if(dt == param->d_ml_step1)
    {
    int i, dir, upd, upd_over;
    TensProd TP;
    long int r, r1, r2;

    // initialyze ml_polycorr_ris_level1 to 1 in case it is needed
    if(param->d_ml_step1==param->d_size[0])
      {
      for(r=0; r<param->d_ml_size; r++)
         {
         one_TensProd(&(GC->ml_polycorr_ris_level1[r]));
         }
      }

    // initialize ml_polycorr_tmp_level1 to 0
    for(r=0; r<param->d_ml_size; r++)
       {
       zero_TensProd(&(GC->ml_polycorr_tmp_level1[r]));
       }

    // perform the update
    for(upd=0; upd< param->d_up_level1; upd++)
       {
       // heatbath
       for(r=0; r<param->d_space_vol; r++) heatbath_w(GC, geo, param, sisp_and_t_to_si(r, t_start, param), 0);
       for(i=1; i<dt; i++)
          {
          for(dir=0; dir<STDIM; dir++)
             {
             for(r=0; r<param->d_space_vol; r++) heatbath_w(GC, geo, param, sisp_and_t_to_si(r, t_start+i, param), dir);
             }
          }

       // overrelaxation
       for(upd_over=0; upd_over<param->d_overrelax; upd_over++)
          {
          for(r=0; r<param->d_space_vol; r++) overrelaxation_w(GC, geo, param, sisp_and_t_to_si(r, t_start, param), 0);
          for(i=1; i<dt; i++)
             {
             for(dir=0; dir<STDIM; dir++)
                {
                for(r=0; r<param->d_space_vol; r++) overrelaxation_w(GC, geo, param, sisp_and_t_to_si(r, t_start+i, param), dir);
                }
             }
          }

       // unitarize
       for(r=0; r<param->d_space_vol; r++) unitarize( &(GC->lattice[sisp_and_t_to_si(r, t_start, param)][0]) );
       for(i=1; i<dt; i++)
          {
          for(dir=0; dir<STDIM; dir++)
             {
             for(r=0; r<param->d_space_vol; r++) unitarize( &(GC->lattice[sisp_and_t_to_si(r, t_start+i, param)][dir]) );
             }
          }

       // compute Polyakov loop restricted to the slice
       GAUGE_GROUP *loc_poly = (GAUGE_GROUP *) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(GAUGE_GROUP));

       compute_local_poly(GC,
                          geo,
                          param,
                          t_start,
                          dt,
                          loc_poly);

       // compute the tensor products
       // and update ml_polycorr_tmp_level2
       for(r=0; r<param->d_space_vol; r++)
          {
          int t_tmp, dir;

          for(dir=1; dir<STDIM; dir++)
             {
             r1=sisp_and_t_to_si(r, 0, param);
             for(i=0; i<param->d_dist_poly; i++) r1=nnp(geo, r1, dir);
             si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1

             TensProd_init(&TP, &(loc_poly[r]), &(loc_poly[r2]) );
             plus_equal_TensProd(&(GC->ml_polycorr_tmp_level1[r+(dir-1)*param->d_space_vol]), &TP);
             }
          }

       free(loc_poly);
       } // end of update

    // normalize tmp_level1
    for(r=0; r<param->d_ml_size; r++)
       {
       times_equal_real_TensProd(&(GC->ml_polycorr_tmp_level1[r]), inv_upd_level1);
       }

    // update ml_ris_level1
    for(r=0; r<param->d_ml_size; r++)
       {
       times_equal_TensProd(&(GC->ml_polycorr_ris_level1[r]), &(GC->ml_polycorr_tmp_level1[r]));
       }
    } // end of level 1
  }




#endif

