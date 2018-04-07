#ifndef GAUGE_CONF_MULTI_C
#define GAUGE_CONF_MULTI_C

#include"../include/macro.h"

#include<complex.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
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
  int i, num_hit;
  long r;

  GAUGE_GROUP matrix;

  if(param->d_dist_poly>1 && param->d_size[1]-param->d_dist_poly>1) // Polyakov loops are separated along the "1" direction
    {
    num_hit=param->d_multihit;
    }
  else
    {
    num_hit=0;
    }

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r, i, matrix)
  #endif
  for(r=0; r<param->d_space_vol/2; r++)
     {
     one(&(loc_poly[r]));
     for(i=0; i<dt; i++)
        {
        multihit(GC,
                 geo,
                 param,
                 sisp_and_t_to_si(r, t_start+i, param),
                 0,
                 num_hit,
                 &matrix);
        times_equal(&(loc_poly[r]), &matrix);
        }
     }

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r, i, matrix)
  #endif
  for(r=param->d_space_vol/2; r<param->d_space_vol; r++)
     {
     one(&(loc_poly[r]));
     for(i=0; i<dt; i++)
        {
        multihit(GC,
                 geo,
                 param,
                 sisp_and_t_to_si(r, t_start+i, param),
                 0,
                 num_hit,
                 &matrix);
        times_equal(&(loc_poly[r]), &matrix);
        }
     }
  }


// single update of a slice
void slice_single_update(Gauge_Conf * GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         int t_start,
                         int dt)
  {
  long r;
  int i, dir, upd_over;

  // heatbath
  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol/2; r++) heatbath(GC, geo, param, sisp_and_t_to_si(r, t_start, param), 0);

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=param->d_space_vol/2; r<param->d_space_vol; r++) heatbath(GC, geo, param, sisp_and_t_to_si(r, t_start, param), 0);

  for(i=1; i<dt; i++)
     {
     for(dir=0; dir<STDIM; dir++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol/2; r++) heatbath(GC, geo, param, sisp_and_t_to_si(r, t_start+i, param), dir);

        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=param->d_space_vol/2; r<param->d_space_vol; r++) heatbath(GC, geo, param, sisp_and_t_to_si(r, t_start+i, param), dir);
        }
     }

  // overrelaxation
  for(upd_over=0; upd_over<param->d_overrelax; upd_over++)
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_space_vol/2; r++) overrelaxation(GC, geo, param, sisp_and_t_to_si(r, t_start, param), 0);

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=param->d_space_vol/2; r<param->d_space_vol; r++) overrelaxation(GC, geo, param, sisp_and_t_to_si(r, t_start, param), 0);

     for(i=1; i<dt; i++)
        {
        for(dir=0; dir<STDIM; dir++)
           {
           #ifdef OPENMP_MODE
           #pragma omp parallel for num_threads(NTHREADS) private(r)
           #endif
           for(r=0; r<param->d_space_vol/2; r++) overrelaxation(GC, geo, param, sisp_and_t_to_si(r, t_start+i, param), dir);

           #ifdef OPENMP_MODE
           #pragma omp parallel for num_threads(NTHREADS) private(r)
           #endif
           for(r=param->d_space_vol/2; r<param->d_space_vol; r++) overrelaxation(GC, geo, param, sisp_and_t_to_si(r, t_start+i, param), dir);
           }
        }
     }

  // unitarize
  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol; r++) unitarize( &(GC->lattice[sisp_and_t_to_si(r, t_start, param)][0]) );

  for(i=1; i<dt; i++)
     {
     for(dir=0; dir<STDIM; dir++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++) unitarize( &(GC->lattice[sisp_and_t_to_si(r, t_start+i, param)][dir]) );
        }
     }
  }


// multilevel
void multilevel(Gauge_Conf * GC,
                Geometry const * const geo,
                GParam const * const param,
                int t_start,
                int dt)
  {
  int level=-2; // initialized just to avoid warnings

  // determine the level to be used
  if(dt>param->d_ml_step[0])
    {
    level=-1;
    }
  else
    {
    int tmp;
    for(tmp=0; tmp<NLEVELS; tmp++)
       {
       if(param->d_ml_step[tmp]==dt)
         {
         level=tmp;
         tmp=NLEVELS+10;
         }
       }
    if(tmp==NLEVELS)
      {
      fprintf(stderr, "Error in the determination of the level in the multilevel (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }

  // LEVEL -1, do not average
  if(level == -1)
    {
    int i;
    long int r;

    // initialyze ml_polycorr_ris[0] to 1
    for(r=0; r<param->d_space_vol; r++)
       {
       one_TensProd(&(GC->ml_polycorr_ris[0][r]));
       }

    // call lower levels
    for(i=0; i<(param->d_size[0])/(param->d_ml_step[0]); i++)
       {
       multilevel(GC,
                  geo,
                  param,
                  t_start+i*param->d_ml_step[0],
                  param->d_ml_step[0]);
       }
    } // end of the outermost level

  else if(level == NLEVELS-1) // INNERMOST LEVEL
    {
    int i, upd;
    TensProd TP;
    long int r, r1, r2;

    // initialize ml_polycorr_tmp[level] to 0
    for(r=0; r<param->d_space_vol; r++)
       {
       zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
       }

    // perform the update
    for(upd=0; upd< param->d_ml_upd[level]; upd++)
       {
       slice_single_update(GC,
                           geo,
                           param,
                           t_start,
                           dt);

       // compute Polyakov loop restricted to the slice
       GAUGE_GROUP *loc_poly = (GAUGE_GROUP *) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(GAUGE_GROUP));

       compute_local_poly(GC,
                          geo,
                          param,
                          t_start,
                          dt,
                          loc_poly);

       // compute the tensor products
       // and update ml_polycorr_tmp[level]
       for(r=0; r<param->d_space_vol; r++)
          {
          int t_tmp, dir=1;

          r1=sisp_and_t_to_si(r, 0, param);
          for(i=0; i<param->d_dist_poly; i++) r1=nnp(geo, r1, dir);
          si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1

          TensProd_init(&TP, &(loc_poly[r]), &(loc_poly[r2]) );
          plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &TP);
          }

       free(loc_poly);
       } // end of update

    // normalize polycorr_tmp
    for(r=0; r<param->d_space_vol; r++)
       {
       times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
       }

    // update polycorr_ris
    for(r=0; r<param->d_space_vol; r++)
       {
       times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
       }
    } // end of innermost level

  #if NLEVELS>1
  else // NOT THE INNERMOST NOT THE OUTERMOST LEVEL
    {
    if(level==-1 || level==NLEVELS-1)
      {
      fprintf(stderr, "Error in the multilevel (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    int i, upd;
    long int r;

    // initialize ml_polycorr_tmp[level] to 0
    for(r=0; r<param->d_space_vol; r++)
       {
       zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
       }

    // perform the update
    for(upd=0; upd< param->d_ml_upd[level]; upd++)
       {
       slice_single_update(GC,
                           geo,
                           param,
                           t_start,
                           dt);

       // initialyze ml_polycorr_ris[level+1] to 1
       for(r=0; r<param->d_space_vol; r++)
          {
          one_TensProd(&(GC->ml_polycorr_ris[level+1][r]));
          }

       // call higher levels
       for(i=0; i<(param->d_ml_step[level])/(param->d_ml_step[level+1]); i++)
          {
          multilevel(GC,
                     geo,
                     param,
                     t_start+i*param->d_ml_step[level+1],
                     param->d_ml_step[level+1]);
          }

       // update polycorr_tmp[level] with polycorr_ris[level+1]
       for(r=0; r<param->d_space_vol; r++)
          {
          plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &(GC->ml_polycorr_ris[level+1][r]));
          }

       } // end of update

    // normalize polycorr_tmp[level]
    for(r=0; r<param->d_space_vol; r++)
       {
       times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
       }

    // update polycorr_ris[level]
    for(r=0; r<param->d_space_vol; r++)
       {
       times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
       }
    } // end of the not innermost not outermost level
  #endif
  } // end of the multilevel


#endif

