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

void multihit(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
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
    calcstaples_wilson(GC, geo, r, dir, &staple);

    for(i=0; i<num_hit; i++)
       {
       single_heatbath(&partial, &staple, param);
       plus_equal(G, &partial);

       unitarize(&partial);
       }
    times_equal_real(G, inv_num_hit);
    }
  else
    {
    equal(G, &(GC->lattice[r][dir]));
    }
  }


// compute polyakov loop on a single slice
void compute_local_poly(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
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
  for(r=0; r<param->d_space_vol; r++)
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

  #ifdef OPENMP_MODE
  if(geo->indexing_type!=0)
    {
    fprintf(stderr, "Wrong indexing used! (indexing_type=%d) (%s, %d)\n", geo->indexing_type, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

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


// multilevel for polyakov QbarQ correlator
void multilevel_pot_QbarQ(Gauge_Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int t_start,
                          int dt)
  {
  int i, upd;
  long int r;
  int level;

  level=-2;
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
    if(level==-2)
      {
      fprintf(stderr, "Error in the determination of the level in the multilevel (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }

  switch(level)
    {
    case -1 :     // LEVEL -1, do not average
      // initialyze ml_polycorr_ris[0] to 1
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         one_TensProd(&(GC->ml_polycorr_ris[0][r]));
         }

      // call lower levels
      for(i=0; i<(param->d_size[0])/(param->d_ml_step[0]); i++)
         {
         multilevel_pot_QbarQ(GC,
                    geo,
                    param,
                    t_start+i*param->d_ml_step[0],
                    param->d_ml_step[0]);
         }
      break;
      // end of the outermost level

    case NLEVELS-1 : // INNERMOST LEVEL
      // initialize ml_polycorr_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
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
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            TensProd TP;
            long r1, r2;
            int t_tmp ;

            r1=sisp_and_t_to_si(r, 0, param);
            for(i=0; i<param->d_dist_poly; i++) r1=nnp(geo, r1, 1);
            si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1

            TensProd_init(&TP, &(loc_poly[r]), &(loc_poly[r2]) );
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &TP);
            }

        free(loc_poly);
        } // end of update

      // normalize polycorr_tmp
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      // update polycorr_ris
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         }

      break;
      // end of innermost level

    default:  // NOT THE INNERMOST NOT THE OUTERMOST LEVEL
      if(level==-1 || level==NLEVELS-1)
        {
        fprintf(stderr, "Error in the multilevel (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }

      // initialize ml_polycorr_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
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
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            one_TensProd(&(GC->ml_polycorr_ris[level+1][r]));
            }

         // call higher levels
         for(i=0; i<(param->d_ml_step[level])/(param->d_ml_step[level+1]); i++)
            {
            multilevel_pot_QbarQ(GC,
                       geo,
                       param,
                       t_start+i*param->d_ml_step[level+1],
                       param->d_ml_step[level+1]);
            }

         // update polycorr_tmp[level] with polycorr_ris[level+1]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &(GC->ml_polycorr_ris[level+1][r]));
            }

         } // end of update

      // normalize polycorr_tmp[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      // update polycorr_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         }
      break;
      // end of the not innermost not outermost level

    } // end of switch
  } // end of multilevel


// multilevel for polyakov QbarQ correlator to be used in long simulations
void multilevel_pot_QbarQ_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int t_start,
                               int dt,
                               int iteration)
  {
  int i, upd;
  long int r;
  int level;

  level=-2;
  // determine the level to be used
  if(dt>param->d_ml_step[0])
    {
    fprintf(stderr, "Value of 'dt' not permitted in multilevel_pot_QbarQ_long (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
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
    if(level==-2)
      {
      fprintf(stderr, "Error in the determination of the level in the multilevel (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }

  if(level==0)
    {
    // initialize ml_polycorr_ris[0] to 1 if needed
    if(t_start==0 && iteration==0)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         one_TensProd(&(GC->ml_polycorr_ris[0][r]));
         }
      }

    if(iteration==0)
      {
      // initialize ml_polycorr_tmp[0] to 0 when needed
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         zero_TensProd(&(GC->ml_polycorr_tmp[0][r]));
         }
      }

    if(NLEVELS==1)
      {
      // perform the update
      for(upd=0; upd< param->d_ml_upd[0]; upd++)
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
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            TensProd TP;
            long r1, r2;
            int t_tmp;

            r1=sisp_and_t_to_si(r, 0, param);
            for(i=0; i<param->d_dist_poly; i++) r1=nnp(geo, r1, 1);
            si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1

            TensProd_init(&TP, &(loc_poly[r]), &(loc_poly[r2]) );
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &TP);
            }

         free(loc_poly);
         } // end of update
      }
    else  // NLEVELS!=1
      {
      // initialyze ml_polycorr_ris[1] to 1
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         one_TensProd(&(GC->ml_polycorr_ris[1][r]));
         }

      // call inner levels
      for(i=0; i<(param->d_size[0])/(param->d_ml_step[1]); i++)
         {
         multilevel_pot_QbarQ(GC,
                    geo,
                    param,
                    t_start+i*param->d_ml_step[1],
                    param->d_ml_step[1]);
         }
      }

    if(iteration==param->d_ml_level0_repeat-1) // iteration starts from zero
      {
      // normalize polycorr_tmp[0]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[0][r]), 1.0/( (double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat));
         }

      // update polycorr_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[0][r]), &(GC->ml_polycorr_tmp[0][r]));
         }
      }
    }
  else // level != 0
    {
    multilevel_pot_QbarQ(GC,
                         geo,
                         param,
                         t_start,
                         dt);
    }
  } // end of multilevel


// compute plaquettes on the t=1 time slice
void compute_plaq_on_slice1(Gauge_Conf const * const GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            double **plaq)
  {
  long r;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol; r++)
     {
     int i, j, tmp;
     long r4;

     tmp=0;
     r4=sisp_and_t_to_si(r, 1, param); // t=1

     i=0;
     for(j=1; j<STDIM; j++)
        {
        plaq[r][tmp]=plaquettep(GC, geo, r4, i, j);
        tmp++;
        }

     for(i=1; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           plaq[r][tmp]=plaquettep(GC, geo, r4, i, j);
           tmp++;
           }
        }
     #ifdef DEBUG
     if(tmp!=STDIM*(STDIM-1)/2)
       {
       fprintf(stderr, "Error in computation of the plaquettes in multilevel (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     #endif
     }
  }


// multilevel for polyakov string width of QbarQ
void multilevel_string_QbarQ(Gauge_Conf * GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             int t_start,
                             int dt)
  {
  const int numplaqs=(STDIM*(STDIM-1))/2;
  int i, upd;
  long int r;
  int level;

  level=-2;
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
    if(level==-2)
      {
      fprintf(stderr, "Error in the determination of the level in the multilevel (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }

  switch(level)
    {
    case -1 :     // LEVEL -1, do not average
      // initialyze ml_polycorr_ris[0] and ml_polyplaq_ris[0] to 1
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         one_TensProd(&(GC->ml_polycorr_ris[0][r]));
         for(i=0; i<numplaqs; i++)
            {
            one_TensProd(&(GC->ml_polyplaq_ris[0][r][i]));
            }
         }

      // call lower levels
      for(i=0; i<(param->d_size[0])/(param->d_ml_step[0]); i++)
         {
         multilevel_string_QbarQ(GC,
                                 geo,
                                 param,
                                 t_start+i*param->d_ml_step[0],
                                 param->d_ml_step[0]);
         }
      break;
      // end of the outermost level

    case NLEVELS-1 : // INNERMOST LEVEL
      // initialize ml_polycorr_tmp[level] and ml_polyplaq_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
         for(i=0; i<numplaqs; i++)
            {
            zero_TensProd(&(GC->ml_polyplaq_tmp[level][r][i]));
            }
         }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         slice_single_update(GC,
                             geo,
                             param,
                             t_start,
                             dt);

         // compute Polyakov loop and plaquettes restricted to the slice
         GAUGE_GROUP *loc_poly = (GAUGE_GROUP *) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(GAUGE_GROUP));

         double **loc_plaq=(double**) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(double*));
         for(r=0; r<param->d_space_vol; r++)
            {
            loc_plaq[r] = (double *) mymalloc(DOUBLE_ALIGN, (unsigned long) numplaqs * sizeof(double) );
            }

         compute_local_poly(GC,
                            geo,
                            param,
                            t_start,
                            dt,
                            loc_poly);

         if((t_start==0 && param->d_ml_step[NLEVELS-1]>1) || (t_start==1 && param->d_ml_step[NLEVELS-1]==1))
           {
           compute_plaq_on_slice1(GC, geo, param, loc_plaq);
           }

         // compute the tensor products
         // and update ml_polycorr_tmp[level], ml_polyplaq_tmp[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r, i)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            TensProd TP, TP2;
            long r1, r2;
            int t_tmp;

            r1=sisp_and_t_to_si(r, 0, param);
            for(i=0; i<param->d_dist_poly; i++) r1=nnp(geo, r1, 1);
            si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1

            TensProd_init(&TP, &(loc_poly[r]), &(loc_poly[r2]) );

            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &TP);

            if((t_start==0 && param->d_ml_step[NLEVELS-1]>1) || (t_start==1 && param->d_ml_step[NLEVELS-1]==1))
              {
              r1=sisp_and_t_to_si(r, 0, param);
              for(i=0; i<param->d_dist_poly/2; i++) r1=nnp(geo, r1, 1);
              for(i=0; i<param->d_trasv_dist; i++) r1=nnp(geo, r1, 2);
              si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1

              for(i=0; i<numplaqs; i++)
                 {
                 equal_TensProd(&TP2, &TP);
                 times_equal_real_TensProd(&TP2, loc_plaq[r2][i]);

                 plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r][i]), &TP2);
                 }
              }
            else
              {
              for(i=0; i<numplaqs; i++)
                 {
                 plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r][i]), &TP);
                 }
              }
            }

        free(loc_poly);
        for(r=0; r<param->d_space_vol; r++)
           {
           free(loc_plaq[r]);
           }
        free(loc_plaq);
        } // end of update

      // normalize polycorr_tmp[level] and polyplaq_tmp[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         for(i=0; i<numplaqs; i++)
            {
            times_equal_real_TensProd(&(GC->ml_polyplaq_tmp[level][r][i]), 1.0/(double) param->d_ml_upd[level]);
            }
         }

      // update polycorr_ris[level] and polyplaq_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         for(i=0; i<numplaqs; i++)
            {
            times_equal_TensProd(&(GC->ml_polyplaq_ris[level][r][i]), &(GC->ml_polyplaq_tmp[level][r][i]));
            }
         }

      break;
      // end of innermost level

    default:  // NOT THE INNERMOST NOT THE OUTERMOST LEVEL
      if(level==-1 || level==NLEVELS-1)
        {
        fprintf(stderr, "Error in the multilevel (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }

      // initialize ml_polycorr_tmp[level] and ml_polyplaq_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
         for(i=0; i<numplaqs; i++)
            {
            zero_TensProd(&(GC->ml_polyplaq_tmp[level][r][i]));
            }
         }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         slice_single_update(GC,
                             geo,
                             param,
                             t_start,
                             dt);

         // initialyze ml_polycorr_ris[level+1] and ml_polyplaq_ris[level+1] to 1
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r, i)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            one_TensProd(&(GC->ml_polycorr_ris[level+1][r]));
            for(i=0; i<numplaqs; i++)
               {
               one_TensProd(&(GC->ml_polyplaq_ris[level+1][r][i]));
               }
            }

         // call higher levels
         for(i=0; i<(param->d_ml_step[level])/(param->d_ml_step[level+1]); i++)
            {
            multilevel_string_QbarQ(GC,
                                    geo,
                                    param,
                                    t_start+i*param->d_ml_step[level+1],
                                    param->d_ml_step[level+1]);
            }

         // update polycorr_tmp[level] with polycorr_ris[level+1]
         // and analously for polyplaq_tmp
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r, i)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &(GC->ml_polycorr_ris[level+1][r]));
            for(i=0; i<numplaqs; i++)
               {
               plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r][i]), &(GC->ml_polyplaq_ris[level+1][r][i]));
               }
            }

         } // end of update

      // normalize polycorr_tmp[level] and polyplaq_tmp[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         for(i=0; i<numplaqs; i++)
            {
            times_equal_real_TensProd(&(GC->ml_polyplaq_tmp[level][r][i]), 1.0/(double) param->d_ml_upd[level]);
            }
         }

      // update polycorr_ris[level] and polyplaq_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         for(i=0; i<numplaqs; i++)
            {
            times_equal_TensProd(&(GC->ml_polyplaq_ris[level][r][i]), &(GC->ml_polyplaq_tmp[level][r][i]));
            }
         }
      break;
      // end of the not innermost not outermost level

    } // end of switch
  } // end of multilevel



#endif

