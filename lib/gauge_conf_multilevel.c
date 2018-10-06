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
    calcstaples_wilson(GC, geo, param, r, dir, &staple);

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
void compute_local_poly(Gauge_Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        int t_start,
                        int dt)
  {
  int num_hit;
  long r;

  if(param->d_dist_poly>1 && param->d_size[1]-param->d_dist_poly>1) // Polyakov loops are separated along the "1" direction
    {
    num_hit=param->d_multihit;
    }
  else
    {
    num_hit=0;
    }

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol; r++)
     {
     int i;
     GAUGE_GROUP matrix;

     one(&(GC->loc_poly[r]));
     for(i=0; i<dt; i++)
        {
        multihit(GC,
                 geo,
                 param,
                 sisp_and_t_to_si(geo, r, t_start+i),
                 0,
                 num_hit,
                 &matrix);
        times_equal(&(GC->loc_poly[r]), &matrix);
        }
     }
  }


// compute the clovers needed to update a slice
// note that for t=t_start only temporal links are updated, thus
// we do not need clovers at t=t_start-1, while we do need those at t=t_start+dt
void slice_compute_clovers(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           int dir,
                           int t_start,
                           int dt)
  {
  GAUGE_GROUP aux;
  long r, rs;
  int t, i, j;

  for(t=t_start; t<t_start+dt; t++)
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r, rs, i, j, aux)
     #endif
     for(rs=0; rs<param->d_space_vol; rs++)
        {
        r=sisp_and_t_to_si(geo, rs, t);
        for(i=0; i<4; i++)
           {
           for(j=i+1; j<4; j++)
              {
              if(i!=dir && j!=dir)
                {
                clover(GC, geo, param, r, i, j, &aux);

                equal(&(GC->clover_array[r][i][j]), &aux);
                minus_equal_dag(&(GC->clover_array[r][i][j]), &aux);  // clover_array[r][i][j]=aux-aux^{dag}

                equal(&(GC->clover_array[r][j][i]), &(GC->clover_array[r][i][j]));
                times_equal_real(&(GC->clover_array[r][j][i]), -1.0); // clover_array[r][j][i]=-clover_array[r][i][j]
                }
              }
           }
        }
     }

  if(t_start+dt<param->d_size[0])
    {
    t=t_start+dt;
    }
  else
    {
    t=0;
    }

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r, rs, i, j, aux)
  #endif
  for(rs=0; rs<param->d_space_vol; rs++)
     {
     r=sisp_and_t_to_si(geo, rs, t);
     for(i=0; i<4; i++)
        {
        for(j=i+1; j<4; j++)
           {
           if(i!=dir && j!=dir)
             {
             clover(GC, geo, param, r, i, j, &aux);

             equal(&(GC->clover_array[r][i][j]), &aux);
             minus_equal_dag(&(GC->clover_array[r][i][j]), &aux);  // clover_array[r][i][j]=aux-aux^{dag}

             equal(&(GC->clover_array[r][j][i]), &(GC->clover_array[r][i][j]));
             times_equal_real(&(GC->clover_array[r][j][i]), -1.0); // clover_array[r][j][i]=-clover_array[r][i][j]
             }
           }
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
  #ifdef THETA_MODE
    slice_compute_clovers(GC, geo, param, 0, t_start, 1);
  #endif

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol/2; r++)
     {
     heatbath(GC, geo, param, sisp_and_t_to_si(geo, r, t_start), 0);
     }

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=param->d_space_vol/2; r<param->d_space_vol; r++)
     {
     heatbath(GC, geo, param, sisp_and_t_to_si(geo, r, t_start), 0);
     }

  for(dir=0; dir<STDIM; dir++)
     {
     #ifdef THETA_MODE
       slice_compute_clovers(GC, geo, param, dir, t_start, dt);
     #endif

     for(i=1; i<dt; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol/2; r++)
           {
           heatbath(GC, geo, param, sisp_and_t_to_si(geo, r, t_start+i), dir);
           }

        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=param->d_space_vol/2; r<param->d_space_vol; r++)
           {
           heatbath(GC, geo, param, sisp_and_t_to_si(geo, r, t_start+i), dir);
           }
        }
     }

  // overrelaxation
  for(upd_over=0; upd_over<param->d_overrelax; upd_over++)
     {
     #ifdef THETA_MODE
       slice_compute_clovers(GC, geo, param, 0, t_start, 1);
     #endif

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_space_vol/2; r++)
        {
        overrelaxation(GC, geo, param, sisp_and_t_to_si(geo, r, t_start), 0);
        }

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=param->d_space_vol/2; r<param->d_space_vol; r++)
        {
        overrelaxation(GC, geo, param, sisp_and_t_to_si(geo, r, t_start), 0);
        }

     for(dir=0; dir<STDIM; dir++)
        {
        #ifdef THETA_MODE
          slice_compute_clovers(GC, geo, param, dir, t_start, dt);
        #endif

        for(i=1; i<dt; i++)
           {
           #ifdef OPENMP_MODE
           #pragma omp parallel for num_threads(NTHREADS) private(r)
           #endif
           for(r=0; r<param->d_space_vol/2; r++)
              {
              overrelaxation(GC, geo, param, sisp_and_t_to_si(geo, r, t_start+i), dir);
              }

           #ifdef OPENMP_MODE
           #pragma omp parallel for num_threads(NTHREADS) private(r)
           #endif
           for(r=param->d_space_vol/2; r<param->d_space_vol; r++)
              {
              overrelaxation(GC, geo, param, sisp_and_t_to_si(geo, r, t_start+i), dir);
              }
           }
        }
     }

  // unitarize
  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol; r++)
     {
     unitarize( &(GC->lattice[sisp_and_t_to_si(geo, r, t_start)][0]) );
     }

  for(i=1; i<dt; i++)
     {
     for(dir=0; dir<STDIM; dir++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++) unitarize( &(GC->lattice[sisp_and_t_to_si(geo, r, t_start+i)][dir]) );
        }
     }
  }


// multilevel for polyakov correlator
void multilevel_polycorr(Gauge_Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int t_start,
                          int dt)
  {
  int i, upd;
  long int r;
  int level;

  // remember that d_ml_step[0]>d_ml_step[1]> d_ml_step[2] ...

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
         multilevel_polycorr(GC,
                              geo,
                              param,
                              t_start+i*param->d_ml_step[0],
                              param->d_ml_step[0]);
         }
      break;
      // end of the outermost level


    case NLEVELS-1 : // INNERMOST LEVEL

      // in case level -1 is never used
      if(level==0 && param->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr_ris[0] to 1
        #ifdef openmp_mode
        #pragma omp parallel for num_threads(nthreads) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           one_TensProd(&(GC->ml_polycorr_ris[0][r]));
           }
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
         slice_single_update(GC, geo, param, t_start, dt);

         // compute Polyakov loop restricted to the slice
         compute_local_poly(GC, geo, param, t_start, dt);

         // compute the tensor products
         // and update ml_polycorr_tmp[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            TensProd TP;
            long r1, r2;
            int j, t_tmp;

            r1=sisp_and_t_to_si(geo, r, 0);
            for(j=0; j<param->d_dist_poly; j++)
               {
               r1=nnp(geo, r1, 1);
               }
            si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

            TensProd_init(&TP, &(GC->loc_poly[r]), &(GC->loc_poly[r2]) );
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &TP);
            }
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

      // in case level -1 is never used
      if(level==0 && param->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr_ris[0] to 1
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           one_TensProd(&(GC->ml_polycorr_ris[0][r]));
           }
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
         slice_single_update(GC, geo, param, t_start, dt);

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
            multilevel_polycorr(GC,
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


// multilevel for polyakov orrelator to be used in long simulations
void multilevel_potlycorr_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int t_start,
                               int dt,
                               int iteration)
  {
  int i, upd;
  long int r;

  if(dt!=param->d_ml_step[0])
    {
    fprintf(stderr, "'dt' has to be equal to ml_step[0] in multilevel_polycorr_long (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

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

  // perform the update
  for(upd=0; upd< param->d_ml_upd[0]; upd++)
     {
     slice_single_update(GC, geo, param, t_start, dt);
     if(NLEVELS==1)
       {
       // compute Polyakov loop restricted to the slice
       compute_local_poly(GC, geo, param, t_start, dt);

       // compute the tensor products
       // and update ml_polycorr_tmp[0]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(r)
       #endif
       for(r=0; r<param->d_space_vol; r++)
          {
          TensProd TP;
          long r1, r2;
          int j, t_tmp;

          r1=sisp_and_t_to_si(geo, r, 0);
          for(j=0; j<param->d_dist_poly; j++)
             {
             r1=nnp(geo, r1, 1);
             }
          si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

          TensProd_init(&TP, &(GC->loc_poly[r]), &(GC->loc_poly[r2]) );
          plus_equal_TensProd(&(GC->ml_polycorr_tmp[0][r]), &TP);
          }
       }
     else  // NLEVELS>1
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
       for(i=0; i<(param->d_ml_step[0])/(param->d_ml_step[1]); i++)
          {
          // important: we have to call the "non long" version in inner levels
          multilevel_polycorr(GC,
                               geo,
                               param,
                               t_start+i*param->d_ml_step[1],
                               param->d_ml_step[1]);
          }

       // update polycorr_tmp[0] with polycorr_ris[1]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(r)
       #endif
       for(r=0; r<param->d_space_vol; r++)
          {
          plus_equal_TensProd(&(GC->ml_polycorr_tmp[0][r]), &(GC->ml_polycorr_ris[1][r]));
          }
       }
    } // end update

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

    // update polycorr_ris[0]
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(r)
    #endif
    for(r=0; r<param->d_space_vol; r++)
       {
       times_equal_TensProd(&(GC->ml_polycorr_ris[0][r]), &(GC->ml_polycorr_tmp[0][r]));
       }
    }
  } // end of multilevel


// compute plaquettes on the t=1 time slice
void compute_plaq_on_slice1(Gauge_Conf *GC,
                            Geometry const * const geo,
                            GParam const * const param)
  {
  long r;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol; r++)
     {
     int j;
     long r4;

     r4=sisp_and_t_to_si(geo, r, 1); // t=1

     // moves to the correct position of the plaquette
     for(j=0; j<param->d_dist_poly/2; j++)
        {
        r4=nnp(geo, r4, 1);
        }
     for(j=0; j<param->d_trasv_dist; j++)
        {
        r4=nnp(geo, r4, 2);
        }

     GC->loc_plaq[r]=plaquettep_complex(GC, geo, param, r4, param->d_plaq_dir[0], param->d_plaq_dir[1]);
     }
  }


// multilevel for polyakov string width
void multilevel_tube_disc(Gauge_Conf * GC,
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

      // initialyze ml_polycorr_ris[0] and ml_polyplaq_ris[0] to 1
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         one_TensProd(&(GC->ml_polycorr_ris[0][r]));
         one_TensProd(&(GC->ml_polyplaq_ris[0][r]));
         }

      // call lower levels
      for(i=0; i<(param->d_size[0])/(param->d_ml_step[0]); i++)
         {
         multilevel_tube_disc(GC,
                                    geo,
                                    param,
                                    t_start+i*param->d_ml_step[0],
                                    param->d_ml_step[0]);
         }
      break;
      // end of the outermost level


    case NLEVELS-1 : // INNERMOST LEVEL

      // in case level -1 is never used
      if(level==0 && param->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr_ris[0] to 1
        #ifdef openmp_mode
        #pragma omp parallel for num_threads(nthreads) private(r, i)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           one_TensProd(&(GC->ml_polycorr_ris[0][r]));
           one_TensProd(&(GC->ml_polyplaq_ris[0][r]));
           }
        }

      // initialize ml_polycorr_tmp[level] and ml_polyplaq_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
         zero_TensProd(&(GC->ml_polyplaq_tmp[level][r]));
         }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         slice_single_update(GC, geo, param, t_start, dt);

         // compute Polyakov loop and plaquettes (when needed) restricted to the slice
         compute_local_poly(GC, geo, param, t_start, dt);
         if(t_start==0)
           {
           compute_plaq_on_slice1(GC, geo, param);
           }

         // compute the tensor products
         // and update ml_polycorr_tmp[level], ml_polyplaq_tmp[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            TensProd TP, TP2;
            long r1, r2;
            int j, t_tmp;

            r1=sisp_and_t_to_si(geo, r, 0);
            for(j=0; j<param->d_dist_poly; j++)
               {
               r1=nnp(geo, r1, 1);
               }
            si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

            TensProd_init(&TP, &(GC->loc_poly[r]), &(GC->loc_poly[r2]) );

            // update ml_polycorr_tmp
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &TP);

            // update ml_polyplaq_tmp when needed
            if(t_start==0)
              {
              equal_TensProd(&TP2, &TP);
              times_equal_complex_TensProd(&TP2, GC->loc_plaq[r]); // the plaquette is computed in such a way that
                                                               // here just "r" is needed
              plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r]), &TP2);
              }
            else
              {
              plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r]), &TP);
              }
            }
         } // end of the for on "upd"

      // normalize polycorr_tmp[level] and polyplaq_tmp[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         times_equal_real_TensProd(&(GC->ml_polyplaq_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      // update polycorr_ris[level] and polyplaq_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         times_equal_TensProd(&(GC->ml_polyplaq_ris[level][r]), &(GC->ml_polyplaq_tmp[level][r]));
         }

      break;
      // end of innermost level


    default:  // NOT THE INNERMOST NOT THE OUTERMOST LEVEL

      if(level==-1 || level==NLEVELS-1)
        {
        fprintf(stderr, "Error in the multilevel (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }

      // in case level -1 is never used
      if(level==0 && param->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr_ris[0] to 1
        #ifdef openmp_mode
        #pragma omp parallel for num_threads(nthreads) private(r, i)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           one_TensProd(&(GC->ml_polycorr_ris[0][r]));
           one_TensProd(&(GC->ml_polyplaq_ris[0][r]));
           }
        }

      // initialize ml_polycorr_tmp[level] and ml_polyplaq_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
         zero_TensProd(&(GC->ml_polyplaq_tmp[level][r]));
         }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         slice_single_update(GC, geo, param, t_start, dt);

         // initialyze ml_polycorr_ris[level+1] and ml_polyplaq_ris[level+1] to 1
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r, i)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            one_TensProd(&(GC->ml_polycorr_ris[level+1][r]));
            one_TensProd(&(GC->ml_polyplaq_ris[level+1][r]));
            }

         // call higher levels
         for(i=0; i<(param->d_ml_step[level])/(param->d_ml_step[level+1]); i++)
            {
            multilevel_tube_disc(GC,
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
            plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r]), &(GC->ml_polyplaq_ris[level+1][r]));
            }

         } // end of update

      // normalize polycorr_tmp[level] and polyplaq_tmp[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         times_equal_real_TensProd(&(GC->ml_polyplaq_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      // update polycorr_ris[level] and polyplaq_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r, i)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         times_equal_TensProd(&(GC->ml_polyplaq_ris[level][r]), &(GC->ml_polyplaq_tmp[level][r]));
         }
      break;
      // end of the not innermost not outermost level

    } // end of switch
  } // end of multilevel


// multilevel for polyakov string width to be used in long simulations
void multilevel_tube_disc_long(Gauge_Conf * GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     int t_start,
                                     int dt,
                                     int iteration)
  {
  int i, upd;
  long int r;

  if(dt!=param->d_ml_step[0])
    {
    fprintf(stderr, "'dt' has to be equal to ml_step[0] in multilevel_polycorr_long (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize ml_polycorr_ris[0] and ml_polyplaq_ris[0] to 1 if needed
  if(t_start==0 && iteration==0)
    {
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(r, i)
    #endif
    for(r=0; r<param->d_space_vol; r++)
       {
       one_TensProd(&(GC->ml_polycorr_ris[0][r]));
       one_TensProd(&(GC->ml_polyplaq_ris[0][r]));
       }
    }

  if(iteration==0)
    {
    // initialize ml_polycorr_tmp[0] and ml_polyplaq_tmp[0] to 0 when needed
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(r, i)
    #endif
    for(r=0; r<param->d_space_vol; r++)
       {
       zero_TensProd(&(GC->ml_polycorr_tmp[0][r]));
       zero_TensProd(&(GC->ml_polyplaq_tmp[0][r]));
       }
    }

  // perform the update
  for(upd=0; upd< param->d_ml_upd[0]; upd++)
     {
     slice_single_update(GC, geo, param, t_start, dt);
     if(NLEVELS==1)
       {
       // compute Polyakov loop and plaquettes (when needed) restricted to the slice
       compute_local_poly(GC, geo, param, t_start, dt);
       if(t_start==0)
         {
         compute_plaq_on_slice1(GC, geo, param);
         }

       // compute the tensor products
       // and update ml_polycorr_tmp[level], ml_polyplaq_tmp[level]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(r)
       #endif
       for(r=0; r<param->d_space_vol; r++)
          {
          TensProd TP, TP2;
          long r1, r2;
          int j, t_tmp;

          r1=sisp_and_t_to_si(geo, r, 0);
          for(j=0; j<param->d_dist_poly; j++)
             {
             r1=nnp(geo, r1, 1);
             }
          si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

          TensProd_init(&TP, &(GC->loc_poly[r]), &(GC->loc_poly[r2]) );

          // update ml_polycorr_tmp
          plus_equal_TensProd(&(GC->ml_polycorr_tmp[0][r]), &TP);

          // update ml_polyplaq_tmp when needed
          if(t_start==0)
            {
            equal_TensProd(&TP2, &TP);
            times_equal_complex_TensProd(&TP2, GC->loc_plaq[r]); // the plaquette is computed in such a way that
                                                             // here just "r" is needed
            plus_equal_TensProd(&(GC->ml_polyplaq_tmp[0][r]), &TP2);
            }
          else
            {
            plus_equal_TensProd(&(GC->ml_polyplaq_tmp[0][r]), &TP);
            }
          }
       }
     else  // NLEVELS>1
       {
       // initialyze ml_polycorr_ris[1] and ml_polyplaq_ris[1] to 1
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(r, i)
       #endif
       for(r=0; r<param->d_space_vol; r++)
          {
          one_TensProd(&(GC->ml_polycorr_ris[1][r]));
          one_TensProd(&(GC->ml_polyplaq_ris[1][r]));
          }

       // call inner levels
       for(i=0; i<(param->d_ml_step[0])/(param->d_ml_step[1]); i++)
          {
          // important: we have to call the "non long" version in inner levels
          multilevel_tube_disc(GC,
                                     geo,
                                     param,
                                     t_start+i*param->d_ml_step[1],
                                     param->d_ml_step[1]);
          }

       // update polycorr_tmp[0] and polyplaq_tmp[0] with polycorr_ris[1] and polyplaq_tmp[1]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(r, i)
       #endif
       for(r=0; r<param->d_space_vol; r++)
          {
          plus_equal_TensProd(&(GC->ml_polycorr_tmp[0][r]), &(GC->ml_polycorr_ris[1][r]));
          plus_equal_TensProd(&(GC->ml_polyplaq_tmp[0][r]), &(GC->ml_polyplaq_ris[1][r]));
          }
       }
    } // end update

  if(iteration==param->d_ml_level0_repeat-1) // iteration starts from zero
    {
    // normalize polycorr_tmp[0] and polyplaq_tmp[0]
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(r, i)
    #endif
    for(r=0; r<param->d_space_vol; r++)
       {
       times_equal_real_TensProd(&(GC->ml_polycorr_tmp[0][r]), 1.0/( (double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat));
       times_equal_real_TensProd(&(GC->ml_polyplaq_tmp[0][r]), 1.0/( (double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat));
       }

    // update polycorr_ris[0] and polyplaq_ris[0]
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(r, i)
    #endif
    for(r=0; r<param->d_space_vol; r++)
       {
       times_equal_TensProd(&(GC->ml_polycorr_ris[0][r]), &(GC->ml_polycorr_tmp[0][r]));
       times_equal_TensProd(&(GC->ml_polyplaq_ris[0][r]), &(GC->ml_polyplaq_ris[0][r]));
       }
    }
  } // end of multilevel


// compute auxlilliary stuff on the first time slice for connected correlator
void compute_poly_polyplaq_and_polyplaqconn_on_slice1(Gauge_Conf *GC,
                                                      Geometry const * const geo,
                                                      GParam const * const param)
  {
  int num_hit;
  long r;

  if(param->d_dist_poly>1 && param->d_size[1]-param->d_dist_poly>1) // Polyakov loops are separated along the "1" direction
    {
    num_hit=param->d_multihit;
    }
  else
    {
    num_hit=0;
    }

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_space_vol; r++)
     {
     int j;
     long r4;
     GAUGE_GROUP matrix;

     #ifdef DEBUG
       long r4start;
     #endif

     r4=sisp_and_t_to_si(geo, r, 0); // t=0 starting point

     multihit(GC, geo, param, r4, 0, num_hit, &matrix);
     equal(&(GC->loc_poly[r]), &matrix);
     equal(&(GC->loc_polyplaqconn[r]), &matrix);

     r4=nnp(geo, r4, 0); // now we are in the t=1 plane
     #ifdef DEBUG
       r4start=r4;
     #endif

     if(param->d_dist_poly/2>=1)
       {
       times_equal(&(GC->loc_polyplaqconn[r]), &GC->lattice[r4][1]);
       r4=nnp(geo, r4, 1);

       for(j=1; j<param->d_dist_poly/2; j++)
          {
          multihit(GC, geo, param, r4, 1, num_hit, &matrix);
          times_equal(&(GC->loc_polyplaqconn[r]), &matrix);
          r4=nnp(geo, r4, 1);
          }
       }

     if(param->d_trasv_dist>=1)
       {
       times_equal(&(GC->loc_polyplaqconn[r]), &GC->lattice[r4][2]);
       r4=nnp(geo, r4, 2);

       for(j=1; j<param->d_trasv_dist; j++)
          {
          multihit(GC, geo, param, r4, 2, num_hit, &matrix);
          times_equal(&(GC->loc_polyplaqconn[r]), &matrix);
          r4=nnp(geo, r4, 2);
          }
       }

     plaquettep_matrix(GC, geo, param, r4, param->d_plaq_dir[0], param->d_plaq_dir[1], &matrix);

     GC->loc_plaq[r]=retr(&matrix)+I*imtr(&matrix);
     times_equal(&(GC->loc_polyplaqconn[r]), &matrix);

     if(param->d_trasv_dist>=1)
       {
       for(j=1; j<param->d_trasv_dist; j++)
          {
          r4=nnm(geo, r4, 2);
          multihit(GC, geo, param, r4, 2, num_hit, &matrix);
          times_equal_dag(&(GC->loc_polyplaqconn[r]), &matrix);
          }

       r4=nnm(geo, r4, 2);
       times_equal_dag(&(GC->loc_polyplaqconn[r]), &GC->lattice[r4][2]);
       }

     if(param->d_dist_poly/2>=1)
       {
       for(j=1; j<param->d_dist_poly/2; j++)
          {
          r4=nnm(geo, r4, 1);
          multihit(GC, geo, param, r4, 1, num_hit, &matrix);
          times_equal_dag(&(GC->loc_polyplaqconn[r]), &matrix);
          }

       r4=nnm(geo, r4, 1);
       times_equal_dag(&(GC->loc_polyplaqconn[r]), &GC->lattice[r4][1]);
       }

     #ifdef DEBUG
     if(r4start!=r4)
       {
       fprintf(stderr, "Loop not closing! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     #endif

     for(j=1; j<param->d_ml_step[NLEVELS-1]; j++)
        {
        multihit(GC, geo, param, r4, 0, num_hit, &matrix);
        times_equal(&(GC->loc_polyplaqconn[r]), &matrix);
        times_equal(&(GC->loc_poly[r]), &matrix);
        r4=nnp(geo, r4, 1);
        }
     }
  }


// multilevel for string width using the connected operator
void multilevel_tube_conn(Gauge_Conf * GC,
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

      // initialyze ml_polycorr_ris[0], ml_polyplaq_ris[0] and ml_polyplaqconn_ris[0] to 1
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         one_TensProd(&(GC->ml_polycorr_ris[0][r]));
         one_TensProd(&(GC->ml_polyplaq_ris[0][r]));
         one_TensProd(&(GC->ml_polyplaqconn_ris[0][r]));
         }

      // call lower levels
      for(i=0; i<(param->d_size[0])/(param->d_ml_step[0]); i++)
         {
         multilevel_tube_conn(GC,
                              geo,
                              param,
                              t_start+i*param->d_ml_step[0],
                              param->d_ml_step[0]);
         }
      break;
      // end of the outermost level


    case NLEVELS-1 : // INNERMOST LEVEL

      // in case level -1 is never used
      if(level==0 && param->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr_ris[0], ml_polyplaq_ris[0] and ml_polyplaqconn_ris[0] to 1
        #ifdef openmp_mode
        #pragma omp parallel for num_threads(nthreads) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           one_TensProd(&(GC->ml_polycorr_ris[0][r]));
           one_TensProd(&(GC->ml_polyplaq_ris[0][r]));
           one_TensProd(&(GC->ml_polyplaqconn_ris[0][r]));
           }
        }

      // initialize ml_polycorr_tmp[level], ml_polyplaq_tmp[level] and ml_polyplaqconn_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
         zero_TensProd(&(GC->ml_polyplaq_tmp[level][r]));
         zero_TensProd(&(GC->ml_polyplaqconn_tmp[level][r]));
         }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         slice_single_update(GC, geo, param, t_start, dt);

         // compute Polyakov loop and plaquettes (when needed) restricted to the slice
         if(t_start==0)
           {
           compute_poly_polyplaq_and_polyplaqconn_on_slice1(GC, geo, param);
           }
         else
           {
           compute_local_poly(GC, geo, param, t_start, dt);
           }

         // compute the tensor products
         // and update ml_polycorr_tmp[level], ml_polyplaq_tmp[level] and ml_polyplaqconn_tmp[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            TensProd TP, TP2;
            long r1, r2;
            int j, t_tmp;

            r1=sisp_and_t_to_si(geo, r, 0);
            for(j=0; j<param->d_dist_poly; j++)
               {
               r1=nnp(geo, r1, 1);
               }
            si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

            TensProd_init(&TP, &(GC->loc_poly[r]), &(GC->loc_poly[r2]) );

            // update ml_polycorr_tmp
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &TP);

            // update ml_polyplaq_tmp and ml_polyplaqconn_tmp when needed
            if(t_start==0)
              {
              equal_TensProd(&TP2, &TP);
              times_equal_complex_TensProd(&TP2, GC->loc_plaq[r]); // the plaquette is computed in such a way that
                                                               // here just "r" is needed
              plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r]), &TP2);

              TensProd_init(&TP, &(GC->loc_polyplaqconn[r]), &(GC->loc_poly[r2]) );
              plus_equal_TensProd(&(GC->ml_polyplaqconn_tmp[level][r]), &TP);
              }
            else
              {
              plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r]), &TP);
              plus_equal_TensProd(&(GC->ml_polyplaqconn_tmp[level][r]), &TP);
              }
            }
         } // end of the for on "upd"

      // normalize polycorr_tmp[level], polyplaq_tmp[level] and ml_polyplaqconn_tmp[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         times_equal_real_TensProd(&(GC->ml_polyplaq_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         times_equal_real_TensProd(&(GC->ml_polyplaqconn_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      // update polycorr_ris[level], polyplaq_ris[level] and ml_polyplaqconn_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         times_equal_TensProd(&(GC->ml_polyplaq_ris[level][r]), &(GC->ml_polyplaq_tmp[level][r]));
         times_equal_TensProd(&(GC->ml_polyplaqconn_ris[level][r]), &(GC->ml_polyplaqconn_tmp[level][r]));
         }

      break;
      // end of innermost level


    default:  // NOT THE INNERMOST NOT THE OUTERMOST LEVEL

      if(level==-1 || level==NLEVELS-1)
        {
        fprintf(stderr, "Error in the multilevel (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }

      // in case level -1 is never used
      if(level==0 && param->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr_ris[0], ml_polyplaq_ris[0] and ml_polyplaqconn_ris[0] to 1
        #ifdef openmp_mode
        #pragma omp parallel for num_threads(nthreads) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           one_TensProd(&(GC->ml_polycorr_ris[0][r]));
           one_TensProd(&(GC->ml_polyplaq_ris[0][r]));
           one_TensProd(&(GC->ml_polyplaqconn_ris[0][r]));
           }
        }

      // initialize ml_polycorr_tmp[level], ml_polyplaq_tmp[level] and ml_polyplaqconn_tmp[level] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         zero_TensProd(&(GC->ml_polycorr_tmp[level][r]));
         zero_TensProd(&(GC->ml_polyplaq_tmp[level][r]));
         zero_TensProd(&(GC->ml_polyplaqconn_tmp[level][r]));
         }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         slice_single_update(GC, geo, param, t_start, dt);

         // initialyze ml_polycorr_ris[level+1], ml_polyplaq_ris[level+1] and ml_polyplaqconn_ris[level+1] to 1
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            one_TensProd(&(GC->ml_polycorr_ris[level+1][r]));
            one_TensProd(&(GC->ml_polyplaq_ris[level+1][r]));
            one_TensProd(&(GC->ml_polyplaqconn_ris[level+1][r]));
            }

         // call higher levels
         for(i=0; i<(param->d_ml_step[level])/(param->d_ml_step[level+1]); i++)
            {
            multilevel_tube_disc(GC,
                                       geo,
                                       param,
                                       t_start+i*param->d_ml_step[level+1],
                                       param->d_ml_step[level+1]);
            }

         // update polycorr_tmp[level] with polycorr_ris[level+1]
         // and analously for polyplaq_tmp and ml_polyplaqconn_ris
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<param->d_space_vol; r++)
            {
            plus_equal_TensProd(&(GC->ml_polycorr_tmp[level][r]), &(GC->ml_polycorr_ris[level+1][r]));
            plus_equal_TensProd(&(GC->ml_polyplaq_tmp[level][r]), &(GC->ml_polyplaq_ris[level+1][r]));
            plus_equal_TensProd(&(GC->ml_polyplaqconn_tmp[level][r]), &(GC->ml_polyplaqconn_ris[level+1][r]));
            }

         } // end of update

      // normalize polycorr_tmp[level], polyplaq_tmp[level] and ml_polyplaqconn_tmp[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_real_TensProd(&(GC->ml_polycorr_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         times_equal_real_TensProd(&(GC->ml_polyplaq_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         times_equal_real_TensProd(&(GC->ml_polyplaqconn_tmp[level][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      // update polycorr_ris[level], polyplaq_ris[level] and ml_polyplaqconn_ris[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr_ris[level][r]), &(GC->ml_polycorr_tmp[level][r]));
         times_equal_TensProd(&(GC->ml_polyplaq_ris[level][r]), &(GC->ml_polyplaq_tmp[level][r]));
         times_equal_TensProd(&(GC->ml_polyplaqconn_ris[level][r]), &(GC->ml_polyplaqconn_tmp[level][r]));
         }
      break;
      // end of the not innermost not outermost level

    } // end of switch
  } // end of multilevel



#endif

