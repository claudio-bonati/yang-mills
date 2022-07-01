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
    GAUGE_GROUP staple, staple2, partial;

    zero(G);
    equal(&partial, &(GC->lattice[r][dir]));

    // gauge part of the staple
    if(fabs(param->d_beta)>MIN_VALUE)
      {
      #ifndef THETA_MODE
        calcstaples_wilson(GC, geo, r, dir, &staple);
      #else
        // compute_clovers in direction "dir" HAS TO BE CALLED BEFORE!
        calcstaples_with_topo(GC, geo, param, r, dir, &staple);
      #endif
      times_equal_real(&staple, param->d_beta);
      }
    else
      {
      zero(&staple);
      }

    // higgs part of the staple
    if(fabs(param->d_higgs_beta)>MIN_VALUE)
      {
      vector_tensor_vector_vecs(&staple2, &(GC->higgs[r]), &(GC->higgs[nnp(geo, r, dir)]));
      times_equal_real(&staple2, param->d_higgs_beta*NHIGGS*NCOLOR); // NCOLOR is needed to compensate the
                                                                     // 1/NCOLOR that is present in the gauge part
      }
    else
      {
      zero(&staple2);
      }

    plus_equal(&staple, &staple2);

    for(i=0; i<num_hit; i++)
       {
       single_heatbath(&partial, &staple);
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


// perform a complete update on the given level
void update_for_multilevel(Gauge_Conf * GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           int level)
   {
   for(int i=0; i<STDIM; i++)
      {
      if(geo->d_size[i]==1)
        {
        fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   long r;
   int j, dir;

   // heatbath
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, dir);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(geo->d_volume)/2; r++)
         {
         int t;
         long rsp;

         si_to_sisp_and_t(&rsp, &t, geo, r);
         if((t % param->d_ml_step[level])!=0 || dir==0)
           {
           heatbath(GC, geo, param, r, dir);
           }
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(geo->d_volume)/2; r<(geo->d_volume); r++)
         {
         int t;
         long rsp;

         si_to_sisp_and_t(&rsp, &t, geo, r);
         if((t % param->d_ml_step[level])!=0 || dir==0)
           {
           heatbath(GC, geo, param, r, dir);
           }
         }
      }

   // overrelax
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, dir);
      #endif

      for(j=0; j<param->d_overrelax; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<(geo->d_volume)/2; r++)
            {
            int t;
            long rsp;

            si_to_sisp_and_t(&rsp, &t, geo, r);
            if((t % param->d_ml_step[level])!=0 || dir==0)
              {
              overrelaxation(GC, geo, param, r, dir);
              }
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=(geo->d_volume)/2; r<(geo->d_volume); r++)
            {
            int t;
            long rsp;

            si_to_sisp_and_t(&rsp, &t, geo, r);
            if((t % param->d_ml_step[level])!=0 || dir==0)
              {
              overrelaxation(GC, geo, param, r, dir);
              }
            }
         }
      }

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }
   }


// perform a complete update of gauge and higgs on the given level
void update_for_multilevel_with_higgs(Gauge_Conf * GC,
                                      Geometry const * const geo,
                                      GParam const * const param,
                                      int level)
   {
   for(int i=0; i<STDIM; i++)
      {
      if(geo->d_size[i]==1)
        {
        fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   long r;
   int j, dir;

   // heatbath
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, dir);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(geo->d_volume)/2; r++)
         {
         int t, acc;
         long rsp;

         si_to_sisp_and_t(&rsp, &t, geo, r);
         if((t % param->d_ml_step[level])!=0 || dir==0)
           {
           heatbath_with_higgs(GC, geo, param, r, dir);
           }
         if((t % param->d_ml_step[level])!=0)
           {
           acc=metropolis_for_higgs(GC, geo, param, r);
           (void) acc; // just to avoid warning at compile time
           }
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(geo->d_volume)/2; r<(geo->d_volume); r++)
         {
         int t, acc;
         long rsp;

         si_to_sisp_and_t(&rsp, &t, geo, r);
         if((t % param->d_ml_step[level])!=0 || dir==0)
           {
           heatbath_with_higgs(GC, geo, param, r, dir);
           }
         if((t % param->d_ml_step[level])!=0)
           {
           acc=metropolis_for_higgs(GC, geo, param, r);
           (void) acc; // just to avoid warning at compile time
           }
         }
      }

   // overrelax
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, dir);
      #endif

      for(j=0; j<param->d_overrelax; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<(geo->d_volume)/2; r++)
            {
            int t;
            long rsp;

            si_to_sisp_and_t(&rsp, &t, geo, r);
            if((t % param->d_ml_step[level])!=0 || dir==0)
              {
              overrelaxation_with_higgs(GC, geo, param, r, dir);
              }
            if((t % param->d_ml_step[level])!=0)
              {
              overrelaxation_for_higgs(GC, geo, r);
              }
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=(geo->d_volume)/2; r<(geo->d_volume); r++)
            {
            int t;
            long rsp;

            si_to_sisp_and_t(&rsp, &t, geo, r);
            if((t % param->d_ml_step[level])!=0 || dir==0)
              {
              overrelaxation_with_higgs(GC, geo, param, r, dir);
              }
            if((t % param->d_ml_step[level])!=0)
              {
              overrelaxation_for_higgs(GC, geo, r);
              }
            }
         }
      }

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(geo->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         normalize_vecs(&(GC->higgs[r]));
         }
      }
   }




// compute polyakov loop on a single slice
void compute_local_poly(Gauge_Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param)
  {
  int num_hit;
  long raux;

  if(param->d_dist_poly >1 && (geo->d_size[1]-param->d_dist_poly) >1) // Polyakov loops are separated along the "1" direction
    {
    num_hit=param->d_multihit;
    }
  else
    {
    num_hit=0;
    }

  #ifdef THETA_MODE
  // clovers are eventually needed by the multihit
  compute_clovers(GC, geo, 0);
  #endif

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(raux)
  #endif
  for(raux=0; raux<(geo->d_space_vol*geo->d_size[0]/param->d_ml_step[NLEVELS-1]); raux++)
     {
     int i, t;
     GAUGE_GROUP matrix;

     const long r = raux/(geo->d_size[0]/param->d_ml_step[NLEVELS-1]);
     const int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[NLEVELS-1]) );

     one(&(GC->loc_poly[slice][r]));
     for(i=0; i<param->d_ml_step[NLEVELS-1]; i++)
        {
        t=slice*(param->d_ml_step[NLEVELS-1])+i;
        multihit(GC,
                 geo,
                 param,
                 sisp_and_t_to_si(geo, r, t),
                 0,
                 num_hit,
                 &matrix);
        times_equal(&(GC->loc_poly[slice][r]), &matrix);
        }
     }
  }


// multilevel for polyakov correlator
void multilevel_polycorr(Gauge_Conf * GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         int dt)
  {
  int upd, level;
  long raux;

  // remember that d_size[0] >= d_ml_step[0]>d_ml_step[1]> d_ml_step[2] ...

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

      // initialyze ml_polycorr[0] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

         zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
         }

      // call lower levels
      multilevel_polycorr(GC,
                          geo,
                          param,
                          param->d_ml_step[0]);

      break;
      // end of the outermost level

    case NLEVELS-1 : // INNERMOST LEVEL

      // in case level -1 is never used
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel(GC, geo, param, level);

         // compute Polyakov loop restricted to the slice
         compute_local_poly(GC, geo, param);

         // compute the tensor products
         // and update ml_polycorr[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            TensProd TP;
            long r1, r2;
            int j, t_tmp;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            r1=sisp_and_t_to_si(geo, r, 0); // r is a 3d index, r1 is the 4d index value of (r,t=0)
            for(j=0; j<param->d_dist_poly; j++)
               {
               r1=nnp(geo, r1, 1);
               }
            si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

            TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);
            }

        } // end of update

      // normalize polycorr
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
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
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel(GC, geo, param, level);

         // initialyze ml_polycorr[level+1] to 0
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level+1]; raux++)
            {
            long r = raux/(geo->d_size[0]/param->d_ml_step[level+1]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level+1]) );

            zero_TensProd(&(GC->ml_polycorr[level+1][slice][r]));
            }

         // call higher levels
         multilevel_polycorr(GC,
                             geo,
                             param,
                             param->d_ml_step[level+1]);

         // update polycorr[level] with polycorr[level+1]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            int j;
            TensProd TP;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            one_TensProd(&TP);
            for(j=0; j<param->d_ml_step[level]/param->d_ml_step[level+1]; j++)
               {
               times_equal_TensProd(&TP, &(GC->ml_polycorr[level+1][slice*param->d_ml_step[level]/param->d_ml_step[level+1]+j][r]));
               }

            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);
            }

         } // end of update

      // normalize polycorr[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      break;
      // end of the not innermost not outermost level

    } // end of switch

  } // end of multilevel


// multilevel for polyakov correlator with higgs
void multilevel_polycorr_with_higgs(Gauge_Conf * GC,
                                    Geometry const * const geo,
                                    GParam const * const param,
                                    int dt)
  {
  int upd, level;
  long raux;

  // remember that d_size[0] >= d_ml_step[0]>d_ml_step[1]> d_ml_step[2] ...

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

      // initialyze ml_polycorr[0] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

         zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
         }

      // call lower levels
      multilevel_polycorr_with_higgs(GC,
                                     geo,
                                     param,
                                     param->d_ml_step[0]);

      break;
      // end of the outermost level

    case NLEVELS-1 : // INNERMOST LEVEL

      // in case level -1 is never used
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel_with_higgs(GC, geo, param, level);

         // compute Polyakov loop restricted to the slice
         compute_local_poly(GC, geo, param);

         // compute the tensor products
         // and update ml_polycorr[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            TensProd TP;
            long r1, r2;
            int j, t_tmp;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            r1=sisp_and_t_to_si(geo, r, 0); // r is a 3d index, r1 is the 4d index value of (r,t=0)
            for(j=0; j<param->d_dist_poly; j++)
               {
               r1=nnp(geo, r1, 1);
               }
            si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

            TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);
            }

        } // end of update

      // normalize polycorr
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
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
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel_with_higgs(GC, geo, param, level);

         // initialyze ml_polycorr[level+1] to 0
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level+1]; raux++)
            {
            long r = raux/(geo->d_size[0]/param->d_ml_step[level+1]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level+1]) );

            zero_TensProd(&(GC->ml_polycorr[level+1][slice][r]));
            }

         // call higher levels
         multilevel_polycorr_with_higgs(GC,
                                        geo,
                                        param,
                                        param->d_ml_step[level+1]);

         // update polycorr[level] with polycorr[level+1]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            int j;
            TensProd TP;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            one_TensProd(&TP);
            for(j=0; j<param->d_ml_step[level]/param->d_ml_step[level+1]; j++)
               {
               times_equal_TensProd(&TP, &(GC->ml_polycorr[level+1][slice*param->d_ml_step[level]/param->d_ml_step[level+1]+j][r]));
               }

            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);
            }

         } // end of update

      // normalize polycorr[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
         }

      break;
      // end of the not innermost not outermost level

    } // end of switch

  } // end of multilevel


// multilevel for polyakov correlator to be used in long simulations
void multilevel_polycorr_long(Gauge_Conf * GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              int dt,
                              int iteration)
  {
  int upd;
  long raux;

  if(dt!=param->d_ml_step[0])
    {
    fprintf(stderr, "'dt' has to be equal to ml_step[0] in multilevel_polycorr_long (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize ml_polycorr[0] to 0 if needed
  if(iteration==0)
    {
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
       }
    }

  // perform the update
  for(upd=0; upd< param->d_ml_upd[0]; upd++)
     {
     // update on level zero
     update_for_multilevel(GC, geo, param, 0);

     #if NLEVELS==1
       // compute Polyakov loop restricted to the slice
       compute_local_poly(GC, geo, param);

       // compute the tensor products
       // and update ml_polycorr[0]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          TensProd TP;
          long r1, r2;
          int j, t_tmp;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          r1=sisp_and_t_to_si(geo, r, 0);
          for(j=0; j<param->d_dist_poly; j++)
             {
             r1=nnp(geo, r1, 1);
             }
          si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

          TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);
          }
     #else  // NLEVELS>1
       // initialyze ml_polycorr[1] to zero
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[1]; raux++)
          {
          long r = raux/(geo->d_size[0]/param->d_ml_step[1]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[1]) );

          zero_TensProd(&(GC->ml_polycorr[1][slice][r]));
          }

       // call inner levels
       // important: we have to call the "non long" version in inner levels
       multilevel_polycorr(GC,
                           geo,
                           param,
                           param->d_ml_step[1]);

       // update polycorr[0] with polycorr[1]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          int j;
          TensProd TP;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          one_TensProd(&TP);
          for(j=0; j<param->d_ml_step[0]/param->d_ml_step[1]; j++)
             {
             times_equal_TensProd(&TP, &(GC->ml_polycorr[1][slice*param->d_ml_step[0]/param->d_ml_step[1]+j][r]));
             }

          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);
          }
     #endif
     } // end update

  if(iteration==param->d_ml_level0_repeat-1) // iteration starts from zero
    {
    // normalize polycorr[0]
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       times_equal_real_TensProd(&(GC->ml_polycorr[0][slice][r]), 1.0/((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat) );
       }
    }
  } // end of multilevel


// multilevel for polyakov correlator to be used in long simulations with higgs fields
void multilevel_polycorr_long_with_higgs(Gauge_Conf * GC,
                                         Geometry const * const geo,
                                         GParam const * const param,
                                         int dt,
                                         int iteration)
  {
  int upd;
  long raux;

  if(dt!=param->d_ml_step[0])
    {
    fprintf(stderr, "'dt' has to be equal to ml_step[0] in multilevel_polycorr_long (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize ml_polycorr[0] to 0 if needed
  if(iteration==0)
    {
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
       }
    }

  // perform the update
  for(upd=0; upd< param->d_ml_upd[0]; upd++)
     {
     // update on level zero
     update_for_multilevel_with_higgs(GC, geo, param, 0);

     #if NLEVELS==1
       // compute Polyakov loop restricted to the slice
       compute_local_poly(GC, geo, param);

       // compute the tensor products
       // and update ml_polycorr[0]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          TensProd TP;
          long r1, r2;
          int j, t_tmp;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          r1=sisp_and_t_to_si(geo, r, 0);
          for(j=0; j<param->d_dist_poly; j++)
             {
             r1=nnp(geo, r1, 1);
             }
          si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

          TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);
          }
     #else  // NLEVELS>1
       // initialyze ml_polycorr[1] to zero
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[1]; raux++)
          {
          long r = raux/(geo->d_size[0]/param->d_ml_step[1]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[1]) );

          zero_TensProd(&(GC->ml_polycorr[1][slice][r]));
          }

       // call inner levels
       // important: we have to call the "non long" version in inner levels
       multilevel_polycorr_with_higgs(GC,
                                      geo,
                                      param,
                                      param->d_ml_step[1]);

       // update polycorr[0] with polycorr[1]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          int j;
          TensProd TP;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          one_TensProd(&TP);
          for(j=0; j<param->d_ml_step[0]/param->d_ml_step[1]; j++)
             {
             times_equal_TensProd(&TP, &(GC->ml_polycorr[1][slice*param->d_ml_step[0]/param->d_ml_step[1]+j][r]));
             }

          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);
          }
     #endif
     } // end update

  if(iteration==param->d_ml_level0_repeat-1) // iteration starts from zero
    {
    // normalize polycorr[0]
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       times_equal_real_TensProd(&(GC->ml_polycorr[0][slice][r]), 1.0/((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat) );
       }
    }
  } // end of multilevel


// compute polyakov loop and plaquette on a single slice
void compute_local_poly_and_plaq(Gauge_Conf *GC,
                                 Geometry const * const geo,
                                 GParam const * const param)
  {
  int num_hit;
  long raux;

  if(param->d_dist_poly>1 && geo->d_size[1]-param->d_dist_poly>1) // Polyakov loops are separated along the "1" direction
    {
    num_hit=param->d_multihit;
    }
  else
    {
    num_hit=0;
    }

  #ifdef THETA_MODE
  // clovers are eventually needed by the multihit
  compute_clovers(GC, geo, 0);
  #endif

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(raux)
  #endif
  for(raux=0; raux<(geo->d_space_vol*geo->d_size[0]/param->d_ml_step[NLEVELS-1]); raux++)
     {
     int i, j, t;
     long r4;
     GAUGE_GROUP matrix;

     const long r = raux/(geo->d_size[0]/param->d_ml_step[NLEVELS-1]);
     const int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[NLEVELS-1]) );

     one(&(GC->loc_poly[slice][r]));
     for(i=0; i<param->d_ml_step[NLEVELS-1]; i++)
        {
        t=slice*(param->d_ml_step[NLEVELS-1])+i;
        multihit(GC,
                 geo,
                 param,
                 sisp_and_t_to_si(geo, r, t),
                 0,
                 num_hit,
                 &matrix);
        times_equal(&(GC->loc_poly[slice][r]), &matrix);
        }

     if(slice==0)
       {
       r4=sisp_and_t_to_si(geo, r, 1); // t=1

       // moves to the correct position of the plaquette
       for(j=0; j<param->d_dist_poly/2; j++) // polyakov loop are separated along direction 1
          {
          r4=nnp(geo, r4, 1);
          }
       for(j=0; j<param->d_trasv_dist; j++) // the transverse direction is direction 2
          {
          r4=nnp(geo, r4, 2);
          }

       GC->loc_plaq[r]=plaquettep_complex(GC, geo, r4, param->d_plaq_dir[0], param->d_plaq_dir[1]);
       }
     }
  }


// multilevel for flux width computation using the disconnected correlator
void multilevel_tube_disc(Gauge_Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int dt)
  {
  int upd, level;
  long raux;

  // remember that d_size[0] >= d_ml_step[0]>d_ml_step[1]> d_ml_step[2] ...

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

      // initialyze ml_polycorr[0] and ml_polyplaq[0] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

         zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
         if(slice==0)
           {
           zero_TensProd(&(GC->ml_polyplaq[0][r]));
           }
         }

      // call lower levels
      multilevel_tube_disc(GC,
                           geo,
                           param,
                           param->d_ml_step[0]);

      break;
      // end of the outermost level

    case NLEVELS-1 : // INNERMOST LEVEL

      // in case level -1 is never used
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0] and ml_polyplaq[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           if(slice==0)
             {
             zero_TensProd(&(GC->ml_polyplaq[0][r]));
             }
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel(GC, geo, param, level);

         // compute Polyakov loop and plaquette restricted to the slice
         compute_local_poly_and_plaq(GC, geo, param);

         // compute the tensor products
         // and update ml_polycorr[level] and ml_polyplaq[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            TensProd TP;
            long r1, r2;
            int j, t_tmp;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            r1=sisp_and_t_to_si(geo, r, 0);
            for(j=0; j<param->d_dist_poly; j++)
               {
               r1=nnp(geo, r1, 1);
               }
            si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

            TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);

            if(slice==0)
              {
              times_equal_complex_TensProd(&TP, GC->loc_plaq[r]);
              plus_equal_TensProd(&(GC->ml_polyplaq[level][r]), &TP);
              }
            }

        } // end of update

      // normalize polycorr and polyplaq
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
         if(slice==0)
           {
           times_equal_real_TensProd(&(GC->ml_polyplaq[level][r]), 1.0/(double) param->d_ml_upd[level]);
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

      // in case level -1 is never used
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0] and ml_polyplaq[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           if(slice==0)
             {
             zero_TensProd(&(GC->ml_polyplaq[0][r]));
             }
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel(GC, geo, param, level);

         // initialyze ml_polycorr[level+1] and ml_polyplaq[level+1] to 0
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level+1]; raux++)
            {
            long r = raux/(geo->d_size[0]/param->d_ml_step[level+1]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level+1]) );

            zero_TensProd(&(GC->ml_polycorr[level+1][slice][r]));
            if(slice==0)
              {
              zero_TensProd(&(GC->ml_polyplaq[level+1][r]));
              }
            }

         // call higher levels
         multilevel_tube_disc(GC,
                              geo,
                              param,
                              param->d_ml_step[level+1]);

         // update polycorr[level] with polycorr[level+1]
         // and analogously for polyplaq
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            int j;
            TensProd TP;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            one_TensProd(&TP);
            for(j=0; j<param->d_ml_step[level]/param->d_ml_step[level+1]; j++)
               {
               times_equal_TensProd(&TP, &(GC->ml_polycorr[level+1][slice*param->d_ml_step[level]/param->d_ml_step[level+1]+j][r]));
               }
            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);

            if(slice==0)
              {
              equal_TensProd(&TP, &(GC->ml_polyplaq[level+1][r]));
              for(j=1; j<param->d_ml_step[level]/param->d_ml_step[level+1]; j++)
                 {
                 times_equal_TensProd(&TP, &(GC->ml_polycorr[level+1][slice*param->d_ml_step[level]/param->d_ml_step[level+1]+j][r]));
                 }
              plus_equal_TensProd(&(GC->ml_polyplaq[level][r]), &TP);
              }
            }

         } // end of update

      // normalize polycorr[level] and polyplaq[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
         if(slice==0)
           {
           times_equal_real_TensProd(&(GC->ml_polyplaq[level][r]), 1.0/(double) param->d_ml_upd[level]);
           }
         }

      break;
      // end of the not innermost not outermost level
    } // end of switch
  } // end of multilevel


// multilevel for flux tube with disconnected correlator to be used in long simulations
void multilevel_tube_disc_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int dt,
                               int iteration)
  {
  int upd;
  long raux;

  if(dt!=param->d_ml_step[0])
    {
    fprintf(stderr, "'dt' has to be equal to ml_step[0] in multilevel_tube_disc_long (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialyze ml_polycorr[0] and ml_polyplaq[0] to 0 if needed
  if(iteration==0)
    {
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
       if(slice==0)
         {
         zero_TensProd(&(GC->ml_polyplaq[0][r]));
         }
       }
    }

  // perform the update
  for(upd=0; upd < param->d_ml_upd[0]; upd++)
     {
     // update on level zero
     update_for_multilevel(GC, geo, param, 0);

     #if NLEVELS==1
       // compute Polyakov loop and plaquette restricted to the slice
       compute_local_poly_and_plaq(GC, geo, param);

       // compute the tensor products
       // and update ml_polycorr[0] and ml_polyplaq[0]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          TensProd TP;
          long r1, r2;
          int j, t_tmp;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          r1=sisp_and_t_to_si(geo, r, 0);
          for(j=0; j<param->d_dist_poly; j++)
             {
             r1=nnp(geo, r1, 1);
             }
          si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

          TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);

          if(slice==0)
            {
            times_equal_complex_TensProd(&TP, GC->loc_plaq[r]);
            plus_equal_TensProd(&(GC->ml_polyplaq[0][r]), &TP);
            }
          }
     #else  // NLEVELS>1
       // initialyze ml_polycorr[1] and ml_polyplaq[1] to 0
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[1]; raux++)
          {
          long r = raux/(geo->d_size[0]/param->d_ml_step[1]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[1]) );

          zero_TensProd(&(GC->ml_polycorr[1][slice][r]));
          if(slice==0)
            {
            zero_TensProd(&(GC->ml_polyplaq[1][r]));
            }
          }

       // call inner levels
       // important: we have to call the "non long" version in inner levels
       multilevel_tube_disc(GC,
                            geo,
                            param,
                            param->d_ml_step[1]);

       // update polycorr[0] with polycorr[1]
       // and analogously for polyplaq
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          int j;
          TensProd TP;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          one_TensProd(&TP);
          for(j=0; j<param->d_ml_step[0]/param->d_ml_step[1]; j++)
             {
             times_equal_TensProd(&TP, &(GC->ml_polycorr[1][slice*param->d_ml_step[0]/param->d_ml_step[1]+j][r]));
             }
          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);

          if(slice==0)
            {
            equal_TensProd(&TP, &(GC->ml_polyplaq[1][r]));
            for(j=1; j<param->d_ml_step[0]/param->d_ml_step[1]; j++)
               {
               times_equal_TensProd(&TP, &(GC->ml_polycorr[1][slice*param->d_ml_step[0]/param->d_ml_step[1]+j][r]));
               }
            plus_equal_TensProd(&(GC->ml_polyplaq[0][r]), &TP);
            }
          }
     #endif
     } // end update

  if(iteration==param->d_ml_level0_repeat-1) // iteration starts from zero
    {
    // normalize polycorr[level] and polyplaq[level]
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       times_equal_real_TensProd(&(GC->ml_polycorr[0][slice][r]), 1.0/((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat) );

       if(slice==0)
         {
         times_equal_real_TensProd(&(GC->ml_polyplaq[0][r]), 1.0/((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat) );
         }
       }
    }
  } // end of multilevel


// compute polyakov loop, plaquette and connected plaquette on a single slice
void compute_local_poly_plaq_and_plaqconn(Gauge_Conf *GC,
                                          Geometry const * const geo,
                                          GParam const * const param)
  {
  int num_hit;
  long raux;

  if(param->d_dist_poly>1 && geo->d_size[1]-param->d_dist_poly>1) // Polyakov loops are separated along the "1" direction
    {
    num_hit=param->d_multihit;
    }
  else
    {
    num_hit=0;
    }

  #ifdef THETA_MODE
  // clovers are eventually needed by the multihit
  compute_clovers(GC, geo, 0);
  #endif

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(raux)
  #endif
  for(raux=0; raux<(geo->d_space_vol*geo->d_size[0]/param->d_ml_step[NLEVELS-1]); raux++)
     {
     int i, j, t;
     long r4;
     GAUGE_GROUP matrix;

     const long r = raux/(geo->d_size[0]/param->d_ml_step[NLEVELS-1]);
     const int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[NLEVELS-1]) );

     one(&(GC->loc_poly[slice][r]));
     for(i=0; i<param->d_ml_step[NLEVELS-1]; i++)
        {
        t=slice*(param->d_ml_step[NLEVELS-1])+i;
        multihit(GC,
                 geo,
                 param,
                 sisp_and_t_to_si(geo, r, t),
                 0,
                 num_hit,
                 &matrix);
        times_equal(&(GC->loc_poly[slice][r]), &matrix);
        }

     if(slice==0)
       {
       r4=sisp_and_t_to_si(geo, r, 0); // t=0 starting point

       multihit(GC, geo, param, r4, 0, num_hit, &matrix);
       equal(&(GC->loc_plaqconn[r]), &matrix);

       r4=nnp(geo, r4, 0); // now we are in the t=1 plane

       #ifdef DEBUG
         long r4start=r4;
       #endif

       // polyakov loop are separated along the "1" direction
       for(j=0; j<param->d_dist_poly/2; j++)
          {
          times_equal(&(GC->loc_plaqconn[r]), &(GC->lattice[r4][1]));
          r4=nnp(geo, r4, 1);
          }

       // the transverse direction is the "2" one
       for(j=0; j<param->d_trasv_dist; j++)
          {
          times_equal(&(GC->loc_plaqconn[r]), &GC->lattice[r4][2]);
          r4=nnp(geo, r4, 2);
          }

       plaquettep_matrix(GC, geo, r4, param->d_plaq_dir[0], param->d_plaq_dir[1], &matrix);

       GC->loc_plaq[r]=retr(&matrix)+I*imtr(&matrix);
       times_equal(&(GC->loc_plaqconn[r]), &matrix);

       // the transverse direction is the "2" one
       for(j=0; j<param->d_trasv_dist; j++)
          {
          r4=nnm(geo, r4, 2);
          times_equal_dag(&(GC->loc_plaqconn[r]), &(GC->lattice[r4][2]));
          }

       // polyakov loop are separated along the "1" direction
       for(j=0; j<param->d_dist_poly/2; j++)
          {
          r4=nnm(geo, r4, 1);
          times_equal_dag(&(GC->loc_plaqconn[r]), &(GC->lattice[r4][1]));
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
          times_equal(&(GC->loc_plaqconn[r]), &matrix);
          r4=nnp(geo, r4, 0);
          }
       }// end of if(slice==0)
     }
  }


// multilevel for flux width computation using the connected correlator
void multilevel_tube_conn(Gauge_Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int dt)
  {
  int upd, level;
  long raux;

  // remember that d_size[0] >= d_ml_step[0]>d_ml_step[1]> d_ml_step[2] ...

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

      // initialyze ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] to 0
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

         zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
         if(slice==0)
           {
           zero_TensProd(&(GC->ml_polyplaq[0][r]));
           zero_TensProd(&(GC->ml_polyplaqconn[0][r]));
           }
         }

      // call lower levels
      multilevel_tube_conn(GC,
                           geo,
                           param,
                           param->d_ml_step[0]);

      break;
      // end of the outermost level

    case NLEVELS-1 : // INNERMOST LEVEL

      // in case level -1 is never used
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           if(slice==0)
             {
             zero_TensProd(&(GC->ml_polyplaq[0][r]));
             zero_TensProd(&(GC->ml_polyplaqconn[0][r]));
             }
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel(GC, geo, param, level);

         // compute Polyakov loop, plaquette and connected plaquette restricted to the slice
         compute_local_poly_plaq_and_plaqconn(GC, geo, param);

         // compute the tensor products
         // and update ml_polycorr[level], ml_polyplaq[level] and ml_polyplaqconn[level]
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            TensProd TP;
            long r1, r2;
            int j, t_tmp;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            r1=sisp_and_t_to_si(geo, r, 0);
            for(j=0; j<param->d_dist_poly; j++)
               {
               r1=nnp(geo, r1, 1);
               }
            si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

            TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);

            if(slice==0)
              {
              times_equal_complex_TensProd(&TP, GC->loc_plaq[r]);
              plus_equal_TensProd(&(GC->ml_polyplaq[level][r]), &TP);

              TensProd_init(&TP, &(GC->loc_plaqconn[r]), &(GC->loc_poly[slice][r2]) );
              plus_equal_TensProd(&(GC->ml_polyplaqconn[level][r]), &TP);
              }
            }
        } // end of update

      // normalize polycorr and polyplaq
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
         if(slice==0)
           {
           times_equal_real_TensProd(&(GC->ml_polyplaq[level][r]), 1.0/(double) param->d_ml_upd[level]);
           times_equal_real_TensProd(&(GC->ml_polyplaqconn[level][r]), 1.0/(double) param->d_ml_upd[level]);
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

      // in case level -1 is never used
      if(level==0 && geo->d_size[0]==param->d_ml_step[0])
        {
        // initialyze ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] to 0
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(raux)
        #endif
        for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
           {
           long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
           int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

           zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
           if(slice==0)
             {
             zero_TensProd(&(GC->ml_polyplaq[0][r]));
             zero_TensProd(&(GC->ml_polyplaqconn[0][r]));
             }
           }
        }

      // perform the update
      for(upd=0; upd< param->d_ml_upd[level]; upd++)
         {
         update_for_multilevel(GC, geo, param, level);

         // initialyze ml_polycorr[level+1], ml_polyplaq[level+1] and ml_polyplaqconn[level+1] to 0
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level+1]; raux++)
            {
            long r = raux/(geo->d_size[0]/param->d_ml_step[level+1]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level+1]) );

            zero_TensProd(&(GC->ml_polycorr[level+1][slice][r]));
            if(slice==0)
              {
              zero_TensProd(&(GC->ml_polyplaq[level+1][r]));
              zero_TensProd(&(GC->ml_polyplaqconn[level+1][r]));
              }
            }

         // call higher levels
         multilevel_tube_conn(GC,
                              geo,
                              param,
                              param->d_ml_step[level+1]);

         // update polycorr[level] with polycorr[level+1]
         // and analogously for polyplaq and polyplaqconn
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(raux)
         #endif
         for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
            {
            int j;
            TensProd TP;

            long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
            int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

            one_TensProd(&TP);
            for(j=0; j<param->d_ml_step[level]/param->d_ml_step[level+1]; j++)
               {
               times_equal_TensProd(&TP, &(GC->ml_polycorr[level+1][slice*param->d_ml_step[level]/param->d_ml_step[level+1]+j][r]));
               }
            plus_equal_TensProd(&(GC->ml_polycorr[level][slice][r]), &TP);

            if(slice==0)
              {
              equal_TensProd(&TP, &(GC->ml_polyplaq[level+1][r]));
              for(j=1; j<param->d_ml_step[level]/param->d_ml_step[level+1]; j++)
                 {
                 times_equal_TensProd(&TP, &(GC->ml_polycorr[level+1][slice*param->d_ml_step[level]/param->d_ml_step[level+1]+j][r]));
                 }
              plus_equal_TensProd(&(GC->ml_polyplaq[level][r]), &TP);

              equal_TensProd(&TP, &(GC->ml_polyplaqconn[level+1][r]));
              for(j=1; j<param->d_ml_step[level]/param->d_ml_step[level+1]; j++)
                 {
                 times_equal_TensProd(&TP, &(GC->ml_polycorr[level+1][slice*param->d_ml_step[level]/param->d_ml_step[level+1]+j][r]));
                 }
              plus_equal_TensProd(&(GC->ml_polyplaqconn[level][r]), &TP);
              }
            }
         } // end of update

      // normalize polycorr[level], polyplaq[level] and polyplaqconn[level]
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(raux)
      #endif
      for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[level]; raux++)
         {
         long r = raux/(geo->d_size[0]/param->d_ml_step[level]);
         int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[level]) );

         times_equal_real_TensProd(&(GC->ml_polycorr[level][slice][r]), 1.0/(double) param->d_ml_upd[level]);
         if(slice==0)
           {
           times_equal_real_TensProd(&(GC->ml_polyplaq[level][r]), 1.0/(double) param->d_ml_upd[level]);
           times_equal_real_TensProd(&(GC->ml_polyplaqconn[level][r]), 1.0/(double) param->d_ml_upd[level]);
           }
         }

      break;
      // end of the not innermost not outermost level
    } // end of switch
  } // end of multilevel



// multilevel for flux tube with connected correlator to be used in long simulations
void multilevel_tube_conn_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int dt,
                               int iteration)
  {
  int upd;
  long raux;

  if(dt!=param->d_ml_step[0])
    {
    fprintf(stderr, "'dt' has to be equal to ml_step[0] in multilevel_tube_conn_long (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialyze ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] to 0 if needed
  if(iteration==0)
    {
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       zero_TensProd(&(GC->ml_polycorr[0][slice][r]));
       if(slice==0)
         {
         zero_TensProd(&(GC->ml_polyplaq[0][r]));
         zero_TensProd(&(GC->ml_polyplaqconn[0][r]));
         }
       }
    }

  // perform the update
  for(upd=0; upd < param->d_ml_upd[0]; upd++)
     {
     // update on level zero
     update_for_multilevel(GC, geo, param, 0);

     #if NLEVELS==1
       // compute Polyakov loop, plaquette and connected plaquette restricted to the slice
       compute_local_poly_plaq_and_plaqconn(GC, geo, param);

       // compute the tensor products
       // and update ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0]
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          TensProd TP;
          long r1, r2;
          int j, t_tmp;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          r1=sisp_and_t_to_si(geo, r, 0);
          for(j=0; j<param->d_dist_poly; j++)
             {
             r1=nnp(geo, r1, 1);
             }
          si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

          TensProd_init(&TP, &(GC->loc_poly[slice][r]), &(GC->loc_poly[slice][r2]) );
          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);

          if(slice==0)
            {
            times_equal_complex_TensProd(&TP, GC->loc_plaq[r]);
            plus_equal_TensProd(&(GC->ml_polyplaq[0][r]), &TP);

            TensProd_init(&TP, &(GC->loc_plaqconn[r]), &(GC->loc_poly[slice][r2]) );
            plus_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &TP);
            }
          }
     #else  // NLEVELS>1
       // initialyze ml_polycorr[1], ml_polyplaq[1] and ml_polyplaqconn[1] to 0
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[1]; raux++)
          {
          long r = raux/(geo->d_size[0]/param->d_ml_step[1]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[1]) );

          zero_TensProd(&(GC->ml_polycorr[1][slice][r]));
          if(slice==0)
            {
            zero_TensProd(&(GC->ml_polyplaq[1][r]));
            zero_TensProd(&(GC->ml_polyplaqconn[1][r]));
            }
          }

       // call inner levels
       // important: we have to call the "non long" version in inner levels
       multilevel_tube_conn(GC,
                            geo,
                            param,
                            param->d_ml_step[1]);

       // update polycorr[0] with polycorr[1]
       // and analogously for polyplaq and polyplaqconn
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(raux)
       #endif
       for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
          {
          int j;
          TensProd TP;

          long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
          int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

          one_TensProd(&TP);
          for(j=0; j<param->d_ml_step[0]/param->d_ml_step[1]; j++)
             {
             times_equal_TensProd(&TP, &(GC->ml_polycorr[1][slice*param->d_ml_step[0]/param->d_ml_step[1]+j][r]));
             }
          plus_equal_TensProd(&(GC->ml_polycorr[0][slice][r]), &TP);

          if(slice==0)
            {
            equal_TensProd(&TP, &(GC->ml_polyplaq[1][r]));
            for(j=1; j<param->d_ml_step[0]/param->d_ml_step[1]; j++)
               {
               times_equal_TensProd(&TP, &(GC->ml_polycorr[1][slice*param->d_ml_step[0]/param->d_ml_step[1]+j][r]));
               }
            plus_equal_TensProd(&(GC->ml_polyplaq[0][r]), &TP);

            equal_TensProd(&TP, &(GC->ml_polyplaqconn[1][r]));
            for(j=1; j<param->d_ml_step[0]/param->d_ml_step[1]; j++)
               {
               times_equal_TensProd(&TP, &(GC->ml_polycorr[1][slice*param->d_ml_step[0]/param->d_ml_step[1]+j][r]));
               }
            plus_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &TP);
            }
          }
     #endif
     } // end update

  if(iteration==param->d_ml_level0_repeat-1) // iteration starts from zero
    {
    // normalize polycorr[level], polyplaq[level] and polyplaqconn[level]
    #ifdef OPENMP_MODE
    #pragma omp parallel for num_threads(NTHREADS) private(raux)
    #endif
    for(raux=0; raux<geo->d_space_vol*geo->d_size[0]/param->d_ml_step[0]; raux++)
       {
       long r = raux/(geo->d_size[0]/param->d_ml_step[0]);
       int slice = (int) (raux % (geo->d_size[0]/param->d_ml_step[0]) );

       times_equal_real_TensProd(&(GC->ml_polycorr[0][slice][r]), 1.0/((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat) );

       if(slice==0)
         {
         times_equal_real_TensProd(&(GC->ml_polyplaq[0][r]), 1.0/((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat) );
         times_equal_real_TensProd(&(GC->ml_polyplaqconn[0][r]), 1.0/((double) param->d_ml_upd[0] * (double) param->d_ml_level0_repeat) );
         }
       }
    }
  } // end of multilevel

#endif

