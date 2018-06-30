#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/mymalloc.h"
#include"../include/tens_prod.h"


// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  long r,
                  int i,
                  int j)
   {
   GAUGE_GROUP matrix;

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix);
   }


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt)
   {
   long r;
   int i, j;
   double ps, pt;

   ps=0.0;
   pt=0.0;

   for(r=0; r<(param->d_volume); r++)
      {
      i=0;
      for(j=1; j<param->d_stdim; j++)
         {
         pt+=plaquettep(GC, geo, r, i, j);
         }
     
      for(i=1; i<param->d_stdim; i++)
         {
         for(j=i+1; j<param->d_stdim; j++)
            {
            ps+=plaquettep(GC, geo, r, i, j);
            }
         }
      }

   if(param->d_stdim>2)
     {
     ps*=param->d_inv_vol;
     ps/=((double) param->d_stdim*(param->d_stdim-1)/2-(param->d_stdim-1));
     }
   else
     {
     ps=0.0;
     }

   pt*=param->d_inv_vol;
   pt/=((double) param->d_stdim-1);

   *plaqs=ps;
   *plaqt=pt;
   }


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly)
   {
   long r, rsp;
   int i;
   double rep, imp;
   GAUGE_GROUP matrix;

   rep=0.0;
   imp=0.0;

   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      r=sisp_and_t_to_si(rsp, 0, param);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

      rep+=retr(&matrix);
      imp+=imtr(&matrix);
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }


// compute the topological charge
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   GAUGE_GROUP aux1, aux2, aux3;
   double ris, real1, real2, loc_charge; 
   const double chnorm=1.0/(128.0*PI*PI);
   long r;
   int i, dir[4][4], sign;

   if(param->d_stdim!=4)
     {
     fprintf(stderr, "Wrong number of dimension! (%d instead of 4) (%s, %d)\n", param->d_stdim, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   dir[1][1] = 1;
   dir[1][2] = 1;
   dir[1][3] = 1;

   dir[2][1] = 2;
   dir[2][2] = 3;
   dir[2][3] = 0;

   dir[3][1] = 3;
   dir[3][2] = 2;
   dir[3][3] = 2;

   dir[0][1] = 0;
   dir[0][2] = 0;
   dir[0][3] = 3;

   ris=0.0;
   for(r=0; r<(param->d_volume); r++)
      {

      sign=1;
      loc_charge=0.0;

      for(i=1; i<4; i++)
         {
         clover(GC, geo, param, r, dir[1][i], dir[2][i], &aux1);
         clover(GC, geo, param, r, dir[3][i], dir[0][i], &aux2);

         times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
         real1=retr(&aux3)*NCOLOR;

         times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
         real2=retr(&aux3)*NCOLOR;
        
         loc_charge+=((double) sign*(real1-real2));
         sign=-sign;
         }
      ris+=(loc_charge*chnorm); 
      }

   return ris;
   }


// compute GParam::d_nummeas values of the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq)
   {
   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

        plaquette(&helperconf, geo, param, &plaqs, &plaqt);
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }

     end_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
     plaquette(GC, geo, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }
     } 
   }


void perform_measures_localobs(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep)
   {
   int i;
   double plaqs, plaqt, polyre, polyim, *charge, *meanplaq, charge_nocooling;

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov(GC, geo, param, &polyre, &polyim);
   charge_nocooling=topcharge(GC, geo, param);

   fprintf(datafilep, "%.12lf %.12lf %.12lf %.12lf %.12lf ", plaqs, plaqt, polyre, polyim, charge_nocooling);

   charge   = (double *) memalign(DOUBLE_ALIGN, (unsigned long) param->d_coolrepeat * sizeof(double));
   meanplaq = (double *) memalign(DOUBLE_ALIGN, (unsigned long) param->d_coolrepeat * sizeof(double));

   topcharge_cooling(GC, geo, param, charge, meanplaq);
   for(i=0; i<param->d_coolrepeat; i++)
      {
      fprintf(datafilep, "%.12f %.12f ", charge[i], meanplaq[i]);
      }
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(charge);
   free(meanplaq);
   }


// to optimize the number of hits to be used in multilevel
void optimize_multihit_polycorr(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp;
  long r, r1, r2;
  double poly_std, poly_average, diff_sec;
  double *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  poly_array = (double *) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(double));

  fprintf(datafilep, "Multihit optimization: ");
  fprintf(datafilep, "the smaller the 'ris' value, the best\n");
  fprintf(datafilep, "mhit  ris\n");
  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop computation
     for(r=0; r<param->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<param->d_size[0]; i++)
           {
           multihit(GC,
                    geo,
                    param,
                    sisp_and_t_to_si(r, i, param),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }
        poly_array[r]=retr(&matrix);
        }

     // average correlator computation
     poly_average=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(r, 0, param);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1
        poly_average+=poly_array[r]*poly_array[r2];
        }
     poly_average*=param->d_inv_space_vol;

     // std computation
     poly_std=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(r, 0, param);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, r1, param); // r2 is the spatial value of r1
        poly_std += (poly_array[r]*poly_array[r2]-poly_average)*(poly_array[r]*poly_array[r2]-poly_average);
        }
     poly_std*=param->d_inv_space_vol;
     poly_std*=param->d_inv_space_vol;

     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.6g  (time:%g)\n", mh, poly_std*mh, diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel
void optimize_multilevel_potQbarQ(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   int i;
   long r;
   double poly_std, poly_average;
   double *poly_array;

   poly_array = (double *) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(double));

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_pot_QbarQ(GC,
                        geo,
                        param,
                        0,
                        param->d_size[0]);

   // polyakov loop correlator
   poly_average=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr_ris[0][r]));
      poly_average+=poly_array[r];
      }
   poly_average*=param->d_inv_space_vol;

   poly_std=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_std+=(poly_average-poly_array[r])*(poly_average-poly_array[r]);
      }
   poly_std*=param->d_inv_space_vol;
   poly_std*=param->d_inv_space_vol;

   for(i=0; i<NLEVELS; i++)
      {
      poly_std*=(double) param->d_ml_upd[i];
      }

   fprintf(datafilep, "%.6g ", poly_std);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_pot_QbarQ(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;

     multilevel_pot_QbarQ(GC,
                geo,
                param,
                0,
                param->d_size[0]);

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr_ris[0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorr(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_potQbarQ(GC, geo, param, datafilep);
   #endif
   }


// print the value of the polyakov loop correlator that has been computed by multilevel
void perform_measures_pot_QbarQ_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   double ris;
   long r;

   ris=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      ris+=retr_TensProd(&(GC->ml_polycorr_ris[0][r]));
      }
   ris*=param->d_inv_space_vol;

   fprintf(datafilep, "%.12g\n", ris);
   fflush(datafilep);
   }


// to optimize the multilevel for stringwidth
void optimize_multilevel_stringQbarQ(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   int i;
   long r;
   double poly_std, poly_average;
   double *poly_array;

   poly_array = (double *) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(double));

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_pot_QbarQ(GC,
                        geo,
                        param,
                        0,
                        param->d_size[0]);

   // polyakov loop correlator
   poly_average=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr_ris[0][r]));
      poly_array[r]-=retr_TensProd(&(GC->ml_polyplaq_ris[0][r][0]));

      poly_average+=poly_array[r];
      }
   poly_average*=param->d_inv_space_vol;

   poly_std=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_std+=(poly_average-poly_array[r])*(poly_average-poly_array[r]);
      }
   poly_std*=param->d_inv_space_vol;
   poly_std*=param->d_inv_space_vol;

   for(i=0; i<NLEVELS; i++)
      {
      poly_std*=(double) param->d_ml_upd[i];
      }

   fprintf(datafilep, "%.6g ", poly_std);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_string_QbarQ(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep)
   {
   #ifdef OPT_MULTILEVEL
     optimize_multilevel_stringQbarQ(GC, geo, param, datafilep);
   #else
     int i;
     const int numplaqs=(param->d_stdim*(param->d_stdim-1))/2;
     double ris;
     long r;

     multilevel_string_QbarQ(GC,
                             geo,
                             param,
                             0,
                             param->d_size[0]);

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr_ris[0][r]));
        }
     ris*=param->d_inv_space_vol;
     fprintf(datafilep, "%.12g ", ris);

     for(i=0; i<numplaqs; i++)
        {
        ris=0.0;
        for(r=0; r<param->d_space_vol; r++)
           {
           ris+=retr_TensProd(&(GC->ml_polyplaq_ris[0][r][i]));
           }
        ris*=param->d_inv_space_vol;
        fprintf(datafilep, "%.12g ", ris);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);
   #endif
   }



#endif
