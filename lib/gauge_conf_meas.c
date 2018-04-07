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

//#define DEBUG

// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const restrict GC,
                  Geometry const * const restrict geo,
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
void plaquette(Gauge_Conf const * const restrict GC,
               Geometry const * const restrict geo,
               GParam const * const restrict param,
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
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, r, i, j);
         }
     
      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, r, i, j);
            }
         }
      }

   if(STDIM>2)
     {
     ps*=param->d_inv_vol;
     ps/=((double) STDIM*(STDIM-1)/2-(STDIM-1));
     }
   else
     {
     ps=0.0;
     }

   pt*=param->d_inv_vol;
   pt/=((double) STDIM-1);

   *plaqs=ps;
   *plaqt=pt;
   }


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const restrict GC,
              Geometry const * const restrict geo,
              GParam const * const restrict param,
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


/*
// compute the topological charge
double topcharge(Gauge_Conf const * const GC,
                 GParam const * const param)
   {
   GAUGE_GROUP aux1, aux2, aux3;
   double ris, real1, real2, loc_charge; 
   const double chnorm=1.0/(128.0*PI*PI);
   long r;
   int i, dir[4][4], sign;

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
         quadrifoglio(GC, r, dir[1][i], dir[2][i], &aux1);
         quadrifoglio(GC, r, dir[3][i], dir[0][i], &aux2);

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


// compute the topological charge density at point "r"
double topchargedens(Gauge_Conf const * const GC,
                     long r)
   {
   GAUGE_GROUP aux1, aux2, aux3;
   double real1, real2, loc_charge; 
   const double chnorm=1.0/(128.0*PI*PI);
   int dir[4][4], sign, i;

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

   sign=1;
   loc_charge=0.0;

   for(i=1; i<4; i++)
      {
      quadrifoglio(GC, r, dir[1][i], dir[2][i], &aux1);
      quadrifoglio(GC, r, dir[3][i], dir[0][i], &aux2);
 
         times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
         real1=retr(&aux3)*NCOLOR;

         times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
         real2=retr(&aux3)*NCOLOR;
        
      loc_charge+=((double) sign*(real1-real2));
      sign=-sign;
      }
   return (loc_charge*chnorm); 
   }



// compute GParam::d_nummeas values of the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf const * const GC,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq)
   {
   if(param->d_cooling>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_nummeas); iter++)
        {
        cooling(&helperconf, param, param->d_cooling);

        ris=topcharge(&helperconf, param);
        charge[iter]=ris;

        plaquette(&helperconf, param, &plaqs, &plaqt);
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }

     end_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, param);
     plaquette(GC, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_nummeas); iter++)
        {
        charge[iter]=ris;
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }
     } 
   }

*/


void perform_measures_localobs(Gauge_Conf const * const restrict GC,
                               Geometry const * const restrict geo,
                               GParam const * const restrict param,
                               FILE *datafilep)
   {
   double plaqs, plaqt, polyre, polyim;

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov(GC, geo, param, &polyre, &polyim);


   fprintf(datafilep, "%.12lf %.12lf %.12lf %.12lf\n", plaqs, plaqt, polyre, polyim);
   fflush(datafilep);
   }


// to optimize the number of hits to be used in multilevel
void optimize_multihit(Gauge_Conf *restrict GC,
                       Geometry const * const restrict geo,
                       GParam const * const restrict param,
                       FILE *datafilep)
  {
  const int max_hit=10;

  int i, mh;
  long r;
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

     poly_average=0.0;
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

     time(&time2);
     diff_sec = difftime(time2, time1);
     fprintf(datafilep, "%d  %.6g  (time:%g)\n", mh, poly_std*mh, diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel
void optimize_multilevel(Gauge_Conf *restrict GC,
                         Geometry const * const restrict geo,
                         GParam const * const restrict param,
                         FILE *datafilep)
   {
   int i;
   long r;
   double poly_std, poly_average;
   double *poly_array;

   poly_array = (double *) mymalloc(DOUBLE_ALIGN, (unsigned long) param->d_space_vol * sizeof(double));

   fprintf(datafilep, "Multilevel optimization for level1 of the multithit: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel(GC,
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


void perform_measures_polycorr_ml(Gauge_Conf *restrict GC,
                                  Geometry const * const restrict geo,
                                  GParam const * const restrict param,
                                  FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;

     multilevel(GC,
                geo,
                param,
                0,
                param->d_size[0]);

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr_ris[0][r]));
        }
     ris/=(double)param->d_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel(GC, geo, param, datafilep);
   #endif
   }




#endif
