#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<malloc.h>
#include<math.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"

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
   int aux_t, aux_s; // used to count the number of temporal and spatial plaquettes
   double ps, pt;

   ps=0.0;
   pt=0.0;

   aux_t=1;
   aux_s=1;

   for(r=0; r<(param->d_volume); r++)
      {
      aux_t=0;
      aux_s=0;

      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, r, i, j);
         aux_t++;
         }
     
      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, r, i, j);
            aux_s++;
            }
         }
      }

   #if STDIM > 2
     ps*=param->d_inv_vol;
     ps/=((double) aux_s);
   #else
     ps=0.0;
   #endif

   pt*=param->d_inv_vol;
   pt/=((double) aux_t);

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


#endif
