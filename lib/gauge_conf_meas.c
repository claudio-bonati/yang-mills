#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"



// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

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


// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double complex plaquettep_complex(Gauge_Conf const * const GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  long r,
                                  int i,
                                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

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

   return retr(&matrix)+I*imtr(&matrix);
   }


// computation of the plaquette (matrix) in position r and positive directions i,j
void plaquettep_matrix(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix)
   {
   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ j
//       |   (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> i
//       r   (1)
//

   equal(matrix, &(GC->lattice[r][i]));
   times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(matrix, &(GC->lattice[r][j]));
   }


// compute the four-leaf clover in position r, in the plane i,j and save it in M
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            GParam const * const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M)
   {
   GAUGE_GROUP aux;
   long k, p;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= STDIM || j >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param;
   #endif

   zero(M);

//
//                   i ^
//                     |
//             (14)    |     (3)
//         +-----<-----++-----<-----+
//         |           ||           |
//         |           ||           |
//   (15)  V      (13) ^V (4)       ^ (2)
//         |           ||           |
//         |   (16)    || r   (1)   |
//    p    +----->-----++----->-----+------>   j
//         +-----<-----++-----<-----+
//         |    (9)    ||   (8)     |
//         |           ||           |
//    (10) V      (12) ^V (5)       ^ (7)
//         |           ||           |
//         |           ||           |
//         +------>----++----->-----+
//              (11)   k      (6)
//
   // avanti-avanti
   equal(&aux, &(GC->lattice[r][i]) );                           // 1
   times_equal(&aux, &(GC->lattice[nnp(geo, r, i)][j]) );        // 2
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, r, j)][i]) );    // 3
   times_equal_dag(&aux, &(GC->lattice[r][j]) );                 // 4
   plus_equal(M, &aux);

   k=nnm(geo, r, j);

   // avanti-indietro
   equal_dag(&aux, &(GC->lattice[k][j]) );                       // 5
   times_equal(&aux, &(GC->lattice[k][i]) );                     // 6
   times_equal(&aux, &(GC->lattice[nnp(geo, k, i)][j]) );        // 7
   times_equal_dag(&aux, &(GC->lattice[r][i]) );                 // 8
   plus_equal(M, &aux);

   p=nnm(geo, r, i);

   // indietro-indietro
   equal_dag(&aux, &(GC->lattice[p][i]) );                       // 9
   times_equal_dag(&aux, &(GC->lattice[nnm(geo, k, i)][j]) );    // 10
   times_equal(&aux, &(GC->lattice[nnm(geo, k, i)][i]) );        // 11
   times_equal(&aux, &(GC->lattice[k][j]) );                     // 12
   plus_equal(M, &aux);

   // indietro-avanti
   equal(&aux, &(GC->lattice[r][j]) );                            // 13
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, p, j)][i]) );     // 14
   times_equal_dag(&aux, &(GC->lattice[p][j]) );                  // 15
   times_equal(&aux, &(GC->lattice[p][i]) );                      // 16
   plus_equal(M, &aux);
   }


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt)
   {
   long r;
   double ps, pt;

   ps=0.0;
   pt=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i, j;
      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, param, r, i, j);
         }
     
      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, param, r, i, j);
            }
         }
      }

   if(STDIM>2)
     {
     ps*=param->d_inv_vol;
     ps/=((double) (STDIM-1)*(STDIM-2)/2);
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


// compute the clover discretization of
// sum_{\mu\nu}  Tr(F_{\mu\nu}F_{\mu\nu})/2
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        double *energy)
  {
  long r;
  double ris;

  ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<param->d_volume; r++)
     {
     int i, j;
     GAUGE_GROUP aux1, aux2;

     for(i=0; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           clover(GC, geo, param, r, i, j, &aux1);

           ta(&aux1);
           equal(&aux2, &aux1);
           times_equal(&aux1, &aux2);
           ris+=-NCOLOR*retr(&aux1)/16.0;
           }
        }
     }

  *energy=ris*param->d_inv_vol;
  }


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly)
   {
   long rsp;
   double rep, imp;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

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
// see readme for more details
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   double ris;
   long r;

   if(STDIM !=4)
     {
     fprintf(stderr, "Wrong number of dimension! (%d instead of 4) (%s, %d)\n", STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      GAUGE_GROUP aux1, aux2, aux3;
      double real1, real2, loc_charge;
      const double chnorm=1.0/(128.0*PI*PI);
      int i, dir[4][3], sign;

      dir[0][0] = 0;
      dir[0][1] = 0;
      dir[0][2] = 0;

      dir[1][0] = 1;
      dir[1][1] = 2;
      dir[1][2] = 3;

      dir[2][0] = 2;
      dir[2][1] = 1;
      dir[2][2] = 1;

      dir[3][0] = 3;
      dir[3][1] = 3;
      dir[3][2] = 2;

      sign=-1;
      loc_charge=0.0;

      for(i=0; i<3; i++)
         {
         clover(GC, geo, param, r, dir[0][i], dir[1][i], &aux1);
         clover(GC, geo, param, r, dir[2][i], dir[3][i], &aux2);

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
        meanplaq[iter]=0.5*(plaqs+plaqt); // this is in 4d: 3 spatial and 3 temporal plaquettes
        }

     free_gauge_conf(&helperconf, param); 
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
   #if STDIM==4
   int i, err;
   double plaqs, plaqt, polyre, polyim, *charge, *meanplaq, charge_nocooling;

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov(GC, geo, param, &polyre, &polyim);
   charge_nocooling=topcharge(GC, geo, param);

   fprintf(datafilep, "%.12g %.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim, charge_nocooling);

   err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   topcharge_cooling(GC, geo, param, charge, meanplaq);
   for(i=0; i<param->d_coolrepeat; i++)
      {
      fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
      }
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(charge);
   free(meanplaq);

   #else

   double plaqs, plaqt, polyre, polyim;

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov(GC, geo, param, &polyre, &polyim);

   fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   #endif
   }


// to optimize the number of hits to be used in multilevel
void optimize_multihit_polycorr(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp, err;
  long r, r1, r2;
  double poly_std, poly_average, diff_sec;
  double complex *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, param, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
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
                    sisp_and_t_to_si(geo, r, i),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }
        poly_array[r]=retr(&matrix)+I*imtr(&matrix);
        }

     // average correlator computation
     poly_average=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
        poly_average+=cabs(poly_array[r]*conj(poly_array[r2]));
        }
     poly_average*=param->d_inv_space_vol;

     // std computation
     poly_std=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
        poly_std += pow(cabs(poly_array[r]*conj(poly_array[r2]))-poly_average, 2.0);
        }
     poly_std*=param->d_inv_space_vol;
     poly_std*=param->d_inv_space_vol;

     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g  (time:%g)\n", mh, poly_average*poly_average*mh, poly_std*mh, diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel
void optimize_multilevel_polycorr(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_std, poly_average;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorr(GC,
                       geo,
                       param,
                       param->d_size[0]);

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<param->d_space_vol; r++)
      {
      int j;
      for(j=1; j<param->d_size[0]/param->d_ml_step[0]; j++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][j][r]) );
         }
      }

   // averages
   poly_average=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      poly_average+=cabs(poly_array[r]);
      }
   poly_average*=param->d_inv_space_vol;

   // fluctuations
   poly_std=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_std+=pow(cabs(poly_array[r])-poly_average, 2.0);
      }
   poly_std*=param->d_inv_space_vol;
   poly_std*=param->d_inv_space_vol;

   // normalizations
   poly_average*=poly_average;
   for(i=0; i<NLEVELS; i++)
      {
      poly_average*=(double) param->d_ml_upd[i];
      }
   for(i=0; i<NLEVELS; i++)
      {
      poly_std*=(double) param->d_ml_upd[i];
      }

   fprintf(datafilep, "%.12g  %.12g ", poly_average, poly_std);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_polycorr(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;

     multilevel_polycorr(GC,
                geo,
                param,
                param->d_size[0]);

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_space_vol; r++)
        {
        int i;
        for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
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
     optimize_multilevel_polycorr(GC, geo, param, datafilep);
   #endif
   }


// to optimize the multilevel
void optimize_multilevel_polycorr_long(Gauge_Conf *GC,
                                       GParam const * const param,
                                       FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_std, poly_average;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<param->d_space_vol; r++)
      {
      int i;
      for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // average
   poly_average=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      poly_average+=cabs(poly_array[r]);
      }
   poly_average*=param->d_inv_space_vol;

   // fluctuation
   poly_std=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_std+=pow(cabs(poly_array[r])-poly_average, 2.0);
      }
   poly_std*=param->d_inv_space_vol;
   poly_std*=param->d_inv_space_vol;

   // normalization
   poly_average*=poly_average;
   for(i=0; i<NLEVELS; i++)
      {
      poly_average*=(double) param->d_ml_upd[i];
      }
   poly_average*=(double) param->d_ml_level0_repeat;

   for(i=0; i<NLEVELS; i++)
      {
      poly_std*=(double) param->d_ml_upd[i];
      }
   poly_std*=(double) param->d_ml_level0_repeat;

   fprintf(datafilep, "%.12g  %.12g ", poly_average, poly_std);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(%d) ", param->d_ml_level0_repeat);
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(poly_array);
   }


// print the value of the polyakov loop correlator that has been computed by multilevel
void perform_measures_polycorr_long(Gauge_Conf *GC,
                                    GParam const * const param,
                                    FILE *datafilep)
   {
   #ifdef OPT_MULTILEVEL
      optimize_multilevel_polycorr_long(GC, param, datafilep);
   #else
     double ris;
     long r;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_space_vol; r++)
        {
        int i;
        for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   }


// perform the computation of the string width with the
// disconnected correlator using the multilevel algorithm
void perform_measures_tube_disc(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;

   multilevel_tube_disc(GC,
                        geo,
                        param,
                        param->d_size[0]);

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<param->d_space_vol; r++)
      {
      int i;
      for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// perform the computation of the string width with the
// connected correlator using the multilevel algorithm
void perform_measures_tube_conn(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;

   multilevel_tube_conn(GC,
                        geo,
                        param,
                        param->d_size[0]);

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<param->d_space_vol; r++)
      {
      int i;
      for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// print the value of the the string width with the
// connected correlator that has been computed by multilevel
void perform_measures_tube_conn_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   double risr, risi;
   long r;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<param->d_space_vol; r++)
      {
      int i;
      for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }




#endif
