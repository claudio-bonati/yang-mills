#ifndef GAUGE_CONF_MEAS_MULTILEVEL_C
#define GAUGE_CONF_MEAS_MULTILEVEL_C

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
  double complex poly_corr;
  double poly_corr_abs, poly_corr_fluct, diff_sec;
  double complex *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) geo->d_space_vol * sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
  fprintf(datafilep, "the smaller the value the better the multihit\n");

  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop computation
     for(r=0; r<geo->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<geo->d_size[0]; i++)
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
     poly_corr=0.0+I*0.0;
     poly_corr_abs=0.0;
     for(r=0; r<geo->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr += poly_array[r]*conj(poly_array[r2]);
        poly_corr_abs += cabs(poly_array[r]*conj(poly_array[r2]));
        }
     poly_corr*=geo->d_inv_space_vol;
     poly_corr_abs*=geo->d_inv_space_vol;

     // fluctuation of the average correlator computation
     poly_corr_fluct=0.0;
     for(r=0; r<geo->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
        poly_corr_fluct+=cabs( poly_array[r]*conj(poly_array[r2]) - poly_corr );
        }
     poly_corr_fluct*=geo->d_inv_space_vol;


     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

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
   double complex poly_corr;
   double poly_corr_abs, poly_corr_fluct;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) geo->d_space_vol * sizeof(double complex));
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
                       geo->d_size[0]);

   for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<geo->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=geo->d_inv_space_vol;
   poly_corr_abs*=geo->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      poly_corr_fluct += cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=geo->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*= sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*= sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g ", poly_corr_abs);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// to optimize the multilevel
void optimize_multilevel_polycorr_with_higgs(Gauge_Conf *GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep)
   {
   int i, err;
   long r;
   double complex poly_corr;
   double poly_corr_abs, poly_corr_fluct;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) geo->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorr_with_higgs(GC,
                                  geo,
                                  param,
                                  geo->d_size[0]);

   for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<geo->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=geo->d_inv_space_vol;
   poly_corr_abs*=geo->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      poly_corr_fluct += cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=geo->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*= sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*= sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g ", poly_corr_abs);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
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
     int i;

     multilevel_polycorr(GC,
                         geo,
                         param,
                         geo->d_size[0]);

     for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<geo->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<geo->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=geo->d_inv_space_vol;

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

// perform the computation of the polyakov loop correlator with the multilevel algorithm
// for the theory with higgs fields
void perform_measures_polycorr_with_higgs(Gauge_Conf *GC,
                                          Geometry const * const geo,
                                          GParam const * const param,
                                          FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;
     int i;

     multilevel_polycorr_with_higgs(GC,
                                    geo,
                                    param,
                                    geo->d_size[0]);

     for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<geo->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<geo->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=geo->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorr(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_polycorr_with_higgs(GC, geo, param, datafilep);
   #endif
   }


// to optimize the multilevel
void optimize_multilevel_polycorr_long(Gauge_Conf *GC,
                                       Geometry const * const geo,
                                       GParam const * const param,
                                       FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_corr_abs, poly_corr_fluct;
   double complex poly_corr;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) geo->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<geo->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // average
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=geo->d_inv_space_vol;
   poly_corr_abs*=geo->d_inv_space_vol;

   // fluctuation
   poly_corr_fluct=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      poly_corr_fluct+=cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=geo->d_inv_space_vol;

   // normalization
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*=sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_ml_level0_repeat);
   poly_corr_fluct*=sqrt(param->d_ml_level0_repeat);

   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g %.12g ", poly_corr_abs, poly_corr_fluct);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "(%d) ", param->d_ml_level0_repeat);
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(poly_array);
   }


// print the value of the polyakov loop correlator that has been computed by multilevel
void perform_measures_polycorr_long(Gauge_Conf *GC,
                                    Geometry const * const geo,
                                    GParam const * const param,
                                    FILE *datafilep)
   {
   #ifdef OPT_MULTILEVEL
      optimize_multilevel_polycorr_long(GC, param, datafilep);
   #else
     double ris;
     long r;
     int i;

     for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<geo->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<geo->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=geo->d_inv_space_vol;

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
   int i;

   multilevel_tube_disc(GC,
                        geo,
                        param,
                        geo->d_size[0]);

   for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<geo->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

// perform the computation of the string width with the
// disconnected correlator that has been computed by multilevel (long version)
void perform_measures_tube_disc_long(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<geo->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
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
   int i;

   multilevel_tube_conn(GC,
                        geo,
                        param,
                        geo->d_size[0]);

   for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<geo->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// print the value of the the string width with the
// connected correlator that has been computed by multilevel
void perform_measures_tube_conn_long(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   for(i=1; i<geo->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<geo->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<geo->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=geo->d_inv_space_vol;
   risi*=geo->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


#endif
