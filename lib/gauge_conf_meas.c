#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/flavour_matrix.h"
#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"
#include"../include/tens_prod_adj.h"
#include"../include/su2_monopoles.h"
#include"../include/sun_monopoles.h"
#include"../include/u1_monopoles.h"


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


// compute the mean plaquettes (spatial, temporal) in fundamental and adjoint rep.
void plaquette_fundadj(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *plaqsf,
                       double *plaqtf,
                       double *plaqsa,
                       double *plaqta)
   {
   #if NCOLOR==1
     fprintf(stderr, "Computing the adjoint plaquette in U(1) is meaningless (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);

     // just to avoid warnings
     (void) GC;
     (void) geo;
     (void) param;
     (void) plaqsf;
     (void) plaqtf;
     (void) plaqsa;
     (void) plaqta;
   #else
     long r;
     double pfs, pft, pas, pat;

     pfs=0.0;
     pft=0.0;
     pas=0.0;
     pat=0.0;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pfs) reduction(+ : pft) reduction(+ : pas) reduction(+ : pat)
     #endif
     for(r=0; r<(param->d_volume); r++)
        {
        GAUGE_GROUP matrix;
        double rtr, itr;
        int i, j;

        i=0;
        for(j=1; j<STDIM; j++)
           {
           plaquettep_matrix(GC, geo, param, r, i, j, &matrix);

           rtr=retr(&matrix);
           itr=imtr(&matrix);

           pft+=rtr;
           pat+=(NCOLOR*NCOLOR*(rtr*rtr+itr*itr)-1)/(NCOLOR*NCOLOR-1);
           }

        for(i=1; i<STDIM; i++)
           {
           for(j=i+1; j<STDIM; j++)
              {
              plaquettep_matrix(GC, geo, param, r, i, j, &matrix);

              rtr=retr(&matrix);
              itr=imtr(&matrix);

              pfs+=rtr;
              pas+=(NCOLOR*NCOLOR*(rtr*rtr+itr*itr)-1)/(NCOLOR*NCOLOR-1);
              }
           }
        }

     if(STDIM>2)
       {
       pfs*=param->d_inv_vol;
       pfs/=((double) (STDIM-1)*(STDIM-2)/2);

       pas*=param->d_inv_vol;
       pas/=((double) (STDIM-1)*(STDIM-2)/2);
       }
     else
       {
       pfs=0.0;
       pas=0.0;
       }

     pft*=param->d_inv_vol;
     pft/=((double) STDIM-1);

     pat*=param->d_inv_vol;
     pat/=((double) STDIM-1);

     *plaqsf=pfs;
     *plaqtf=pft;

     *plaqsa=pas;
     *plaqta=pat;
   #endif
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


// compute the mean Polyakov loop in the adjoint representation (the trace of)
void polyakov_adj(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  double *repoly,
                  double *impoly)
   {
   long rsp;
   double rep, imp;
   double complex tr;

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
      tr=NCOLOR*retr(&matrix)+NCOLOR*imtr(&matrix)*I;

      #if NCOLOR==1
        (void) tr;
        rep+=0.0;
      #else
        rep+=(cabs(tr)*cabs(tr)-1)/(NCOLOR*NCOLOR-1);
      #endif

      imp+=0.0;
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }


// compute the mean Polyakov loop and its powers (trace of) in the presence of trace deformation
void polyakov_for_tracedef(Gauge_Conf const * const GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            double *repoly,
                            double *impoly)
   {
   long rsp;
   double **rep, **imp;
   int j, err;
   long i;

   for(j=0;j<(int)floor(NCOLOR/2);j++)
      {
      repoly[j]=0.0;
      impoly[j]=0.0;
      }

   err=posix_memalign((void**)&rep, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   err=posix_memalign((void**)&imp, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(i=0; i<param->d_space_vol; i++)
      {
      err=posix_memalign((void**)&(rep[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      err=posix_memalign((void**)&(imp[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      for(j=0; j<(int)floor(NCOLOR/2); j++)
         {
         rep[i][j] = 0.0;
         imp[i][j] = 0.0;
         }
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int k;
      GAUGE_GROUP matrix, matrix2;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(k=0; k<param->d_size[0]; k++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

       rep[rsp][0] = retr(&matrix);
       imp[rsp][0] = imtr(&matrix);

       equal(&matrix2, &matrix);

      for(k=1; k<(int)floor(NCOLOR/2.0); k++)
         {
         times_equal(&matrix2, &matrix);
         rep[rsp][k] = retr(&matrix2);
         imp[rsp][k] = imtr(&matrix2);
         }
      }

    for(j=0; j<(int)floor(NCOLOR/2); j++)
       {
       for(i=0; i<param->d_space_vol; i++)
          {
          repoly[j] += rep[i][j];
          impoly[j] += imp[i][j];
          }
       }

   for(j=0; j<(int)floor(NCOLOR/2.0); j++)
      {
      repoly[j] *= param->d_inv_space_vol;
      impoly[j] *= param->d_inv_space_vol;
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      free(rep[i]);
      free(imp[i]);
      }
   free(rep);
   free(imp);
   }


// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     long r)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     (void) GC;
     (void) geo;
     (void) param;
     (void) r;
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris=0.0; // initialized just to avoid compiler warnings

   #if (STDIM==4 && NCOLOR>1)
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
     ris=(loc_charge*chnorm);
   #endif

   #if (STDIM==2 && NCOLOR==1)
     GAUGE_GROUP u1matrix;
     double angle;

     plaquettep_matrix(GC, geo, param, r, 0, 1, &u1matrix);
     angle=atan2(cimag(u1matrix.comp), creal(u1matrix.comp))/PI2;

     ris=angle;
   #endif

   return ris;
   }


// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris;
   long r;

   ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      ris+=loc_topcharge(GC, geo, param, r);
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
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

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
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
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
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }


void perform_measures_localobs(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep,
                               FILE *monofilep)
   {
   double plaqs, plaqt, polyre, polyim;

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov(GC, geo, param, &polyre, &polyim);

   fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);

   // topological observables
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int i, err;
     double*charge, *meanplaq, charge_nocooling;

     charge_nocooling=topcharge(GC, geo, param);
     fprintf(datafilep, " %.12g ", charge_nocooling);

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

     free(charge);
     free(meanplaq);
   #else
     fprintf(datafilep, "\n");
   #endif
   fflush(datafilep);

   // monopole observables
   if(param->d_mon_meas == 1)
     {
     #if STDIM==4
     int subg, subgnum;
     Gauge_Conf helperconf;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     alloc_diag_proj_stuff(&helperconf, param);

     // MAG gauge fixing
     max_abelian_gauge_fix(&helperconf, geo, param);
 
     //diagonal projection
     diag_projection(&helperconf, param);
   
     //loop on all the U(1) subgroups
     if(NCOLOR>1)
       {
       subgnum=NCOLOR-1;
       }
     else
       {
       subgnum=1;
       }
     for(subg=0; subg<subgnum; subg++)
        {
        // extract the abelian component subg and save it to GC->u1_subg
        U1_extract(&helperconf, param, subg);

        // compute monopole observables
        monopoles_obs(&helperconf, geo, param, subg, monofilep);
        }

     free_diag_proj_stuff(&helperconf, param);
     free_gauge_conf(&helperconf, param);

     fflush(monofilep);
     #else
     (void) monofilep;
     #endif
     }
   }


// perform local observables in the case of trace deformation, it computes all the order parameters
void perform_measures_localobs_with_tracedef(Gauge_Conf const * const GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep,
                                             FILE *monofilep)
   {
   int i;
   double plaqs, plaqt, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1]; // +1 just to avoid warning if NCOLOR=1

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov_for_tracedef(GC, geo, param, polyre, polyim);

   fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

   for(i=0; i<(int)floor(NCOLOR/2); i++)
      {
      fprintf(datafilep, "%.12g %.12g ", polyre[i], polyim[i]);
      }

   // topological observables
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int err;
     double *charge, *meanplaq, charge_nocooling;

     charge_nocooling=topcharge(GC, geo, param);

     fprintf(datafilep, "%.12g ", charge_nocooling);

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

     free(charge);
     free(meanplaq);
   #else
     fprintf(datafilep, "\n");
   #endif

   fflush(datafilep);

   // monopole observables
   if(param->d_mon_meas == 1)
     {
     #if(STDIM==4)
     Gauge_Conf helperconf;
     int subg, subgnum;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     alloc_diag_proj_stuff(&helperconf, param);

     // MAG gauge fixing
     max_abelian_gauge_fix(&helperconf, geo, param);

     //diagonal projection
     diag_projection(&helperconf, param);

     //loop on all the U(1) subgroups
     if(NCOLOR>1)
       {
       subgnum=NCOLOR-1;
       }
     else
       {
       subgnum=1;
       }
     for(subg=0; subg<subgnum; subg++)
        {
        // extract the abelian component subg and save it to GC->u1_subg
        U1_extract(&helperconf, param, subg);

        // compute monopole observables
        monopoles_obs(&helperconf, geo, param, subg, monofilep);
        }

     free_diag_proj_stuff(&helperconf, param);
     free_gauge_conf(&helperconf, param);

     fflush(monofilep);
     #else
     (void) monofilep;
     #endif
     }
   }


void perform_measures_localobs_fundadj(Gauge_Conf const * const GC,
                                       Geometry const * const geo,
                                       GParam const * const param,
                                       FILE *datafilep)
   {
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int i, err;
     double plaqfs, plaqft, plaqas, plaqat, polyre, polyim, *charge, *meanplaq, charge_nocooling;

     plaquette_fundadj(GC, geo, param, &plaqfs, &plaqft, &plaqas, &plaqat);
     polyakov(GC, geo, param, &polyre, &polyim);
     charge_nocooling=topcharge(GC, geo, param);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g %.12g %.12g %.12g ", plaqfs, plaqft, plaqas, plaqat, polyre, polyim, charge_nocooling);

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
     double plaqfs, plaqft, plaqas, plaqat, polyre, polyim;

     plaquette_fundadj(GC, geo, param, &plaqfs, &plaqft, &plaqas, &plaqat);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g %.12g %.12g", plaqfs, plaqft, plaqas, plaqat, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);
   #endif
   }


// compute the average value of \sum_{flavours} Re(H_x U_{x,mu} H_{x+mu})
void higgs_interaction(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *he)
  {
  long r;
  double ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<(param->d_volume); r++)
     {
     int i;
     double aux=0.0;
     GAUGE_VECS v1;
     GAUGE_GROUP matrix;

     for(i=0; i<STDIM; i++)
        {
        equal(&matrix, &(GC->lattice[r][i]));

        matrix_times_vector_all_vecs(&v1, &matrix, &(GC->higgs[nnp(geo, r, i)]));
        aux+=re_scal_prod_vecs(&(GC->higgs[r]), &v1);
        }

     ris+=aux;
     }

  ris/=(double) STDIM;
  ris*=param->d_inv_vol;

  *he=ris;
  }


// compute flavour related observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this function
//
// tildeG0=ReTr[(\sum_x Q_x)(\sum_y Q_y)]/volume/NHIGGS
// tildeGminp=ReTr[(\sum_x Q_xe^{ipx})(\sum_y Q_ye^{-ipy)]/volume/NHIGGS
//
// tildeG0 is susceptibility/NHIGGS, tildeGminp is used to compute the 2nd momentum correlation function
//
// tildeD0=conj(\sum_x D_x) (\sum_y D_y) / volume
// tildeDminp=(\sum_x D_x e^{ipx}) conj(\sum_y D_y e^{ipy}) /volume
//
// tildeD0 is a U1 susceptibility, tildeDminp is used to compute the 2nd momentum correlation function
void compute_flavour_observables(Gauge_Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp,
                                 double *tildeD0,
                                 double *tildeDminp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)param->d_size[1];
  double complex D, Dp;
  FMatrix Q, Qp, Qmp, tmp1, tmp2;

  // Q =sum_x Q_x
  // Qp=sum_x e^{ipx}Q_x
  // Qmp=sum_x e^{-ipx}Q_x
  //
  // D, Dp and are the analogous of Q and Qp for D

  D=0.0+0.0*I;
  Dp=0.0+0.0*I;

  zero_FMatrix(&Q);
  zero_FMatrix(&Qp);
  zero_FMatrix(&Qmp);
  for(r=0; r<(param->d_volume); r++)
     {
     equal_FMatrix(&tmp1, &(GC->Qh[r]));
     equal_FMatrix(&tmp2, &tmp1);

     plus_equal_FMatrix(&Q, &tmp1);
     D+=(GC->Dh[r]);

     si_to_cart(coord, r, param);

     times_equal_complex_FMatrix(&tmp1, cexp(I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qp, &tmp1);
     Dp+=((GC->Dh[r]) * cexp(I*((double)coord[1])*p) );

     times_equal_complex_FMatrix(&tmp2, cexp(-I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qmp, &tmp2);
     }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0=retr_FMatrix(&tmp1)*param->d_inv_vol;
  *tildeD0=creal(conj(D)*D)*param->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp=retr_FMatrix(&tmp1)*param->d_inv_vol;
  *tildeDminp=creal(Dp*conj(Dp))*param->d_inv_vol;
  }


// compute correlators of flavour observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this function
//
// corrQQ is the correlato ReTr[Q_x Q_{x+d}]/N_higgs
// corr0string0 is the correlator \sum_f Re[hf^{dag} U_{x,1}U_{x+1,1}....Q_{x+d-1,1} hf], where hf is the f-th flavour
// corr0string1 is the correlator Re[h0^{dag} U_{x,1}U_{x+1,1}....U_{x+d-1,1} h1], where h1 is the second flavour
void compute_flavour_observables_corr(Gauge_Conf const * const GC,
                                      Geometry const * const geo,
                                      GParam const * const param,
                                      double *corrQQ,
                                      double *corr0string0,
                                      double *corr0string1)
  {
  int dist;
  long r;
  double accumulator1, accumulator2;

  for(dist=0; dist<param->d_size[1]; dist++)
     {
     accumulator1=0.0;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : accumulator1)
     #endif
     for(r=0; r<(param->d_volume); r++)
        {
        int i;
        long r1;
        FMatrix tmp1;

        equal_FMatrix(&tmp1, &(GC->Qh[r]));
        r1=r;
        for(i=0; i<dist; i++)
           {
           r1=nnp(geo, r1, 1);
           }
        times_equal_FMatrix(&tmp1, &(GC->Qh[r1]));
        accumulator1+=retr_FMatrix(&tmp1);
        }
     accumulator1*=param->d_inv_vol;
     corrQQ[dist]=accumulator1;

     accumulator1=0.0;
     accumulator2=0.0;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : accumulator1) reduction(+ : accumulator2)
     #endif
     for(r=0; r<(param->d_volume); r++)
        {
        int i;
        long r1;
        GAUGE_VECS phi1, phi2;
        GAUGE_GROUP U;

        equal_vecs(&phi1, &(GC->higgs[r]));
        r1=r;
        one(&U);
        for(i=0; i<dist; i++)
           {
           times_equal(&U, &(GC->lattice[r1][1]));
           r1=nnp(geo, r1, 1);
           }
        matrix_times_vector_all_vecs(&phi2, &U, &(GC->higgs[r1]));
        accumulator1+=re_scal_prod_vecs(&phi1, &phi2);
        #if NHIGGS >1
         accumulator2+=re_scal_prod_single_vecs(&phi1, &phi2, 0, 1);
        #else
         accumulator2+=0.0;
        #endif
        }
     accumulator1*=param->d_inv_vol;
     accumulator2*=param->d_inv_vol;

     corr0string0[dist]=accumulator1;
     corr0string1[dist]=accumulator2;
     }
  }


void perform_measures_higgs(Gauge_Conf *GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            FILE *datafilep)
   {
   double plaqs, plaqt, polyre, polyim, he, tildeG0, tildeGminp, tildeD0, tildeDminp;
   long r;

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov(GC, geo, param, &polyre, &polyim);
   higgs_interaction(GC, geo, param, &he);

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      init_FMatrix_vecs(&(GC->Qh[r]), &(GC->higgs[r]));
      GC->Dh[r] = HiggsU1Obs_vecs(&(GC->higgs[r]));
      }

   compute_flavour_observables(GC,
                               param,
                               &tildeG0,
                               &tildeGminp,
                               &tildeD0,
                               &tildeDminp);

   fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
   fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
   fprintf(datafilep, "%.12g ", he);
   fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);
   fprintf(datafilep, "%.12g %.12g ", tildeD0, tildeDminp);

   /*
   // for correlators

   int err, i;
   double *corrQQ, *corr0string0, *corr0string1;
   err=posix_memalign((void**) &(corrQQ), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[1] * sizeof(double));
   err+=posix_memalign((void**) &(corr0string0), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[1] * sizeof(double));
   err+=posix_memalign((void**) &(corr0string1), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[1] * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating the correlators! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   compute_flavour_observables_corr(GC,
                                    geo,
                                    param,
                                    corrQQ,
                                    corr0string0,
                                    corr0string1);
   for(i=0; i<param->d_size[1]; i++)
      {
      fprintf(datafilep, "%.12g ", corrQQ[i]);
      }
   for(i=0; i<param->d_size[1]; i++)
      {
      fprintf(datafilep, "%.12g ", corr0string0[i]);
      }
    for(i=0; i<param->d_size[1]; i++)
      {
      fprintf(datafilep, "%.12g ", corr0string1[i]);
      }

   free(corrQQ);
   free(corr0string0);
   free(corr0string1);
   */

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


// this is a function to be used just to test some fine points
// most notably TrP^2=TrQ^2+1/NHIGGS
void perform_measures_higgs_for_testing(Gauge_Conf *GC,
                                        Geometry const * const geo,
                                        GParam const * const param,
                                        FILE *datafilep)
   {
   double plaqs, plaqt, polyre, polyim, he, p2, tildeG0, tildeGminp, tildeD0, tildeDminp;
   long r;

   plaquette(GC, geo, param, &plaqs, &plaqt);
   polyakov(GC, geo, param, &polyre, &polyim);
   higgs_interaction(GC, geo, param, &he);

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      init_FMatrix_vecs(&(GC->Qh[r]), &(GC->higgs[r]));
      GC->Dh[r] = HiggsU1Obs_vecs(&(GC->higgs[r]));
      }

   fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
   fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
   fprintf(datafilep, "%.12g ", he);

   compute_flavour_observables(GC,
                               param,
                               &tildeG0,
                               &tildeGminp,
                               &tildeD0,
                               &tildeDminp);

   fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);
   fprintf(datafilep, "%.12g %.12g ", tildeD0, tildeDminp);

   p2=0.0;
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+: p2)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      FMatrix tmp1, tmp2;
      equal_FMatrix(&tmp1, &(GC->Qh[r]));
      equal_FMatrix(&tmp2, &tmp1);
      times_equal_FMatrix(&tmp1, &tmp2);
      p2+=retr_FMatrix(&tmp1)*NHIGGS+1./NHIGGS;
      }
   p2*=param->d_inv_vol;

   fprintf(datafilep, "%.12g ", p2);

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }



// fix maximal abelian gauge
// following the procedure described in
// C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]
void max_abelian_gauge_fix(Gauge_Conf *GC,
                           Geometry const * const geo,
                           GParam const * const param)
   {
   int i, dir;
   long r;
   double lambda[NCOLOR];
   const double overrelaxparam=1.85; // 1.0 means no overrelaxation
   const double target=1.0e-8;
   double nondiag, nondiagaux;

   // inizialize the matrix lambda = diag((N-1)/2, (N-1)/2-1, ..., -(N-1)/2)
   for(i=0; i<NCOLOR; i++)
      {
      lambda[i] = ( (double) NCOLOR -1.)/2. - (double) i;
      }

   nondiag=1;
   while(nondiag > target)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
        #endif
        for(r=0; r<param->d_volume/2; r++)
           {
           GAUGE_GROUP G_mag, help, X_links[2*STDIM];   // X_links contains the 2*STDIM links used in the computation of X(n)

           // initialize X_links[2*STDIM] with the 2*STDIM links surrounding the point r
           // links 0 to (STDIM-1) are forward, while links STDIM to (2*STDIM-1) are backwards.
           for(dir=0; dir<STDIM; dir++)
              {
              equal(&(X_links[dir]), &(GC->lattice[r][dir]));
              equal(&(X_links[dir+STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
              }

           comp_MAG_gauge_transformation(X_links, lambda, overrelaxparam, &G_mag);
 
           // apply the gauge transformation
           for(dir=0; dir<STDIM; dir++)
              {
              times(&help, &G_mag, &(GC->lattice[r][dir]));
              equal(&(GC->lattice[r][dir]), &help);

              times_equal_dag(&(GC->lattice[nnm(geo, r, dir)][dir]), &G_mag);
              }
           }

        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
        #endif
        for(r=param->d_volume/2; r<param->d_volume; r++)
           {
           GAUGE_GROUP G_mag, help, X_links[2*STDIM];   // X_links contains the 2*STDIM links used in the computation of X(n)

           // initialize X_links[2*STDIM] with the 2*STDIM links surrounding the point r
           // links 0 to (STDIM-1) are forward, while links STDIM to (2*STDIM-1) are backwards.
           for(dir=0; dir<STDIM; dir++)
              {
              equal(&(X_links[dir]), &(GC->lattice[r][dir]));
              equal(&(X_links[dir+STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
              }

           comp_MAG_gauge_transformation(X_links, lambda, overrelaxparam, &G_mag);

           // apply the gauge transformation
           for(dir=0; dir<STDIM; dir++)
              {
              times(&help, &G_mag, &(GC->lattice[r][dir]));
              equal(&(GC->lattice[r][dir]), &help);

              times_equal_dag(&(GC->lattice[nnm(geo, r, dir)][dir]), &G_mag);
              }
           }

        // nondiagaux is the sum of the squares of the out-diagonal terms
        nondiagaux=0;

        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r, dir)  reduction(+ : nondiagaux)
        #endif
        for(r=0; r<param->d_volume; r++)
           {
           GAUGE_GROUP X_links[2*STDIM];   // X_links contains the 2*STDIM links used in the computation of X(n)
           double counter;

           for(dir=0; dir<STDIM; dir++)
              {
              equal(&(X_links[dir]), &(GC->lattice[r][dir]));
              equal(&(X_links[dir+STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
              }
           comp_outdiagnorm_of_X(X_links, lambda, &counter);
           nondiagaux += counter;
           }
     
        nondiag = nondiagaux * param->d_inv_vol / (double)NCOLOR / (double) NCOLOR;

        // printf("%g  %g\n", nondiag, nondiag/target);
        // fflush(stdout);
        }

   // unitarize all the links
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }
   } 


// extract the diagonal part of the links after gauge fixing.
// the phases are saved in GC->diag_proj but these are NOT the monopole phases (see U1_extract)
void diag_projection(Gauge_Conf *GC,
                     GParam const * const param)
   {
   int dir; 
   long r;

   for(r=0;r<param->d_volume;r++)
      {
      for(dir=0;dir<STDIM;dir++)
         {
         diag_projection_single_site(GC, &(GC->lattice[r][dir]), r, dir);
         }
      }
   }


// extract the abelian components of the link
// following the procedure described in
// Bonati, D'Elia https://arxiv.org/abs/1308.0302
// and save them in GC->u1_subg
//
// also intialize GC->uflag to zero
void U1_extract(Gauge_Conf *GC, 
                GParam const * const param,
                int subg)
   { 
   int dir, i;
   long r;    

   for(r=0;r<param->d_volume;r++)
      {
      for(dir=0;dir<STDIM;dir++)
         {
         GC->u1_subg[r][dir] = 0.0;
         for(i=0;i<=subg;i++)
            {
            GC->u1_subg[r][dir] += GC->diag_proj[r][dir][i];
            }

         GC->uflag[r][dir] = 0;
         }
      }
   }


// Compute the forward derivative of the abelian part of the plaquette Fjk in direction i.
// the angle is chosen in between -pi and pi.
void Di_Fjk(Gauge_Conf *GC,
            Geometry const * const geo,
            long r,
            int idir,
            int jdir,
            int kdir,
            double *DiFjk)

   {
   double prpi, pr; // pr -> plaquette at site r, prpi plaquette at site r+idir

//
//       ^ k
//       |
//       +---<---+
//       |       |
//       V       ^         pr
//       |       |
//       +--->---+---> j
//       r
//

   pr  = GC->u1_subg[r][jdir] - GC->u1_subg[r][kdir];
   pr += GC->u1_subg[nnp(geo, r, jdir)][kdir] - GC->u1_subg[nnp(geo, r, kdir)][jdir];


//
//       ^ k
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)        prpi
//       |       |
//       +--->---+---> j
//      r+i    (4)
//

   r=nnp(geo, r, idir);

   prpi  = GC->u1_subg[r][jdir] - GC->u1_subg[r][kdir];
   prpi += GC->u1_subg[nnp(geo, r, jdir)][kdir] - GC->u1_subg[nnp(geo, r, kdir)][jdir];
 
   *DiFjk = 2.0*(atan(tan(prpi/2.0)) - atan(tan(pr/2.0)));
   }


// compute the DeGrand-DeTar currents
int DeGrand_current(Gauge_Conf *GC,
                    Geometry const * const geo,
                    long r,
                    int dir)
   {
   if(STDIM!=4)
     {
     fprintf(stderr, "Wrong number of dimensions! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double der1, der2, der3;
   int ris;

   if(dir == 0)
     {
     Di_Fjk(GC, geo, r, 1, 2, 3, &der1);
     Di_Fjk(GC, geo, r, 3, 1, 2, &der2);
     Di_Fjk(GC, geo, r, 2, 1, 3, &der3);
   
     ris = (int) round( ((der1 + der2 - der3)/PI2) );
     }
   else if(dir ==1)
          {
          Di_Fjk(GC, geo, r, 3, 2, 0, &der1);
          Di_Fjk(GC, geo, r, 0, 3, 2, &der2);
          Di_Fjk(GC, geo, r, 2, 3, 0, &der3);

          ris = (int) round( ((der1 + der2 - der3)/PI2) );
          }
   else if(dir == 2)
          {
          Di_Fjk(GC, geo, r, 3, 0, 1, &der1);
          Di_Fjk(GC, geo, r, 0, 1, 3, &der2);
          Di_Fjk(GC, geo, r, 1, 0, 3, &der3);

          ris = (int) round( ((der1 + der2 - der3)/PI2) );
          }
   else
     {
     Di_Fjk(GC, geo, r, 0,2,1, &der1);
     Di_Fjk(GC, geo, r, 2,1,0, &der2);
     Di_Fjk(GC, geo, r, 1,2,0, &der3);

     ris = (int) round( ((der1 + der2 - der3)/PI2) );
     }

   return ris;
   } 


// search for monopole wrappings passing from r_tback
// this function can be invoked in two different ways
//
// or r=nnp(geo, r_tback, 0) and DeGrand_current(GC, geo, r_tback, 0)!=0 (forward case)
// or r=nnm(geo, r_tback, 0) and DeGrand_current(GC, geo, r, 0)!=0  (backward case)
//
// GC->uflag[][] is initialized in monopole_obs
//
// nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are associated to uflag[r][dir]
//
// num_wrap = number of wrappings
void wrap_search(Gauge_Conf *GC,
                 Geometry const * const geo,
                 GParam const * const param,
                 long r,
                 long r_tback,
                 int *num_wrap)
   {
   #if STDIM!=4
     fprintf(stderr, "Wrong number of dimensions! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
   #endif

   int dir, n_mu;

   if(r == r_tback)
     {
     return;
     }
   else
     {
     // forward case
     for(dir=0; dir<STDIM; dir++)
        {
        n_mu=DeGrand_current(GC, geo, nnp(geo, r, dir), dir);

        // if not all the monopole currents have been followed
        if(n_mu > GC->uflag[r][dir])
          {
          GC->uflag[r][dir] += 1;

          if( (geo->d_timeslice[r] == param->d_size[0]-1) && (dir == 0) )
            {
            *num_wrap += 1;
            }

          wrap_search(GC, geo, param, nnp(geo, r, dir), r_tback, num_wrap);

          return;
          }
        }

     //backward case
     for(dir=0;dir<STDIM;dir++)
        {
        n_mu=DeGrand_current(GC, geo, r, dir);

        if(n_mu < GC->uflag[nnm(geo, r, dir)][dir])
          {
          GC->uflag[nnm(geo, r, dir)][dir] -= 1;

          if( (geo->d_timeslice[r] == 0) && (dir == 0) )
            {
            *num_wrap -= 1;
            }

          wrap_search(GC, geo, param, nnm(geo, r, dir), r_tback, num_wrap);

          return;
          }
        }
     }
   }


// GC->uflag[][] has to be initialized to zero before calling this function
// (when GC->uflag is allocated it is also initialized to zero)
void monopoles_obs(Gauge_Conf *GC, 
                   Geometry const * const geo,
                   GParam const * const param, 
                   int subg, 
                   FILE* monofilep)
   {
   double mean_wrap;
   long r, rsp, r_tback, r_tbackback;
   int n_mu, num_wrap, mono_charge;
   int cartcoord[4];

   mean_wrap = 0.0;     // mean value of monopole wraps for unit volume

   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      r = sisp_and_t_to_si(geo, rsp, 1);                             // t=1 slice
      r_tback = sisp_and_t_to_si(geo, rsp, 0);                       // t=0 slice
      r_tbackback = sisp_and_t_to_si(geo, rsp, param->d_size[0]-1);  // t=T-1 slice

      // check the t=1 temporal slice to find monopoles currents
      n_mu=DeGrand_current(GC, geo, r, 0);

      // start following monopole charge in forward direction. Maximum lattice charge is +2 so we try twice
      for(mono_charge = 0; mono_charge<2; mono_charge++)
         {
         // nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are associated to uflag[r][dir]
         if(n_mu > GC->uflag[r_tback][0])
           {
           GC->uflag[r_tback][0] += 1;

           num_wrap = 0;
           wrap_search(GC, geo, param, r, r_tback, &num_wrap);

           mean_wrap += abs(num_wrap);

           lexeo_to_cart(cartcoord, r_tback, param);
           if(n_mu == 1)
             {
             fprintf(monofilep, "%ld ", GC->update_index);

             for(int k = 0; k< 4; k++)
                {
                fprintf(monofilep, "%d ", cartcoord[k]);
                }
             fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
             }
           else if(GC->uflag[r_tback][0] == 1) // this is to print only once monopole of charge +2
                  {
                  fprintf(monofilep, "%ld ", GC->update_index);

                  for(int k = 0; k<4; k++)
                     {
                     fprintf(monofilep, "%d ", cartcoord[k]);
                     }
                  fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
                  }
           }
         }

      n_mu=DeGrand_current(GC, geo, r_tback, 0);

      // start following monopole charge in backward direction. Maximum lattice charge is +2 so we try twice
      for(mono_charge = 0; mono_charge<2; mono_charge++)
         {
         // nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are associated to uflag[r][dir]
         if(n_mu < GC->uflag[r_tbackback][0])
           {
           GC->uflag[r_tbackback][0] -= 1;

           num_wrap = -1;
           wrap_search(GC, geo, param, r_tbackback, r_tback, &num_wrap);

           lexeo_to_cart(cartcoord, r_tback, param);
           if(n_mu == -1)
             {
             fprintf(monofilep, "%ld ", GC->update_index);

             for(int k = 0; k<4; k++)
                {
                fprintf(monofilep, "%d ", cartcoord[k]);
                }
             fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
             }
           else if(GC->uflag[r][0] == -1)  // this is to print only once monopole of charge +2
                  {
                  fprintf(monofilep, "%ld ", GC->update_index);

                  for(int k=0; k<4; k++)
                     {
                     fprintf(monofilep, "%d ", cartcoord[k]);
                     }
                  fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
                  }
           }
         }
      }
   }


#endif












