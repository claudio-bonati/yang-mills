#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<math.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

//#define DEBUG

// compute the staple in position r, direction i and save it in M
void calcstaples_w(Gauge_Conf const * const restrict GC,
                   Geometry const * const restrict geo,
                   long r,
                   int i,
                   GAUGE_GROUP * restrict M)
   {
   int j, l;
   long k;
   GAUGE_GROUP link1, link2, link3, link12, stap;

   zero(M); // M=0

   for(l=i+1; l<i+STDIM; l++)
      {
      j = (l % STDIM);

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//         |           |
//         |           V (2)
//         |           |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

      equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
      equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
      equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

      times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
      times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

      plus_equal(M, &stap);

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |           |
//     (2) V           |
//         |           |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

      k=nnm(geo, r, j);

      equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
      equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
      equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

      times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
      times(&stap, &link12, &link3);        // stap=link12*link3

      plus_equal(M, &stap);
      }
    }


// perform an update with heatbath
void heatbath_w(Gauge_Conf * restrict GC,
                Geometry const * const restrict geo,
                GParam const * const restrict param,
                long r,
                int i)
   {
   #ifdef DEBUG
   if(r > param->d_volume)
     {
     fprintf(stderr, "r too large: %ld > %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i > STDIM)
     {
     fprintf(stderr, "i too large: i=%d > %s (%s, %d)\n", i, QUOTEME(STDIM), __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   GAUGE_GROUP stap;

   calcstaples_w(GC, geo, r, i, &stap);
   single_heatbath(&(GC->lattice[r][i]), &stap, param);
   }


// perform an update with overrelaxation
void overrelaxation_w(Gauge_Conf * restrict GC,
                      Geometry const * const restrict geo,
                      GParam const * const restrict param,
                      long r,
                      int i)
   {
   #ifdef DEBUG
   if(r > param->d_volume)
     {
     fprintf(stderr, "r too large: %ld > %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i > STDIM)
     {
     fprintf(stderr, "i too large: i=%d > %s (%s, %d)\n", i, QUOTEME(STDIM), __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif
   (void) param; // just to avoid wornings

   GAUGE_GROUP stap;

   calcstaples_w(GC, geo, r, i, &stap);
   single_overrelaxation(&(GC->lattice[r][i]), &stap);
   }


/*
// compute the four-leaf clover in position r, in the plane i,j and save it in M
void quadrifoglio(Gauge_Conf const * const GC,
                  long r,
                  int j,
                  int i,
                  GAUGE_GROUP *M)
   {
   GAUGE_GROUP aux; 
   long k, p;

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
   equal(&aux, &(GC->lattice[r][j]) );                               // 1
   times_equal(&aux, &(GC->lattice[nnp(GC->geo, r, j)][i]) );        // 2
   times_equal_dag(&aux, &(GC->lattice[nnp(GC->geo, r, i)][j]) );    // 3
   times_equal_dag(&aux, &(GC->lattice[r][i]) );                     // 4
   plus_equal(M, &aux);
 
   k=nnm(GC->geo, r, i);

   // avanti-indietro
   equal_dag(&aux, &(GC->lattice[k][i]) );                       // 5
   times_equal(&aux, &(GC->lattice[k][j]) );                     // 6
   times_equal(&aux, &(GC->lattice[nnp(GC->geo, k, j)][i]) );    // 7
   times_equal_dag(&aux, &(GC->lattice[r][j]) );                 // 8
   plus_equal(M, &aux);

   p=nnm(GC->geo, r, j);

   // indietro-indietro
   equal_dag(&aux, &(GC->lattice[p][j]) );                           // 9
   times_equal_dag(&aux, &(GC->lattice[nnm(GC->geo, k, j)][i]) );    // 10
   times_equal(&aux, &(GC->lattice[nnm(GC->geo, k, j)][j]) );        // 11
   times_equal(&aux, &(GC->lattice[k][i]) );                         // 12
   plus_equal(M, &aux);

   // indietro-avanti
   equal(&aux, &(GC->lattice[r][i]) );                                // 13
   times_equal_dag(&aux, &(GC->lattice[nnp(GC->geo, p, i)][j]) );     // 14
   times_equal_dag(&aux, &(GC->lattice[p][i]) );                      // 15
   times_equal(&aux, &(GC->lattice[p][j]) );                          // 16
   plus_equal(M, &aux);
   }


// compute the antisymmetric part of the four-leaf clover in position r, in the plane i,j and
// save it in  GC->quadri (------IMPORTANT------: i>j)
void compute_quadri(Gauge_Conf  *GC,
            long r,
            int i,
            int j)
   {
   #ifdef THETA_MODE
   GAUGE_GROUP aux;

   quadrifoglio(GC, r, i, j, &aux);

   equal(&(GC->quadri[r][i][j]), &aux);            
   minus_equal_dag(&(GC->quadri[r][i][j]), &aux);  // quadri[r][i][j]=aux-(aux)^{dag}

   equal(&(GC->quadri[r][j][i]), &(GC->quadri[r][i][j]));
   times_equal_real(&(GC->quadri[r][j][i]),-1.0);  // quadri[r][j][i]=-quadri[r][i][j]
   #else
   // to avoid compile time warning on unused variables
   (void) GC;
   (void) r;
   (void) i;
   (void) j;
   #endif
   }
*/

// perform a complete update
void update_w(Gauge_Conf * restrict GC,
              Geometry const * const restrict geo,
              GParam const * const restrict param)
   {
   long r;
   int j, dir;

   // heatbath
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif 
      for(r=0; r<(param->d_volume)/2; r++)
         {
         heatbath_w(GC, geo, param, r, dir);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif 
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         heatbath_w(GC, geo, param, r, dir);
         } 
      }

   // overrelax
   for(dir=0; dir<STDIM; dir++)
      {
      for(j=0; j<param->d_overrelax; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=0; r<(param->d_volume)/2; r++)
            {
            overrelaxation_w(GC, geo, param, r, dir);
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            overrelaxation_w(GC, geo, param, r, dir);
            }
         }
      }
   
   // final unitarization
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

   GC->update_index++;
   }

/*
// perform n cooling steps minimizing the action at theta=0
void cooling(Gauge_Conf *GC,
             GParam const * const param,
             int n)
   {
   long r;
   int i, k;

   for(k=0; k<n; k++)
      {
      // cooling
      for(i=0; i<4; i++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=0; r<(param->d_volume)/2; r++)
            {
            GAUGE_GROUP staple;
            calcstaples_notopo(GC, r, i, &staple);
            cool(&(GC->lattice[r][i]), &staple);  
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            GAUGE_GROUP staple;
            calcstaples_notopo(GC, r, i, &staple);
            cool(&(GC->lattice[r][i]), &staple);  
            }
         }
      }

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, i)
   #endif 
   for(r=0; r<(param->d_volume); r++)
      {
      for(i=0; i<4; i++)
         {
         unitarize(&(GC->lattice[r][i]));
         } 
      }
   }
*/

#endif
