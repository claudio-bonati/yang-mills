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

// compute the staple in position r, direction i and save it in M
void calcstaples_wilson(Gauge_Conf const * const restrict GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i,
                        GAUGE_GROUP * restrict M)
   {
   int j, l;
   long k;
   GAUGE_GROUP link1, link2, link3, link12, stap;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= param->d_stdim)
     {
     fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, param->d_stdim, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   zero(M); // M=0

   for(l=i+1; l< i + param->d_stdim; l++)
      {
      j = (l % param->d_stdim);

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
void heatbath(Gauge_Conf *GC,
                Geometry const * const geo,
                GParam const * const param,
                long r,
                int i)
   {
   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= param->d_stdim)
     {
     fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, param->d_stdim, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   GAUGE_GROUP stap;

   calcstaples_wilson(GC, geo, param, r, i, &stap);
   single_heatbath(&(GC->lattice[r][i]), &stap, param);
   }


// perform an update with overrelaxation
void overrelaxation(Gauge_Conf *GC,
                      Geometry const * const geo,
                      GParam const * const param,
                      long r,
                      int i)
   {
   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= param->d_stdim)
     {
     fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, param->d_stdim, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif
   (void) param; // just to avoid wornings

   GAUGE_GROUP stap;

   calcstaples_wilson(GC, geo, param, r, i, &stap);
   single_overrelaxation(&(GC->lattice[r][i]), &stap);
   }

/*
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
void update(Gauge_Conf * GC,
            Geometry const * const geo,
            GParam const * const param)
   {
   long r;
   int j, dir;

   #ifdef OPENMP_MODE
   if(geo->indexing_type!=0)
     {
     fprintf(stderr, "Wrong indexing used! (indexing_type=%d) (%s, %d)\n", geo->indexing_type, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   // heatbath
   for(dir=0; dir<param->d_stdim; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif 
      for(r=0; r<(param->d_volume)/2; r++)
         {
         heatbath(GC, geo, param, r, dir);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif 
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         heatbath(GC, geo, param, r, dir);
         } 
      }

   // overrelax
   for(dir=0; dir<param->d_stdim; dir++)
      {
      for(j=0; j<param->d_overrelax; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=0; r<(param->d_volume)/2; r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }
         }
      }
   
   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif 
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<param->d_stdim; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         } 
      }

   GC->update_index++;
   }


// perform n cooling steps minimizing the action at theta=0
void cooling(Gauge_Conf *GC,
             Geometry const * const geo,
             GParam const * const param,
             int n)
   {
   long r;
   int i, k;

   #ifdef OPENMP_MODE
   if(geo->indexing_type!=0)
     {
     fprintf(stderr, "Wrong indexing used! (indexing_type=%d) (%s, %d)\n", geo->indexing_type, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   for(k=0; k<n; k++)
      {
      // cooling
      for(i=0; i<param->d_stdim; i++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=0; r<(param->d_volume)/2; r++)
            {
            GAUGE_GROUP staple;
            calcstaples_wilson(GC, geo, param, r, i, &staple);
            cool(&(GC->lattice[r][i]), &staple);  
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            GAUGE_GROUP staple;
            calcstaples_wilson(GC, geo, param, r, i, &staple);
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
      for(i=0; i<param->d_stdim; i++)
         {
         unitarize(&(GC->lattice[r][i]));
         } 
      }
   }


// perform a single step of the Runge Kutta integrator for the Wilson flow
// as described in Luscher arXiv:1006.4518 app. C
void gradflow_RKstep(Gauge_Conf *GC,
                     Gauge_Conf *helper1,
                     Gauge_Conf *helper2,
                     Geometry const * const geo,
                     GParam const *const param,
                     double dt)
   {
   long r;
   int dir;

   #ifdef OPENMP_MODE
   if(geo->indexing_type!=0)
     {
     fprintf(stderr, "Wrong indexing used! (indexing_type=%d) (%s, %d)\n", geo->indexing_type, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   // initialize
   for(dir=0; dir<param->d_stdim; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_volume; r++)
         {
         equal(&(helper1->lattice[r][dir]), &(GC->lattice[r][dir]));
         equal(&(helper2->lattice[r][dir]), &(GC->lattice[r][dir]));
         }
      }


   for(dir=0; dir<param->d_stdim; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_volume; r++)
         {
         GAUGE_GROUP staple, aux, link;

         calcstaples_wilson(helper1, geo, param, r, dir, &staple);
         equal(&link, &(helper1->lattice[r][dir]));
         times(&aux, &link, &staple);                // aux=link*staple
         times_equal_real(&aux, -dt/4.0);
         equal(&(helper2->lattice[r][dir]), &aux);    // helper2=aux
         taexp(&aux);
         times(&(GC->lattice[r][dir]), &aux, &link); // GC=aux*link
         }
      }

   // now helper1=W_0, helper2=(1/4)Z_0 and GC=W_1

   for(dir=0; dir<param->d_stdim; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_volume; r++)
         {
         GAUGE_GROUP staple, aux, link;

         calcstaples_wilson(GC, geo, param, r, dir, &staple);
         equal(&link, &(GC->lattice[r][dir]));
         times(&aux, &link, &staple);               // aux=link*staple
         times_equal_real(&aux, -dt*8.0/9.0);
         minus_equal_times_real(&aux, &(helper2->lattice[r][dir]), 17.0/9.0); // 1/4 was in Z_0
         equal(&(helper2->lattice[r][dir]), &aux);
         taexp(&aux);
         times(&(helper1->lattice[r][dir]), &aux, &link); // helper1=aux*link
         }
      }

   // now helper1=W_2, helper2=(8/9)Z_1-(17/36)Z_0 and GC=W_1

   for(dir=0; dir<param->d_stdim; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_volume; r++)
         {
         GAUGE_GROUP staple, aux, link;

         calcstaples_wilson(helper1, geo, param, r, dir, &staple);
         equal(&link, &(helper1->lattice[r][dir]));
         times(&aux, &link, &staple);                   // aux=link*staple
         times_equal_real(&aux, -dt*3.0/4.0);
         minus_equal(&aux, &(helper2->lattice[r][dir])); // aux=(3/4)Z_2-(8/9)Z_1+(17/36)Z_0
         taexp(&aux);
         times(&(GC->lattice[r][dir]), &aux, &link);    // GC=aux*link
         }
      }

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i;
      for(i=0; i<param->d_stdim; i++)
         {
         unitarize(&(GC->lattice[r][i]));
         }
      }
   }



// n step of ape smearing with parameter alpha
// new=Proj[old + alpha *staple ]
void ape_smearing(Gauge_Conf *GC,
                  Geometry const * const geo,
                  GParam const *const param,
                  double alpha,
                  int n)
   {
   Gauge_Conf helper1;
   long r;
   int dir, count;

   #ifdef OPENMP_MODE
   if(geo->indexing_type!=0)
     {
     fprintf(stderr, "Wrong indexing used! (indexing_type=%d) (%s, %d)\n", geo->indexing_type, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   init_gauge_conf_from_gauge_conf(&helper1, GC, param); //helper1=GC

   for(count=0; count<n; count++)
      {
      if(count%2==0) // smear(helper1)->GC
        {
        for(dir=0; dir<param->d_stdim; dir++)
           {
           #ifdef OPENMP_MODE
           #pragma omp parallel for num_threads(NTHREADS) private(r)
           #endif
           for(r=0; r<param->d_volume; r++)
              {
              GAUGE_GROUP staple, link;

              calcstaples_wilson(&helper1, geo, param, r, dir, &staple);
              equal(&link, &(helper1.lattice[r][dir]));
              times_equal_real(&link, 1-alpha);
              times_equal_real(&staple, alpha/6.0);
              plus_equal_dag(&link, &staple);
              unitarize(&link);
              equal(&(GC->lattice[r][dir]), &link);
              }
           }
        }
      else // smear(GC)->helper1
        {
        for(dir=0; dir<param->d_stdim; dir++)
           {
           #ifdef OPENMP_MODE
           #pragma omp parallel for num_threads(NTHREADS) private(r)
           #endif
           for(r=0; r<param->d_volume; r++)
              {
              GAUGE_GROUP staple, link;

              calcstaples_wilson(GC, geo, param, r, dir, &staple);
              equal(&link, &(GC->lattice[r][dir]));
              times_equal_real(&link, 1-alpha);
              times_equal_real(&staple, alpha/6.0);
              plus_equal_dag(&link, &staple);
              unitarize(&link);
              equal(&(helper1.lattice[r][dir]), &link);
              }
           }
        }
      }

   if(n>0 && n%2==0) // GC=helper1
     {
     for(dir=0; dir<param->d_stdim; dir++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_volume; r++)
           {
           equal(&(GC->lattice[r][dir]), &(helper1.lattice[r][dir]));
           }
        }
     }

   end_gauge_conf(&helper1, param);
   }




#endif
