#ifndef RANDOM_C
#define RANDOM_C

#include"../include/macro.h"

#include<limits.h>
#include<math.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<time.h>

#include"../include/dSFMT.h"
#include"../include/random.h"


// this is the status of the random number generator
dsfmt_t rng_status[NTHREADS];


// random number in (0,1)
double casuale(void)
 {
 double ris;
 int thread_id=0;
 #ifdef OPENMP_MODE
   thread_id=omp_get_thread_num();
 #endif

 ris=dsfmt_genrand_open_open(&(rng_status[thread_id]));

 return ris;
 }


// initialize random generator
void initrand(unsigned int s)
  {
  int thread_id;
  unsigned int seed=s;

  if(s==0)
    {
    seed=((unsigned int) time(NULL)+10) % UINT_MAX;
    if(seed==0)
      {
      seed=1;
      }
    }

  for(thread_id=0; thread_id<NTHREADS; thread_id++)
     {
     dsfmt_init_gen_rand(&(rng_status[thread_id]), seed);
     seed=(seed+10) % UINT_MAX;
     if(seed==0)
       {
       seed=1;
       }
     }
  }


// normal gaussian random number generator (polar method, knuth vol 2, p. 117)
double gauss1()
   {
   double v1, v2, s, ris;

   do
     {
     v1=1.0-2.0*casuale();
     v2=1.0-2.0*casuale();
     s=v1*v1+v2*v2;
     }
   while(s >= 1);

   ris=v1*sqrt(-2*log(s)/s);
   return ris;
   }


// normal gaussian random number generator (polar method, knuth vol 2, p. 117)
void gauss2(double *ris1, double *ris2)
   {
   double v1, v2, s;

   do
     {
     v1=-1.0+2.0*casuale();
     v2=-1.0+2.0*casuale();
     s=v1*v1+v2*v2;
     }
   while(s >= 1);

   *ris1=v1*sqrt(-2*log(s)/s);
   *ris2=v2*sqrt(-2*log(s)/s);
   }


#endif
