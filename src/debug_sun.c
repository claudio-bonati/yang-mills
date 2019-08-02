#ifndef DEBUG_SUN_C
#define DEBUG_SUN_C

#include<math.h>
#include<stdlib.h>

#include"../include/macro.h"
#include"../include/random.h"
#include"../include/sun.h"
#include"../include/sun_upd.h"


int main(void)
  {
  unsigned int seme=0;
  double energy, beta;
  
  SuN M, N, L, T;
   
  // initialize random seed
  initrand(seme);

  // fix a value for d_beta
  beta=6.0;
    
  printf("\n*******************************\n");
  printf("PROGRAM FOR THE DEBUG OF SU(N)\n");
  printf("*******************************\n\n");


  printf("N=%s", QUOTEME(NCOLOR));


  printf("\n\n");
  printf("VERIFY THAT THE RANDOM MATRIX IS IN SU(N)\n\n");
  printf("  random matrix ....");
  rand_matrix_SuN(&M);
  if(scheck_SuN(&M) == 0)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }

  printf("\n\n");
  printf("VERIFY THAT UPDATE SU(N)->SU(N)\n\n");
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

  // heatbath
  times_equal_real_SuN(&N, beta);
  single_heatbath_SuN(&M, &N);
  printf("  Heatbath ...");
  if(scheck_SuN(&M) == 0)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }

  // overrelaxation
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple) */
  single_overrelaxation_SuN(&M, &N);
  printf("  Overrelaxation ...");
  if(scheck_SuN(&M) == 0)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }



  printf("\n\n");
  printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

  times_SuN(&T, &M, &N);  // T=M*N
  energy=retr_SuN(&T);    // initial energy
  single_overrelaxation_SuN(&M, &N);
  times_SuN(&T, &M, &N);  // T=M*N
  energy-=retr_SuN(&T);    // -=final energy
  if(fabs(energy)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!   DeltaE=%g\n", energy);
    return EXIT_FAILURE;
    }



  printf("\n\n");
  printf("VERIFY THAT COOLING DECREASES THE ENERGY ...");
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

  times_SuN(&T, &M, &N);  // T=M*N
  energy=retr_SuN(&T);    // initial energy

  cool_SuN(&M, &N);

  times_SuN(&T, &M, &N);  // T=M*N
  energy-=retr_SuN(&T);   // -=final energy
  if(energy<MIN_VALUE)

    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!   DeltaE=%g\n", energy);
    return EXIT_FAILURE;
    }

  printf("\nTEST PASSED\n\n");

  return EXIT_SUCCESS;
  }

#endif
