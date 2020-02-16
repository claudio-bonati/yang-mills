#ifndef DEBUG_SON_C
#define DEBUG_SON_C

#include<math.h>
#include<stdlib.h>

#include"../include/macro.h"
#include"../include/random.h"
#include"../include/son.h"
#include"../include/son_upd.h"


int main(void)
  {
  unsigned int seme=0;
  double energy, beta;

  SoN M, N, L, T;

  // initialize random seed
  initrand(seme);

  // fix a value for d_beta
  beta=6.0;

  printf("\n*******************************\n");
  printf("PROGRAM FOR THE DEBUG OF SO(N)\n");
  printf("*******************************\n\n");


  printf("N=%s", QUOTEME(NCOLOR));


  printf("\n\n");
  printf("VERIFY THAT THE RANDOM MATRIX IS IN SO(N)\n\n");
  printf("  random matrix ....");
  rand_matrix_SoN(&M);
  if(scheck_SoN(&M) == 0)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }

  printf("\n\n");
  printf("VERIFY THAT UPDATE SO(N)->SO(N)\n\n");
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

  // heatbath
  times_equal_real_SoN(&N, beta);
  single_heatbath_SoN(&M, &N);
  printf("  Heatbath ...");
  if(scheck_SoN(&M) == 0)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }

  // overrelaxation
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple) */
  single_overrelaxation_SoN(&M, &N);
  printf("  Overrelaxation ...");
  if(scheck_SoN(&M) == 0)
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
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

  times_SoN(&T, &M, &N);  // T=M*N
  energy=retr_SoN(&T);    // initial energy
  single_overrelaxation_SoN(&M, &N);
  times_SoN(&T, &M, &N);  // T=M*N
  energy-=retr_SoN(&T);    // -=final energy
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
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); // N+=L,  M in SU(N), N no   (M=link, N=staple)

  times_SoN(&T, &M, &N);  // T=M*N
  energy=retr_SoN(&T);    // initial energy

  cool_SoN(&M, &N);

  times_SoN(&T, &M, &N);  // T=M*N
  energy-=retr_SoN(&T);   // -=final energy
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

