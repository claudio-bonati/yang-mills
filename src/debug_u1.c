#ifndef DEBUG_U1_C
#define DEBUG_U1_C

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/macro.h"
#include"../include/random.h"
#include"../include/u1.h"
#include"../include/u1_upd.h"

int main(void)
  {
  unsigned int seme=0;
  double energy, beta;

  U1 M, N, L, T, mI;

  // initialize random seed
  initrand(seme);

  // fix a value for d_beta
  beta=0.9;

  printf("\n******************************\n");
  printf("PROGRAM FOR THE DEBUG OF U(1)\n");
  printf("******************************\n");

  printf("\n");
  printf("VERIFY THAT THE RANDOM MATRIX IS IN U(1)\n\n");
  printf("  random matrix ....");
  rand_matrix_U1(&M);
  one_U1(&mI);
  times_equal_real_U1(&mI, -1.0);  // mI=-1.0
  times_dag2_U1(&T, &M, &M);
  plus_equal_U1(&T, &mI);
  if(norm_U1(&T) <=MIN_VALUE)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }



  printf("\n\n");
  printf("VERIFY THAT UPDATE U(1)->U(1)\n\n");
  rand_matrix_U1(&M);
  rand_matrix_U1(&N);
  rand_matrix_U1(&L);
  plus_equal_U1(&N, &L); // N+=L,  M in U(1), N no   (M=link, N=staple)

  // heatbath
  times_equal_real_U1(&N, beta);
  single_heatbath_U1(&M, &N);
  printf("  Heatbath ...");
  times_dag2_U1(&T, &M, &M);
  plus_equal_U1(&T, &mI);
  if(norm_U1(&T) <=MIN_VALUE)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }

  // overrelaxation
  rand_matrix_U1(&M);
  rand_matrix_U1(&N);
  rand_matrix_U1(&L);
  plus_equal_U1(&N, &L); // N+=L,  M in U(1), N no   (M=link, N=staple)
  single_overrelaxation_U1(&M, &N);
  printf("  Overrelaxation ...");
  times_dag2_U1(&T, &M, &M);
  plus_equal_U1(&T, &mI);
  if(norm_U1(&T) <=MIN_VALUE)
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
  rand_matrix_U1(&M);
  rand_matrix_U1(&N);
  rand_matrix_U1(&L);
  plus_equal_U1(&N, &L); // N+=L,  M in U(1), N no   (M=link, N=staple)

  times_U1(&T, &M, &N);  // T=M*N
  energy=retr_U1(&T);    // initial energy
  single_overrelaxation_U1(&M, &N);
  times_U1(&T, &M, &N);  // T=M*N
  energy-=retr_U1(&T);    // -=final energy
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
  rand_matrix_U1(&M);
  rand_matrix_U1(&N);
  rand_matrix_U1(&L);
  plus_equal_U1(&N, &L); // N+=L,  M in U(1), N no   (M=link, N=staple)

  times_U1(&T, &M, &N);  // T=M*N
  energy=retr_U1(&T);    // initial energy

  cool_U1(&M, &N);

  times_U1(&T, &M, &N);  // T=M*N
  energy-=retr_U1(&T);   // -=final energy
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
