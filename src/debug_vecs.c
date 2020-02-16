#ifndef DEBUG_VECS_C
#define DEBUG_VECS_C

#include<complex.h>
#include<math.h>
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/macro.h"
#include"../include/random.h"


int main(void)
  {
  unsigned int seme=0;
  int i;
  double energy;

  GAUGE_VECS M, N, L;
  GAUGE_GROUP matrix;

  // initialize random seed
  initrand(seme);

  printf("\n***************************************************\n");
  printf("PROGRAM FOR THE DEBUG OF VECS [for U1, Su2 and SuN] \n");
  printf("***************************************************\n\n");

  printf("GGROUP=%s\n", QUOTEME(GGROUP));
  printf("NCOLOR=%s\n", QUOTEME(NCOLOR));
  printf("NHIGGS=%s\n\n", QUOTEME(NHIGGS));


  printf("VERIFY THAT THE RANDOM VECTOR IS CORRECTLY NORMALIZED ...");
  rand_vecs(&M);
  equal_vecs(&N, &M);
  energy=re_scal_prod_vecs(&M,&N);
  if(fabs(energy-1.0) < MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");

  printf("VERIFY THAT THE RANDOM ROTATION DOES NOT CHANGE THE NORMALIZATION ...");
  rand_vecs(&M);
  rand_matrix(&matrix);
  i=(int) (NHIGGS*casuale()-MIN_VALUE);
  matrix_times_vector_single_vecs(&L, &matrix, &M, i);
  equal_vecs(&N, &L);
  energy=re_scal_prod_vecs(&L,&N);
  if(fabs(energy-1) < MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");

  printf("VERIFY THAT THE RANDOM PHASE MULTIPLICATION DOES NOT CHANGE THE NORMALIZATION ...");
  rand_vecs(&L);
  i=(int) (NHIGGS*casuale()-MIN_VALUE);
  times_equal_complex_single_vecs(&L, cexp(PI2*casuale()*I), i);
  equal_vecs(&N, &L);
  energy=re_scal_prod_vecs(&L, &N);
  if(fabs(energy-1) < MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");

  printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
  rand_vecs(&M);
  equal_vecs(&L, &M);
  rand_vecs(&N);
  times_equal_real_vecs(&N, casuale());

  energy=re_scal_prod_vecs(&M, &N);

  single_overrelaxation_vecs(&M, &N);

  energy-=re_scal_prod_vecs(&M, &N);

  if(fabs(energy)<MIN_VALUE)
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

