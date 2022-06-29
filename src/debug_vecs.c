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
  int i, k;
  double energy;

  GAUGE_VECS M, N, L;
  GAUGE_GROUP matrix;

  // initialize random seed
  initrand(seme);

  printf("\n***************************************************\n");
  printf("PROGRAM FOR THE DEBUG OF VECS \n");
  printf("***************************************************\n\n");

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

  if(NHIGGS>1)
    {
    printf("VERIFY THAT THE TWO COMPONENT ROTATION DOES NOT CHANGE THE NORMALIZATION ...");
    rand_vecs(&M);
    i=(int) (NHIGGS*casuale()-MIN_VALUE);
    k=(i+1 + (int)((NHIGGS-1)*casuale()*(1.0 - MIN_VALUE)) )% NHIGGS;
    rotate_two_components_vecs(&L, &M, i, k, PI2*casuale());
    energy=norm_vecs(&L);
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
    }

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

  printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE ENERGY AND NORMALIZATION ...");
  rand_vecs(&M);
  equal_vecs(&L, &M);
  rand_vecs(&N);
  times_equal_real_vecs(&N, casuale());

  energy=re_scal_prod_vecs(&M, &N);

  single_overrelaxation_vecs(&M, &N);

  energy-=re_scal_prod_vecs(&M, &N);

  if(fabs(energy)<MIN_VALUE)
    {
    printf("  OK");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!   DeltaE=%g  ", energy);
    return EXIT_FAILURE;
    }

  energy=norm_vecs(&M);
  if(fabs(energy-1)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");

  printf("VERIFY THAT OVERRELAXATION IS REVERSIBLE ...");
  rand_vecs(&M);
  equal_vecs(&L, &M);
  rand_vecs(&N);
  times_equal_real_vecs(&N, casuale());

  single_overrelaxation_vecs(&M, &N);
  single_overrelaxation_vecs(&M, &N);

  times_equal_real_vecs(&M, -1);
  plus_equal_vecs(&M, &L);

  energy=norm_vecs(&M);
  if(fabs(energy)<MIN_VALUE)
    {
    printf("  OK");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!   Delta=%g  ", energy);
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("\nTEST PASSED\n\n");

  return EXIT_SUCCESS;
  }

#endif

