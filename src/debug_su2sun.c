#ifndef DEBUG_SU2SUN_C
#define DEBUG_SU2SUN_C

#include<math.h>
#include<stdlib.h>

#include"../include/flavour_matrix.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/su2.h"
#include"../include/sun.h"


void convert_Su2Vecs_to_SuNVecs(SuNVecs * restrict vN, Su2Vecs const * const restrict v2)
  {
  int i;

  for(i=0; i<NHIGGS; i++)
     {
     vN->comp[2*i+0]=v2->comp[4*i+2]+I*v2->comp[4*i+1];
     vN->comp[2*i+1]=v2->comp[4*i+0]-I*v2->comp[4*i+3];
     }
  }

void convert_Su2_to_SuN(SuN * restrict MN, Su2 const * const restrict M2)
  {
  #if NCOLOR!=2
    (void) MN;
    (void) M2;

    printf("This code require N_c=2\n\n");
    exit(EXIT_FAILURE);
  #else
    MN->comp[m(0,0)]=M2->comp[0]+I*M2->comp[3];
    MN->comp[m(1,1)]=M2->comp[0]-I*M2->comp[3];
  
    MN->comp[m(0,1)]=+M2->comp[2]+I*M2->comp[1];
    MN->comp[m(1,0)]=-M2->comp[2]+I*M2->comp[1];
  #endif
  }


int main(void)
  {
  unsigned int seme=0;

  int i;
  double x;
  double complex y;
  Su2 M2, U2;
  SuN MN, UN;

  Su2Vecs v2, w2;
  SuNVecs vN, wN;

  FMatrix fm2, fmN, fmaux;

  // initialize random seed
  initrand(seme);

  printf("\n************************************************\n");
  printf("PROGRAM TO VERIFY CONSISTECY BETWEEN Su2 and SuN\n");
  printf("************************************************\n\n");

  #if NCOLOR!=2
    printf("This code require N_c=2\n\n");
    return EXIT_FAILURE;
  #else
    printf("NCOLOR=%s\n", QUOTEME(NCOLOR));
  #endif
  printf("NHIGGS=%s\n\n", QUOTEME(NHIGGS));


  printf("VERIFY THAT Su2->SuN IS UNITARY ...");
  rand_matrix_Su2(&M2);
  convert_Su2_to_SuN(&MN, &M2);
  i=scheck_SuN(&MN);
  if(i==0)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("VERIFY CONSISTENCY BETWEEN Su2*Su2 and SuN*SuN ...");
  rand_matrix_Su2(&M2);
  rand_matrix_Su2(&U2);
  convert_Su2_to_SuN(&MN, &M2);
  convert_Su2_to_SuN(&UN, &U2);

  times_equal_Su2(&U2, &M2);
  times_equal_SuN(&UN, &MN);

  convert_Su2_to_SuN(&MN, &U2);
  minus_equal_SuN(&UN, &MN);
  x=norm_SuN(&UN);
  if(fabs(x)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("VERIFY THAT Su2Vecs->SuNVecs IS UNITARY ...");
  rand_vec_Su2Vecs(&v2);
  convert_Su2Vecs_to_SuNVecs(&vN, &v2);
  x=norm_SuNVecs(&vN);
  if(fabs(x-1)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("VERIFY CONSISTENCY BETWEEN Su2Vecs -> FMatrix and SuNVecs -> FMatrix ...");
  rand_vec_Su2Vecs(&v2);
  convert_Su2Vecs_to_SuNVecs(&vN, &v2);

  init_FMatrix_Su2Vecs(&fm2, &v2);
  init_FMatrix_SuNVecs(&fmN, &vN);

  minus_equal_FMatrix(&fm2, &fmN);
  x=norm_FMatrix(&fm2);
  if(fabs(x)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("VERIFY THAT FMatrix IS TRACELESS HERMITIAN ...");
  rand_vec_Su2Vecs(&v2);
  rand_vec_SuNVecs(&vN);

  init_FMatrix_Su2Vecs(&fm2, &v2);
  init_FMatrix_SuNVecs(&fmN, &vN);

  x=fabs(retr_FMatrix(&fm2))+fabs(retr_FMatrix(&fmN));

  equal_FMatrix(&fmaux, &fm2);
  minus_equal_dag_FMatrix(&fmaux, &fm2);
  x+=norm_FMatrix(&fmaux);

  equal_FMatrix(&fmaux, &fmN);
  minus_equal_dag_FMatrix(&fmaux, &fmN);
  x+=norm_FMatrix(&fmaux);

  if(fabs(x/4.0)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("VERIFY CONSISTENCY BETWEEN Su2*Su2Vecs and SuN*SuNVecs ...");
  rand_vec_Su2Vecs(&v2);
  convert_Su2Vecs_to_SuNVecs(&vN, &v2);

  rand_matrix_Su2(&M2);
  convert_Su2_to_SuN(&MN, &M2);

  matrix_times_vector_all_Su2Vecs(&w2, &M2, &v2);
  matrix_times_vector_all_SuNVecs(&wN, &MN, &vN);

  convert_Su2Vecs_to_SuNVecs(&vN, &w2);
  minus_equal_SuNVecs(&vN, &wN);
  x=norm_SuNVecs(&vN);
  if(fabs(x)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("VERIFY CONSISTENCY BETWEEN vector_tensor_vector_Su2Vecs and vector_tensor_vector_SuNVecs ...");
  rand_vec_Su2Vecs(&v2);
  rand_vec_Su2Vecs(&w2);
  convert_Su2Vecs_to_SuNVecs(&vN, &v2);
  convert_Su2Vecs_to_SuNVecs(&wN, &w2);

  vector_tensor_vector_Su2Vecs(&M2, &v2, &w2);
  vector_tensor_vector_SuNVecs(&MN, &vN, &wN);

  rand_matrix_Su2(&U2);
  convert_Su2_to_SuN(&UN, &U2);

  times_equal_Su2(&U2, &M2);
  times_equal_SuN(&UN, &MN);

  x=retr_Su2(&U2)-retr_SuN(&UN);
  if(fabs(x)<MIN_VALUE)
    {
    printf("  OK\n");
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");


  printf("VERIFY CONSISTENCY BETWEEN Su2Vecs -> HiggsU1Obs and SuNVecs -> HiggsU1Obs ...");
  rand_vec_Su2Vecs(&v2);
  convert_Su2Vecs_to_SuNVecs(&vN, &v2);

  y=HiggsU1Obs_Su2Vecs(&v2);
  y-=HiggsU1Obs_SuNVecs(&vN);

  x=cabs(y);
  if(fabs(x)<MIN_VALUE)
    {
    if(NHIGGS<NCOLOR)
      {
      printf("  OK but trivial since NHIGGS<NCOLOR\n");
      }
    else
      {
      printf("  OK\n");
      }
    }
  else
    {
    printf("  ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }
  printf("\n");

  printf("\nTEST PASSED\n\n");

  return EXIT_SUCCESS;
  }

#endif
