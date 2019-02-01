#ifndef DEBUG_SU2_ADJ_C
#define DEBUG_SU2_ADJ_C

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/macro.h"
#include"../include/random.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"

#define m2adj(X,Y) ((X)*3 + (Y))

int main()
  {
  unsigned int seme=0;

  Su2 B1, B2;
  Su2Adj A1;

  int i;
  double risth, ris;

  // initialize random seed
  initrand(seme);

  printf("\n*********************************\n");
  printf("PROGRAM FOR THE DEBUG OF SU(2)Adj\n");
  printf("*********************************\n");

  printf("\n");
  printf("VERIFY THAT Tr(Adj)=|Tr(Fund)|^2-1\n\n");
  printf("  ....");

  //random matrix in the fundamental repr.
  rand_matrix_Su2(&B1);
  rand_matrix_Su2(&B2);

  //compute matrices (B1,B2) in adjoint representation (A1,A2)
  fund_to_adj_Su2(&A1, &B1);

  // shoud be trace
  risth=pow(2.0*retr_Su2(&B1),2.0)-1;

  ris=0.0;
  for(i=0; i<2*2-1; i++)
     {
     ris+=A1.comp[m2adj(i,i)];
     }

  if(fabs(risth-ris) <=MIN_VALUE)
    {
    printf("    OK\n");
    }
  else
    {
    printf("    ERROR!!!!!!!!!!!\n");
    return EXIT_FAILURE;
    }

  printf("\nTEST PASSED\n\n");

  return EXIT_SUCCESS;
  }

#undef m2adj

#endif
