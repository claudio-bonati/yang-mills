#ifndef DEBUG_RNG_C
#define DEBUG_RNG_C

#include<stdio.h>
#include<stdlib.h>

#include"../include/macro.h"
#include"../include/random.h"

int main(void)
   {
   int i;
   unsigned int seme;

   printf("\n");

   seme=1;
   initrand(seme);
   printf("seed=%u\n", seme);
   for(i=0; i<5; i++)
      {
      printf("  random[%d] = %.16g\n", i, casuale());
      }

   printf("\n");

   seme=2;
   initrand(seme);
   printf("seed=%u\n", seme);
   for(i=0; i<5; i++)
      {
      printf("  random[%d] = %.16g\n", i, casuale());
      }

   printf("\n");

   seme=1;
   initrand(seme);
   printf("seed=%u\n", seme);
   for(i=0; i<5; i++)
      {
      printf("  random[%d] = %.16g\n", i, casuale());
      }

   printf("\n");

   seme=0;
   initrand(seme);
   printf("seed=time()\n");
   for(i=0; i<5; i++)
      {
      printf("  random[%d] = %.16g\n", i, casuale());
      }

   printf("\n");

   return EXIT_SUCCESS;
   }
#endif
