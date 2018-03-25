#ifndef ENDIANNESS_C
#define ENDIANNESS_C

#include<stdlib.h>

#include"../include/endianness.h"

int endian(void)
   {
   int i = 1;
   char *p = (char *)&i;

   if (p[0] == 1)
        return 0; // LITTLE ENDIAN
   else
        return 1; // BIG ENDIAN
   }


void SwapBytesInt(void *pv)
    {
    char *p = pv;
    size_t lo, hi;
    char tmp;

    for(lo=0, hi=sizeof(int)-1; hi>lo; lo++, hi--)
       {
       tmp=p[lo];
       p[lo] = p[hi];
       p[hi] = tmp;
       }
    }

void SwapBytesFloat(void *pv)
    {
    char *p = pv;
    size_t lo, hi;
    char tmp;

    for(lo=0, hi=sizeof(float)-1; hi>lo; lo++, hi--)
       {
       tmp=p[lo];
       p[lo] = p[hi];
       p[hi] = tmp;
       }
    }

void SwapBytesDouble(void *pv)
    {
    char *p = pv;
    size_t lo, hi;
    char tmp;

    for(lo=0, hi=sizeof(double)-1; hi>lo; lo++, hi--)
       {
       tmp=p[lo];
       p[lo] = p[hi];
       p[hi] = tmp;
       }
    }


#endif
