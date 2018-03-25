#ifndef MYMALLOC_C
#define MYMALLOC_C

#include"../include/macro.h"

#ifdef MEMALIGN_MODE
  #include<malloc.h>
#endif
#include<stdlib.h>

void * mymalloc(size_t alignement, size_t size)
  {
  #ifdef MEMALIGN_MODE
  void *ptr=memalign(alignement, size);
  return ptr;
  #else
  (void) alignement; // just to avoid warning at compile time
  void *ptr=malloc(size);
  return ptr;
  #endif
  }

#endif

