#ifndef MACRO_H
#define MACRO_H

#include"../config.h"

#if GGROUP == 0
  #if NCOLOR == 1
    #define GAUGE_GROUP     U1
    #define GAUGE_GROUP_ADJ U1Adj
    #define GAUGE_VECS      U1Vecs
  #elif NCOLOR == 2
    #define GAUGE_GROUP     Su2
    #define GAUGE_GROUP_ADJ Su2Adj
    #define GAUGE_VECS      Su2Vecs
  #else
    #define GAUGE_GROUP     SuN
    #define GAUGE_GROUP_ADJ SuNAdj
    #define GAUGE_VECS      SuNVecs
  #endif
#elif GGROUP == 1
  #define GAUGE_GROUP     SoN
  #define GAUGE_GROUP_ADJ SoNAdj
  #define GAUGE_VECS      SoNVecs

  #if NCOLOR == 1
    #error N_c==1 incompatible with SoN
  #endif
#endif

// function to access matrix elements
#define m(X,Y) ((X)*NCOLOR + (Y))  // for gauge group
#define mf(X,Y) ((X)*NHIGGS + (Y)) // for higgs field
#if GGROUP == 0
  #define madj(X,Y) ((X)*(NCOLOR*NCOLOR -1) + (Y))      // for the adjoint rep of SuN
#elif GGROUP == 1
  #define madj(X,Y) ((X)*(NCOLOR*(NCOLOR -1)/2) + (Y))  // for the adjoint rep of SoN
#endif


#define MIN_VALUE 1.0e-13

#define INT_ALIGN 16
#define DOUBLE_ALIGN 32

static const double PI=3.141592653589793238462643383279502884197169399375105820974944;
static const double PI2=6.283185307179586476925286766559005768394338798750211641949889;
static const double HALF_PI=1.570796326794896619231321691639751442098584699687552910487472;

#define STD_STRING_LENGTH 50 // standarg lenght of unknown strings

// way to print a macro: if
// #define val1 val2
// then QUOTEME(val1) give the string "val2"
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

// to activate posix_memalign in stdlib.h
#ifndef _POSIX_C_SOURCE 
#define _POSIX_C_SOURCE 200809L
#endif

#endif
