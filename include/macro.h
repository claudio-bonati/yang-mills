#ifndef MACRO_H
#define MACRO_H

#include"../config.h"

#if NCOLOR == 2
  #define GAUGE_GROUP     Su2
  #define GAUGE_VECS      Su2Vecs
#else
  #define GAUGE_GROUP     SuN
  #define GAUGE_VECS      SuNVecs
#endif

// function to access matrix elements
#define m(X,Y) ((X)*NCOLOR + (Y))  // for gauge group
#define mf(X,Y) ((X)*NHIGGS + (Y)) // for higgs field

#define MIN_VALUE 1.0e-13

#define INT_ALIGN 16
#define DOUBLE_ALIGN 32

static const double PI=3.141592653589793238462643383279502884197169399375105820974944;
static const double PI2=6.283185307179586476925286766559005768394338798750211641949889;
static const double HALF_PI=1.570796326794896619231321691639751442098584699687552910487472;
static const double PI2_N = PI2 / (double)NCOLOR;


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

#define OPENSSL_SUPPRESS_DEPRECATED

#endif
