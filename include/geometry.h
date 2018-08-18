#ifndef GEOMETRY_H
#define GEOMETRY_H

#include"macro.h"
#include"gparam.h"

typedef struct Geometry {
   long **d_nnp;      // d_nnp_loc[r][i] = next neighbour (on the local lattice) in dir.  i of the site r
   long **d_nnm;      // d_nnm_loc[r][i] = next neighbour (on the local lattice) in dir. -i of the site r
   int *d_timeslice;  // d_timeslice[r]  = time component of r
   long *d_spacecomp; // d_spacecomp[r]  = space component of r
   long **d_tsp;      // d_tsp[t][rsp] = r such that d_timeslice[r]=t and d_spacecomp[r]=rsp
} Geometry;

// these are the functions to be used in shwitching between different indices
long (*cart_to_si)(int const * const cartcoord, GParam const * const param); // cartesian coordinates -> single index
void (*si_to_cart)(int *cartcoord, long si, GParam const * const param);     // single index -> cartesian coordinates
long (*lex_to_si)(long lex, GParam const * const param);          // lexicographic -> single index
long (*si_to_lex)(long si, GParam const * const param);           // lexicographic -> single index
long (*sisp_and_t_to_si_compute)(long sisp, int t, GParam const * const param);            // single index spatial and time -> single index tot
void (*si_to_sisp_and_t_compute)(long *sisp, int *t, long si, GParam const * const param); // single index tot -> single index spatial and time

// general functions
void init_indexing_lexeo(void); // has to be called before init_geometry
void init_geometry(Geometry *geo, GParam const * const param);
void free_geometry(Geometry *geo, GParam const * const param);

// next neighbour in + direction
inline long nnp(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnp[r][i];
  }

// next neighbour in - direction
inline long nnm(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnm[r][i];
  }

// single index spatial and time -> single index tot
inline long sisp_and_t_to_si(Geometry const * const geo, long sisp, int t)
  {
  return geo->d_tsp[t][sisp];
  }

// single index tot -> single index spatial and time
inline void si_to_sisp_and_t(long *sisp, int *t, Geometry const * const geo, long si)
  {
  *sisp=geo->d_spacecomp[si];
  *t=geo->d_timeslice[si];
  }

// for debug
void test_geometry(Geometry const * const geo, GParam const * const param);

//------------ these are not to be used outside geometry.c ----------------

long cart_to_lex(int const * const cartcoord, GParam const * const param);   // cartesian coordinates -> lexicographic index
void lex_to_cart(int *cartcoord, long lex, GParam const * const param);      // lexicographic index -> cartesian coordinates

long cart_to_lexeo(int const * const cartcoord, GParam const * const param); // cartesian coordinates -> lexicographic eo index
void lexeo_to_cart(int *cartcoord, long lexeo, GParam const * const param);  // lexicographic eo index -> cartesian coordinates

long lex_to_lexeo(long lex, GParam const * const param);                     //  lexicographic index -> lexicographic eo index
long lexeo_to_lex(long lexeo, GParam const * const param);                   //  lexicographic eo index -> lexicographic index

long cartsp_to_lexsp(int const * const ccsp, GParam const * const param); // spatial cartesian coordinates -> spatial lexicographic index
void lexsp_to_cartsp(int *ccsp, long lexsp, GParam const * const param);  // spatial lexicographic index -> spatial cartesian coordinates

long cartsp_to_lexeosp(int const * const ccsp, GParam const * const param);  // spatial cartesian coordinates -> spatial lexicographic eo index
void lexeosp_to_cartsp(int *ccsp, long lexeosp, GParam const * const param); // spatial lexicographic eo index -> spatial cartesian coordinates

long lexsp_to_lexeosp(long lexsp, GParam const * const param);     //  spatial lexicographic index -> spatial lexicographic eo index
long lexeosp_to_lexsp(long lexeosp, GParam const * const param);   //  spatial lexicographic eo index -> spatial lexicographic index

long lexeosp_and_t_to_lexeo(long lexeosp, int t, GParam const * const param);    // lexicographic eo spatial and time -> lexicographic eo index
void lexeo_to_lexeosp_and_t(long *lexeosp, int *t, long lexeo, GParam const * const param); // lex. eo index -> lex. eo spatial and t

#endif
