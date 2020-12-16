#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<math.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

// compute the staple in position r, direction i and save it in M
void calcstaples_wilson(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i,
                        GAUGE_GROUP *M)
  {
  int j, l;
  long k;
  GAUGE_GROUP link1, link2, link3, link12, stap;

  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #else
  (void) param; // just to avoid warnings
  #endif

  zero(M); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

     equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
     equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

     times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
     times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

     plus_equal(M, &stap);

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |
//     (2) V           |
//         |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

     k=nnm(geo, r, j);

     equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
     equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

     times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
     times(&stap, &link12, &link3);        // stap=link12*link3

     plus_equal(M, &stap);
     }
   }


// compute the components of the staple in position r, direction i and save it in M[2*(STDIM-1)+1]
// position 0 of M is not used. It is used in simulations at imaginary theta values
void calcstaples_wilson_nosum(Gauge_Conf const * const GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              long r,
                              int i,
                              GAUGE_GROUP *M)
  {
  int j, l, count;
  long k;
  GAUGE_GROUP link1, link2, link3, link12, stap;

  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #else
  (void) param; // just to avoid warnings
  #endif

  count=0;
  zero(&M[count]); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     count++;
     zero(&M[count]); // M=0

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

     equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
     equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

     times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
     times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

     equal(&(M[count]), &stap);

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |
//     (2) V           |
//         |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

     count++;
     zero(&M[count]); // M=0

     k=nnm(geo, r, j);

     equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
     equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

     times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
     times(&stap, &link12, &link3);        // stap=link12*link3

     equal(&(M[count]), &stap);
     }

   #ifdef DEBUG
   GAUGE_GROUP helper;
   int m;

   calcstaples_wilson(GC, geo, param, r, i, &helper);
   for(m=0; m<2*(STDIM-1)+1; m++)
      {
      minus_equal(&helper, &(M[m]));
      }
   if(norm(&helper)>MIN_VALUE)
     {
     fprintf(stderr, "Problems in calcstaples_wilson_nosum (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif
   }


// compute the staple for the trace deformed theory:
// in practice a Polyakov loop without a link
void calcstaples_tracedef(Gauge_Conf const * const GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          long r,
                          int i,
                          GAUGE_GROUP * M)
  {
  if(i!=0)
    {
    zero(M);
    #ifdef DEBUG
    fprintf(stderr, "Using calcstaples_tracedef for a non-temporal link (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    #endif
    }
  else
    {
    int j;
    long int rnext;
    GAUGE_GROUP aux;

    one(&aux);

    rnext=r;
    for(j=1; j<param->d_size[0]; j++)
       {
       rnext=nnp(geo, rnext, 0);
       times_equal(&aux, &(GC->lattice[rnext][0]));
       }

    equal(M, &aux);
    }
  }


// compute all the clovers in directions ortogonal to "dir"
void compute_clovers(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     int dir)
  {
  long r;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<param->d_volume; r++)
     {
     GAUGE_GROUP aux;
     int i, j;

     for(i=0; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           if(i!=dir && j!=dir)
             {
             clover(GC, geo, param, r, i, j, &aux);

             equal(&(GC->clover_array[r][i][j]), &aux);
             minus_equal_dag(&(GC->clover_array[r][i][j]), &aux);  // clover_array[r][i][j]=aux-aux^{dag}

             equal(&(GC->clover_array[r][j][i]), &(GC->clover_array[r][i][j]));
             times_equal_real(&(GC->clover_array[r][j][i]), -1.0); // clover_array[r][j][i]=-clover_array[r][i][j]
             }
           }
        }
     }
  }


// compute the staple in position r, direction i and save it in M
// when an imaginary theta term is present
void calcstaples_with_topo(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           long r,
                           int i,
                           GAUGE_GROUP *M)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifndef THETA_MODE
  calcstaples_wilson(GC, geo, param, r, i, M);
  #else

  if(STDIM!=4)
    {
    fprintf(stderr, "Error: imaginary theta term can be used only in 4 dimensions! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  GAUGE_GROUP link1, link2, link3, link12, stap, topo_stap, aux;
  const double coeff=(param->d_theta)*((double) NCOLOR)/(param->d_beta*128.0*PI*PI);
  long k;
  int j, l;
  int i0, j0;
  int sood1[4][4], sood2[4][4]; // signed ordered orthogonal directions

  zero(M); // M=0
  zero(&topo_stap); // topo_stap=0

  // the theta term is written as
  // theta/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood1[\mu][\nu] sood2[\mu][\nu]} )
  // the independent permutations are 0123 0231 0312

  sood1[0][1] = 2;
  sood2[0][1] = 3;
  sood1[1][0] = 3;
  sood2[1][0] = 2;

  sood1[0][2] = 3;
  sood2[0][2] = 1;
  sood1[2][0] = 1;
  sood2[2][0] = 3;

  sood1[0][3] = 1;
  sood2[0][3] = 2;
  sood1[3][0] = 2;
  sood2[3][0] = 1;

  sood1[1][2] = 0;
  sood2[1][2] = 3;
  sood1[2][1] = 3;
  sood2[2][1] = 0;

  sood1[1][3] = 2;
  sood2[1][3] = 0;
  sood1[3][1] = 0;
  sood2[3][1] = 2;

  sood1[2][3] = 0;
  sood2[2][3] = 1;
  sood1[3][2] = 1;
  sood2[3][2] = 0;

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     i0=sood1[i][j];
     j0=sood2[i][j];

//
//       i ^
//         |   (1)
//     (b) +----->-----+ (c)
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//     (a) +-----<-----+-->   j
//       r     (3)    (d)
//

     equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
     equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

     times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
     times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

     plus_equal(M, &stap);

     // clover insertion in (a)
     times(&aux, &stap, &(GC->clover_array[r][i0][j0])); // stap*clover
     plus_equal(&topo_stap, &aux);

     // clover insertion in (b)
     times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap);  // clover*stap
     plus_equal(&topo_stap, &aux);

     // clover insertion in (c)
     times(&aux, &link1, &(GC->clover_array[nnp(geo, nnp(geo, r, i), j)][i0][j0]));  // link1*clover
     times_equal_dag(&aux, &link2);       // *=link2^{dag}
     times_equal_dag(&aux, &link3);       // *=link3^{dag}
     plus_equal(&topo_stap, &aux);

     // clover insertion in (d)
     times(&aux, &link12, &(GC->clover_array[nnp(geo, r, j)][i0][j0]));  // link1*link2*quadri
     times_equal_dag(&aux, &link3);          // *=link3^{dag}

     plus_equal(&topo_stap, &aux);

//
//       i ^
//         |   (1)
//     (d) +----<------+ (a)
//         |           |
//         |
//     (2) V           |
//         |
//         |           | (b)
//     (c) +------>----+--->j
//        k     (3)    r
//

     k=nnm(geo, r, j);

     equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
     equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

     times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
     times(&stap, &link12, &link3);        // stap=link12*link3

     plus_equal(M, &stap);

     // clover insertion in (a)
     times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap); // clover*stap
     minus_equal(&topo_stap, &aux);

     // clover insertion in (b)
     times(&aux, &stap, &(GC->clover_array[r][i0][j0])); // stap*clover
     minus_equal(&topo_stap, &aux);

     // clover insertion in (c)
     times(&aux, &link12, &(GC->clover_array[k][i0][j0])); // link1^{dag}*link2^{dag}*clover
     times_equal(&aux, &link3);                            // *=link3
     minus_equal(&topo_stap, &aux);

     // clover insertion in (d)
     times_dag1(&aux, &link1, &(GC->clover_array[nnp(geo, k, i)][i0][j0]));  // link1^{dag}*clover
     times_equal_dag(&aux, &link2);            // *=link2^{dag}
     times_equal(&aux, &link3);                // *=link3

     minus_equal(&topo_stap, &aux);
     }

  times_equal_real(&topo_stap, coeff);
  plus_equal(M, &topo_stap);
  #endif
  }


// compute the components of the staple in position r, direction i and save it in M[2*(STDIM-1)+1]
// position 0 of M is the topological staple
void calcstaples_with_topo_nosum(Gauge_Conf const * const GC,
                                 Geometry const * const geo,
                                 GParam const * const param,
                                 long r,
                                 int i,
                                 GAUGE_GROUP *M)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifndef THETA_MODE
  calcstaples_wilson_nosum(GC, geo, param, r, i, M);
  #else

  if(STDIM!=4)
    {
    fprintf(stderr, "Error: imaginary theta term can be used only in 4 dimensions! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  GAUGE_GROUP link1, link2, link3, link12, stap, topo_stap, aux;
  const double coeff=(param->d_theta)*((double) NCOLOR)/(param->d_beta*128.0*PI*PI);
  long k;
  int j, l, count;
  int i0, j0;
  int sood1[4][4], sood2[4][4]; // signed ordered orthogonal directions

  // the theta term is written as
  // theta/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood1[\mu][\nu] sood2[\mu][\nu]} )
  // the independent permutations are 0123 0231 0312

  sood1[0][1] = 2;
  sood2[0][1] = 3;
  sood1[1][0] = 3;
  sood2[1][0] = 2;

  sood1[0][2] = 3;
  sood2[0][2] = 1;
  sood1[2][0] = 1;
  sood2[2][0] = 3;

  sood1[0][3] = 1;
  sood2[0][3] = 2;
  sood1[3][0] = 2;
  sood2[3][0] = 1;

  sood1[1][2] = 0;
  sood2[1][2] = 3;
  sood1[2][1] = 3;
  sood2[2][1] = 0;

  sood1[1][3] = 2;
  sood2[1][3] = 0;
  sood1[3][1] = 0;
  sood2[3][1] = 2;

  sood1[2][3] = 0;
  sood2[2][3] = 1;
  sood1[3][2] = 1;
  sood2[3][2] = 0;

  zero(&topo_stap); // topo_stap=0

  count=0;

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     count++;
     zero(&(M[count]));

     i0=sood1[i][j];
     j0=sood2[i][j];

//
//       i ^
//         |   (1)
//     (b) +----->-----+ (c)
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//     (a) +-----<-----+-->   j
//       r     (3)    (d)
//

     equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
     equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

     times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
     times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

     equal(&(M[count]), &stap);

     // clover insertion in (a)
     times(&aux, &stap, &(GC->clover_array[r][i0][j0])); // stap*clover
     plus_equal(&topo_stap, &aux);

     // clover insertion in (b)
     times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap);  // clover*stap
     plus_equal(&topo_stap, &aux);

     // clover insertion in (c)
     times(&aux, &link1, &(GC->clover_array[nnp(geo, nnp(geo, r, i), j)][i0][j0]));  // link1*clover
     times_equal_dag(&aux, &link2);       // *=link2^{dag}
     times_equal_dag(&aux, &link3);       // *=link3^{dag}
     plus_equal(&topo_stap, &aux);

     // clover insertion in (d)
     times(&aux, &link12, &(GC->clover_array[nnp(geo, r, j)][i0][j0]));  // link1*link2*quadri
     times_equal_dag(&aux, &link3);          // *=link3^{dag}

     plus_equal(&topo_stap, &aux);

//
//       i ^
//         |   (1)
//     (d) +----<------+ (a)
//         |           |
//         |
//     (2) V           |
//         |
//         |           | (b)
//     (c) +------>----+--->j
//        k     (3)    r
//

     count++;
     zero(&(M[count]));

     k=nnm(geo, r, j);

     equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
     equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

     times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
     times(&stap, &link12, &link3);        // stap=link12*link3

     equal(&(M[count]), &stap);

     // clover insertion in (a)
     times(&aux, &(GC->clover_array[nnp(geo, r, i)][i0][j0]), &stap); // clover*stap
     minus_equal(&topo_stap, &aux);

     // clover insertion in (b)
     times(&aux, &stap, &(GC->clover_array[r][i0][j0])); // stap*clover
     minus_equal(&topo_stap, &aux);

     // clover insertion in (c)
     times(&aux, &link12, &(GC->clover_array[k][i0][j0])); // link1^{dag}*link2^{dag}*clover
     times_equal(&aux, &link3);                            // *=link3
     minus_equal(&topo_stap, &aux);

     // clover insertion in (d)
     times_dag1(&aux, &link1, &(GC->clover_array[nnp(geo, k, i)][i0][j0]));  // link1^{dag}*clover
     times_equal_dag(&aux, &link2);            // *=link2^{dag}
     times_equal(&aux, &link3);                // *=link3

     minus_equal(&topo_stap, &aux);
     }

  times_equal_real(&topo_stap, coeff);
  equal(&(M[0]), &topo_stap);
  #endif


  #ifdef DEBUG
  GAUGE_GROUP helper;
  int m;

  calcstaples_with_topo(GC, geo, param, r, i, &helper);
  for(m=0; m<2*(STDIM-1)+1; m++)
     {
     minus_equal(&helper, &(M[m]));
     }
  if(norm(&helper)>MIN_VALUE)
    {
    fprintf(stderr, "Problems in calcstaples_with_topo_nosum (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif
  }


// perform an update with heatbath
void heatbath(Gauge_Conf *GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_GROUP stap;

  #ifndef THETA_MODE
    calcstaples_wilson(GC, geo, param, r, i, &stap);
  #else
    calcstaples_with_topo(GC, geo, param, r, i, &stap);
  #endif

  times_equal_real(&stap, param->d_beta);
  single_heatbath(&(GC->lattice[r][i]), &stap);
  }


// perform an update with overrelaxation
void overrelaxation(Gauge_Conf *GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif
  (void) param; // just to avoid wornings

  GAUGE_GROUP stap;

  #ifndef THETA_MODE
    calcstaples_wilson(GC, geo, param, r, i, &stap);
  #else
    calcstaples_with_topo(GC, geo, param, r, i, &stap);
  #endif

  single_overrelaxation(&(GC->lattice[r][i]), &stap);
  }


// perform an update with metropolis
// return 1 if the proposed update is accepted
int metropolis(Gauge_Conf *GC,
               Geometry const * const geo,
               GParam const * const param,
               long r,
               int i,
               int numhits)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_GROUP stap, new_link, tmp_matrix, rnd_matrix;
  double action_new, action_old;
  int acc, hits;

  #ifndef THETA_MODE
    calcstaples_wilson(GC, geo, param, r, i, &stap);
  #else
    calcstaples_with_topo(GC, geo, param, r, i, &stap);
  #endif

  acc=0;

  for(hits=0; hits<numhits; hits++)
     {
     // compute old action
     times(&tmp_matrix, &(GC->lattice[r][i]), &stap);
     action_old=param->d_beta*(1.0-retr(&tmp_matrix));

     // compute the new link
     one(&tmp_matrix);
     rand_matrix(&rnd_matrix);
     times_equal_real(&rnd_matrix, param->d_epsilon_metro);
     plus_equal(&rnd_matrix, &tmp_matrix);
     unitarize(&rnd_matrix);   // rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
     if(casuale()<0.5)
       {
       times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }
     else
       {
       times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }

     // new action
     times(&tmp_matrix, &new_link, &stap);
     action_new=param->d_beta*(1.0-retr(&tmp_matrix));

     if(casuale()< exp(action_old-action_new))
       {
       equal(&(GC->lattice[r][i]), &new_link);
       acc+=1;
       }
     else
       {
       acc+=0;
       }
     }

  return acc;
  }


// perform an update with metropolis with trace deformations
// return 1 if the proposed update is accepted
int metropolis_with_tracedef(Gauge_Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             long r,
                             int i,
                             int numhits)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_GROUP stap_w, stap_td, new_link, tmp_matrix, rnd_matrix, poly;
  double action_new, action_old;
  double rpart, ipart;
  int j, acc, hits;

  #ifndef THETA_MODE
    calcstaples_wilson(GC, geo, param, r, i, &stap_w);
  #else
    calcstaples_with_topo(GC, geo, param, r, i, &stap_w);
  #endif

  acc=0;

  for(hits=0; hits<numhits; hits++)
     {
     // compute old action
     times(&tmp_matrix, &(GC->lattice[r][i]), &stap_w);
     action_old=param->d_beta*(1.0-retr(&tmp_matrix));
     if(i==0) // just if we are updating a temporal link
       {
       // "staple" for trace deformation
       calcstaples_tracedef(GC, geo, param, r, i, &stap_td);

       // trace deformation contribution to action_old
       times(&poly, &(GC->lattice[r][i]), &stap_td);
       one(&tmp_matrix);
       for(j=0; j<(int)floor(NCOLOR/2.0); j++)
          {
          times_equal(&tmp_matrix, &poly);
          rpart=NCOLOR*retr(&tmp_matrix);
          ipart=NCOLOR*imtr(&tmp_matrix);
          action_old += param->d_h[j]*(rpart*rpart+ipart*ipart);
          }
       }

     // compute the update to be proposed
     one(&tmp_matrix);
     rand_matrix(&rnd_matrix);
     times_equal_real(&rnd_matrix, param->d_epsilon_metro);
     plus_equal(&rnd_matrix, &tmp_matrix);
     unitarize(&rnd_matrix);   // rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
     if(casuale()<0.5)
       {
       times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }
     else
       {
       times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }

     // compute the new action
     times(&tmp_matrix, &new_link, &stap_w);
     action_new=param->d_beta*(1.0-retr(&tmp_matrix));
     if(i==0) // just if we are updating a temporal link
       {
       // trace deformation contribution to action_new
       times(&poly, &new_link, &stap_td);
       one(&tmp_matrix);
       for(j=0; j<(int)floor(NCOLOR/2.0); j++)
          {
          times_equal(&tmp_matrix, &poly);
          rpart=NCOLOR*retr(&tmp_matrix);
          ipart=NCOLOR*imtr(&tmp_matrix);
          action_new += param->d_h[j]*(rpart*rpart+ipart*ipart);
          }
       }

     if(casuale()< exp(action_old-action_new))
       {
       equal(&(GC->lattice[r][i]), &new_link);
       acc+=1;
       }
     else
       {
       acc+=0;
       }
     }

  return acc;
  }


// perform an update with metropolis
// return 1 if the proposed update is accepted
int metropolis_fundadj(Gauge_Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int numhits)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_GROUP stap[2*(STDIM-1)+1], new_link, tmp_matrix, rnd_matrix;
  double action_new, action_old, re, im;
  int acc, count, hits;

  #ifndef THETA_MODE
    calcstaples_wilson_nosum(GC, geo, param, r, i, stap);
  #else
    calcstaples_with_topo_nosum(GC, geo, param, r, i, stap);
  #endif

  acc=0;

  for(hits=0; hits<numhits; hits++)
     {
     // compute old action
     times(&tmp_matrix, &(GC->lattice[r][i]), &(stap[0]));
     action_old=param->d_beta*(1.0-retr(&tmp_matrix));  // count=0 corresponds to theta term
     for(count=1; count<2*(STDIM-1)+1; count++)
        {
        times(&tmp_matrix, &(GC->lattice[r][i]), &(stap[count]));
        re=retr(&tmp_matrix);
        im=imtr(&tmp_matrix);

        action_old+=param->d_beta*(1.0-re);
        #if NCOLOR!=1
          action_old+=param->d_adjbeta*(1.0-(NCOLOR*NCOLOR*(re*re+im*im)-1.0)/(NCOLOR*NCOLOR-1));
        #else
          (void) im; // just to avoid warnings
        #endif
        }

     // compute the new link
     one(&tmp_matrix);
     rand_matrix(&rnd_matrix);
     times_equal_real(&rnd_matrix, param->d_epsilon_metro);
     plus_equal(&rnd_matrix, &tmp_matrix);
     unitarize(&rnd_matrix);   // rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
     if(casuale()<0.5)
       {
       times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }
     else
       {
       times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }

     // new action
     times(&tmp_matrix, &new_link, &(stap[0]));
     action_new=param->d_beta*(1.0-retr(&tmp_matrix));   // count=0 corresponds to theta term
     for(count=1; count<2*(STDIM-1)+1; count++)
        {
        times(&tmp_matrix, &new_link, &(stap[count]));
        re=retr(&tmp_matrix);
        im=imtr(&tmp_matrix);

        action_new+=param->d_beta*(1.0-re);
        #if NCOLOR!=1
          action_new+=param->d_adjbeta*(1.0-(NCOLOR*NCOLOR*(re*re+im*im)-1.0)/(NCOLOR*NCOLOR-1));
        #else
          (void) im; // just to avoid warnings
        #endif
        }

     if(casuale()< exp(action_old-action_new))
       {
       equal(&(GC->lattice[r][i]), &new_link);
       acc+=1;
       }
     else
       {
       acc+=0;
       }
     }

  return acc;
  }


// perform a complete update
void update(Gauge_Conf * GC,
            Geometry const * const geo,
            GParam const * const param)
   {
   for(int i=0; i<STDIM; i++)
      {
      if(param->d_size[i]==1)
        {
        fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   long r;
   int j, dir;

   // heatbath
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, dir);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif 
      for(r=0; r<(param->d_volume)/2; r++)
         {
         heatbath(GC, geo, param, r, dir);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif 
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         heatbath(GC, geo, param, r, dir);
         } 
      }

   // overrelax
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, dir);
      #endif

      for(j=0; j<param->d_overrelax; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=0; r<(param->d_volume)/2; r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif 
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }
         }
      }
   
   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif 
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         } 
      }

   GC->update_index++;
   }


// perform a complete update with trace deformation
void update_with_trace_def(Gauge_Conf * GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           double *acc)
   {
   int err, *a;
   long r, asum;
   int j, dir, t;

   const int maxhits=5;

   err=posix_memalign((void**)&a, (size_t)INT_ALIGN, (size_t) param->d_space_vol * sizeof(int));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(r=0; r<param->d_space_vol; r++)
      {
      a[r]=0;
      }

   // heatbath on spatial links
   for(dir=1; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, dir);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume)/2; r++)
         {
         heatbath(GC, geo, param, r, dir);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         heatbath(GC, geo, param, r, dir);
         }
      }

   // metropolis on temporal links
   for(t=0; t<param->d_size[0]; t++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, 0);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_space_vol)/2; r++)
         {
         long r4=sisp_and_t_to_si(geo, r, t);
         a[r]+=metropolis_with_tracedef(GC, geo, param, r4, 0, maxhits);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_space_vol)/2; r<(param->d_space_vol); r++)
         {
         long r4=sisp_and_t_to_si(geo, r, t);
         a[r]+=metropolis_with_tracedef(GC, geo, param, r4, 0, maxhits);
         }
      }

   asum=0;
   #ifdef OPENMP_MODE
   #pragma omp parallel for reduction(+:asum) private(r)
   #endif
   for(r=0; r<param->d_space_vol; r++)
      {
      asum+=(long)a[r];
      }

   *acc=((double)asum)*param->d_inv_vol/(double)maxhits;

   // overrelax spatial links
   for(dir=1; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, dir);
      #endif

      for(j=0; j<param->d_overrelax; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<(param->d_volume)/2; r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }
         }
      }

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }

   free(a);

   GC->update_index++;
   }


// perform a complete update with fundamental-adjoint action
void update_fundadj(Gauge_Conf * GC,
                   Geometry const * const geo,
                   GParam const * const param,
                   double *acc)
   {
   int err, *a;
   long r, asum;
   int dir;
   const int maxhits=5;

   err=posix_memalign((void**)&a, (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(int));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(r=0; r<param->d_volume; r++)
      {
      a[r]=0;
      }

   // metropolis
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, dir);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume)/2; r++)
         {
         a[r]+=metropolis_fundadj(GC, geo, param, r, dir, maxhits);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         a[r]+=metropolis_fundadj(GC, geo, param, r, dir, maxhits);
         }
      }

   asum=0;
   #ifdef OPENMP_MODE
   #pragma omp parallel for reduction(+:asum) private(r)
   #endif
   for(r=0; r<param->d_volume; r++)
      {
      asum+=(long)a[r];
      }

   *acc=((double)asum)*param->d_inv_vol/(double) STDIM/(double)maxhits;

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }

   free(a);

   GC->update_index++;
   }


// perform n cooling steps minimizing the action at theta=0
void cooling(Gauge_Conf *GC,
             Geometry const * const geo,
             GParam const * const param,
             int n)
  {
  long r;
  int i, k;

  for(k=0; k<n; k++)
     {
     // cooling
     for(i=0; i<STDIM; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<(param->d_volume)/2; r++)
           {
           GAUGE_GROUP staple;
           calcstaples_wilson(GC, geo, param, r, i, &staple);
           cool(&(GC->lattice[r][i]), &staple);
           }

        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=(param->d_volume)/2; r<(param->d_volume); r++)
           {
           GAUGE_GROUP staple;
           calcstaples_wilson(GC, geo, param, r, i, &staple);
           cool(&(GC->lattice[r][i]), &staple);
           }
        }
     }

  // final unitarization
  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r, i)
  #endif
  for(r=0; r<(param->d_volume); r++)
     {
     for(i=0; i<STDIM; i++)
        {
        unitarize(&(GC->lattice[r][i]));
        }
     }
  }


// perform a single step of the Runge Kutta integrator for the Wilson flow
// as described in Luscher arXiv:1006.4518 app. C
void gradflow_RKstep(Gauge_Conf *GC,
                     Gauge_Conf *helper1,
                     Gauge_Conf *helper2,
                     Geometry const * const geo,
                     GParam const *const param,
                     double dt)
  {
  long r;
  int dir;

  // initialize
  for(dir=0; dir<STDIM; dir++)
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        equal(&(helper1->lattice[r][dir]), &(GC->lattice[r][dir]));
        equal(&(helper2->lattice[r][dir]), &(GC->lattice[r][dir]));
        }
     }


  for(dir=0; dir<STDIM; dir++)
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        GAUGE_GROUP staple, aux, link;

        calcstaples_wilson(helper1, geo, param, r, dir, &staple);
        equal(&link, &(helper1->lattice[r][dir]));
        times(&aux, &link, &staple);                // aux=link*staple
        times_equal_real(&aux, -dt/4.0);
        equal(&(helper2->lattice[r][dir]), &aux);    // helper2=aux
        taexp(&aux);
        times(&(GC->lattice[r][dir]), &aux, &link); // GC=aux*link
        }
     }

  // now helper1=W_0, helper2=(1/4)Z_0 and GC=W_1

  for(dir=0; dir<STDIM; dir++)
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        GAUGE_GROUP staple, aux, link;

        calcstaples_wilson(GC, geo, param, r, dir, &staple);
        equal(&link, &(GC->lattice[r][dir]));
        times(&aux, &link, &staple);               // aux=link*staple
        times_equal_real(&aux, -dt*8.0/9.0);
        minus_equal_times_real(&aux, &(helper2->lattice[r][dir]), 17.0/9.0); // 1/4 was in Z_0
        equal(&(helper2->lattice[r][dir]), &aux);
        taexp(&aux);
        times(&(helper1->lattice[r][dir]), &aux, &link); // helper1=aux*link
        }
     }

  // now helper1=W_2, helper2=(8/9)Z_1-(17/36)Z_0 and GC=W_1

  for(dir=0; dir<STDIM; dir++)
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        GAUGE_GROUP staple, aux, link;

        calcstaples_wilson(helper1, geo, param, r, dir, &staple);
        equal(&link, &(helper1->lattice[r][dir]));
        times(&aux, &link, &staple);                   // aux=link*staple
        times_equal_real(&aux, -dt*3.0/4.0);
        minus_equal(&aux, &(helper2->lattice[r][dir])); // aux=(3/4)Z_2-(8/9)Z_1+(17/36)Z_0
        taexp(&aux);
        times(&(GC->lattice[r][dir]), &aux, &link);    // GC=aux*link
        }
     }

  // final unitarization
  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r)
  #endif
  for(r=0; r<(param->d_volume); r++)
     {
     int i;
     for(i=0; i<STDIM; i++)
        {
        unitarize(&(GC->lattice[r][i]));
        }
     }
  }



// n step of ape smearing with parameter alpha
// new=Proj[old + alpha *staple ]
void ape_smearing(Gauge_Conf *GC,
                  Geometry const * const geo,
                  GParam const *const param,
                  double alpha,
                  int n)
  {
  Gauge_Conf helper1;
  long r;
  int dir, count;

  init_gauge_conf_from_gauge_conf(&helper1, GC, param); //helper1=GC

  for(count=0; count<n; count++)
     {
     if(count%2==0) // smear(helper1)->GC
       {
       for(dir=0; dir<STDIM; dir++)
          {
          #ifdef OPENMP_MODE
          #pragma omp parallel for num_threads(NTHREADS) private(r)
          #endif
          for(r=0; r<param->d_volume; r++)
             {
             GAUGE_GROUP staple, link;

             calcstaples_wilson(&helper1, geo, param, r, dir, &staple);
             equal(&link, &(helper1.lattice[r][dir]));
             times_equal_real(&link, 1-alpha);
             times_equal_real(&staple, alpha/6.0);
             plus_equal_dag(&link, &staple);
             unitarize(&link);
             equal(&(GC->lattice[r][dir]), &link);
             }
          }
       }
     else // smear(GC)->helper1
       {
       for(dir=0; dir<STDIM; dir++)
          {
          #ifdef OPENMP_MODE
          #pragma omp parallel for num_threads(NTHREADS) private(r)
          #endif
          for(r=0; r<param->d_volume; r++)
             {
             GAUGE_GROUP staple, link;

             calcstaples_wilson(GC, geo, param, r, dir, &staple);
             equal(&link, &(GC->lattice[r][dir]));
             times_equal_real(&link, 1-alpha);
             times_equal_real(&staple, alpha/6.0);
             plus_equal_dag(&link, &staple);
             unitarize(&link);
             equal(&(helper1.lattice[r][dir]), &link);
             }
          }
       }
     }

  if(n>0 && n%2==0) // GC=helper1
    {
    for(dir=0; dir<STDIM; dir++)
       {
       #ifdef OPENMP_MODE
       #pragma omp parallel for num_threads(NTHREADS) private(r)
       #endif
       for(r=0; r<param->d_volume; r++)
          {
          equal(&(GC->lattice[r][dir]), &(helper1.lattice[r][dir]));
          }
       }
    }

  free_gauge_conf(&helper1, param);
  }


// perform an update with heatbath in the higgs theory
void heatbath_with_higgs(Gauge_Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r,
                         int i)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_GROUP stap1, stap2;

  if(fabs(param->d_beta)>MIN_VALUE)
    {
    calcstaples_wilson(GC, geo, param, r, i, &stap1);
    times_equal_real(&stap1, param->d_beta);
    }
  else
    {
    zero(&stap1);
    }

  if(fabs(param->d_higgs_beta)>MIN_VALUE)
    {
    vector_tensor_vector_vecs(&stap2, &(GC->higgs[r]), &(GC->higgs[nnp(geo, r, i)]));
    times_equal_real(&stap2, param->d_higgs_beta*NHIGGS*NCOLOR); // NCOLOR is needed to compensate the
                                                                 // 1/NCOLOR that is present in the gauge part
    }
  else
    {
    zero(&stap2);
    }

  plus_equal(&stap1, &stap2);

  single_heatbath(&(GC->lattice[r][i]), &stap1);
  }


// perform an update with overrelaxation in the higgs theory
void overrelaxation_with_higgs(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               long r,
                               int i)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_GROUP stap1, stap2;

  if(fabs(param->d_beta)>MIN_VALUE)
    {
    calcstaples_wilson(GC, geo, param, r, i, &stap1);
    times_equal_real(&stap1, param->d_beta);
    }
  else
    {
    zero(&stap1);
    }

  if(fabs(param->d_higgs_beta)>MIN_VALUE)
    {
    vector_tensor_vector_vecs(&stap2, &(GC->higgs[r]), &(GC->higgs[nnp(geo, r, i)]));
    times_equal_real(&stap2, param->d_higgs_beta*NHIGGS*NCOLOR); // NCOLOR is needed to compensate the
                                                                 // 1/NCOLOR that is present in the gauge part
    }
  else
    {
    zero(&stap2);
    }

  plus_equal(&stap1, &stap2);

  single_overrelaxation(&(GC->lattice[r][i]), &stap1);
  }


// compute the staple for the higgs field
void calcstaples_for_higgs(Gauge_Conf *GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           long r,
                           GAUGE_VECS *staple)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #else
  (void) param; // just to avoid warnings
  #endif

  int i;
  GAUGE_VECS aux;
  GAUGE_GROUP matrix;

  zero_vecs(staple);

  // forward
  for(i=0; i<STDIM; i++)
     {
     equal(&matrix, &(GC->lattice[r][i]) );
     matrix_times_vector_all_vecs(&aux, &matrix, &(GC->higgs[nnp(geo, r, i)]) );
     plus_equal_vecs(staple, &aux);
     }

  // backward
  for(i=0; i<STDIM; i++)
     {
     equal_dag(&matrix, &(GC->lattice[nnm(geo, r, i)][i]) );
     matrix_times_vector_all_vecs(&aux, &matrix, &(GC->higgs[nnm(geo, r, i)]) );
     plus_equal_vecs(staple, &aux);
     }
  }


// perform an update of the higgs field with overrelaxation
void overrelaxation_for_higgs(Gauge_Conf *GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              long r)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_VECS staple;

  calcstaples_for_higgs(GC, geo, param, r, &staple);

  single_overrelaxation_vecs(&(GC->higgs[r]), &staple);

  normalize_vecs(&(GC->higgs[r]));
  }


// perform an update of the higgs field with metropolis
// return the number of accepted moves (from 0 to NHIGGS)
int metropolis_for_higgs(Gauge_Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j, acc;
  double old_energy, new_energy;
  GAUGE_VECS staple, new_vector;
  GAUGE_GROUP identity_matrix, matrix, rnd_matrix;

  calcstaples_for_higgs(GC, geo, param, r, &staple);

  one(&identity_matrix);

  acc=0;
  for(i=0; i<NHIGGS; i++)
     {
     old_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&(GC->higgs[r]), &staple);

     j=(int)(NHIGGS*casuale()*(1.0 - MIN_VALUE));

     rand_matrix(&matrix);
     times_equal_real(&matrix, param->d_epsilon_metro);
     plus_equal(&matrix, &identity_matrix);
     unitarize(&matrix);   // matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
     if(casuale()<0.5)
       {
       equal(&rnd_matrix, &matrix);
       matrix_times_vector_single_vecs(&new_vector, &rnd_matrix, &(GC->higgs[r]), j);
       }
     else
       {
       equal_dag(&rnd_matrix, &matrix);
       matrix_times_vector_single_vecs(&new_vector, &rnd_matrix, &(GC->higgs[r]), j);
       }

     new_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&new_vector, &staple);

     if(casuale()< exp(old_energy-new_energy))
       {
       equal_vecs(&(GC->higgs[r]), &new_vector);
       acc+=1;
       }
     }
  normalize_vecs(&(GC->higgs[r]));

  #if GGROUP == 0 // only for SuN groups
  for(i=0; i<NHIGGS; i++)
     {
     old_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&(GC->higgs[r]), &staple);

     equal_vecs(&new_vector, &(GC->higgs[r]));
     times_equal_complex_single_vecs(&new_vector, cexp(param->d_epsilon_metro*PI*(2.0*casuale()-1)*I), i);
     new_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&new_vector, &staple);

     if(casuale()< exp(old_energy-new_energy))
       {
       equal_vecs(&(GC->higgs[r]), &new_vector);
       acc+=1;
       }
     }
  #else  // SoN groups
  for(i=0; i<NHIGGS; i++)
     {
     if(casuale()<0.5)
       {
       old_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&(GC->higgs[r]), &staple);
  
       equal_vecs(&new_vector, &(GC->higgs[r]));
       times_equal_real_single_vecs(&new_vector, -1.0, i);
       new_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&new_vector, &staple);
  
       if(casuale()< exp(old_energy-new_energy))
         {
         equal_vecs(&(GC->higgs[r]), &new_vector);
         acc+=1;
         }
       }
     }
  #endif

  if(NHIGGS>1)
    {
    int k;
    double angle;

    for(i=0; i<NHIGGS; i++)
       {
       old_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&(GC->higgs[r]), &staple);

       j=(int)(NHIGGS*casuale()*(1.0 - MIN_VALUE));
       k=(j+1 + (int)((NHIGGS-1)*casuale()*(1.0 - MIN_VALUE)) )% NHIGGS;

       angle=PI*param->d_epsilon_metro*(2.0*casuale()-1.0);

       rotate_two_components_vecs(&new_vector, &(GC->higgs[r]), j, k, angle);

       new_energy=-NHIGGS*param->d_higgs_beta*re_scal_prod_vecs(&new_vector, &staple);

       if(casuale()< exp(old_energy-new_energy))
         {
         equal_vecs(&(GC->higgs[r]), &new_vector);
         acc+=1;
         }
       }
    normalize_vecs(&(GC->higgs[r]));
    }


  return acc/3;
  }


// perform a complete update with higgs field
void update_with_higgs(Gauge_Conf * GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *acc)
   {
   #ifdef THETA_MODE
    fprintf(stderr, "THETA_MODE not yet implemented in the higgs case, check everything (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
   #endif

   int err, *a;
   long r, asum;
   int j, dir;

   err=posix_memalign((void**)&a, (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(int));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(r=0; r<param->d_volume; r++)
      {
      a[r]=0;
      }

   // heatbath on links
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume)/2; r++)
         {
         heatbath_with_higgs(GC, geo, param, r, dir);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         heatbath_with_higgs(GC, geo, param, r, dir);
         }
      }

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }

   // overrelax links and higgs
   for(j=0; j<param->d_overrelax; j++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<(param->d_volume)/2; r++)
            {
            overrelaxation_with_higgs(GC, geo, param, r, dir);
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            overrelaxation_with_higgs(GC, geo, param, r, dir);
            }

         // final unitarization
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
         #endif
         for(r=0; r<(param->d_volume); r++)
            {
            for(dir=0; dir<STDIM; dir++)
               {
               unitarize(&(GC->lattice[r][dir]));
               }
            }
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume)/2; r++)
         {
         overrelaxation_for_higgs(GC, geo, param, r);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         overrelaxation_for_higgs(GC, geo, param, r);
         }

      // normalization for higgs is included in the update functions
      }

   // metropolis on higgs
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume)/2; r++)
      {
      a[r]+=metropolis_for_higgs(GC, geo, param, r);
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=(param->d_volume)/2; r<(param->d_volume); r++)
      {
      a[r]+=metropolis_for_higgs(GC, geo, param, r);
      }

   // normalization for higgs is included in the update functions

   // acceptance computation
   asum=0;
   #ifdef OPENMP_MODE
   #pragma omp parallel for reduction(+:asum) private(r)
   #endif
   for(r=0; r<param->d_volume; r++)
      {
      asum+=(long)a[r];
      }

   *acc=((double)asum)*param->d_inv_vol/(double)NHIGGS;

   free(a);

   GC->update_index++;
   }






#endif
