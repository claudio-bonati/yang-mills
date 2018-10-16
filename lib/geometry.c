#ifndef GEOMETRY_C
#define GEOMETRY_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>

#include"../include/geometry.h"
#include"../include/gparam.h"

// single index 4d = even/odd lexicographic index 4d
// single index 3d = even/odd lexicographic index 3d
void init_indexing_lexeo(void)
  {
  cart_to_si = &cart_to_lexeo;
  si_to_cart = &lexeo_to_cart;
  lex_to_si = &lex_to_lexeo;
  si_to_lex = &lexeo_to_lex;
  sisp_and_t_to_si_compute=&lexeosp_and_t_to_lexeo;
  si_to_sisp_and_t_compute=&lexeo_to_lexeosp_and_t;
  }


// initialize geometry
void init_geometry(Geometry *geo, GParam const * const param)
  {
  int i, value, valuep, valuem, err;
  long r, rm, rp;
  int cartcoord[STDIM];

  // allocate memory
  err=posix_memalign((void**)&(geo->d_nnp), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_nnm), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(geo->d_nnp[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     err=posix_memalign((void**)&(geo->d_nnm[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  err=posix_memalign((void**)&(geo->d_timeslice), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(int));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_spacecomp), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_tsp), (size_t)INT_ALIGN, (size_t) param->d_size[0] * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<param->d_size[0]; r++)
     {
     err=posix_memalign((void**)&(geo->d_tsp[r]), (size_t)INT_ALIGN, (size_t) param->d_space_vol * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // INITIALIZE
  for(r=0; r<param->d_volume; r++)
     {
     si_to_cart(cartcoord, r, param);

     for(i=0; i<STDIM; i++)
        {
        value=cartcoord[i];

        valuep=value+1;
        if(valuep >= param->d_size[i])
          {
          valuep-=param->d_size[i];
          }
        cartcoord[i]=valuep;
        rp=cart_to_si(cartcoord, param);
        geo->d_nnp[r][i]=rp;

        valuem=value-1;
        if(valuem<0)
          {
          valuem+=param->d_size[i];
          }
        cartcoord[i]=valuem;
        rm=cart_to_si(cartcoord, param);
        geo->d_nnm[r][i]=rm;

        cartcoord[i]=value;
        }
     } // end of loop on r

  for(r=0; r<param->d_volume; r++)
     {
     si_to_sisp_and_t_compute(&rp, &value, r, param);
     geo->d_spacecomp[r]=rp;
     geo->d_timeslice[r]=value;
     geo->d_tsp[value][rp]=r;
     }

  #ifdef DEBUG
    test_geometry(geo, param);
  #endif
  }  


// free memory
void free_geometry(Geometry *geo, GParam const * const param)
  {
  long r;

  for(r=0; r<param->d_volume; r++)
     {
     free(geo->d_nnp[r]);
     free(geo->d_nnm[r]);
     }
  free(geo->d_nnp);
  free(geo->d_nnm);

  free(geo->d_timeslice);
  free(geo->d_spacecomp);
  for(r=0; r<param->d_size[0]; r++)
     {
     free(geo->d_tsp[r]);
     }
  free(geo->d_tsp);
  }


long nnp(Geometry const * const geo, long r, int i);


long nnm(Geometry const * const geo, long r, int i);


long sisp_and_t_to_si(Geometry const * const geo, long sisp, int t);


void si_to_sisp_and_t(long *sisp, int *t, Geometry const * const geo, long si);


void test_geometry(Geometry const * const geo, GParam const * const param)
  {
  long si, ris_test, si_bis, sisp;
  int dir, cart[STDIM], cartsp[STDIM-1], t;

  // test of lex_to_cart <-> cart_to_lex
  for(si=0; si < param->d_volume; si++)
     {
     lex_to_cart(cart, si, param);
     ris_test=cart_to_lex(cart, param);

     if(si != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of lexeo_to_cart <-> cart_to_lexeo
  for(si=0; si < param->d_volume; si++)
     {
     lexeo_to_cart(cart, si, param);
     ris_test=cart_to_lexeo(cart, param);

     if(si != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of nnp <-> nnm
  for(si=0; si < param->d_volume; si++)
     {
     for(dir=0; dir<STDIM; dir++)
        {
        si_bis=nnp(geo, si, dir);
        ris_test=nnm(geo, si_bis, dir);

        if(si != ris_test)
          {
          fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }

  // test of lexsp_to_cartsp <-> cartsp_to_lexsp
  for(sisp=0; sisp < param->d_space_vol; sisp++)
     {
     lexsp_to_cartsp(cartsp, sisp, param);
     ris_test=cartsp_to_lexsp(cartsp, param);

     if(sisp != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of lexeosp_to_cartsp <-> cartsp_to_lexeosp
  for(sisp=0; sisp < param->d_space_vol; sisp++)
     {
     lexeosp_to_cartsp(cartsp, sisp, param);
     ris_test=cartsp_to_lexeosp(cartsp, param);

     if(sisp != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

   // test of lexeosp_and_t_to_lexeo <-> lexeo_to_lexeosp_and_t
   for(si=0; si<param->d_volume; si++)
      {
      lexeo_to_lexeosp_and_t(&sisp, &t, si, param);
      ris_test=lexeosp_and_t_to_lexeo(sisp, t, param);

      if(si != ris_test)
        {
        fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }
  }



//------------ these are not to be used outside geometry.c ----------------

// cartesian coordinates -> lexicographic index
long cart_to_lex(int const * const cartcoord, GParam const * const param)
  {
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=0; i<STDIM; i++)
     {
     ris+=cartcoord[i]*aux;
     aux*=param->d_size[i];
     }

  // ris = cartcoord[0]
  //      +cartcoord[1]*size[0]
  //      +cartcoord[2]*size[0]*size[1]
  //      +...
  //      +cartcoord[STDIM-1]*size[0]*size[1]*...*size[STDIM-2]

  return ris;

  /*
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=STDIM-1; i>=0; i--)
     {
     ris+=cartcoord[i]*aux;
     aux*=param->d_size[i];
     }

  // ris = cartcoord[STDIM-1] +
  //      +cartcoord[STDIM-2]*size[STDIM-1]+
  //      +cartcoord[STDIM-3]*size[STDIM-1]*size[STDIM-2]+
  //      +...
  //      +cartcoord[0]*size[STDIM-1]*size[STDIM-2]*...*size[1]

  return ris;
  */
  }


// lexicographic index -> cartesian coordinates
void lex_to_cart(int *cartcoord, long lex, GParam const * const param)
  {
  int i;
  long aux[STDIM];

  aux[0]=1;
  for(i=1; i<STDIM; i++)
     {
     aux[i]=aux[i-1]*param->d_size[i-1];
     }
  // aux[0]=1
  // aux[1]=size[0]
  // aux[2]=size[0]*size[1]
  // ...
  // aux[STDIM-1]=size[0]*size[1]*...*size[STDIM-2]

  for(i=STDIM-1; i>=0; i--)
     {
     cartcoord[i]=(int) (lex/aux[i]);
     lex-=aux[i]*cartcoord[i];
     }

  /*
  int i;
  long aux[STDIM];

  aux[STDIM-1]=1;
  for(i=STDIM-2; i>=0; i--)
     {
     aux[i]=aux[i+1]*param->d_size[i+1];
     }
  // aux[STDIM-1] = 1
  // aux[STDIM-2] = size[STDIM-1]
  // aux[STDIM-3] = size[STDIM-1]*size[STDIM-2]
  // ...
  // aux[0]       = size[STDIM-1]*size[STDIM-2]*...*size[1]

  for(i=0; i<STDIM; i++)
     {
     cartcoord[i]=(int) (lex/aux[i]);
     lex-=aux[i]*cartcoord[i];
     }
  */
  }


// cartesian coordinates -> lexicographic eo index
long cart_to_lexeo(int const * const cartcoord, GParam const * const param)
  {
  long lex;
  int i, eo;

  lex=cart_to_lex(cartcoord, param);

  eo=0;
  for(i=0; i<STDIM; i++)
     {
     eo+=cartcoord[i];
     }

  if(eo % 2==0)
    {
    return lex/2;
    }
  else
    {
    return (lex + param->d_volume)/2;
    }
  // even sites are written first
  }


// lexicographic eo index -> cartesian coordinates
void lexeo_to_cart(int *cartcoord, long lexeo, GParam const * const param)
  {
  long lex;
  int i, eo;

  if(param->d_volume % 2 == 0)
    {
    if(lexeo < param->d_volume/2)
      {
      lex=2*lexeo;
      }
    else
      {
      lex=2*(lexeo-param->d_volume/2);
      }
    lex_to_cart(cartcoord, lex, param);

    eo=0;
    for(i=0; i<STDIM; i++)
       {
       eo+=cartcoord[i];
       }
    eo = eo % 2;

    // this is to take care of the case d_volume is even but not
    // all the lattice extents are even
    if( (eo == 0 && lexeo >= param->d_volume/2) ||
        (eo == 1 && lexeo < param->d_volume/2) )
      {
      lex+=1;
      lex_to_cart(cartcoord, lex, param);
      }
    }
  else
    {
    if(lexeo <= param->d_volume/2)
      {
      lex=2*lexeo;
      }
    else
      {
      lex=2*(lexeo-param->d_volume/2)-1;
      }
    lex_to_cart(cartcoord, lex, param);
    }
  }


//  lexicographic index -> lexicographic eo index
long lex_to_lexeo(long lex, GParam const * const param)
  {
  int cartcoord[STDIM];

  lex_to_cart(cartcoord, lex, param);

  return cart_to_lexeo(cartcoord, param);
  }


//  lexicographic eo index -> lexicographic index
long lexeo_to_lex(long lexeo, GParam const * const param)
  {
  int cartcoord[STDIM];

  lexeo_to_cart(cartcoord, lexeo, param);

  return cart_to_lex(cartcoord, param);
  }


// spatial cartesian coordinates -> spatial lexicographic index
long cartsp_to_lexsp(int const * const ccsp, GParam const * const param)
  {
  // the index for the spatial cartesian coord. goes from 0 to STDIM-2 hence ccsp[STDIM-1]
  // cc   = t x1 x2 ... x_{STDIM-1}
  // ccsp =   x1 x2     x_{STDIM-1}
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=0; i<STDIM-1; i++)
     {
     ris+=ccsp[i]*aux;
     aux*=param->d_size[i+1];
     }

  // ris = ccsp[0]
  //      +ccsp[1]*size[1]
  //      +ccsp[2]*size[1]*size[2]
  //      +...
  //      +ccsp[STDIM-2]*size[1]*size[2]*...*size[STDIM-2]

  return ris;

  /*
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=STDIM-2; i>=0; i--)
     {
     ris+=ccsp[i]*aux;
     aux*=param->d_size[i+1];
     }

  // ris = ccsp[STDIM-2] +
  //      +ccsp[STDIM-3]*size[STDIM-1]+
  //      +ccsp[STDIM-4]*size[STDIM-1]*size[STDIM-2]+
  //      +ccsp
  //      +ccsp[0]*size[STDIM-1]*size[STDIM-2]*...*size[2]

  return ris;
  */
  }


// spatial lexicographic index -> spatial cartesian coordinates
void lexsp_to_cartsp(int *ccsp, long lexsp, GParam const * const param)
  {
  // the index for the spatial cartesian coord. goes from 0 to STDIM-2 hence ccsp[STDIM-1]
  // cc   = t x1 x2 ... x_{STDIM-1}
  // ccsp =   x1 x2     x_{STDIM-1}

  int i;
  long aux[STDIM-1];

  aux[0]=1;
  for(i=1; i<STDIM-1; i++)
     {
     aux[i]=aux[i-1]*param->d_size[i];
     }
  // aux[0]=1
  // aux[1]=size[1]
  // aux[2]=size[1]*size[2]
  // ...
  // aux[STDIM-2]=size[1]*size[2]*...*size[STDIM-2]

  for(i=STDIM-2; i>=0; i--)
     {
     ccsp[i]=(int) (lexsp/aux[i]);
     lexsp-=aux[i]*ccsp[i];
     }

  /*
  int i;
  long aux[STDIM-1];

  aux[STDIM-2]=1;
  for(i=STDIM-3; i>=0; i--)
     {
     aux[i]=aux[i+1]*param->d_size[i+2];
     }
  // aux[STDIM-2] = 1
  // aux[STDIM-3] = size[STDIM-1]
  // aux[STDIM-4] = size[STDIM-1]*size[STDIM-2]
  // ...
  // aux[0]       = size[STDIM-1]*size[STDIM-2]*...*size[2]

  for(i=0; i<STDIM-1; i++)
     {
     ccsp[i]=(int) (lexsp/aux[i]);
     lexsp-=aux[i]*ccsp[i];
     }
  */
  }


// spatial cartesian coordinates -> spatial lexicographic eo index
long cartsp_to_lexeosp(int const * const ccsp, GParam const * const param)
  {
  long lexsp;
  int i, eo;

  lexsp=cartsp_to_lexsp(ccsp, param);

  eo=0;
  for(i=0; i<STDIM-1; i++)
     {
     eo+=ccsp[i];
     }

  if(eo % 2==0)
    {
    return lexsp/2;
    }
  else
    {
    return (lexsp + param->d_space_vol)/2;
    }
  }


// spatial lexicographic eo index -> spatial cartesian coordinates
void lexeosp_to_cartsp(int *ccsp, long lexeosp, GParam const * const param)
  {
  long lexsp;
  int i, eo;

  if(param->d_space_vol % 2 == 0)
    {
    if(lexeosp < param->d_space_vol/2)
      {
      lexsp=2*lexeosp;
      }
    else
      {
      lexsp=2*(lexeosp - param->d_space_vol/2);
      }
    lexsp_to_cartsp(ccsp, lexsp, param);

    eo=0;
    for(i=0; i<STDIM-1; i++)
       {
       eo+=ccsp[i];
       }
    eo = eo % 2;

    if( (eo == 0 && lexeosp >= param->d_space_vol/2) ||
        (eo == 1 && lexeosp < param->d_space_vol/2) )
      {
      lexsp+=1;
      lexsp_to_cartsp(ccsp, lexsp, param);
      }
    }
  else
    {
    if(lexeosp <= param->d_space_vol/2)
      {
      lexsp=2*lexeosp;
      }
    else
      {
      lexsp=2*(lexeosp - param->d_space_vol/2)-1;
      }
    lexsp_to_cartsp(ccsp, lexsp, param);
    }
  }


// spatial lexicographic index -> spatial lexicographic eo index
long lexsp_to_lexeosp(long lexsp, GParam const * const param)
  {
  int ccsp[STDIM];

  lexsp_to_cartsp(ccsp, lexsp, param);

  return cartsp_to_lexeosp(ccsp, param);
  }


//  spatial lexicographic eo index -> spatial lexicographic index
long lexeosp_to_lexsp(long lexeosp, GParam const * const param)
  {
  int ccsp[STDIM];

  lexeosp_to_cartsp(ccsp, lexeosp, param);

  return cartsp_to_lexsp(ccsp, param);
  }


// lexicographic eo spatial and time -> lexicographic eo index
long lexeosp_and_t_to_lexeo(long lexeosp, int t, GParam const * const param)
  {
  int cc[STDIM];

  lexeosp_to_cartsp(cc+1, lexeosp, param);
  cc[0]=t;

  return cart_to_lexeo(cc, param);
  }


// lexicographic eo index -> lexicographic eo spatial and time
void lexeo_to_lexeosp_and_t(long *lexeosp, int *t, long lexeo, GParam const * const param)
  {
  int i, cc[STDIM], ccsp[STDIM-1];

  lexeo_to_cart(cc, lexeo, param);

  *t=cc[0];

  for(i=0; i<STDIM-1; i++)
     {
     ccsp[i]=cc[i+1];
     }

  *lexeosp=cartsp_to_lexeosp(ccsp, param);
  }


#endif
