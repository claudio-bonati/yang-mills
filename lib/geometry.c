#ifndef GEOMETRY_C
#define GEOMETRY_C

#include<malloc.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/macro.h"
#include"../include/mymalloc.h"

//#define DEBUG

void init_indexing_lexeo(void)
  {
  cart_to_si = &cart_to_lexeo;
  si_to_cart = &lexeo_to_cart;
  lex_to_si = &lex_to_lexeo;
  si_to_lex = &lexeo_to_lex;
  sisp_and_t_to_si=&lexsp_and_t_to_lexeo;
  }


// initialize geometry
void init_geometry(Geometry *geo, GParam const * const param)
  {
  int i, value, valuep, valuem;
  long r, rm, rp;
  int cartcoord[STDIM];

  // allocate memory
  geo->d_nnp = (long **) mymalloc(INT_ALIGN, (unsigned long) param->d_volume * sizeof(long *));
  if(geo->d_nnp==NULL)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  geo->d_nnm = (long **) mymalloc(INT_ALIGN, (unsigned long) param->d_volume * sizeof(long *));
  if(geo->d_nnm==NULL)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     geo->d_nnp[r] = (long *) mymalloc(INT_ALIGN, STDIM * sizeof(long));
     if(geo->d_nnp[r]==NULL)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     geo->d_nnm[r] = (long *) mymalloc(INT_ALIGN, STDIM * sizeof(long));
     if(geo->d_nnm[r]==NULL)
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
  }


long nnp(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnp[r][i];
  }


long nnm(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnm[r][i];
  }


void test_geometry(Geometry const * const geo, GParam const * const param)
  {
  long si, ris_test, si_bis;
  int dir, cart[STDIM];

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
  // ris=cartcoord[0] + param->d_size[0]*cartcoord[1] + param->d_size[0]*param->d_size[1]*cartcoord[2] + ...

  return ris;
  }


// lexicographic index -> cartesian coordinates
void lex_to_cart(int *cartcoord, long lex, GParam const * const param)
  {
  int i;
  long aux[STDIM];

  aux[0]=1;
  for(i=1; i<STDIM; i++)
     {
     aux[i]=aux[i-1]*(param->d_size[i-1]);
     }

  for(i=STDIM-1; i>=0; i--)
     {
     cartcoord[i]=(int) (lex/aux[i]);
     lex-=aux[i]*cartcoord[i];
     }
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


// lexicographic spatial and time -> lexicographic index
long lexsp_and_t_to_lexeo(long lexsp, int t, GParam const * const param)
  {
  return t+param->d_size[0]*lexsp;
  }


#endif
