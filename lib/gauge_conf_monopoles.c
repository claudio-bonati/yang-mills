#ifndef GAUGE_CONF_MONO_C
#define GAUGE_CONF_MONO_C

#include"../include/macro.h"

#include<complex.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"


void max_abelian_gauge(Gauge_Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param)
  {
  // Indici
  long r;
  int dir;
  

  //params
  double OverrelaxParameter = 1.85;     // The parameter for the Overrelaxation is the una used in the SU(2) case
  double precision = 0.0000000001;      // This is used to check that the diagonal elements of X(n) are near zero
  
  // variabili dell'algoritmo
  double funzionale_old;                // unico dubbio, mai veramente usato
  double *lambda;   //  L = diag((N-1)/2, (N-1)/2 -1, ..., -(N-1)/2)
  double non_diag_contr = 0; 
  GAUGE_GROUP X, G_mag;
  
 

  //alloco la matrice lambda che è un array lungo NCOLOR
  err=posix_memalign((void**)&lambda, (size_t)DOUBLE_ALIGN, (size_t) NCOLOR * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
 
  // Inizializzo la matrice lambda 
  compute_lambda_matrix(lambda);

  // Adesso nel codice fortran vengo unitarizzati tutti i link del reticolo, qui lo sono già


  //
  //      ORA COMINCIA IL GAUGE FIXING VERO E PROPRIO
  //

  while(non_diag_contr > precision)
       {
       for(r=0; r<(param->d_volume); r++) // Ciclo su tutti i siti del reticolo
          {
          // Prima viene calcolato l'operatore X(n) nel sito r
          comp_operator_X(Gc, geo, param, r, lambda, &X);
     
          // Viene ora calcolata la matrice di gauge G_mag massimizzando X
          maximize_operator_X(GC, geo, param, r, &X, &G_mag);
     
          // Ora si applica la trasformazione di gauge G_mag nel punto r
          for(dir=0; dir<STDIM; dir++)
             {
             apply_mag_single_site(GC, geo, param, &G_mag);
             }
          }
     
       // Viene calcolata la media del contributo degli elementi non diagonali di X(n)
       comp_non_diagonal_contribution(GC, &non_diag_contr)
      }
      
      // Faccio una unitarizzazione della conf una volta applicata la gauge
      for(r=0; r<(param->d_volume); r++)
         {
         for(dir=0; dir<STDIM; dir++)
            {
            unitarize(&(GC->lattice[r][dir]))
            }
         }
  }







