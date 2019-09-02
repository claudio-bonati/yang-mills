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
                       GParam const * const param,
                       FILE *monofilep)
  {
  // Indici
  long r;
  int i, dir, err;
  

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
 
  // Inizializzo la matrice lambda L= diag((N-1)/2, (N-1)/2-1, ..., -(N-1)/2)
    for(i=0; i<NCOLOR; i++)
       {
       lambda[i] = (double) (NCOLOR -1)/2 - i;
       }

  // Adesso nel codice fortran vengo unitarizzati tutti i link del reticolo, qui lo sono già

  // metto a uno la matrice dell'operatore X(n) e G_mag
  zero(&X);
  one(&G_mag);

  //
  //      ORA COMINCIA IL GAUGE FIXING VERO E PROPRIO
  //
  comp_operator_X(GC, geo, param, r=0, lambda, &X);
/*
  while(non_diag_contr > precision)
       {
       for(r=0; r<(param->d_volume); r++) // Ciclo su tutti i siti del reticolo
          {
          // Prima viene calcolato l'operatore X(n) nel sito r
          comp_operator_X(GC, geo, param, r, lambda, &X);
     
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
         }*/
  free(lambda);
  }



// questa funzione calcola l'operatore X(n) che deve essere massimizzato
//  X(n) = sum_alldir { U_{\mu}(n)*lambda*Udag_mu(n) + Udag_{\mu}(n - \mu)*lambda*U_mu(n - \mu)}
void comp_operator_X (Gauge_Conf const * const GC,
		      Geometry const * const geo,
                      GParam const * const param,
                      long r,
                      double *lambda,
                      GAUGE_GROUP *X)
  {
  int dir, i;
  GAUGE_GROUP aux1, aux2, aux3;

  for(dir=0; dir<NCOLOR; dir++) 
     {
     //calcolo lambda*Udag_mu(n) e lambda*U_mu(n- mu)
     diag_matrix_times_dag(&aux1, lambda, &(GC->lattice[r][dir]));
     diag_matrix_times(&aux2, lambda, &(GC->lattice[nnm(geo, r, dir)][dir]));
     
     //Calcolo la matrice X(n)
     times(&aux3, &(GC->lattice[r][dir]), &aux1);
     plus_equal(X, &aux3);
     times_dag1(&aux3, &(GC->lattice[nnm(geo, r, dir)][dir]), &aux2);
     plus_equal(X, &aux3);
     }   
  }




#endif




