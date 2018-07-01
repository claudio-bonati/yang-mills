#ifndef SUN_C
#define SUN_C

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/sun.h"
#include"../include/sun_upd.h"
#include"../include/tens_prod.h"

// A=1
void one_SuN(SuN * restrict A)
  {
  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0+0.0*I;
     }

  for(i=0; i<NCOLOR; i++)
     {
     A->comp[m(i,i)]=1.0+0.0*I;
     }
  }


// A=0
void zero_SuN(SuN * restrict A)
  {
  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0+0.0*I;
     }
  }


// A=B
void equal_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A=B^{dag}
void equal_dag_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=conj(B->comp[m(j,i)]);
        }
     }
  }


// A+=B
void plus_equal_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A+=B^{dag}
void plus_equal_dag_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]+=conj(B->comp[m(j,i)]);
        }
     }
  }


// A-=B
void minus_equal_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A-=(r*B)
void minus_equal_times_real_SuN(SuN * restrict A, SuN const * restrict B, double r)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=(r*B->comp[i]);
     }
  }


// A-=B^{dag}
void minus_equal_dag_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]-=conj(B->comp[m(j,i)]);
        }
     }
  }


// A=b*B+c*C
void lin_comb_SuN(SuN * restrict A,
                  double b, SuN const * restrict B,
                  double c, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=b*(B->comp[i])+c*(C->comp[i]);
     }
  }


// A=b*B^{dag}+c*C
void lin_comb_dag1_SuN(SuN * restrict A,
                       double b, SuN const * restrict B,
                       double c, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*conj(B->comp[m(j,i)])+c*(C->comp[m(i,j)]);
        }
     }
  }


// A=b*B+c*C^{dag}
void lin_comb_dag2_SuN(SuN * restrict A,
                       double b, SuN const * restrict B,
                       double c, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*(B->comp[m(i,j)])+c*conj(C->comp[m(j,i)]);
        }
     }
  }


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_SuN(SuN * restrict A,
                        double b, SuN const * restrict B,
                        double c, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*conj(B->comp[m(j,i)])+c*conj(C->comp[m(j,i)]);
        }
     }
  }


// A*=r
void times_equal_real_SuN(SuN * restrict A, double r)
  {
  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=r
void times_equal_complex_SuN(SuN * restrict A, double complex r)
  {
  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=B
void times_equal_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  int i, j, k;
  double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }
 
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*(B->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A*=B^{dag}
void times_equal_dag_SuN(SuN * restrict A, SuN const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  int i, j, k;
  double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*conj(B->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B*C
void times_SuN(SuN * restrict A, SuN const * restrict B, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B^{dag}*C
void times_dag1_SuN(SuN * restrict A, SuN const * restrict B, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j, k;
  double complex sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[m(k,i)])*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     } 
  }


// A=B*C^{dag}
void times_dag2_SuN(SuN * restrict A, SuN const * restrict B, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*conj(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     } 
  }


// A=B^{dag}*C^{dag}
void times_dag12_SuN(SuN * restrict A, SuN const * restrict B, SuN const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j, k;
  double complex sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[m(k,i)])*conj(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     } 
  }


// SU(N) random matrix
// generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices
void rand_matrix_SuN(SuN * restrict A)
  {
  int i, j, k;
  double p0, p1, p2, p3, p;
  double complex aux00, aux01, aux10, aux11, temp0, temp1;

  one_SuN(A);

  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        // SU(2) random components
        p=2.0;
        while(p>1.0)
             { 
             p0=1.0-2.0*casuale();
             p1=1.0-2.0*casuale();
             p2=1.0-2.0*casuale();
             p3=1.0-2.0*casuale();
             p=sqrt(p0*p0+p1*p1+p2*p2+p3*p3);
             }

        p0/=p;
        p1/=p;
        p2/=p;
        p3/=p;

        aux00= p0 + p3*I;
        aux01= p2 + p1*I;
        aux10=-p2 + p1*I;
        aux11= p0 - p3*I;

        for(k=0; k<NCOLOR; k++)
           {
           temp0=A->comp[m(k,i)]*aux00 + A->comp[m(k,j)]*aux10;
           temp1=A->comp[m(k,i)]*aux01 + A->comp[m(k,j)]*aux11;
           A->comp[m(k,i)]=temp0;
           A->comp[m(k,j)]=temp1;
           }
        }
     }
  }


// generate a matrix in the algebra of SuN with gaussian
// random components in the base T_i such that Tr(T_iT_j)=delta_{ij}
void rand_algebra_gauss_matrix_SuN(SuN * restrict A)
  {
  #if NCOLOR==1
    (void) A; // just to avoid warnings
  #else
  int i, j;
  double d1, d2, dd[NCOLOR-1];
  const double nfactor=sqrt(2.0/(double)(NCOLOR*NCOLOR-NCOLOR));

  zero_SuN(A);

  // out of diagonal elements
  for(i=0; i<NCOLOR; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        gauss2(&d1, &d2);
        A->comp[m(i,j)]=d1-d2*I;
        A->comp[m(j,i)]=d1+d2*I;
        }
     }

  // random numbes to be used in the diagonal
  for(i=0; i<NCOLOR-1; i++)
     {
     dd[i]=gauss1();
     }

  // diagonal
  if(NCOLOR==2)
    {
    A->comp[m(0,0)]=dd[0];
    A->comp[m(1,1)]=-dd[0];
    }
  else
    {
    for(i=0; i<NCOLOR-2; i++)
       {
       A->comp[m(i,i)]+=dd[i];
       A->comp[m(i+1,i+1)]-=dd[i];
       }
    for(i=0; i<NCOLOR-1; i++)
       {
       A->comp[m(i,i)]+=nfactor*dd[NCOLOR-2];
       }
    A->comp[m(NCOLOR-1,NCOLOR-1)]=nfactor*(1.0-(double)NCOLOR)*dd[NCOLOR-2];
    }

  times_equal_real_SuN(A, 1./sqrt(2.0));
  #endif
  }


// l2 norm of the matrix
double norm_SuN(SuN const * restrict A)
  {
  int i;
  double aux, ris;

  ris=0.0;
  for(i=0; i<NCOLOR*NCOLOR; i++)
     { 
     aux=cabs(A->comp[i]);
     ris+=aux*aux;
     }
  return sqrt(ris);
  }


// real part of the trace /N
double retr_SuN(SuN const * restrict A)
  {
  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[m(i,i)];
     }
  ris=creal(tr)/(double)NCOLOR;
  return ris;
  }


// imaginary part of the trace /N
double imtr_SuN(SuN const * restrict A)
  {
  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[m(i,i)];
     }
  ris=cimag(tr)/(double)NCOLOR;
  return ris;
  }


// LU decomposition with partial pivoting
//   from Numerical Recipes in C, pag 46
void LU_SuN(SuN const * restrict A, SuN * restrict ris, int * restrict sign)
  {
  int i, imax, j, k;
  double  big, temp;
  double complex sum, dum;
  double vv[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));

  imax=0;
  equal_SuN(ris, A);

  (*sign)=1;
  for(i=0; i<NCOLOR; i++)
     {
     big=0.0;
     for(j=0; j<NCOLOR; j++)
        {
        temp = cabs(ris->comp[m(i,j)]);
        if( temp>big ) big=temp;
        }
     vv[i]=1.0/big;
     }

  for(j=0; j<NCOLOR; j++)
     {
     for(i=0; i<j; i++)
        {
        sum=ris->comp[m(i,j)];
        for(k=0; k<i; k++)
           { 
           sum-=(ris->comp[m(i,k)])*(ris->comp[m(k,j)]);
           }
        ris->comp[m(i,j)]=sum;
        }
     
     big=0.0;
     for(i=j; i<NCOLOR; i++)
        {
        sum=ris->comp[m(i,j)];
        for(k=0; k<j; k++)
           {
           sum-=(ris->comp[m(i,k)])*(ris->comp[m(k,j)]);
           }
        ris->comp[m(i,j)]=sum;
       
        temp = vv[i]*cabs(sum);
        if(temp >= big)
          {
          big=temp;
          imax=i;
          }
        }

     if(j!= imax)
       {
       for(k=0; k<NCOLOR; k++)
          {
          dum=ris->comp[m(imax,k)];
          ris->comp[m(imax,k)]=ris->comp[m(j,k)];
          ris->comp[m(j,k)]=dum;
          }
       (*sign)*=(-1);
       vv[imax]=vv[j];
       }
 
     if(j!= NCOLOR-1)
       {
       dum=(1.0+0.0*I)/(ris->comp[m(j,j)]);
       for(i=j+1; i<NCOLOR; i++)
          { 
          (ris->comp[m(i,j)])*=dum;
          }
       }
     }
  }


// determinant
complex double det_SuN(SuN const * restrict A)
  {
  #if NCOLOR==3
    complex double ris=0.0+0.0*I;
 
    ris+=(A->comp[m(0,0)])*(A->comp[m(1,1)])*(A->comp[m(2,2)]);
    ris+=(A->comp[m(1,0)])*(A->comp[m(2,1)])*(A->comp[m(0,2)]);
    ris+=(A->comp[m(2,0)])*(A->comp[m(0,1)])*(A->comp[m(1,2)]);
    ris-=(A->comp[m(2,0)])*(A->comp[m(1,1)])*(A->comp[m(0,2)]);
    ris-=(A->comp[m(1,0)])*(A->comp[m(0,1)])*(A->comp[m(2,2)]);
    ris-=(A->comp[m(0,0)])*(A->comp[m(2,1)])*(A->comp[m(1,2)]);

    return ris;
  #else
    int i;
    double complex ris;
    SuN lu; 

    LU_SuN(A, &lu, &i); 

    if(i>0)
      {
      ris=1.0+0.0*I;
      }
    else
     {
     ris=-1.0+0.0*I;
     }

    for(i=0; i<NCOLOR; i++)
       {
       ris*=(lu.comp[m(i,i)]);
       }

    return ris;
  #endif
  }


// gives 0 if the matrix is in SU(N) and 1 otherwise
int scheck_SuN(SuN const * restrict A)
  {
  int i, j, k, ris;
  double complex aux;

  ris=0;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           aux+=(A->comp[m(i,k)])*conj(A->comp[m(j,k)]);
           }
        if(i==j) aux-=(1.0+0.0*I);
        if(cabs(aux)>MIN_VALUE) ris=1;
        }
     }

  if(ris==0)
    {
    if(cabs(det_SuN(A)-1)>MIN_VALUE)
      {
      ris=1;
      }
    }

  return ris;
  }


// sunitarize
void unitarize_SuN(SuN * restrict A)
  {
  const double beta_aux=1.0e+20;
  double check;
  SuN force, guess, guess_old, helper, helper2;

  if(scheck_SuN(A)!=0)
    {
    one_SuN(&guess);
    check=1.0;

    equal_dag_SuN(&force, A);

    while(check>MIN_VALUE)
         {
         equal_SuN(&guess_old, &guess);

         // heatbath
         single_heatbath_aux_SuN(&guess, &force, beta_aux); // minimize Tr(force*guess)

         // compute the check
         equal_SuN(&helper, &guess);
         minus_equal_SuN(&helper, &guess_old);
         times_SuN(&helper2, &helper, &helper);
         check=sqrt(fabs(retr_SuN(&helper2))/(double) NCOLOR);

         //printf("aux: %g\n", check);
         }
    equal_SuN(A, &guess);
    }
  }


// eponential of the traceless antihermitian part
void taexp_SuN(SuN * restrict A)
  {
  SuN aux, aux1, uno;
  complex double trace;
  int i;

  equal_SuN(&aux, A);
  equal_dag_SuN(&aux1, A);
  minus_equal_SuN(&aux, &aux1);
  times_equal_real_SuN(&aux, 0.5); // now aux=(A-A^{dag})/2

  trace=aux.comp[m(0,0)];
  for(i=1; i<NCOLOR; i++)
     {
     trace+=aux.comp[m(i,i)];
     }
  trace/=(double)NCOLOR;
  for(i=0; i<NCOLOR; i++)
     {
     aux.comp[m(i,i)]-=trace;
     }

  one_SuN(&uno);

  // now aux is the traceless antihermitian part of the initial matrix
  // and we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/5.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/4.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/3.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/2.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  plus_equal_SuN(A, &uno);

  unitarize_SuN(A);
  }


// takes the traceless antihermitian part
void ta_SuN(SuN * restrict A)
  {
  SuN aux, aux1;
  double complex trace;
  int i;

  equal_SuN(&aux, A);
  equal_dag_SuN(&aux1, A);
  minus_equal_SuN(&aux, &aux1);
  times_equal_real_SuN(&aux, 0.5); // now aux=(A-A^{dag})/2

  trace=aux.comp[m(0,0)];
  for(i=1; i<NCOLOR; i++)
     {
     trace+=aux.comp[m(i,i)];
     }
  trace/=(double)NCOLOR;

  for(i=0; i<NCOLOR; i++)
     {
     aux.comp[m(i,i)]-=trace;
     }

  equal_SuN(A, &aux);
  }


// return 0 if matrix is traceless antihermitian, 1 otherwise
int ta_check_SuN(const SuN * restrict A)
  {
  double complex aux;
  int i, j, ris;

  ris=0;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux+=(A->comp[m(i,j)]+conj(A->comp[m(j,i)]));
        }
     }
  if(cabs(aux)>MIN_VALUE) ris=1;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     aux+=A->comp[m(i,i)];
     }
  if(cabs(aux)>MIN_VALUE) ris=1;

  return ris;
  }


// exponential of a TA matrix
void exp_of_ta_SuN(SuN * restrict A)
  {
  // we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  #ifdef DEBUG
  if(ta_check_SuN(A)!=0)
    {
    fprintf(stderr, "Trying to exp. a non TA matrix! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  SuN aux, uno;

  one_SuN(&uno);
  equal_SuN(&aux, A); // in aux the initial matrix is stored

  times_equal_real_SuN(A, 1.0/5.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/4.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/3.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/2.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  plus_equal_SuN(A, &uno);

  unitarize_SuN(A);
  }


void print_on_screen_SuN(SuN const * const A)
  {
  int i, j;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        printf("%.16f %.16f ", creal(A->comp[m(i,j)]), cimag(A->comp[m(i,j)]));
        }
     }
  printf("\n");
  }


void print_on_file_SuN(FILE *fp, SuN const * const A)
  {
  int i, j, err;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fprintf(fp, "%.16f %.16f ", creal(A->comp[m(i,j)]), cimag(A->comp[m(i,j)]));
        if(err<0)
          {
          fprintf(stderr, "Problem in writing on file a SuN matrix (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }
  fprintf(fp, "\n");
  }


// print on binary file without changing endiannes
void print_on_binary_file_noswap_SuN(FILE *fp, SuN const * const A)
  {
  int i, j;
  size_t err;
  double re, im;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        re=creal(A->comp[m(i,j)]);
        im=cimag(A->comp[m(i,j)]);

        err=fwrite(&re, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SuN matrix (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        err=fwrite(&im, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SuN matrix (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }
  }


// print on binary file changing endiannes
void print_on_binary_file_swap_SuN(FILE *fp, SuN const * const A)
  {
  int i, j;
  size_t err;
  double re, im;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        re=creal(A->comp[m(i,j)]);
        im=cimag(A->comp[m(i,j)]);

        SwapBytesDouble(&re);
        SwapBytesDouble(&im);

        err=fwrite(&re, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SuN matrix (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        err=fwrite(&im, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SuN matrix (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }
  }


// print on binary file in bigendian
void print_on_binary_file_bigen_SuN(FILE *fp, SuN const * const A)
  {
  if(endian()==0) // little endian machine
    {
    print_on_binary_file_swap_SuN(fp, A);
    }
  else
    {
     print_on_binary_file_noswap_SuN(fp, A);
    }
  }


void read_from_file_SuN(FILE *fp, SuN *A)
  {
  int i, j, err;
  double re, im;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fscanf(fp, "%lg %lg", &re, &im);
        if(err!=2)
          {
          fprintf(stderr, "Problems reading SuN matrix from file (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        A->comp[m(i,j)]=re+im*I; 
        }
     }
  }


// read from binary file without changing endiannes
void read_from_binary_file_noswap_SuN(FILE *fp, SuN *A)
  {
  size_t err;
  int i, j;
  double re, im;
  double aux[2];
 
  err=0;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err+=fread(&re, sizeof(double), 1, fp);
        err+=fread(&im, sizeof(double), 1, fp);
        aux[0]=re;
        aux[1]=im;

        memcpy((void *)&(A->comp[m(i,j)]), (void*)aux, sizeof(aux));
        //equivalent to A->comp[m(i,j)]=re+im*I;
        }
     }

  if(err!=2*NCOLOR*NCOLOR)
    {
    fprintf(stderr, "Problems reading SuN matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }
 

// read from binary file changing endianness
void read_from_binary_file_swap_SuN(FILE *fp, SuN *A)
  {
  int i, j;
  size_t err;
  double re, im;
  double aux[2];
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=0;
        err+=fread(&re, sizeof(double), 1, fp);
        err+=fread(&im, sizeof(double), 1, fp);
        if(err!=2)
         {
         fprintf(stderr, "Problems reading SuN matrix from file (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }

        SwapBytesDouble(&re);
        SwapBytesDouble(&im);
        aux[0]=re;
        aux[1]=im;

        memcpy((void *)&(A->comp[m(i,j)]), (void*)aux, sizeof(aux));
        // equivalent to A->comp[m(i,j)]=re+im*I;
        }
     }
  }


// read from binary file written in bigendian
void read_from_binary_file_bigen_SuN(FILE *fp, SuN *A)
  {
  if(endian()==0) // little endian machine
    {
    read_from_binary_file_swap_SuN(fp, A);
    }
  else
    {
    read_from_binary_file_noswap_SuN(fp, A);
    }
  }


void TensProd_init_SuN(TensProd * restrict TP, SuN const * restrict A1, SuN const * restrict A2)
  {
  int i, j, k, l;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              TP->comp[i][j][k][l]=conj(A1->comp[m(i,j)])*A2->comp[m(k,l)];
              }
           }
        }
     }
  }


#endif
