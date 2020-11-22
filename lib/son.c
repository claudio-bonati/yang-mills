#ifndef SON_C
#define SON_C

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/flavour_matrix.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/son.h"
#include"../include/son_upd.h"
#include"../include/tens_prod.h"
#include"../include/tens_prod_adj.h"


// ***************** for SoN



// A=1
void one_SoN(SoN *A);


// A=0
void zero_SoN(SoN *A);


// A=B
void equal_SoN(SoN *A, SoN const * const B);


// A=B^{dag}
void equal_dag_SoN(SoN *A, SoN const * const B);


// A+=B
void plus_equal_SoN(SoN *A, SoN const * const B);


// A+=B^{dag}
void plus_equal_dag_SoN(SoN *A, SoN const * const B);


// A-=B
void minus_equal_SoN(SoN *A, SoN const * const B);


// A-=(r*B)
void minus_equal_times_real_SoN(SoN *A, SoN const * const B, double r);


// A-=B^{dag}
void minus_equal_dag_SoN(SoN *A, SoN const * const B);


// A=b*B+c*C
void lin_comb_SoN(SoN *A,
                  double b, SoN const * const B,
                  double c, SoN const * const C);


// A=b*B^{dag}+c*C
void lin_comb_dag1_SoN(SoN *A,
                       double b, SoN const * const B,
                       double c, SoN const * const C);


// A=b*B+c*C^{dag}
void lin_comb_dag2_SoN(SoN *A,
                       double b, SoN const * const B,
                       double c, SoN const * const C);


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_SoN(SoN *A,
                        double b, SoN const * const B,
                        double c, SoN const * const C);


// A*=r
void times_equal_real_SoN(SoN *A, double r);


// A*=r
void times_equal_complex_SoN(SoN *A, double complex r);


// A*=B
void times_equal_SoN(SoN *A, SoN const * const B);


// A*=B^{dag}
void times_equal_dag_SoN(SoN *A, SoN const * const B);


// A=B*C
void times_SoN(SoN *A, SoN const * const B, SoN const * const C);


// A=B^{dag}*C
void times_dag1_SoN(SoN *A, SoN const * const B, SoN const * const C);


// A=B*C^{dag}
void times_dag2_SoN(SoN *A, SoN const * const B, SoN const * const C);


// A=B^{dag}*C^{dag}
void times_dag12_SoN(SoN *A, SoN const * const B, SoN const * const C);


// A = lambda*B with lambda diagonal marix
void diag_matrix_times_SoN(SoN * restrict A, double const lambda[NCOLOR], SoN const * const restrict B);


// A=lambda*B^{dag} with lambda diagonal matrix
void diag_matrix_times_dag_SoN(SoN * restrict A, double const lambda[NCOLOR], SoN const * const restrict B);


// SU(N) random matrix
// generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices
void rand_matrix_SoN(SoN *A)
  {
  int i, j, k;
  double p;
  double aux00, aux01, aux10, aux11, temp0, temp1;

  one_SoN(A);

  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        // SO(2) random components
        p=2*PI*casuale();

        aux00= cos(p);
        aux01= sin(p);
        aux10=-sin(p);
        aux11= cos(p);

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


// l2 norm of the matrix
double norm_SoN(SoN const * const A);


// real part of the trace /N
double retr_SoN(SoN const * const A);


// imaginary part of the trace /N
double imtr_SoN(SoN const * const A);


// LU decomposition with partial pivoting
//   from Numerical Recipes in C, pag 46
void LU_SoN(SoN const * const restrict A, SoN * restrict ris, int * restrict sign)
  {
  int i, imax, j, k;
  double big, temp;
  double sum, dum;
  double vv[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));

  imax=0;
  equal_SoN(ris, A);

  (*sign)=1;
  for(i=0; i<NCOLOR; i++)
     {
     big=0.0;
     for(j=0; j<NCOLOR; j++)
        {
        temp = fabs(ris->comp[m(i,j)]);
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

        temp = vv[i]*fabs(sum);
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
       dum=1.0/(ris->comp[m(j,j)]);
       for(i=j+1; i<NCOLOR; i++)
          {
          (ris->comp[m(i,j)])*=dum;
          }
       }
     }
  }


// determinant
double det_SoN(SoN const * const A);


// gives 0 if the matrix is in SU(N) and 1 otherwise
int scheck_SoN(SoN const * const restrict A)
  {
  int i, j, k, ris;
  double aux;

  ris=0;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           aux+=(A->comp[m(i,k)])*(A->comp[m(j,k)]);
           }
        if(i==j) aux-=(1.0);
        if(fabs(aux)>MIN_VALUE) ris=1;
        }
     }

  if(ris==0)
    {
    if(fabs(det_SoN(A)-1)>MIN_VALUE)
      {
      ris=1;
      }
    }

  return ris;
  }


// sunitarize
void unitarize_SoN(SoN * restrict A)
  {
  double check;
  SoN force, guess, guess_old, helper, helper1, helper2;

  if(scheck_SoN(A)!=0)
    {
    one_SoN(&guess);
    check=1.0;

    equal_dag_SoN(&force, A);

    while(check>MIN_VALUE)
         {
         equal_SoN(&guess_old, &guess);

         // maximize Tr(force*guess)
         cool_SoN(&guess, &force);

         // compute the check
         equal_SoN(&helper, &guess);
         minus_equal_SoN(&helper, &guess_old);
         equal_SoN(&helper1, &helper);
         times_SoN(&helper2, &helper, &helper1);
         check=sqrt(fabs(retr_SoN(&helper2))/(double) NCOLOR);

         //printf("aux: %g\n", check);
         }
    equal_SoN(A, &guess);
    }
  }


// takes the traceless antihermitian part
void ta_SoN(SoN *A);


// eponential of the traceless antihermitian part
void taexp_SoN(SoN *A);


// return 0 if matrix is traceless antihermitian, 1 otherwise
int ta_check_SoN(const SoN * const A);


// exponential of a TA matrix
void exp_of_ta_SoN(SoN *A);


// print on screen
void print_on_screen_SoN(SoN const * const A)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        printf("%.16f ", A->comp[m(i,j)]);
        }
     }
  printf("\n");
  }


// print on file
int print_on_file_SoN(FILE *fp, SoN const * const A)
  {
  int i, j, err;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fprintf(fp, "%.16f ", A->comp[m(i,j)]);
        if(err<0)
          {
          fprintf(stderr, "Problem in writing on file a SoN matrix (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        }
     }
  fprintf(fp, "\n");

  return 0;
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_SoN(FILE *fp, SoN const * const A)
  {
  int i, j;
  size_t err;
  double re;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        re=A->comp[m(i,j)];

        err=fwrite(&re, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SoN matrix (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        }
     }

  return 0;
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_SoN(FILE *fp, SoN const * const A)
  {
  int i, j;
  size_t err;
  double re;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        re=A->comp[m(i,j)];

        SwapBytesDouble(&re);

        err=fwrite(&re, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SoN matrix (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        }
     }

  return 0;
  }


// print on binary file in bigendian
int print_on_binary_file_bigen_SoN(FILE *fp, SoN const * const A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_SoN(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_SoN(fp, A);
    }

  return err;
  }


// read from file
int read_from_file_SoN(FILE *fp, SoN *A)
  {
  int i, j, err;
  double re;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fscanf(fp, "%lg", &re);
        if(err!=1)
          {
          fprintf(stderr, "Problems reading SoN matrix from file (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        A->comp[m(i,j)]=re;
        }
     }
  return 0;
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_SoN(FILE *fp, SoN *A)
  {
  size_t err;
  int i, j;
  double re;

  err=0;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err+=fread(&re, sizeof(double), 1, fp);

        A->comp[m(i,j)]=re;
        }
     }

  if(err!=NCOLOR*NCOLOR)
    {
    fprintf(stderr, "Problems reading SoN matrix from file (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }

  return 0;
  }


// read from binary file changing endianness
int read_from_binary_file_swap_SoN(FILE *fp, SoN *A)
  {
  int i, j;
  size_t err;
  double re;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=0;
        err+=fread(&re, sizeof(double), 1, fp);
        if(err!=1)
         {
         fprintf(stderr, "Problems reading SoN matrix from file (%s, %d)\n", __FILE__, __LINE__);
         return 1;
         }

        SwapBytesDouble(&re);

        A->comp[m(i,j)]=re;
        }
     }
  return 0;
  }


// read from binary file written in bigendian
int read_from_binary_file_bigen_SoN(FILE *fp, SoN *A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_SoN(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_SoN(fp, A);
    }

  return err;
  }


// initialize tensor product
void TensProd_init_SoN(TensProd *TP, SoN const * const A1, SoN const * const A2);



// ***************** for SoNAdj



// convert the fundamental representation matrix B to the adjoint representation matrix A
void fund_to_adj_SoN(SoNAdj * restrict A, SoN const * const restrict B);


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
void TensProdAdj_init_SoN(TensProdAdj * restrict TP, SoN const * const restrict A1, SoN const * const restrict A2);


// initialize tensor product in the adjoint representation
// using two matrices in the adjoint representation
void TensProdAdj_init_SoNAdj(TensProdAdj * restrict TP, SoNAdj const * const restrict A1, SoNAdj const * const restrict A2);


// A=1
void one_SoNAdj(SoNAdj * restrict A);


// A=0
void zero_SoNAdj(SoNAdj * restrict A);


// A+=B
void plus_equal_SoNAdj(SoNAdj * restrict A, SoNAdj const * const restrict B);


// A*=r
void times_equal_real_SoNAdj(SoNAdj * restrict A, double r);


// A*=B
void times_equal_SoNAdj(SoNAdj * restrict A, SoNAdj const * const restrict B);


// trace in the adjoint rep.
double retr_SoNAdj(SoNAdj * restrict A);


// ***************** for SoNVecs


// A=1
void one_SoNVecs(SoNVecs * restrict A);


// A=0
void zero_SoNVecs(SoNVecs * restrict A);


// A=B
void equal_SoNVecs(SoNVecs * restrict A, SoNVecs const * const restrict B);


// A -> A^{\dag}
void conjugate_SoNVecs(SoNVecs * restrict A);


// A-=B
void minus_equal_SoNVecs(SoNVecs * restrict A, SoNVecs const * const restrict B);


// A+=B
void plus_equal_SoNVecs(SoNVecs * restrict A, SoNVecs const * const restrict B);


// *= with real number
void times_equal_real_SoNVecs(SoNVecs * restrict A, double r);

// *= with double for a single component
void times_equal_real_single_SoNVecs(SoNVecs * restrict A, double r, int j);

// *= with complex number for a single component
void times_equal_complex_single_SoNVecs(SoNVecs * restrict A, double complex r, int j);


// norm
double norm_SoNVecs(SoNVecs const * const restrict A);


// normalize
void normalize_SoNVecs(SoNVecs * restrict A);


// random vector (normalized)
void rand_vec_SoNVecs(SoNVecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     A->comp[i] = 1.0-2.0*casuale();
     }

  normalize_SoNVecs(A);
  }

// real part of the scalar product re(v_1^{\dag}v_2)
double re_scal_prod_SoNVecs(SoNVecs const * const restrict v1, SoNVecs const * const restrict v2);


// real part of the scalar product re(v_1[a]^{\dag}v_2[b]) with a, b flavour indices
double re_scal_prod_single_SoNVecs(SoNVecs const * const restrict v1, SoNVecs const * const restrict v2, int a, int b);


// the i-th component of v2 is multiplied by "matrix"
// v1=matrix*v2
void matrix_times_vector_single_SoNVecs(SoNVecs * restrict v1,
                                        SoN const * const restrict matrix,
                                        SoNVecs const * const restrict v2,
                                        int i);


// all the components of v2 are multiplied by "matrix"
// v1=matrix*v2
void matrix_times_vector_all_SoNVecs(SoNVecs * restrict v1,
                                     SoN const * const restrict matrix,
                                     SoNVecs const * const restrict v2);


// rotate two components of the vector
void rotate_two_components_SoNVecs(SoNVecs * restrict v1,
                                   SoNVecs const * const restrict v2,
                                   int i,
                                   int j,
                                   double angle);


// tensor product of two vectors
// Re(v1^{\dag} * aux * v2) = ReTr(aux * matrix)
void vector_tensor_vector_SoNVecs(SoN * restrict matrix,
                                  SoNVecs const * const restrict v1,
                                  SoNVecs const * const restrict v2);


// initialize the flavour matrix with a vector
// FM[mf(i,j)]=\sum_{on_gauge}conj(v1[i])v1[j] - delta^{ij}/N
// i, j are the flavour indices
void init_FMatrix_SoNVecs(FMatrix * restrict fmatrix, SoNVecs const * const restrict v1);


// return a double coumplex number to check the fate of U(1) flavour symmetry
double complex HiggsU1Obs_SoNVecs(SoNVecs const * const restrict v1);


// print on file
int print_on_file_SoNVecs(FILE *fp, SoNVecs const * const A)
  {
  int i, j, err;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fprintf(fp, "%.16f ", A->comp[NCOLOR*i+j]);
        if(err<0)
          {
          fprintf(stderr, "Problem in writing on file a SoN vector (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        }
     }
  fprintf(fp, "\n");

  return 0;
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_SoNVecs(FILE *fp, SoNVecs const * const A)
  {
  int i, j;
  size_t err;
  double re;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        re=A->comp[NCOLOR*i+j];

        err=fwrite(&re, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SoN vector (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        }
     }

  return 0;
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_SoNVecs(FILE *fp, SoNVecs const * const A)
  {
  int i, j;
  size_t err;
  double re;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        re=A->comp[NCOLOR*i+j];

        SwapBytesDouble(&re);

        err=fwrite(&re, sizeof(double), 1, fp);
        if(err!=1)
          {
          fprintf(stderr, "Problem in binary writing on file a SoN vector (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        }
     }

  return 0;
  }


// print on binary file in bigendian
int print_on_binary_file_bigen_SoNVecs(FILE *fp, SoNVecs const * const A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_SoNVecs(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_SoNVecs(fp, A);
    }

  return err;
  }


// read from file
int read_from_file_SoNVecs(FILE *fp, SoNVecs *A)
  {
  int i, j, err;
  double re;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fscanf(fp, "%lg", &re);
        if(err!=1)
          {
          fprintf(stderr, "Problems reading SoN vector from file (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }
        A->comp[NCOLOR*i+j]=re;
        }
     }
  return 0;
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_SoNVecs(FILE *fp, SoNVecs *A)
  {
  size_t err;
  int i, j;
  double re;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fread(&re, sizeof(double), 1, fp);

        if(err!=1)
          {
          fprintf(stderr, "Problems reading SoN vector from file (%s, %d)\n", __FILE__, __LINE__);
          return 1;
          }

        A->comp[NCOLOR*i+j]=re;
        }
     }

  return 0;
  }


// read from binary file changing endianness
int read_from_binary_file_swap_SoNVecs(FILE *fp, SoNVecs *A)
  {
  int i, j;
  size_t err;
  double re;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fread(&re, sizeof(double), 1, fp);

        if(err!=1)
         {
         fprintf(stderr, "Problems reading SoN vector from file (%s, %d)\n", __FILE__, __LINE__);
         return 1;
         }

        SwapBytesDouble(&re);

        A->comp[NCOLOR*i+j]=re;
        }
     }
  return 0;
  }


// read from binary file written in bigendian
int read_from_binary_file_bigen_SoNVecs(FILE *fp, SoNVecs *A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_SoNVecs(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_SoNVecs(fp, A);
    }

  return err;
  }


#endif

