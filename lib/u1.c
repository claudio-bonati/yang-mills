#ifndef U1_C
#define U1_C

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/flavour_matrix.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/tens_prod.h"
#include"../include/u1.h"


// ***************** for U1


// initialize
void init_U1(U1 *A, double complex vec);


// A=1
void one_U1(U1 *A);


// A=0
void zero_U1(U1 *A);


// A=B
void equal_U1(U1 *A, U1 const * const B);


// A=B^{dag}
void equal_dag_U1(U1 *A, U1 const * const B);


// A+=B
void plus_equal_U1(U1 *A, U1 const * const B);


// A+=B^{dag}
void plus_equal_dag_U1(U1 *A, U1 const * const B);


// A-=B
void minus_equal_U1(U1 *A, U1 const * const B);


// A-=(r*B)
void minus_equal_times_real_U1(U1 *A, U1 const * const B, double r);


// A-=B^{dag}
void minus_equal_dag_U1(U1 *A, U1 const * const B);


// A=b*B+c*C
void lin_comb_U1(U1 *A,
                  double b, U1 const * const B,
                  double c, U1 const * const C);


// A=b*B^{dag}+c*C
void lin_comb_dag1_U1(U1 *A,
                       double b, U1 const * const B,
                       double c, U1 const * const C);


// A=b*B+c*C^{dag}
void lin_comb_dag2_U1(U1 *A,
                       double b, U1 const * const B,
                       double c, U1 const * const C);


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_U1(U1 *A,
                        double b, U1 const * const B,
                        double c, U1 const * const C);


// A*=r
void times_equal_real_U1(U1 *A, double r);


// A*=r
void times_equal_complex_U1(U1 *A, double complex r);


// A*=B
void times_equal_U1(U1 *A, U1 const * const B);


// A*=B^{dag}
void times_equal_dag_U1(U1 *A, U1 const * const B);


// A=B*C
void times_U1(U1 *A, U1 const * const B, U1 const * const C);


// A=lambda*B with lambda diagonal matrix
void diag_matrix_times_U1(U1 * restrict A, double const lambda[1], U1 const * const restrict B);


// A=lambda*B^{dag} with lambda diagonal matrix
void diag_matrix_times_dag_U1(U1 * restrict A, double const lambda[1], U1 const * const restrict B);


// A=B^{dag}*C
void times_dag1_U1(U1 *A, U1 const * const B, U1 const * const C);


// A=B*C^{dag}
void times_dag2_U1(U1 *A, U1 const * const B, U1 const * const C);


// A=B^{dag}*C^{dag}
void times_dag12_U1(U1 *A, U1 const * const B, U1 const * const C);


// random matrix
void rand_matrix_U1(U1 * restrict A)
  {
  double p0, p1, p;

  do
    {
    p0=1.0-2.0*casuale();
    p1=1.0-2.0*casuale();

    p=sqrt(p0*p0+p1*p1);
    }
  while(p<MIN_VALUE);

  p0/=p;
  p1/=p;

  A->comp = p0 + p1*I;
  }


// l2 norm of the matrix
double norm_U1(U1 const * const A);


// real part of the trace
double retr_U1(U1 const * const A);


// imaginary part of the trace
double imtr_U1(U1 const * const A);


// unitarize the matrix
void unitarize_U1(U1 *A);


// antihermitian part (NO TRACELESS!)
void ta_U1(U1 *A);


// exponential of the antihermitian part (NO TRACELESS!)
void taexp_U1(U1 *A);


// print on screen
void print_on_screen_U1(U1 const * const A)
  {
  printf("%.16lf %.16lf\n", creal(A->comp), cimag(A->comp));
  }


// print on file
int print_on_file_U1(FILE *fp, U1 const * const A)
  {
  int err;

  err=fprintf(fp, "%.16lf %.16lf\n", creal(A->comp), cimag(A->comp));
  if(err<0)
    {
    fprintf(stderr, "Problem in writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }

  return 0;
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_U1(FILE *fp, U1 const * const A)
  {
  size_t err=0;
  double re, im;

  re=creal(A->comp);
  im=cimag(A->comp);

  err=fwrite(&re, sizeof(double), 1, fp);
  if(err!=1)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }
  err=fwrite(&im, sizeof(double), 1, fp);
  if(err!=1)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }

  return 0;
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_U1(FILE *fp, U1 const * const A)
  {
  double tmp;
  size_t err=0;

  tmp=creal(A->comp);
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);

  tmp=cimag(A->comp);
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);

  if(err!=2)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }

  return 0;
  }


// print on binary file in big endian format
int print_on_binary_file_bigen_U1(FILE *fp, const U1 * const A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_U1(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_U1(fp, A);
    }

  return err;
  }


// read from file
int read_from_file_U1(FILE *fp, U1 *A)
  {
  double re, im;
  int err;

  err=fscanf(fp, "%lg %lg", &re, &im);
  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }
  A->comp=re+im*I;

  return 0;
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_U1(FILE *fp, U1 *A)
  {
  double re, im;
  size_t err=0;
  err+=fread(&re, sizeof(double), 1, fp);
  err+=fread(&im, sizeof(double), 1, fp);
  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }
  A->comp=re+I*im;

  return 0;
  }


// read from binary file changing endiannes
int read_from_binary_file_swap_U1(FILE *fp, U1 *A)
  {
  double re, im;
  size_t err=0;

  err+=fread(&re, sizeof(double), 1, fp);
  err+=fread(&im, sizeof(double), 1, fp);

  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }

  SwapBytesDouble((void *)&re);
  SwapBytesDouble((void *)&im);

  A->comp=re+I*im;

  return 0;
  }


// read from binary file written in big endian
int read_from_binary_file_bigen_U1(FILE *fp, U1 *A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_U1(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_U1(fp, A);
    }

  return err;
  }


// initialize tensor product
void TensProd_init_U1(TensProd *TP, U1 const * const A1, U1 const * const A2);



// ***************** for U1Adj

// convert the fundamental representation matrix B to the adjoint representation matrix A
void fund_to_adj_U1(U1Adj * restrict A, U1 const * const restrict B);


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
void TensProdAdj_init_U1(TensProdAdj * restrict TP, U1 const * const restrict A1, U1 const * const restrict A2);


// initialize tensor product in the adjoint representation
// using two matrices in the adjoint representation
void TensProdAdj_init_U1Adj(TensProdAdj * restrict TP, U1Adj const * const restrict A1, U1Adj const * const restrict A2);


// A=1
void one_U1Adj(U1Adj * restrict A);


// A=0
void zero_U1Adj(U1Adj * restrict A);


// A+=B
void plus_equal_U1Adj(U1Adj * restrict A, U1Adj const * const restrict B);


// A*=r
void times_equal_real_U1Adj(U1Adj * restrict A, double r);


// A*=B
void times_equal_U1Adj(U1Adj * restrict A, U1Adj const * const restrict B);


// trace in the adjoint rep.
double retr_U1Adj(U1Adj * restrict A);


// ***************** for U1Vecs


// A=1
void one_U1Vecs(U1Vecs * restrict A);


// A=0
void zero_U1Vecs(U1Vecs * restrict A);


// A=B
void equal_U1Vecs(U1Vecs * restrict A, U1Vecs const * const restrict B);


// A -> A^{dag}
void conjugate_U1Vecs(U1Vecs * restrict A);


// A-=B
void minus_equal_U1Vecs(U1Vecs * restrict A, U1Vecs const * const restrict B);


// A+=B
void plus_equal_U1Vecs(U1Vecs * restrict A, U1Vecs const * const restrict B);


// *= with real number
void times_equal_real_U1Vecs(U1Vecs * restrict A, double r);


// *= with real for a single component
void times_equal_real_single_U1Vecs(U1Vecs * restrict A, double r, int j);


// *= with complex number for a single component
void times_equal_complex_single_U1Vecs(U1Vecs * restrict A, double complex r, int j);


// norm
double norm_U1Vecs(U1Vecs const * const restrict A);


// normalize
void normalize_U1Vecs(U1Vecs * restrict A);


// random vector (normalized)
void rand_vec_U1Vecs(U1Vecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double p0, p1;

  for(i=0; i<NHIGGS; i++)
     {
     p0=1.0-2.0*casuale();
     p1=1.0-2.0*casuale();

     A->comp[i] = p0 + p1*I;
     }

  normalize_U1Vecs(A);
  }


// real part of the scalar product re(v_1^{\dag}v_2)
double re_scal_prod_U1Vecs(U1Vecs const * const restrict v1, U1Vecs const * const restrict v2);


// real part of the scalar product re(v_1[a]^{\dag}v_2[b]) with a, b flavour indices
double re_scal_prod_single_U1Vecs(U1Vecs const * const restrict v1, U1Vecs const * const restrict v2, int a, int b);


// the i-th component of v2 is multiplied by "matrix"
// v1=matrix*v2
void matrix_times_vector_single_U1Vecs(U1Vecs * restrict v1, U1 const * const restrict matrix, U1Vecs const * const restrict v2, int i);


// all the components of v2 are multiplied by "matrix"
// v1=matrix*v2
void matrix_times_vector_all_U1Vecs(U1Vecs * restrict v1, U1 const * const restrict matrix, U1Vecs const * const restrict v2);


// tensor product of two vectors
// Re(v1^{\dag} * aux * v2) = ReTr(aux * matrix)
void vector_tensor_vector_U1Vecs(U1 * restrict matrix, U1Vecs const * const restrict v1, U1Vecs const * const restrict v2);


// initialize the flavour matrix with a vector
// FM[mf(i,j)]=conj(v1[i])v1[j]-delta^{ij}/N
void init_FMatrix_U1Vecs(FMatrix * restrict fmatrix, U1Vecs const * const restrict v1);


// return a double coumplex number to check the fate of U(1) flavour symmetry
double complex HiggsU1Obs_U1Vecs(U1Vecs const * const restrict v1);


// print on file
int print_on_file_U1Vecs(FILE *fp, U1Vecs const * const A)
  {
  int i, err;

  for(i=0; i<NHIGGS; i++)
     {
     err=fprintf(fp, "%.16lf %.16lf\n", creal(A->comp[i]), cimag(A->comp[i]));
     if(err<0)
       {
       fprintf(stderr, "Problem in writing on a file a U1 vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_U1Vecs(FILE *fp, U1Vecs const * const A)
  {
  int i;
  size_t err=0;
  double re, im;

  for(i=0; i<NHIGGS; i++)
     {
     re=creal(A->comp[i]);
     im=cimag(A->comp[i]);

     err=fwrite(&re, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on a file a U1 vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     err=fwrite(&im, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on a file a U1 vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_U1Vecs(FILE *fp, U1Vecs const * const A)
  {
  int i;
  double tmp;
  size_t err;

  for(i=0;i<NHIGGS; i++)
     {
     tmp=creal(A->comp[i]);
     SwapBytesDouble(&tmp);
     err=fwrite(&(tmp), sizeof(double), 1, fp);

     tmp=cimag(A->comp[i]);
     SwapBytesDouble(&tmp);
     err+=fwrite(&(tmp), sizeof(double), 1, fp);

     if(err!=2)
       {
       fprintf(stderr, "Problem in binary writing on a file a U1 vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file in big endian format
int print_on_binary_file_bigen_U1Vecs(FILE *fp, const U1Vecs * const A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_U1Vecs(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_U1Vecs(fp, A);
    }

  return err;
  }


// read from file
int read_from_file_U1Vecs(FILE *fp, U1Vecs *A)
  {
  double re, im;
  int i, err;

  for(i=0; i<NHIGGS; i++)
     {
     err=fscanf(fp, "%lg %lg", &re, &im);
     if(err!=2)
       {
       fprintf(stderr, "Problems reading U1 vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     A->comp[i]=re+im*I;
     }

  return 0;
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_U1Vecs(FILE *fp, U1Vecs *A)
  {
  int i;
  double re, im;
  size_t err=0;

  for(i=0; i<NHIGGS; i++)
     {
     err+=fread(&re, sizeof(double), 1, fp);
     err+=fread(&im, sizeof(double), 1, fp);
     if(err!=2)
       {
       fprintf(stderr, "Problems reading U1 vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     A->comp[i]=re+I*im;
     }

  return 0;
  }


// read from binary file changing endiannes
int read_from_binary_file_swap_U1Vecs(FILE *fp, U1Vecs *A)
  {
  int i;
  double re, im;
  size_t err;

  for(i=0; i<NHIGGS; i++)
     {
     err=fread(&re, sizeof(double), 1, fp);
     err+=fread(&im, sizeof(double), 1, fp);

     if(err!=2)
       {
       fprintf(stderr, "Problems reading U1 vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     SwapBytesDouble((void *)&re);
     SwapBytesDouble((void *)&im);

     A->comp[i]=re+I*im;
     }

  return 0;
  }


// read from binary file written in big endian
int read_from_binary_file_bigen_U1Vecs(FILE *fp, U1Vecs *A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_U1Vecs(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_U1Vecs(fp, A);
    }

  return err;
  }


#endif
