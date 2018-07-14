#ifndef U1_C
#define U1_C

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/tens_prod.h"
#include"../include/u1.h"


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
void print_on_file_U1(FILE *fp, U1 const * const A)
  {
  int err;
  err=fprintf(fp, "%.16lf %.16lf\n", creal(A->comp), cimag(A->comp));
  if(err<0)
    {
    fprintf(stderr, "Problem in writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file without changing endiannes
void print_on_binary_file_noswap_U1(FILE *fp, U1 const * const A)
  {
  size_t err=0;
  double re, im;

  re=creal(A->comp);
  im=cimag(A->comp);

  err=fwrite(&re, sizeof(double), 1, fp);
  if(err!=1)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=fwrite(&im, sizeof(double), 1, fp);
  if(err!=1)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file changing endiannes
void print_on_binary_file_swap_U1(FILE *fp, U1 const * const A)
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
    exit(EXIT_FAILURE);
    }
  }


// print on binary file in big endian format
void print_on_binary_file_bigen_U1(FILE *fp, const U1 * const A)
  {
  if(endian()==0) // little endian machine
    {
    print_on_binary_file_swap_U1(fp, A);
    }
  else
    {
    print_on_binary_file_noswap_U1(fp, A);
    }
  }


// read from file
void read_from_file_U1(FILE *fp, U1 *A)
  {
  double re, im;
  int err;

  err=fscanf(fp, "%lg %lg", &re, &im);
  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  A->comp=re+im*I;
  }


// read from binary file without changing endiannes
void read_from_binary_file_noswap_U1(FILE *fp, U1 *A)
  {
  double re, im;
  size_t err=0;
  err+=fread(&re, sizeof(double), 1, fp);
  err+=fread(&im, sizeof(double), 1, fp);
  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  A->comp=re+I*im;
  }


// read from binary file changing endiannes
void read_from_binary_file_swap_U1(FILE *fp, U1 *A)
  {
  double re, im;
  size_t err=0;

  err+=fread(&re, sizeof(double), 1, fp);
  err+=fread(&im, sizeof(double), 1, fp);

  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  SwapBytesDouble((void *)&re);
  SwapBytesDouble((void *)&im);

  A->comp=re+I*im;
  }


// read from binary file written in big endian
void read_from_binary_file_bigen_U1(FILE *fp, U1 *A)
  {
  if(endian()==0) // little endian machine
    {
    read_from_binary_file_swap_U1(fp, A);
    }
  else
    {
    read_from_binary_file_noswap_U1(fp, A);
    }
  }


// initialize tensor product
void TensProd_init_U1(TensProd *TP, U1 const * const A1, U1 const * const A2);



#endif
