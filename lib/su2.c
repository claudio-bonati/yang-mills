#ifndef SU2_C
#define SU2_C

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/su2.h"

//
// An Su2 matrix is represented as comp[0]+i\sum_{j=1}^3 comp[j]\sigma_j where
// sigma_j are Pauli matrices, comp[j] are real and \sum_{j=0}^3 comp[j]^2=1
//

void init_Su2(Su2 *A, double vec[4]);

// A=1
void one_Su2(Su2 *A);

// A=0
void zero_Su2(Su2 *A);

// A=B
void equal_Su2(Su2 *A, Su2 const * const B);

// A=B^{dag}
void equal_dag_Su2(Su2 *A, Su2 const * const B);

// A+=B
void plus_equal_Su2(Su2 *A, Su2 const * const B);


// A+=B^{dag}
void plus_equal_dag_Su2(Su2 *A, Su2 const * const B);


// A-=B
void minus_equal_Su2(Su2 *A, Su2 const * const B);


// A-=(r*B)
void minus_equal_times_real_Su2(Su2 *A, Su2 const * const B, double r);


// A-=B^{dag}
void minus_equal_dag_Su2(Su2 *A, Su2 const * const B);

// A=b*B+c*C
void lin_comb_Su2(Su2 *A,
                  double b, Su2 const * const B,
                  double c, Su2 const * const C);

// A=b*B^{dag}+c*C
void lin_comb_dag1_Su2(Su2 *A,
                       double b, Su2 const * const B,
                       double c, Su2 const * const C);


// A=b*B+c*C^{dag}
void lin_comb_dag2_Su2(Su2 *A,
                       double b, Su2 const * const B,
                       double c, Su2 const * const C);


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_Su2(Su2 *A,
                        double b, Su2 const * const B,
                        double c, Su2 const * const C);

// A*=r
void times_equal_real_Su2(Su2 *A, double r);

// A*=r
void times_equal_complex_Su2(Su2 *A, double complex r);


// A*=B
void times_equal_Su2(Su2 *A, Su2 const * const B);


// A*=B^{dag}
void times_equal_dag_Su2(Su2 *A, Su2 const * const B);


// A=B*C
void times_Su2(Su2 *A, Su2 const * const B, Su2 const * const C);


// A=B^{dag}*C
void times_dag1_Su2(Su2 *A, Su2 const * const B, Su2 const * const C);


// A=B*C^{dag}
void times_dag2_Su2(Su2 *A, Su2 const * const B, Su2 const * const C);


// A=B^{dag}*C^{dag}
void times_dag12_Su2(Su2 *A, Su2 const * const B, Su2 const * const C);


// random SU(2) matrix
void rand_matrix_Su2(Su2 * restrict A)
  {
  register double p0, p1, p2, p3, p;

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

  A->comp[0]=p0;
  A->comp[1]=p1;
  A->comp[2]=p2;
  A->comp[3]=p3;
  }


// random SU(2) matrix with p0 given (used in the update)
void rand_matrix_p0_Su2(double p0, Su2 * restrict A)
  {
  register double p1, p2, p3, p;

  p=2.0;
  while(p>1.0)
       {
       p1=1.0-2.0*casuale();
       p2=1.0-2.0*casuale();
       p3=1.0-2.0*casuale();
       p=p1*p1+p2*p2+p3*p3;
       }

  p/=(1.0-p0*p0);
  p=sqrt(p);

  p1/=p;
  p2/=p;
  p3/=p;

  A->comp[0]=p0;
  A->comp[1]=p1;
  A->comp[2]=p2;
  A->comp[3]=p3;
  }


// sqrt of the determinant
double sqrtdet_Su2(Su2 const * const A);


// l2 norm of the matrix
double norm_Su2(Su2 const * const A);


// real part of the trace /2
double retr_Su2(Su2 const * const A);


// imaginary part of the trace /2
double imtr_Su2(Su2 const * const A);


// unitarize the matrix
void unitarize_Su2(Su2 *A);


// traceless antihermitian part
void ta_Su2(Su2 *A);


// exponential of the traceless antihermitian part
void taexp_Su2(Su2 *A);


// print on screen
void print_on_screen_Su2(Su2 const * const restrict A)
  {
  printf("%.16lf %.16lf %.16lf %.16lf\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);
  }


// print on file
void print_on_file_Su2(FILE *fp, Su2 const * const restrict A)
  {
  int err;
  err=fprintf(fp, "%.16lf %.16lf %.16lf %.16lf\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);
  if(err<0)
    {
    fprintf(stderr, "Problem in writing on a file a Su2 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file without changing endiannes
void print_on_binary_file_noswap_Su2(FILE *fp, Su2 const * const restrict A)
  {
  size_t err=0;
  err+=fwrite(&(A->comp[0]), sizeof(double), 1, fp);
  err+=fwrite(&(A->comp[1]), sizeof(double), 1, fp);
  err+=fwrite(&(A->comp[2]), sizeof(double), 1, fp);
  err+=fwrite(&(A->comp[3]), sizeof(double), 1, fp);
  if(err!=4)
    {
    fprintf(stderr, "Problem in binary writing on a file a Su2 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file changing endiannes
void print_on_binary_file_swap_Su2(FILE *fp, Su2 const * const restrict A)
  {
  double tmp;
  size_t err=0;

  tmp=A->comp[0];
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);

  tmp=A->comp[1];
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);

  tmp=A->comp[2];
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);

  tmp=A->comp[3];
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);
  if(err!=4)
    {
    fprintf(stderr, "Problem in binary writing on a file a Su2 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file in big endian format
void print_on_binary_file_bigen_Su2(FILE *fp, Su2 const * const restrict A)
  {
  if(endian()==0) // little endian machine
    {
    print_on_binary_file_swap_Su2(fp, A);
    }
  else
    {
    print_on_binary_file_noswap_Su2(fp, A);
    }
  }


// read from file
void read_from_file_Su2(FILE *fp, Su2 * restrict A)
  {
  int err=fscanf(fp, "%lg %lg %lg %lg", &(A->comp[0]), &(A->comp[1]), &(A->comp[2]), &(A->comp[3]));

  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// read from binary file without changing endiannes
void read_from_binary_file_noswap_Su2(FILE *fp, Su2 * restrict A)
  {
  size_t err=0;
  err+=fread(&(A->comp[0]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[1]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[2]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[3]), sizeof(double), 1, fp);
  
  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// read from binary file changing endiannes
void read_from_binary_file_swap_Su2(FILE *fp, Su2 * restrict A)
  {
  size_t err=0;
  err+=fread(&(A->comp[0]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[1]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[2]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[3]), sizeof(double), 1, fp);
  
  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  SwapBytesDouble((void *)&(A->comp[0]));
  SwapBytesDouble((void *)&(A->comp[1]));
  SwapBytesDouble((void *)&(A->comp[2]));
  SwapBytesDouble((void *)&(A->comp[3]));
  }


// read from binary file written in big endian
void read_from_binary_file_bigen_Su2(FILE *fp, Su2 * restrict A)
  {
  if(endian()==0) // little endian machine
    {
    read_from_binary_file_swap_Su2(fp, A);
    }
  else
    {
    read_from_binary_file_noswap_Su2(fp, A);
    }
  }


// initialize tensor product
void TensProd_init_Su2(TensProd *TP, Su2 const * const A1, Su2 const * const A2);

#endif
