#ifndef SU2_C
#define SU2_C

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/flavour_matrix.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/su2.h"

//
// An Su2 matrix is represented as comp[0]+i\sum_{j=1}^3 comp[j]\sigma_j where
// sigma_j are Pauli matrices, comp[j] are real and \sum_{j=0}^3 comp[j]^2=1
//


// ***************** for Su2


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


// A = lambda*B with lambda diagonal marix
void diag_matrix_times_Su2(Su2 * restrict A, double const lambda[2], Su2 const * const restrict B);

// A=lambda*B^{dag} with lambda diagonal matrix
void diag_matrix_times_dag_Su2(Su2 * restrict A, double const lambda[2], Su2 const * const restrict B);


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
int print_on_file_Su2(FILE *fp, Su2 const * const restrict A)
  {
  int err;

  err=fprintf(fp, "%.16lf %.16lf %.16lf %.16lf\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);

  if(err<0)
    {
    fprintf(stderr, "Problem in writing on a file a Su2 matrix (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }
  else
    {
    return 0;
    }
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_Su2(FILE *fp, Su2 const * const restrict A)
  {
  size_t err=0;
  err+=fwrite(&(A->comp[0]), sizeof(double), 1, fp);
  err+=fwrite(&(A->comp[1]), sizeof(double), 1, fp);
  err+=fwrite(&(A->comp[2]), sizeof(double), 1, fp);
  err+=fwrite(&(A->comp[3]), sizeof(double), 1, fp);

  if(err!=4)
    {
    fprintf(stderr, "Problem in binary writing on a file a Su2 matrix (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }
  else
    {
    return 0;
    }
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_Su2(FILE *fp, Su2 const * const restrict A)
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
    return 1;
    }
  else
    {
    return 0;
    }
  }


// print on binary file in big endian format
int print_on_binary_file_bigen_Su2(FILE *fp, Su2 const * const restrict A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_Su2(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_Su2(fp, A);
    }
  return err;
  }


// read from file
int read_from_file_Su2(FILE *fp, Su2 * restrict A)
  {
  int err=fscanf(fp, "%lg %lg %lg %lg", &(A->comp[0]), &(A->comp[1]), &(A->comp[2]), &(A->comp[3]));

  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }
  else
    {
    return 0;
    }
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_Su2(FILE *fp, Su2 * restrict A)
  {
  size_t err=0;
  err+=fread(&(A->comp[0]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[1]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[2]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[3]), sizeof(double), 1, fp);
  
  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }
  else
    {
    return 0;
    }
  }


// read from binary file changing endiannes
int read_from_binary_file_swap_Su2(FILE *fp, Su2 * restrict A)
  {
  size_t err=0;
  err+=fread(&(A->comp[0]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[1]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[2]), sizeof(double), 1, fp);
  err+=fread(&(A->comp[3]), sizeof(double), 1, fp);
  
  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    return 1;
    }

  SwapBytesDouble((void *)&(A->comp[0]));
  SwapBytesDouble((void *)&(A->comp[1]));
  SwapBytesDouble((void *)&(A->comp[2]));
  SwapBytesDouble((void *)&(A->comp[3]));

  return 0;
  }


// read from binary file written in big endian
int read_from_binary_file_bigen_Su2(FILE *fp, Su2 * restrict A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_Su2(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_Su2(fp, A);
    }

  return err;
  }


// initialize tensor product
void TensProd_init_Su2(TensProd *TP, Su2 const * const A1, Su2 const * const A2);



// ***************** for Su2Vecs


// A=1
void one_Su2Vecs(Su2Vecs * restrict A);


// A=0
void zero_Su2Vecs(Su2Vecs * restrict A);


// A=B
void equal_Su2Vecs(Su2Vecs * restrict A, Su2Vecs const * const restrict B);


// A -> A^{dag}
void conjugate_Su2Vecs(Su2Vecs * restrict A);


// A-=B
void minus_equal_Su2Vecs(Su2Vecs * restrict A, Su2Vecs const * const restrict B);


// A+=B
void plus_equal_Su2Vecs(Su2Vecs * restrict A, Su2Vecs const * const restrict B);


// *= with real number
void times_equal_real_Su2Vecs(Su2Vecs * restrict A, double r);


// *= with real for a single component
void times_equal_real_single_Su2Vecs(Su2Vecs * restrict A, double r, int j);


// *= with complex number for a single component
void times_equal_complex_single_Su2Vecs(Su2Vecs * restrict A, double complex r, int j);


// norm
double norm_Su2Vecs(Su2Vecs const * const restrict A);


// normalize
void normalize_Su2Vecs(Su2Vecs * restrict A);


// random vector (normalized)
void rand_vec_Su2Vecs(Su2Vecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double p0;

  for(i=0; i<4*NHIGGS; i++)
     {
     p0=1.0-2.0*casuale();
     A->comp[i] = p0;
     }

  normalize_Su2Vecs(A);
  }


// real part of the scalar product re(v_1^{\dag}v_2)
double re_scal_prod_Su2Vecs(Su2Vecs const * const restrict v1, Su2Vecs const * const restrict v2);


// real part of the scalar product re(v_1[a]^{\dag}v_2[b]) with a, b flavour indices
double re_scal_prod_single_Su2Vecs(Su2Vecs const * const restrict v1, Su2Vecs const * const restrict v2, int a, int b);


// the i-th component of v2 is multiplied by "matrix"
// v1=matrix*v2
void matrix_times_vector_single_Su2Vecs(Su2Vecs * restrict v1, Su2 const * const restrict matrix, Su2Vecs const * const restrict v2, int i);


// all the components of v2 are multiplied by "matrix"
// v1=matrix*v2
void matrix_times_vector_all_Su2Vecs(Su2Vecs * restrict v1, Su2 const * const restrict matrix, Su2Vecs const * const restrict v2);


// rotate two components of the vector
void rotate_two_components_Su2Vecs(Su2Vecs * restrict v1,
                                   Su2Vecs const * const restrict v2,
                                   int i,
                                   int j,
                                   double angle);


// tensor product of two vectors
// Re(v1^{\dag} * aux * v2) = ReTr(aux * matrix)
void vector_tensor_vector_Su2Vecs(Su2 * restrict matrix, Su2Vecs const * const restrict v1, Su2Vecs const * const restrict v2);


// initialize the flavour matrix with a vector
// FM[mf(i,j)]=\sum_{on_gauge}conj(v1[i])v1[j] - delta^{ij}/N
// i, j are the flavour indices
void init_FMatrix_Su2Vecs(FMatrix * restrict fmatrix, Su2Vecs const * const restrict v1);


// return a double coumplex number to check the fate of U(1) flavour symmetry
double complex HiggsU1Obs_Su2Vecs(Su2Vecs const * const restrict v1);


// print on file
int print_on_file_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A)
  {
  int i, err;

  for(i=0; i<NHIGGS; i++)
     {
     err=fprintf(fp, "%.16lf %.16lf %.16lf %.16lf\n", A->comp[4*i+0], A->comp[4*i+1], A->comp[4*i+2], A->comp[4*i+3]);

     if(err<0)
       {
       fprintf(stderr, "Problem in writing on a file a Su2 vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A)
  {
  int i;
  size_t err;

  for(i=0; i<NHIGGS; i++)
     {
     err =fwrite(&(A->comp[4*i+0]), sizeof(double), 1, fp);
     err+=fwrite(&(A->comp[4*i+1]), sizeof(double), 1, fp);
     err+=fwrite(&(A->comp[4*i+2]), sizeof(double), 1, fp);
     err+=fwrite(&(A->comp[4*i+3]), sizeof(double), 1, fp);

     if(err!=4)
       {
       fprintf(stderr, "Problem in binary writing on a file a Su2 vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A)
  {
  int i;
  double tmp;
  size_t err;

  for(i=0; i<NHIGGS; i++)
     {
     tmp=A->comp[4*i+0];
     SwapBytesDouble(&tmp);
     err=fwrite(&(tmp), sizeof(double), 1, fp);

     tmp=A->comp[4*i+1];
     SwapBytesDouble(&tmp);
     err+=fwrite(&(tmp), sizeof(double), 1, fp);

     tmp=A->comp[4*i+2];
     SwapBytesDouble(&tmp);
     err+=fwrite(&(tmp), sizeof(double), 1, fp);

     tmp=A->comp[4*i+3];
     SwapBytesDouble(&tmp);
     err+=fwrite(&(tmp), sizeof(double), 1, fp);

     if(err!=4)
       {
       fprintf(stderr, "Problem in binary writing on a file a Su2 vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file in big endian format
int print_on_binary_file_bigen_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_Su2Vecs(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_Su2Vecs(fp, A);
    }
  return err;
  }


// read from file
int read_from_file_Su2Vecs(FILE *fp, Su2Vecs * restrict A)
  {
  int i, err;

  for(i=0; i<NHIGGS; i++)
     {
     err=fscanf(fp, "%lg %lg %lg %lg", &(A->comp[4*i+0]), &(A->comp[4*i+1]), &(A->comp[4*i+2]), &(A->comp[4*i+3]));

     if(err!=4)
       {
       fprintf(stderr, "Problems reading Su2 vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_Su2Vecs(FILE *fp, Su2Vecs * restrict A)
  {
  int i;
  size_t err;

  for(i=0; i<NHIGGS; i++)
     {
     err =fread(&(A->comp[4*i+0]), sizeof(double), 1, fp);
     err+=fread(&(A->comp[4*i+1]), sizeof(double), 1, fp);
     err+=fread(&(A->comp[4*i+2]), sizeof(double), 1, fp);
     err+=fread(&(A->comp[4*i+3]), sizeof(double), 1, fp);

     if(err!=4)
       {
       fprintf(stderr, "Problems reading Su2 vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// read from binary file changing endiannes
int read_from_binary_file_swap_Su2Vecs(FILE *fp, Su2Vecs * restrict A)
  {
  int i;
  size_t err;

  for(i=0; i<NHIGGS; i++)
     {
     err =fread(&(A->comp[4*i+0]), sizeof(double), 1, fp);
     err+=fread(&(A->comp[4*i+1]), sizeof(double), 1, fp);
     err+=fread(&(A->comp[4*i+2]), sizeof(double), 1, fp);
     err+=fread(&(A->comp[4*i+3]), sizeof(double), 1, fp);

     if(err!=4)
       {
       fprintf(stderr, "Problems reading Su2 vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     SwapBytesDouble((void *)&(A->comp[4*i+0]));
     SwapBytesDouble((void *)&(A->comp[4*i+1]));
     SwapBytesDouble((void *)&(A->comp[4*i+2]));
     SwapBytesDouble((void *)&(A->comp[4*i+3]));
     }

  return 0;
  }


// read from binary file written in big endian
int read_from_binary_file_bigen_Su2Vecs(FILE *fp, Su2Vecs * restrict A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_Su2Vecs(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_Su2Vecs(fp, A);
    }

  return err;
  }


#endif
