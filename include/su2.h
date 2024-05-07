#ifndef SU2_H
#define SU2_H

#include<complex.h>
#include<math.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"macro.h"
#include"tens_prod.h"

//
// An Su2 matrix is represented as comp[0]+i\sum_{j=1}^3 comp[j]\sigma_j where
// sigma_j are Pauli matrices, comp[j] are real and \sum_{j=0}^3 comp[j]^2=1
//
typedef struct Su2 {
     double comp[4] __attribute__((aligned(DOUBLE_ALIGN)));
} Su2;


typedef struct Su2Adj {
     double comp[9] __attribute__((aligned(DOUBLE_ALIGN)));
} Su2Adj;

//
// a complex 2-vector \Phi = (\phi_1+i\phi_2, \phi_3+i\phi_4) can be conveniently rewritten as a 2x2 matrix \varphi
// by using the auxilliary vector \tilde{\Phi}=i\sigma_2\Phi* and (see Montvay Muenster eq. 6.1)
// \varphi=( \tilde{\Phi} | \Phi ) = \phi_3 Id + i\sigma_1\phi_2 +i\sigma_2\phi_1-i\sigma_3\phi_4
// It is not difficult to show that
// \tilde{U\Phi} = U \tilde{\Phi}
// U\varphi=(\tilde{U\Phi}| U\Phi )
// 2 Re (\Phi_A^{\dag} U \Phi_B) = Tr (\varphi_A^{\dag} U \varphi_B)
//
// the 2-dim complex form of the flavour "i" vector is
// (v[4*i+2]+I v[4*i+1], v[4*i+0]-I v[4*i+3])
//
typedef struct Su2Vecs {
     double comp[4*NHIGGS] __attribute__((aligned(DOUBLE_ALIGN)));
} Su2Vecs;


// ***************** for Su2


inline void init_Su2(Su2 * restrict A, double vec[4])
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=vec[0];
  A->comp[1]=vec[1];
  A->comp[2]=vec[2];
  A->comp[3]=vec[3];
  }


// A=1
inline void one_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=1.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=0
inline void zero_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=0.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=B
inline void equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=B->comp[0];
  A->comp[1]=B->comp[1];
  A->comp[2]=B->comp[2];
  A->comp[3]=B->comp[3];
  }


// A=B^{dag}
inline void equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= B->comp[0];
  A->comp[1]=-B->comp[1];
  A->comp[2]=-B->comp[2];
  A->comp[3]=-B->comp[3];
  }


// A+=B
inline void plus_equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]+=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A+=B^{dag}
inline void plus_equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]+=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=B
inline void minus_equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]-=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=(r*B)
inline void minus_equal_times_real_Su2(Su2 * restrict A, Su2 const * const restrict B, double r)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]-=(r*B->comp[0]);
  A->comp[1]-=(r*B->comp[1]);
  A->comp[2]-=(r*B->comp[2]);
  A->comp[3]-=(r*B->comp[3]);
  }


// A-=B^{dag}
inline void minus_equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]-=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A=b*B+c*C
inline void lin_comb_Su2(Su2 * restrict A,
                  double b, Su2 const * const restrict B,
                  double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] + c*C->comp[1];
  A->comp[2]= b*B->comp[2] + c*C->comp[2];
  A->comp[3]= b*B->comp[3] + c*C->comp[3];
  }


// A=b*B^{dag}+c*C
inline void lin_comb_dag1_Su2(Su2 * restrict A,
                       double b, Su2 const * const restrict B,
                       double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] + c*C->comp[1];
  A->comp[2]= -b*B->comp[2] + c*C->comp[2];
  A->comp[3]= -b*B->comp[3] + c*C->comp[3];
  }


// A=b*B+c*C^{dag}
inline void lin_comb_dag2_Su2(Su2 * restrict A,
                       double b, Su2 const * const restrict B,
                       double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] - c*C->comp[1];
  A->comp[2]= b*B->comp[2] - c*C->comp[2];
  A->comp[3]= b*B->comp[3] - c*C->comp[3];
  }


// A=b*B^{dag}+c*C^{dag}
inline void lin_comb_dag12_Su2(Su2 * restrict A,
                        double b, Su2 const * const restrict B,
                        double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] - c*C->comp[1];
  A->comp[2]= -b*B->comp[2] - c*C->comp[2];
  A->comp[3]= -b*B->comp[3] - c*C->comp[3];
  }


// A*=r
inline void times_equal_real_Su2(Su2 * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]*=r;
  A->comp[1]*=r;
  A->comp[2]*=r;
  A->comp[3]*=r;
  }


// A*=r
inline void times_equal_complex_Su2(Su2 * restrict A, double complex r)
  {
  #ifdef DEBUG
  if(fabs(cimag(r))>MIN_VALUE)
    {
    fprintf(stderr, "Trying to multiply SU(2) matrix by a non-real number (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]*=creal(r);
  A->comp[1]*=creal(r);
  A->comp[2]*=creal(r);
  A->comp[3]*=creal(r);
  }


// A*=B
inline void times_equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  register double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];

  A->comp[0]= a0*B->comp[0] - a1*B->comp[1] - a2*B->comp[2] - a3*B->comp[3];
  A->comp[1]= a0*B->comp[1] + a1*B->comp[0] - a2*B->comp[3] + a3*B->comp[2];
  A->comp[2]= a0*B->comp[2] + a2*B->comp[0] + a1*B->comp[3] - a3*B->comp[1];
  A->comp[3]= a0*B->comp[3] + a3*B->comp[0] - a1*B->comp[2] + a2*B->comp[1];
  }


// A*=B^{dag}
inline void times_equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  register double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];

  A->comp[0]=  a0*B->comp[0] + a1*B->comp[1] + a2*B->comp[2] + a3*B->comp[3];
  A->comp[1]= -a0*B->comp[1] + a1*B->comp[0] + a2*B->comp[3] - a3*B->comp[2];
  A->comp[2]= -a0*B->comp[2] + a2*B->comp[0] - a1*B->comp[3] + a3*B->comp[1];
  A->comp[3]= -a0*B->comp[3] + a3*B->comp[0] + a1*B->comp[2] - a2*B->comp[1];
  }


// A=B*C
inline void times_Su2(Su2 * restrict A,
                      Su2 const * const restrict B,
                      Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// A=lambda*B with lambda diagonal matrix
inline void diag_matrix_times_Su2(Su2 * restrict A, 
                                  double const lambda[2],
                                  Su2 const * const restrict B)
  {
  fprintf(stderr, "The function diag_matrix_times_Su2 cannot be defined using the Pauli representation of SU(2) ");
  fprintf(stderr, "(%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);

  (void) A;
  (void) lambda;
  (void) B;
  }


// A=lambda*B^{dag} with lambda diagonal matrix
inline void diag_matrix_times_dag_Su2(Su2 * restrict A,
                                      double const lambda[2],
                                      Su2 const * const restrict B)
  {
  fprintf(stderr, "The function diag_matrix_times_dag_Su2 cannot be defined using the Pauli representation of SU(2) ");
  fprintf(stderr, "(%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);

  (void) A;
  (void) lambda;
  (void) B;
  }


// A=B^{dag}*C
inline void times_dag1_Su2(Su2 * restrict A,
                    Su2 const * const restrict B,
                    Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B*C^{dag}
inline void times_dag2_Su2(Su2 * restrict A,
                    Su2 const * const restrict B,
                    Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B^{dag}*C^{dag}
inline void times_dag12_Su2(Su2 * restrict A,
                     Su2 const * const restrict B,
                     Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// random SU(2) matrix
void rand_matrix_Su2(Su2 *A);


// random SU(2) matrix with p0 given (used in the update)
void rand_matrix_p0_Su2(double p0, Su2 *A);


// sqrt of the determinant
inline double sqrtdet_Su2(Su2 const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  return sqrt(A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1]
             + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3]);
  }


// l2 norm of the matrix
inline double norm_Su2(Su2 const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  return sqrtdet_Su2(A);
  }


// real part of the trace /2
inline double retr_Su2(Su2 const * const restrict A)
  {
  return A->comp[0];
  }


// imaginary part of the trace /2
inline double imtr_Su2(Su2 const * const restrict A)
  {
  (void) A; // to suppress compilation warning of unused variable
  return 0.0;
  }


// unitarize the matrix
inline void unitarize_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double p;

  p=A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1] + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3];
  p=1.0/sqrt(p);

  A->comp[0]*=p;
  A->comp[1]*=p;
  A->comp[2]*=p;
  A->comp[3]*=p;
  }


// traceless antihermitian part
inline void ta_Su2(Su2 * restrict A)
  {
  A->comp[0]=0;
  }


// exponential of the traceless antihermitian part
inline void taexp_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double v1, v2, v3, norm, s;

  // comp[0] is neglected since we consider the ta part
  norm=(A->comp[1]*A->comp[1]);
  norm+=(A->comp[2]*A->comp[2]);
  norm+=(A->comp[3]*A->comp[3]);
  norm=sqrt(norm);

  v1=A->comp[1]/norm;
  v2=A->comp[2]/norm;
  v3=A->comp[3]/norm;

  s=sin(norm);

  A->comp[0]=cos(norm);
  A->comp[1]=v1*s;
  A->comp[2]=v2*s;
  A->comp[3]=v3*s;
  }


// print on screen
void print_on_screen_Su2(Su2 const * const A);


// print on file
int print_on_file_Su2(FILE *fp, Su2 const * const A);


// print on binary file without changing endiannes
int print_on_binary_file_noswap_Su2(FILE *fp, Su2 const * const A);


// print on binary file changing endiannes
int print_on_binary_file_swap_Su2(FILE *fp, Su2 const * const A);


// print on binary file in big endian format
int print_on_binary_file_bigen_Su2(FILE *fp, Su2 const * const A);


// read from file
int read_from_file_Su2(FILE *fp, Su2 *A);


// read from binary file without changing endiannes
int read_from_binary_file_noswap_Su2(FILE *fp, Su2 *A);


// read from binary file changing endiannes
int read_from_binary_file_swap_Su2(FILE *fp, Su2 *A);


// read from binary file written in big endian
int read_from_binary_file_bigen_Su2(FILE *fp, Su2 *A);


// initialize tensor product
inline void TensProd_init_Su2(TensProd * restrict TP, Su2 const * const restrict A1, Su2 const * const restrict A2)
  {
  #ifdef DEBUG
  if(A1==A2)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A1->comp), DOUBLE_ALIGN);
    __assume_aligned(&(A2->comp), DOUBLE_ALIGN);
    __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k, l;
  double complex aux1[4] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex aux2[4] __attribute__((aligned(DOUBLE_ALIGN)));

  #define m2(X,Y) ((X)*2 + (Y))

  // reconstruct the complex form of the matrices
  aux1[m2(0,0)] = A1->comp[0]+I*A1->comp[3];
  aux1[m2(0,1)] = A1->comp[2]+I*A1->comp[1];
  aux1[m2(1,0)] =-A1->comp[2]+I*A1->comp[1];
  aux1[m2(1,1)] = A1->comp[0]-I*A1->comp[3];

  aux2[m2(0,0)] = A2->comp[0]+I*A2->comp[3];
  aux2[m2(0,1)] = A2->comp[2]+I*A2->comp[1];
  aux2[m2(1,0)] =-A2->comp[2]+I*A2->comp[1];
  aux2[m2(1,1)] = A2->comp[0]-I*A2->comp[3];

  for(i=0; i<2; i++)
     {
     for(j=0; j<2; j++)
        {
        for(k=0; k<2; k++)
           {
           for(l=0; l<2; l++)
              {
              TP->comp[i][j][k][l]=conj(aux1[m2(i,j)])*aux2[m2(k,l)];
              }
           }
        }
     }

  #undef m2
  }



// ***************** for Su2Vecs

// A=1
inline void one_Su2Vecs(Su2Vecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<4*NHIGGS; i++)
     {
     A->comp[i]=1.0;
     }
  }


// A=0
inline void zero_Su2Vecs(Su2Vecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<4*NHIGGS; i++)
     {
     A->comp[i]=0.0;
     }
  }


// A=B
inline void equal_Su2Vecs(Su2Vecs * restrict A, Su2Vecs const * const restrict B)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<4*NHIGGS; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A -> A^{dag}
inline void conjugate_Su2Vecs(Su2Vecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NHIGGS; i++)
     {
     A->comp[4*i+1]=-A->comp[4*i+1];
     A->comp[4*i+3]=-A->comp[4*i+3];
     }
  }


// A-=B
inline void minus_equal_Su2Vecs(Su2Vecs * restrict A, Su2Vecs const * const restrict B)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<4*NHIGGS; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A+=B
inline void plus_equal_Su2Vecs(Su2Vecs * restrict A, Su2Vecs const * const restrict B)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<4*NHIGGS; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// *= with real number
inline void times_equal_real_Su2Vecs(Su2Vecs * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<4*NHIGGS; i++)
     {
     A->comp[i]*=r;
     }
  }


// *= with real for a single component
inline void times_equal_real_single_Su2Vecs(Su2Vecs * restrict A, double r, int j)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[4*j+0]*=r;
  A->comp[4*j+1]*=r;
  A->comp[4*j+2]*=r;
  A->comp[4*j+3]*=r;
  }


// *= with complex number for a single component
inline void times_equal_complex_single_Su2Vecs(Su2Vecs * restrict A, double complex r, int j)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double complex a, b;

  a=A->comp[4*j+2]+I*A->comp[4*j+1];
  b=A->comp[4*j+0]-I*A->comp[4*j+3];

  a*=r;
  b*=r;

  A->comp[4*j+2]=creal(a);
  A->comp[4*j+1]=cimag(a);
  A->comp[4*j+0]=creal(b);
  A->comp[4*j+3]=-cimag(b);
  }


// norm
inline double norm_Su2Vecs(Su2Vecs const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris=0.0;

  for(i=0; i<4*NHIGGS; i++)
     {
     ris+=fabs(A->comp[i])*fabs(A->comp[i]);
     }

  return sqrt(ris);
  }


// normalize
inline void normalize_Su2Vecs(Su2Vecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double norm=norm_Su2Vecs(A);

  times_equal_real_Su2Vecs(A, 1.0/norm);
  }


// random vector (normalized)
void rand_vec_Su2Vecs(Su2Vecs * restrict A);


// real part of the scalar product re(v_1^{\dag}v_2)
inline double re_scal_prod_Su2Vecs(Su2Vecs const * const restrict v1, Su2Vecs const * const restrict v2)
  {
   #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int i, j;
  Su2 aux1, aux2 __attribute__((aligned(DOUBLE_ALIGN)));
  double ris=0.0;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<4; j++)
        {
        aux1.comp[j]=v1->comp[4*i+j];
        aux2.comp[j]=v2->comp[4*i+j];
        }

     times_equal_dag_Su2(&aux2, &aux1);
     ris+=retr_Su2(&aux2);  // the 0.5 that is due to the different normalization of matrices and vectors
                            // (see text above the definition of Su2Vecs) is not needed since retr = trace/2
     }

  return ris;
  }


// real part of the scalar product re(v_1[a]^{\dag}v_2[b]) with a, b flavour indices
inline double re_scal_prod_single_Su2Vecs(Su2Vecs const * const restrict v1, Su2Vecs const * const restrict v2, int a, int b)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int j;
  Su2 aux1, aux2 __attribute__((aligned(DOUBLE_ALIGN)));
  double ris;

  for(j=0; j<4; j++)
     {
     aux1.comp[j]=v1->comp[4*a+j];
     aux2.comp[j]=v2->comp[4*b+j];
     }

  times_equal_dag_Su2(&aux2, &aux1);
  ris=retr_Su2(&aux2);  // the 0.5 that is due to the different normalization of matrices and vectors
                        // (see text above the definition of Su2Vecs) is not needed since retr = trace/2

  return ris;
  }


// the i-th component of v2 is multiplied by "matrix"
// v1=matrix*v2
inline void matrix_times_vector_single_Su2Vecs(Su2Vecs * restrict v1, Su2 const * const restrict matrix, Su2Vecs const * const restrict v2, int i)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(matrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int j;
  Su2 aux1, aux2 __attribute__((aligned(DOUBLE_ALIGN)));

  equal_Su2Vecs(v1, v2);

  for(j=0; j<4; j++)
     {
     aux2.comp[j]=v2->comp[4*i+j];
     }

  times_Su2(&aux1, matrix, &aux2);

  for(j=0; j<4; j++)
     {
     v1->comp[4*i+j]=aux1.comp[j];
     }
  }


// all the components of v2 are multiplied by "matrix"
// v1=matrix*v2
inline void matrix_times_vector_all_Su2Vecs(Su2Vecs * restrict v1, Su2 const * const restrict matrix, Su2Vecs const * const restrict v2)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(matrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int i, j;
  Su2 aux1, aux2 __attribute__((aligned(DOUBLE_ALIGN)));

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<4; j++)
        {
        aux2.comp[j]=v2->comp[4*i+j];
        }

     times_Su2(&aux1, matrix, &aux2);

     for(j=0; j<4; j++)
        {
        v1->comp[4*i+j]=aux1.comp[j];
        }
     }
  }


// rotate two components of the vector
inline void rotate_two_components_Su2Vecs(Su2Vecs * restrict v1,
                                          Su2Vecs const * const restrict v2,
                                          int i,
                                          int j,
                                          double angle)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int k;

  equal_Su2Vecs(v1, v2);

  for(k=0; k<4; k++)
     {
     v1->comp[4*i+k]= cos(angle)*v2->comp[4*i+k] + sin(angle)*v2->comp[4*j+k];
     v1->comp[4*j+k]=-sin(angle)*v2->comp[4*i+k] + cos(angle)*v2->comp[4*j+k];
     }
  }


// tensor product of two vectors
// Re(v1^{\dag} * aux * v2) = ReTr(aux * matrix)
inline void vector_tensor_vector_Su2Vecs(Su2 * restrict matrix, Su2Vecs const * const restrict v1, Su2Vecs const * const restrict v2)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(matrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int i, j;
  Su2 aux1, aux2, aux3 __attribute__((aligned(DOUBLE_ALIGN)));

  zero_Su2(matrix);

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<4; j++)
        {
        aux1.comp[j]=v1->comp[4*i+j];
        aux2.comp[j]=v2->comp[4*i+j];
        }
     times_dag2_Su2(&aux3, &aux2, &aux1);

     plus_equal_Su2(matrix, &aux3);
     }

  times_equal_real_Su2(matrix, 0.5);   // 0.5 is due to the different normalization of matrices and vectors
                                       // see text above the definition of Su2Vecs
  }


// initialize the flavour matrix with a vector
// FM[mf(i,j)]=\sum_{on_gauge}conj(v1[i])v1[j] - delta^{ij}/N
// i, j are the flavour indices
inline void init_FMatrix_Su2Vecs(FMatrix * restrict fmatrix, Su2Vecs const * const restrict v1)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(fmatrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  zero_FMatrix(fmatrix);

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NHIGGS; j++)
        {
        fmatrix->comp[mf(i,j)]+=(v1->comp[4*i+0])*(v1->comp[4*j+0]);
        fmatrix->comp[mf(i,j)]+=(v1->comp[4*i+1])*(v1->comp[4*j+1]);
        fmatrix->comp[mf(i,j)]+=(v1->comp[4*i+2])*(v1->comp[4*j+2]);
        fmatrix->comp[mf(i,j)]+=(v1->comp[4*i+3])*(v1->comp[4*j+3]);

        fmatrix->comp[mf(i,j)]+=(v1->comp[4*i+2])*(v1->comp[4*j+1])*I;
        fmatrix->comp[mf(i,j)]-=(v1->comp[4*i+1])*(v1->comp[4*j+2])*I;

        fmatrix->comp[mf(i,j)]+=(v1->comp[4*i+3])*(v1->comp[4*j+0])*I;
        fmatrix->comp[mf(i,j)]-=(v1->comp[4*i+0])*(v1->comp[4*j+3])*I;
        }
     }

  for(i=0; i<NHIGGS; i++)
     {
     fmatrix->comp[mf(i,i)] -= 1.0/(double)NHIGGS;
     }
  }


// return a double coumplex number to check the fate of U(1) flavour symmetry
inline double complex HiggsU1Obs_Su2Vecs(Su2Vecs const * const restrict v1)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  #endif

  #if NHIGGS >=2
    double complex a, b, c, d;

    a = v1->comp[2] + (v1->comp[1])*I;
    b = v1->comp[0] - (v1->comp[3])*I;  // first flavour = (a,b)

    c = v1->comp[4+2] + (v1->comp[4+1])*I;
    d = v1->comp[4+0] - (v1->comp[4+3])*I;  // second flavour = (c, d)

    return a*d-b*c;   // = det[ (a,b), (c,d) ]
  #else
    (void) v1;
    return 0.0 + 0.0*I;
  #endif
  }

// print on file
int print_on_file_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A);


// print on binary file without changing endiannes
int print_on_binary_file_noswap_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A);


// print on binary file changing endiannes
int print_on_binary_file_swap_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A);


// print on binary file in big endian format
int print_on_binary_file_bigen_Su2Vecs(FILE *fp, Su2Vecs const * const restrict A);


// read from file
int read_from_file_Su2Vecs(FILE *fp, Su2Vecs * restrict A);


// read from binary file without changing endiannes
int read_from_binary_file_noswap_Su2Vecs(FILE *fp, Su2Vecs * restrict A);


// read from binary file changing endiannes
int read_from_binary_file_swap_Su2Vecs(FILE *fp, Su2Vecs * restrict A);


// read from binary file written in big endian
int read_from_binary_file_bigen_Su2Vecs(FILE *fp, Su2Vecs * restrict A);



#endif
