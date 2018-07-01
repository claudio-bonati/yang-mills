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

void init_Su2(Su2 * restrict A, double vec[4])
  {
  A->comp[0]=vec[0];
  A->comp[1]=vec[1];
  A->comp[2]=vec[2];
  A->comp[3]=vec[3];
  }
  

// A=1
void one_Su2(Su2 * restrict A)
  {
  A->comp[0]=1.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=0
void zero_Su2(Su2 * restrict A)
  {
  A->comp[0]=0.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=B
void equal_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]=B->comp[0];
  A->comp[1]=B->comp[1];
  A->comp[2]=B->comp[2];
  A->comp[3]=B->comp[3];
  }


// A=B^{dag}
void equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]= B->comp[0];
  A->comp[1]=-B->comp[1];
  A->comp[2]=-B->comp[2];
  A->comp[3]=-B->comp[3];
  }


// A+=B
void plus_equal_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]+=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A+=B^{dag}
void plus_equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]+=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=B
void minus_equal_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]-=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=(r*B)
void minus_equal_times_real_Su2(Su2 * restrict A, Su2 const * restrict B, double r)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]-=(r*B->comp[0]);
  A->comp[1]-=(r*B->comp[1]);
  A->comp[2]-=(r*B->comp[2]);
  A->comp[3]-=(r*B->comp[3]);
  }


// A-=B^{dag}
void minus_equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]-=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A=b*B+c*C
void lin_comb_Su2(Su2 * restrict A,
                  double b, Su2 const * restrict B,
                  double c, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] + c*C->comp[1];
  A->comp[2]= b*B->comp[2] + c*C->comp[2];
  A->comp[3]= b*B->comp[3] + c*C->comp[3];
  }


// A=b*B^{dag}+c*C
void lin_comb_dag1_Su2(Su2 * restrict A,
                       double b, Su2 const * restrict B,
                       double c, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] + c*C->comp[1];
  A->comp[2]= -b*B->comp[2] + c*C->comp[2];
  A->comp[3]= -b*B->comp[3] + c*C->comp[3];
  }


// A=b*B+c*C^{dag}
void lin_comb_dag2_Su2(Su2 * restrict A,
                       double b, Su2 const * restrict B,
                       double c, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] - c*C->comp[1];
  A->comp[2]= b*B->comp[2] - c*C->comp[2];
  A->comp[3]= b*B->comp[3] - c*C->comp[3];
  }


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_Su2(Su2 * restrict A,
                        double b, Su2 const * restrict B,
                        double c, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] - c*C->comp[1];
  A->comp[2]= -b*B->comp[2] - c*C->comp[2];
  A->comp[3]= -b*B->comp[3] - c*C->comp[3];
  }


// A*=r
void times_equal_real_Su2(Su2 * restrict A, double r)
  {

  A->comp[0]*=r;
  A->comp[1]*=r;
  A->comp[2]*=r;
  A->comp[3]*=r;
  }


// A*=B
void times_equal_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
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
void times_equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
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
void times_Su2(Su2 * restrict A, Su2 const * restrict B, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]= B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// A=B^{dag}*C
void times_dag1_Su2(Su2 * restrict A, Su2 const * restrict B, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]= B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B*C^{dag}
void times_dag2_Su2(Su2 * restrict A, Su2 const * restrict B, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]=  B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B^{dag}*C^{dag}
void times_dag12_Su2(Su2 * restrict A, Su2 const * restrict B, Su2 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  A->comp[0]=  B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// random SU(2) matrix
void rand_matrix_Su2(Su2 * restrict A)
  {
  double p0, p1, p2, p3, p;

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
double sqrtdet_Su2(Su2 const * restrict A)
  {
  return sqrt(A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1]
             + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3]);
  }


// l2 norm of the matrix
double norm_Su2(Su2 const * restrict A)
  {
  return sqrtdet_Su2(A);
  }


// real part of the trace /2
double retr_Su2(Su2 const * restrict A)
  {
  return A->comp[0];
  }


// imaginary part of the trace /2
double imtr_Su2(Su2 const * restrict A)
  {
  (void) A; // to suppress compilation warning of unused variable
  return 0.0;
  }


// unitarize the matrix
void unitarize_Su2(Su2 * restrict A)
  {
  double p;

  p=A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1] + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3];
  p=1.0/sqrt(p);

  A->comp[0]*=p;
  A->comp[1]*=p;
  A->comp[2]*=p;
  A->comp[3]*=p;
  }


// traceless antihermitian part
void ta_Su2(Su2 * restrict A)
  {
  A->comp[0]=0;
  }


// exponential of the traceless antihermitian part
void taexp_Su2(Su2 * restrict A)
  {
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
void print_on_screen_Su2(Su2 const * const A)
  {
  printf("%.16lf %.16lf %.16lf %.16lf\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);
  }


// print on file
void print_on_file_Su2(FILE *fp, Su2 const * const A)
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
void print_on_binary_file_noswap_Su2(FILE *fp, Su2 const * const A)
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
void print_on_binary_file_swap_Su2(FILE *fp, Su2 const * const A)
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
void print_on_binary_file_bigen_Su2(FILE *fp, const Su2 * const A)
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
void read_from_file_Su2(FILE *fp, Su2 *A)
  {
  int err=fscanf(fp, "%lg %lg %lg %lg", &(A->comp[0]), &(A->comp[1]), &(A->comp[2]), &(A->comp[3]));

  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// read from binary file without changing endiannes
void read_from_binary_file_noswap_Su2(FILE *fp, Su2 *A)
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
void read_from_binary_file_swap_Su2(FILE *fp, Su2 *A)
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
void read_from_binary_file_bigen_Su2(FILE *fp, Su2 *A)
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


void TensProd_init_Su2(TensProd * restrict TP, Su2 const * restrict A1, Su2 const * restrict A2)
  {
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


#endif
