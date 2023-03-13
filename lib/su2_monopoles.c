#ifndef SU2_MONOPOLES_C
#define SU2_MONOPOLES_C

#include"../include/macro.h"

#include<complex.h>
#include<stdlib.h>
#include<stdio.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/su2_monopoles.h"
#include"../include/sun_monopoles.h"


// This function compute the MAG gauge transformation
void comp_MAG_gauge_transformation_Su2(Su2 X_links[2*STDIM],
                                       double const lambda[2],
                                       double overrelaxparam,
                                       Su2 *G_mag)

   {
   int dir, i;
   Su2 Ur, Urm;
   double complex X[4], Ur_matrix[4], Urm_matrix[4], aux1[4], aux2[4];
   double vec_x[3];

   // the general 2x2 matrix A is written in the form
   // a[0]  a[1]
   // a[2]  a[3]

   for(i=0; i<4; i++)
      {
      X[i] = 0.0 + 0.0*I;
      }

   for(dir=0; dir<STDIM; dir++)
      {
      equal_Su2(&Ur, &(X_links[dir]));          // U_dir(r)
      equal_Su2(&Urm, &(X_links[dir+STDIM]));   // U_dir(r-dir)

      // the matrices are written as 2x2 complex matrices (no Pauli)
      // since the result of the product with lambda cannot be represented in the Pauli form

      Ur_matrix[0]  =  Ur.comp[0] + Ur.comp[3]*I;
      Ur_matrix[1]  =  Ur.comp[2] + Ur.comp[1]*I;
      Ur_matrix[2]  = -Ur.comp[2] + Ur.comp[1]*I;
      Ur_matrix[3]  =  Ur.comp[0] - Ur.comp[3]*I;

      Urm_matrix[0]  =  Urm.comp[0] + Urm.comp[3]*I;
      Urm_matrix[1]  =  Urm.comp[2] + Urm.comp[1]*I;
      Urm_matrix[2]  = -Urm.comp[2] + Urm.comp[1]*I;
      Urm_matrix[3]  =  Urm.comp[0] - Urm.comp[3]*I;
    
      // aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n-dir)
      aux1[0] = lambda[0]*conj(Ur_matrix[0]);
      aux1[1] = lambda[0]*conj(Ur_matrix[2]);
      aux1[2] = lambda[1]*conj(Ur_matrix[1]);
      aux1[3] = lambda[1]*conj(Ur_matrix[3]);

      aux2[0] = lambda[0]*Urm_matrix[0];
      aux2[1] = lambda[0]*Urm_matrix[1];
      aux2[2] = lambda[1]*Urm_matrix[2];
      aux2[3] = lambda[1]*Urm_matrix[3];

      // X = U_{dir}(n)*lambda*U^{dag}_{dir}(n) + U^{dag}_{dir}(n-dir)*lambda*U_{dir}(n-dir)
      // as in C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]

      X[0] += Ur_matrix[0]*aux1[0] + Ur_matrix[1]*aux1[2];
      X[1] += Ur_matrix[0]*aux1[1] + Ur_matrix[1]*aux1[3];
      X[2] += Ur_matrix[2]*aux1[0] + Ur_matrix[3]*aux1[2];
      X[3] += Ur_matrix[2]*aux1[1] + Ur_matrix[3]*aux1[3];

      X[0] += conj(Urm_matrix[0])*aux2[0] + conj(Urm_matrix[2])*aux2[2];
      X[1] += conj(Urm_matrix[0])*aux2[1] + conj(Urm_matrix[2])*aux2[3];
      X[2] += conj(Urm_matrix[1])*aux2[0] + conj(Urm_matrix[3])*aux2[2];
      X[3] += conj(Urm_matrix[1])*aux2[1] + conj(Urm_matrix[3])*aux2[3];
      }

   #ifdef DEBUG
   // check if X is traceless hermitian (as it should!)
   double herm_real_check=0, herm_imag_check=0;
   double complex trace_check = 0.0 + 0.0*I;

   trace_check = creal(X[0]) + creal(X[3]) + (cimag(X[0]) + cimag(X[3]))*I;
   herm_real_check = creal(X[1]) - creal(X[2]);
   herm_imag_check = cimag(X[1]) + cimag(X[2]);

   if(fabs(herm_real_check) > 10*MIN_VALUE || fabs(herm_imag_check) > 10*MIN_VALUE || cabs(trace_check) > 10*MIN_VALUE )
     {
     fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g ", herm_real_check, herm_imag_check);
     fprintf(stderr, "trace_check_real = %.12g trace_check_imag = %.12g ", creal(trace_check), cimag(trace_check));
     fprintf(stderr, "(%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif
 
   // write X as X = vec_x[0]*sigma1 + vec_x[1]*sigma2 + vec_x[2]*sigma3
   vec_x[0] = creal(X[1]);
   vec_x[1] = cimag(X[2]);
   vec_x[2] = creal(X[0]);

   // find the gauge transformation that makes X diagonal
   diagonalize_X_Su2_aux(overrelaxparam, vec_x, G_mag);
   }



// diagonalize the X matrix defined in
// C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]
// and returns the gauge transformation that implement the diagonalization
void diagonalize_X_Su2_aux(double overrelaxparam,
                           double X[3],
                           Su2 *G)
   {
   double alpha_max=0, twoalpha_max=0, X_mod=0;

   X_mod = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
   if(X_mod > MIN_VALUE)
     {
     X[0] = X[0]/X_mod;
     X[1] = X[1]/X_mod;
     X[2] = X[2]/X_mod;

     X_mod = sqrt(X[0]*X[0] + X[1]*X[1]);

     if(X_mod > MIN_VALUE)
       {
       X[0] = X[0]/X_mod; // =sin(phi_max) in the paper
       X[1] = X[1]/X_mod; // -cos(phi_max) in the paper

       twoalpha_max = acos(X[2]);

       alpha_max = overrelaxparam*(twoalpha_max/2.0);

       G->comp[0] =  cos(alpha_max);
       G->comp[1] = -sin(alpha_max)*X[1];
       G->comp[2] =  sin(alpha_max)*X[0];
       G->comp[3] = 0.0;
       }
     else
       {
       one_Su2(G);
       }
     }
   else
     {
     one_Su2(G);
     }
 }


// compute the square norm of the out-of-diagonal elements of X
void comp_outdiagnorm_of_X_Su2(Su2 X_links[2*STDIM],
                               double const lambda[2],
                               double *non_diag_contr)
   {
   int dir, i;
   Su2 Ur, Urm;
   complex double X[4], Ur_matrix[4], Urm_matrix[4], aux1[4], aux2[4];

   // the general 2x2 matrix A is written in the form
   // a[0]  a[1]
   // a[2]  a[3]

   for(i=0; i<4; i++)
      {
      X[i] = 0.0 + 0.0*I;
      }

   for(dir=0; dir<STDIM; dir++)
      {
      equal_Su2(&Ur, &(X_links[dir]));
      equal_Su2(&Urm, &(X_links[dir+STDIM]));

      // the matrices are written as 2x2 complex matrices (no Pauli)
      // since the result of the product with lambda cannot be represented in the Pauli form

      Ur_matrix[0]  =  Ur.comp[0] + Ur.comp[3]*I;
      Ur_matrix[1]  =  Ur.comp[2] + Ur.comp[1]*I;
      Ur_matrix[2]  = -Ur.comp[2] + Ur.comp[1]*I;
      Ur_matrix[3]  =  Ur.comp[0] - Ur.comp[3]*I;

      Urm_matrix[0]  =  Urm.comp[0] + Urm.comp[3]*I;
      Urm_matrix[1]  =  Urm.comp[2] + Urm.comp[1]*I;
      Urm_matrix[2]  = -Urm.comp[2] + Urm.comp[1]*I;
      Urm_matrix[3]  =  Urm.comp[0] - Urm.comp[3]*I;

      // aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n-dir)
      aux1[0] = lambda[0]*conj(Ur_matrix[0]);
      aux1[1] = lambda[0]*conj(Ur_matrix[2]);
      aux1[2] = lambda[1]*conj(Ur_matrix[1]);
      aux1[3] = lambda[1]*conj(Ur_matrix[3]);

      aux2[0] = lambda[0]*Urm_matrix[0];
      aux2[1] = lambda[0]*Urm_matrix[1];
      aux2[2] = lambda[1]*Urm_matrix[2];
      aux2[3] = lambda[1]*Urm_matrix[3];

      // X = U_{dir}(n)*lambda*U^{dag}_{dir}(n) + U^{dag}_{dir}(n-dir)*lambda*U_{dir}(n-dir)
      // as in C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]

      X[0] += Ur_matrix[0]*aux1[0] + Ur_matrix[1]*aux1[2];
      X[1] += Ur_matrix[0]*aux1[1] + Ur_matrix[1]*aux1[3];
      X[2] += Ur_matrix[2]*aux1[0] + Ur_matrix[3]*aux1[2];
      X[3] += Ur_matrix[2]*aux1[1] + Ur_matrix[3]*aux1[3];

      X[0] += conj(Urm_matrix[0])*aux2[0] + conj(Urm_matrix[2])*aux2[2];
      X[1] += conj(Urm_matrix[0])*aux2[1] + conj(Urm_matrix[2])*aux2[3];
      X[2] += conj(Urm_matrix[1])*aux2[0] + conj(Urm_matrix[3])*aux2[2];
      X[3] += conj(Urm_matrix[1])*aux2[1] + conj(Urm_matrix[3])*aux2[3];
      }

   #ifdef DEBUG
   // check if X is traceless hermitian (as it should!)
   double herm_real_check=0, herm_imag_check=0;
   double complex trace_check = 0.0 + 0.0*I;

   trace_check = creal(X[0]) + creal(X[3]) + (cimag(X[0]) + cimag(X[3]))*I;
   herm_real_check = creal(X[1]) - creal(X[2]);
   herm_imag_check = cimag(X[1]) + cimag(X[2]);

   if(herm_real_check > MIN_VALUE*10 || herm_imag_check > MIN_VALUE*10 || creal(trace_check) > MIN_VALUE*10 || cimag(trace_check) > MIN_VALUE*10)
     {
     fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g ", herm_real_check, herm_imag_check);
     fprintf(stderr, "trace_check_real = %.12g trace_check_imag = %.12g", creal(trace_check), cimag(trace_check));
     fprintf(stderr, "(%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

   *non_diag_contr =  creal(X[1])*creal(X[1]) + cimag(X[1])*cimag(X[1])  + creal(X[2])*creal(X[2]) + cimag(X[2])*cimag(X[2]);
   }


// compute the value functional to be maximized in MAG
// as in C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]
void comp_functional_fmag_Su2(Su2 X_links[2*STDIM], 
                              double const lambda[2],
                              double *fmag)
   {
   int dir;
   double sum;
   Su2 U, Udag;
   complex double U_matrix[4], Udag_matrix[4], aux1[4], aux2[4];

   // the general 2x2 matrix A is written in the form
   // a[0]  a[1]
   // a[2]  a[3]

   sum=0.0;

   for(dir=0; dir<STDIM; dir++)
      {
      equal_Su2(&U, &(X_links[dir]));
      equal_dag_Su2(&Udag, &(X_links[dir]));

      // the matrices are written as 2x2 complex matrices (no Pauli)
      // since the result of the product with lambda cannot be represented in the Pauli form

      U_matrix[0]  =  U.comp[0] + U.comp[3]*I;
      U_matrix[1]  =  U.comp[2] + U.comp[1]*I;
      U_matrix[2]  = -U.comp[2] + U.comp[1]*I;
      U_matrix[3]  =  U.comp[0] - U.comp[3]*I;
  
      Udag_matrix[0]  =  Udag.comp[0] + Udag.comp[3]*I;
      Udag_matrix[1]  =  Udag.comp[2] + Udag.comp[1]*I;
      Udag_matrix[2]  = -Udag.comp[2] + Udag.comp[1]*I;
      Udag_matrix[3]  =  Udag.comp[0] - Udag.comp[3]*I;
    
      // aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n)
      aux1[0] = lambda[0]*Udag_matrix[0];
      aux1[1] = lambda[0]*Udag_matrix[2];
      aux1[2] = lambda[1]*Udag_matrix[1];
      aux1[3] = lambda[1]*Udag_matrix[3];

      aux2[0] = lambda[0]*U_matrix[0];
      aux2[1] = lambda[0]*U_matrix[1];
      aux2[2] = lambda[1]*U_matrix[2];
      aux2[3] = lambda[1]*U_matrix[3];

      // sum += tr(aux1*aux2) (note that this trace is a real number, thus the creal)
      sum += creal( aux1[0]*aux2[0] + aux1[1]*aux2[2] + aux1[2]*aux2[1] + aux1[3]*aux1[3] );
      }
   
   *fmag = sum;
   }


// extract the diagonal part of the link and save it in GC->diag_proj
void diag_projection_single_site_Su2(Gauge_Conf *GC,
                                     Su2 *link, 
                                     long r,
                                     int dir)
   {
   double phi0, phi1;
   double complex U_matrix0, U_matrix3;
 
   U_matrix0  =  link->comp[0] + link->comp[3]*I;
   U_matrix3  =  link->comp[0] - link->comp[3]*I;
 
   phi0 = atan2(cimag(U_matrix0), creal(U_matrix0));
   phi1 = atan2(cimag(U_matrix3), creal(U_matrix3));

   GC->diag_proj[r][dir][0] = phi0;
   GC->diag_proj[r][dir][1] = phi1;
   }



#endif 
