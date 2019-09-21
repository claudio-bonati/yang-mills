#ifndef SU2_MONOPOLES_C
#define SU2_MONOPOLES_C

#include"../include/macro.h"

#include<complex.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/su2_monopoles.h"
#include"../include/sun_monopoles.h"


// This function compute the MAG gauge transformation
void comp_MAG_gauge_transformation_Su2 (Su2 X_links[2*STDIM],
                                        double lambda[NCOLOR],
                                        double OverRelaxParam,
                                        Su2 *G_mag)

   {
   // The general 2x2 matrix A is written in the form
   // 
   //             a[0]   a[1]
   //             a[2]   a[3]
   //

   int dir, i;
   Su2 Ur, Umr;
   complex double X[4], Ur_matrix[4], Umr_matrix[4], aux1[4], aux2[4];
   double vec_x[3];  
 
   for(i=0; i<4; i++)
     {
     X[i] = 0.0 + 0.0*I;
     }

   for(dir=0; dir<STDIM; dir++)
     {
     for(i=0; i<4; i++)
       {
       aux1[i] = 0.0 + 0.0*I;
       aux2[i] = 0.0 + 0.0*I;
       }
     // I need to change the format of the Su2 variable
     equal_Su2(&Umr, &(X_links[dir+STDIM]));
     equal_Su2(&Ur, &(X_links[dir]));
 
     Ur_matrix[0]  =  Ur.comp[0] + Ur.comp[3]*I;     
     Ur_matrix[1]  =  Ur.comp[2] + Ur.comp[1]*I;     
     Ur_matrix[2]  = -Ur.comp[2] + Ur.comp[1]*I;     
     Ur_matrix[3]  =  Ur.comp[0] - Ur.comp[3]*I;     
  
     Umr_matrix[0]  =  Umr.comp[0] + Umr.comp[3]*I;     
     Umr_matrix[1]  =  Umr.comp[2] + Umr.comp[1]*I;     
     Umr_matrix[2]  = -Umr.comp[2] + Umr.comp[1]*I;     
     Umr_matrix[3]  =  Umr.comp[0] - Umr.comp[3]*I;     
    
     // I compute aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n-dir)
     aux1[0] = lambda[0]*conj(Ur_matrix[0]);
     aux1[1] = lambda[0]*conj(Ur_matrix[2]);
     aux1[2] = lambda[1]*conj(Ur_matrix[1]);
     aux1[3] = lambda[1]*conj(Ur_matrix[3]);
     
     aux2[0] = lambda[0]*Umr_matrix[0];
     aux2[1] = lambda[0]*Umr_matrix[1];
     aux2[2] = lambda[1]*Umr_matrix[2];
     aux2[3] = lambda[1]*Umr_matrix[3];
     
     // I compute  X = U_{dir}(n)*lambda*U^{dag}_{dir}(n) and U^{dag}_{dir}(n-dir)*lambda*U_{dir}(n-dir)
    
     X[0] += Ur_matrix[0]*aux1[0] + Ur_matrix[1]*aux1[2]; 
     X[1] += Ur_matrix[0]*aux1[1] + Ur_matrix[1]*aux1[3]; 
     X[2] += Ur_matrix[2]*aux1[0] + Ur_matrix[3]*aux1[2]; 
     X[3] += Ur_matrix[2]*aux1[1] + Ur_matrix[3]*aux1[3];

     X[0] += conj(Umr_matrix[0])*aux2[0] + conj(Umr_matrix[2])*aux2[2]; 
     X[1] += conj(Umr_matrix[0])*aux2[1] + conj(Umr_matrix[2])*aux2[3]; 
     X[2] += conj(Umr_matrix[1])*aux2[0] + conj(Umr_matrix[3])*aux2[2]; 
     X[3] += conj(Umr_matrix[1])*aux2[1] + conj(Umr_matrix[3])*aux2[3]; 
     }
  
     #ifdef DEBUG
     double herm_real_check=0, herm_imag_check=0;
     complex trace_check = 0.0 + 0.0*I;

     trace_check = creal(X[0]) + creal(X[3]) + (cimag(X[0]) + cimag(X[3]))*I;
     herm_real_check = creal(X[1]) - creal(X[2]); 
     herm_imag_check = cimag(X[1]) + cimag(X[2]); 
    
     if(herm_real_check > MIN_VALUE*10 || herm_imag_check > MIN_VALUE*10 || creal(trace_check) > MIN_VALUE*10 || cimag(trace_check) > MIN_VALUE*10)
      {
      fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g trace_check_real = %.12g trace_check_imag = %.12g\n", herm_real_check, herm_imag_check, creal(trace_check), cimag(trace_check));
      }
     #endif   
 
     // Convert the matrix X in the form X = \vec{x}*\vec{sigma} 
     vec_x[0] = creal(X[1]);
     vec_x[1] = cimag(X[2]);
     vec_x[2] = creal(X[0]);
    
 
     max_X_comp_G_Su2_aux(OverRelaxParam, vec_x, G_mag);
   }







// ref C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 
void max_X_comp_G_Su2_aux (double OverRelaxParam,
		           double X[3],
                           Su2 *G)
   {
   double alpha_max=0, duealpha_max=0, X_mod=0; 
   X_mod = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
   
   X[0] = X[0]/X_mod;
   X[1] = X[1]/X_mod;
   X[2] = X[2]/X_mod;
  

   X_mod = sqrt(X[0]*X[0] + X[1]*X[1]);
  
   X[0] = X[0]/X_mod; // sin(phi_max) 
   X[1] = X[1]/X_mod; // -cos(phi_max)
   
   //printf("sinphi = %.12g, -cosphi =%.12g ", X[0], X[1]);
 
   if(X_mod > MIN_VALUE)
    {
    if((1.0-X[2]) > 0.2)
     {
     duealpha_max = acos(X[2]);
     //printf("primaarcococose = %.12g ", X[2]);
     }
    else
      {
      duealpha_max = X_mod + pow(X_mod, 3)/6.0;
      }

    alpha_max = OverRelaxParam*(duealpha_max/2.0);

    //printf("%.12g\n", alpha_max);

    if(alpha_max < 0.01)
     {
     G->comp[0] =  1.0 - alpha_max*alpha_max/2.0 + pow(alpha_max, 4)/24.0;
     G->comp[1] = -(alpha_max - pow(alpha_max,3)/6.0)*X[1];
     G->comp[2] =  (alpha_max - pow(alpha_max,3)/6.0)*X[0];
     G->comp[3] = 0.0;
     }
    else
      {
      G->comp[0] =  cos(alpha_max);
      G->comp[1] = -sin(alpha_max)*X[1];
      G->comp[2] =  sin(alpha_max)*X[0];
      G->comp[3] = 0.0;
      }
    //printf("COMPONENTS OF G MATRIX = %.12g %.12g %.12g %.12g\n", G->comp[0], G->comp[1], G->comp[2], G->comp[3]);
    }
 }


 
void comp_non_diagonal_contribution_Su2 (Su2 X_links[2*STDIM], 
                                         double lambda[2],
                                         double *non_diag_contr)
   {
   int dir, i;
   Su2 Ur, Umr;
   complex double X[4], Ur_matrix[4], Umr_matrix[4], aux1[4], aux2[4];
 
   for(i=0; i<4; i++)
     {
     X[i] = 0.0 + 0.0*I;
     }

   for(dir=0; dir<STDIM; dir++)
     {
     for(i=0; i<4; i++)
       {
       aux1[i] = 0.0 + 0.0*I;
       aux2[i] = 0.0 + 0.0*I;
       }
     // I need to change the format of the Su2 variable
     equal_Su2(&Umr, &(X_links[dir+STDIM]));
     equal_Su2(&Ur, &(X_links[dir]));
 
     Ur_matrix[0]  =  Ur.comp[0] + Ur.comp[3]*I;     
     Ur_matrix[1]  =  Ur.comp[2] + Ur.comp[1]*I;     
     Ur_matrix[2]  = -Ur.comp[2] + Ur.comp[1]*I;     
     Ur_matrix[3]  =  Ur.comp[0] - Ur.comp[3]*I;     
  
     Umr_matrix[0]  =  Umr.comp[0] + Umr.comp[3]*I;     
     Umr_matrix[1]  =  Umr.comp[2] + Umr.comp[1]*I;     
     Umr_matrix[2]  = -Umr.comp[2] + Umr.comp[1]*I;     
     Umr_matrix[3]  =  Umr.comp[0] - Umr.comp[3]*I;     
    
     // I compute aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n-dir)
     aux1[0] = lambda[0]*conj(Ur_matrix[0]);
     aux1[1] = lambda[0]*conj(Ur_matrix[2]);
     aux1[2] = lambda[1]*conj(Ur_matrix[1]);
     aux1[3] = lambda[1]*conj(Ur_matrix[3]);
     
     aux2[0] = lambda[0]*Umr_matrix[0];
     aux2[1] = lambda[0]*Umr_matrix[1];
     aux2[2] = lambda[1]*Umr_matrix[2];
     aux2[3] = lambda[1]*Umr_matrix[3];
     
     // I compute  X = U_{dir}(n)*lambda*U^{dag}_{dir}(n) and U^{dag}_{dir}(n-dir)*lambda*U_{dir}(n-dir)
    
     X[0] += Ur_matrix[0]*aux1[0] + Ur_matrix[1]*aux1[2]; 
     X[1] += Ur_matrix[0]*aux1[1] + Ur_matrix[1]*aux1[3]; 
     X[2] += Ur_matrix[2]*aux1[0] + Ur_matrix[3]*aux1[2]; 
     X[3] += Ur_matrix[2]*aux1[1] + Ur_matrix[3]*aux1[3];

     X[0] += conj(Umr_matrix[0])*aux2[0] + conj(Umr_matrix[2])*aux2[2]; 
     X[1] += conj(Umr_matrix[0])*aux2[1] + conj(Umr_matrix[2])*aux2[3]; 
     X[2] += conj(Umr_matrix[1])*aux2[0] + conj(Umr_matrix[3])*aux2[2]; 
     X[3] += conj(Umr_matrix[1])*aux2[1] + conj(Umr_matrix[3])*aux2[3]; 
     }
  
     #ifdef DEBUG
     double herm_real_check=0, herm_imag_check=0;
     complex trace_check = 0.0 + 0.0*I;

     trace_check = creal(X[0]) + creal(X[3]) + (cimag(X[0]) + cimag(X[3]))*I;
     herm_real_check = creal(X[1]) - creal(X[2]); 
     herm_imag_check = cimag(X[1]) + cimag(X[2]); 
    
     if(herm_real_check > MIN_VALUE*10 || herm_imag_check > MIN_VALUE*10 || creal(trace_check) > MIN_VALUE*10 || cimag(trace_check) > MIN_VALUE*10)
      {
      fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g trace_check_real = %.12g trace_check_imag = %.12g\n", herm_real_check, herm_imag_check, creal(trace_check), cimag(trace_check));
      }
     #endif

     *non_diag_contr =  creal(X[1])*creal(X[1]) + cimag(X[1])*cimag(X[1])  + creal(X[2])*creal(X[2]) + cimag(X[2])*cimag(X[2]);  
   }


void comp_functional_fmag_Su2(Su2 X_links[2*STDIM], 
                              double lambda[2],
                              double *fmag)
   {
   int dir, i;
   Su2 U, Udag;
   complex double X[4], U_matrix[4], Udag_matrix[4], aux1[4], aux2[4];  
 
   for(i=0; i<4; i++)
     {
     X[i] = 0.0 + 0.0*I;
     }

   for(dir=0; dir<STDIM; dir++)
     {
     for(i=0; i<4; i++)
       {
       aux1[i] = 0.0 + 0.0*I;
       aux2[i] = 0.0 + 0.0*I;
       }

     // I need to change the format of the Su2 variable
     equal_Su2(&U, &(X_links[dir]));
     equal_dag_Su2(&Udag, &(X_links[dir]));
 
     U_matrix[0]  =  U.comp[0] + U.comp[3]*I;     
     U_matrix[1]  =  U.comp[2] + U.comp[1]*I;     
     U_matrix[2]  = -U.comp[2] + U.comp[1]*I;     
     U_matrix[3]  =  U.comp[0] - U.comp[3]*I;     
  
     Udag_matrix[0]  =  Udag.comp[0] + Udag.comp[3]*I;     
     Udag_matrix[1]  =  Udag.comp[2] + Udag.comp[1]*I;     
     Udag_matrix[2]  = -Udag.comp[2] + Udag.comp[1]*I;     
     Udag_matrix[3]  =  Udag.comp[0] - Udag.comp[3]*I;     
    
     // I compute aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n-dir)
     aux1[0] = lambda[0]*Udag_matrix[0];
     aux1[1] = lambda[0]*Udag_matrix[2];
     aux1[2] = lambda[1]*Udag_matrix[1];
     aux1[3] = lambda[1]*Udag_matrix[3];
     
     aux2[0] = lambda[0]*U_matrix[0];
     aux2[1] = lambda[0]*U_matrix[1];
     aux2[2] = lambda[1]*U_matrix[2];
     aux2[3] = lambda[1]*U_matrix[3];
     
     // I compute  X = U_{dir}(n)*lambda*U^{dag}_{dir}(n) and U^{dag}_{dir}(n-dir)*lambda*U_{dir}(n-dir)
    
     X[0] += aux1[0]*aux2[0]; 
     X[1] += aux1[0]*aux2[1]; 
     X[2] += aux1[2]*aux2[0]; 
     X[3] += aux1[2]*aux2[1]; 
     }
   
   *fmag = creal(X[0]) + cimag(X[0]) + creal(X[3]) + cimag(X[3]);
   }
















#endif 
