#ifndef SUN_MONOPOLES_C
#define SUN_MONOPOLES_C

#include"../include/macro.h"

#include<complex.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/sun_monopoles.h"
#include"../include/su2_monopoles.h"
#include"../include/sun.h"
#include"../include/sun_aux.h"


// This function compute the MAG gauge transformation in the SU(N) case
void comp_MAG_gauge_transformation_SuN (SuN X_links[2*STDIM],
                                        double lambda[NCOLOR],
                                        double OverRelaxParam,
                                        SuN *G_mag)

   {
   int i, j, dir;
   double vec_x[3];
   double X_mod;
   SuN X, aux1, aux2, aux3, aux_g;
   Su2 G_mag_su2;   
   
   zero_SuN(&X);
   one_SuN(G_mag);


   for(dir=0; dir<STDIM; dir++) 
      {
      //calcolo lambda*Udag_mu(n) e lambda*U_mu(n- mu)
      diag_matrix_times_dag_SuN(&aux1, lambda, &(X_links[dir]));
      diag_matrix_times_SuN(&aux2, lambda, &(X_links[dir+STDIM]));
      
      //Calcolo la matrice X(n)
      times_SuN(&aux3, &(X_links[dir]), &aux1);
      plus_equal_SuN(&X, &aux3);
      times_dag1_SuN(&aux3, &(X_links[dir+STDIM]), &aux2);
      plus_equal_SuN(&X, &aux3);
      }   
   
     #ifdef DEBUG
     double herm_real_check=0, herm_imag_check=0;
     complex trace_check = 0.0 + 0.0*I;
     
     for(i=0;i<NCOLOR-1; i++)
       {
       trace_check += creal(X.comp[m(i,i)])  + cimag(X.comp[m(i,i)])*I;
       for(j=i; j<NCOLOR; j++)
         { 
         herm_real_check += creal(X.comp[m(i,j)]) - creal(X.comp[m(j,i)]); 
         herm_imag_check += cimag(X.comp[m(i,j)]) + cimag(X.comp[m(j,i)]);
         } 
       }    
 
     if(herm_real_check > MIN_VALUE*10 || herm_imag_check > MIN_VALUE*10 || creal(trace_check) > MIN_VALUE*1000 || cimag(trace_check) > MIN_VALUE*10)
      {
      fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g trace_check_real = %.12g trace_check_imag = %.12g\n", herm_real_check, herm_imag_check, creal(trace_check)-24, cimag(trace_check));
      }
     #endif
     
 
      // Cycle on all the SU(2) subgroups
      // REF: M.D'Elia and C. Bonati Nuc. PHys B 877 (2012) 233-259
      for(i=0;i<NCOLOR-1;i++)
        {
        for(j=i+1;j<NCOLOR;j++)
          {
          vec_x[0] = creal(X.comp[m(i,j)]);
          vec_x[1] = cimag(X.comp[m(j,i)]);
          vec_x[2] = 0.5*(creal(X.comp[m(i,i)]) - creal(X.comp[m(j,j)]));
            
          //printf("SOTTOGRUPPO %d %d : %.12g %.12g %.12g\n", i, j, vec_x[0], vec_x[1], vec_x[2]);
 
          X_mod = sqrt(vec_x[0]*vec_x[0] + vec_x[1]*vec_x[1] + vec_x[2]*vec_x[2]);

          if(X_mod<MIN_VALUE)
           {
           }
          else
           { 
           diagonalize_X_Su2_aux(OverRelaxParam, vec_x, &G_mag_su2);
           duetoenne(&G_mag_su2, i, j, &aux_g);
         
           //print_on_screen(&X);
           // Update the X(n) matrix: X->GXGdag
           times_dag2_SuN(&aux3, &X, &aux_g);
           times_SuN(&aux2, &aux_g, &aux3);
           equal_SuN(&X, &aux2);
           // printf("%d %d\n", i, j);
          // print_on_screen(&X);
          
          
           // Build the gauge transformation matrix
           times_SuN(&aux2, &aux_g, G_mag);
           equal_SuN(G_mag, &aux2);

          //printf("SOTTOGRUPPO %d %d\n", i, j);
          //print_on_screen(G_mag);
           }
         }
        }
 
     #ifdef DEBUG
     trace_check = 0.0 + 0.0*I;
     herm_real_check=0;
     herm_imag_check=0; 
     for(i=0;i<NCOLOR; i++)
       {
       trace_check += creal(X.comp[m(i,i)])  + cimag(X.comp[m(i,i)])*I;
       for(j=i; j<NCOLOR; j++)
         {
         herm_real_check += creal(X.comp[m(i,j)]) - creal(X.comp[m(j,i)]); 
         herm_imag_check += cimag(X.comp[m(i,j)]) + cimag(X.comp[m(j,i)]);
         } 
       }
     if(herm_real_check > MIN_VALUE*10 || herm_imag_check > MIN_VALUE*10 || creal(trace_check) > MIN_VALUE*1000 || cimag(trace_check) > MIN_VALUE*10)
       {
       fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g trace_check_real = %.12g trace_check_imag = %.12g\n", herm_real_check, herm_imag_check, creal(trace_check), cimag(trace_check));
       }
     #endif     
}

//questa funzione calcola la media del modulo quadro dei contributi fuori diagonale dell'operatre X(n)
void comp_outdiagnorm_of_X_SuN (SuN X_links[2*STDIM],
                                         double *lambda,
                                         double *counter)
   {
   int i, j, dir;
   double aux=0;
   SuN X, aux1, aux2, aux3;
   
   zero_SuN(&X);


   for(dir=0; dir<STDIM; dir++) 
      {
      //calcolo lambda*Udag_mu(n) e lambda*U_mu(n- mu)
      diag_matrix_times_dag_SuN(&aux1, lambda, &(X_links[dir]));
      diag_matrix_times_SuN(&aux2, lambda, &(X_links[dir+STDIM]));
      
      //Calcolo la matrice X(n)
      times_SuN(&aux3, &(X_links[dir]), &aux1);
      plus_equal_SuN(&X, &aux3);
      times_dag1_SuN(&aux3, &(X_links[dir+STDIM]), &aux2);
      plus_equal_SuN(&X, &aux3);
      }   
  
     #ifdef DEBUG
     double herm_real_check=0, herm_imag_check=0;
     complex trace_check = 0.0 + 0.0*I;
     
     for(i=0;i<NCOLOR-1; i++)
       {
       trace_check += creal(X.comp[m(i,i)])  + cimag(X.comp[m(i,i)])*I;
       for(j=i+1; j<NCOLOR; j++)
         {
         herm_real_check += creal(X.comp[m(i,j)]) - creal(X.comp[m(j,i)]); 
         herm_imag_check += cimag(X.comp[m(i,j)]) + cimag(X.comp[m(j,i)]);
         } 
       }    
 
     if(herm_real_check > MIN_VALUE*10 || herm_imag_check > MIN_VALUE*10 || creal(trace_check) > MIN_VALUE*1000 || cimag(trace_check) > MIN_VALUE*10)
      {
      fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g trace_check_real = %.12g trace_check_imag = %.12g\n", herm_real_check, herm_imag_check, creal(trace_check), cimag(trace_check));
      }
     #endif 
      
      //compute the elements outside the diagonal
      for(i=0;i<NCOLOR-1;i++)
        {
        for(j=i+1;j<NCOLOR;j++)
           {
           aux += creal(X.comp[m(i,j)])*creal(X.comp[m(i,j)]) + cimag(X.comp[m(i,j)])*cimag(X.comp[m(i,j)]);
           }
        }

   *counter = aux;
   }



void comp_functional_fmag_SuN(SuN X_links[2*STDIM],
                              double lambda[NCOLOR],
                              double *fmag)
   {
   int dir, i;
   SuN F_mag, aux1, aux2, aux3;
   double fmag_aux=0;


   zero_SuN(&F_mag);
    
   for(dir=0;dir<STDIM;dir++)
     { 
     diag_matrix_times_dag_SuN(&aux1, lambda, &X_links[dir]);
     diag_matrix_times_SuN(&aux2, lambda, &X_links[dir]);
     times_SuN(&aux3, &aux1, &aux2);
     plus_equal_SuN(&F_mag, &aux3);
     }
   
   for(i=0;i<NCOLOR;i++)
     {
     fmag_aux += creal(F_mag.comp[m(i,i)]) + cimag(F_mag.comp[m(i,i)]);
     }
   *fmag=fmag_aux;
   }


void diag_projection_single_site_SuN(Gauge_Conf *GC,
                                     SuN *link, 
                                     long r,
                                     int dir)
   {
   int subg;
   double phi[NCOLOR], inv_rho[NCOLOR], dphi, inv_rho_sum, phi_aux;
   complex double u_ii; 


   dphi=0;
   inv_rho_sum=0;
   for(subg=0;subg<NCOLOR;subg++)
     {
     u_ii = link->comp[m(subg,subg)];
     phi[subg] = atan2(cimag(u_ii), creal(u_ii));

   //printf("INIZIO site %ld dir %d sub %d : %.12lf\n", r, dir, subg, phi[subg]);
     

     //D = diag(exp(i phi1),exp(i phi2),exp(i phi3), ... . exp(i phiN))
     //condition for maximum Re(Tr(U*D^dag)) is
     //|U_11| sin(phi1 - theta1) = |U_22| sin(phi2 - theta2) = |U_22| sin(phi3 - theta3) = ...
     //assuming phases are already almost factorized, we linearize "sin", then
     //(phii - thetai) \propto 1/|U_ii| 
     // the general formula is phidiag_i = phi_i - dphi *(|U_ii|^{-1}/sum_j|U_jj|^{-1})       


     inv_rho[subg] = 1.0/(sqrt(creal(u_ii)*creal(u_ii) + cimag(u_ii)*cimag(u_ii)));
     inv_rho_sum += inv_rho[subg];
     
     dphi += phi[subg]; 
     }
     if(dphi > PI)
      {
      dphi = dphi - PI2;
      }
     else if(dphi < -PI)
      {
      dphi = dphi + PI2;
      }
     

     for(subg=0;subg<NCOLOR;subg++)
       {      
       phi_aux =  phi[subg] - dphi*inv_rho[subg]/inv_rho_sum;
       if(phi_aux > PI)
        {
        phi_aux = phi_aux - PI2;
        }
       else if (phi_aux < -PI)
        {
        phi_aux = phi_aux + PI2;
        } 
       
       (GC->diag_proj[r][dir][subg]) = phi_aux;
       //printf("final angles site %ld dir %d sub %d : %.12lf\n", r, dir, subg, GC->diag_proj[r][dir][subg]);
       }
   }









#endif




