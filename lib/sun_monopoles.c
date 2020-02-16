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

// This function compute the MAG gauge transformation in the SU(N) case
void comp_MAG_gauge_transformation_SuN(SuN X_links[2*STDIM],
                                       double const lambda[NCOLOR],
                                       double overrelaxparam,
                                       SuN *G_mag)

   {
   int i, j, dir, k;
   double vec_x[3], X_mod;
   complex double gii, gij, gji, gjj, tmp0, tmp1;
   SuN X, aux1, aux2, aux3;
   Su2 G2;
   
   zero_SuN(&X);
   one_SuN(G_mag);

   for(dir=0; dir<STDIM; dir++) 
      {
      // aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n-dir)
      diag_matrix_times_dag_SuN(&aux1, lambda, &(X_links[dir]));
      diag_matrix_times_SuN(&aux2, lambda, &(X_links[dir+STDIM]));
      
      // X = U_{dir}(n)*lambda*U^{dag}_{dir}(n) + U^{dag}_{dir}(n-dir)*lambda*U_{dir}(n-dir)
      // as in C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]
      times_SuN(&aux3, &(X_links[dir]), &aux1);
      plus_equal_SuN(&X, &aux3);
      times_dag1_SuN(&aux3, &(X_links[dir+STDIM]), &aux2);
      plus_equal_SuN(&X, &aux3);
      }   
   
   #ifdef DEBUG
   double herm_real_check=0.0, herm_imag_check=0.0;
   double complex trace_check = 0.0+0.0*I;

   for(i=0; i<NCOLOR; i++)
      {
      trace_check += X.comp[m(i,i)];
      for(j=i+1; j<NCOLOR; j++)
         {
         herm_real_check += creal(X.comp[m(i,j)]) - creal(X.comp[m(j,i)]);
         herm_imag_check += cimag(X.comp[m(i,j)]) + cimag(X.comp[m(j,i)]);
         }
      }

   if(fabs(herm_real_check) > MIN_VALUE*10 || fabs(herm_imag_check) > MIN_VALUE*10 || cabs(trace_check) > MIN_VALUE*10 )
     {
     fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g ", herm_real_check, herm_imag_check);
     fprintf(stderr, "trace_check_real = %.12g trace_check_imag = %.12g ", creal(trace_check), cimag(trace_check));
     fprintf(stderr, "(%s, %d)\n", __FILE__, __LINE__);
     }
   #endif
 
   // cycle on all the SU(2) subgroups
   // see M.D'Elia and C. Bonati Nuc. PHys B 877 (2012) 233-259
   for(i=0; i<NCOLOR-1; i++)
      {
      for(j=i+1; j<NCOLOR; j++)
         {
         vec_x[0] = creal(X.comp[m(i,j)]);
         vec_x[1] = cimag(X.comp[m(j,i)]);
         vec_x[2] = 0.5*(creal(X.comp[m(i,i)]) - creal(X.comp[m(j,j)]));

         X_mod = sqrt(vec_x[0]*vec_x[0] + vec_x[1]*vec_x[1] + vec_x[2]*vec_x[2]);

         if(X_mod>MIN_VALUE)
           {
           diagonalize_X_Su2_aux(overrelaxparam, vec_x, &G2);

           gii= G2.comp[0] + (G2.comp[3])*I;
           gij= G2.comp[2] + (G2.comp[1])*I;
           gji=-G2.comp[2] + (G2.comp[1])*I;
           gjj= G2.comp[0] - (G2.comp[3])*I;

           // update the X matrix: X->G2*X*G2dag
           // 1st step X->X*G2dag
           equal_SuN(&aux1, &X);
           for(k=0; k<NCOLOR; k++)
              {
              tmp0=aux1.comp[m(k,i)]*conj(gii)+aux1.comp[m(k,j)]*conj(gij);
              tmp1=aux1.comp[m(k,i)]*conj(gji)+aux1.comp[m(k,j)]*conj(gjj);
              X.comp[m(k,i)]=tmp0;
              X.comp[m(k,j)]=tmp1;
              }

           // 2nd step X->G2*X
           for(k=0; k<NCOLOR; k++)
              {
              tmp0=gii*X.comp[m(i,k)]+gij*X.comp[m(j,k)];
              tmp1=gji*X.comp[m(i,k)]+gjj*X.comp[m(j,k)];
              X.comp[m(i,k)]=tmp0;
              X.comp[m(j,k)]=tmp1;
              }

           // build the gauge transformation matrix G_mag->G2*G_mag
           for(k=0; k<NCOLOR; k++)
              {
              tmp0=gii*G_mag->comp[m(i,k)]+gij*G_mag->comp[m(j,k)];
              tmp1=gji*G_mag->comp[m(i,k)]+gjj*G_mag->comp[m(j,k)];
              G_mag->comp[m(i,k)]=tmp0;
              G_mag->comp[m(j,k)]=tmp1;
              }
           }
         }
      }
   }


// compute the squared absolute values of out-of-diagonal terms of X(n)
void comp_outdiagnorm_of_X_SuN(SuN X_links[2*STDIM],
                               double const lambda[NCOLOR],
                               double *outdiagnorm2)
   {
   int i, j, dir;
   double aux=0;
   SuN X, aux1, aux2, aux3;
   
   zero_SuN(&X);

   for(dir=0; dir<STDIM; dir++) 
      {
      // aux1 = lambda*U^{dag}_{dir}(n) and aux2 = lambda*U_{dir}(n-dir)
      diag_matrix_times_dag_SuN(&aux1, lambda, &(X_links[dir]));
      diag_matrix_times_SuN(&aux2, lambda, &(X_links[dir+STDIM]));
      
      // X = U_{dir}(n)*lambda*U^{dag}_{dir}(n) + U^{dag}_{dir}(n-dir)*lambda*U_{dir}(n-dir)
      // as in C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]
      times_SuN(&aux3, &(X_links[dir]), &aux1);
      plus_equal_SuN(&X, &aux3);
      times_dag1_SuN(&aux3, &(X_links[dir+STDIM]), &aux2);
      plus_equal_SuN(&X, &aux3);
      }   
  
   #ifdef DEBUG
   double herm_real_check=0.0, herm_imag_check=0.0;
   double complex trace_check = 0.0+0.0*I;

   for(i=0; i<NCOLOR; i++)
      {
      trace_check += X.comp[m(i,i)];
      for(j=i+1; j<NCOLOR; j++)
         {
         herm_real_check += creal(X.comp[m(i,j)]) - creal(X.comp[m(j,i)]);
         herm_imag_check += cimag(X.comp[m(i,j)]) + cimag(X.comp[m(j,i)]);
         }
      }

   if(fabs(herm_real_check) > MIN_VALUE*10 || fabs(herm_imag_check) > MIN_VALUE*10 || cabs(trace_check) > MIN_VALUE*10 )
     {
     fprintf(stderr, "herm_real_check = %.12g herm_imag_check = %.12g ", herm_real_check, herm_imag_check);
     fprintf(stderr, "trace_check_real = %.12g trace_check_imag = %.12g ", creal(trace_check), cimag(trace_check));
     fprintf(stderr, "(%s, %d)\n", __FILE__, __LINE__);
     }
   #endif

   //compute the out-of-diag elements
    for(i=0; i<NCOLOR; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          aux += cabs(X.comp[m(i,j)])*cabs(X.comp[m(i,j)]);
          }
       }

   *outdiagnorm2 = aux;
   }


// compute the value of the functional to be maximized in MAG
void comp_functional_fmag_SuN(SuN X_links[2*STDIM],
                              double const lambda[NCOLOR],
                              double *fmag)
   {
   int dir, i;
   SuN F_mag, aux1, aux2, aux3;
   double fmag_aux=0.0;

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
      fmag_aux += creal(F_mag.comp[m(i,i)]);
      }

   *fmag=fmag_aux;
   }

// extract the diagonal part of the link and save it in GC->diag_proj
void diag_projection_single_site_SuN(Gauge_Conf *GC,
                                     SuN *link, 
                                     long r,
                                     int dir)
   {
   int subg;
   double phi[NCOLOR], inv_rho[NCOLOR], dphi, inv_rho_sum, phi_aux;
   double complex u_ii;

   dphi=0;
   inv_rho_sum=0;

   // cycle on subgroups
   for(subg=0; subg<NCOLOR; subg++)
      {
      u_ii = link->comp[m(subg,subg)];
      phi[subg] = atan2(cimag(u_ii), creal(u_ii));

      // see C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ] app. B
      // for an explanation of the following procedure

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
     
   for(subg=0; subg<NCOLOR; subg++)
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
       
      GC->diag_proj[r][dir][subg] = phi_aux;
      }
   }


#endif




