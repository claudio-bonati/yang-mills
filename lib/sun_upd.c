#ifndef SUN_UPD_C
#define SUN_UPD_C

#include<complex.h>
#include<math.h>
#include<stdio.h>

#include"../include/gparam.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"
#include"../include/sun.h"

// given the matrix NCOLOR*NCOLOR "in" extracts the i, j lines and column and
// returns "xi" real number and "u" in SU(2)
// 4 xi^2 = redet2[s-s^(dag)+1*tr(s^(dag))]
// u = [s-s^(dag)+1*tr(s^(dag))]/2/xi
// (see Kennedy, Pendleton Phys. Lett. B 156, 393 (1985))
void ennetodue_SuN(SuN const * const in, int i, int j, double *xi, Su2 *u)
   {
   double s[2][2][2], auxr[2][2], auxi[2][2];
   double p;

   s[0][0][0]=creal(in->comp[m(i,i)]);
   s[0][0][1]=cimag(in->comp[m(i,i)]);

   s[0][1][0]=creal(in->comp[m(i,j)]);
   s[0][1][1]=cimag(in->comp[m(i,j)]);

   s[1][0][0]=creal(in->comp[m(j,i)]);
   s[1][0][1]=cimag(in->comp[m(j,i)]);

   s[1][1][0]=creal(in->comp[m(j,j)]);
   s[1][1][1]=cimag(in->comp[m(j,j)]);

   auxr[0][0]=s[0][0][0]+s[1][1][0];
   auxi[0][0]=s[0][0][1]-s[1][1][1];

   auxr[0][1]=s[0][1][0]-s[1][0][0];
   auxi[0][1]=s[0][1][1]+s[1][0][1];

   auxr[1][0]=s[1][0][0]-s[0][1][0];
   auxi[1][0]=s[1][0][1]+s[0][1][1];

   auxr[1][1]=s[0][0][0]+s[1][1][0];
   auxi[1][1]=s[1][1][1]+s[1][1][1]-s[0][0][1]-s[1][1][1];

   p=auxr[0][0]*auxr[1][1]-auxi[0][0]*auxi[1][1]-auxr[0][1]*auxr[1][0]+auxi[0][1]*auxi[1][0];
   p=sqrt(p);

   (*xi)=p/2.0;

   if(*xi>MIN_VALUE)
     {
     auxr[0][0]/=p;
     auxi[0][1]/=p;
     auxr[0][1]/=p;
     auxi[0][0]/=p;
     }
   u->comp[0]=auxr[0][0];
   u->comp[1]=auxi[0][1];
   u->comp[2]=auxr[0][1];
   u->comp[3]=auxi[0][0];
   }


// given a 2*2 matrix extend to NCOLOR*NCOLOR with 1 on the diagonal
void duetoenne_SuN(Su2 const * const in, int i, int j, SuN *out)
   {
   one_SuN(out);

   out->comp[m(i,i)]= in->comp[0] + (in->comp[3])*I;
   out->comp[m(i,j)]= in->comp[2] + (in->comp[1])*I;
   out->comp[m(j,i)]=-in->comp[2] + (in->comp[1])*I;
   out->comp[m(j,j)]= in->comp[0] - (in->comp[3])*I;
   }


// Pseudo-heatbath by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
// Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985))
//
// generate link according to the distibution \exp[-(1/N_c)ReTr(link*staple)]
// the coupling is inside the staple
//
void single_heatbath_SuN(SuN *link, SuN const * const staple)
    {
    SuN aux;
    Su2 u, v, w;
    double xi, p0; 
    double complex temp0, temp1;
    double complex fii, fij, fji, fjj;
    int i, j, k;

    equal_SuN(&aux, staple);     // aux=staple
    times_equal_SuN(&aux, link); // aux=staple*link

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          ennetodue_SuN(&aux, i, j, &xi, &u);

          xi*=2.0/((double) NCOLOR);

          if(xi>MIN_VALUE)
            {
            randheat_Su2(xi, &p0);

            equal_dag_Su2(&w, &u);   // w=u^{dag}
            rand_matrix_p0_Su2(p0, &v);
            times_equal_Su2(&w, &v); // w*=v
           
            fii= w.comp[0] + (w.comp[3])*I;
            fij= w.comp[2] + (w.comp[1])*I;
            fji=-w.comp[2] + (w.comp[1])*I;
            fjj= w.comp[0] - (w.comp[3])*I;

            // link*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=link->comp[m(k,i)]*fii + link->comp[m(k,j)]*fji;
               temp1=link->comp[m(k,i)]*fij + link->comp[m(k,j)]*fjj;
               link->comp[m(k,i)]=temp0;
               link->comp[m(k,j)]=temp1;
               }

            // aux*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=aux.comp[m(k,i)]*fii + aux.comp[m(k,j)]*fji;
               temp1=aux.comp[m(k,i)]*fij + aux.comp[m(k,j)]*fjj;
               aux.comp[m(k,i)]=temp0;
               aux.comp[m(k,j)]=temp1;
               }
            }
          else
            {
            rand_matrix_Su2(&w);

            fii= w.comp[0] + (w.comp[3])*I;
            fij= w.comp[2] + (w.comp[1])*I;
            fji=-w.comp[2] + (w.comp[1])*I;
            fjj= w.comp[0] - (w.comp[3])*I;

            // link*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=link->comp[m(k,i)]*fii + link->comp[m(k,j)]*fji;
               temp1=link->comp[m(k,i)]*fij + link->comp[m(k,j)]*fjj;
               link->comp[m(k,i)]=temp0;
               link->comp[m(k,j)]=temp1;
               }

            // aux*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=aux.comp[m(k,i)]*fii + aux.comp[m(k,j)]*fji;
               temp1=aux.comp[m(k,i)]*fij + aux.comp[m(k,j)]*fjj;
               aux.comp[m(k,i)]=temp0;
               aux.comp[m(k,j)]=temp1;
               }
            }
          }
       }
    }


// Pseudo-overrelaxation by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
// Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985))
void single_overrelaxation_SuN(SuN *link, SuN const * const staple)
    {
    SuN aux;
    Su2 u,v;
    double xi; 
    double complex temp0, temp1;
    double complex fii, fij, fji, fjj;
    int i, j, k;

    equal_SuN(&aux, staple);     // aux=staple
    times_equal_SuN(&aux, link); // aux=staple*link

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          { 
          ennetodue_SuN(&aux, i, j, &xi, &u);

          if(xi>MIN_VALUE)
            {
            equal_dag_Su2(&v, &u);       // v=u^{dag}
            times_equal_dag_Su2(&v, &u); // v=(u^{dag})^2

            fii= v.comp[0] + (v.comp[3])*I;
            fij= v.comp[2] + (v.comp[1])*I;
            fji=-v.comp[2] + (v.comp[1])*I;
            fjj= v.comp[0] - (v.comp[3])*I;

            //link*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=link->comp[m(k,i)]*fii + link->comp[m(k,j)]*fji;
               temp1=link->comp[m(k,i)]*fij + link->comp[m(k,j)]*fjj;
               link->comp[m(k,i)]=temp0;
               link->comp[m(k,j)]=temp1;
               }

            // aux*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=aux.comp[m(k,i)]*fii + aux.comp[m(k,j)]*fji;
               temp1=aux.comp[m(k,i)]*fij + aux.comp[m(k,j)]*fjj;
               aux.comp[m(k,i)]=temp0;
               aux.comp[m(k,j)]=temp1;
               }
            }
          else
            {
            rand_matrix_Su2(&u);

            fii= u.comp[0] + (u.comp[3])*I;
            fij= u.comp[2] + (u.comp[1])*I;
            fji=-u.comp[2] + (u.comp[1])*I;
            fjj= u.comp[0] - (u.comp[3])*I;

            // link*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=link->comp[m(k,i)]*fii + link->comp[m(k,j)]*fji;
               temp1=link->comp[m(k,i)]*fij + link->comp[m(k,j)]*fjj;
               link->comp[m(k,i)]=temp0;
               link->comp[m(k,j)]=temp1;
               }

            // aux*=final
            for(k=0; k<NCOLOR; k++)
               {
               temp0=aux.comp[m(k,i)]*fii + aux.comp[m(k,j)]*fji;
               temp1=aux.comp[m(k,i)]*fij + aux.comp[m(k,j)]*fjj;
               aux.comp[m(k,i)]=temp0;
               aux.comp[m(k,j)]=temp1;
               }
            }
          }
       }
    }


// cooling
void cool_SuN(SuN *link, SuN const * const staple)
  {
  SuN prod;
  Su2 u, udag;
  double complex temp0, temp1;
  double complex fii, fij, fji, fjj;
  double aux; 
  int i, j, k;

  equal_SuN(&prod, staple);         // prod=staple
  times_equal_SuN(&prod, link);     // prod=staple*link
 
  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        ennetodue_SuN(&prod, i, j, &aux, &u); // aux not used in the following

        equal_dag_Su2(&udag, &u);   // udag = (u)^{dag}

        fii= udag.comp[0] + (udag.comp[3])*I;
        fij= udag.comp[2] + (udag.comp[1])*I;
        fji=-udag.comp[2] + (udag.comp[1])*I;
        fjj= udag.comp[0] - (udag.comp[3])*I;

        // link*=final
        for(k=0; k<NCOLOR; k++)
           {
           temp0=link->comp[m(k,i)]*fii + link->comp[m(k,j)]*fji;
           temp1=link->comp[m(k,i)]*fij + link->comp[m(k,j)]*fjj;
           link->comp[m(k,i)]=temp0;
           link->comp[m(k,j)]=temp1;
           }

        // prod*=final
        for(k=0; k<NCOLOR; k++)
           {
           temp0=prod.comp[m(k,i)]*fii + prod.comp[m(k,j)]*fji;
           temp1=prod.comp[m(k,i)]*fij + prod.comp[m(k,j)]*fjj;
           prod.comp[m(k,i)]=temp0;
           prod.comp[m(k,j)]=temp1;
           }
        }
     }

  // multiply by the center element that maximizes ReTr[staple * link]
  aux = argtr_SuN(&prod);                         // aux = phase of ReTr[staple * link]
  aux = round(aux / PI2_N ) * PI2_N;              // round aux to nearest center phase 
  times_equal_complex_SuN(link, cexp(-I * aux));  // link *= exp(-i * aux)
  }


void single_overrelaxation_SuNVecs(SuNVecs *restrict link, SuNVecs const * const staple)
  {
  SuNVecs newlink;
  double norm;
  double complex scalprod;

  norm=norm_SuNVecs(staple);
  if(norm>MIN_VALUE)
    {
    scalprod=complex_scal_prod_SuNVecs(staple, link);

    equal_SuNVecs(&newlink, staple);
    times_equal_complex_SuNVecs(&newlink, 2.0*scalprod/norm/norm);
    minus_equal_SuNVecs(&newlink, link);

    equal_SuNVecs(link, &newlink);
    }
  else
    {
    rand_vec_SuNVecs(link);
    }
  }


#endif
