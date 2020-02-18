#ifndef SON_UPD_C
#define SON_UPD_C

#include<complex.h>
#include<math.h>
#include<stdio.h>

#include"../include/gparam.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/son.h"


// given the matrix NCOLOR*NCOLOR "in" extracts the i, j lines and column and
// returns "xi" real number and "u" in SO(2) parametrized by
// u = rotation by an angle "theta"
//
// see e.g. Richard Lau "SO(N) gauge theories in 2+1 dimensions" PhD thesis (Oxford)
//
// "u" si such that for any 2x2 matrix m we have "tr(m*in)=xi*tr(m*u)
void ennetodue_SoN(SoN const * const in, int i, int j, double *xi, double *theta)
   {
   double m00, m01, m10, m11;

   m00=(in->comp[m(i,i)]+in->comp[m(j,j)])/2.0;
   m11=m00;
   m01=(in->comp[m(i, j)]-in->comp[m(j, i)])/2.0;
   m10=-m01;

   *theta=atan2(m01, m00);

   *xi=sqrt(m00*m11-m01*m10);
   }


// given a 2*2 matrix parametrized by the angle theta it extends it to NCOLOR*NCOLOR with 1 on the diagonal
void duetoenne_SoN(double theta, int i, int j, SoN *out)
   {
   one_SoN(out);

   out->comp[m(i,i)]= cos(theta);
   out->comp[m(i,j)]= sin(theta);
   out->comp[m(j,i)]=-sin(theta);
   out->comp[m(j,j)]= cos(theta);
   }


// Pseudo-heatbath by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982))
//
// generate link according to the distibution \exp[-(1/N_c)ReTr(link*staple)]
// the coupling is inside the staple
//
void single_heatbath_SoN(SoN *link, SoN const * const staple)
    {
    SoN aux;
    double xi, theta_old, theta_new, y1, y2, prob;
    double fii, fij, fji, fjj, temp0, temp1;
    int i, j, k;

    equal_SoN(&aux, staple);     // aux=staple
    times_equal_SoN(&aux, link); // aux=staple*link

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          ennetodue_SoN(&aux, i, j, &xi, &theta_old);

          xi*=(2.0/(double)NCOLOR);

          if(xi>MIN_VALUE)
            {
            y1=casuale();
            y2=casuale();

            theta_new=sqrt( -2.0/xi*log(1 - y1*(1-exp(-xi*PI*PI/2.0)) ) )*cos(PI2*(y2-0.5));
            prob=exp( xi * (cos(theta_new) -theta_old*theta_old/2.0 -cos(theta_old) +theta_new*theta_new/2.0 ));

            if(casuale()<prob)
              {
              fii= cos(theta_new-theta_old);
              fij= sin(theta_new-theta_old);
              fji=-sin(theta_new-theta_old);
              fjj= cos(theta_new-theta_old);

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
          else
            {
            theta_new=2*PI*casuale();

            fii= cos(theta_new);
            fij= sin(theta_new);
            fji=-sin(theta_new);
            fjj= cos(theta_new);

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
void single_overrelaxation_SoN(SoN *link, SoN const * const staple)
    {
    SoN aux;
    double xi, theta, theta_new;
    double fii, fij, fji, fjj, temp0, temp1;
    int i, j, k;

    equal_SoN(&aux, staple);     // aux=staple
    times_equal_SoN(&aux, link); // aux=staple*link

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          ennetodue_SoN(&aux, i, j, &xi, &theta);

          if(xi>MIN_VALUE)
            {
            fii= cos(-2.0*theta);
            fij= sin(-2.0*theta);
            fji=-sin(-2.0*theta);
            fjj= cos(-2.0*theta);

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
            theta_new=2*PI*casuale();

            fii= cos(theta_new);
            fij= sin(theta_new);
            fji=-sin(theta_new);
            fjj= cos(theta_new);

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
void cool_SoN(SoN *link, SoN const * const staple)
  {
  SoN prod;
  double fii, fij, fji, fjj, temp0, temp1;
  double xi, theta; // xi not used
  int i, j, k;

  equal_SoN(&prod, staple);         // prod=staple
  times_equal_SoN(&prod, link);     // prod=staple*link

  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        ennetodue_SoN(&prod, i, j, &xi, &theta);

        fii= cos(-theta);
        fij= sin(-theta);
        fji=-sin(-theta);
        fjj= cos(-theta);

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
  }


void single_overrelaxation_SoNVecs(SoNVecs *restrict link, SoNVecs const * const staple)
  {
  SoNVecs newlink;
  double norm, scalprod;

  norm=norm_SoNVecs(staple);
  if(norm>MIN_VALUE)
    {
    scalprod=re_scal_prod_SoNVecs(link, staple);

    equal_SoNVecs(&newlink, staple);
    times_equal_real_SoNVecs(&newlink, 2.0*scalprod/norm/norm);
    minus_equal_SoNVecs(&newlink, link);

    equal_SoNVecs(link, &newlink);
    }
  else
    {
    rand_vec_SoNVecs(link);
    }
  }

#endif
