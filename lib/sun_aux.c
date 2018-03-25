#ifndef SUN_AUX_C
#define SUN_AUX_C

#include<complex.h>
#include<math.h>

#include"../include/macro.h"
#include"../include/su2.h"
#include"../include/sun.h"

// given the matrix NCOLOR*NCOLOR "in" extracts the i, j lines and column and
// gives "xi" real number and "u" in SU(2)
// 4 xi^2 = redet2[s-s^(dag)+1*tr(s^(dag))]
// u = [s-s^(dag)+1*tr(s^(dag))]/2/xi
// (see Kennedy, Pendleton Phys. Lett. B 156, 393 (1985))
void ennetodue(SuN const * const in, int i, int j, double *xi, Su2 *u)
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
void duetoenne(Su2 const * const in, int i, int j, SuN *out)
   {
   one_SuN(out);

   out->comp[m(i,i)]= in->comp[0] + (in->comp[3])*I;
   out->comp[m(i,j)]= in->comp[2] + (in->comp[1])*I;
   out->comp[m(j,i)]=-in->comp[2] + (in->comp[1])*I;
   out->comp[m(j,j)]= in->comp[0] - (in->comp[3])*I;
   }

#endif
