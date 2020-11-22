#ifndef SU2_UPD_C
#define SU2_UPD_C

#include<math.h>
#include<stdio.h>

#include"../include/gparam.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"

// random number generator for heabath (see Kennedy, Pendleton Phys. Lett. B 156, 393 (1985))
// given "k" return "a" in [-1,1] with P(a) = sqrt(1-a*a)*exp(k*a)
void randheat_Su2(double k, double *out)
    {
    double r, r1, r2, c, r3;

    if(k>1.6847)
      {
      // kennedy pendleton
      while(1==1)
        {
        r=casuale();
        r1=casuale();

        r=-log(r);  
        r/=k;
        r1=-log(r1); 
        r1/=k;
 
        r2=casuale();
        c=cos(r2*PI2);
        c*=c;

        r*=c;  
        r1+=r;  

        r3=casuale();
        r=1-r1*0.5-r3*r3;
        if(r>0)
          {
          *out=1-r1;
          break;
          }
        }
      }
    else
      {
      // creutz
      r1=exp(-2.0*k);
      while(1==1)
        {
        do
        r=casuale();
        while(r<r1);
        r=1.0+log(r)/k;
 
        r2=casuale();
        if(r*r<r2*(2.0-r2))
          {
          *out=r;
          break;
          }
        }
      }
    }


// heatbath
//
// generate link according to the distibution \exp[-(1/N_c)ReTr(link*staple)]
// the coupling is inside the staple
//
void single_heatbath_Su2(Su2 *link, Su2 const * const staple)
    {
    double p, p0;
    Su2 matrix1, matrix2;

    equal_Su2(&matrix1, staple);                    // matrix1=staple

    p=sqrtdet_Su2(&matrix1);
    if(p>MIN_VALUE)
      {
      times_equal_real_Su2(&matrix1, 1.0/p);  // matrix1 *= 1.0/p
      equal_dag_Su2(&matrix2, &matrix1);      // matrix2 = matrix1^{dag}

      randheat_Su2(p, &p0);
 
      rand_matrix_p0_Su2(p0, link);
      times_equal_Su2(link, &matrix2);        // link*=matrix2
      }
    else
      {
      rand_matrix_Su2(link);
      }
    }


// overrelaxation
void single_overrelaxation_Su2(Su2 *link, Su2 const * const staple)
    {
    double p;
    Su2 matrix1, matrix2;

    equal_Su2(&matrix1, staple); // matrix1=staple

    p=sqrtdet_Su2(&matrix1);   
    if(p>MIN_VALUE)

      {
      times_equal_real_Su2(&matrix1, 1.0/p);     // matrix1*=(1/p)
      equal_dag_Su2(&matrix2, &matrix1);         // matrix2=matrix1^{dag}

      times_dag1_Su2(&matrix1, &matrix2, link);  // matrix1=matrix2*link^{dag}
      times_Su2(link, &matrix1, &matrix2);       // link=matrix1*matrix2
      }
    else
      {
      rand_matrix_Su2(link);
      }
    }


// cooling
void cool_Su2(Su2 *link, Su2 const * const staple)
    {
    Su2 matrix1;
    
    equal_Su2(&matrix1, staple);   // matrix1=staple
    unitarize_Su2(&matrix1);

    equal_dag_Su2(link, &matrix1);  // link=matrix1^{dag}
    }


void single_overrelaxation_Su2Vecs(Su2Vecs *restrict link, Su2Vecs const * const staple)
  {
  Su2Vecs newlink;
  double norm, scalprod;

  norm=norm_Su2Vecs(staple);
  if(norm>MIN_VALUE)
    {
    scalprod=re_scal_prod_Su2Vecs(link, staple);

    equal_Su2Vecs(&newlink, staple);
    times_equal_real_Su2Vecs(&newlink, 2.0*scalprod/norm/norm);
    minus_equal_Su2Vecs(&newlink, link);

    equal_Su2Vecs(link, &newlink);
    }
  else
    {
    rand_vec_Su2Vecs(link);
    }
  }


#endif
