#ifndef U1_UPD_C
#define U1_UPD_C

#include<complex.h>
#include<math.h>
#include<stdio.h>

#include"../include/macro.h"
#include"../include/random.h"
#include"../include/u1.h"
#include"../include/u1_upd.h"

// heatbath ( in fact it is not exactly an heatbath, but an approximate heatbath
//            corrected by a Metropolis step )
//
// generate link according to the distibution \exp[ReTr(link*staple)]
// the coupling is inside the staple
//
// based on Giovanni Iannelli master thesis (october 2018)
// In practice it approximates cos(x)~1-x^2/2 and then correct
// with a Metropolis step
//
void single_heatbath_U1(U1 * restrict link, U1 const * const restrict staple)
  {
  double k, xold, y1, y2, xnew, prob;
  U1 aux;

  k=norm_U1(staple);

  if(k>MIN_VALUE)
    {
    times_U1(&aux, link, staple);
    xold=atan2(cimag(aux.comp), creal(aux.comp));

    y1=casuale();
    y2=casuale();

    xnew=sqrt( -2./k*log(1 - y1*(1-exp(-k*PI*PI/2.0)) ) )*cos(PI2*(y2-0.5));
    prob=exp( k * (cos(xnew) -xold*xold/2.0 -cos(xold) +xnew*xnew/2.0 ));

    if(casuale()<prob)
      {
      link->comp = cos(xnew)+sin(xnew)*I;
      times_equal_dag_U1(link, staple);
      times_equal_real_U1(link, 1./k);
      }
    }
  else
    {
    rand_matrix_U1(link);
    }
  }


/*
// heatbath by Moriarty Phys. Rev. D 25, p2185 (1982)
//
// very inefficient for large values of the lattice coupling
//
void single_heatbath_U1(U1 * restrict link, U1 const * const restrict staple)
  {
  double alpha, theta, theta_staple, q, q_max, r1, r2, x;

  alpha=norm_U1(staple);

  if(alpha>MIN_VALUE)
    {
    q_max=exp(0.210513662353018684327769*alpha);  // see later for this number
    theta_staple=atan2(-cimag(staple->comp), creal(staple->comp));

    do
      {
      r1=casuale();
      r2=casuale();
      x = -1.0 + log(1.0 + (exp(2.0*alpha) - 1.0)*r1 )/alpha;
      q=exp(alpha * ( cos(HALF_PI*(1.0-x)) - x) );
      }
    while(q<=q_max*r2);

    theta = (1.0 - x)*HALF_PI;

    r1=casuale();
    if(r1>0.5)
      {
      theta=-theta;
      }
    theta += theta_staple;

    link->comp = cos(theta)+sin(theta)*I;
    }
  else
    {
    rand_matrix_U1(link);
    }

  //  WorkingPrecision->1000;
  //  Q[x_]:=Exp[Cos[Pi/2*(1-x)]-x]
  //  FindRoot[D[Q[x],x]==0, {x,1}, WorkingPrecision->100]
  //  Out[52]= {x -> 0.5606641805798867176366776048997096707812104519411362714885751166519976969907076829844764496162691569}
  //  N[Log[Q[x]/.%41], 100]
  //  Out[53]= 0.210513662353018684327769435155832317434879346989632455087165428289411464536813734686515674410131331
  }
*/


void single_overrelaxation_U1(U1 * restrict link, U1 const * const restrict staple)
  {
  double p;
  U1 newlink, helper;

  p=norm_U1(staple);

  if(p>MIN_VALUE)
    {
    equal_U1(&helper, staple);
    times_equal_real_U1(&helper, 1.0/p);

    times_dag12_U1(&newlink, &helper, link);
    times_equal_dag_U1(&newlink, &helper);

    equal_U1(link, &newlink);
    }
  else
    {
    rand_matrix_U1(link);
    }
  }


void cool_U1(U1 * restrict link, U1 const * const restrict staple)
  {
  U1 aux;

  equal_U1(&aux, staple);
  unitarize_U1(&aux);
  equal_dag_U1(link, &aux);
  }


void single_overrelaxation_U1Vecs(U1Vecs *restrict link, U1Vecs const * const staple)
  {
  U1Vecs newlink;
  double norm;
  double complex scalprod;

  norm=norm_U1Vecs(staple);
  if(norm>MIN_VALUE)
    {
    scalprod=complex_scal_prod_U1Vecs(staple, link);

    equal_U1Vecs(&newlink, staple);
    times_equal_complex_U1Vecs(&newlink, 2.0*scalprod/norm/norm);
    minus_equal_U1Vecs(&newlink, link);

    equal_U1Vecs(link, &newlink);
    }
  else
    {
    rand_vec_U1Vecs(link);
    }
  }





#endif

