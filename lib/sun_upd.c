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
#include"../include/sun_aux.h"

// Pseudo-heatbath by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
// Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985))
void single_heatbath_SuN(SuN *link, SuN const * const staple, GParam const * const param)
    {
    SuN aux;
    Su2 u, v, w;
    double xi, p0; 
    double complex temp0, temp1;
    double complex fii, fij, fji, fjj;
    int i, j, k;

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          equal_SuN(&aux, staple);     // aux=staple
          times_equal_SuN(&aux, link); // aux=staple*link
          ennetodue(&aux, i, j, &xi, &u);

          xi*=(param->d_beta)*2.0/((double) NCOLOR);

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
            }
          }
       }
    }


// Pseudo-heatbath by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
// Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985))
void single_heatbath_aux_SuN(SuN *link, SuN const * const staple, double beta)
    {
    SuN aux;
    Su2 u, v, w;
    double xi, p0;
    double complex temp0, temp1;
    double complex fii, fij, fji, fjj;
    int i, j, k;

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          equal_SuN(&aux, staple);     // aux=staple
          times_equal_SuN(&aux, link); // aux=staple*link
          ennetodue(&aux, i, j, &xi, &u);

          xi*=beta*2.0/((double) NCOLOR);

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

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          { 
          equal_SuN(&aux, staple);     // aux=staple
          times_equal_SuN(&aux, link); // aux=staple*link
          ennetodue(&aux, i, j, &xi, &u);

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
            }
          }
       }
    }


// cooling
void cool_SuN(SuN *link, SuN const * const staple)
  {
  SuN prod;
  Su2 aux2, aux2bis;
  double complex temp0, temp1;
  double complex fii, fij, fji, fjj;
  double xi; // not used
  int i, j, k;

  equal_SuN(&prod, staple);         // prod=staple
  times_equal_SuN(&prod, link);     // prod=staple*link
 
  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        ennetodue(&prod, i, j, &xi, &aux2);

        equal_dag_Su2(&aux2bis, &aux2);   // aux2bis = (aux2)^{dag}

        fii= aux2bis.comp[0] + (aux2bis.comp[3])*I;
        fij= aux2bis.comp[2] + (aux2bis.comp[1])*I;
        fji=-aux2bis.comp[2] + (aux2bis.comp[1])*I;
        fjj= aux2bis.comp[0] - (aux2bis.comp[3])*I;

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

#endif
