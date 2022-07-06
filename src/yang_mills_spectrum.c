#ifndef YM_SPECTRUM_C
#define YM_SPECTRUM_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"


// compute blocked link for a given site
void spatialblocking_singlesite(Gauge_Conf const * const GC,
                                Geometry const * const geo,
                                long r,
                                int i,
                                double blockcoeff,
                                GAUGE_GROUP *M)
  {
  int j, l;
  long k, k1;

  GAUGE_GROUP link1, link2, link3, link4, staptmp, stap;

  #if DEBUG
  if(r >= geo->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i == 0)
    {
    fprintf(stderr, "time direction selected: i=%d (%s, %d)\n", i, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  zero(&stap); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     if(j!=0)
       {
//
//       i ^
//         |    (1)
//         +----->-----+
//         |           |
//         |           V (2)
//         |           |
//         |           |
//       k +-----------+
//         |           |
//         |           |
//         |           V (3)
//         |           |
//         +-----<-----+-->   j
//       r     (4)
//

       k=nnp(geo, r, i);
       equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[nnp(geo, k, j)][i]));  // link2 = (2)
       equal(&link3, &(GC->lattice[nnp(geo, r, j)][i]));  // link3 = (3)
       equal(&link4, &(GC->lattice[r][j]));               // link3 = (4)

       times_dag2(&staptmp, &link1, &link2);   // staptmp=link1*link2^{dag}
       times_equal_dag(&staptmp, &link3);      // staptmp*=link3^{dag}
       times_equal_dag(&staptmp, &link4);      // staptmp*=link4^{dag}

       plus_equal(&stap, &staptmp);

//
//       i ^
//         |   (1)
//         +----<------+
//         |           |
//     (2) V           |
//         |           |
//      k1 +-----------+
//         |           |
//     (3) V           |
//         |           |
//         +------>----+--->j
//        k     (4)    r
//

       k=nnm(geo, r, j);
       k1=nnp(geo, k, i);

       equal(&link1, &(GC->lattice[nnp(geo, k1, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[k1][i]));               // link2 = (2)
       equal(&link3, &(GC->lattice[k][i]));                // link3 = (3)
       equal(&link4, &(GC->lattice[k][j]));                // link4 = (4)

       times_dag12(&staptmp, &link1, &link2); // stap=link1^{dag}*link2^{dag}
       times_equal_dag(&staptmp, &link3);     // stap*=link3^{dag}
       times_equal(&staptmp, &link4);         // stap*=link4

       plus_equal(&stap, &staptmp);
       }
     }

   equal(M, &(GC->lattice[r][i]));
   times_equal(M, &(GC->lattice[nnp(geo, r, i)][i]));

   times_equal_real(&stap, blockcoeff);
   plus_equal_dag(M, &stap);
   unitarize(M);
   }


// create spatially blocked configurations (i.e. L_spatial->L_spatial/2)
void init_spatial_blocked_conf(Gauge_Conf *blockGC,
                               Geometry const * const blockgeo,
                               Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               double blockcoeff)
  {

  long rb, r;
  int blockcart[STDIM], cart[STDIM];
  int i, mu, err;
  GAUGE_GROUP M;

  for(i=1; i<STDIM; i++)
     {
     if(geo->d_size[i] % 2 != 0)
       {
       fprintf(stderr, "Problem with spatial size not even: %d ! (%s, %d)\n", geo->d_size[i], __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // allocate the lattice
  err=posix_memalign((void**)&(blockGC->lattice), (size_t) DOUBLE_ALIGN, (size_t) blockgeo->d_volume * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(blockgeo->d_volume); r++)
     {
     err=posix_memalign((void**)&(blockGC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t) STDIM * sizeof(GAUGE_GROUP));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // initialize GC
  for(rb=0; rb<(blockgeo->d_volume); rb++)
     {
     si_to_cart(blockcart, rb, blockgeo);
     cart[0]=blockcart[0];
     for(i=1; i<STDIM; i++)
        {
        cart[i]=2*blockcart[i];
        }
     r=cart_to_si(cart, geo);

     equal(&(blockGC->lattice[rb][0]), &(GC->lattice[r][0]) );

     for(mu=1; mu<STDIM; mu++)
        {
        // this is the real point where the blocking is performed
        spatialblocking_singlesite(GC, geo, r, mu, blockcoeff, &M);
        equal(&(blockGC->lattice[rb][mu]), &M);
        }
     }

  blockGC->update_index=GC->update_index;
  }


// compute real part of the plaquette
void plaquette_re(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  long r,
                  int i,
                  int j,
                  double *replaq)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= geo->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   *replaq =  retr(&matrix);
   }


// compute spatial polyakov loop
void space_polyakov(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    long r,
                    double complex *poly)
   {
   int i;
   GAUGE_GROUP matrix;

   one(&matrix);
   for(i=0; i<geo->d_size[1]; i++) //Calcolo il loop spaziale lungo la direzione 1
      {
      times_equal(&matrix, &(GC->lattice[r][1]));
      r=nnp(geo, r, 1);
      }
         
   *poly = retr(&matrix)+imtr(&matrix);
   }


// compute staples with no time direction
void calcstaples_wilson_no_time(Gauge_Conf const * const GC,
                                Geometry const * const geo,
                                long r,
                                int i,
                                GAUGE_GROUP *M)
  {
  int j, l;
  long k;
  GAUGE_GROUP link1, link2, link3, link12, stap;

  #if DEBUG
  if(r >= geo->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i == 0)
    {
    fprintf(stderr, "time direction selected: i=%d (%s, %d)\n", i, __FILE__, __LINE__); 
    exit(EXIT_FAILURE);
    }
  #endif

  zero(M); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     if(j!=0)
       { 

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

       equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
       equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

       times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
       times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

       plus_equal(M, &stap);

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |
//     (2) V           |
//         |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

       k=nnm(geo, r, j);

       equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
       equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

       times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
       times(&stap, &link12, &link3);        // stap=link12*link3

       plus_equal(M, &stap);
       }
     }
   }


// perform smearing on spatial links
void spatial_smearing(Gauge_Conf const * const GC,
                      Geometry const * const geo,
                      double alpha,
                      int smearing_steps)
  {
  int i, step;
  long r;
  GAUGE_GROUP M;
  Gauge_Conf staple_GC;
  
  init_gauge_conf_from_gauge_conf_noclover(&staple_GC, GC, geo);

  for(step=0; step<smearing_steps; step++)
     {
     for(r = 0; r < geo->d_volume; r++)
        {
        for(i = 1; i < STDIM; i++)
           {
           calcstaples_wilson_no_time(GC, geo, r, i, &M);
           equal(&(staple_GC.lattice[r][i]), &M);
           }
        }

     for(r = 0; r < geo->d_volume; r++)
        {
        for(i = 1; i < STDIM; i++)
           {
           times_equal_real(&(staple_GC.lattice[r][i]), alpha);
           plus_equal_dag(&GC->lattice[r][i], &(staple_GC.lattice[r][i]));
           unitarize(&GC->lattice[r][i]); 
           }
        }
     }

  free_gauge_conf_noclover(&staple_GC, geo);
  }


// compute p=0 averages of the plaquette and polyakov operators
void compute_operators_for_spectrum(Gauge_Conf const * const block_GC,
                                    Geometry const * const blockgeo,
                                    double *plaq,
                                    double *av_plaq,
                                    double complex *poly)
  {
  int i, j, t;
  long rsp, r;
  double plaq_loc;
  double complex poly_loc;

  for(t=0; t<blockgeo->d_size[0]; t++)
     {
     plaq[t]=0.0;
     poly[t]=0.0+I*0.0;
     }

  //Calcolo la somma delle placchette spaziali al tempo t
  for(t = 0; t<blockgeo->d_size[0]; t++)
     {
     for(rsp = 0; rsp < blockgeo->d_space_vol; rsp++)//rsp è l'indice spaziale
        {
        r = sisp_and_t_to_si(blockgeo, rsp, t); //Dall'indice spaziale e da quello temporale ottengo l'indice 4D
        for(i = 1; i < STDIM; i++)
           {
           for(j = i + 1; j < STDIM; j++)
              {
              plaquette_re(block_GC, blockgeo, r, i, j, &plaq_loc);
              plaq[t] += plaq_loc;  //somma placchette al tempo t
              }
           }
        }
     plaq[t]/=(double) blockgeo->d_space_vol;
     plaq[t]/= (double) ((STDIM-1)*(STDIM-2)/2);
     }

  *av_plaq = 0.0;
  for(t = 0; t < blockgeo->d_size[0]; t++)
     {
     *av_plaq += plaq[t];
     }
  *av_plaq /= (double)blockgeo->d_size[0];

  //toreloni
  for(t = 0; t<blockgeo->d_size[0]; t++)
     {
     for(rsp = 0; rsp < blockgeo->d_space_vol; rsp++)
        {
        r = sisp_and_t_to_si(blockgeo, rsp, t);
        space_polyakov(block_GC, blockgeo, r, &poly_loc);
        poly[t] += poly_loc;
        }
     poly[t] /= (double complex) (blockgeo->d_space_vol);
     }
  }


// print correlators of the observables
void compute_and_print_corr_for_spectrum(double **plaq,
                                         double *av_plaq,
                                         double complex **poly,
                                         Geometry const * const geo,
                                         int numblocks,
                                         FILE *datafilep)
  {
  int i, j, t, t1, t2;
  double corr_plaq, corr_poly;

  // Faccio la media dei correlatori della parte reale delle plaquette a distanza fissa t e la stampo
  for(t = 0; t<geo->d_size[0]/2; t++)
     {
     for(i=0; i<numblocks+1; i++)
        {
        for(j=i; j<numblocks+1; j++)
           {
           corr_plaq = 0.0;

           for(t1 = 0; t1<geo->d_size[0]; t1++)
              {
              t2=(t1+t) % geo->d_size[0];
              corr_plaq += plaq[i][t2]*plaq[j][t1];
              }
           corr_plaq/=(double) geo->d_size[0];

           fprintf(datafilep, " %.12f", corr_plaq);
           }
        }
     }

  // Faccio la media dei correlatori dei polyakov spaziali a distanza fissa t e la stampo
  for(t = 0; t<geo->d_size[0]/2; t++)
     {
     for(i=0; i<numblocks+1; i++)
        {
        for(j=i; j<numblocks+1; j++)
           {
           corr_poly = 0.0;

           for(t1 = 0; t1<geo->d_size[0]; t1++)
              {
              t2=(t1+t) % geo->d_size[0];
              corr_poly += creal(conj(poly[i][t2])*poly[j][t1]);
              }
           corr_poly/=(double) geo->d_size[0];

           fprintf(datafilep, " %.12f", corr_poly);
           }
        }
     }

  // stampo le medie
  for(i=0; i<numblocks+1; i++)
     {
     fprintf(datafilep, " %.12f", av_plaq[i]);
     }

  fprintf(datafilep, "\n");
  fflush(datafilep);
  }



void perform_measures_spectrum(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               FILE *datafilep)
  {
  // number of blocking steps, must be >1
  #define NUMBLOCK 2

  // smearing parameter
  const double alpha=0.5;

  // smearing steps
  const int smearing_steps = 3;

  int i, j, sizes[STDIM];
  double complex *poly[NUMBLOCK+1];
  double *plaq[NUMBLOCK+1], av_plaq[NUMBLOCK+1];

  Gauge_Conf smeared_GC, block_GC[NUMBLOCK];
  Geometry blockgeo[NUMBLOCK];

  for(i=0; i<NUMBLOCK+1; i++)
     {
     plaq[i] = (double *) malloc((long unsigned int)(geo->d_size[0])*sizeof(double));
     poly[i] = (double complex *) malloc((long unsigned int)(geo->d_size[0])*sizeof(double complex));
     }

  // copy of the initial configuration
  init_gauge_conf_from_gauge_conf_noclover(&smeared_GC, GC, geo);

  // perform smearing
  spatial_smearing(&smeared_GC, geo, alpha, smearing_steps);

  compute_operators_for_spectrum(&smeared_GC,
                                 geo,
                                 plaq[0],
                                 &(av_plaq[0]),
                                 poly[0]);

  // geometry for blocking
  sizes[0]=geo->d_size[0];
  for(i=1; i<STDIM; i++)
     {
     sizes[i]=geo->d_size[i]/2;
     }
  init_geometry(&(blockgeo[0]), sizes);

  // blocking
  init_spatial_blocked_conf(&(block_GC[0]),
                            &(blockgeo[0]),
                            &smeared_GC,
                            geo,
                            alpha);

  compute_operators_for_spectrum(&block_GC[0],
                                 &blockgeo[0],
                                 plaq[1],
                                 &(av_plaq[1]),
                                 poly[1]);

  for(i=1; i<NUMBLOCK; i++)
     {
     // geometry for blocking
     for(j=1; j<STDIM; j++)
        {
        sizes[j]/=2;
        }
     init_geometry(&(blockgeo[i]), sizes);

     // blocking
     init_spatial_blocked_conf(&(block_GC[i]),
                               &(blockgeo[i]),
                               &(block_GC[i-1]),
                               &(blockgeo[i-1]),
                               alpha);

     compute_operators_for_spectrum(&block_GC[i],
                                    &blockgeo[i],
                                    plaq[i+1],
                                    &(av_plaq[i+1]),
                                    poly[i+1]);
     }

  compute_and_print_corr_for_spectrum(plaq,
                                      av_plaq,
                                      poly,
                                      geo,
                                      NUMBLOCK,
                                      datafilep);

  free_gauge_conf_noclover(&smeared_GC, geo);
  for(i=NUMBLOCK-1; i>=0; i--)
     {
     free_gauge_conf_noclover(&(block_GC[i]), &(blockgeo[i]));
     free_geometry(&(blockgeo[i]));
     }

  for(i=0; i<NUMBLOCK+1; i++)
     {
     free(plaq[i]);
     free(poly[i]);
     }

  #undef NUMBLOCK
  }


void real_main(char *in_file)
    {
    Gauge_Conf GC;
    Geometry geo;
    GParam param;

    int count;
    FILE *datafilep; 
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize geometry
    init_geometry(&geo, param.d_sizeg);

    // initialize gauge configuration
    init_gauge_conf(&GC, &geo, &param);

    // montecarlo
    time(&time1);
    // count starts from 1 to avoid problems using %
    for(count=1; count < param.d_sample + 1; count++)
       {
       update(&GC, &geo, &param);

       if(count % param.d_measevery ==0 && count >= param.d_thermal)
         {
         perform_measures_spectrum(&GC, &geo, datafilep);
         }

       // save configuration for backup
       if(param.d_saveconf_back_every!=0)
         {
         if(count % param.d_saveconf_back_every == 0 )
           {
           // simple
           write_conf_on_file(&GC, &geo, &param);

           // backup copy
           write_conf_on_file_back(&GC, &geo, &param);
           }
         }
      }

    time(&time2);
    // montecarlo end

    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      write_conf_on_file(&GC, &geo, &param);
      }

    // print simulation details
    print_parameters_spectrum(&param, time1, time2);

    // free gauge configuration
    free_gauge_conf(&GC, &geo);

    // free geometry
    free_geometry(&geo);
    }


void print_template_input(void)
  {
  FILE *fp;

  fp=fopen("template_input.in", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "size 4 4 4 4\n");
    fprintf(fp,"\n");
    fprintf(fp, "beta 5.705\n");
    fprintf(fp, "theta 1.5\n");
    fprintf(fp,"\n");
    fprintf(fp, "sample    10\n");
    fprintf(fp, "thermal   0\n");
    fprintf(fp, "overrelax 5\n");
    fprintf(fp, "measevery 1\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                   0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every     5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp,"\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    #(0=time)\n");
    fclose(fp);
    }
  }


int main (int argc, char **argv)
    {
    char in_file[50];

    if(argc != 2)
      {
      printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tGGROUP: %s\n", QUOTEME(GGROUP));
      printf("\tN_c (number of colors): %d\n", NCOLOR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\tNum_levels (number of levels): %d\n", NLEVELS);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      #ifdef THETA_MODE
        printf("\n\tusing imaginary theta\n");
      #endif

      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      print_template_input();


      return EXIT_SUCCESS;
      }
    else
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
        }
      else
        {
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    
    } 

#endif
