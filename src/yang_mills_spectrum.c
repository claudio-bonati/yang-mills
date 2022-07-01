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
//         |----<------+
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
       k1=nnp(geo, r, i);

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

   times_equal_real(&stap, blockcoeff/(2.0*(double)(STDIM-2)));
   times_equal_real(M, 1.0-blockcoeff);
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


void space_polyakov(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    long r,
                    double *repoly)
   {
   int i;
   GAUGE_GROUP matrix;

   one(&matrix);
   for(i=0; i<geo->d_size[1]; i++) //Calcolo il loop spaziale lungo la direzione 1
      {
      times_equal(&matrix, &(GC->lattice[r][1]));
      r=nnp(geo, r, 1);
      }
         
   *repoly = retr(&matrix);  
   }


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


void spatial_ape_smearing(Gauge_Conf const * const GC,
                          Geometry const * const geo,
                          double alpha,
                          int smearing_steps)
  {
  int i, step;
  long r;
  GAUGE_GROUP M;
  Gauge_Conf staple_GC;
  
  init_gauge_conf_from_gauge_conf(&staple_GC, GC, geo);

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
           times_equal_real(&(staple_GC.lattice[r][i]), alpha/(2.0*(double)(STDIM-2)));
           times_equal_real(&GC->lattice[r][i], 1.0-alpha);
           plus_equal_dag(&GC->lattice[r][i], &(staple_GC.lattice[r][i]));
           unitarize(&GC->lattice[r][i]); 
           }
        }
     }

  free_gauge_conf(&staple_GC, geo);
  }


void compute_and_print_corr_for_spectrum(Gauge_Conf const * const block_GC,
                                         Geometry const * const blockgeo,
                                         FILE *datafilep)
  {
  int i, j, t, t1, t2;
  long rsp, r;
  double *plaq_re, *poly_re;
  double replaq, av_plaq_re, repoly, av_poly_re;
  double corr_plaq_re_re, corr_poly_re_re;

  plaq_re = (double*) malloc((long unsigned int)(blockgeo->d_size[0])*sizeof(double));
  poly_re = (double*) malloc((long unsigned int)(blockgeo->d_size[0])*sizeof(double));

  for(t=0; t<blockgeo->d_size[0]; t++)
     {
     plaq_re[t]=0.0;
     poly_re[t]=0.0;
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
              plaquette_re(block_GC, blockgeo, r, i, j, &replaq);
              plaq_re[t] += replaq;  //somma placchette al tempo t
              }
           }
        }
     plaq_re[t]/=(double) blockgeo->d_space_vol;
     plaq_re[t]/= (double) ((STDIM-1)*(STDIM-2)/2);
     }

  av_plaq_re = 0.0;
  for(t = 0; t < blockgeo->d_size[0]; t++)
     {
     av_plaq_re += plaq_re[t];
     }
  av_plaq_re /= (double)blockgeo->d_size[0];

  //toreloni
  for(t = 0; t<blockgeo->d_size[0]; t++)
     {
     for(rsp = 0; rsp < blockgeo->d_space_vol/blockgeo->d_size[1]; rsp++)
        {
        r = sisp_and_t_to_si(blockgeo, rsp, t);
        space_polyakov(block_GC, blockgeo, r, &repoly);
        poly_re[t] += repoly;
        }
     poly_re[t] /= (double) (blockgeo->d_space_vol/blockgeo->d_size[1]);
     }

  av_poly_re = 0.0;
  for(t = 0; t < blockgeo->d_size[0]; t++)
     {
     av_poly_re += poly_re[t];
     }
  av_poly_re /= (double)blockgeo->d_size[0];

  // Faccio la media dei correlatori della parte reale delle plaquette a distanza fissa t e la stampo
  for(t = 0; t<blockgeo->d_size[0]/2; t++)
     {
     corr_plaq_re_re = 0.0;

     for(t1 = 0; t1<blockgeo->d_size[0]; t1++)
        {
        t2=(t1+t) % blockgeo->d_size[0];
        corr_plaq_re_re += plaq_re[t2]*plaq_re[t1];
        }
     corr_plaq_re_re/=(double) blockgeo->d_size[0];

     fprintf(datafilep, " %.12f", corr_plaq_re_re);
     }

  // Faccio la media dei correlatori della parte reale dei polyakov spaziali a distanza fissa t e la stampo
  for(t = 0; t<blockgeo->d_size[0]/2; t++)
     {
     corr_poly_re_re = 0.0;

     for(t1 = 0; t1<blockgeo->d_size[0]; t1++)
        {
        t2=(t1+t) % blockgeo->d_size[0];
        corr_poly_re_re += poly_re[t2]*poly_re[t1];
        }
     corr_poly_re_re/=(double) blockgeo->d_size[0];

     fprintf(datafilep, " %.12f", corr_poly_re_re);
     }

  // stampo le medie
  fprintf(datafilep, " %.12f", av_plaq_re);
  fprintf(datafilep, " %.12f", av_poly_re);

  fprintf(datafilep, "\n");
  fflush(datafilep);

  free(plaq_re);
  free(poly_re);
  }


void perform_measures_spectrum(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               FILE *datafilep)
  {
  int i, sizes[STDIM];
  Gauge_Conf smeared_GC, block_GC;
  Geometry blockgeo;

  // smearing parameter
  const double alpha=0.2;

  // smearing steps
  const int smearing_steps = 2;

  // copy of the initial configuration
  init_gauge_conf_from_gauge_conf(&smeared_GC, GC, geo);

  // perform smearing
  spatial_ape_smearing(&smeared_GC, geo, alpha, smearing_steps);

  // geometry for blocking
  sizes[0]=geo->d_size[0];
  for(i=1; i<STDIM; i++)
     {
     sizes[i]=geo->d_size[i]/2;
     }
  init_geometry(&blockgeo, sizes);

  // blocking
  init_spatial_blocked_conf(&block_GC,
                            &blockgeo,
                            &smeared_GC,
                            geo,
                            alpha);

  compute_and_print_corr_for_spectrum(&block_GC,
                                      &blockgeo,
                                      datafilep);

  free_gauge_conf(&smeared_GC, geo);
  free_gauge_conf(&block_GC, &blockgeo);
  free_geometry(&blockgeo);
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
       printf("%d ", count);
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
    print_parameters_local(&param, time1, time2);

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
