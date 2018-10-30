bin_PROGRAMS = conf_check \
               debug_rng \
               debug_su2 \
               debug_sun \
               debug_u1 \
               yang_mills_local \
               yang_mills_polycorr \
               yang_mills_polycorr_long \
               yang_mills_t0 \
               yang_mills_tube_conn \
               yang_mills_tube_disc \
               yang_mills_tracedef 

dSFMT = ./include/dSFMT-common.h \
        ./include/dSFMT.h \
        ./include/dSFMT-params11213.h \
        ./include/dSFMT-params1279.h \
        ./include/dSFMT-params132049.h \
        ./include/dSFMT-params19937.h \
        ./include/dSFMT-params216091.h \
        ./include/dSFMT-params2203.h \
        ./include/dSFMT-params4253.h \
        ./include/dSFMT-params44497.h \
        ./include/dSFMT-params521.h \
        ./include/dSFMT-params86243.h \
        ./include/dSFMT-params.h \
        ./lib/dSFMT.c

endian = ./include/endianness.h  \
         ./lib/endianness.c 

funcpoint = ./include/function_pointers.h \
            ./lib/function_pointers.c

gauge_conf = ./include/gauge_conf.h \
             ./lib/gauge_conf_def.c \
             ./lib/gauge_conf_meas.c \
             ./lib/gauge_conf_multilevel.c \
             ./lib/gauge_conf_upd.c 

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c 

su2 = ./include/aligncheck.h \
      ./include/su2.h \
      ./include/su2_upd.h \
      ./lib/su2.c \
      ./lib/su2_upd.c 

sun = ./include/aligncheck.h \
      ./include/sun_aux.h \
      ./include/sun.h \
      ./include/sun_upd.h \
      ./lib/sun_aux.c \
      ./lib/sun.c \
      ./lib/sun_upd.c 

tens_prod = ./include/aligncheck.h \
            ./include/tens_prod.h \
            ./lib/tens_prod.c

u1 = ./include/aligncheck.h \
     ./include/u1.h \
     ./include/u1_upd.h \
     ./lib/u1.c \
     ./lib/u1_upd.c

conf_check_SOURCES=./src/conf_check.c \
                   $(dSFMT) \
                   $(endian) \
                   $(funcpoint) \
                   $(gparam) \
                   $(random) \
                   $(su2) \
                   $(sun) \
                   $(u1) 

debug_rng_SOURCES=./src/debug_rng.c \
                  $(dSFMT) \
                  $(random)

debug_su2_SOURCES=./src/debug_su2.c \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(su2)

debug_sun_SOURCES=./src/debug_sun.c \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(su2) \
                  $(sun)

debug_u1_SOURCES=./src/debug_u1.c \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(u1)

yang_mills_local_SOURCES=./src/yang_mills_local.c \
                          $(dSFMT)  \
                          $(endian) \
                          $(funcpoint) \
                          $(gauge_conf) \
                          $(geometry) \
                          $(gparam) \
                          ./include/macro.h \
                          $(random) \
                          $(su2) \
                          $(sun) \
                          $(tens_prod) \
                          $(u1) 

yang_mills_polycorr_SOURCES=./src/yang_mills_polycorr.c \
                             $(dSFMT)  \
                             $(endian) \
                             $(funcpoint) \
                             $(gauge_conf) \
                             $(geometry) \
                             $(gparam) \
                             ./include/macro.h \
                             $(random) \
                             $(su2) \
                             $(sun) \
                             $(tens_prod) \
                             $(u1) 

yang_mills_polycorr_long_SOURCES=./src/yang_mills_polycorr_long.c \
                                  $(dSFMT)  \
                                  $(endian) \
                                  $(funcpoint) \
                                  $(gauge_conf) \
                                  $(geometry) \
                                  $(gparam) \
                                  ./include/macro.h \
                                  $(random) \
                                  $(su2) \
                                  $(sun) \
                                  $(tens_prod) \
                                  $(u1) 

yang_mills_t0_SOURCES=./src/yang_mills_t0.c \
                       $(dSFMT)  \
                       $(endian) \
                       $(funcpoint) \
                       $(gauge_conf) \
                       $(geometry) \
                       $(gparam) \
                       ./include/macro.h \
                       $(random) \
                       $(su2) \
                       $(sun) \
                       $(tens_prod) \
                       $(u1) 

yang_mills_tube_conn_SOURCES=./src/yang_mills_tube_conn.c \
                             $(dSFMT)  \
                             $(endian) \
                             $(funcpoint) \
                             $(gauge_conf) \
                             $(geometry) \
                             $(gparam) \
                             ./include/macro.h \
                             $(random) \
                             $(su2) \
                             $(sun) \
                             $(tens_prod) \
                             $(u1) 


yang_mills_tube_disc_SOURCES=./src/yang_mills_tube_disc.c \
                             $(dSFMT)  \
                             $(endian) \
                             $(funcpoint) \
                             $(gauge_conf) \
                             $(geometry) \
                             $(gparam) \
                             ./include/macro.h \
                             $(random) \
                             $(su2) \
                             $(sun) \
                             $(tens_prod) \
                             $(u1) 

yang_mills_tracedef_SOURCES=./src/yang_mills_tracedef.c \
                            $(dSFMT)  \
                            $(endian) \
                            $(funcpoint) \
                            $(gauge_conf) \
                            $(geometry) \
                            $(gparam) \
                            ./include/macro.h \
                            $(random) \
                            $(su2) \
                            $(sun) \
                            $(tens_prod) \
                            $(u1) 

AM_CFLAGS = -std=c99 -O2 -Wall -Wextra -Werror -pedantic -Wconversion @OPENMP_CFLAGS@

CLEANFILES = dati.dat log.dat dati.log conf.dat* ml.dat template_input.in 
