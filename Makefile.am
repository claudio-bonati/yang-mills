bin_PROGRAMS = debug_rng \
               debug_su2 \
               debug_sun \
               yang_mills

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
             ./lib/gauge_conf_upd.c 

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

mymalloc = ./include/mymalloc.h \
           ./lib/mymalloc.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c 

su2 = ./include/su2.h \
      ./include/su2_upd.h \
      ./lib/su2.c \
      ./lib/su2_upd.c 

sun = ./include/sun_aux.h \
      ./include/sun.h \
      ./include/sun_upd.h \
      ./lib/sun_aux.c \
      ./lib/sun.c \
      ./lib/sun_upd.c 


if ONE_FILE_COMPILATION
  yang_mills_SOURCES=./src/yang_mills.c \
                     $(dSFMT)
  EXTRA_DIST = ./src/input.in \
               $(endian) \
               $(funcpoint) \
               $(gauge_conf) \
               $(geometry) \
               $(gparam) \
               ./include/macro.h \
               $(mymalloc) \
               $(random) \
               $(su2) \
               $(sun) 

  DISTCHECK_CONFIGURE_FLAGS = --enable-one-file-compilation
else
  yang_mills_SOURCES=./src/yang_mills.c \
                     $(dSFMT)  \
                     $(endian) \
                     $(funcpoint) \
                     $(gauge_conf) \
                     $(geometry) \
                     $(gparam) \
                     ./include/macro.h \
                     $(mymalloc) \
                     $(random) \
                     $(su2) \
                     $(sun) 
  EXTRA_DIST = ./src/input.in 
endif

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

AM_CFLAGS = -std=c99 -O3 -Wall -Wextra -pedantic -Werror -Wconversion @OPENMP_CFLAGS@

CLEANFILES = dati*.dat log*.dat dati*.log conf*.dat* input_check.in 
