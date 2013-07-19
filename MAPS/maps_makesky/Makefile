FFLAGS := -ffixed-line-length-132 -Wall
INCS=-I$(SIM_INC) 
CFLAGS=-g -Wall $(INCS) -D_FILE_OFFSET_BITS=64

all: maps_makesky diffuse_sky construct_sky

maps_makesky: maps_makesky.c
	cc $(CFLAGS) -o maps_makesky maps_makesky.c -I$(CFITSIO_INC) -I$(WCS_INC) -L$(CFITSIO_LIB) -L$(WCS_LIB) -lwcs -lcfitsio -lm

diffuse_sky: gsm/gsm.f
	cd gsm;$(MAKE) all; cd ../ 

construct_sky: scripts/construct_sky.pl
	cd scripts;$(MAKE) all; cd ../

install: maps_makesky diffuse_sky construct_sky
	mv maps_makesky $(SIM_BIN)
	mv gsm/diffuse_sky $(SIM_BIN)	
	mv gsm/convert_hpx $(SIM_BIN)
	cp scripts/*.pl $(SIM_BIN)

clean:
	rm -f *.o gsm/*.o

test_wcs: test_wcs.c
	cc $(CFLAGS) -o test_wcs test_wcs.c -I$(CFITSIO_INC) -I/usr/local/include/wcslib -L$(CFITSIO_LIB) -lcfitsio -lwcs -lm

