INCS = -I$(CFITSIO_INC)
LIBS = -L$(CFITSIO_LIB)
CFLAGS= -Wall -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE

MAPS_im2uv: MAPS_im2uv.c
#	cc $(CFLAGS) -o MAPS_im2uv MAPS_im2uv.c $(INCS) $(LIBS) -lfftw3_threads -lfftw3 -lcfitsio -lm -lpthread 
	cc $(CFLAGS) -o MAPS_im2uv MAPS_im2uv.c $(INCS) $(LIBS) -lfftw3 -lcfitsio -lm

install: MAPS_im2uv
	mv MAPS_im2uv $(SIM_BIN)

clean:
	rm -f *.o MAPS_im2uv
