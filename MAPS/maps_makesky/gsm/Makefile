GSM_FFLAGS := -ffixed-line-length-132 -Wall
GSM_CFLAGS =-g -Wall $(CFLAGS) -I$(SIM_INC) -D_FILE_OFFSET_BITS=64 -I$(CFITSIO_INC) -I$(WCS_INC)   
GSM_LFLAGS = $(LDFLAGS) -L$(WCS_LIB) -L$(SIM_LIB) -L$(CFITSIO_LIB)

GSM_CPPFLAGS := $(GSM_CFLAGS)
GSM_CFLAGS := -std=gnu99 $(GSM_CFLAGS)

GSM_LIBRARY = libutils.a

GSM_LIB_OBJ = get_pix_area.o \
	  convert_to_Iv.o \
	  convert_to_Jy.o \
	  hpxcarre.o \
	  dat2carre.o \
	  hpxwrap.o \
	  npix2nside.o \
	  pix2ang_ring.o \
	  wcs_utils.o \
	  simple_poly.o

ifeq ($(HAVE_CBLAS),TRUE)	  
  GSM_CPPFLAGS := $(GSM_CPPFLAGS) -I$(CBLAS_INC)
  GSM_PROGS:=diffuse_sky convert_hpx convert_dat HPXcvt wmap_pol trace_pol build_pol.pl \
  		fix_sun
else
  GSM_PROGS:=diffuse_sky convert_hpx convert_dat HPXcvt wmap_pol build_pol.pl
endif

all: $(GSM_LIBRARY) $(GSM_PROGS) 
	
$(GSM_LIBRARY): $(GSM_LIB_OBJ)
	$(AR) $(ARFLAGS) $(GSM_LIBRARY) $(GSM_LIB_OBJ)


.SUFFIXES:
.SUFFIXES:.o .c .cpp

%.o:%.c
	$(CC) $(GSM_CFLAGS) $(GSM_INC) -c $<
%.o:%.cpp
	g++ $(GSM_CPPFLAGS) $(INC) -c $<
convert_dat: convert_dat.o $(LIBRARY)
	g++ $(GSM_CPPFLAGS) -o convert_dat convert_dat.o $(GSM_LFLAGS) $(GSM_LIBRARY) -lcfitsio $(GSM_LIBRARY)
convert_hpx: convert_hpx.o $(GSM_LIBRARY)
	g++ $(GSM_CPPFLAGS) -o convert_hpx convert_hpx.o $(GSM_LFLAGS) -L$(WCS_LIB) $(GSM_LIBRARY) -lcfitsio -lwcs -lnovas_c $(PGPLOT_LDFLAGS)
HPXcvt: HPXcvt.o
	g++ $(GSM_CPPFLAGS) -o HPXcvt HPXcvt.o  $(GSM_LFLAGS) -lcfitsio 
diffuse_sky: gsm.f
	gfortran $(GSM_FFLAGS) gsm.f -o diffuse_sky
wmap_pol: wmap_pol.o
	g++ $(GSM_CPPFLAGS) -o wmap_pol wmap_pol.o $(GSM_LFLAGS) -lcfitsio
trace_pol: trace_pol.o
	g++ $(GSM_CPPFLAGS) -o trace_pol trace_pol.o $(GSM_LFLAGS) -L$(CBLAS_LIB) -lcfitsio -lcblas -lm

fix_sun: fix_sun.o
	g++ $(GSM_CPPFLAGS) -o fix_sun fix_sun.o $(GSM_LFLAGS) -L$(CBLAS_LIB) -lcfitsio -lcblas -lm
clean:
	rm -f *.o
install:
	cp $(GSM_PROGS) $(SIM_BIN)
