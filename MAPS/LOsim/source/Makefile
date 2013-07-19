# Makefile for LOsim	RB, 25 Nov 2003. Updated RBW. Jan 2007

# Modify the environments below for your system
CFLAGS= -D_FILE_OFFSET_BITS=64 -g -O
LOSIM_INC = $(SIM)/LOsim/include

.c.o:
	gcc -c -I$(LOSIM_INC) $(CFLAGS) $<

SOURCES = LOsim.c SumSources.c SourceListRead.c ScreenIO.c \
          RowTransformBrightness.c ColumnTransformBrightness.c OOBSources.c \
          GlobalVariables.c DescriptionFileRead.c \
          ObservationRead.c StationListRead.c ImagesRead.c ImageAdd.c \
          sidereal_time.c julian_date.c \
          RefreshRow.c \
          SinesAndTangents.c BitReverse.c Butterfly.c \
          GetColumnBlock.c LargeFiles.c  SimFilename.c Utilities.c \
          DFT.c FITSread.c FITSwrite.c ComplexArithmetic.c

OBJECTS = LOsim.o SumSources.o SourceListRead.o ScreenIO.o \
          RowTransformBrightness.o ColumnTransformBrightness.o OOBSources.o \
          GlobalVariables.o DescriptionFileRead.o \
          ObservationRead.o StationListRead.o ImagesRead.o ImageAdd.o \
          sidereal_time.o julian_date.o \
          RefreshRow.o \
          SinesAndTangents.o BitReverse.o Butterfly.o \
          GetColumnBlock.o LargeFiles.o  SimFilename.o Utilities.o \
          DFT.o FITSread.o FITSwrite.o ComplexArithmetic.o


OBJECTS_PATCHTEST = open_gridfile.o get_patch.o \
                    GlobalVariables.o LargeFiles.o SimFilename.o DescriptionFileRead.o Utilities.o \
                    ObservationRead.o

LIBRARIES = -lm
.PHONY: all
all: LOsim 
# Vistest

LOsim: $(OBJECTS)
#	gcc -static -o LOsim $(OBJECTS) $(LIBRARIES) # link command creates executable LOsim
	gcc $(CFLAGS) -o LOsim $(OBJECTS) $(LIBRARIES) # link command creates executable LOsim

install: LOsim
	mv LOsim $(SIM_BIN)

clean:
	rm -f  $(OBJECTS) $(OBJECTS_VV) $(OBJECTS_PATCHTEST)

