#
# Makefile for brigen
# L. Benkevitch Jan 16 2011
# Modified 15 Feb 2011 (thanks to Jeoffrey B. Crew)
#
CFLAGS = -g -Wall
CPPFLAGS = -I/usr/include/ -I../brigen/ -I../parser/ -I../visgen/  \
	   -I../utilities/ -I${SIM_INC}
LDLIBS = -lcfitsio -lm

CC = gcc
YACC = bison
LEX = flex

PARSER = ${SIM}/parser

brigen: brigen.o sources_to_image.o add_extimages.o image_to_fits.o \
	$(PARSER)/read_prj.tab.o $(PARSER)/lex.yy.o 

$(PARSER)/read_prj.tab.o $(PARSER)/lex.yy.o:
	echo "Compiling the parser"; \
	pushd $(PARSER) && $(MAKE) install && popd

#read_prj.tab.c read_prj.tab.h: read_prj.y
#	$(YACC) -d read_prj.y

#lex.yy.c: read_prj.l
#	$(LEX) read_prj.l

install: brigen
	 mv brigen $(SIM_BIN)

clean:
	rm -f *.o 
#	rm -f *.o lex.yy.c read_prj.tab.c read_prj.tab.h read_prj.output
