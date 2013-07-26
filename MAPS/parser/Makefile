CFLAGS = -g -Wall
CPPFLAGS = -I/usr/include/ -I../brigen/ -I../visgen/ -I../utilities/  \
	-I${SIM_INC}
LDLIBS = -lm

CC = gcc
YACC = bison
LEX = flex

.PHONY: install
install: lex.yy.o read_prj.tab.o

lex.yy.o read_prj.tab.o : lex.yy.c read_prj.tab.c

lex.yy.c read_prj.tab.c read_prj.tab.h: read_prj.y read_prj.l
	$(YACC) -d -v read_prj.y
	$(LEX) read_prj.l

clean:
	rm -f *.o lex.yy.c read_prj.tab.c read_prj.tab.h read_prj.output

