#!/bin/bash
#
# Top level Makefile for MAPS                   DAM Sep 27 2007
#

BASEDIR = $(PWD)
SUBDIRS = parser novas-c201 utilities LOsim/source maps_makesky maps_im2uv \
	visgen maps2uvfits brigen
BINEXECUTABLES = convert_hpx diffuse_sky LOsim maps2uvfits MAPS_im2uv \
	maps_makesky visgen brigen
GSMEXECUTABLES = HPXcvt wmap_pol convert_dat

install:
	@for subdir in $(SUBDIRS); do \
	  echo "Installing all in $$subdir"; \
	  cd $$subdir && $(MAKE) install && cd $(BASEDIR); \
	done

clean:
	@for subdir in $(SUBDIRS); do \
	  echo "Cleaning object files in $$subdir"; \
	  cd $$subdir && $(MAKE) clean && cd $(BASEDIR); \
	done

purge:
	@echo 'Cleaning remporary emacs *~ and vi #*# files'; \
	rm -f `find ./ -name '*~'` && rm -f `find ./ -name '#*#'`; \
	echo 'Cleaning executables in /bin'; \
	cd $(BASEDIR)/bin && rm -f $(BINEXECUTABLES)  && cd $(BASEDIR); \
	echo 'Cleaning libraries in lib'; \
	cd $(BASEDIR)/lib && rm -f *.a  && cd $(BASEDIR); \
	echo 'Cleaning executables and libraries in maps_makesky/gsm'; \
	cd $(BASEDIR)/maps_makesky/gsm && rm -f *.a; \
	    rm -f $(GSMEXECUTABLES) &&  cd $(BASEDIR); \

unversion:
	@echo 'Remove all the Subversion (.svn) direcrories recursively';
	rm -rf `find . -type d -name .svn`
