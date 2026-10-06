CC = cc
PREFIX ?= /opt/homebrew
CPPFLAGS = -I$(PREFIX)/include
CFLAGS = -O2 -fcommon -Wall
LDLIBS = -L$(PREFIX)/lib -lfftw3 -lgsl -lgslcblas -lm
SOURCES = spino5.c get_input_spino2.c init_confv2.c evolve12.c eigen_strain.c calc_bn.c out_conf.c gsl_support.c
binorder: $(SOURCES) binary.h gsl_support.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) $(LDLIBS) -o $@
clean:
	rm -f binorder
