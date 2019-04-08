SHELL	=	/bin/sh
#
# Local definitions
#
#C      =       cc4 -O3 -ffast-math -mfpmath=387 -mieee-fp -ftree-vectorize
#C      =       cc4 -O3 -ffast-math -mfpmath=387 -mieee-fp -funroll-all-loops -ftree-vectorize
#C	=	mpicc4 -O3 
#C	=	icc -O2
#C	=	mpicc -O2 -mcmodel=large  -shared-intel -L/opt/intel/fce/10.0.023/lib
C	=	mpicc -O2 -mcmodel=large  -shared-intel -L/opt/intel/fce/15.0.090/lib
#C	=	mpicxx -g -O3
#C	=	icc -O2 -mcmodel=medium -shared-intel
#
# TARGET RULES
#
#
.f.o:;   
#
FLAGS  = -mpitrace
FLAGS  =
#path for source files
SRCDIR = ./src
#path for source files
OBJDIR = ./obj
#path for header files	
VPATH = ./src		

meam_SYSLIB = -lmeam -lifcore -lsvml -limf
meam_SYSPATH = -L/opt/intel/fce/15.0.090/lib
#
#LIB= -llammps -lmeam -lfftw -lmpich -lmpl -lpthread -lstdc++ -lgsl -lgslcblas -lm ${meam_SYSLIB}
LIB= -llammps -lpthread -lstdc++ -lgsl -lgslcblas -lm ${meam_SYSLIB}
INC =  -I/home-3/pyi7@jhu.edu/usr/include
HEADERS = distributions.h ensembles.h forcefield.h history.h init.h io.h lists.h motion.h position.h random.h \
rebridge.h roots.h sample.h units.h varbridge.h vector.h input.h header.h globals.h types.h

OBJECTS = distributions.o ensembles.o forcefield.o history.o init.o io.o lists.o motion.o position.o random.o \
rebridge.o roots.o sample.o units.o varbridge.o vector.o
	
#

#

main: $(OBJECTS) main.o
	$(C) -o main  $(FLAGS) $(LIB) $(OBJECTS) main.o;
#	mv main ~/bin/;
#	rm $(OBJECTS) main.o;

block:	$(OBJECTS) block.o
	$(C) -o block $(FLAGS) $(LIB) block.o $(OBJECTS);
	rm $(OBJECTS) block.o

builder: $(OBJECTS) builder.o
	$(C) -o builder $(FLAGS) $(LIB) builder.o $(OBJECTS);
	rm $(OBJECTS) builder.o

conf2car: $(OBJECTS) conf2car.o
	$(C) -o conf2car $(FLAGS) $(LIB) conf2car.o $(OBJECTS);
	mv conf2car ~/bin/;
#	rm $(OBJECTS) conf2car.o

conf2lammps: $(OBJECTS) conf2lammps.o
	$(C) -o conf2lammps $(FLAGS) $(LIB) conf2lammps.o $(OBJECTS);
	rm $(OBJECTS) conf2lammps.o

conf2DL_POLY: $(OBJECTS) conf2DL_POLY.o
	$(C) -o conf2DL_POLY $(FLAGS) $(LIB) conf2DL_POLY.o $(OBJECTS);
	rm $(OBJECTS) conf2DL_POLY.o

conf2gro: $(OBJECTS) conf2gro.o
	$(C) -o conf2gro $(FLAGS) $(LIB) conf2gro.o $(OBJECTS);
	rm $(OBJECTS) conf2gro.o

config: $(OBJECTS) config.o
	$(C) -o config $(FLAGS) $(LIB) config.o $(OBJECTS);
	rm $(OBJECTS) config.o;

dumphst: $(OBJECTS) dumphst.o
	$(C) -o dumphst $(FLAGS) $(LIB) dumphst.o $(OBJECTS);
	rm $(OBJECTS) dumphst.o;

embed:	$(OBJECTS) embed.o
	$(C) -o embed $(FLAGS) $(LIB) embed.o $(OBJECTS);
	rm $(OBJECTS) embed.o

hist:	$(OBJECTS) histogram.o
	$(C) -o hist $(FLAGS) $(LIB) histogram.o $(OBJECTS);
	rm $(OBJECTS) histogram.o;

hst:	$(OBJECTS) hst.o
	$(C) -o hst $(FLAGS) $(LIB) hst.o $(OBJECTS);
	rm $(OBJECTS) hst.o

hstcomb:	$(OBJECTS) hstcomb.o
	$(C) -o hstcomb $(FLAGS) $(LIB) hstcomb.o $(OBJECTS);
	rm $(OBJECTS) hstcomb.o

lammps2car: $(OBJECTS) lammps2car.o
	$(C) -o lammps2car $(FLAGS) $(LIB) lammps2car.o $(OBJECTS);
	rm $(OBJECTS) lammps2car.o

lammps2hst: $(OBJECTS) lammps2hst.o
	$(C) -o lammps2hst $(FLAGS) $(LIB) lammps2hst.o $(OBJECTS);
	#mv lammps2hst ~/bin/lammps2hst;
	#mv lammps2hst ~/bin/MPIlammps2hst;
	mv lammps2hst MPIlammps2hst;
#	rm $(OBJECTS) lammps2hst.o

lammpstest: $(OBJECTS) lammpstest.o
	$(C) -o lammpstest $(FLAGS) $(LIB) lammpstest.o $(OBJECTS);
	rm $(OBJECTS) lammpstest.o

lmp_conf: $(OBJECTS) lmp_conf.o
	$(C) -o lmp_conf $(FLAGS) $(LIB) lmp_conf.o $(OBJECTS);
	rm $(OBJECTS) lmp_conf.o

metal2hst: $(OBJECTS) metal2hst.o
	$(C) -o metal2hst $(FLAGS) metal2hst.o $(OBJECTS) $(LIB);

extractlammps: $(OBJECTS) extractlammps.o
	$(C) -o extractlammps $(FLAGS) $(LIB) extractlammps.o $(OBJECTS);
	rm $(OBJECTS) extractlammps.o

test: test.o
	$(C) -o testprog $(FLAGS) $(LIB) test.o;
	rm test.o
test.o: test.c $(HEADERS)
	$(C) -c $(INC) test.c
#
#
main.o: main.c $(HEADERS)
	$(C) -c $(SRCDIR)/main.c
builder.o: builder.c $(HEADERS)
	$(C) -c $(SRCDIR)/builder.c
block.o:  block.c $(HEADERS)
	$(C) -c $(SRCDIR)/block.c
conf2car.o: conf2car.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2car.c
conf2lammps.o: conf2lammps.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2lammps.c
conf2DL_POLY.o: conf2DL_POLY.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2DL_POLY.c
conf2gro.o: conf2gro.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2gro.c
config.o: config.c $(HEADERS)
	$(C) -c $(SRCDIR)/config.c
dumphst.o: dumphst.c  $(HEADERS)
	$(C) -c $(SRCDIR)/dumphst.c
embed.o: embed.c $(HEADERS)
	$(C) -c $(SRCDIR)/embed.c
histogram.o:	histogram.c $(HEADERS)
	$(C) -c $(SRCDIR)/histogram.c
hst.o:	hst.c $(HEADERS)
	$(C) -c $(SRCDIR)/hst.c
hstcomb.o:	hstcomb.c $(HEADERS)
	$(C) -c $(SRCDIR)/hstcomb.c
lammps2car.o: lammps2car.c $(HEADERS)
	$(C) -c $(SRCDIR)/lammps2car.c
lammps2hst.o: lammps2hst.c $(HEADERS)
	$(C) -c $(SRCDIR)/lammps2hst.c
lmp_conf.o: lmp_conf.c $(HEADERS)
	$(C) -c $(SRCDIR)/lmp_conf.c
metal2hst.o: metal2hst.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/metal2hst.c
extractlammps.o: extractlammps.c $(HEADERS)
	$(C) -c $(SRCDIR)/extractlammps.c


distributions.o: distributions.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/distributions.c
ensembles.o: ensembles.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/ensembles.c
forcefield.o: forcefield.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/forcefield.c
history.o: history.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/history.c
init.o:	 init.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/init.c
io.o:	io.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/io.c
lists.o: lists.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/lists.c
motion.o: motion.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/motion.c
position.o: position.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/position.c
random.o:  random.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/random.c
rebridge.o: rebridge.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/rebridge.c
roots.o: roots.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/roots.c
sample.o: sample.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/sample.c
units.o:  units.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/units.c
varbridge.o: varbridge.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/varbridge.c
vector.o: vector.c $(HEADERS)
	$(C) -c $(INC) $(SRCDIR)/vector.c

clean:
	rm *.o *~ output

