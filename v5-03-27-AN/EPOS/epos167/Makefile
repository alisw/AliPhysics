SOURCE_DIR=
LIBDIR=$(OBJ)/$(shell uname)/
OUTDIR=outputs/
VER=165
FC=$(shell root-config --f77)
CPP=$(shell root-config --cxx)
LFLAGS= 
FFLAGS= 
UNAME=$(shell uname)
CFLAGS=$(shell root-config --cflags)
LIBS=$(shell root-config --libs) -lg2c
ifeq ($(UNAME),Linux)
  FC=$(shell root-config --f77)
##  FC=gfortran
  CPP=$(shell root-config --cxx)
##   ---     -pg for using gprof
  SFFLAGS= -g -Wall -fno-automatic -finit-local-zeros -fbounds-check
  SLFLAGS= -g 
  URQMD=
##   URQMD= -g ../urqmd23/obj/Linux/urqmd/*.o 	   
endif
ifeq ($(UNAME),SunOS)
  SFFLAGS= 
  SLFLAGS=
endif
ifeq ($(UNAME),AIX)
  SFFLAGS= 
  SLFLAGS=
endif

ifeq ($(UNAME),HP-UX)
  SFFLAGS= 
  SLFLAGS=+FPVZO
endif

FILES=  \
epos-app-158.f	\
epos-bas-168.f  epos-con-161.f                 \
epos-dky-165.f  epos-dro-168.f epos-ems-165.f  \
epos-fra-163.f  epos-ico-158.f epos-yyy-168.f  \
epos-ids-165.f  epos-int-168.f epos-jps-164.f  \
epos-omg-160.f  epos-par-129.f epos-sem-165.f  \
epos-rsh-165.f  epos-qsh-165.f epos-tim-155.f  \
epos-uti-165.f  epos-xan-154.f epos-xpr-165.f  \
                                           \
eposm.f                             \
                                     \
eposu.f      \
#eposurqmd.f                          \
 
CFILES=epostrapfpe.c epostimer.c main.c

IFILES=epos.inc epos.incems epos.incico epos.incpar epos.incsem epos.param \
	epos.inirj epos.inidi epos.ini1b epos.iniev epos.inirj epos.initl \
	Makefile epos.inics .epos epos cmake isigma.optns epos_all add.pl \
	iems.optns models.f

HEADERS=eposproc.h

X: dirs Xepos
	echo $(SOURCE_DIR)
Xepos:  $(CFILES:%.c=$(LIBDIR)%.o) $(FILES:%.f=$(LIBDIR)%.o)
	$(CPP) $(SLFLAGS) $(LFLAGS) $(URQMD) $(LIBS) -o $(LIBDIR)$@ $^  
$(FILES:%.f=$(LIBDIR)%.o) : $(LIBDIR)%.o : $(SOURCE_DIR)%.f  
	$(FC) $(SFFLAGS) $(FFLAGS) -o $@ -c $<
$(CFILES:%.c=$(LIBDIR)%.o) : $(LIBDIR)%.o : $(SOURCE_DIR)%.c $(HEADERS)
	$(CPP) $(CFLAGS) -o $@ -c $< 
tar: $(FILES) $(CFILES) $(IFILES)
	tar czvf epos$(VER).tgz $^
cors: $(FILES) epos.inc epos.incems epos.param epos.inirj epos.inics
	scp $^ $(HOME)/corsika/corsika.svn/epos/
copy: $(FILES) epos.inc epos.incems epos.param Makefile
	scp $^ iklx62.fzk.de:nexus/epos/epos121
cpcc: $(FILES) epos.inc epos.param epos.inirj epos.inics
	scp $^ ccali.in2p3.fr:epos115
dirs:
	@if [ ! -d $(LIBDIR)  ]  ;then                       \
	set -x; mkdir -p $(LIBDIR); set +x;                  \
	fi     
clean:
	/bin/rm -r $(LIBDIR)
