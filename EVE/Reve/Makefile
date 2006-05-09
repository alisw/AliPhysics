# $Header$

# MODULE ~include sub-directory, namespace
# NAME   ~ name of the library

MODULE   = Reve
NAME     = Reve
LIBNAME  = lib${NAME}.${LDEXT}
DICTSTEM = ${NAME}_Dict

HDRS    = $(filter-out %LinkDef.h %_Dict.h, $(wildcard *.h))
SRCS    = $(wildcard *.cxx)
OBJS    = $(SRCS:.cxx=.o)
LIBOBJS = ${OBJS} ${DICTSTEM}.o
ALLOBJS = ${LIBOBJS} reve_main/reve_main.o

LIBS    = ${LIBNAME}
EXES    = reve

ROOTLIBS  = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

include ../make_base.inc

${LIBNAME}: ${LIBOBJS}
	${LD} ${LDFLAGS} -o $@ $^

${DICTSTEM}.cc: ${HDRS} LinkDef.h
	${ROOTCINT} -f $@ -c ${ROOTCINTFLAGS} $^

reve: reve_main/reve_main.o ${LIBNAME}
	${CXX} ${CPPFLAGS} -o $@ ${CXXFLAGS} $< -L. -lReve \
	   ${ROOTGLIBS} -lTreePlayer -lEG -lGeom -lGeomPainter -lGed -lRGL
