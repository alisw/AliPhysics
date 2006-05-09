# $Header$

# MODULE ~include sub-directory, namespace
# NAME   ~ name of the library

MODULE   = Alieve
NAME     = Alieve
LIBNAME  = lib${NAME}.${LDEXT}
DICTSTEM = ${NAME}_Dict

HDRS = $(filter-out %LinkDef.h %_Dict.h, $(wildcard *.h))
SRCS = $(wildcard *.cxx)
OBJS = $(SRCS:.cxx=.o)
ALLOBJS = ${OBJS} ${DICTSTEM}.o

LIBS    = ${LIBNAME}

include ../make_base.inc

ALICE_INCLUDES = -I${ALICE_ROOT}/include -I${ALICE_ROOT}/TPC -I${ALICE_ROOT}/ITS

CPPFLAGS += ${ALICE_INCLUDES}

${DICTSTEM}.cc: ${HDRS} LinkDef.h
	${ROOTCINT} -f $@ -c ${ROOTCINTFLAGS} ${ALICE_INCLUDES} $^
