# $Header$

ifeq (${REVESYS},)
REVESYS := ${CURDIR}/..
endif
MKDEPINC := make_depend.inc

LD       = ${CXX}
LDEXT    = so

MAKEDEPEND = ${ROOTSYS}/bin/rmkdepend -Y -w 666 -f-
ROOTCINT   = ${ROOTSYS}/bin/rootcint

CPPFLAGS = -I. -I${REVESYS} -I${ROOTSYS}/include
CXXFLAGS = -g -fPIC -Wall -Woverloaded-virtual -Wno-long-long -pedantic-errors -ansi
LDFLAGS  = -shared

ROOTCINTFLAGS = -I. -I${REVESYS}

.SUFFIXES:  .cxx .cc .h .o
.SECONDARY: ${ALLOBJS} 

.SILENT: ${ALLOBJS} 

%.o: %.cc
	@echo --- Compiling $<
	${CXX} ${CPPFLAGS} -c -o $@ ${CXXFLAGS} $<

%.o: %.cxx
	@echo --- Compiling $<
	${CXX} ${CPPFLAGS} -c -o $@ ${CXXFLAGS} $<

lib%.so: ${OBJS} %_Dict.o
	${LD} ${LDFLAGS} -o $@ $^ ${ROOTLIBS} ${EXTLIBS}

#--------------------------------------------------------------

all: ${LIBS} ${EXES} .inst_libs .inst_exes

#--------------------------------------------------------------

.inst_libs: ${LIBS}
	@mkdir -p ${REVESYS}/lib
	@if [ -n '${LIBS}' ]; then rsync -a ${LIBS} ${REVESYS}/lib; fi
	@touch .inst_libs

.inst_exes: ${EXES}
	@mkdir -p ${REVESYS}/bin
	@if [ -n '${EXES}' ]; then rsync -a ${EXES} ${REVESYS}/bin; fi
	@touch .inst_exes

#--------------------------------------------------------------

depend:
	@rm -f ${MKDEPINC}

${MKDEPINC}:
	@echo "Creating dependencies for module=${MODULE}, name=${NAME}."
	@${MAKEDEPEND} -- -- ${SRCS} > ${MKDEPINC} 2>/dev/null

clean:
	@rm -f ${ALLOBJS} *_Dict.* ${LIBS} ${EXES}

distclean: clean
	@rm -f ${MKDEPINC} .inst_*
	@(cd ${REVESYS}/lib;    rm -f $(notdir ${LIBS}))
	@(cd ${REVESYS}/bin;    rm -f $(notdir ${EXES}))

#--------------------------------------------------------------

echo_module:
	@echo ${MODULE}

echo_hdrs:
	@echo ${HDRS}

echo_srcs:
	@echo ${SRCS}

echo_objs:
	@echo ${OBJS}

echo_libs:
	@echo ${LIBS}

echo_exes:
	@echo ${EXES}

#--------------------------------------------------------------

-include ${MKDEPINC}
