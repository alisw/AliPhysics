include ../Makefile.inc

CLEANFILES = *.$(ObjSuf) *.$(DllSuf) *Dict.* *~ core
noinst_HEADERS = TPCCommonLinkDef.h
pkginclude_HEADERS = Common.h
BUILT_SOURCES = TPCCommonDict.h TPCCommonDict.cxx

bin_PROGRAMS = libAmoreTPCCommon.$(DllSuf)
libAmoreTPCCommon_CPPFLAGS = $(CPPFLAGS)
libAmoreTPCCommon_CXXFLAGS = $(CXXFLAGS)
libAmoreTPCCommon_LIBS = $(LIBS)
libAmoreTPCCommon_LDFLAGS = $(LDFLAGS) -shared -Wl,-soname,$(bin_PROGRAMS) -o $(bin_PROGRAMS)
libAmoreTPCCommon_SOURCES = TPCCommonDict.cxx Common.cxx
libAmoreTPCCommon_OBJECTS = TPCCommonDict.$(ObjSuf) $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(wildcard *.$(SrcSuf)))

%.$(ObjSuf):	%.$(SrcSuf) %.$(HdrSuf)
	@echo $<
	@$(CXX) $(INCLUDES) $(libAmoreTPCCommon_CPPFLAGS) $(libAmoreTPCCommon_CXXFLAGS) -c $<

default: all

all: $(bin_PROGRAMS)

$(bin_PROGRAMS): $(libAmoreTPCCommon_OBJECTS)
	@$(LD) $(libAmoreTPCCommon_LDFLAGS) $(libAmoreTPCCommon_LIBS) $^
#$(OutPutOpt)$@

install: all
	@cp $(bin_PROGRAMS) ${AMORE_SITE}/lib;
	
uninstall:
	@rm -rf ${AMORE_SITE}/lib/$(bin_PROGRAMS);

clean:
	@rm -rf $(CLEANFILES)
	
%Dict.$(SrcSuf) %Dict.$(HdrSuf): $(pkginclude_HEADERS)
	@echo Creating dictionary...
	@rootcint -f $*Dict.$(SrcSuf) -c -p $(INCLUDES) $(libAmoreTPCCommon_CPPFLAGS) $(libAmoreTPCCommon_CXXFLAGS) $(pkginclude_HEADERS) $*LinkDef.h
