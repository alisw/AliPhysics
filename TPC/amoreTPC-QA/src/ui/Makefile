include ../Makefile.inc

CLEANFILES = *.$(ObjSuf) *.$(DllSuf) *Dict.* *~ core
noinst_HEADERS = TPCUILinkDef.h
pkginclude_HEADERS = UIQA.h
BUILT_SOURCES = TPCUIDict.h TPCUIDict.cxx

bin_PROGRAMS = libAmoreTPCUI.$(DllSuf)
libAmoreTPCUI_CPPFLAGS = $(CPPFLAGS) -I../common
libAmoreTPCUI_CXXFLAGS = $(CXXFLAGS)
libAmoreTPCUI_LIBS = $(LIBS) -lAmoreTPCCommon
libAmoreTPCUI_LDFLAGS = $(LDFLAGS) -L../common -L$(AMORE_SITE)/lib -shared -Wl,-soname,$(bin_PROGRAMS) -o $(bin_PROGRAMS)
libAmoreTPCUI_SOURCES = TPCUIDict.cxx UIQA.cxx
libAmoreTPCUI_OBJECTS = TPCUIDict.$(ObjSuf) $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(wildcard *.$(SrcSuf)))

%.$(ObjSuf):	%.$(SrcSuf) %.$(HdrSuf)
	@echo $<
	@$(CXX) $(INCLUDES) $(libAmoreTPCUI_CPPFLAGS) $(libAmoreTPCUI_CXXFLAGS) -c $<

default: all

all: $(bin_PROGRAMS)

$(bin_PROGRAMS): $(libAmoreTPCUI_OBJECTS)
	@$(LD) $(libAmoreTPCUI_LDFLAGS) $(libAmoreTPCUI_LIBS) $^
#$(OutPutOpt)$@

install: all
	@cp $(bin_PROGRAMS) ${AMORE_SITE}/lib;
	
uninstall:
	@rm -rf ${AMORE_SITE}/lib/$(bin_PROGRAMS);

clean:
	@rm -rf $(CLEANFILES)

%Dict.$(SrcSuf) %Dict.$(HdrSuf): $(pkginclude_HEADERS)
	@echo Creating dictionary...
	@rootcint -f $*Dict.$(SrcSuf) -c -p $(INCLUDES) $(libAmoreTPCUI_CPPFLAGS) $(libAmoreTPCUI_CXXFLAGS) $(pkginclude_HEADERS) $*LinkDef.h
