include ../Makefile.inc

CLEANFILES = *.$(ObjSuf) *.$(DllSuf) *Dict.* *~ core
noinst_HEADERS = TPCSubscriberLinkDef.h
pkginclude_HEADERS = SubscriberQA.h
BUILT_SOURCES = TPCSubscriberDict.h TPCSubscriberDict.cxx

bin_PROGRAMS = libAmoreTPCSubscriber.$(DllSuf)
libAmoreTPCSubscriber_CPPFLAGS = $(CPPFLAGS) -I../common
libAmoreTPCSubscriber_CXXFLAGS = $(CXXFLAGS)
libAmoreTPCSubscriber_LIBS = $(LIBS) -lAmoreTPCCommon -lProof
libAmoreTPCSubscriber_LDFLAGS = $(LDFLAGS) -L../common -L$(AMORE_SITE)/lib -shared -Wl,-soname,$(bin_PROGRAMS) -o $(bin_PROGRAMS)
libAmoreTPCSubscriber_SOURCES = TPCSubscriberDict.cxx SubscriberQA.cxx
libAmoreTPCSubscriber_OBJECTS = TPCSubscriberDict.$(ObjSuf) $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(wildcard *.$(SrcSuf)))

%.$(ObjSuf):	%.$(SrcSuf) %.$(HdrSuf)
	@echo $<
	@$(CXX) $(INCLUDES) $(libAmoreTPCSubscriber_CPPFLAGS) $(libAmoreTPCSubscriber_CXXFLAGS) -c $<

default: all

all: $(bin_PROGRAMS)

$(bin_PROGRAMS): $(libAmoreTPCSubscriber_OBJECTS)
	@$(LD) $(libAmoreTPCSubscriber_LDFLAGS) $(libAmoreTPCSubscriber_LIBS) $^
#$(OutPutOpt)$@

install: all
	@cp $(bin_PROGRAMS) ${AMORE_SITE}/lib;
	
uninstall:
	@rm -rf ${AMORE_SITE}/lib/$(bin_PROGRAMS);

clean:
	@rm -rf $(CLEANFILES)

%Dict.$(SrcSuf) %Dict.$(HdrSuf): $(pkginclude_HEADERS)
	@echo Creating dictionary...
	@rootcint -f $*Dict.$(SrcSuf) -c -p $(INCLUDES) $(libAmoreTPCSubscriber_CPPFLAGS) $(libAmoreTPCSubscriber_CXXFLAGS) $(pkginclude_HEADERS) $*LinkDef.h
