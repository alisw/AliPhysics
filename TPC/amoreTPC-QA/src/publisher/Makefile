include ../Makefile.inc

CLEANFILES = *.$(ObjSuf) *.$(DllSuf) *Dict.* *~ core
noinst_HEADERS = TPCPublisherLinkDef.h
pkginclude_HEADERS = PublisherQA.h
BUILT_SOURCES = TPCPublisherDict.h TPCPublisherDict.cxx

bin_PROGRAMS = libAmoreTPCPublisher.$(DllSuf)
libAmoreTPCPublisher_CPPFLAGS = $(CPPFLAGS) -I../common -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/RAW -I$(ALICE_ROOT)/TPC
libAmoreTPCPublisher_CXXFLAGS = $(CXXFLAGS)
libAmoreTPCPublisher_LIBS = $(LIBS) -lAmoreTPCCommon -lProof -lEG -lVMC -lGeom -lMinuit -lRAWDatabase -lSTEERBase -lSTEER -lESD -lRAWDatarec -lCDB -lTPCbase -lTPCrec -lTPCsim
libAmoreTPCPublisher_LDFLAGS = $(LDFLAGS) -L../common -L$(AMORE_SITE)/lib -L$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET) -shared -Wl,-soname,$(bin_PROGRAMS) -o $(bin_PROGRAMS)
libAmoreTPCPublisher_SOURCES = TPCPublisherDict.cxx PublisherQA.cxx
libAmoreTPCPublisher_OBJECTS = TPCPublisherDict.$(ObjSuf) $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(wildcard *.$(SrcSuf)))

%.$(ObjSuf):	%.$(SrcSuf) %.$(HdrSuf)
	@echo $<
	@$(CXX) $(INCLUDES) $(libAmoreTPCPublisher_CPPFLAGS) $(libAmoreTPCPublisher_CXXFLAGS) -c $<

default: all

all: $(bin_PROGRAMS)

$(bin_PROGRAMS): $(libAmoreTPCPublisher_OBJECTS)
	@$(LD) $(libAmoreTPCPublisher_LDFLAGS) $(libAmoreTPCPublisher_LIBS) $^
#$(OutPutOpt)$@

install: all
	@cp $(bin_PROGRAMS) ${AMORE_SITE}/lib;
	
uninstall:
	@rm -rf ${AMORE_SITE}/lib/$(bin_PROGRAMS);

clean:
	@rm -rf $(CLEANFILES)

%Dict.$(SrcSuf) %Dict.$(HdrSuf): $(pkginclude_HEADERS)
	@echo Creating dictionary...
	@rootcint -f $*Dict.$(SrcSuf) -c -p $(INCLUDES) $(libAmoreTPCPublisher_CPPFLAGS) $(libAmoreTPCPublisher_CXXFLAGS) $(pkginclude_HEADERS) $*LinkDef.h
