include Makefile.arch

ROOTLIBS         = `root-config --libs`  -lVMC -lGeom -lEG
ALIROOTLIBS      = -L$(ALICE_ROOT)/lib/tgt_linux  -lESD -lRAWDatarec -lRAWDatabase -lSTEERBase


#LIBS		 = $(LIBS) $(ROOTLIBS) $(ALIROOTLIBS)
LIBS		 = $(ROOTLIBS) $(ALIROOTLIBS)

ROOTINCLUDE      = -I$(ROOTSYS)/include
ALIROOTINCLUDE   = -I$(ALICE_ROOT)/include  -I$(ALICE_ROOT)/RAW
INCLUDES	 =  $(ROOTINCLUDE) $(ALIROOTINCLUDE)

#MAINO       = testDecoder.$(ObjSuf) AliHLTAltroData.$(ObjSuf)  AliHLTAltroBunch.$(ObjSuf) AliHLTDDLDecoder.$(ObjSuf)
#MAINS       = testDecoder.$(SrcSuf) AliHLTAltroData.$(SrcSuf)  AliHLTAltroBunch.$(SrcSuf) AliHLTDDLDecoder.$(SrcSuf)

MAINO       = testDecoder.$(ObjSuf)
MAINS       = testDecoder.$(SrcSuf) 
MAIN        = testDecoder$(ExeSuf)

OBJS          = $(MAINO) 
PROGRAMS      = $(MAIN) 


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)


$(MAIN):	$(MAINO)
		$(LD) $(LDFLAGS) $(INCLUDES) $^ $(LIBS) $(GLIBS) -lTreePlayer -lMathCore -lThread -lpthread $(OutPutOpt)$(MAIN)
		@echo "$@ done"


clean:
		@rm -f $(OBJS) core

distclean:      clean
		-@mv -f linearIO.root linearIO.roott
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations
		@rm -rf cxx_repository
		-@mv -f linearIO.roott linearIO.root
		-@cd RootShower && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)


#AliHLTAltroData.$(ObjSuf): AliHLTAltroData.h
#AliHLTDDLDecoder.$(ObjSuf): AliHLTDDLDecoder.h
#testDecoder.$(ObjSuf): AliHLTAltroData.h AliHLTDDLDecoder.h
testDecoder.$(ObjSuf): 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<
