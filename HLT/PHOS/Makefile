include Makefile.arch


MAINO       = testDecoder.$(ObjSuf) AliHLTAltroData.$(ObjSuf)  AliHLTAltroBunch.$(ObjSuf) AliHLTDDLDecoder.$(ObjSuf)

MAINS       = testDecoder.$(SrcSuf) AliHLTAltroData.$(SrcSuf)  AliHLTAltroBunch.$(SrcSuf) AliHLTDDLDecoder.$(SrcSuf)

MAIN        = testDecoder$(ExeSuf)


OBJS          = $(MAINO) 
PROGRAMS      = $(MAIN) 


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)


$(MAIN):	$(MAINO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(GLIBS) -lTreePlayer -lMathCore -lThread -lpthread $(OutPutOpt)$(MAIN)
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
	
AliHLTAltroData.$(ObjSuf): AliHLTAltroData.h
AliHLTDDLDecoder.$(ObjSuf): AliHLTDDLDecoder.h
testDecoder.$(ObjSuf): AliHLTAltroData.h AliHLTDDLDecoder.h

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
