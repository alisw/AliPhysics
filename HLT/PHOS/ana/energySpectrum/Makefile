include Makefile.arch


MAINO       = energySpectrum.$(ObjSuf)

MAINS       = energySpectrum.$(SrcSuf)

MAIN        = energySpectrum$(ExeSuf)


OBJS          = $(MAINO) 
PROGRAMS      = $(MAIN) 


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)


$(MAIN):	$(MAINO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(GLIBS) $(ALILIBS) -lTreePlayer -lMathCore -lThread -lpthread $(HLTLIBS) -lHLTbase -lAliHLTUtil -lAliHLTPHOS $(OutPutOpt)$(MAIN)
		@echo "$@ done"
#		$(LD) $(LDFLAGS) $^ $(LIBS) $(GLIBS) -lTreePlayer -lMathCore -lThread -lpthread $(OutPutOpt)$(MAIN)



clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
