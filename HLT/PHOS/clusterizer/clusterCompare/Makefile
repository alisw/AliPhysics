include Makefile.arch


MAINO       = clusterCompare.$(ObjSuf)

MAINS       = clusterCompare.$(SrcSuf)

MAIN        = clusterCompare$(ExeSuf)


OBJS          = $(MAINO) 
PROGRAMS      = $(MAIN) 


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)


$(MAIN):	$(MAINO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(GLIBS) $(ALILIBS) -lTreePlayer -lMathCore -lThread -lpthread -L/home/odjuvsla/Workspace/alice/AliRoot/lib/tgt_linux -lRAWDatabase -lRAWDatarec -lRAWDatasim -lPHOSbase -lPHOSrec -lPHOSsim -lSTEER -lSTEERBase -lTPCbase -lAliHLTUtil -lAliHLTPHOS -lHLTbase -lHLTinterface $(OutPutOpt)$(MAIN)
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
