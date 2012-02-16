# this make file is used to quickly compile locally the code (to check
# for errors) and to produce par files + tags

OBJECTS =  AliAnalysisTaskPerformanceStrange.o #AliAnalysisLambdaK0HistoManager.o AliAnalysisTaskLambdaK0AOD.o AliAnalysisLambdaK0Cuts.o
HEADERS = $(OBJECTS:.o=.h) 
MACROS = 

 #DrawRatiosVsPt_C.so
	  #Get_mass_histos_C.so LoadHistoAndRecorrect_C.so FitPt2Par_C.so
 #Get_pt_histos_debug_C.so \
MOREHEADERS = 
BASENAME = BG

#############################################################
TARGET   = lib$(BASENAME).so
DICT     = lib$(BASENAME)Dict

CXX = g++
CXXFLAGS = -O -Wall -fPIC -D_REENTRANT -I$(ROOTSYS)/include \
	-I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/obj/include -I$(ALICE_ROOT)/PWG0 -I$(ALICE_ROOT)/ANALYSIS	-I$(ALICE_ROOT)/PWG1/background \
	-c -g
CINTFLAGS = -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/obj/include -I$(ALICE_ROOT)/PWG0 -I$(ALICE_ROOT)/ANALYSIS 

all:  $(MACROS) $(TARGET)
#############################################################
%.o: %.cxx %.h $(HEADERS)
	$(CXX) $(CXXFLAGS) $<

%_C.so: %.C 
	echo "{" > tmp.C
	echo "gROOT->LoadMacro(\"compile.C\");" >> tmp.C
	echo "compile(\"$<\");" >> tmp.C
	echo "}" >> tmp.C
	echo ".x tmp.C" | root -b -q

%_cxx.so: %.cxx
	echo "{" > tmp.C
	echo "gROOT->LoadMacro(\"compile.C\");" >> tmp.C
	echo "compile(\"$<\");" >> tmp.C
	echo "}" >> tmp.C
	echo ".x tmp.C" | root -b
#	echo "gROOT->LoadMacro(\"IncludeMichele.C\");" >> tmp.C 
#	echo "IncludeMichele();" >> tmp.C 	
#	echo "gROOT->LoadMacro(\"$<+g\");" >> tmp.C
#	echo "}" >> tmp.C
#	echo ".x tmp.C" | aliroot -b


$(DICT).cxx: $(HEADERS) LinkDef.h
	rootcint  -f $@ -c $(CINTFLAGS) -p $(MOREHEADERS) $^

$(DICT).o: $(DICT).cxx
	$(CXX) $(CXXFLAGS) $<


$(TARGET): $(OBJECTS) $(DICT).o 
	g++ -shared -O $^ -o $@ 
#	cd .. ; make

tags:
	cd $(ALICE_ROOT); \
	find . -iname "*.cxx  " | grep -v target | grep -v tgt | etags - ;\
	find . -iname "*.h"   | grep -v target | grep -v tgt | grep -v include | etags --append - 	
#	find . -iname "*.cxx  " | grep -v target | grep -v tgt | etags \
#	--include=$(ALICE_ROOT)/TAGS --include=$(ROOTSYS)/TAGS  -
#	find . -iname "*.h"   | grep -v target | grep -v tgt | grep -v include | etags --append -

pars:
	cd $(ALICE_ROOT)/obj ; \
	make STEERBase.par; \
	make ESD.par ;\
	make AOD.par ;\
	make ANALYSIS.par ;\
	make ANALYSISalice.par ;\
	make CORRFW.par	;\
	make PWG0base.par

.PHONY: clean

clean:
	@rm *.o *.so *Dict.cxx

