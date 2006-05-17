ObjSuf		= o
SrcSuf		= cxx
HeadSuf		= h

ROOTLIBS	= -L$(ROOTSYS)/lib -lNew -lCint \
		-lGraf -lGraf3d -lHist -lHtml -lMatrix -lMinuit\
		-lPostscript -lTree -lASImage \
		-lCore -lEG -lFoam -lFumili -lGedOld -lGed \
		-lGeom -lGpad -lGraf -lGuiBld -lGui -lMatrix\
		-lRint -lRootAuth  -lThread -lVMC


PEAKFINDERTESTO		= PeakFinderTest.$(ObjSuf)
PEAKFINDERTESTS		= PeakFinderTest.$(SrcSuf)

ALIPHOSPULSEGENERATORO	= AliPHOSPulseGenerator.$(ObjSuf)
ALIPHOSPULSEGENERATORS	= AliPHOSPulseGenerator.$(SrcSuf)
ALIPHOSPULSEGENERATORH	= AliPHOSPulseGenerator.$(HeadSuf)

ALIPHOSFITTERO		= AliPHOSFitter.$(ObjSuf)
ALIPHOSFITTERS		= AliPHOSFitter.$(SrcSuf)
ALIPHOSFITTERH		= AliPHOSFitter.$(HeadSuf)

OBJS	= $(PEAKFINDERTESTO)  $(ALIPHOSPULSEGENERATORO) $(ALIPHOSFITTERO)

pftest:$(OBJS)
	g++ $(ROOTLIBS) $(ROOTGLIBS) $(OBJS) -o pftest

$(PEAKFINDERTESTO):$(PEAKFINDERTESTS)
	g++ -I$(ROOTSYS)/include -c $(PEAKFINDERTESTS)

$(ALIPHOSPULSEGENERATORO):$(ALIPHOSPULSEGENERATORS) $(ALIPHOSPULSEGENERATORH)
	g++ -I$(ROOTSYS)/include -c $(ALIPHOSPULSEGENERATORS)

$(ALIPHOSFITTERO):$(ALIPHOSFITTERS) $(ALIPHOSFITTERH)
	g++ -I$(ROOTSYS)/include -c $(ALIPHOSFITTERS)

clean:
	@rm -f *o *~ *# pftest
