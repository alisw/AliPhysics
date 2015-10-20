# $Id$

PACKAGE = TestShuttle

include $(ROOTSYS)/etc/Makefile.arch
include $(PACKAGE).pkg

ALICEINC      = -I.
ifneq ($(ALICE_ROOT),)
  ALICEINC += -I$(ALICE_ROOT)/include
endif

CXXFLAGS     += $(ALICEINC) -g
LDFLAGS      += `root-config --ldflags --libs` -L$(ALICE_ROOT)/lib -lACORDEbase -lADbase -lSTEER -lCDB -lSTEERBase -lESD
SRCS         += dict.cxx
OBJS          = $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) -o $@
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) core *~ *.so dict.* $(PARFILE)

dict.cxx dict.h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	@[ -f dict ] || rm -f $@
	rootcint -f $@ -c $(ALICEINC) $^
