include $(ROOTSYS)/test/Makefile.arch
include libPWG0base.pkg

ALICEINC = -I$(ALICE_ROOT)/include
CXXFLAGS += $(ALICEINC)
SRCS += dict.cxx
OBJS          = $(SRCS:.cxx=.o)
PROGRAM       = libPWG0.so

all:          $(PROGRAM)

$(PROGRAM):     $(OBJS)
	@echo "Linking" $(PROGRAM) ...
	@/bin/rm -f $(PROGRAM)
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $(PROGRAM)
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) -o $(PROGRAM)
endif
	@chmod a+x $(PROGRAM)
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) core *~ *.so dict.*

dict.cxx dict.h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	@[ -f dict ] || rm -f $@
	rootcint -f $@ -c $(ALICEINC) $^
