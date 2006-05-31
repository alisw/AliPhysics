PACKAGE = invalid-only-for-proof

include $(ROOTSYS)/test/Makefile.arch
include lib$(PACKAGE).pkg

SRCS         += dict.cxx
OBJS          = $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) -o $@
endif
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) core *~ *.so dict.*

dict.cxx dict.h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	@[ -f dict ] || rm -f $@
	rootcint -f $@ -c $(ALICEINC) $^

