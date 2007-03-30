PACKAGE = invalid-only-for-proof

include $(ROOTSYS)/test/Makefile.arch
include lib$(PACKAGE).pkg


ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
   ifneq ($(ANALYSIS_NEW_INCLUDE),)
       ALICEINC += -I../$(ANALYSIS_NEW_INCLUDE)
    endif
   ifneq ($(ANALYSIS_INCLUDE),)
       ALICEINC += -I../$(ANALYSIS_INCLUDE)
    endif
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/include
     endif
endif

CXXFLAGS     += $(ALICEINC) -g


SRCS         += G__$(PACKAGE).cxx
OBJS          = $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) *.so G__$(PACKAGE).*

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

