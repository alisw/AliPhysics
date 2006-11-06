
include ../../../build/Makefile.$(ALICE_TARGET)


CLASSES = TestServer.cxx 
DICT = DictTest.cxx

HEADERS := $(patsubst %.cxx,%.h,$(CLASSES))

SRCS = $(CLASSES) $(DICT)

OBJS := $(patsubst %.cxx,.obj/%.o,$(SRCS))
DEPS := $(patsubst .obj/%.o,.dep/%.d,$(OBJS))

INCDIR = -I$(shell root-config --incdir) -I../ -I../../../include 
CXXFLAGS += $(INCDIR)
LIBDIR = $(shell root-config --libdir)


SHARED = -shared

.PHONY: all clean distclean dep

all: $(DEP) libTest.so

$(DICT): $(HEADERS) LinkDef.h
	rootcint -f $@ -c $(INCDIR) $(HEADERS) LinkDef.h 


libTest.so: $(OBJS)
	$(LD) $(SHARED) $(LDFLAGS) $^ -L$(LIBDIR) $(LIBS) -o $@


dep: $(DICT) $(DEPS)

clean:
	rm -rf DictTest.h
	rm -rf DictTest.cxx
	rm -rf .dep
	rm -rf .obj
	rm -rf libTest.so

distclean: clean 

ifeq ($(filter %clean,$(MAKECMDGOALS)),)
DEPS_MAGIC := $(shell mkdir .dep > /dev/null 2>&1 || :)
endif
ifeq ($(strip $(filter %clean,$(MAKECMDGOALS)) $(filter dep,$(MAKECMDGOALS))),)
OBJS_MAGIC := $(shell mkdir .obj > /dev/null 2>&1 || :)
-include $(DEPS)
endif

.dep/%.d: %.cxx
	set -e; $(CXX) -MM $(CXXFLAGS) $< | \
	awk '/\.o/ {print ".dep/$(*F).d .obj/"$$0} ! /\.o/ {print $$0}' > .dep/$(*F).d

.obj/%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o .obj/$(*F).o 

