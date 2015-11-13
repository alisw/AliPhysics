
include Makefile.arch


CLASSES = AliHLTCreateGRP.cxx 
	
DICT = AliHLTCreateGRPDict.cxx

HEADERS := $(patsubst %.$(SrcSuf),%.h,$(CLASSES))

SRCS = $(CLASSES) $(DICT)

INCLUDES := $(ALICE_ROOT)/include -I$(DIMDIR)/dim
LIBDIR := $(shell root-config --libdir) -L$(DIMDIR)/linux
LIBS += -ldim
SHARED = -shared

CXXFLAGS += -I$(INCLUDES)

OBJS := $(patsubst %.$(SrcSuf),.obj/%.$(ObjSuf),$(SRCS))
DEPS := $(patsubst .obj/%.$(ObjSuf),.dep/%.d,$(OBJS))


.PHONY: all clean distclean dep

all: $(DEP) libAliHLTCreateGRP.so 

$(DICT): $(HEADERS) AliHLTCreateGRP.h
	rootcint -f $@ -c -I$(INCLUDES) $^



libAliHLTCreateGRP.so: $(OBJS)
	$(LD) $(SHARED) $(LDFLAGS) $^ -L$(LIBDIR) $(LIBS) -o $@

dep: $(DICT) $(DEPS)

clean:
	rm -rf AliHLTCreateGRPDict.h
	rm -rf AliHLTCreateGRPDict.cxx
	rm -rf libAliHLTCreateGRP.so
	rm -rf .dep
	rm -rf .obj


distclean: clean ;

ifeq ($(filter %clean,$(MAKECMDGOALS)),)
DEPS_MAGIC := $(shell mkdir .dep > /dev/null 2>&1 || :)
endif
ifeq ($(strip $(filter %clean,$(MAKECMDGOALS)) $(filter dep,$(MAKECMDGOALS))),)
OBJS_MAGIC := $(shell mkdir .obj > /dev/null 2>&1 || :)
-include $(DEPS)
endif

.dep/%.d: %.$(SrcSuf)
	set -e; $(CXX) -MM $(CXXFLAGS) $< | \
	awk '/\.o/ {print ".dep/$(*F).d .obj/"$$0} ! /\.o/ {print $$0}' > .dep/$(*F).d

.obj/%.$(ObjSuf): %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $< -o .obj/$(*F).o 

