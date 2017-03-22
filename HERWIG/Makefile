# Makefile to compile Herwig6 as a shared library for
# ROOT. Needs the Herwig6 routines and the F77 to C++
# interface routines

# General Definitions for linux compilers

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = HERWIG

# C++ sources

SRCS          =

# C sources

CSRCS       = main.c  herwig6_address.c 

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

##### MACROS #####

FSRCS	= HERWIG6100.f herwig6_common_block_address.f herwig6_called_from_cc.f 

FOBJS	= $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

SRCS	= $(FSRCS) $(CSRCS)
OBJS	= $(FOBJS) $(COBJS)

DSRCS	= 

DOBJS	= $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(DSRCS))

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS)

# C compilation flags

CFLAGS      = $(COPT) $(CLIBCOPT) $(CLIBDEFS)

# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS)



##### TARGETS #####

# Target

SLIBRARY	= $(LIBDIR)/libherwig6.$(SL)
ALIBRARY	= $(LIBDIR)/libherwig6.a

default:	$(SLIBRARY)

$(LIBDIR)/libherwig6.$(SL):  $(OBJS)

depend:	          $(CSRCS) $(DSRCS)

TOCLEAN		= $(OBJS) $(DOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 





























ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
