############################### GEANT321 Makefile #############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = GEANT321

GDIRS = gbase gcons gdraw geocad ggeom gheisha ghits ghrout ghutils giface \
        giopa gkine gparal gphys gscan gstrag gtrak matx55 miface miguti \
        neutron peanut fiface cgpack fluka block comad

XDIRS = gxint cdf

# C++ sources

SRCS          =
##### MACROS #####

GFSRC	= $(wildcard $(patsubst %,%/*.F,$(GDIRS))) gcinit.F
GCSRC	= $(wildcard $(patsubst %,%/*.c,$(GDIRS))) 
GCXXSRC	= $(wildcard $(patsubst %,%/*.cxx,$(GDIRS))) 

GSRC	= $(GFSRC) $(GCSRC)
GFOBJ	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(GFSRC))
GCOBJ	= $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(GCSRC))
GCXXOBJ	= $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(GCXXSRC))
GOBJ	= $(GFOBJ) $(GCOBJ) $(GCXXOBJ)

XFSRC	= $(wildcard $(patsubst %,%/*.F,$(XDIRS))) 
XCSRC	= $(wildcard $(patsubst %,%/*.c,$(XDIRS))) 

XSRC	= $(XFSRC) $(XCSRC) $(XCXXSRC)
XFOBJ	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(XFSRC))
XCOBJ	= $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(XCSRC))
XOBJ	= $(XFOBJ) $(XCOBJ) $(XCXXOBJ)

SRCS	= $(GSRC) $(XSRC)
OBJS	= $(GOBJ) $(XOBJ)

# Global Macros

ALLSRCS       = $(SRCS) $(HDRS) $(FSRCS) $(CSRCS) $(DICT)
ALLOBJS       = $(OBJS) $(FOBJS) $(COBJS) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS) -I.
 
# C compilation flags

CFLAGS      = $(COPT) $(CLIBCOPT) $(CLIBDEFS) -I.
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS) -I.
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libgeant321.$(SL) $(LIBDIR)/libgxint321.$(SL)
ALIBRARY	= $(LIBDIR)/libgeant321.a $(LIBDIR)/libgxint321.a

ifeq ($(PLATFORM),alpha)
        default:	$(ALIBRARY)
else
        default:	$(SLIBRARY)
endif

$(LIBDIR)/libgeant321.$(SL) $(LIBDIR)/libgeant321.a:  $(GOBJ)
$(LIBDIR)/libgxint321.$(SL) $(LIBDIR)/libgxint321.$(SL):  $(XOBJ)

depend:			$(SRCS)

TOCLEAN		= $(GOBJ) $(XOBJ)

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
