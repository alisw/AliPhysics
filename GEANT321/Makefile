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
GFOBJ	= $(GFSRC:.F=.o)
GCOBJ	= $(GCSRC:.c=.o)
GCXXOBJ	= $(GCXXSRC:.cxx=.o)
GOBJ	= $(GFOBJ) $(GCOBJ) $(GCXXOBJ)

XFSRC	= $(wildcard $(patsubst %,%/*.F,$(XDIRS))) 
XCSRC	= $(wildcard $(patsubst %,%/*.c,$(XDIRS))) 

XSRC	= $(XFSRC) $(XCSRC) $(XCXXSRC)
XFOBJ	= $(XFSRC:.F=.o)
XCOBJ	= $(XCSRC:.c=.o)
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

SLIBRARY	= $(ALICE_ROOT)/lib/libgeant321.$(SL) $(ALICE_ROOT)/lib/libgxint321.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libgeant321.a $(ALICE_ROOT)/lib/libgxint321.a

ifeq ($(PLATFORM),alpha)
	default:	$(ALIBRARY)
else
	default:	$(SLIBRARY)
endif

$(ALICE_ROOT)/lib/libgeant321.$(SL) $(ALICE_ROOT)/lib/libgeant321.a:  $(GOBJ)
$(ALICE_ROOT)/lib/libgxint321.$(SL) $(ALICE_ROOT)/lib/libgxint321.$(SL):  $(XOBJ)

depend:			$(SRCS)

TOCLEAN		= $(GOBJ) $(XOBJ)

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
