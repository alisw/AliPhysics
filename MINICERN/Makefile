############################### MINICERN Makefile #############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = minicern

INC	= -Ipacklib/zebra \
          -Ipacklib/kernlib/kerngen \
 	  -Imathlib/gen \
	  -Ipacklib/kernlib/kernnum

 
# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(INC) $(CLIBDEFS)
 
# C compilation flags

CFLAGS      = $(COPT) $(CLIBCOPT) $(INC) $(CLIBDEFS)
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(INC) $(CLIBDEFS)
 
FSRCS	= dummies.F hbcdf.F fzdummies.F \
	  mathlib/gen/b/atg.F \
	  mathlib/gen/c/gamma.F \
	  mathlib/gen/c/gamma128.F \
	  mathlib/gen/c/gamma64.F \
	  mathlib/gen/d/gauss.F \
	  mathlib/gen/d/gauss128.F \
	  mathlib/gen/d/gauss64.F \
	  mathlib/gen/g/prob.F \
	  mathlib/gen/n/mtlprt.F \
	  mathlib/gen/n/mtlset.F \
	  mathlib/gen/v/ranlux.F \
	  mathlib/gen/v/rnorml.F \
	  mathlib/gen/v/rnormx.F \
	  mathlib/gen/v/rnpssn.F \
	  packlib/kernlib/kerngen/ccgen/getenvf.F \
	  packlib/kernlib/kerngen/ccgencf/cfopen.F \
	  packlib/kernlib/kerngen/ccgencf/cfstat.F \
	  packlib/kernlib/kerngen/tcgen/cfill.F \
	  packlib/kernlib/kerngen/tcgen/cleft.F \
	  packlib/kernlib/kerngen/tcgen/cltou.F \
	  packlib/kernlib/kerngen/tcgen/cross.F \
	  packlib/kernlib/kerngen/tcgen/cutol.F \
	  packlib/kernlib/kerngen/tcgen/flpsor.F \
	  packlib/kernlib/kerngen/tcgen/gshowf.F \
	  packlib/kernlib/kerngen/tcgen/icluns.F \
	  packlib/kernlib/kerngen/tcgen/iumode.F \
	  packlib/kernlib/kerngen/tcgen/kerngt.F \
	  packlib/kernlib/kerngen/tcgen/lenocc.F \
	  packlib/kernlib/kerngen/tcgen/lorenf.F \
	  packlib/kernlib/kerngen/tcgen/ranmar.F \
	  packlib/kernlib/kerngen/tcgen/traceq.F \
	  packlib/kernlib/kerngen/tcgen/ubits.F \
	  packlib/kernlib/kerngen/tcgen/uoptc.F \
	  packlib/kernlib/kerngen/tcgen/utrans.F \
	  packlib/kernlib/kerngen/tcgenc/datimh.F \
	  packlib/kernlib/kerngen/tcgenc/intrac.F \
	  packlib/kernlib/kerngen/tcgenc/timed.F \
	  packlib/kernlib/kerngen/tcgenc/timel.F \
	  packlib/kernlib/kerngen/tcgens/iucomp.F \
	  packlib/kernlib/kerngen/tcgens/iufind.F \
	  packlib/kernlib/kerngen/tcgens/iusame.F \
	  packlib/kernlib/kerngen/tcgens/jbit.F \
	  packlib/kernlib/kerngen/tcgens/jbyt.F \
	  packlib/kernlib/kerngen/tcgens/jrsbyt.F \
	  packlib/kernlib/kerngen/tcgens/sbit0.F \
	  packlib/kernlib/kerngen/tcgens/sbit1.F \
	  packlib/kernlib/kerngen/tcgens/sbyt.F \
	  packlib/kernlib/kerngen/tcgens/ublow.F \
	  packlib/kernlib/kerngen/tcgens/ubunch.F \
	  packlib/kernlib/kerngen/tcgens/ucopiv.F \
	  packlib/kernlib/kerngen/tcgens/ucopy.F \
	  packlib/kernlib/kerngen/tcgens/ucopy2.F \
	  packlib/kernlib/kerngen/tcgens/uctoh.F \
	  packlib/kernlib/kerngen/tcgens/uctoh1.F \
	  packlib/kernlib/kerngen/tcgens/ufill.F \
	  packlib/kernlib/kerngen/tcgens/uhtoc.F \
	  packlib/kernlib/kerngen/tcgens/upkbyt.F \
	  packlib/kernlib/kerngen/tcgens/vblank.F \
	  packlib/kernlib/kerngen/tcgens/vfill.F \
	  packlib/kernlib/kerngen/tcgens/vmax.F \
	  packlib/kernlib/kerngen/tcgens/vmod.F \
	  packlib/kernlib/kerngen/tcgens/vscale.F \
	  packlib/kernlib/kerngen/tcgens/vunit.F \
	  packlib/kernlib/kerngen/tcgens/vzero.F \
	  packlib/kernlib/kernnum/e104fort/fint.F \
	  packlib/kernlib/kernnum/e106fort/locatf.F \
	  packlib/kernlib/kernnum/f010fort/f010pr.F \
	  packlib/kernlib/kernnum/f010fort/rinv.F \
	  packlib/kernlib/kernnum/f011fort/rfact.F \
	  packlib/kernlib/kernnum/f011fort/rfinv.F \
	  packlib/kernlib/kernnum/f011fort/tmprnt.F \
	  packlib/kernlib/kernnum/n001fort/kerset.F \
	  packlib/zebra/dzebra/dzare1.F \
	  packlib/zebra/dzebra/dzbkdv.F \
	  packlib/zebra/dzebra/dzbkhd.F \
	  packlib/zebra/dzebra/dzbkup.F \
	  packlib/zebra/dzebra/dzbkxr.F \
	  packlib/zebra/dzebra/dzdata.F \
	  packlib/zebra/dzebra/dziopd.F \
	  packlib/zebra/dzebra/dzioph.F \
	  packlib/zebra/dzebra/dziopr.F \
	  packlib/zebra/dzebra/dzmap.F \
	  packlib/zebra/dzebra/dzopt.F \
	  packlib/zebra/dzebra/dzshpr.F \
	  packlib/zebra/dzebra/dzsnap.F \
	  packlib/zebra/dzebra/dzswap.F \
	  packlib/zebra/dzebra/dztext.F \
	  packlib/zebra/dzebra/dztyp.F \
	  packlib/zebra/dzebra/nzleft.F \
	  packlib/zebra/dzebra/zpage.F \
	  packlib/zebra/mq/mzbook.F \
	  packlib/zebra/mq/mzdiv.F \
	  packlib/zebra/mq/mzdrop.F \
	  packlib/zebra/mq/mzebra.F \
	  packlib/zebra/mq/mzend.F \
	  packlib/zebra/mq/mzflag.F \
	  packlib/zebra/mq/mzinco.F \
	  packlib/zebra/mq/mzinqd.F \
	  packlib/zebra/mq/mzinqs.F \
	  packlib/zebra/mq/mzlift.F \
	  packlib/zebra/mq/mzlink.F \
	  packlib/zebra/mq/mzlint.F \
	  packlib/zebra/mq/mzlogl.F \
	  packlib/zebra/mq/mzmaps.F \
	  packlib/zebra/mq/mzneed.F \
	  packlib/zebra/mq/mzpudx.F \
	  packlib/zebra/mq/mzpush.F \
	  packlib/zebra/mq/mzstor.F \
	  packlib/zebra/mq/mzvers.F \
	  packlib/zebra/mq/mzwork.F \
	  packlib/zebra/mqg/mzfgap.F \
	  packlib/zebra/mqg/mzgar1.F \
	  packlib/zebra/mqg/mzgarb.F \
	  packlib/zebra/mqg/mzgsta.F \
	  packlib/zebra/mqg/mzmove.F \
	  packlib/zebra/mqg/mzrelb.F \
	  packlib/zebra/mqg/mzrell.F \
	  packlib/zebra/mqg/mzrelx.F \
	  packlib/zebra/mqg/mztabc.F \
	  packlib/zebra/mqg/mztabf.F \
	  packlib/zebra/mqg/mztabh.F \
	  packlib/zebra/mqg/mztabm.F \
	  packlib/zebra/mqg/mztabr.F \
	  packlib/zebra/mqg/mztabs.F \
	  packlib/zebra/mqg/mztabx.F \
	  packlib/zebra/mqg/mzwipe.F \
	  packlib/zebra/mqs/mzchln.F \
	  packlib/zebra/mqs/mzchls.F \
	  packlib/zebra/mqs/mzchnb.F \
	  packlib/zebra/mqs/mzdvac.F \
	  packlib/zebra/mqs/mzfdiv.F \
	  packlib/zebra/mqs/mzform.F \
	  packlib/zebra/mqs/mziocf.F \
	  packlib/zebra/mqs/mzioch.F \
	  packlib/zebra/mqs/mziocr.F \
	  packlib/zebra/mqs/mzixco.F \
	  packlib/zebra/mqs/mzresv.F \
	  packlib/zebra/mqs/mzsdiv.F \
	  packlib/zebra/mqs/mzxref.F \
	  packlib/zebra/mqs/mzxrup.F \
	  packlib/zebra/mqv/zvautx.F \
	  packlib/zebra/qend/zabend.F \
	  packlib/zebra/qend/zend.F \
	  packlib/zebra/qend/zfatal.F \
	  packlib/zebra/qend/zfatam.F \
	  packlib/zebra/qend/zphase.F \
	  packlib/zebra/qend/zpostm.F \
	  packlib/zebra/qend/ztell.F \
	  packlib/zebra/qend/ztelus.F \
	  packlib/zebra/qutil/izbcd.F \
	  packlib/zebra/qutil/izbcdt.F \
	  packlib/zebra/qutil/izhnum.F \
	  packlib/zebra/qutil/qnext.F \
	  packlib/zebra/qutil/zhtoi.F  \
	  packlib/zebra/qutil/zitoh.F \
	  packlib/zebra/qutil/zshunt.F \
	  packlib/zebra/rz/rzallo.F \
	  packlib/zebra/rz/rzcdir.F \
	  packlib/zebra/rz/rzdate.F \
	  packlib/zebra/rz/rzdelk.F \
	  packlib/zebra/rz/rzdlok.F \
	  packlib/zebra/rz/rzend.F \
	  packlib/zebra/rz/rzexpd.F \
	  packlib/zebra/rz/rzfdir.F \
	  packlib/zebra/rz/rzfile.F \
	  packlib/zebra/rz/rzfree.F \
	  packlib/zebra/rz/rzin.F \
	  packlib/zebra/rz/rzink.F \
	  packlib/zebra/rz/rzins.F \
	  packlib/zebra/rz/rziodo.F \
	  packlib/zebra/rz/rzllok.F \
	  packlib/zebra/rz/rzlock.F \
	  packlib/zebra/rz/rzmake.F \
	  packlib/zebra/rz/rzmdir.F \
	  packlib/zebra/rz/rzmods.F \
	  packlib/zebra/rz/rzobkn.F \
	  packlib/zebra/rz/rzopen.F \
	  packlib/zebra/rz/rzout.F \
	  packlib/zebra/rz/rzpaff.F \
	  packlib/zebra/rz/rzpath.F \
	  packlib/zebra/rz/rzprnk.F \
	  packlib/zebra/rz/rzpurf.F \
	  packlib/zebra/rz/rzread.F \
	  packlib/zebra/rz/rzrtop.F \
	  packlib/zebra/rz/rzsame.F \
	  packlib/zebra/rz/rzsave.F \
	  packlib/zebra/rz/rzstrip.F \
	  packlib/zebra/rz/rzstrir.F \
	  packlib/zebra/rz/rzused.F \
	  packlib/zebra/rz/rzvcyc.F \
	  packlib/zebra/rz/rzver1.F \
	  packlib/zebra/rz/rzver2.F \
	  packlib/zebra/rz/rzveri.F \
	  packlib/zebra/rz/rzwrt.F \
 
 
ifeq ($(PLATFORM),linux)
	  FSRC = $(FSRCS) \
	  packlib/kernlib/kerngen/unix/lnxgs/rdmin.F
endif
ifeq ($(PLATFORM),hpux)
	  FSRC = $(FSRCS)
endif
ifeq ($(PLATFORM),solaris)
	  FSRC = $(FSRCS) \
          packlib/kernlib/kerngen/tcgen/uset.F
endif
ifeq ($(PLATFORM),alpha)
	  FSRC = $(FSRCS) 
endif

CSRCS	= \
	  packlib/kernlib/kerngen/ccgen/abend.c \
	  packlib/kernlib/kerngen/ccgen/geteni.c \
	  packlib/kernlib/kerngen/ccgen/jumptn.c \
	  packlib/kernlib/kerngen/ccgen/jumpxn.c \
	  packlib/kernlib/kerngen/ccgen/locf.c \
	  packlib/kernlib/kerngen/ccgen/qnexte.c \
	  packlib/kernlib/kerngen/ccgen/vxinvb.c \
	  packlib/kernlib/kerngen/ccgencf/cfclos.c \
	  packlib/kernlib/kerngen/ccgencf/cfget.c \
	  packlib/kernlib/kerngen/ccgencf/cfopei.c \
	  packlib/kernlib/kerngen/ccgencf/cfput.c \
	  packlib/kernlib/kerngen/ccgencf/cfseek.c \
	  packlib/kernlib/kerngen/ccgencf/cfstati.c \
	  packlib/kernlib/kerngen/ccgenu/fchput.c \
	  packlib/kernlib/kerngen/ccgenu/fchtak.c \
 	  packlib/kernlib/kerngen/ccgen/datime.c \
 
ifeq ($(PLATFORM),linux)
	  CSRC = $(CSRCS) \
	  packlib/kernlib/kerngen/unix/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),hpux)
	  CSRC = $(CSRCS) \
	  packlib/kernlib/kerngen/unix/hpxgs/traceqc.c \
          packlib/kernlib/kerngen/ccgen/lnblnk.c
endif
ifeq ($(PLATFORM),solaris)
	  CSRC = $(CSRCS) 
endif
ifeq ($(PLATFORM),alpha)
	  CSRC = $(CSRCS) 
endif
	

CXXSRC	=

FOBJ	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRC))
COBJ	= $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRC))
CXXOBJ	= $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(CXXSRC))

SRCS	= $(FSRC) $(CSRC) $(CXXSRC)
OBJS	= $(FOBJ) $(COBJ) $(CXXOBJ)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libminicern.$(SL)
ALIBRARY	= $(LIBDIR)/libminicern.a

ifeq ($(PLATFORM),alpha)
        default:	$(ALIBRARY) $(SLIBRARY)
else
        default:	$(SLIBRARY)
endif

$(LIBDIR)/libminicern.$(SL):  $(OBJS)
$(LIBDIR)/libminicern.a:  $(OBJS)

depend:			$(SRCS)

TOCLEAN		= $(OBJS)

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
