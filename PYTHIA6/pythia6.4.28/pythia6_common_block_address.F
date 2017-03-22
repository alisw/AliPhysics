#if defined (__ia64) || defined (__x86_64)
      integer*8 
#else
      integer
#endif
     + function pythia6_common_block_address(common_block_name)
c-----------------------------------------------------------------------
c  revision history:
c  -----------------
c *0001 Jun 29 1998 P.Murat: created
c *0002 Apr 30 1999 P.Murat: force loading in the BLOCK DATA
c-----------------------------------------------------------------------
#include "pythia6.inc"
c#include "stdhep.inc"
c
      character*(*) common_block_name
      external      pydata
c
c-----------------------------------------------------------------------
#if defined (__ia64) || defined (__x86_64)
      integer*8
#else
      integer
#endif
     +     pythia6_addressc, pythia6_addressi,
     +     pythia6_addressd, pythia6_addressf,
     +     common_block_address
      external      pythia6_addressc, pythia6_addressi
      external      pythia6_addressd, pythia6_addressf
c-----------------------------------------------------------------------
      common_block_address = 0
c-----------------------------------------------------------------------
c  GEANT common-blocks
c-----------------------------------------------------------------------
      if     (common_block_name(1:6).eq."HEPEVT") then
        common_block_address = pythia6_addressi(NEVHEP)
      elseif (common_block_name(1:6).eq."PYJETS") then
        common_block_address = pythia6_addressi(N)
      elseif (common_block_name(1:6).eq."PYDAT1") then
        common_block_address = pythia6_addressi(MSTU)
      elseif (common_block_name(1:6).eq."PYDAT2") then
        common_block_address = pythia6_addressi(KCHG)
      elseif (common_block_name(1:6).eq."PYDAT3") then
        common_block_address = pythia6_addressi(MDCY)
      elseif (common_block_name(1:6).eq."PYDAT4") then
        common_block_address = pythia6_addressc(CHAF)
      elseif (common_block_name(1:6).eq."PYDATR") then
        common_block_address = pythia6_addressi(MRPY)
      elseif (common_block_name(1:6).eq."PYSUBS") then
        common_block_address = pythia6_addressi(MSEL)
      elseif (common_block_name(1:6).eq."PYPARS") then
        common_block_address = pythia6_addressi(MSTP)
      elseif (common_block_name(1:6).eq."PYINT1") then
        common_block_address = pythia6_addressi(MINT)
      elseif (common_block_name(1:6).eq."PYINT2") then
        common_block_address = pythia6_addressi(ISET)
      elseif (common_block_name(1:6).eq."PYINT3") then
        common_block_address = pythia6_addressd(XSFX)
      elseif (common_block_name(1:6).eq."PYINT4") then
        common_block_address = pythia6_addressi(MWID)
      elseif (common_block_name(1:6).eq."PYINT5" ) then
        common_block_address = pythia6_addressi(NGENPD)
      elseif (common_block_name(1:6).eq."PYINT6" ) then
        common_block_address = pythia6_addressc(PROC)
      elseif (common_block_name(1:6).eq."PYINT7") then
        common_block_address = pythia6_addressd(SIGT)
      elseif (common_block_name(1:6).eq."PYINT8") then
        common_block_address = pythia6_addressd(XPVMD)
      elseif (common_block_name(1:6).eq."PYINT9") then
        common_block_address = pythia6_addressd(VXPVMD)
      elseif (common_block_name(1:6).eq."PYUPPR") then
        common_block_address = pythia6_addressi(NUP)
      elseif (common_block_name(1:6).eq."PYMSSM") then
        common_block_address = pythia6_addressi(IMSS)
      elseif (common_block_name(1:6).eq."PYSSMT") then
        common_block_address = pythia6_addressd(ZMIX)
      elseif (common_block_name(1:6).eq."PYINTS") then
        common_block_address = pythia6_addressd(XXM(1))
      elseif (common_block_name(1:6).eq."PYBINS") then
        common_block_address = pythia6_addressi(IHIST(1))
      endif
c
      pythia6_common_block_address = common_block_address
      end
