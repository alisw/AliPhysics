      integer function herwig6_common_block_address(common_block_name)
c-----------------------------------------------------------------------
      include "herwig65.inc"

      character*(*) common_block_name
      external      HWUDAT
c-----------------------------------------------------------------------
      integer       herwig6_addressc, herwig6_addressi
      integer       herwig6_addressd, herwig6_addressf
      integer       herwig6_addressl, herwig6_addressdc
      external      herwig6_addressc, herwig6_addressi
      external      herwig6_addressd, herwig6_addressf
      external      herwig6_addressl, herwig6_addressdc
      integer       common_block_address
c-----------------------------------------------------------------------
      common_block_address = 0
c-----------------------------------------------------------------------
      if     (common_block_name.eq."HEPEVT") then
        common_block_address = herwig6_addressi(NEVHEP)
      elseif (common_block_name.eq."HWBEAM") then
        common_block_address = herwig6_addressi(IPART1)
      elseif (common_block_name.eq."HWBMCH") then
        common_block_address = herwig6_addressc(PART1)
      elseif (common_block_name.eq."HWPROC") then
        common_block_address = herwig6_addressd(EBEAM1)
      elseif (common_block_name.eq."HWPRAM") then
        common_block_address = herwig6_addressd(AFCH)
      elseif (common_block_name.eq."HWPRCH") then
        common_block_address = herwig6_addressc(AUTPDF)
      elseif (common_block_name.eq."HWPART") then
        common_block_address = herwig6_addressi(NEVPAR)
      elseif (common_block_name.eq."HWPARP") then
        common_block_address = herwig6_addressd(DECPAR)
      elseif (common_block_name.eq."HWBOSC") then
        common_block_address = herwig6_addressd(ALPFAC)
      elseif (common_block_name.eq."HWPARC") then
        common_block_address = herwig6_addressi(JCOPAR)
      elseif (common_block_name.eq."HWBRCH") then
        common_block_address = herwig6_addressd(ANOMSC)
      elseif (common_block_name.eq."HWEVNT") then
        common_block_address = herwig6_addressd(AVWGT)
      elseif (common_block_name.eq."HWHARD") then
        common_block_address = herwig6_addressd(ASFIXD)
      elseif (common_block_name.eq."HWPROP") then
        common_block_address = herwig6_addressd(RLTIM)
      elseif (common_block_name.eq."HWUNAM") then
        common_block_address = herwig6_addressc(RNAME)
      elseif (common_block_name.eq."HWUPDT") then
        common_block_address = herwig6_addressd(BRFRAC)
      elseif (common_block_name.eq."HWUWTS") then
        common_block_address = herwig6_addressd(REPWT)
      elseif (common_block_name.eq."HWUCLU") then
        common_block_address = herwig6_addressd(CLDKWT)
      elseif (common_block_name.eq."HWDIST") then
        common_block_address = herwig6_addressd(EXAG)
      elseif (common_block_name.eq."HWQDKS") then
        common_block_address = herwig6_addressd(VTXQDK)
      elseif (common_block_name.eq."HWUSUD") then
        common_block_address = herwig6_addressd(ACCUR)
      elseif (common_block_name.eq."HWSUSY") then
        common_block_address = herwig6_addressd(TANB)
      elseif (common_block_name.eq."HWRPAR") then
        common_block_address = herwig6_addressd(LAMDA1)
      elseif (common_block_name.eq."HWMINB") then
        common_block_address = herwig6_addressd(PMBN1)
      elseif (common_block_name.eq."HWCLUS") then
        common_block_address = herwig6_addressd(PPCL)
      elseif (common_block_name.eq."HWGRAV") then
        common_block_address = herwig6_addressd(GRVLAM)
      elseif (common_block_name.eq."HW6202") then
        common_block_address = herwig6_addressd(VIPWID)
      elseif (common_block_name.eq."HW6203") then
        common_block_address = herwig6_addressd(ABWGT)
      elseif (common_block_name.eq."HW6300") then
        common_block_address = herwig6_addressd(MJJMIN)
      elseif (common_block_name.eq."HWPMRS") then
        common_block_address = herwig6_addressd(FMRS)
      elseif (common_block_name.eq."HWCIRC") then
        common_block_address = herwig6_addressi(CIRCOP)
      elseif (common_block_name.eq."HWDSPB") then
        common_block_address = herwig6_addressd(ABMODE)
      elseif (common_block_name.eq."HWDSP2") then
        common_block_address = herwig6_addressd(A2MODE)
      elseif (common_block_name.eq."HWDSP3") then
        common_block_address = herwig6_addressd(A3MODE)
      elseif (common_block_name.eq."HWDSP4") then
        common_block_address = herwig6_addressd(A4MODE)
      elseif (common_block_name.eq."HWDSPN") then
        common_block_address = herwig6_addressi(NDECSY)
      elseif (common_block_name.eq."HWSPIN") then
        common_block_address = herwig6_addressdc(MESPN)
      elseif (common_block_name.eq."HWSTAU") then
        common_block_address = herwig6_addressi(JAK1)
      elseif (common_block_name.eq."HWGUPR") then
        common_block_address = herwig6_addressd(LHWGT)
      elseif (common_block_name.eq."HW6500") then
        common_block_address = herwig6_addressl(PRESPL)
      elseif (common_block_name.eq."HW6504") then
        common_block_address = herwig6_addressi(ITOPRD)
      elseif (common_block_name.eq."HW6506") then
        common_block_address = herwig6_addressd(PDFX0)
      endif
c-----------------------------------------------------------------------
      herwig6_common_block_address = common_block_address
      end
