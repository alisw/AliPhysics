*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:09  mclareni
* Zebra
*
*
*
* fzicvd1.inc
*
#if defined(CERNLIB_QMCDC)
      PARAMETER    (JCHMIN = 31,  JCHMAX = 92)
      DIMENSION    IPKVU(5), IPKVP(5), MTAB(93)
#endif
*
* fzicvd1.inc
*
#if defined(CERNLIB_QMCRY)
      PARAMETER    (JEXMIN = -1023,  JEXMAX=1023)
      PARAMETER    (IOVPMS = 0 60000 40000000 00000000 B)
      PARAMETER    (IOVNMS = 1 60000 40000000 00000000 B)
      PARAMETER    (IBADMS = IOVPMS)
      PARAMETER    (IBADMD = IBADMS, IBADME = 0)
      PARAMETER    (IOVPMD = IOVPMS, IOVNMD = IOVNMS)
      PARAMETER    (IBL32  = 4010020040 B)
#endif
*
* fzicvd1.inc
*
#if (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      PARAMETER    (IBL32 = '20202020'X)
      PARAMETER    (MSKB32= '80000000'X)
      PARAMETER    (MSKL32= 'FFFFFFFF00000000'X)
#endif
*
* fzicvd1.inc
*
#if defined(CERNLIB_QMIBM)
      PARAMETER    (JCHMIN = 31,  JCHMAX = 96)
      DIMENSION    IPKVEC(5), MTAB(97)
#endif
*
* fzicvd1.inc
*
#if defined(CERNLIB_QMIBMD)
      PARAMETER    (JCHMIN = 31,  JCHMAX = 96)
      DIMENSION    IPKVEC(5), MTAB(97)
#endif
