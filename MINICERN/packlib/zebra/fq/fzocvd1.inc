*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*
* fzocvd1.inc
*
#if defined(CERNLIB_QMCDC)
      DIMENSION    IPKVU(5), IPKVP(5), MTAB(64)
#endif
*
* fzocvd1.inc
*
#if defined(CERNLIB_QMIBM)
      PARAMETER    (JCHMIN = 63,  JCHMAX = 187)
      DIMENSION    IPKVEC(5), MTAB(188)
#endif
*
* fzocvd1.inc
*
#if defined(CERNLIB_QMIBMD)
      PARAMETER    (JCHMIN = 63,  JCHMAX = 187)
      DIMENSION    IPKVEC(5), MTAB(188)
#endif
