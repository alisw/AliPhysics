*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:09  mclareni
* Zebra
*
*
*   no longer needed
*
* fzicxd1.inc
*
#if defined(CERNLIB_QMCDC)
      PARAMETER    (JEXMIN = -975,  JEXMAX=1023)
      PARAMETER    (IBADMS = O"17770000000000000000")
      PARAMETER    (IOVPMS = O"37774000000000000000")
      PARAMETER    (IOVNMS = O"40003777777777777777")
      PARAMETER    (IBADMD = IBADMS, IBADME = 0)
      PARAMETER    (IOVPMD = IOVPMS, IOVNMD = IOVNMS)
#endif
*   no longer needed
*
* fzicxd1.inc
*
#if defined(CERNLIB_QMIBM)
      PARAMETER    (JEXMIN = -259,  JEXMAX=251)
C     PARAMETER    (IOVPMS =   Z 7FFFFFF0)
      PARAMETER    (IOVPMS = 2147 483 632)
C     PARAMETER    (IOVNMS =   Z FFFFFFF0)
      PARAMETER    (IOVNMS = -16)
      PARAMETER    (IBADMD = IBADMS, IBADME = 0)
      PARAMETER    (IOVPMD = IOVPMS, IOVNMD = IOVNMS)
#endif
*   no longer needed
*
* fzicxd1.inc
*
#if defined(CERNLIB_QMIBMD)
      PARAMETER    (JEXMIN = -259,  JEXMAX=251)
C     PARAMETER    (IOVPMS =   Z 7FFFFFF0)
      PARAMETER    (IOVPMS = 2147 483 632)
C     PARAMETER    (IOVNMS =   Z FFFFFFF0)
      PARAMETER    (IOVNMS = -16)
      PARAMETER    (IBADMD = IBADMS, IBADME = 0)
      PARAMETER    (IOVPMD = IOVPMS, IOVNMD = IOVNMS)
#endif
*   no longer needed
*
* fzicxd1.inc
*
#if defined(CERNLIB_QMND3)
      PARAMETER    (JEXMIN = -255,  JEXMAX=253)
      PARAMETER    (IOVPMS = 17777777770 B)
      PARAMETER    (IOVNMS = 37777777770 B)
      PARAMETER    (IBADMD = IBADMS, IBADME = 0)
      PARAMETER    (IOVPMD = IOVPMS, IOVNMD = IOVNMS)
#endif
*    use default CALL IE3TOS for input single prec.
#ifndef CERNLIB_FQIE3TSC
#endif
*    use default CALL IE3TOD for input double prec.
#ifndef CERNLIB_FQIE3TDC
#endif
*    --------------  sequences for FZMACH  ---------------------
*   no longer needed
*
* fzicxd1.inc
*
#if defined(CERNLIB_QMVAX)
      PARAMETER    (JEXMIN = -127,  JEXMAX=125)
      PARAMETER    (IOVPMS = '00007F81'X)
      PARAMETER    (IOVNMS = '0000FF81'X)
      PARAMETER    (IBADMS = '00008001'X)
      PARAMETER    (IBADMD = IBADMS, IBADME = 0)
      PARAMETER    (IOVPMD = IOVPMS, IOVNMD = IOVNMS)
#endif
*    use default CALL IE3TOS for input single prec.
#ifndef CERNLIB_FQIE3TSC
#endif
*    use default CALL IE3TOD for input double prec.
#ifndef CERNLIB_FQIE3TDC
#endif
