*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:09  mclareni
* Zebra
*
*
*
* fznwmach.inc
*
#if defined(CERNLIB_QMCDC)
      NWMIN = (32*NWMIN) / 60
#endif
*
* fznwmach.inc
*
#if defined(CERNLIB_QMCRY)
      NWMIN = NWMIN / 2
#endif
*
* fznwmach.inc
*
#if (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      NWMIN = NWMIN / 2
#endif
*
* fznwmach.inc
*
#if defined(CERNLIB_QMIBMD)
      NWMIN = NWMIN / 2
#endif
*
* fznwmach.inc
*
#if defined(CERNLIB_QMIRTD)
      NWMIN = NWMIN / 2
#endif
