*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:20:20  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:46  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
*
* mzeunit.inc
*
#if defined(CERNLIB_ALLIANT)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_APOLLO)
      IQREAD = 5
      IQPRNT = 6
      IQPNCH = 3

#elif defined(CERNLIB_CDC)
      IQREAD = L"INPUT"
      IQPRNT = L"OUTPUT"
      IQPNCH = L"PUNCH"

#elif (defined(CERNLIB_CRAY))&&!defined(CERNLIB_QMCRU)
      IQREAD = 3L$IN
      IQPRNT = 4L$OUT
      IQPNCH = 6L$PUNCH

#elif defined(CERNLIB_DGE)
      IQREAD = 1
      IQPRNT = 2
      IQPNCH = 3

#elif defined(CERNLIB_MSDOS) || defined(CERNLIB_WINNT)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_HPUX)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_IBMAIX)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_IBMRT)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_IBMRTD)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_LINUX)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4


#elif defined(CERNLIB_MACMPW)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4


#elif defined(CERNLIB_NEXT)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_OS9)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_PDP)
      IQREAD = 1
      IQPRNT = 3
      IQPNCH = 20

#elif defined(CERNLIB_SGI)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_SUN)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_TMO)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_VAXVMS)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_DECS)
      IQREAD = 2
      IQPRNT = 3
      IQPNCH = 4

#elif defined(CERNLIB_QMVDS)
      IQREAD = 5
      IQPRNT = 6
      IQPNCH = 3

#elif 1
      IQREAD = 5
      IQPRNT = 6
      IQPNCH = 7

#endif
