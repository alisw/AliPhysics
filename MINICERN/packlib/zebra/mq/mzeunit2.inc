*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:20:20  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:47  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
*
* mzeunit2.inc
*
#if defined(CERNLIB_ALLIANT)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_APOLLO)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_CONVEX)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_DGE)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_MSDOS) || defined(CERNLIB_WINNT)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_HPUX)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_IBMAIX)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_IBMRT)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_IBMRTD)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_LINUX)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_MACMPW)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_NEXT)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_NORD)
      IQTTIN = 1
      IQTYPE = 1

#elif defined(CERNLIB_OS9)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_PDP)
      IQTTIN = 5
      IQTYPE = 5

#elif defined(CERNLIB_SGI)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_SUN)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_TMO)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_VAXVMS)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_DECS)
      IQTTIN = 5
      IQTYPE = 6

#elif defined(CERNLIB_QMVDS)
      IQTTIN = 1
      IQTYPE = 1

#elif 1
      IQTTIN = 0
      IQTYPE = 0

#endif
