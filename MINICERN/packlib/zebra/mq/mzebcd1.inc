*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:20:19  mclareni
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
* mzebcd1.inc
*
#if defined(CERNLIB_ALLIANT)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_CDC)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;:'

#elif defined(CERNLIB_MSDOS) || defined(CERNLIB_WINNT)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_HPUX)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_IBMAIX)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_IBMRT)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'

#elif defined(CERNLIB_IBMRTD)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'

#elif defined(CERNLIB_LINUX)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_MACMPW)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_NEXT)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_OS9)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_SGI)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_SUN)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_TMO)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_DECS)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif defined(CERNLIB_QMVDS)
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(61:61) = CHAR(92)

#elif 1
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'

#endif
