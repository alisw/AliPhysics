*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*              EoF on direct-access file
*
* fzdaeof.inc
*
#if defined(CERNLIB_QMAPO)
      PARAMETER    (IOSEOF = 84148232)
#elif defined(CERNLIB_QMCRY)
      PARAMETER    (IOSEOF = 5)
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      PARAMETER    (IOSEOF = 134)
#elif (defined(CERNLIB_QMCV32))&&(!defined(CERNLIB_QIEEE))
      PARAMETER    (IOSEOF = 134)
#elif (defined(CERNLIB_QMCV32))&&(defined(CERNLIB_QIEEE))
      PARAMETER    (IOSEOF = 134)
#elif defined(CERNLIB_QMHPX)
      PARAMETER    (IOSEOF = 922)
#elif defined(CERNLIB_QMIRT)
***   PARAMETER    (IOSEOF = 10)     !*  for xlf 2.2
***   PARAMETER    (IOSEOF = 122)    !*  for xlf 2.3
      PARAMETER    (IOSEOF = 1)      !*  for xlf 3.2
#elif defined(CERNLIB_QMIRTD)
***   PARAMETER    (IOSEOF = 10)     !*  for xlf 2.2
***   PARAMETER    (IOSEOF = 122)    !*  for xlf 2.3
      PARAMETER    (IOSEOF = 1)      !*  for xlf 3.2
#elif defined(CERNLIB_QMLNX)
      PARAMETER    (IOSEOF = 36)
#elif defined(CERNLIB_QMOS9)
      PARAMETER    (IOSEOF = 326)
#elif defined(CERNLIB_QMSGI)
      PARAMETER    (IOSEOF = 27)
#elif defined(CERNLIB_QMVAX)
      PARAMETER    (IOSEOF = 36)
#elif defined(CERNLIB_QMVMI)
      PARAMETER    (IOSEOF = 36)
#elif defined(CERNLIB_WINDOWS)
*    only valid for dec fortran (e.g. not f2c + microsoft c)
      PARAMETER    (IOSEOF = 36)
#elif 1
      PARAMETER    (IOSEOF = -1)
#endif
