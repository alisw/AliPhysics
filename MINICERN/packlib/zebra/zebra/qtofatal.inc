*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:54  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_QTOFATAL_INC
#define CERNLIB_ZEBRA_QTOFATAL_INC
*
*
* qtofatal.inc
*
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      IQUEST(10)= NAMESR
      CALL ZFATAL
#elif defined(CERNLIB_QTRHOLL)
      IQUEST(9) = NAMESR(1)
      IQUEST(10)= NAMESR(2)
      CALL ZFATAL
#elif defined(CERNLIB_A6M)
      CALL UCTOH (NAMESR,IQUEST(10),6,6)
      CALL ZFATAL
#elif 1 
      CALL UCTOH (NAMESR,IQUEST(9),4,8)
      CALL ZFATAL
*  -------------   for ZEBRA FZ   ------------------------------

#endif
#endif
