*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:55  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_JZTIMIN_INC
#define CERNLIB_ZEBRA_JZTIMIN_INC
*
*
* jztimin.inc
*
#if defined(CERNLIB_JZTIME)
      IF (JQTIME.NE.0)  THEN
          CALL TIMED (DELTAT)
          Q(KQS+LQSV+5) = Q(KQS+LQSV+5) + DELTAT
        ENDIF
#endif

#endif
