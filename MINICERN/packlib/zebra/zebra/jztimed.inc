*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:55  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_JZTIMED_INC
#define CERNLIB_ZEBRA_JZTIMED_INC
*
*
* jztimed.inc
*
#if defined(CERNLIB_JZTIME)
      IF (JQTIME.NE.0)  CALL TIMED (DELTAT)
#endif

#endif
