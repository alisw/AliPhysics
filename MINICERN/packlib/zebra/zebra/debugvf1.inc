*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:56  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_DEBUGVF1_INC
#define CERNLIB_ZEBRA_DEBUGVF1_INC
*
*   VFORT (Put just in front of SUBROUTINE st.)
*
* debugvf1.inc
*
#if defined(CERNLIB_QMIBM)
      DEBUG TRACE,SUBTRACE,INIT
      AT 1001
      TRACE ON
      END DEBUG
#endif

#endif
