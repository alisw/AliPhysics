*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:57  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_QTOABEND_INC
#define CERNLIB_ZEBRA_QTOABEND_INC
*
* qtoabend.inc
*
#if defined(CERNLIB_IBM)||defined(CERNLIB_IBMD)
*             FORCE LAST LINE TO BE OUT OF THE BUFFER
*
      WRITE (IQLOG,9701)
 9701 FORMAT ('0GOING TO ABEND.'/ 20(1X,20(1H*)/) )
      CALL ABEND
#elif defined(CERNLIB_PDP)
      WRITE (IQLOG, 9701)
 9701 FORMAT ('0GOING TO ABEND.', 4X, 20('*'))
      STOP '********************'
#elif 1
      CALL ABEND
#endif

#endif
