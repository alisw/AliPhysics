*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZNDIR(CHPATH,CHOPT)
*
************************************************************************
*
*      To set,read or print the naming directory
* Input:
*   *CHPATH* Character variable  specifying the complete pathname  of the
*            naming directory ('S' option)
*   CHOPT    Character variable specifying the option
*            'R'   Read the naming directory pathname into CHPATH
*            'P'   Print the naming directory
*            ' '   Set  the naming  directory to  the  path specified  in
*                  CHPATH
* Output:
*   *CHPATH* Character variable  containing the complete pathname  of the
*            naming directory (default option).
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 01.04.86
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
      DIMENSION    IOPTV(2)
      EQUIVALENCE (IOPTR,IOPTV(1)), (IOPTP,IOPTV(2))
      CHARACTER*(*) CHPATH,CHOPT
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
      CALL UOPTC (CHOPT,'RP',IOPTV)
*
*            Read Naming directory
*
      IF(IOPTR.NE.0)THEN
         CALL RZPAFF(CHNDIR,NLNDIR,CHPATH)
         GO TO 99
      ENDIF
*
*            Print Naming directory
*
      IF(IOPTP.NE.0)THEN
         CALL RZPAFF(CHNDIR,NLNDIR,CHL)
         WRITE(IQPRNT,1000)CHL(1:90)
 1000 FORMAT(' Naming Directory = ',A)
         GO TO 99
      ENDIF
*
*            Set Naming directory
*
      CALL RZPATH(CHPATH)
      NLNDIR= NLPAT
      DO 10 I=1,NLPAT
         CHNDIR(I)=CHPAT(I)
  10  CONTINUE
*
  99  RETURN
      END
