*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXLUNF(LUN,ICASE,IFREE)
*
*             If ICASE=1 check if LUN is free
*             If ICASE=2 check if LUN is used
*             If ICASE=3 return first free unit starting at LUN
*
#include "geant321/gcxlun.inc"
#include "geant321/gcunit.inc"
      CHARACTER*8 CHCASE(10)
      SAVE CHCASE
      DATA CHCASE/'FZ      ','FZOUT   ','HBOOK   ','HIGZ    '
     +           ,'RZ      ','KUIP    ','COMIS   ','GKS     '
     +           ,'FORTRAN ','USER    '/
*________________________________________________________
*
      IF(LUN.LE.0.OR.LUN.GT.128)THEN
         PRINT *,' Invalid Logical unit number'
         IFREE=6
         IF(ICASE.NE.1)IFREE=0
         GO TO 99
      ENDIF
      IFREE=IABS(LUNIT(LUN))
      IF(ICASE.EQ.1)THEN
         IF(IFREE.EQ.0)GO TO 99
         IFR=MIN(IFREE,10)
         WRITE(CHMAIL,1000) LUN,CHCASE(IFR)
         CALL GMAIL(0,0)
      ELSEIF(ICASE.EQ.2)THEN
         IF(IFREE.NE.0)GO TO 99
         WRITE(CHMAIL,2000) LUN
         CALL GMAIL(0,0)
      ELSE
         IFREE=0
         DO 10 I=LUN,99
            IF(LUNIT(I).EQ.0)THEN
               IFREE=I
               GO TO 99
            ENDIF
  10     CONTINUE
         WRITE(CHMAIL,3000) LUN
         CALL GMAIL(0,0)
      ENDIF
*
  99  RETURN
 1000 FORMAT(' Unit ',I3,' already used by ',A)
 2000 FORMAT(' Unit ',I3,' is not active')
 3000 FORMAT(' No more free logical units above ',I3)
      END
