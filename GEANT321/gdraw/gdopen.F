*
* $Id$
*
* $Log$
* Revision 1.2  1997/01/24 14:16:06  mclareni
* NZLEFT 2nd argument 0 instead of blank for Msoft
*
* Revision 1.1.1.1.2.1  1997/01/21 11:24:20  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDOPEN(IVIEW)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Open view bank IVIEW for storing; if IVIEW=0 then        *
C.    *       none view bank is opened and drawing is done             *
C.    *       on the screen.                                           *
C.    *                                                                *
C.    *       Q(JDRAW+IVIEW) contains a control word that is :         *
C.    *                                                                *
C.    *       1 for empty banks (created just to avoid gaps)           *
C.    *         or for deleted banks                                   *
C.    *                                                                *
C.    *       2 for opened banks (all banks that have been opened)     *
C.    *                                                                *
C.    *       3 for protected banks (all banks that can't be deleted   *
C.    *         by the user)                                           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Authors : R.Brun, P.Zanarini   *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
      COMMON/GCLOCA/NLOCAL(2),JLOCA(19),JV
C.
C.    ------------------------------------------------------------------
C.
C
C             Book JDRAW structure for view banks if not there
C
      IF (JDRAW.LE.0) CALL MZBOOK(IXCONS,JDRAW,JDRAW,1,'DRAW',0,0,0,3,0)
C
C             Previous view number must be 0 (screen)
C
      IF (IDVIEW.NE.0) GO TO 60
C
C             If want open a view bigger than current number of views
C             then do a pushing of view banks
C
      IF(IVIEW.LE.0)GO TO 999
      IF (IVIEW.GT.NKVIEW) GO TO 20
C
C             Cannot open a view already used
C
      IF (LQ(JDRAW-IVIEW).NE.0) GO TO 40
      IDVIEW=IVIEW
      IF (IDVIEW.EQ.0) GO TO 999
C
C             Open view number IDVIEW
C
   10 Q(JDRAW+IDVIEW)=2.
      CALL MZBOOK(IXCONS,JV,JDRAW,-IDVIEW,'DRAW',6,6,22,3,0)
C
C             Set number of words for creation and for pushing of banks
C             graphic unit, graphic segment, text unit, text segment
C
C
      NLEFT=NZLEFT(IXCONS,0)
      NL=NLEFT/10
      MAXGS=MIN(NL,10000)
      MAXGU=MAXGS/10
      MORGS=MAXGS/2
      MORGU=MAXGU/2
      MAXTU=10
      MORTU=MAXTU
      MAXTS=100
      MORTS=MAXTS
C
      Q(JV+2)=MAXGU
      Q(JV+3)=MORGU
      Q(JV+5)=MAXGS
      Q(JV+6)=MORGS
      Q(JV+8)=MAXTU
      Q(JV+9)=MORTU
      Q(JV+11)=MAXTS
      Q(JV+12)=MORTS
C
      Q(JV+16)=PLTRNX/2.
      Q(JV+17)=PLTRNY/2.
C
C           Set ICUT of the view bank to 0
C
      Q(JV+22)=0.
C
C             Graphic unit, graphic segment, text unit, text segment
C             pointers are reset to zero
C
      IGU=0
      IGS=0
      ITU=0
      ITS=0
C
C             Create view banks
C
      ND1=Q(JV+2)
      ND2=ND1
      ND3=Q(JV+8)
      ND4=Q(JV+5)
      ND5=ND4
      ND6=Q(JV+11)
      CALL MZBOOK(IXCONS,LBANK,JV,-1,'DRAW',0,0,ND1,3,0)
      CALL MZBOOK(IXCONS,LBANK,JV,-2,'DRAW',0,0,ND2,3,0)
      CALL MZBOOK(IXCONS,LBANK,JV,-3,'DRAW',0,0,ND3,3,0)
      CALL MZBOOK(IXCONS,LBANK,JV,-4,'DRAW',0,0,ND4,3,0)
      CALL MZBOOK(IXCONS,LBANK,JV,-5,'DRAW',0,0,ND5,3,0)
      CALL MZBOOK(IXCONS,LBANK,JV,-6,'DRAW',0,0,ND6,3,0)
      GO TO 999
C
C             Push view banks
C
   20 IDVIEW=IVIEW
      NDL=IDVIEW-NKVIEW
      CALL MZPUSH(IXCONS,JDRAW,NDL,NDL,'I')
      IF(IEOTRI.NE.0)GO TO 50
      DO 30 I=NKVIEW+1,IDVIEW-1
         Q(JDRAW+I)=1.
   30 CONTINUE
      NKVIEW=IDVIEW
      GO TO 10
C
   40 WRITE (CHMAIL,1000)
      CALL GMAIL(0,0)
      GO TO 999
C
   50 WRITE (CHMAIL,1100)
      CALL GMAIL(0,0)
      GO TO 999
C
   60 WRITE (CHMAIL,1200)
      CALL GMAIL(0,0)
C
 1000 FORMAT (' ERROR IN GDOPEN: VIEW ALREADY EXISTING')
 1100 FORMAT (' ERROR IN GDOPEN: MEMORY OVERFLOW IN PUSHING A BANK')
 1200 FORMAT (' ERROR IN GDOPEN: ANOTHER VIEW IS ALREADY OPENED')
  999 RETURN
      END
