*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GGORDQ (IVO)
C.
C.    *****************************************************************
C.    *                                                               *
C.    *    Find and order the boundaries of the contents of the       *
C.    *    IVOth volume, with respect to coordinate IAX :             *
C.    *           IAX = 1    X Axis                                   *
C.    *           IAX = 2    Y Axis                                   *
C.    *           IAX = 3    Z Axis                                   *
C.    *           IAX = 4    Rxy                                      *
C.    *           IAX = 5    Rxyz                                     *
C.    *           IAX = 6    PHI   (PHI=0 => X axis)                  *
C.    *           IAX = 7    THETA (THETA=0 => Z axis)                *
C.    *    All values of IAX will be tried and then that value is     *
C.    *    chosen, that results in the smallest number of volumes per *
C.    *    division.                                                  *
C.    *    Called by : GGCLOS                                         *
C.    *    Author: Stephan Egli (large parts are copies of GGORD)     *
C.    *                                                               *
C.    *****************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
*
 
      DIMENSION CLOW(500),CHIGH(500),CORD(1000),ITYPE(1000),
     +ICONT(500),ICON(1000),ICONS(500)
      EQUIVALENCE (CLOW(1),WS(1)),(CHIGH(1),WS(501))
      EQUIVALENCE (CORD(1),WS(1001)),(ITYPE(1),WS(2001))
      EQUIVALENCE (ICONT(1),WS(3001)),(ICON(1),WS(3501))
      EQUIVALENCE (ICONS(1),WS(4501))
C
      CHARACTER*4 NAME
 
C.    ------------------------------------------------------------------
*
      JVO = LQ(JVOLUM-IVO)
      CALL UHTOC(IQ(JVOLUM+IVO),4,NAME,4)
      NIN = Q(JVO+3)
      IAXNOW = 0
      IF(IQ(JVO-2).EQ.NIN+2) THEN
*
* *** This is to allow re-entry in the routine from the interactive
* *** version.
         JNEAR = LQ(JVO-NIN-1)
         IF(JNEAR.GT.0) THEN
            JSB = LQ(JNEAR)
            IF(JSB.GT.0) THEN
               IAXNOW = Q(JSB+1)
            ENDIF
         ENDIF
      ENDIF
      IF(IAXNOW.EQ.0) THEN
         IAXNOW=-Q(JVO+1)
      ENDIF
*   assume that ordering can not be done unless proven otherwise
      Q(JVO+1)=0.
      RBEST=1.E9
 
* try all possible axes
 
      DO 1 IAX=1,7
 
*   count number of additional words needed and total number of volumes
*   in all divisions
 
      NCOALL=0
*
* *** Find the upper and lower coordinates of each content
*
      DO 50 IN = 1,NIN
         CALL GFCLIM (JVO, IN, IAX, CLOW(IN), CHIGH(IN), IERR)
         IF (IERR.NE.0) GOTO 1
   50 CONTINUE
*
* *** Order the coordinate limits, keeping track of the associated
*           content number
*
      CALL GFCORD (NIN, CLOW, CHIGH, CORD, ITYPE, ICON)
      NC = NIN*2
*
*  **   Count and load up the distinct boundaries
*
      IBO = 0
      DO 60 IC = 1,NC
         IBO = IBO +1
         IF (IBO.EQ.1) GO TO 60
         IF (CORD(IC)-CORD(IC-1).LT.1.E-4) IBO = IBO -1
   60 CONTINUE
      NDIV  = IBO -1
      IF (IAX.EQ.6) NDIV = IBO
 
*   *   Load up number of contents in each section
*
      IDIV    = 0
      NCONT   = 1
      ICONT(1)= ICON(1)
      IF (IAX.NE.6) GO TO 70
      NCONT   = 0
      NSTOR   = 0
      ICONT(1)= 0
      DO 65 IN = 1,NIN
         IF (CHIGH(IN).GT.CLOW(IN)) GO TO 65
*           (this content straddles PHI=0.)
         NSTOR = NSTOR +1
         ICONS(NSTOR) = IN
         IF (ICON(1).EQ.IN) GO TO 65
*           (IN is in 1st division as well)
         NCONT = NCONT +1
         ICONT(NCONT) = IN
   65 CONTINUE
*
      IF (ITYPE(1).EQ.2) GO TO 70
*            (first boundary is a low, add the new content)
      NCONT = NCONT +1
      ICONT(NCONT) = ICON(1)
*
   70 CONTINUE
*
      DO 130 IC = 2,NC
         IDIV = IDIV +1
         IF (CORD(IC)-CORD(IC-1).LT.1.E-4) GO TO 90
*
*          New division, load up last division
*
         IF (NCONT.LE.0) GO TO 100
         NCOALL=NCOALL+NCONT
         GO TO 100
   90    CONTINUE
         IDIV = IDIV -1
*
  100    CONTINUE
*
*         Update contents of current division
*
         IF (ITYPE(IC).EQ.1) GO TO 120
*
*         This boundary was a high, so one less content
*
         ICP = 0
         DO 110 ICNT = 1,NCONT
            IF (ICONT(ICNT).EQ.ICON(IC)) ICP=1
         IF (ICP.EQ.1) ICONT(ICNT) = ICONT(ICNT+1)
  110    CONTINUE
         NCONT = NCONT -1
         GO TO 130
*
  120    CONTINUE
*
*          This boundary was a low, so one extra content
*
         NCONT = NCONT +1
         ICONT(NCONT) = ICON(IC)
*
  130 CONTINUE
*
      IF(IAX.EQ.6) NCOALL = NCOALL+NSTOR
      RNOW=FLOAT(NCOALL)/NDIV
      IF(RNOW.LT.RBEST)THEN
        IAXOPT=IAX
        RBEST=RNOW
        NDIVB=NDIV
      ENDIF
 
* end of loop over IAX
 
1     CONTINUE
 
 
* now the best axis is selected - compare with axis requested by CALL
* to GSORD (if any)
 
      IF(IAXNOW.GT.0)THEN
 
        WRITE (CHMAIL,1002) NAME,NIN,IAXOPT,NDIVB,RBEST,IAXNOW
        CALL  GMAIL (0, 0)
 1002   FORMAT(' GGORDQ : Volume ',A4,2X,'NIN=',I4,' IAX=',I2,2X,
     +    'NDIV=',I3,2X,'NVOL/DIV=',F5.1,2X,'IAX wanted by user:',I2)
 
      ELSE
 
        WRITE (CHMAIL,1003) NAME,NIN,IAXOPT,NDIVB,RBEST
        CALL  GMAIL (0, 0)
 1003   FORMAT(' GGORDQ : Volume ',A4,2X,'NIN=',I4,' IAX=',I2,2X,
     +    'NDIV=',I3,2X,'NVOL/DIV=',F5.1)
 
      ENDIF
 
* overwrite old axis and store sorting information for new axis
 
      Q(JVO+1)=-IAXOPT
      CALL GGORD(IVO)
 
      END
