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
      SUBROUTINE GGORD (IVO)
C.
C.    *****************************************************************
C.    *                                                               *
C.    *    Find and order the boundaries of the contents of the       *
C.    *    IVOth volume, with respect to coordinate IAX :             *
C.    *           IAX = 1    X Axis                                   *
C.    *           IAX = 2    Y Axis                                   *
C.    *           IAX = 3    Z Axis                                   *
C.    *           IAX = 4    Rxy (static ordering only  -> GTMEDI)    *
C.    *           IAX = 14   Rxy (also dynamic ordering -> GTNEXT)    *
C.    *           IAX = 5    Rxyz (static ordering only -> GTMEDI)    *
C.    *           IAX = 15   Rxyz (also dynamic ordering -> GTNEXT)   *
C.    *           IAX = 6    PHI   (PHI=0 => X axis)                  *
C.    *           IAX = 7    THETA (THETA=0 => Z axis)                *
C.    *                                                               *
C.    *    Called by : GGCLOS                                         *
C.    *    Authors   : R.Brun, F.Bruyant, A.McPherson  *********      *
C.    *                                                               *
C.    *****************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      DIMENSION CLOW(500),CHIGH(500),CORD(1000),ITYPE(1000),
     +ICONT(500),ICON(1000),ICONS(500)
      DIMENSION P(100)
      EQUIVALENCE (CLOW(1),WS(1)),(CHIGH(1),WS(501))
      EQUIVALENCE (CORD(1),WS(1001)),(ITYPE(1),WS(2001))
      EQUIVALENCE (ICONT(1),WS(3001)),(ICON(1),WS(3501))
      EQUIVALENCE (ICONS(1),WS(4501)),(P(1),WS(5001))
C.
C.    ------------------------------------------------------------------
*
      JVO = LQ(JVOLUM-IVO)
      NIN = Q(JVO+3)
      IAX = -Q(JVO+1)
*
* *** Inhibit dynamic ordering in Rxy and R unless explicitly required
*     by user (ordering axis 14 or 15)
*
      IF(IAX.EQ.4.OR.IAX.EQ.5)THEN
*
* *** Only static ordering allowed, dynamic inhibited (GNEXT,GTNEXT)
*
         Q(JVO+1)=-1.
      ELSE
         Q(JVO+1)=-2.
      ENDIF
      IAX=MOD(IAX,10)
*
* *** Find the upper and lower coordinates of each content
*
      DO 50 IN = 1,NIN
         CALL GFCLIM (JVO, IN, IAX, CLOW(IN), CHIGH(IN), IERR)
         IF (IERR.NE.0) THEN
            JIN = LQ(JVO-IN)
            IVOT  = Q(JIN+2)
            IROTT = Q(JIN+4)
            INUMT = Q(JIN+3)
            Q(JVO+1) = 0.
            WRITE (CHMAIL, 1001) IQ(JVOLUM+IVOT), INUMT, IQ(JVOLUM+IVO),
     +                           IAX, IROTT
            CALL GMAIL (0, 0)
            GO TO 999
         ENDIF
   50 CONTINUE
*
* *** Order the coordinate limits, keeping track of the associated
*           content number
*
      CALL GFCORD (NIN, CLOW, CHIGH, CORD, ITYPE, ICON)
*
* *** Book the JSB bank to store the boundaries
*
      JNEAR = LQ(JVO-NIN-1)
      IF(LQ(JNEAR).EQ.0) THEN
         CALL MZBOOK (IXCONS,JSB,JNEAR,0,'VOBO',0,0,NIN*2+2,3,0)
         JVO = LQ(JVOLUM-IVO)
      ELSE
         JSB = LQ(JNEAR)
      ENDIF
      IQ(JSB-5) = IVO
*
      Q(JSB+1) = IAX
      NC = NIN*2
*
*  **   Count and load up the distinct boundaries
*
      IBO = 0
      DO 60 IC = 1,NC
         IBO = IBO +1
         Q(JSB+IBO+2) = CORD(IC)
         IF (IBO.EQ.1) GO TO 60
         IF (CORD(IC)-CORD(IC-1).LT.1.E-4) IBO = IBO -1
   60 CONTINUE
      Q(JSB+2) = IBO
      NDIV  = IBO -1
      IF (IAX.EQ.6) NDIV = IBO
*
*  **   Book the JSC0 bank to store the number of contents in each
*        section (between neighbouring boundaries)
*
      JSC0 = LQ(JVO-NIN-2)
      IF(JSC0.GT.0) THEN
         CALL MZDROP(IXCONS,JSC0,'L')
      ENDIF
      CALL MZBOOK (IXCONS,JSC0,JVO,-NIN-2,'VOBC',NDIV,NDIV,NDIV,2,0)
      IQ(JSC0-5) = IVO
*
*   *   Load up number of contents in each section and when greater
*        than 0 book and load bank of contents
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
         IQ(JSC0+IDIV) = NCONT
         IF (NCONT.LE.0) GO TO 100
*
*          Book bank for contents
*
         CALL MZBOOK (IXCONS,JSCV,JSC0,-IDIV,'VODC',0,0,NCONT,2,0)
         JVO = LQ(JVOLUM-IVO)
         JSC0= LQ(JVO-NIN-2)
*
*          Load up contents
*
         DO 80 ICNT = 1,NCONT
            IQ(JSCV+ICNT) = ICONT(ICNT)
   80    CONTINUE
*
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
      IF (IAX.NE.6) GO TO 150
      IQ(JSC0+NDIV) = NSTOR
      IF (NSTOR.EQ.0) GO TO 150
      CALL MZBOOK (IXCONS,JSCV,JSC0,-NDIV,'VOID',0,0,NSTOR,2,0)
*
      DO 140 IS = 1,NSTOR
         IQ(JSCV+IS) = ICONS(IS)
  140 CONTINUE
      JVO  = LQ(JVOLUM-IVO)
      JSC0 = LQ(JVO-NIN-2)
  150 CONTINUE
      DO 159 I = 1,NDIV
         IF (IQ(JSC0+I).GT.1) GO TO 999
  159 CONTINUE
      IQ(JSC0) = IBSET(IQ(JSC0),0)
*
 1001 FORMAT (' GGORD : Error in GFCLIM for content ',A4,I7,' in ',A4,
     +         ' along axis',I5,' IROT= ',I5)
*                                                              END GGORD
  999 END
 
