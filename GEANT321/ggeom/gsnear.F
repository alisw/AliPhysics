*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSNEAR (CHMOTH,IN,NLIST,LIST)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Gives an ordered LIST of NLIST daughter volumes to search  *
C.    *     for when leaving the INth daughter of the volume CHMOTH    *
C.    *     (LIST(1) = 0 means back in mother at exit of INth content) *
C.    *                                                                *
C.    *     If IN = -1, the mother is regarded to have no contents     *
C.    *         contiguous to its boundaries (LIST and NLIST unused)   *
C.    *        (Bit 4 set in CHMOTH volume bank for action in GGCLOS)  *
C.    *                                                                *
C.    *     Default facility :                                         *
C.    *         If IN = 0, for each content in turn, GSNEAR builds a   *
C.    *         LIST  limited to 1 element LIST(1)=0                   *
C.    *                                                                *
C.    *     Called by : 'USER'                                         *
C.    *     Author    : F.Bruyant  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      CHARACTER*4 CHMOTH
      DIMENSION LIST(*)
C.    ------------------------------------------------------------------
*
      IVO = 0
      IF (JVOLUM.GT.0) CALL GLOOK (CHMOTH, IQ(JVOLUM+1), NVOLUM, IVO)
      IF (IVO.EQ.0) THEN
         WRITE (CHMAIL, 1001)
         CALL GMAIL (0, 0)
         GO TO 999
      ENDIF
*
* *** Check that volume CHMOTH has contents positioned by GSPOS or
*      GSPOSP, and not ordered by User
*
      JVO    = LQ(JVOLUM-IVO)
      ISEARC = Q(JVO+1)
      NIN    = Q(JVO+3)
      IF (NIN.LE.0.OR.ISEARC.GT.0) THEN
         WRITE (CHMAIL, 1002)
         CALL GMAIL (0, 0)
         GO TO 999
      ENDIF
*
      IF (IN.GT.0) THEN
*
* ***    Lift a JNEAR bank for current content
*
         JIN = LQ(JVO-IN)
         CALL MZBOOK (IXCONS, JNEAR,JIN,-1,'VONE', 0,0,NLIST+1, 2,0)
         IQ(JNEAR-5) = 100*IVO +IN
         IQ(JNEAR+1) = NLIST
         DO 19 I = 1,NLIST
            IQ(JNEAR+I+1) = LIST(I)
   19    CONTINUE
*
      ELSE IF (IN.EQ.0) THEN
*
* ***    Default option for all contents
*
         DO 49 I = 1,NIN
            JIN = LQ(JVO-I)
            CALL MZBOOK (IXCONS, JNEAR,JIN,-1,'VONE',0,0,2, 2,0)
            IQ(JNEAR-5) = 100*IVO+I
            JVO = LQ(JVOLUM-IVO)
            IQ(JNEAR+1)  = 1
   49    CONTINUE
*
      ELSE
*
         IQ(JVO) = IBSET(IQ(JVO),3)
*
      ENDIF
*
 1001 FORMAT (' GSNEAR : Volume',A4,' does not exist  *****')
 1002 FORMAT (' GSNEAR : Volume',A4,' not a candidate for GSNEAR *****')
*                                                             END GSNEAR
  999 END
