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
      SUBROUTINE GSNEXT (MOTHER,IN,NLIST,LIST)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Gives an ordered LIST of NLIST daughter volumes to search  *
C.    *     for when leaving the INth daughter of the volume MOTHER    *
C.    *                                                                *
C.    *     Default facility :                                         *
C.    *         If IN = 0, for each content in turn, GSNEXT builds a   *
C.    *         LIST  limited to next (if any), previous (if any) and  *
C.    *         current content                                        *
C.    *                                                                *
C.    *     Called by : 'USER'                                         *
C.    *     Author    : F.Bruyant  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      CHARACTER*4 MOTHER
      DIMENSION LIST(*)
C.    ------------------------------------------------------------------
*
      IVO = 0
      IF (JVOLUM.GT.0) CALL GLOOK (MOTHER, IQ(JVOLUM+1), NVOLUM, IVO)
      IF (IVO.EQ.0) THEN
         WRITE (CHMAIL, 1001)
         CALL GMAIL (0, 0)
         GO TO 999
      ENDIF
*
* *** Check that volume MOTHER has contents positioned by GSPOS or
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
         NVOLST = 0
         DO 19 I = 1,NLIST
            IF(LIST(I).GT.NIN.OR.LIST(I).LT.0) THEN
               WRITE(CHMAIL,1003) IQ(JVOLUM+IVO),I,LIST(I),NIN
               CALL GMAIL(0,0)
            ELSE
               DO 20 J=1, NVOLST
                  IF(IQ(JNEAR+1+J).EQ.LIST(I)) THEN
                     WRITE(CHMAIL,1004) IQ(JVOLUM+IVO), LIST(I)
                     CALL GMAIL(0,0)
                     GOTO 19
                  ENDIF
   20          CONTINUE
               NVOLST = NVOLST+1
               IQ(JNEAR+NVOLST+1) = LIST(I)
            ENDIF
   19    CONTINUE
         IQ(JNEAR+1) = NVOLST
         IF(NVOLST.LT.NLIST) THEN
             CALL MZPUSH (IXCONS, JNEAR, 0, NVOLST-NLIST, 'R')
         ENDIF
*
      ELSE IF (IN.EQ.0) THEN
*
* ***    Default option for all contents
*
         DO 49 I = 1,NIN
            JIN = LQ(JVO-I)
            IF (LQ(JIN-1).NE.0) GO TO 49
            N = 3
            IF (I.EQ.1.OR.I.EQ.NIN) N = 2
            CALL MZBOOK (IXCONS, JNEAR,JIN,-1,'VONE',0,0,N+1, 2,0)
            IQ(JNEAR-5) = 100*IVO+I
            JVO = LQ(JVOLUM-IVO)
            IQ(JNEAR+1)  = N
            IQ(JNEAR+N+1) = I
            IF (I.NE.1) GO TO 41
            IQ(JNEAR+N) = I +1
            GO TO 49
   41       IF (I.NE.NIN) GO TO 42
            IQ(JNEAR+N) = I -1
            GO TO 49
   42       IQ(JNEAR+N)   = I +1
            IQ(JNEAR+N-1) = I -1
   49    CONTINUE
*
      ENDIF
*
 1001 FORMAT (' GSNEXT : Volume',A4,' does not exist  *****')
 1002 FORMAT (' GSNEXT : Volume',A4,' not a candidate for GSNEXT *****')
 1003 FORMAT (' GSNEXT : Volume ',A4,': LIST(',I3,')=',I3,
     +        ' NIN=',I3,' value ignored')
 1004 FORMAT (' GSNEXT : Volume',A4,': repetition of content ',
     +        I3,' ignored')
*                                                             END GSNEXT
  999 END
