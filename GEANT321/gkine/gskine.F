*
* $Id$
*
* $Log$
* Revision 1.2  1998/02/10 16:04:19  japost
*   Comments are correctly redistributed.
*
* Revision 1.1.1.1  1995/10/24 10:21:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 28/03/94  01.30.59  by  S.Giani
*-- Author :
      SUBROUTINE GSKINE(PLAB,IPART,NV,BUF,NWBUF,NT)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Store kinematics of track NT into data structure         *
C     *       Track is coming from vertex NV                           *
C.    *                                                                *
C.    *                                                                *
C.    * In detail:                                                     *
C.    *                                                                *
C.    *  Stores long life track parameters.                            *
C.    *                                                                *
C.    *  Input:                                                        *
C.    *    PLAB      components of momentum                            *
C.    *    IPART     type of particle (see GSPART)                     *
C.    *    NV        vertex number origin of track                     *
C.    *    UBUF      array of NUBUF floating point user parameters     *
C.    *    NUBUF                                                       *
C.    *                                                                *
C.    *  Output:                                                       *
C.    *    NT        track number (if=0 error).                        *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUKINE,GUSTEP                        *
C.    *       Authors    R.Brun, F.Carena, M.Hansroul  *********       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcpush.inc"
      DIMENSION PLAB(3),BUF(1)
      SAVE JV
C.
C.    ------------------------------------------------------------------
C.
      IF(NV.LE.0)GO TO 10
      IF(NV.GT.NVERTX)GO TO 90
      JV = LQ(JVERTX- NV)
      IF(JV.LE.0)GO TO 90
C
   10 NT     = NTRACK + 1
      IF (JKINE.EQ.0)THEN
         CALL MZBOOK(IXDIV,JKINE,JKINE,1,'KINE',NCKINE,NCKINE,1,2,0)
         IQ(JKINE-5)=0
      ENDIF
      IF(NT.GT.IQ(JKINE-2)) CALL MZPUSH(IXDIV,JKINE,NPKINE,0,'I')
C
      CALL MZBOOK(IXDIV,JK,JKINE,-NT,'KINE',1,1,8,3,0)
C
      IF (JPART.LE.0) GO TO 90
      IF (IPART.GT.IQ(JPART-2)) GO TO 90
      JPA = LQ(JPART- IPART)
      IF (JPA.LE.0) GO TO 90
      P2     = PLAB(1)**2 + PLAB(2)**2 + PLAB(3)**2
      AMASS = Q(JPA + 7)
      E      = SQRT(P2 + AMASS*AMASS)
C
C            Store kinematics in JKINE
C
      Q(JK + 1) = PLAB(1)
      Q(JK + 2) = PLAB(2)
      Q(JK + 3) = PLAB(3)
      Q(JK + 4) = E
      Q(JK + 5) = IPART
      Q(JK + 6) = NV
C
C            Copy user words
C
      IF(NWBUF.GT.0)THEN
         CALL MZBOOK(IXDIV,JUK,JK,-1,'KINU',0,0,NWBUF,3,0)
         IQ(JUK-5)=NT
C
         DO 15 I=1,NWBUF
            Q(JUK+I)=BUF(I)
   15    CONTINUE
      ENDIF
C
C
C            Connect track NT to vertex NV
C
      IF(NV.GT.0)THEN
         JV = LQ(JVERTX- NV)
         NTG = Q(JV + 7)
         NFREE=IQ(JV-1)-7-NTG
         IF(NFREE.LE.0)CALL MZPUSH(IXDIV,JV,0,2,'I')
         Q(JV + NTG + 8) = NT
         Q(JV + 7) = NTG + 1
      ENDIF
C
      NTRACK = NT
      IQ(JKINE+1)=NTRACK
      GO TO 99
C
C            Error
C
  90  NT     = 0
  99  RETURN
      END
