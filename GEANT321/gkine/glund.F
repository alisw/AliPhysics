*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GLUND
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        GEANT INTERFACE TO THE LUND MONTE-CARLO                 *
C.    *                                                                *
C.    *  FOR TECHNICAL DETAILS OF THE LUND MONTE CARLO,SEE             *
C.    *  1. THE LUND MONTE CARLO FOR JET FRAGMENTATION - T.SJOSTRAND   *
C.    *     LU TP 82-3  and COMPUTER PHYSICS COMMUNICATION 27(1982)243 *
C.    *  2. THE LUND MONTE CARLO FOR E+E- JET PHYSICS  - T.SJOSTRAND   *
C.    *     LU TP 82-7 and COMPUTER PHYSICS COMMUNICATION 28(1983)229  *
C.    *                                                                *
C.    *    The parameters to use the LUND interface are set in GLUNDI  *
C.    *    They can be changed by data cards (see GLUNDI)              *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUKINE                               *
C.    *       Author    R.Hemingway  *********                         *
C.    *       modified by R.Brun  18/11/85 **                          *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION VERT(3)
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcpush.inc"
#include "geant321/lujets.inc"
#include "geant321/gclund.inc"
      PARAMETER (NSEL=44)
      DIMENSION IPLUND(NSEL)
      REAL      PGEANT(3),BUF(1)
      EQUIVALENCE (BUF,K)
      DATA IPLUND/
     +    22,   -11,    11,    12,   -13,    13,   111,   211,
     +  -211,   130,   321,  -321,  2112,  2212, -2212,   310,
     +   221,  3122,  3222,  3212,  3112,  3322,  3312,  3334,
     + -2112, -3122, -3112, -3212, -3222, -3322, -3312, -3334,
     +   -15,    15,   411,  -411,   421,  -421,   431,  -431,
     +  4122,    24,   -24,    23/
*.
*.    ------------------------------------------------------------------
*.
      IF (IFINIT(3).EQ.0)  CALL GLUNDI
*
      CALL LUEEVT(IFLUND,ECLUND)
*
      NCVERT = 10
      NCKINE = N+10
      NCJXYZ = NCKINE
      NPVERT =  5
      NPKINE = 10
      NPJXYZ = 10
*
      VERT(1)   = 0.
      VERT(2)   = 0.
      VERT(3)   = 0.
      NWBUF=2+10*N
      CALL GSVERT(V,0,0,BUF,NWBUF,NVTX)
      IF (NVTX.GT.0) THEN
         DO 20 I = 1,N
            IF(K(I,1).GE.1.AND.K(I,1).LE.10) THEN
               DO 10 IGE=1, NSEL
                  IF(IPLUND(IGE).EQ.K(I,2)) THEN
                     PGEANT(1)=P(I,1)
                     PGEANT(2)=P(I,2)
                     PGEANT(3)=P(I,3)
                     CALL GSKINE(PGEANT,IGE,NVTX,BUF,0,NT)
                     GOTO 20
                  ENDIF
   10          CONTINUE
            ENDIF
   20    CONTINUE
*
         JV = LQ(JVERTX-1)
         JUV = LQ(JV-1)
         Q(JUV+1)= 527.0
         Q(JUV+2)= N
*
*     Copy the LUND record into a bank of real variables.
*
         DO 30 NP=1, N
            Q(JUV+2    +NP) = K(NP,1)
            Q(JUV+2+  N+NP) = K(NP,2)
            Q(JUV+2+2*N+NP) = K(NP,3)
            Q(JUV+2+3*N+NP) = K(NP,4)
            Q(JUV+2+4*N+NP) = K(NP,5)
            Q(JUV+2+5*N+NP) = P(NP,1)
            Q(JUV+2+6*N+NP) = P(NP,2)
            Q(JUV+2+7*N+NP) = P(NP,3)
            Q(JUV+2+8*N+NP) = P(NP,4)
            Q(JUV+2+9*N+NP) = P(NP,5)
   30    CONTINUE
*
      ENDIF
      END
