*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE XSECN3(KM,KE,RHO,IN,IDICTS,LDICT,ISIGTS,LSIGT,BUF,
     +IBUF,TCS,LIM,LAST)
C       THIS ROUTINE CREATES MACROSCOPIC TOTAL CROSS SECTIONS
C       AND THEN MIXES AND THINS THESE CROSS SECTIONS ACCORDING
C       TO THE MIXING TABLE
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mpoint.inc"
#include "geant321/mmicab.inc"
      DIMENSION BUF(*),IBUF(*),KM(*),KE(*),RHO(*),IN(*),
     +IDICTS(NNR,NNUC),LDICT(NNR,NNUC),ISIGTS(*),LSIGT(*),TCS(*)
C       ASSIGN THE INITIAL VALUES
C       LST EQUALS THE LAST POSITION USED IN THE BUF ARRAY
C   (I.E. (BUF(LST) = D(LAST)))
C       LEN EQUALS THE CORE SPACE AVAILABLE
      LST=0
      LEN=LIM-LAST
      TOL = 1.0
C       READ IN TWO CROSS SECTION ARRAYS AND CREATE
C       MACROSCOPIC CROSS SECTIONS
      DO 160 J=1,MEDIA
         JI=0
         K=0
C       READ IN THE FIRST ARRAY
         DO 140 IJ=1,NMIX
            IF(KM(IJ).NE.J)GO TO 140
            JI=JI+1
            K=K+1
            II=IN(IJ)
            TOL = AMIN1(TCS(LFP210+II-1)/5.,TOL)
            IF(JI.EQ.2)GO TO 20
            LZ=LDICT(1,II)
            ISLZ=IDICTS(1,II)+LMOX2
            N=LZ
            IF(LEN.LT.N)GO TO 180
            NP=LZ/2
            DO 10 M=1,NP
               BUF(LST+2*M-1)=TCS(ISLZ+2*(M-1))
               BUF(LST+2*M)=TCS(ISLZ+2*M-1)*RHO(IJ)
   10       CONTINUE
            GO TO 140
   20       CONTINUE
C       READ IN THE SECOND ARRAY
            LZ2=LZ+1
            LZ1=LZ
            LZ=LDICT(1,II)
            ISLZ=IDICTS(1,II)+LMOX2
            N=2*(LZ+LZ1)
            IF(N.GE.LEN)GO TO 180
            NP=LZ/2
            DO 30 M=1,NP
               BUF(LST+LZ1+2*M-1)=TCS(ISLZ+2*(M-1))
               BUF(LST+LZ1+2*M)=TCS(ISLZ+2*M-1)*RHO(IJ)
   30       CONTINUE
            GO TO 40
C       MIX THE TWO ARRAYS
   40       K=2
            L=2
            IF(BUF(LST+1).NE.1.E-5)GO TO 170
            IF(BUF(LST+LZ2).NE.1.E-5)GO TO 170
            NXSEC=1
            BUF(LST+LZ1+LZ+1)=1.E-5
            BUF(LST+LZ1+LZ+2)=BUF(LST+2)+BUF(LST+LZ2+1)
C       DETERMINE THE NEXT ENERGY POINT
   50       IF(BUF(LST+1+K).EQ.BUF(LST+LZ2+L))GO TO 90
            IF(BUF(LST+1+K).LT.BUF(LST+LZ2+L))GO TO 70
C       DETERMINE THE CROSS SECTION AT ENERGY POINT BUF(LST+LZ2+L)
            CALL CTERP(BUF(LST+K-1),BUF(LST+K+1),BUF(LST+LZ2+L),
     +                BUF(LST+K), BUF(LST+K+2),SIGMA)
            NXSEC=NXSEC+1
            LP=LZ1+LZ+1+2*(NXSEC-1)
            BUF(LST+LP)=BUF(LST+LZ2+L)
            BUF(LST+LP+1)=BUF(LST+LZ2+L+1)+SIGMA
            L=L+2
            IF(L.LT.LZ)GO TO 50
C       ALL THE POINTS IN THE SECOND ARRAY HAVE NOW BEEN USED
   60       NXSEC=NXSEC+1
            LP=LZ1+LZ+1+2*(NXSEC-1)
            BUF(LST+LP)=BUF(LST+1+K)
            BUF(LST+LP+1)=BUF(LST+2+K)
            K=K+2
            IF(K.LT.LZ1)GO TO 60
            GO TO 100
C       DETERMINE THE CROSS SECTION AT ENERGY POINT BUF(LST+1+K)
   70       CALL CTERP(BUF(LST+LZ2+L-2),BUF(LST+LZ2+L),BUF(LST+1+K),
     +      BUF(LST+LZ2+L-1),BUF(LST+LZ2+L+1),SIGMA)
            NXSEC=NXSEC+1
            LP=LZ1+LZ+1+2*(NXSEC-1)
            BUF(LST+LP)=BUF(LST+1+K)
            BUF(LST+LP+1)=BUF(LST+K+2)+SIGMA
            K=K+2
            IF(K.LT.LZ1)GO TO 50
C       ALL THE POINTS IN THE FIRST ARRAY HAVE NOW BEEN USED
   80       NXSEC=NXSEC+1
            LP=LZ1+LZ+2*NXSEC-1
            BUF(LST+LP)=BUF(LST+LZ2+L)
            BUF(LST+LP+1)=BUF(LST+LZ2+L+1)
            L=L+2
            IF(L.LT.LZ)GO TO 80
            GO TO 100
C       THE ENERGY POINTS COINCIDE
   90       NXSEC=NXSEC+1
            LP=LZ1+LZ+1+2*(NXSEC-1)
            BUF(LST+LP)=BUF(LST+LZ2+L)
            BUF(LST+LP+1)=BUF(LST+2+K)+BUF(LST+LZ2+L+1)
            L=L+2
            K=K+2
            IF((L.LT.LZ).AND.(K.LT.LZ1))GO TO 50
            IF((L.GT.LZ).AND.(K.LT.LZ1))GO TO 60
            IF((L.LT.LZ).AND.(K.GT.LZ1))GO TO 80
C       FINISHED MIXING NOW THIN
  100       L=1
            NXSEC2=1
            LP=LZ1+LZ
            BUF(LST+NXSEC2)=BUF(LST+LP+L)
            BUF(LST+NXSEC2+1)=BUF(LST+LP+L+1)
            KI=0
  110       L=L+2
            KI=KI+1
C       CHECK TO SEE IF AT END OF CROSS SECTION ARRAY
            L2=L+2
            N=2*NXSEC
            IF(L2.LT.N)GO TO 120
C       FINISHED THINING
            NXSEC2=NXSEC2+1
            N=2*(NXSEC2-1)
            BUF(LST+1+N)=BUF(LST+LP+L)
            BUF(LST+2+N)=BUF(LST+LP+L+1)
            LZ=2*NXSEC2
            JI=1
            GO TO 140
  120       DO 130 I=1,KI
C       ESTIMATE THE CROSS SECTION AT KI NODES
               CALL CTERP(BUF(LST+LP+L-2*KI),BUF(LST+LP+L2),
     +                   BUF(LST+LP+L-2*I+2),BUF(LST+LP+L-2*KI+1),
     +                   BUF(LST+LP+L2+1),SIGMA)
               ER=ABS(SIGMA-BUF(LST+LP+L-2*I+3))
C       IF ERROR IS WITHIN ALLOWABLE TOLERANCE, CHECK NEXT POINT
               ERMAX=BUF(LST+LP+L-2*I+3)*TOL
               IF(ER.LE.ERMAX)GO TO 130
C       NOT WITHIN ALLOWABLE TOLERANCE, MUST ADD NODE L-2 TO MESH
               IF(L.GT.3.AND.KI.GT.1) L = L - 2
               NXSEC2=NXSEC2+1
               N=2*(NXSEC2-1)
               BUF(LST+1+N)=BUF(LST+LP+L)
               BUF(LST+2+N)=BUF(LST+LP+L+1)
               KI = 0
               GO TO 110
  130       CONTINUE
C       ALL KI POINTS ARE WITHIN ALLOWABLE TOLERANCE
C       CHECK THE NEXT POINT
            GO TO 110
  140    CONTINUE
C       FINISHED WITH MEDIUM J, NOW STORE IN CORE
         N=2*NXSEC2
         IF(K.EQ.1)N=LZ
         LSIGT(J)=N
         ISIGTS(J)=LAST+1-LMOX3
  150    CONTINUE
         LAST=LAST+N
         LST=LST+N
C       FINISHED MIXING AND THINING
  160 CONTINUE
      GO TO 200
  170 WRITE(IOUT,10000)BUF(LST+1),BUF(LST+LZ2)
10000 FORMAT(' MICAP: ERROR-BEGINNING ENERGY DOES NOT START AT 1.-5',
     +1P2E12.4)
      GOTO 190
  180 CONTINUE
      L=LEN
      WRITE(IOUT,10100)L,N
10100 FORMAT(' MICAP: NOT ENOUGH ROOM TO MIX CROSS SECTIONS',/,5X,
     +'SPACE AVAILABLE=',I10,/,5X,'SPACE NEEDED=',I10)
  190 PRINT '('' CALOR: ERROR in XSECN3 ====> STOP'')'
      STOP
  200 RETURN
      END
