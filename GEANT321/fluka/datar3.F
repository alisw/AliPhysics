*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE DATAR3.FOR
*COPY DATAR3
*
*=== datar3 ===========================================================*
*
      SUBROUTINE DATAR3
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/bamjcm.inc"
#include "geant321/inpdat.inc"
      LOGICAL LSWTCH
      COMMON /FKVALI/ LSWTCH
      DIMENSION IVV(36),IP(36),IB(126),IBB(126),IA(126),IAA(126)
C     DEFINE THE FIELDS FOR PARTICLE CLASSIFICATION
C     IMPS=PSEUDO SCALAR MESONS (SPIN=0)
C     IMVE=VECTOR MESONS (SPIN=1)
C     IB08(IA08)=BARYONS (ANTIBARYONS) (SPIN=1/2)
C     IB10(IA10)=BARYONS (ANTIBARYONS) (SPIN=3/2)
      DATA IP/
     *23,14,16,116,0,0,13,23,25,117,0,0,15,24,31,120,0,0,119,118,121,
     *122,14*0/
      LSWTCH = .FALSE.
      L=0
C     PRINT 7769
C7769  FORMAT (' ***********************************************'/
C    1' *********************************************************'/
C    1'         ETA SUPRESSED IN DATAR3 AND HKLASS   SEPT 85 J.RANFT'/
C    1'     ETA (31)--->RHO0(33);  ETA*(95)--->OMEG(35)          '/
C    1'*************************************************************'/
C    1'*************************************************************')
      DO 1 I=1,6
      DO 2 J=1,6
      L=L+1
      IMPS(I,J)=IP(L)
    2 CONTINUE
    1 CONTINUE
      DATA IVV/
     *33,34,38,123,0,0,32,33,39,124,0,0,36,37,96,127,0,0,126,125,128,
     *129,14*0/
      L=0
      DO 3 I=1,6
      DO 4 J=1,6
      L=L+1
      IMVE(I,J)=IVV(L)
    4 CONTINUE
    3 CONTINUE
      DATA IB/
     *0,1,21,140,0,0,8,22,137,0,0,97,138,0,0,146,5*0,
     *1,8,22,137,0,0,0,20,142,0,0,98,139,0,0,147,5*0,
     *21,22,97,138,0,0,20,98,139,0,0,0,145,0,0,148,5*0,
     *140,137,138,146,0,0,142,139,147,0,0,145,148,50*0/
      L=0
      DO 5 I=1,6
      DO 6 J=1,21
      L=L+1
      IB08(I,J)=IB(L)
    6 CONTINUE
    5 CONTINUE
      DATA IBB/
     *53,54,104,161,0,0,55,105,162,0,0,107,164,0,0,167,5*0,
     *54,55,105,162,0,0,56,106,163,0,0,108,165,0,0,168,5*0,
     *104,105,107,164,0,0,106,108,165,0,0,109,166,0,0,169,5*0,
     *161,162,164,167,0,0,163,165,168,0,0,166,169,0,0,170,47*0/
      L=0
      DO 7 I=1,6
      DO 8 J=1,21
      L=L+1
      IB10(I,J)=IBB(L)
    8 CONTINUE
    7 CONTINUE
      DATA IA/
     *0,2,99,152,0,0,9,100,149,0,0,102,150,0,0,158,5*0,
     *2,9,100,149,0,0,0,101,154,0,0,103,151,0,0,159,5*0,
     *99,100,102,150,0,0,101,103,151,0,0,0,157,0,0,160,5*0,
     *152,149,150,158,0,0,154,151,159,0,0,157,160,50*0/
      L=0
      DO 9 I=1,6
      DO 10 J=1,21
      L=L+1
      IA08(I,J)=IA(L)
   10 CONTINUE
    9 CONTINUE
      DATA IAA/
     *67,68,110,171,0,0,69,111,172,0,0,113,174,0,0,177,5*0,
     *68,69,111,172,0,0,70,112,173,0,0,114,175,0,0,178,5*0,
     *110,111,113,174,0,0,112,114,175,0,0,115,176,0,0,179,5*0,
     *171,172,174,177,0,0,173,175,178,0,0,176,179,0,0,180,47*0/
      L=0
      DO 11 I=1,6
      DO 12 J=1,21
      L=L+1
      IA10(I,J)=IAA(L)
   12 CONTINUE
   11 CONTINUE
C     DEFINE THE FREE PARAMETERS FOR THE MONTE-CARLO PROGRAMMES BAMJET
C     PARJET HAPAQ
      A1=0.88D0
      B3=6.D0
      B1=8.D0
      B2=8.D0
      ISU=4
      BET=8.D0
      AS=0.25D0
      AME=0.93D0
      LT=0
      LE=0
      B8=0.33D0
      DIQ=0.375D0
C
C     BAMJCM INITIALIZATION
      DO 13 J = 1,KMXJCM
         IV (J) = 0
         RE (J) = 0.D+00
         KFR1 (J) = 0
         KFR2 (J) = 0
         RPX (J) = 0.D+00
         RPY (J) = 0.D+00
   13 CONTINUE
      RETURN
      END
