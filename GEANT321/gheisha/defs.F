*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:07  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
C
C *** VARIOUS VECTOR OPERATIONS ***
C
C COPIED FROM F14BLO.PAMLIB 23.4.82
C AUTHOR: V.BLOBEL (UNIVERSITY OF HAMBURG)
C DESYLIB
C
C *** BLANK COMMON REPLACED BY /VECUTY/ TO MATCH GEANT/GHEISHA CODE ***
C *** NOTE THAT P(10,100) HAS BECOME PV(10,200) DUE TO THIS         ***
C
C UN-USED ENTRIES REMOVED :
C     "PCOP" "PEXC" "PZER" "PWRT" "DOT4" "IMPU" "IMPULI" "ADD3"
C     "SUB3" "CROSS" "DOT" "SMUL" "NORZ" "PARPER" "PUNIT" "TRAP"
C
C *** ALL ENTRIES RE-WRITTEN AS SUBROUTINES USING ONLY NECESSARY ***
C *** "DOUBLE PRECISION" STMTS. AND ALL SPECIFIC FUNCTIONS HAVE  ***
C *** BEEN CHANGED TO THEIR GENERIC EQUIVALENCES                 ***
C *** NVE 29-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (22-JUNE-1984)
C
      SUBROUTINE DEFS(K,L,M)
C
C --- PV-ARRAY ---
#include "geant321/mxgkgh.inc"
#include "geant321/s_blankp.inc"
C
      DOUBLE PRECISION A,B
C
      MX=M
      MY=M+1
      MZ=M+2
      DO 52 I=1,3
      F=PV(I,K)
      H=PV(I,L)
      PV(I,MY)=F
      PV(I,MZ)=H
 52   CONTINUE
      A=PV(2,MY)*PV(3,MZ)
      B=PV(3,MY)*PV(2,MZ)
      PV(1,MX)=A-B
      A=PV(3,MY)*PV(1,MZ)
      B=PV(1,MY)*PV(3,MZ)
      PV(2,MX)=A-B
      A=PV(1,MY)*PV(2,MZ)
      B=PV(2,MY)*PV(1,MZ)
      PV(3,MX)=A-B
      A=PV(2,MZ)*PV(3,MX)
      B=PV(3,MZ)*PV(2,MX)
      PV(1,MY)=A-B
      A=PV(3,MZ)*PV(1,MX)
      B=PV(1,MZ)*PV(3,MX)
      PV(2,MY)=A-B
      A=PV(1,MZ)*PV(2,MX)
      B=PV(2,MZ)*PV(1,MX)
      PV(3,MY)=A-B
      DO 58 J=MX,MZ
      A=0.0
      DO 54 I=1,3
      A=A+PV(I,J)*PV(I,J)
 54   CONTINUE
      B=SQRT(A)
      IF (B .NE. 0.0) B=1.0/B
      DO 56 I=1,3
      PV(I,J)=B*PV(I,J)
 56   CONTINUE
 58   CONTINUE
      RETURN
      END
