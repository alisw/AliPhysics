*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE MXMAD(A,B,C,I,J,K)
C
C CERN PROGLIB# F110    MXMAD           .VERSION KERNFOR  1.0   650809
C ORIG. 01/01/64 RKB
C
      DIMENSION A(*),B(*),C(*)
C
      IANDJ=21
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMAD1
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMAD1  (A,B,C,I,J,K)
#endif
C
      IANDJ=22
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMAD2
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMAD2  (A,B,C,I,J,K)
#endif
C
      IANDJ=23
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMAD3
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMAD3  (A,B,C,I,J,K)
#endif
C
      IANDJ=24
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY   (A,B,C,I,J,K)
#endif
C
      IANDJ=11
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY1
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY1  (A,B,C,I,J,K)
#endif
C
      IANDJ=12
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY2
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY2  (A,B,C,I,J,K)
#endif
C
      IANDJ=13
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY3
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMPY3  (A,B,C,I,J,K)
#endif
C
      IANDJ=14
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB   (A,B,C,I,J,K)
#endif
C
      IANDJ=31
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB1
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB1  (A,B,C,I,J,K)
#endif
C
      IANDJ=32
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB2
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB2  (A,B,C,I,J,K)
#endif
C
      IANDJ=33
      GO TO 1
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB3
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMUB3  (A,B,C,I,J,K)
#endif
C
      IANDJ=34
C
    1 N1=IANDJ/10
      N2=IANDJ-10*N1
      IF((I.EQ.0).OR.(K.EQ.0)) GO TO 100
      GO TO  (3,3,5,5),N2
C
    3 IIA=1
      IOA=J
      GO TO 7
    5 IIA=I
      IOA=1
    7 GO TO (9,11,9,11),N2
    9 IIB=K
      IOB=1
      GO TO 13
   11 IIB=1
      IOB=J
   13 IA=1
      IC=1
      DO 25 L=1,I
      IB =1
      DO 23 M=1,K
      GO TO (15,18,16),N1
   15 C(IC)=0.
      GO TO 18
   16 C(IC)=-C(IC)
   18 IF (J.EQ.0) GO TO 23
      JA=IA
      JB=IB
      DO 20 N=1,J
      C(IC)=C(IC)+A(JA)*B(JB)
      JA=JA+IIA
   20 JB=JB+IIB
      IB=IB+IOB
   23 IC=IC+1
   25 IA=IA+IOA
  100 RETURN
      END
