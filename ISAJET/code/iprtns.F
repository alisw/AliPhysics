#include "isajet/pilot.h"
      SUBROUTINE IPRTNS(NPRTNS,PRTNS,IDQ)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-      Fill PINITS common block
C-   Inputs  : 
C-     IDQ(2)= id's of partons starting reaction
C-
C-   Created  10-OCT-1991   Serban D. Protopopescu
C-   Renamed from IPINIT to avoid name clash with Cern Library
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      INTEGER NPRTNS,IDQ(2)
      REAL    PRTNS(4,NPRTNS)
#include "isajet/jetpar.inc"
#include "isajet/pinits.inc"
      REAL AMASS, AM1SQ,AM2SQ,ROOT,QPL,QMN,P1PL,P1MN,P2PL,P2MN
      INTEGER I
C----------------------------------------------------------------------
C          sum P+ and P-, shat
C          assumes sum of transverse momenta is zero
      QPL=0
      QMN=0
      DO 1 I=1,NPRTNS
        QPL=QPL+PRTNS(4,I)+PRTNS(3,I)
        QMN=QMN+PRTNS(4,I)-PRTNS(3,I)
    1 CONTINUE
      SHAT=QPL*QMN
C
C          fill PINITS
      DO 2 I=1,2
        IDINIT(I)=IDQ(I)
        PINITS(5,I)=AMASS(IDQ(I))
        PINITS(1,I)=0.
        PINITS(2,I)=0.
    2 CONTINUE
C          and solve initial kinematics
      AM1SQ=PINITS(5,1)**2
      AM2SQ=PINITS(5,2)**2
      ROOT=SQRT((QPL*QMN-AM1SQ-AM2SQ)**2-4.*AM1SQ*AM2SQ)
      P1PL=(QPL*QMN+AM1SQ-AM2SQ+ROOT)/(2.*QMN)
      P1MN=AM1SQ/P1PL
      P2MN=(QPL*QMN+AM2SQ-AM1SQ+ROOT)/(2.*QPL)
      P2PL=AM2SQ/P2MN
      PINITS(3,1)=.5*(P1PL-P1MN)
      PINITS(4,1)=.5*(P1PL+P1MN)
      PINITS(3,2)=.5*(P2PL-P2MN)
      PINITS(4,2)=.5*(P2PL+P2MN)
  999 RETURN
      END
