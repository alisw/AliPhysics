#include "isajet/pilot.h"
      REAL FUNCTION SSPOLE(MGMS,MUSQ,AS)
C*********************************************************************
C* Computes the on-shell (pole) gluino mass for given running (MSbar)*
C* gluino mass, defined at scale MUSQ, and given alpha_s (AS). The   *
C* squark masses are stored in the SQUARK COMMON block.              *
C* This function needs the complex functions B0 and B1.              *
C* Contributed by M. Drees; modified by H. Baer                      *
C                                                                    *
C  Version 7.30: Cast COMPLEX*16 to REAL*8 in standard way. :-(      *
C*********************************************************************
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
#include "isajet/ssinf.inc"
      REAL MGMS,MUSQ,AS,MGSQ,FAC
      DOUBLE PRECISIONDMUSQ,DFAC
      COMPLEX*16 SSB1
      DMUSQ=MUSQ
      XLAM = LOG(DMUSQ)
      MGSQ = MGMS*MGMS
C
C          Cast COMPLEX*16 to REAL*8:
C
      DFAC = -8*SSB1(MGSQ,0.,AMULSS) - (SSB1(MGSQ,AMTP,AMT1SS)+
     $SSB1(MGSQ,AMTP,AMT2SS)+SSB1(MGSQ,4.0,AMB1SS)+
     $SSB1(MGSQ,4.0,AMB2SS))
      DFAC = DFAC + 12.D0 + 9.D0*LOG(DMUSQ/MGSQ)
      FAC=DFAC
      SSPOLE = MGMS*(1.0 + .0796*AS*FAC )
      RETURN
      END
