#include "isajet/pilot.h"
        REAL FUNCTION SSGZG2(XARG)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + gl
C          Baer's FUNI1- removed masses to simplify
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL XARG
        DOUBLE PRECISION TOP,BOT,XX,F,MQ,MX,MS,MG
C
        MQ=TMP(1)
        MX=TMP(2)
        MS=TMP(3)
        MG=AMGLSS
C
        XX=XARG
        TOP=-MG**2*XX*(1.D0-XX)+MS**2*XX+MQ**2*(1.D0-XX)
        BOT=-MX**2*XX*(1.D0-XX)+MS**2*XX+MQ**2*(1.D0-XX)
        F=DLOG(TOP/BOT)
        SSGZG2=F
        RETURN
        END
