#include "isajet/pilot.h"
        REAL FUNCTION SSGZT(E)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + tp + tb
C          Baer's TOPINT
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL E
      DOUBLE PRECISION MG,MT,MS,MZ,SSDLAM,P,PSI,XLOG,PHI,SN,C1,C2
      DOUBLE PRECISION XT,MUT,MUZ,XMIN,XMAX,EMIN,EMAX
C          Convert to double precision
      MG=AMGLSS
      C1=TMP(1)
      C2=TMP(2)
      MS=TMP(3)
      MT=TMP(4)
      MZ=TMP(5)
      SN=TMP(6)
C
      XT=2*E/MG
      MUT=(MT/MG)**2
      MUZ=(MZ/MG)**2
      XMIN=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)-DSQRT((XT**2-4*MUT)*
     $SSDLAM((1.D0+MUT-XT),MUT,MUZ)))/2.D0/(1.D0-XT+MUT)
      XMAX=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)+DSQRT((XT**2-4*MUT)*
     $SSDLAM((1.D0+MUT-XT),MUT,MUZ)))/2.D0/(1.D0-XT+MUT)
      EMIN=XMIN*MG/2.D0
      EMAX=XMAX*MG/2.D0
      P=SQRT(E**2-MT**2)
      PSI=P*E*(MG**2-MZ**2-2*MG*E)*
     $DSQRT(SSDLAM((MG**2+MT**2-2*MG*E),MZ**2,MT**2))/MG/
     $(MG**2+MT**2-2*MG*E)/(MG**2+MT**2-2*MG*E-MS**2)**2
      XLOG=DLOG((MG**2+MT**2-2*MG*EMAX-MS**2)/
     $(MG**2+MT**2-2*MG*EMIN-MS**2))
      PHI=MZ*(-(EMAX-EMIN)-(2*E*MG+MZ**2-MT**2-MS**2)*
     $XLOG/2.D0/MG)/
     $2.D0/MG/(MG**2+MT**2-MS**2-2*MG*E)
      SSGZT=C1*PSI+SN*C2*PHI
      RETURN
      END
