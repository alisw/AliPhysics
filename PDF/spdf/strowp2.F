*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:48  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:05  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C*********************************************************************

      SUBROUTINE STROWP2(X,SCALE,UPV,DNV,SEA,STR,CHM,GL)
C ::::::::::::::  OWENS SET 2 PION STRUCTURE FUNCTION  :::::::::::::::
#include "pdf/impdp.inc"
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION DGAMMA
#endif
#include "pdf/expdp.inc"
     +       COW(3,5,4),TS(6),XQ(9)

C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1),IS=1,5),IP=1,3)/
     1  4.0000D-01,  6.2800D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -5.9090D-02,  6.4360D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     3 -6.5240D-03,  1.4510D-02,  0.0000D+00,  0.0000D+00,  0.0000D+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2),IS=1,5),IP=1,3)/
     1  7.9400D-01,  0.0000D+00,  2.8900D+00,  6.0000D+00,  0.0000D+00,
     2 -9.1440D-01, -1.2370D+00,  5.9660D-01, -3.6710D+00, -8.1910D+00,
     3  5.9660D-01,  6.5820D-01, -2.5500D-01, -2.3040D+00,  7.7580D+00/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3),IS=1,5),IP=1,3)/
     1  9.0000D-01,  0.0000D+00,  5.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -1.4170D-01, -1.6970D-01, -2.4740D+00, -2.5340D+00,  5.6210D-01,
     3 -1.7400D-01, -9.6230D-02,  1.5750D+00,  1.3780D+00, -2.7010D-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4),IS=1,5),IP=1,3)/
     1  0.0000D+00, -8.8200D-02,  1.9240D+00,  0.0000D+00,  0.0000D+00,
     2  6.2290D-02, -2.8920D-01,  2.4240D-01, -4.4630D+00, -8.3670D-01,
     3 -4.0990D-02, -1.0820D-01,  2.0360D+00,  5.2090D+00, -4.8400D-02/

       DATA ZEROD/0.D0/, ONED/1.D0/, SIXD/6.D0/
       DATA ALAM/0.4D0/, Q02/4.D0/, QMAX2/2.D3/
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

C...Determine set, Lambda and s expansion variable.
        Q2 = SCALE*SCALE
        Q2IN = MIN( QMAX2,MAX( Q02,Q2))
        SD = LOG( LOG( Q2IN/ALAM**2)/ LOG( Q02/ALAM**2))

C...Calculate structure functions.
        DO 10 KFL=1,4
        DO 20 IS=1,5
   20   TS(IS)=COW(1,IS,KFL)+COW(2,IS,KFL)*SD+
     &  COW(3,IS,KFL)*SD*SD
        IF(KFL.EQ.1) THEN
#if defined(CERNLIB_SINGLE)
          DENOM = GAMMA(TS(1))*GAMMA(TS(2)+ONED)/GAMMA(TS(1)+TS(2)+ONED)
#endif
#if defined(CERNLIB_DOUBLE)
          DENOM = DGAMMA(TS(1))*DGAMMA(TS(2)+ONED)/
     +                                          DGAMMA(TS(1)+TS(2)+ONED)
#endif
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/DENOM
        ELSE
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2)
        ENDIF
   10   CONTINUE

C...output
        UPV = XQ(1)
        DNV = XQ(1)
        SEA = XQ(3)/SIXD
        STR = XQ(3)/SIXD
        CHM = XQ(4)
        GL  = XQ(2)
C
        RETURN
        END
