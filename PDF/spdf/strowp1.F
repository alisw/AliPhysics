*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:47  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:05  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C*********************************************************************

      SUBROUTINE STROWP1(X,SCALE,UPV,DNV,SEA,STR,CHM,GL)
C ::::::::::::::  OWENS SET 1 PION STRUCTURE FUNCTION  :::::::::::::::
#include "pdf/impdp.inc"
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION DGAMMA
#endif
#include "pdf/expdp.inc"
     +       COW(3,5,4),TS(6),XQ(9)

C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1),IS=1,5),IP=1,3)/
     1  4.0000D-01,  7.0000D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -6.2120D-02,  6.4780D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     3 -7.1090D-03,  1.3350D-02,  0.0000D+00,  0.0000D+00,  0.0000D+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2),IS=1,5),IP=1,3)/
     1  8.8800D-01,  0.0000D+00,  3.1100D+00,  6.0000D+00,  0.0000D+00,
     2 -1.8020D+00, -1.5760D+00, -1.3170D-01,  2.8010D+00, -1.7280D+01,
     3  1.8120D+00,  1.2000D+00,  5.0680D-01, -1.2160D+01,  2.0490D+01/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3),IS=1,5),IP=1,3)/
     1  9.0000D-01,  0.0000D+00,  5.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -2.4280D-01, -2.1200D-01,  8.6730D-01,  1.2660D+00,  2.3820D+00,
     3  1.3860D-01,  3.6710D-03,  4.7470D-02, -2.2150D+00,  3.4820D-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4),IS=1,5),IP=1,3)/
     1  0.0000D+00, -2.2120D-02,  2.8940D+00,  0.0000D+00,  0.0000D+00,
     2  7.9280D-02, -3.7850D-01,  9.4330D+00,  5.2480D+00,  8.3880D+00,
     3 -6.1340D-02, -1.0880D-01, -1.0852D+01, -7.1870D+00, -1.1610D+01/

       DATA ZEROD/0.D0/, ONED/1.D0/, SIXD/6.D0/
       DATA ALAM/0.2D0/, Q02/4.D0/, QMAX2/2.D3/
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

C...Determine set, Lambda and s expansion variable.
        Q2 = SCALE*SCALE
        Q2IN = MIN( QMAX2,MAX( Q02,Q2))
        SD = LOG( LOG( Q2IN/ALAM**2)/ LOG( Q02/ALAM**2))

C...Calculate structure functions.
        DO 240 KFL=1,4
        DO 230 IS=1,5
  230   TS(IS)=COW(1,IS,KFL)+COW(2,IS,KFL)*SD+
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
  240   CONTINUE

C...Put into output arrays.
        UPV = XQ(1)
        DNV = XQ(1)
        SEA = XQ(3)/SIXD
        STR = XQ(3)/SIXD
        CHM = XQ(4)
        BOT = ZEROD
        TOP = ZEROD
        GL  = XQ(2)
C
        RETURN
        END
