*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:16  cernlib
* Version 7.04
*
* Revision 1.2  1996/10/28 15:13:56  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:48  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE STRBEBC(DX,DQ,DUPV,DDNV,DSEA,DSTR,DGL)
C     SUBROUTINE STRBEBC(DX,DQ,DUDV,DUDS,DSEA,DSTR,DGL)
*********************************************************************
*                                                                   *
*     P.Gorbunov   /ITEP                                            *
*     A.Maslennikov/ITEP                             May   1990     *
*     coded after Z. Phys. C 36 (1987) 1, BEBC WA59                 *
*                                                                   *
*     Original publication:                                         *
*     A.J. Buras and K.J.F. Gaemers, Nucl. Phys. B132 (1978) 249    *
*                                                                   *
*                                                                   *
*     adapted and corrected by H.Plothow-Besch/PPE   April 1991     *
*                                                                   *
*********************************************************************
C.
#include "pdf/expdp.inc"
     +       DX,DQ,DUPV,DDNV,DSEA,DSTR,DGL
      REAL    X, Q, UPV, DNV, SEA, STR, GL, UDVAL, UDSEA
      PARAMETER (BVQ02=0.77, BSQ02=5.3, BGQ02=4., Q02=0.5)
      PARAMETER (EPS=0.5)
C
      EXTERNAL VALDIS,V3
C
      COMMON/W5051I1/
     +           GPSI(2:3),GAA(2:3),GPSIA(2:3),GAPSI(2:3),GMINUS(2:3),
     1           AL(2:3),BET(2:3),GPLUS(2:3),
     2           VNQ02(2:3),SNQ02(2:3),GNQ02(2:3)
C
      PARAMETER (IX=112,IQ=40,NARG=2,KFUN=2)
      REAL ARG(NARG),ENT(IX+IQ),XQDUM(IX,IQ,KFUN),Q2P(IQ),XCO(IX)
      INTEGER NENT(NARG)
C
      COMMON/W5051I2/ Q2L
      DATA ZERO/0./,HALF/0.5/,ONE/1./,TWO/2./,THRE/3./
      DATA INIT/0/
C
      DATA XCO   /
     +0.5000E-04,0.5545E-04,0.6150E-04,0.6820E-04,0.7564E-04,0.8389E-04,
     +0.9303E-04,0.1032E-03,0.1144E-03,0.1269E-03,0.1407E-03,0.1561E-03,
     +0.1731E-03,0.1920E-03,0.2129E-03,0.2361E-03,0.2619E-03,0.2904E-03,
     +0.3221E-03,0.3572E-03,0.3961E-03,0.4393E-03,0.4872E-03,0.5404E-03,
     +0.5993E-03,0.6646E-03,0.7371E-03,0.8174E-03,0.9066E-03,0.1005E-02,
     +0.1115E-02,0.1237E-02,0.1371E-02,0.1521E-02,0.1687E-02,0.1871E-02,
     +0.2075E-02,0.2301E-02,0.2552E-02,0.2830E-02,0.3139E-02,0.3481E-02,
     +0.3860E-02,0.4281E-02,0.4748E-02,0.5266E-02,0.5840E-02,0.6476E-02,
     +0.7182E-02,0.7966E-02,0.8834E-02,0.9797E-02,0.1087E-01,0.1205E-01,
     +0.1336E-01,0.1482E-01,0.1644E-01,0.1823E-01,0.2022E-01,0.2242E-01,
     +0.2487E-01,0.2758E-01,0.3058E-01,0.3392E-01,0.3762E-01,0.4172E-01,
     +0.4627E-01,0.5131E-01,0.5690E-01,0.6311E-01,0.6999E-01,0.7762E-01,
     +0.8608E-01,0.9547E-01,0.1059E+00,0.1174E+00,0.1302E+00,0.1444E+00,
     +0.1602E+00,0.1776E+00,0.1970E+00,0.2185E+00,0.2423E+00,0.2687E+00,
     +0.2980E+00,0.3305E+00,0.3666E+00,0.4065E+00,0.4508E+00,0.5000E+00,
     +0.5321E+00,0.5643E+00,0.5964E+00,0.6286E+00,0.6607E+00,0.6929E+00,
     +0.7250E+00,0.7571E+00,0.7893E+00,0.8214E+00,0.8536E+00,0.8857E+00,
     +0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/
C
      SAVE NENT, ENT, XQDUM
C
      IF(INIT.EQ.0) THEN
C--------- initialize Buras-Gaemers parametrization
C          BSMALL=(33.-8.)/3.
C... in Buras-Gaemers paper (Nucl.Phys.132 (1978) 249) this is like:
           BSMALL=(33.-8.)/4.
           S=0.
           DO N=2,3
              S=S+1./N
              GPSI  (N)=8./6./BSMALL* (1.-2./N/(N+1)+4.*S)
              GAA   (N)=3./BSMALL* (1./3.-4./N/(N-1)-4./(N+1)/(N+2) +
     1                  4.*S)  + 8./3./BSMALL
              GPSIA (N)=2./BSMALL* (8./(N+2)+16./N/(N+1)/(N+2))
              GAPSI (N)=4./3./BSMALL* (1./(N+1)+2./N/(N-1))
              GMINUS(N)=0.5*(GPSI(N)+GAA(N) -
     1                  SQRT((GPSI(N)-GAA(N))**2+4.*GPSIA(N)*GAPSI(N)))
              GPLUS (N)=0.5*(GPSI(N)+GAA(N) +
     1                  SQRT((GPSI(N)-GAA(N))**2+4.*GPSIA(N)*GAPSI(N)))
              AL    (N)=GPSIA(N)*GAPSI(N)/
     1              (GPSIA(N)*GAPSI(N)+(GMINUS(N)-GPSI(N))**2)
              BET   (N)=0.5*AL(N)*(GPSI(N)-GMINUS(N))/GAPSI(N)
           ENDDO

           Q2L     =Q02
C
           VNQ02(2)= GAUSS(VALDIS,ZERO,ONE,EPS)
           VNQ02(3)= GAUSS(V3    ,ZERO,ONE,EPS)
           SNQ02(2)=VNQ02(2)*(ONE-BVQ02)/BVQ02
           SNQ02(3)=SNQ02(2)/(BSQ02+TWO)
           GNQ02(2)=ONE-SNQ02(2)-VNQ02(2)
           GNQ02(3)=GNQ02(2)/(BGQ02+TWO)

*--------- initialize interpolation tables
           DO I=1,IQ
              Q2L=I-HALF
              ENT(IX+I)=LOG10(Q2L)
              DO J=1,IX
                 ENT(J)=LOG10(XCO(J))
                 XX = XCO(J)
                 XQDUM(J,I,1)=VALDIS(XX)
                 XQDUM(J,I,2)=SEADIS(XX)
              ENDDO
           ENDDO
           NENT(1)=IX
           NENT(2)=IQ

           INIT = 1
      ENDIF
      X = DX
      Q = DQ
      Q2IN  = Q*Q
      ARG(1)= LOG10(X)
      ARG(2)= LOG10( MAX( Q02,Q2IN))
      UDVAL = MAX(ZERO,FINT(NARG,ARG,NENT,ENT,XQDUM(1,1,1)))
      UPV = UDVAL*TWO/THRE
      DNV = UPV/TWO
      UDSEA = MAX(ZERO,FINT(NARG,ARG,NENT,ENT,XQDUM(1,1,2))) /
     +        (TWO+HALF)
      SEA  = UDSEA*HALF
      STR  = SEA*HALF
      GL  = 0.
C
      DUPV = UPV
      DDNV = DNV
      DSEA = SEA
      DSTR = STR
      DGL  = GL
      RETURN
      END
