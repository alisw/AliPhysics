*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:40  cernlib
* Version 7.04
*
* Revision 1.2  1996/10/28 15:02:20  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       SUBROUTINE GRVGALO (ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZGL)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*      G R V - P H O T O N - P A R A M E T R I Z A T I O N S      *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE :                *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/31             *
*                                                                 *
*    THE OUTPUT IS ALWAYS   1./ ALPHA(EM) * X * PARTON DENSITY    *
*    output modified by HPB to be always    X * PARTON DENSITY    *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE PARTON DISTRIBUTIONS   *
*   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   *
*   / HO) AND  1.E6 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   *
*                                                                 *
*              HEAVY QUARK THRESHOLDS  Q(H) = M(H) :              *
*         M(C)  =  1.5,  M(B)  =  4.5,  M(T)  =  100  GEV         *
*                                                                 *
*      CORRESPONDING LAMBDA(F) VALUES FOR F ACTIVE FLAVOURS :     *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,   LAMBDA(6)  =  0.082  GEV     *
*      HO :   LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131,   LAMBDA(6)  =  0.053  GEV     *
*                                                                 *
*      HO DISTRIBUTIONS REFER TO THE DIS(GAMMA) SCHEME, SEE :     *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/26             *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
       IMPLICIT REAL (A - Y)
#include "pdf/expdp.inc"
     +        ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZGL
       REAL  X, Q
       DATA ALPHEM/7.29927D-3/
       X = ZX
       Q = ZQ
       MU2  = 0.25
       LAM2 = 0.232 * 0.232
       Q2 = Q*Q
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       SS = SQRT (S)
       S2 = S * S
C...X * U = X * UBAR :
       AL =  1.717
       BE =  0.641
       AK =  0.500 - 0.176 * S
       BK = 15.00  - 5.687 * SS - 0.552 * S2
       AG =  0.235 + 0.046 * SS
       BG =  0.082 - 0.051 * S  + 0.168 * S2
       C  =   0.0  + 0.459 * S
       D  =  0.354 - 0.061 * S
       E  =  4.899 + 1.678 * S
       ES =  2.046 + 1.389 * S
       UL =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZUV = UL * ALPHEM
       ZUB = ZUV
C...X * D = X * DBAR :
       AL =  1.549
       BE =  0.782
       AK =  0.496 + 0.026 * S
       BK =  0.685 - 0.580 * SS + 0.608 * S2
       AG =  0.233 + 0.302 * S
       BG =   0.0  - 0.818 * S  + 0.198 * S2
       C  =  0.114 + 0.154 * S
       D  =  0.405 - 0.195 * S  + 0.046 * S2
       E  =  4.807 + 1.226 * S
       ES =  2.166 + 0.664 * S
       DL  =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZDV = DL * ALPHEM
       ZDB = ZDV
C...X * G :
       AL =  0.676
       BE =  1.089
       AK =  0.462 - 0.524 * SS
       BK =  5.451              - 0.804 * S2
       AG =  0.535 - 0.504 * SS + 0.288 * S2
       BG =  0.364 - 0.520 * S
       C  = -0.323              + 0.115 * S2
       D  =  0.233 + 0.790 * S  - 0.139 * S2
       E  =  0.893 + 1.968 * S
       ES =  3.432 + 0.392 * S
       GL =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZGL = GL * ALPHEM
C...X * S = X * SBAR :
       SF =   0.0
       AL =  1.609
       BE =  0.962
       AK =  0.470              - 0.099 * S2
       BK =  3.246
       AG =  0.121 - 0.068 * SS
       BG = -0.090 + 0.074 * S
       C  =  0.062 + 0.034 * S
       D  =   0.0  + 0.226 * S  - 0.060 * S2
       E  =  4.288 + 1.707 * S
       ES =  2.122 + 0.656 * S
       SL =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZSB = SL * ALPHEM
C...X * C = X * CBAR :
       SF =  0.888
       AL =  0.970
       BE =  0.545
       AK =  1.254 - 0.251 * S
       BK =  3.932              - 0.327 * S2
       AG =  0.658 + 0.202 * S
       BG = -0.699
       C  =  0.965
       D  =   0.0  + 0.141 * S  - 0.027 * S2
       E  =  4.911 + 0.969 * S
       ES =  2.796 + 0.952 * S
       CL =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZCB = CL * ALPHEM
C...X * B = X * BBAR :
       SF =  1.351
       AL =  1.016
       BE =  0.338
       AK =  1.961 - 0.370 * S
       BK =  0.923 + 0.119 * S
       AG =  0.815 + 0.207 * S
       BG = -2.275
       C  =  1.480
       D  = -0.223 + 0.173 * S
       E  =  5.426 + 0.623 * S
       ES =  3.819 + 0.901 * S
       BL =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZBB = BL * ALPHEM
C
       RETURN
       END
