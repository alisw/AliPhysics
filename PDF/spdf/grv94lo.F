*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:36  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*    G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S     *
*                                                                 *
*                         1994 UPDATE                             *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE                  *
*                   M. GLUECK, E.REYA, A.VOGT :                   *
*                   DO-TH 94/24  =  DESY 94-206                   *
*                    (TO APPEAR IN Z. PHYS. C)                    *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE EVOLVED PARTONS FOR    *
*        Q**2 / GEV**2  BETWEEN   0.4   AND  1.E6                 *
*             X         BETWEEN  1.E-5  AND   1.                  *
*   LARGE-X REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION   *
*   IS NEGLIGIBLY SMALL, WERE EXCLUDED FROM THE FIT.              *
*                                                                 *
*   HEAVY QUARK THRESHOLDS  Q(H) = M(H)  IN THE BETA FUNCTION :   *
*                   M(C)  =  1.5,  M(B)  =  4.5                   *
*   CORRESPONDING LAMBDA(F) VALUES IN GEV FOR  Q**2 > M(H)**2 :   *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,                                *
*      NLO :  LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131.                                *
*   THE NUMBER OF ACTIVE QUARK FLAVOURS IS  NF = 3  EVERYWHERE    *
*   EXCEPT IN THE BETA FUNCTION, I.E. THE HEAVY QUARKS C,B,...    *
*   ARE NOT PRESENT AS PARTONS IN THE Q2-EVOLUTION.               *
*   IF NEEDED, HEAVY QUARK DENSITIES CAN BE TAKEN FROM THE 1991   *
*   GRV PARAMETRIZATION.                                          *
*                                                                 *
*   NLO DISTRIBUTIONS ARE GIVEN IN MS-BAR FACTORIZATION SCHEME    *
*   (SUBROUTINE GRV94HO) AS WELL AS IN THE DIS SCHEME (GRV94DI),  *
*   THE LEADING ORDER PARAMETRIZATION IS PROVIDED BY "GRV94LO".   *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS :
*
*    X   = MOMENTUM FRACTION
*    Q2  = SCALE Q**2 IN GEV**2
*
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION) :
*
*    UV  = U(VAL) = U - U(BAR)
*    DV  = D(VAL) = D - D(BAR)
*    DEL = D(BAR) - U(BAR)
*    UDB = U(BAR) + D(BAR)
*    SB  = S = S(BAR)
*    GL  = GLUON
*
*...LO PARAMETRIZATION :
*
       SUBROUTINE GRV94LO (ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZTB,ZGL)
       IMPLICIT REAL (A - Y)
C+SEQ, IMPDP.
#include "pdf/expdp.inc"
     +      ZX,ZQ,ZUV,ZDV,ZDEL,ZUDB,ZUB,ZDB,ZSB,ZCB,ZBB,ZTB,ZGL
       REAL X, Q
       X = ZX
       Q = ZQ
       MU2  = 0.23
       LAM2 = 0.2322   * 0.2322
       Q2 = Q*Q
       ARG1 = MAX(MU2,Q2)/LAM2
       ARG2 = MU2/LAM2
       S  =  LOG ( LOG( ARG1  ) /  LOG( ARG2   ))
       DS =  SQRT (S)
       S2 = S * S
       S3 = S2 * S
*...UV :
       NU  =  2.284 + 0.802 * S + 0.055 * S2
       AKU =  0.590 - 0.024 * S
       BKU =  0.131 + 0.063 * S
       AU  = -0.449 - 0.138 * S - 0.076 * S2
       BU  =  0.213 + 2.669 * S - 0.728 * S2
       CU  =  8.854 - 9.135 * S + 1.979 * S2
       DU  =  2.997 + 0.753 * S - 0.076 * S2
       UV  = GR94FV (X, NU, AKU, BKU, AU, BU, CU, DU)
       ZUV = UV
*...DV :
       ND  =  0.371 + 0.083 * S + 0.039 * S2
       AKD =  0.376
       BKD =  0.486 + 0.062 * S
       AD  = -0.509 + 3.310 * S - 1.248 * S2
       BD  =  12.41 - 10.52 * S + 2.267 * S2
       CD  =  6.373 - 6.208 * S + 1.418 * S2
       DD  =  3.691 + 0.799 * S - 0.071 * S2
       DV  = GR94FV (X, ND, AKD, BKD, AD, BD, CD, DD)
       ZDV = DV
*...DEL :
       NE  =  0.082 + 0.014 * S + 0.008 * S2
       AKE =  0.409 - 0.005 * S
       BKE =  0.799 + 0.071 * S
       AE  = -38.07 + 36.13 * S - 0.656 * S2
       BE  =  90.31 - 74.15 * S + 7.645 * S2
       CE  =  0.0
       DE  =  7.486 + 1.217 * S - 0.159 * S2
       DEL = GR94FV (X, NE, AKE, BKE, AE, BE, CE, DE)
       ZDEL= DEL
*...UDB :
       ALX =  1.451
       BEX =  0.271
       AKX =  0.410 - 0.232 * S
       BKX =  0.534 - 0.457 * S
       AGX =  0.890 - 0.140 * S
       BGX = -0.981
       CX  =  0.320 + 0.683 * S
       DX  =  4.752 + 1.164 * S + 0.286 * S2
       EX  =  4.119 + 1.713 * S
       ESX =  0.682 + 2.978 * S
       UDB = GR94FW(X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
       ZUDB= UDB
*...UB
       ZUB = (ZUDB - ZDEL)/2.D0
*...DB
       ZDB = (ZUDB + ZDEL)/2.D0
*...SB :
       ALS =  0.914
       BES =  0.577
       AKS =  1.798 - 0.596 * S
       AS  = -5.548 + 3.669 * DS - 0.616 * S
       BS  =  18.92 - 16.73 * DS + 5.168 * S
       DST =  6.379 - 0.350 * S  + 0.142 * S2
       EST =  3.981 + 1.638 * S
       ESS =  6.402
       SB  = GR94FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
       ZSB = SB
*...GL :
       ALG =  0.524
       BEG =  1.088
       AKG =  1.742 - 0.930 * S
       BKG =        - 0.399 * S2
       AG  =  7.486 - 2.185 * S
       BG  =  16.69 - 22.74 * S  + 5.779 * S2
       CG  = -25.59 + 29.71 * S  - 7.296 * S2
       DG  =  2.792 + 2.215 * S  + 0.422 * S2 - 0.104 * S3
       EG  =  0.807 + 2.005 * S
       ESG =  3.841 + 0.316 * S
       GL  = GR94FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       ZGL = GL
*...
       ZCB = 0.D0
       ZBB = 0.D0
       ZTB = 0.D0
*...
       RETURN
       END
