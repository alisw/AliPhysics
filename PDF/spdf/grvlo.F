*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:43  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:24  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C. * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C. *                                                                 *
C. *            G R V  -  P A R A M E T R I Z A T I O N S            *
C. *                                                                 *
C. *   adapted by H. Plothow-Besch            August 1991            *
C. *                                                                 *
C. *                                                                 *
C. *                 FOR A DETAILED EXPLANATION SEE :                *
C. *              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/07             *
C. *                                                                 *
C. *   THE PARAMETRIZATIONS ARE FITTED TO OUR PARTON DISTRIBUTIONS   *
C. *   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   *
C. *   / HO) AND  1.E8 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   *
C. *   REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION IS NEG-   *
C. *   LIGIBLE, I.E. BELOW ABOUT 1.E-4, WERE EXCLUDED FROM THE FIT.  *
C. *                                                                 *
C. *              HEAVY QUARK THRESHOLDS  Q(H) = M(H) :              *
C. *         M(C)  =  1.5,  M(B)  =  4.5,  M(T)  =  100  GEV         *
C. *                                                                 *
C. *      CORRESPONDING LAMBDA(F) VALUES FOR F ACTIVE FLAVOURS :     *
C. *      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
C. *             LAMBDA(5)  =  0.153,   LAMBDA(6)  =  0.082  GEV     *
C. *      HO :   LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
C. *             LAMBDA(5)  =  0.131,   LAMBDA(6)  =  0.053  GEV     *
C. *                                                                 *
C. *   HO DISTRIBUTION REFER TO THE MS-BAR SCHEME OF BARDEEN ET AL.  *
C. *                                                                 *
C. * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C...LO PARAMETRIZATIONS :
       SUBROUTINE GRVLO (ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL)
       IMPLICIT REAL (A - Y)
#include "pdf/expdp.inc"
     +       ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL
       REAL  X, Q
       X = ZX
       Q = ZQ
       MU2  = 0.25
       LAM2 = 0.232 * 0.232
       Q2 = Q*Q
       S  = ALOG (ALOG( MAX( MU2,Q2)/LAM2) / ALOG(MU2/LAM2))
       S2 = S * S
       S3 = S2 * S
C...X * (UV + DV) :
       NUD  = 0.663 + 0.191 * S - 0.041 * S2 + 0.031 * S3
       AKUD = 0.326
       AGUD = -1.97 +  6.74 * S -  1.96 * S2
       BUD  =  24.4 -  20.7 * S +  4.08 * S2
       DUD  =  2.86 +  0.70 * S -  0.02 * S2
       UDV  = GRVFV (X, NUD, AKUD, AGUD, BUD, DUD)
C...X * DV :
       ND  = 0.579 + 0.283 * S + 0.047 * S2
       AKD = 0.523 - 0.015 * S
       AGD =  2.22 -  0.59 * S -  0.27 * S2
       BD  =  5.95 -  6.19 * S +  1.55 * S2
       DD  =  3.57 +  0.94 * S -  0.16 * S2
       DV  = GRVFV (X, ND, AKD, AGD, BD, DD)
       ZDV = DV
C...X * UV :
       UV = UDV - DV
       ZUV = UV
C...X * G :
       ALG =  0.558
       BEG =  1.218
       AKG =   1.00 -  0.17 * S
       BKG =   0.0
       AGG =   0.0  + 4.879 * S - 1.383 * S2
       BGG =  25.92 - 28.97 * S + 5.596 * S2
       CG  = -25.69 + 23.68 * S - 1.975 * S2
       DG  =  2.537 + 1.718 * S + 0.353 * S2
       EG  =  0.595 + 2.138 * S
       ESG =  4.066
       GL  = GRVFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
       ZGL = GL
C...X * UBAR = X * DBAR :
       ALU =  1.396
       BEU =  1.331
       AKU =  0.412 - 0.171 * S
       BKU =  0.566 - 0.496 * S
       AGU =  0.363
       BGU = -1.196
       CU  =  1.029 + 1.785 * S - 0.459 * S2
       DU  =  4.696 + 2.109 * S
       EU  =  3.838 + 1.944 * S
       ESU =  2.845
       UDB = GRVFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
       ZUDB = UDB
C...X * SBAR = X * S :
       SS  =   0.0
       ALS =  0.803
       BES =  0.563
       AKS =  2.082 - 0.577 * S
       AGS = -3.055 + 1.024 * S **  0.67
       BS  =   27.4 -  20.0 * S ** 0.154
       DS  =   6.22
       EST =   4.33 + 1.408 * S
       ESS =   8.27 - 0.437 * S
       SB  = GRVFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
       ZSB = SB
C...X * CBAR = X * C :
       SC  =  0.888
       ALC =   1.01
       BEC =   0.37
       AKC =   0.0
       AGC =   0.0
       BC  =   4.24 - 0.804 * S
       DC  =   3.46 + 1.076 * S
       EC  =   4.61 + 1.490 * S
       ESC =  2.555 + 1.961 * S
       CB  = GRVFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
       ZCB = CB
C...X * BBAR = X * B :
       SBO =  1.351
       ALB =   1.00
       BEB =   0.51
       AKB =   0.0
       AGB =   0.0
       BBO =  1.848
       DB  =  2.929 + 1.396 * S
       EB  =   4.71 + 1.514 * S
       ESB =   4.02 + 1.239 * S
       BB  = GRVFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
       ZBB = BB
C...X * TBAR = X * T :
       TB=0.
       ZTB = TB
       RETURN
       END
