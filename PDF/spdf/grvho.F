*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:42  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
C...HO PARAMATRIZATIONS :
       SUBROUTINE GRVHO (ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL)
       IMPLICIT REAL (A - Y)
#include "pdf/expdp.inc"
     +        ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL
       REAL  X, Q
       X = ZX
       Q = ZQ
       MU2  = 0.3
       LAM2 = 0.248 * 0.248
       Q2 = Q*Q
       S  = ALOG (ALOG( MAX(MU2,Q2)/LAM2) / ALOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
       S3 = S2 * S
C...X * (UV + DV) :
       NUD  = 0.330 + 0.151 * S - 0.059 * S2 + 0.027 * S3
       AKUD = 0.285
       AGUD = -2.28 + 15.73 * S -  4.58 * S2
       BUD  =  56.7 -  53.6 * S + 11.21 * S2
       DUD  =  3.17 +  1.17 * S -  0.47 * S2 +  0.09 * S3
       UDV  = GRVFV (X, NUD, AKUD, AGUD, BUD, DUD)
C...X * DV :
       ND  = 0.459 + 0.315 * DS + 0.515 * S
       AKD = 0.624              - 0.031 * S
       AGD =  8.13 -  6.77 * DS +  0.46 * S
       BD  =  6.59 - 12.83 * DS +  5.65 * S
       DD  =  3.98              +  1.04 * S  -  0.34 * S2
       DV  = GRVFV (X, ND, AKD, AGD, BD, DD)
       ZDV = DV
C...X * UV :
       UV = UDV - DV
       ZUV = UV
C...X * G :
       ALG =  1.128
       BEG =  1.575
       AKG =  0.323 + 1.653 * S
       BKG =  0.811 + 2.044 * S
       AGG =   0.0  + 1.963 * S - 0.519 * S2
       BGG =  0.078 +  6.24 * S
       CG  =  30.77 - 24.19 * S
       DG  =  3.188 + 0.720 * S
       EG  = -0.881 + 2.687 * S
       ESG =  2.466
       GL  = GRVFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
       ZGL = GL
C...X * UBAR = X * DBAR :
       ALU =  0.594
       BEU =  0.614
       AKU =  0.636 - 0.084 * S
       BKU =   0.0
       AGU =  1.121 - 0.193 * S
       BGU =  0.751 - 0.785 * S
       CU  =   8.57 - 1.763 * S
       DU  =  10.22 + 0.668 * S
       EU  =  3.784 + 1.280 * S
       ESU =  1.808 + 0.980 * S
       UDB = GRVFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
       ZUDB = UDB
C...X * SBAR = X * S :
       SS  =   0.0
       ALS =  0.756
       BES =  0.101
       AKS =  2.942 - 1.016 * S
       AGS =  -4.60 + 1.167 * S
       BS  =   9.31 - 1.324 * S
       DS  =  11.49 - 1.198 * S + 0.053 * S2
       EST =  2.630 + 1.729 * S
       ESS =   8.12
       SB  = GRVFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
       ZSB = SB
C...X * CBAR = X * C :
       SC  =  0.820
       ALC =   0.98
       BEC =   0.0
       AKC = -0.625 - 0.523 * S
       AGC =   0.0
       BC  =  1.896 + 1.616 * S
       DC  =   4.12 + 0.683 * S
       EC  =   4.36 + 1.328 * S
       ESC =  0.677 + 0.679 * S
       CB  = GRVFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
       ZCB = CB
C...X * BBAR = X * B :
       SBO =  1.297
       ALB =   0.99
       BEB =   0.0
       AKB =   0.0  - 0.193 * S
       AGB =   0.0
       BBO =   0.0
       DB  =  3.447 + 0.927 * S
       EB  =   4.68 + 1.259 * S
       ESB =  1.892 + 2.199 * S
       BB  = GRVFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
       ZBB = BB
C...X * TBAR = X * T :
       TB=0.
       ZTB = TB
       RETURN
       END
