*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      FUNCTION GBRSGM(Z,T,BCUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates cross-section in current material                  *
C.    *  for discrete(hard) muon BREMSSTRAHLUNG.                       *
C.    *         (SIG in barn/atom)                                     *
C.    *                                                                *
C.    *    ==>Called by : GBRSGA                                       *
C.    *       Author  L.Urban  *********                               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION C(52)
      SAVE C
      PARAMETER (AKSI=2.30,ALFA=1.06,GAM=0.63)
C
      PARAMETER (CONMX=0.130653)
C       CONMX=0.75*SQRT(2.7182...)*EMMU
      DATA C/ 0.949313E-07,-0.819600E-07, 0.529075E-07,-0.832023E-08
     +      , 0.539299E-09,-0.127042E-10,-0.165784E-08,-0.307788E-07
     +      , 0.977905E-08,-0.113658E-08, 0.574481E-10,-0.106221E-11
     +      , 0.968339E-09,-0.108640E-08,-0.177634E-09, 0.889497E-10
     +      ,-0.876878E-11, 0.264303E-12, 0.216263E-08,-0.152680E-08
     +      , 0.380989E-09,-0.455274E-10, 0.264172E-11,-0.596016E-13
     +      , 0.444927E-09,-0.272978E-09, 0.645634E-10,-0.748783E-11
     +      , 0.424890E-12,-0.940837E-14, 0.162289E-10,-0.362486E-11
     +      ,-0.576652E-12, 0.211269E-12,-0.185482E-13, 0.522065E-15
     +      ,-0.215590E-09, 0.112204E-09,-0.819133E-11, 0.145128E-12
     +      ,-0.206029E-09, 0.559940E-10,-0.483350E-11, 0.134252E-12
     +      ,-0.368469E-10, 0.999457E-11,-0.904967E-12, 0.272717E-13
     +      ,-0.303446E-11, 0.853429E-12,-0.785466E-13, 0.236435E-14/
C.
C.    ------------------------------------------------------------------
C.
      GBRSGM=0.
      IF(BCUT.GE.T) GOTO 99
C
      E=T+EMMU
      ECMAX=E-CONMX*Z**0.333333
      IF(BCUT.GE.ECMAX) GOTO 99
      X=LOG(E/EMMU)
      Y=LOG(BCUT/EMMU)
C
      S=0.
      YY=1.
      DO 30 I=1,6
         XX=1.
         DO 20 J=1,6
            K=6*I+J-6
            S=S+C(K)*XX*YY
            XX=XX*X
  20     CONTINUE
         YY=YY*Y
  30  CONTINUE
      SS=0.
      YY=1.
      DO 50 I=1,4
         XX=1.
         DO 40 J=1,4
            K=4*I+J+32
            SS=SS+C(K)*XX*YY
            XX=XX*X
  40     CONTINUE
         YY=YY*Y
  50  CONTINUE
      S=S+Z*SS
      IF(S.LE.0.) GOTO 99
C
      FAC=LOG(ECMAX/BCUT)
      IF(FAC.LE.0.) GOTO 99
      FAC=Z*(Z+AKSI*(1.+GAM*LOG(Z)))*FAC**ALFA
      GBRSGM=FAC*S
C
  99  RETURN
      END
