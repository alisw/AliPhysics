*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:59  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:18  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      PROGRAM PDFUPVO
C   define loop parameters (= maximal number of available PDF sets)
C+SEQ, W5051P2.
      PARAMETER (MODEMX = 281)
C   define HBOOK settings
      PARAMETER (NHBMEM = 500000)
      COMMON /PAWC/ HMEM(NHBMEM)
      PARAMETER (NB=900, ID=  0)
C   define DOUBLE PRECISION variables for calling sequences
#include "pdf/expdp.inc"
     +       DX,DSCALE,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL,
     +       ALF,ALPHAS2
      REAL    X, SCALE, UPV, DNV, USEA, DSEA, STR, CHM, BOT, TOP, GL
#include "pdf/w50516.inc"
      CHARACTER*20 PARM(20)
#include "pdf/expdp.inc"
     +       VAL(20)
C
      DATA SCALE/80.140/
      DATA X00/0.002/, DX0/0.001/
C
      CALL HLIMIT(NHBMEM)
C
      DSCALE=SCALE
C   loop over all existing sets of proton structure functions (SF)
      DO 20 ISTRUC=0,MODEMX
C   skip obsolete (or non existing) sets of SF
      IF(ISTRUC.EQ. 19                    ) GOTO 20
      IF(ISTRUC.GE. 48 .AND. ISTRUC.LE. 49) GOTO 20
      IF(ISTRUC.GE. 61 .AND. ISTRUC.LE. 69) GOTO 20
      IF(ISTRUC.EQ. 74                    ) GOTO 20
      IF(ISTRUC.GE. 76 .AND. ISTRUC.LE. 79) GOTO 20
      IF(ISTRUC.GE. 81 .AND. ISTRUC.LE.105) GOTO 20
      IF(ISTRUC.GE.108 .AND. ISTRUC.LE.119) GOTO 20
      IF(ISTRUC.GE.123 .AND. ISTRUC.LE.171) GOTO 20
      IF(ISTRUC.GE.174 .AND. ISTRUC.LE.179) GOTO 20
      IF(ISTRUC.GE.183 .AND. ISTRUC.LE.205) GOTO 20
      IF(ISTRUC.GE.208 .AND. ISTRUC.LE.209) GOTO 20
      IF(ISTRUC.GE.214 .AND. ISTRUC.LE.219) GOTO 20
      IF(ISTRUC.GE.223 .AND. ISTRUC.LE.229) GOTO 20
      IF(ISTRUC.GE.233 .AND. ISTRUC.LE.270) GOTO 20
      IF(ISTRUC.GE.274 .AND. ISTRUC.LE.279) GOTO 20
C   book histograms for each set of SF separately
      CALL HBOOK1(ID+ISTRUC,'U Valence$',NB,0.0015,0.9015,0.)
C   force label printing for each set of SF (not only 1st)
      FIRST = .TRUE.
C   select and set parameters
      PARM(1) = 'Mode'
      VAL(1)  = ISTRUC
      CALL PDFSET(PARM,VAL)
C   loop over all x bins
      DO 10 I=1,NB
        X= X00 +  (I-1)*DX0
        DX=X
        IF(X.LT.0.0015 .OR. X.GT.0.9015) GOTO 10
       CALL STRUCTM(DX,DSCALE,
     +                    DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
        UPV=DUPV
        IF(X.GT.0.499 .AND. X.LE.0.500) WRITE(6,1000) X,SCALE,UPV
 1000   FORMAT(/,'  X= ',F6.4,'  Q= ',F6.3,'  UPV= ',F8.4)
        CALL HF1(ID+ISTRUC,X,UPV)
 10   CONTINUE
C   get alpha(s) for selected set of SF at Q = SCALE
      ALF = ALPHAS2(DSCALE)
      WRITE(6,1001) ISTRUC,DSCALE,ALF
 1001 FORMAT(/,'  Mode = ',I3,'   SCALE = ',F8.4,'   alpha(s) = ',F6.4)
C   get error summary for each set of SF
      CALL PDFSTA
 20   CONTINUE
C
      STOP
      END
