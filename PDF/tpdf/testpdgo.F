*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:56  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:18  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      PROGRAM TESTPDGO
C
      PARAMETER (LFLPRT = 2)
C+SEQ, W5051P2.
      PARAMETER (MODEMX = 281)
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50515.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
C
#include "pdf/expdp.inc"
     +       DX,DQ,DXPDF(-6:6),
     +       ALF,ALPHAS2
      REAL    X, Q, UPV, DNV, USEA, DSEA, SSEA, CSEA, BOT, TOP, GLUE
      DIMENSION XX(3),Q2(3)
      CHARACTER*20 PARM(20)
#include "pdf/expdp.inc"
     +       VAL(20)
C
      DATA NXX,XX/3, .002,.140,.700/
      DATA NQ2,Q2/3, 1.30,24.0,33.5/
C
      N6 = L6
      IFLPRT = LFLPRT
C
      CH1=4./9.
      CH2=1./9.
C
      SCAL=SQRT(Q2(2))
C
      DO 300 ISTRUC=-3,MODEMX
      FIRST=.TRUE.
      IF(ISTRUC.LE.-2) IFLSET = 0
C... test different MODE writings
      IF(ISTRUC .LT. 0) GOTO 310
      PARM(1) = 'MODE'
      IF(ISTRUC.GT.20) PARM(1) = 'Mode'
      IF(ISTRUC.GT.30) PARM(1) = 'mode'
      IF(ISTRUC.GT.40) PARM(1) = 'MOde'
      IF(ISTRUC.GT.50) PARM(1) = 'moDE'
      VAL(1)  = ISTRUC
      IF(ISTRUC.EQ. 19                    ) GOTO 300
      IF(ISTRUC.GE. 48 .AND. ISTRUC.LE. 49) GOTO 300
      IF(ISTRUC.GE. 61 .AND. ISTRUC.LE. 69) GOTO 300
      IF(ISTRUC.EQ. 74                    ) GOTO 300
      IF(ISTRUC.GE. 76 .AND. ISTRUC.LE. 79) GOTO 300
      IF(ISTRUC.GE. 81 .AND. ISTRUC.LE.105) GOTO 300
      IF(ISTRUC.GE.108 .AND. ISTRUC.LE.129) GOTO 300
      IF(ISTRUC.GE.133 .AND. ISTRUC.LE.171) GOTO 300
      IF(ISTRUC.GE.174 .AND. ISTRUC.LE.179) GOTO 300
      IF(ISTRUC.GE.183 .AND. ISTRUC.LE.205) GOTO 300
      IF(ISTRUC.GE.208 .AND. ISTRUC.LE.209) GOTO 300
      IF(ISTRUC.GE.214 .AND. ISTRUC.LE.219) GOTO 300
      IF(ISTRUC.GE.223 .AND. ISTRUC.LE.229) GOTO 300
      IF(ISTRUC.GE.233 .AND. ISTRUC.LE.270) GOTO 300
      IF(ISTRUC.GE.274 .AND. ISTRUC.LE.279) GOTO 300
      IF(ISTRUC.GE.282                    ) GOTO 300
C
C... test different MODE writings and INTEGER settings
  310 CONTINUE
      IF(ISTRUC.GE.0 .AND. ISTRUC.LE.MODEMX) THEN
         CALL PDFSET(PARM,VAL)
      ELSEIF (ISTRUC.EQ.-3) THEN
         GOTO 5
      ELSEIF (ISTRUC.EQ.-2) THEN
         PARM(1) = 'MODE'
         VAL(1)  = 45.D0
         CALL PDFSET(PARM,VAL)
      ELSEIF (ISTRUC.EQ.-1) THEN
         PARM(1) = 'mode'
         VAL(1)  = 45.D0
         CALL PDFSET(PARM,VAL)
      ENDIF
C
      DQ = SCAL
      ALF = ALPHAS2(DQ)
      WRITE (N6,2000) MODE,QCDL4,QCDL5,DQ,ALF
 2000 FORMAT(/,'  Mode = ',I3,'  QCDL4, QCDL5 = ',2F7.4,
     +         '  SCALE = ',F7.4,'  alpha(s) = ',F8.5)
C
   5  CONTINUE
      DO 200 IX=1,NXX
      X=XX(IX)
      Q=SQRT(Q2(IX))
        DX=X
        DQ=Q
        CALL PFTOPDG(DX,DQ,DXPDF)
        DNV=DXPDF(1)
        UPV=DXPDF(2)
        DSEA=DXPDF(-1)
        USEA=DXPDF(-2)
        SSEA=DXPDF(3)
        CSEA=DXPDF(4)
        GLUE=DXPDF(0)
C
      FFPU = UPV
      FFPUB = USEA
      FFPD = DNV
      FFPDB = DSEA
      FFPS = SSEA
      FFPSB = SSEA
      FFPC = CSEA
      FFPCB = CSEA
      FFNU = FFPD
      FFNUB = FFPDB
      FFND = FFPU
      FFNDB = FFPUB
      FFNS = FFPS
      FFNSB = FFPSB
      FFNC = FFPC
      FFNCB = FFPCB
      FFPGL = GLUE
      SFFP=CH1*(FFPU+FFPUB+FFPC+FFPCB)+CH2*(FFPD+FFPDB+FFPS+FFPSB)
      SFFN=CH1*(FFNU+FFNUB+FFNC+FFNCB)+CH2*(FFND+FFNDB+FFNS+FFNSB)
      FACNP = 0.
      IF (SFFP .NE. 0.) FACNP=SFFN/SFFP
      UP = FFPU
      UV = FFPU - FFPUB
      US = FFPUB
      DP = FFPD
      DV = FFPD - FFPDB
      DS = FFPDB
      GL = FFPGL
      F2P = SFFP
C
      WRITE(N6,1000)X,Q2(IX),FACNP,UP,UV,US,DP,DV,GL,F2P
 1000 FORMAT(1X,' X:',F6.3,2X,'Q2:',F5.1,
     1    3X,'N/P:',F7.4,2X,'UP:',F6.3,
     1    1X,'UV:',F6.3,1X,'US:',F6.3,1X,'DP:',F6.3,1X,'DV:',F6.3,
     2    1X,'GL:',F6.3,2X,'F2P:',F7.4)
C
 200  CONTINUE
C
      IF(ISTRUC.EQ.-3) THEN
        DQ = SCAL
        ALF = ALPHAS2(DQ)
        WRITE (N6,2000) MODE,QCDL4,QCDL5,DQ,ALF
      ENDIF
      CALL PDFSTA
C
 300  CONTINUE
C
      STOP
      END
