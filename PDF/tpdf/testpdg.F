*
* $Id$
*
* $Log$
* Revision 1.1.1.3  2000/04/19 13:30:02  mclareni
* Import version 8.04
*
* Revision 1.2  2000/04/09 10:06:56  plothow
* new version 8.00 - MRS 98 and 99, CTEQ5 and GRV98 added
*
* Revision 1.1.1.1  1996/04/12 15:30:17  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      PROGRAM TESTPDG
C
      PARAMETER (LFLPRT = 2)
#include "pdf/w5051p2.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50515.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
#include "pdf/w505110.inc"
#include "pdf/w505120.inc"
#include "pdf/w505122.inc"
C
#include "pdf/expdp.inc"
     +       DX,DQ,DXPDF(-6:6),
     +       ALF,ALPHAS2
      REAL    X, Q, UPV, DNV, USEA, DSEA, SSEA, CSEA, BOT, TOP, GLUE
      DIMENSION XX(5),Q2(5)
      CHARACTER*20 PARM(NCHDIM)
#include "pdf/expdp.inc"
     +       VAL(NCHDIM)
C
      DATA NXX,XX/5, .002,.140,.700,0.01,0.05/
      DATA NQ2,Q2/5, 1.30,24.0,33.5,200.0,200.0/
C
      N6 = L6
      IFLPRT = LFLPRT
C
      CH1=4./9.
      CH2=1./9.
C
      SCAL=SQRT(Q2(2))
C
C First call to PDFSET to initialize COMMON/W505120/
      PARM(1) = 'Init0'
      VAL(1)  = 0.D0
      CALL PDFSET(PARM,VAL)
C.
C Test new PDFLIB versions (4.0 or bigger)
      DO KPTY = 1,NPTYMX
         DO KGRO = 1,NGRMAX
            IF(NPGSMX(KPTY,KGRO) .EQ. 0) GOTO 400
            DO ISET=-9,NPGSMX(KPTY,KGRO)
               IF((KPTY*KGRO).NE.1 .AND. ISET.LE.0) GOTO 300
               WRITE (N6,*) ' ==========================================
     +============================================'
               FIRST=.TRUE.
               IF(ISET.LE.-8) IFLSET = 0
               IF(ISET.LT.0) GOTO 50
               PARM(1) = 'Nptype'
               VAL(1)  = KPTY
               PARM(2) = 'Ngroup'
               VAL(2)  = KGRO
               PARM(3) = 'Nset'
               VAL(3)  = ISET
C... test different NSET writings
               IF(ISET.GT. 5) PARM(3) = 'NSET'
               IF(ISET.GT.10) PARM(3) = 'nset'
               IF(ISET.GT.15) PARM(3) = 'NSet'
               IF(ISET.GT.20) PARM(3) = 'nsET'
C
C... test different NSET writings and INVALID parameter settings
   50          CONTINUE
               IF(ISET.GE.0) THEN
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-9) THEN
                  GOTO 100
               ELSEIF (ISET.EQ.-8) THEN
                  PARM(1) = 'dflm'
                  VAL(1) = 8.D0
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-7) THEN
                  PARM(1) = 'GRV'
                  VAL(1) = 3.D0
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-6) THEN
                  PARM(1) = 'abfkw'
                  VAL(1) = 1.D0
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-5) THEN
                  PARM(1) = 'GRVph'
                  VAL(1) = 2.D0
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-4) THEN
                  PARM(1) = 'Nptype'
                  VAL(1)  = -1
                  PARM(2) = 'Ngroup'
                  VAL(2)  = KGRO
                  PARM(3) = 'Nset'
                  VAL(3)  = ISET
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-3) THEN
                  PARM(1) = 'Nptype'
                  VAL(1)  = -1
                  PARM(2) = 'Ngroup'
                  VAL(2)  = KGRO
                  PARM(3) = 'Nset'
                  VAL(3)  =  0
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-2) THEN
                  PARM(1) = 'Nptype'
                  VAL(1)  = KPTY
                  PARM(2) = 'Ngroup'
                  VAL(2)  = -1
                  PARM(3) = 'Nset'
                  VAL(3)  = ISET
                  CALL PDFSET(PARM,VAL)
               ELSEIF (ISET.EQ.-1) THEN
                  PARM(1) = 'Nptype'
                  VAL(1)  = KPTY
                  PARM(2) = 'Ngroup'
                  VAL(2)  = KGRO
                  PARM(3) = 'Nset'
                  VAL(3)  = 80
                  CALL PDFSET(PARM,VAL)
               ENDIF
C
               DQ = SCAL
               ALF = ALPHAS2(DQ)
               WRITE(N6,4000) NPTYPE,NGROUP,NSET,
     +         SFNAME(NPTYPE,NGROUP,NSET),MODECR(NPTYPE,NGROUP,NSET)
               WRITE(N6,4001) QCDL4,QCDL5,DQ,ALF
 4000          FORMAT(/,'  Nptype = ',I1,'  Ngroup = ',I2,
     +                  '  Nset = ',I3,'  Name  = "',A8,
     +                  '"  CrMode = ',I3)
 4001          FORMAT(1H ,' QCDL4, QCDL5 = ',2F7.4,
     +                    '    Scale = ',F8.4,'  alpha(s) = ',F7.5)
C
  100          CONTINUE
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
               WRITE(N6,3000)X,Q2(IX),FACNP,UP,UV,US,DP,DV,GL,F2P
 3000          FORMAT(1X,' X:',F6.3,2X,'Q2:',F5.1,
     1         3X,'N/P:',F7.4,2X,'UP:',F6.3,
     1         1X,'UV:',F6.3,1X,'US:',F6.3,1X,'DP:',F6.3,1X,'DV:',F6.3,
     2         1X,'GL:',F6.3,2X,'F2P:',F7.4)
C
 200           CONTINUE
C
               IF(ISET.EQ.-9) THEN
                  DQ = SCAL
                  ALF = ALPHAS2(DQ)
                  WRITE(N6,4000) NPTYPE,NGROUP,NSET,
     +            SFNAME(NPTYPE,NGROUP,NSET),MODECR(NPTYPE,NGROUP,NSET)
                  WRITE(N6,4001) QCDL4,QCDL5,DQ,ALF
               ENDIF
               CALL PDFSTA
C
 300           CONTINUE
             ENDDO
 400         CONTINUE
          ENDDO
       ENDDO
C
      STOP
      END
