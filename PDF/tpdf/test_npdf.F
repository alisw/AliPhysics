#include "pdf/pilot.h"
      PROGRAM TEST_NPDF
C
      PARAMETER (LFLPRT = 2)
#include "pdf/w5051p2.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50511a.inc"
#include "pdf/w50512.inc"
#include "pdf/w50515.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
#include "pdf/w505110.inc"
#include "pdf/w505120.inc"
#include "pdf/w505122.inc"
C
#include "pdf/expdp.inc"
     +       DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL,
     +             DUV,DDV,DUS,DDS,DST,DCH,DB,DT,DG,
     +             RUV,RDV,RUB,RDB,RS,RC,RB,RT,RG,
     +       DQ2,DP2,DANO,ALF,ALPHAS2
      REAL   X, Q, UPV, DNV, USEA, DSEA, SSEA, CSEA, BOT, TOP, GLUE
      REAL   P2, ANO
      DIMENSION XX(6),Q2(6),ANO(3)
      CHARACTER*20 PARM(NCHDIM)
#include "pdf/expdp.inc"
     +       VAL(NCHDIM)
C
      DATA NXX,XX/6, .0001,.001,.01,.1,.3,.700/
      DATA NQ2,Q2/6, 2.25,10000.,39.9,2.25,39.9,39.9/
      DATA P2/5./
      DATA IP2/0/
      DATA ANO/12.,117.,208./
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
      KPTY=1
         DO KGRO = 1,NGRMAX
            IF(NPGSMX(KPTY,KGRO) .EQ. 0) GOTO 400
            DO ISET=0,NPGSMX(KPTY,KGRO)
               IF((KPTY*KGRO).NE.1 .AND. ISET.LE.0) GOTO 300
               WRITE (N6,*) ' ==========================================
     +============================================'
               FIRST=.TRUE.
               PARM(1) = 'Nptype'
               VAL(1)  = KPTY
               PARM(2) = 'Ngroup'
               VAL(2)  = KGRO
               PARM(3) = 'Nset'
               VAL(3)  = ISET
C now do the NPDF's
                  PARM(4) = 'NAtype'
                  VAL(4)  = 4
                  PARM(5) = 'NAgroup'
                  VAL(5)  = 1
                  PARM(6) = 'NAset'
                  VAL(6)  = 1
C... test different NSET writings
               IF(ISET.GT. 5) PARM(3) = 'NSET'
               IF(ISET.GT.10) PARM(3) = 'nset'
               IF(ISET.GT.15) PARM(3) = 'NSet'
               IF(ISET.GT.20) PARM(3) = 'nsET'
C
                  CALL PDFSET(PARM,VAL)
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
               WRITE(N6,4002) NATYPE,NAGROUP,NASET,
     +                        SFNAME(NATYPE,NAGROUP,NASET)
 4002          FORMAT(/,'  NAtype = ',I1,'  NAgroup = ',I2,
     +                  '  NAset = ',I3,'  Name  = "',A8,'"')
C
  100          CONTINUE
               DO IA=1,3
               IAFL=0
               DO 200 IX=1,NXX
               X=XX(IX)
               Q=SQRT(Q2(IX))
                 DX=X
                 DQ=Q
                 DANO=ANO(IA)
                    CALL STRUCTA(DX,DQ,DANO,
     +                    DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
                 UPV=DUPV
                 DNV=DDNV
                 USEA=DUSEA
                 DSEA=DDSEA
                 SSEA=DSTR
                 CSEA=DCHM
                 GLUE=DGL
C
               FFPU = UPV + USEA
               FFPUB = USEA
               FFPD = DNV + DSEA
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
               IF (IAFL.EQ.0) WRITE(N6,2999) ANO(IA)
 2999          FORMAT(1X,' A:',F5.0)
               IAFL=1
                IF(Q2(IX).LT.100.) THEN
               WRITE(N6,3000) X,Q2(IX),FACNP,UP,UV,US,DP,DV,GL,F2P
                ELSE
               WRITE(N6,3001) X,Q2(IX),FACNP,UP,UV,US,DP,DV,GL,F2P
                ENDIF
 3000          FORMAT(1X,' X:',F6.3,2X,'Q2:',F6.1,
     1         3X,'N/P:',F7.4,2X,'UP:',F6.3,
     1         1X,'UV:',F6.3,1X,'US:',F6.3,1X,'DP:',F6.3,1X,'DV:',F6.3,
     2         1X,'GL:',F6.3,2X,'F2P:',F7.4)
 3001           FORMAT(1X,' X:',F6.3,2X,'Q2:',F6.0,
     1         3X,'N/P:',F7.4,2X,'UP:',F6.3,
     1         1X,'UV:',F6.3,1X,'US:',F6.3,1X,'DP:',F6.3,1X,'DV:',F6.3,
     2         1X,'GL:',F6.3,2X,'F2P:',F7.4)
C
 200           CONTINUE
               WRITE(N6,*) ' '
               ENDDO
C
               CALL PDFSTA
C
 300           CONTINUE
             ENDDO
 400         CONTINUE
          ENDDO
C       ENDDO
C
      STOP
      END
