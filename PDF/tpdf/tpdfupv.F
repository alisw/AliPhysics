*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:57  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:18  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      PROGRAM PDFUPV
C
      PARAMETER (LFLPRT = 2)
#include "pdf/w5051p2.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w505110.inc"
#include "pdf/w505120.inc"
#include "pdf/w505122.inc"
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
      CHARACTER*20 PARM(NCHDIM)
#include "pdf/expdp.inc"
     +       VAL(NCHDIM)
C
      DATA SCALE/80.140/
      DATA X00/0.002/, DX0/0.001/
C
      N6 = L6
      IFLPRT = LFLPRT
C
      CALL HLIMIT(NHBMEM)
C
      DSCALE=SCALE
C
C First call to PDFSET to initialize COMMON/W505120/
      PARM(1) = 'Init0'
      VAL(1)  = 0.D0
      CALL PDFSET(PARM,VAL)
C   loop over all existing sets of Nucleon structure functions (SF)
      NPTYPE = 1
      NSF = 0
      DO IGR = 1,NGRMAX
         IF(NPGSMX(NPTYPE,IGR).EQ.0) GOTO 30
         DO ISET = 1,NPGSMX(NPTYPE,IGR)
            WRITE (N6,*) ' =============================================
     +================================================================='
            NSF = NSF + 1
C   book histograms for each set of SF separately
            CALL HBOOK1(ID+NSF,'U Valence$',NB,0.0015,0.9015,0.)
C   force label printing for each set of SF (not only 1st)
            FIRST = .TRUE.
C   select and set parameters
            PARM(1) = 'Nptype'
            VAL(1)  = NPTYPE
            PARM(2) = 'Ngroup'
            VAL(2)  = IGR
            PARM(3) = 'Nset'
            VAL(3)  = ISET
            CALL PDFSET(PARM,VAL)
C   loop over all x bins
            DO 10 I=1,NB
               X= X00 +  (I-1)*DX0
               DX=X
               IF(X.LT.0.0015 .OR. X.GT.0.9015) GOTO 10
               CALL STRUCTM(DX,DSCALE,
     +                    DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
               UPV=DUPV
              IF(X.GT.0.499 .AND. X.LE.0.500) WRITE(N6,1000) X,SCALE,UPV
 1000          FORMAT(/,'  X= ',F6.4,'  Q= ',F6.3,'  UPV= ',F8.4)
               CALL HF1(ID+NSF,X,UPV)
 10         CONTINUE
C   get alpha(s) for selected set of SF at Q = SCALE
            ALF = ALPHAS2(DSCALE)
           WRITE(N6,4000) NPTYPE,NGROUP,NSET,SFNAME(NPTYPE,NGROUP,NSET),
     +                    MODECR(NPTYPE,NGROUP,NSET),NSF
            WRITE(N6,4001) DSCALE,ALF
 4000       FORMAT(/,'  Nptype = ',I1,'  Ngroup = ',I1,
     +               '  Nset = ',I2,'  Name  = "',A8,
     +               '"  CrMode = ',I3,'  HBId = ',I3)
 4001       FORMAT(1H ,' Scale  = ',F8.4,'  alpha(s) = ',F6.4)
C   get error summary for each set of SF
            CALL PDFSTA
   20       CONTINUE
         ENDDO
   30    CONTINUE
      ENDDO
C
      STOP
      END
