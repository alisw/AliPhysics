*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.39  by  S.Giani
*-- Author :
      SUBROUTINE PHASP
C
C *** NVE 29-MAR-1988 CERN GENEVA ***
C
C CALLED BY : NUCREC TWOCLU GENXPT
C ORIGIN : H.FESEFELDT (02-DEC-1986)
C
#include "geant321/s_prntfl.inc"
#include "geant321/s_genio.inc"
#include "geant321/limits.inc"
C
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION WTMAX,WTMAXQ,WTFC,TWGT,ONE,TEXPXL,TEXPXU
      PARAMETER (ONE=1.D0)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (ONE=1.0)
#endif
      LOGICAL LZERO
      DIMENSION EMM(18)
      DIMENSION RNO(50)
      DIMENSION EM(18),PD(18),EMS(18),SM(18),FFQ(18),PCM1(90)
      EQUIVALENCE (NT,NPG),(AMASS(1),EM(1)),(PCM1(1),PCM(1,1))
      SAVE KNT
C
      DATA  FFQ/0.,3.141592, 19.73921, 62.01255, 129.8788, 204.0131,
     $                       256.3704, 268.4705, 240.9780, 189.2637,
     $                       132.1308,  83.0202,  47.4210,  24.8295,
     $                        12.0006,   5.3858,   2.2560,   0.8859/
      DATA  KNT , TWOPI /  1 , 6.2831853073 /
C
C --- Initialise local arrays and the result array PCM ---
      DO 10 JZERO=1,18
         PCM(1,JZERO)=0.
         PCM(2,JZERO)=0.
         PCM(3,JZERO)=0.
         PCM(4,JZERO)=0.
         PCM(5,JZERO)=0.
         EMM(JZERO)  =0.
         PD(JZERO)   =0.
         EMS(JZERO)  =0.
         SM(JZERO)   =0.
  10  CONTINUE
C
      KNT = KNT + 1
      IF (.NOT.NPRT(3).AND..NOT.NPRT(4)) GOTO 100
      WRITE(NEWBCD,1200) NPG,TECM,(AMASS(JK),JK=1,NPG)
  100 CONTINUE
  150 IF (NT .LT. 2)  GO TO 1001
      IF (NT .GT. 18)  GO TO 1002
      NTM1=NT-1
      NTM2=NT-2
      NTNM4 = 3*NT - 4
      EMM(1)=EM(1)
      TM=0.0
      DO 200 I=1,NT
      EMS(I)=EM(I)**2
      TM=TM+EM(I)
 200  SM(I)=TM
      WGT=1.
 210  TECMTM=TECM-TM
      IF (TECMTM .LE. 0.0)  GO TO 1000
      EMM(NT)=TECM
      IF (KGENEV.GT.1) GO TO 400
      EMMAX=TECMTM+EM(1)
      EMMIN=0.0
C
C          For weight calculation, form sum of log's of terms
C          instead of product of terms. Note that thereby WTMAX
C          and WTMAXQ are changed in their contents; they are
C          currently not used outside the range from here to
C          label 531. We also need to check for zero factors now.
C          Negative values cannot appear as GPDK always returns a
C          nonnegative number. As coded, even the exotic cases
C          NT<2 (first loop not executed) and NTM1<1 (second loop
C          not executed) should be safe and give the same result
C          for WTG in the end as the old code.
C
      WTMAX=0.0
      LZERO=.TRUE.
      DO 350 I=2,NT
      EMMIN=EMMIN+EM(I-1)
      EMMAX=EMMAX+EM(I)
      WTFC=GPDK(EMMAX,EMMIN,EM(I))
      IF(WTFC.LE.0.) THEN
      LZERO=.FALSE.
      GOTO 351
      ENDIF
      WTMAX=WTMAX+LOG(WTFC)
 350  CONTINUE
 351  WTMAXQ= EXPXU
      IF(LZERO) WTMAXQ= -WTMAX
      GO TO 455
  400 WTMAXQ=LOG(ONE*TECMTM**NTM2*FFQ(NT) / TECM)
  455 CONTINUE
      CALL GRNDM(RNO,NTNM4)
      IF(NTM2) 900,509,460
  460 CONTINUE
      CALL FLPSOR(RNO,NTM2)
      DO 508 J=2,NTM1
  508 EMM(J)=RNO(J-1)*TECMTM+SM(J)
  509 TWGT=WTMAXQ
      IR=NTM2
      LZERO=.TRUE.
      DO 530 I=1,NTM1
      PD(I)=GPDK(EMM(I+1),EMM(I),EM(I+1))
      IF(PD(I).LE.0.0) THEN
      LZERO=.FALSE.
      ELSE
      TWGT=TWGT+LOG(ONE*PD(I))
      ENDIF
  530 CONTINUE
  531 WGT=0.0
      IF(LZERO) THEN
      TEXPXU=EXPXU
      TEXPXL=EXPXL
      WGT=EXP(MAX(MIN(TWGT,TEXPXU),TEXPXL))
      ENDIF
      PCM(1,1)=0.0
      PCM(2,1)=PD(1)
      PCM(3,1)=0.0
      DO 570 I=2,NT
      PCM(1,I)=0.0
      PCM(2,I) = -PD(I-1)
      PCM(3,I)=0.0
      IR=IR+1
      BANG=TWOPI*RNO(IR)
      CB=COS(BANG)
      SB=SIN(BANG)
      IR=IR+1
      C=2.0*RNO(IR)-1.0
      S=SQRT(ABS(1.0-C*C))
      IF(I.EQ.NT) GO TO 1567
      ESYS=SQRT(PD(I)**2+EMM(I)**2)
      BETA=PD(I)/ESYS
      GAMA=ESYS/EMM(I)
      DO 568 J=1,I
      NDX = 5*J - 5
      AA= PCM1(NDX+1)**2 + PCM1(NDX+2)**2 + PCM1(NDX+3)**2
      PCM1(NDX+5) = SQRT(AA)
      PCM1(NDX+4) = SQRT(AA+EMS(J))
      CALL ROTES2(C,S,CB,SB,PCM,J)
      PSAVE = GAMA*(PCM(2,J)+BETA*PCM(4,J))
  568 PCM(2,J)=PSAVE
      GO TO 570
 1567 DO 1568 J=1,I
      AA=PCM(1,J)**2 + PCM(2,J)**2 + PCM(3,J)**2
      PCM(5,J)=SQRT(AA)
      PCM(4,J)=SQRT(AA+EMS(J))
      CALL ROTES2(C,S,CB,SB,PCM,J)
 1568 CONTINUE
  570 CONTINUE
  900 CONTINUE
      RETURN
 1000 DO 212 I=1,NPG
      PCM(1,I)=0.
      PCM(2,I)=0.
      PCM(3,I)=0.
      PCM(4,I)=AMASS(I)
  212 PCM(5,I)=AMASS(I)
      WGT=0.
      RETURN
 1001 IF(NPRT(3).OR.NPRT(4)) WRITE(NEWBCD,1101)
      GO TO 1050
 1002 WRITE(NEWBCD,1102)
 1050 WRITE(NEWBCD,1150) KNT
      WRITE(NEWBCD,1200) NPG,TECM,(AMASS(JK),JK=1,NPG)
      RETURN
 1100 FORMAT (1H0,'AVAILABLE ENERGY NEGATIVE')
 1101 FORMAT (1H0,'LESS THAN 2 OUTGOING PARTICLES')
 1102 FORMAT (1H0,'MORE THAN 18 OUTGOING PARTICLES')
 1150 FORMAT (1H0,'ABOVE ERROR DETECTED IN PHASP AT CALL NUMBER',I7)
 1200 FORMAT (1H0,'INPUT DATA TO PHASP.         NPG= ' ,I6/
     +2X,9H   TECM=  ,D15.7   ,18H  PARTICLE MASSES=,5D15.7/(42X,5D15.7)
     +)
      END
