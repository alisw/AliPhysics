*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
      SUBROUTINE ISOTPE(D,LD,KM,RHO,IN,IDICTS,LDICT,E,TSIG,NMED,
     +                  IIN,IIM)
C       THIS ROUTINE DETERMINES WHICH ISOTOPE HAS BEEN STRUCK
C       IN MEDIA NMED
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mmicab.inc"
C
      DIMENSION D(*),LD(*),KM(*),RHO(*),IN(*),IDICTS(NNR,NNUC),
     +          LDICT(NNR,NNUC)
      SAVE
C
      R=FLTRNF(0)
      NOA=0
      SUM=0.
#if defined(CERNLIB_MDEBUG)
      DO 10 K=1,NMIX
         IF(KM(K).NE.NMED)GO TO 10
C       DETERMINE ISOTOPE NUMBER
         K1=IN(K)
         K2=K
C       DETERMINE TOTAL CROSS SECTION FOR THIS ISOTOPE
         LS1=IDICTS(1,K1)+LMOX2
         L1=LDICT(1,K1)
         LEN=L1/2
         CALL TBSPLT(D(LS1),E,LEN,X)
         SUM=SUM+X*RHO(K)
         PRINT *,' ISOTPE: K=',K,' RHO=',RHO(K),' Sig=',X*RHO(K),
     +       ' SUM=',SUM,' TSIG=',TSIG,' R=',R
   10 CONTINUE
      SUM = 0.0
#endif
   20 DO 30 K=1,NMIX
         IF(KM(K).NE.NMED)GO TO 30
C       DETERMINE ISOTOPE NUMBER
         K1=IN(K)
         K2=K
C       DETERMINE TOTAL CROSS SECTION FOR THIS ISOTOPE
         LS1=IDICTS(1,K1)+LMOX2
         L1=LDICT(1,K1)
         LEN=L1/2
         CALL TBSPLT(D(LS1),E,LEN,X)
         SUM=SUM+X*RHO(K)
C       CHECK TO SEE IF THIS ISOTOPE WAS HIT
         IF(R.LE.SUM/TSIG)GO TO 40
   30 CONTINUE
C       AN ISOTOPE WAS NOT CHOSEN, TRY AGAIN
      NOA=NOA+1
      IF(NOA.GT.5)GO TO 50
      SUM=0.0
      R=FLTRNF(0)
      GO TO 20
   40 IIN=K1
      IIM=K2
#if defined(CERNLIB_MDEBUG)
      PRINT *,' Isotope chosen : K=',K
#endif
      RETURN
   50 WRITE(IOUT,10000)NMED,TSIG
10000 FORMAT(' MICAP: AN ISOTOPE WAS NOT CHOSEN IN 5 ATTEMPTS IN ',
     +'ROUTINE ISOTPE',/,3X,'MEDIUM=',I5,5X,'MACROSCOPIC XSEC=',
     +1PE12.4)
      WRITE(IOUT,10100)R,SUM,TSIG,X,E,RHO(K2),NMED,K1,K2
10100 FORMAT('0',1X,1P6E12.4,3I10)
      WRITE(6,*) ' CALOR: ERROR in ISOTPE =====> STOP '
      STOP
      END
