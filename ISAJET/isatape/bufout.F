#include "isajet/pilot.h"
      SUBROUTINE BUFOUT(IL)
C
C          INVERSE OF BUFIN.
C          IF CDCPACK IS USED, PACK TWO ZEVEL WORDS INTO EACH ZVOUT WORD
C          BY CALL PAIRPAK AND BUFFER OUT ZVOUT ARRAY.
C          OTHERWISE WRITE OUT ZEVEL.
C
C          NRECS INCREMENTED BY 1
C          IL SET TO 3


#include "isajet/itapes.inc"
#include "isajet/final.inc"
#include "isajet/idrun.inc"
#include "isajet/zevel.inc"
#include "isajet/zvout.inc"
      DIMENSION W(2),IW(2)
      EQUIVALENCE(W(1),IW(1))
C          THESE ARE NOW INITIALIZED IN BLOCK DATA
      ITA=IABS(ITEVT)
#if defined(CERNLIB_CDCPACK)
C          USE CDC ASSEMBLY LANGUAGE ROUTINE PAIRPAK TO PACK 2 ZEVEL
C          WORDS INTO 1 ZVOUT WORD AND THEN CALL CDC BUFFER OUT.
      IW(1)=IZEVEL(1)
      IZEVEL(2)=IL
      IW(2)=IL
      CALL PAIRPAK(W(1),W(2),WOUT,IFL)
      ZVOUT(1)=WOUT
      NW=IL/2+MOD(IL,2)
      DO 1 I=2,NW
      II=2*I-1
      CALL MOVLEV(IZEVEL(II),W,2)
      CALL PAIRPAK(W(1),W(2),WOUT,IFL)
      ZVOUT(I)=WOUT
    1 CONTINUE
      IL=3
      NRECS=NRECS+1
      BUFFER OUT(ITA,1) (ZVOUT(1),ZVOUT(NW))
      IF(UNIT(ITA,ZVOUT(1),ZVOUT(NW))) 10,10,11
#endif
#if defined(CERNLIB_STDIO)
C          STANDARD FORTRAN 77 WRITE.
      IZEVEL(2)=IL
      WRITE(ITA,ERR=11) (ZEVEL(I),I=1,IL)
      IL=3
      NRECS=NRECS+1
#endif
   10 RETURN
   11 WRITE(ITLIS,200) IEVT,IDG
  200 FORMAT(//'  BAD WRITE, EVENT NO.',I10,5X,'EVENT ID',2I10)
      RETURN
      END
