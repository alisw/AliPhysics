*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ASSNDX(MODE,A,N,M,IDA,K,SUM,IW,IDW)
      LOGICAL LSW
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'ASSNDX')

      DIMENSION A(IDA,*),K(*),IW(IDW,*)

      IF(N .LT. 1 .OR. M .LT. 1) THEN
       WRITE(ERRTXT,101) N,M
       CALL MTLPRT(NAME,'H301.1',ERRTXT)
       RETURN
      ENDIF
      IMAX=MAX(N,M)
      IMIN=MIN(N,M)
      SUM=0
      IF(N .LE. M) THEN
       DO 1 I = 1,N
       RMIN=A(I,1)
       DO 2 J = 1,M
    2  RMIN=MIN(RMIN,A(I,J))
       SUM=SUM+RMIN
       DO 3 J = 1,M
    3  A(I,J)=A(I,J)-RMIN
    1  CONTINUE
      ENDIF
      IF(N .GE. M) THEN
       DO 4 J = 1,M
       RMIN=A(1,J)
       DO 5 I = 1,N
    5  RMIN=MIN(RMIN,A(I,J))
       SUM=SUM+RMIN
       DO 7 I = 1,N
    7  A(I,J)=A(I,J)-RMIN
    4  CONTINUE
      ENDIF
      DO 8 I = 1,IMAX
      K(I)=0
    8 IW(I,1)=0

      DO 12 I = 1,N
      DO 13 J = 1,M
      IF(A(I,J)+IW(J,1) .EQ. 0) THEN
       K(I)=J
       IW(J,1)=I
       GO TO 12
      ENDIF
   13 CONTINUE
   12 CONTINUE

   10 IFLAG=N
      IRL=0
      ICL=0
      IRS=1

      DO 11 I = 1,N
      IW(I,5)=0
      IF(K(I) .EQ. 0) THEN
       IRL=IRL+1
       IW(IRL,6)=I
       IW(I,5)=-1
       IFLAG=IFLAG-1
      ENDIF
   11 CONTINUE
      IF(IFLAG .EQ. IMIN) THEN
       IF(MODE .EQ. 2) THEN
        DO 70 I = 1,IMAX
   70   K(I)=IW(I,1)
       ENDIF
       RETURN
      ENDIF

      DO 14 J = 1,M
   14 IW(J,4)=0

   30 I=IW(IRS,6)
      IRS=IRS+1
      DO 31 J = 1,M
      IF(A(I,J)+IW(J,4) .EQ. 0) THEN
       IW(J,4)=I
       ICL=ICL+1
       IW(ICL,2)=J
       NEW=IW(J,1)
       IF(NEW .EQ. 0) THEN
        J1=J
   61   IW(J1,1)=IW(J1,4)
        I=IW(J1,4)
        IF(K(I) .EQ. 0) THEN
         K(I)=J1
         GO TO 10
        ENDIF
        JSV=J1
        J1=K(I)
        K(I)=JSV
        GO TO 61
       ENDIF
       IRL=IRL+1
       IW(IRL,6)=NEW
       IW(NEW,5)=I
      ENDIF
   31 CONTINUE
      IF(IRS .LE. IRL) GO TO 30

      LSW=.TRUE.
      ICL0=ICL
      ICBL=0
      DO 41 J = 1,M
      IF(IW(J,4) .EQ. 0) THEN
       ICBL=ICBL+1
       IW(ICBL,3)=J
      ENDIF
   41 CONTINUE
      RMIN=A(IW(1,6),IW(1,3))
      DO 42 I = 1,IRL
      DO 42 J = 1,ICBL
   42 RMIN=MIN(RMIN,A(IW(I,6),IW(J,3)))
      SUM=SUM+RMIN*(IRL+ICBL-IMAX)

      DO 44 I = 1,N
      IF(IW(I,5) .EQ. 0) THEN
       DO 49 IPP = 1,ICL0
   49  A(I,IW(IPP,2))=A(I,IW(IPP,2))+RMIN
       GO TO 44
      ENDIF
      DO 45 IPP = 1,ICBL
      NEW=IW(IPP,3)
      A(I,NEW)=A(I,NEW)-RMIN
      IF(LSW .AND. A(I,NEW)+IW(NEW,4) .EQ. 0) THEN
       IW(NEW,4)=I
       IF(IW(NEW,1) .EQ. 0) THEN
        J1=NEW
        LSW=.FALSE.
       ELSE
        ICL=ICL+1
        IW(ICL,2)=NEW
        IRL=IRL+1
        IW(IRL,6)=IW(NEW,1)
       END IF
      END IF
   45 CONTINUE
   44 CONTINUE
      IF(LSW) THEN
       DO 51 I = ICL0+1,ICL
   51  IW(IW(IW(I,2),1),5)=IW(I,2)
       GO TO 30
      ELSE
   62  IW(J1,1)=IW(J1,4)
       I=IW(J1,4)
       IF(K(I) .EQ. 0) THEN
        K(I)=J1
        GO TO 10
       ENDIF
       JSV=J1
       J1=K(I)
       K(I)=JSV
       GO TO 62
      ENDIF
  101 FORMAT('N = ',I5,' < 1   OR   M = ',I5,' < 1')
      END
