*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:55  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION CLGAMA(Z)
#include "gen/defc64.inc"
     +     CLGAMA
     +    ,CLOGAM
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION WLGAMA(Z)
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +     WLGAMA
     +    ,WLOGAM
#endif
#include "gen/defc64.inc"
C    +     Z,W,U,V,H,P,R,GCONJG,GCMPLX
     +     Z,  U,V,H,P,R
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CLGAMA')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CLGAMA/WLGAMA')
#endif
      DIMENSION C(10)

      PARAMETER (Z1 = 1, HF = Z1/2)

#include "gen/gcmpfun.inc"

      DATA PI /3.14159 26535 89793 24D+0/
      DATA C1 /9.18938 53320 46727 42D-1/
      DATA C2 /1.14472 98858 49400 17D+0/

      DATA C( 1) / 8.33333 33333 33333 33D-2/
      DATA C( 2) /-2.77777 77777 77777 78D-3/
      DATA C( 3) / 7.93650 79365 07936 51D-4/
      DATA C( 4) /-5.95238 09523 80952 38D-4/
      DATA C( 5) / 8.41750 84175 08417 51D-4/
      DATA C( 6) /-1.91752 69175 26917 53D-3/
      DATA C( 7) / 6.41025 64102 56410 26D-3/
      DATA C( 8) /-2.95506 53594 77124 18D-2/
      DATA C( 9) / 1.79644 37236 88305 73D-1/
      DATA C(10) /-1.39243 22169 05901 12D+0/
C     GREAL(U)=DREAL(U)
C     GIMAG(U)=DIMAG(U)
C     GCONJG(U)=DCONJG(U)
C     GCMPLX(X,Y)=DCMPLX(X,Y)

#if !defined(CERNLIB_DOUBLE)
      ENTRY CLOGAM(Z)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY WLOGAM(Z)
#endif

      X=Z
      Y=GIMAG(Z)
      IF(Y .EQ. 0 .AND. -ABS(X) .EQ. INT(X)) THEN
       H=0
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C306.1',ERRTXT)
      ELSE
       YA=ABS(Y)
       U=GCMPLX(X,YA)
       IF(X .LT. 0) U=1-U
       H=0
       UR=U
       IF(UR .LT. 7) THEN
        UI=GIMAG(U)
        A=ATAN2(UI,UR)
        H=U
        DO 1 I = 1,6-INT(UR)
        UR=UR+1
        U=GCMPLX(UR,UI)
        H=H*U
    1   A=A+ATAN2(UI,UR)
        H=GCMPLX(HF*LOG(GREAL(H)**2+GIMAG(H)**2),A)
        U=U+1
       ENDIF
       R=1/U**2
       P=R*C(10)
       DO 2 I = 9,2,-1
    2  P=R*(C(I)+P)
       H=C1+(U-HF)*LOG(U)-U+(C(1)+P)/U-H
       IF(X .LT. 0) THEN
        UR=INT(X)-1
        UI=PI*(X-UR)
        X=PI*YA
        T=EXP(-X-X)
        A=SIN(UI)
        T=X+HF*LOG(T*A**2+(HF*(1-T))**2)
        A=ATAN2(COS(UI)*TANH(X),A)-UR*PI
        H=C2-GCMPLX(T,A)-H
       ENDIF
       IF(Y .LT. 0) H=GCONJG(H)
      ENDIF
#if defined(CERNLIB_DOUBLE)
      WLGAMA=H
#endif
#if !defined(CERNLIB_DOUBLE)
      CLGAMA=H
#endif
      RETURN
  101 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END
