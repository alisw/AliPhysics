*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D501SF (VERS,SUB,MODE,M,A,N,F,DF,K,NX,X,Y,SY,W,NERROR)

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE COMPUTES FUNCTION VALUES AND DERIVATIVES DEPENDING ON
*   THE VALUE OF THE PARAMETER  VERS.
*
*************************************************************************

#include "gen/imp64.inc"
      DIMENSION A(*),F(*),DF(N,*),X(*),Y(*),SY(*),W(*)
      CHARACTER VERS*6

      NERROR=0

      IF (VERS .EQ. 'DSUMSQ') THEN
       CALL SUB (M,A,N,F,DF,MODE,NERROR)
       IF (NERROR .NE. 0) NERROR=3
       RETURN
      ENDIF

      IF (VERS .EQ. 'DFUNFT') THEN
       IX=1
       DO 20 I=1,N
        CALL SUB (K,X(IX),M,A,SF,W,MODE,NERROR)
        IF (NERROR .NE. 0) THEN
         NERROR=3
         RETURN
        ENDIF
        F(I)=(Y(I)-SF)/SY(I)
        IX=IX+NX
       IF (MODE .EQ. 0) GOTO 20
       DO 10 J=1,M
   10  DF(I,J)=-W(J)/SY(I)
   20  CONTINUE
       RETURN
      ENDIF

      END



