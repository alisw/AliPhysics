*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDSCAL(U0,V0)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw an unitary length segment at position (U0,V0)       *
C.    *       correspondent to current scale factor GSCU               *
C.    *       and zoom factor GZUA                                     *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDSPEC                      *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
      DIMENSION U(2),V(2)
      CHARACTER*8 IUMC
C.
C.    ------------------------------------------------------------------
C.
      IF (GSCU*GZUA.LE.0.03) THEN
         IUMC='100< CM$'
         UML=100.
      ELSE IF (GSCU*GZUA.LE.0.3) THEN
         IUMC='10 <CM $'
         UML=10.
      ELSE IF (GSCU*GZUA.LE.3.) THEN
         IUMC=' 1 <CM $'
         UML=1.
      ELSE IF (GSCU*GZUA.LE.30.) THEN
         IUMC='0.1<CM $'
         UML=0.1
      ELSE IF (GSCU*GZUA.LE.300.) THEN
         IUMC='100< [M$'
         UML=0.01
      ELSE IF (GSCU*GZUA.LE.3000.) THEN
         IUMC='10 <[M $'
         UML=0.001
      ELSE
         IUMC=' 1<[M $'
         UML=0.0001
      ENDIF
C
      SLEN=UML*GSCU*GZUA
C
C             Normalize to PLTRNX,PLTRNY
C
      SLEN=SLEN*MIN(PLTRNX,PLTRNY)/20.
C
      DX=SLEN/2.
      DY=DX/4.
      U(1)=U0-DX
      U(2)=U0+DX
      V(1)=V0
      V(2)=V0
      CALL GDRAWV(U,V,2)
      U(2)=U0-DX
      V(1)=V0+DY
      V(2)=V0-DY
      CALL GDRAWV(U,V,2)
      U(1)=U0+DX
      U(2)=U0+DX
      CALL GDRAWV(U,V,2)
C
      CSIZ=SLEN/4.
      VV0=V0+3.*DY
      CALL GDRAWT(U0,VV0,IUMC,CSIZ,0.,1,0)
C
      END
