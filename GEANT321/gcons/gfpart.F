*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:14  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GFPART(IPART,NAPART,ITRTYP,AMASS,CHARGE,TLIFE
     +          ,UBUF,NWBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Return parameters for particle of type IPART             *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GTREVE                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      CHARACTER*20 NAPART
      DIMENSION UBUF(10)
C.
C.    ------------------------------------------------------------------
C.
      ITRTYP=0
      IF (JPART.LE.0)GO TO 99
      IF (IPART.LE.0)GO TO 99
      IF (IPART.GT.NPART)GO TO 99
      JPA = LQ(JPART- IPART)
      IF (JPA.LE.0)GO TO 99
C
      CALL UHTOC(IQ(JPA+1),4,NAPART,20)
      ITRTYP = Q(JPA + 6)
      AMASS = Q(JPA + 7)
      CHARGE = Q(JPA + 8)
      TLIFE = Q(JPA + 9)
      NWBUF = IQ(JPA-1) - 9
      IF(NWBUF.LE.0) GO TO 99
      DO 20 I=1,NWBUF
  20  UBUF(I)=Q(JPA+9+I)
C
  99  RETURN
      END
