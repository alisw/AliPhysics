*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE BKLASS.FOR
*COPY BKLASS
*
*=== bklass ===========================================================*
*
      SUBROUTINE BKLASS(I,J,K,I8,I10)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C*** I,J,K QUARK FLAVOURS U,D,S,C=1,2,3,4
C*** AQ = -Q
C*** I8,I10 BARYON INDICES
#include "geant321/inpdat.inc"
      IF (I) 10,10,100
C*** BARYON
100   CONTINUE
      CALL INDEX2(J,K,IND)
      I8=IB08(I,IND)
      I10=IB10(I,IND)
         IF (I8.LE.0) I8=I10
      RETURN
10    CONTINUE
C*** ANTIBARYONS
      II=IABS(I)
      JJ=IABS(J)
      KK=IABS(K)
      CALL INDEX2(JJ,KK,IND)
      I8=IA08(II,IND)
      I10=IA10(II,IND)
         IF (I8.LE.0) I8=I10
      RETURN
      END
