*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GFDETH(IUSET,IUDET,NH,NAMESH,NBITSH,ORIG,FACT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Returns hit parameters for detector IUDET of set IUSET.      *
C.    *                                                                *
C.    *             Input parameters                                   *
C.    * IUSET     user set identifier                                  *
C.    * IUDET     user detector identifier                             *
C.    *                                                                *
C.    *             Output parameters                                  *
C.    * NH        number of elements per hit                           *
C.    * NAMESH    the NH variable names for the hit elements           *
C.    * NBITSH    the NH bit numbers for packing the variable values   *
C.    * ORIG      The quantity packed in the structure JHITS for the   *
C.    *           Ith variable is a  positive integer with NBITSH(I)   *
C.    *           bits and such that                                   *
C.    * FACT      IVAR(I) = (VAR(I)+ORIG(I))*FACT(I)                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun ,M.Maire *********                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      DIMENSION NBITSH(1),ORIG(1),FACT(1)
      CHARACTER*4 IUSET,IUDET,NAMESH(1)
C.
C.    ------------------------------------------------------------------
C.
      IF(JSET.LE.0)GO TO 90
      NSET=IQ(JSET-1)
      IF(NSET.LE.0)GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF(ISET.LE.0)GO TO 90
      JS=LQ(JSET-ISET)
      NDET=IQ(JS-1)
      IF(NDET.LE.0)GO TO 90
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF(IDET.LE.0)GO TO 95
      JD=LQ(JS-IDET)
      JDH=LQ(JD-1)
      NH=IQ(JD+4)
C
      IF(NH.GT.0)THEN
         DO 10 I=1,NH
            CALL UHTOC(IQ(JDH+4*I-3),4,NAMESH(I),4)
            NBITSH(I)=IQ(JDH+4*I-2)
            ORIG(I)  = Q(JDH+4*I-1)
            FACT(I)  = Q(JDH+4*I )
  10     CONTINUE
      ENDIF
      GO TO 99
C
C              Error
C
   90 WRITE (CHMAIL, 1000) IUSET
      CALL GMAIL(0,0)
 1000 FORMAT (' ***** GFDETH ERROR SET ',A4,' NOT FOUND')
      GO TO 99
   95 WRITE (CHMAIL, 2000) IUSET,IUDET
      CALL GMAIL(0,0)
 2000 FORMAT (' ***** GFDETH ERROR FOR SET ',A4,
     + ' DETECTOR ',A4,' NOT FOUND')
C
  99  RETURN
      END
