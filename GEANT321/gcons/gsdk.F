*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 02/08/94  18.13.20  by  S.Ravndal
*-- Author :
      SUBROUTINE GSDK(IPART,BRATIO,MODE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Defines branching ratios and decay modes for standard         *
C.    *  GEANT particles.                                              *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPART                                *
C.    *       Author    G.Patrick  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
      DIMENSION BRATIO(6)
      DIMENSION MODE(6)
C.
C.    ------------------------------------------------------------------
C.
      IF (IPART.LE.0)                              GO TO 99
      IF (IPART.GT.NPART)                          GO TO 99
C
C     Particle pointer.
C
      JPA = LQ(JPART-IPART)
      IF (JPA.LE.0)                                GO TO 99
C
C     Book decay links and data banks.
C
      JDK1=LQ(JPA-1)
      JDK2=LQ(JPA-2)
      IF(JDK1+JDK2.NE.0) THEN
         WRITE(CHMAIL, 10000)
         CALL GMAIL(1,0)
         CALL GPPART(IPART)
         CALL MZDROP(IXCONS,LQ(JPA-1),' ')
         CALL MZDROP(IXCONS,LQ(JPA-2),' ')
      ENDIF
      CALL MZBOOK(IXCONS,JDK1,JPA,-1,'PABR',0,0,6,3,0)
      JPA=LQ(JPART-IPART)
      CALL MZBOOK(IXCONS,JDK2,JPA,-2,'PAMO',0,0,6,2,0)
      JPA=LQ(JPART-IPART)
      JDK1=LQ(JPA-1)
      IQ(JDK1-5)=IPART
      IQ(JDK2-5)=IPART
C
C     Store branching ratios & decay modes.
C
      DO 20 I=1,6
         Q(JDK1+I) = BRATIO(I)
         IQ(JDK2+I) = MODE(I)
   20 CONTINUE
C
  99  RETURN
10000 FORMAT(' *** GSDK ***: Warning, redefinition of decay ',
     +       'for particle:')
      END
