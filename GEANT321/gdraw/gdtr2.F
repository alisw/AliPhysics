*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR2(IROOT,JMAX,LEVMAX,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Scan all JVOLUM structure                                *
C.    *                                                                *
C.    *       IROOT = root of the tree to be built up (input)          *
C.    *       JMAX = total number of nodes (output)                    *
C.    *       LEVMAX = maximum number of levels to be scanned          *
C.    *       IER  = error flag to detect nodes overflow (output)      *
C.    *                                                                *
C.    *    ==>Called by : GDTR0                                        *
C.    *       Authors : P.Zanarini, S.Banerjee   *********             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
C.
C.    ------------------------------------------------------------------
C.
      IER=0
C
      IQ(JNAM+1)=IROOT
      IQ(JMOT+1)=0
      IQ(JBRO+1)=0
      J=1
      JMAX=1
      LNOD=1
      LMAX=1
   10 CONTINUE
      IF (IQ(JSCA+J).EQ.0) THEN
         CALL GDTR1(J,JMAX,JNEW,IER)
         IF (IER.NE.0) GO TO 999
C
C             Fill array IQ(JDVM+I) for multiplicity by divisioning
C
         JJ=IQ(JMOT+J)
         IF (JJ.EQ.0) THEN
            IQ(JDVM+J)=0
         ELSE
            IVO=IUCOMP(IQ(JNAM+JJ),IQ(JVOLUM+1),NVOLUM)
            JVO=LQ(JVOLUM-IVO)
            NIN=Q(JVO+3)
            IF (NIN.LT.0) THEN
               JIN=LQ(JVO-1)
               IQ(JDVM+J)=ABS(Q(JIN+3))
            ELSE
               IQ(JDVM+J)=0
            ENDIF
         ENDIF
C
         JMAX=MAX(JMAX,JNEW)
         IF (JMAX.GT.MAXNOD) GO TO 20
      ENDIF
      J=J+1
*
***   IF (J.LE.JMAX) GO TO 10
*
      IF (J.LE.LMAX) THEN
         GO TO 10
      ELSE
         LNOD=LNOD+1
         IF (LNOD.GT.LEVMAX.AND.LEVMAX.GT.0) GO TO 999
         LMAX=JMAX
         IF (J.LE.JMAX) GO TO 10
      ENDIF
*
      GO TO 999
   20 WRITE (CHMAIL,1000) MAXNOD
      CALL GMAIL(0,0)
      IER=1
 1000 FORMAT (' GDTR2 : MAXNOD = ',I5,' - TREE OVERFLOW')
  999 RETURN
      END
