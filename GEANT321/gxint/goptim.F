*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GOPTIM(IOPT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Changes the geometry optimization level.                    *
C.    *    Reprocess the geometry initialization phase.                *
C.    *                                                                *
C.    *    Called by : <GXINT>                                         *
C.    *    Authors   : R.Brun  *********                               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcopti.inc"
C.    ------------------------------------------------------------------
*
      IOPTIM=IOPT
      IF(JVOLUM.LE.0)GO TO 999
      JN=LQ(JVOLUM)
      IF(JN.GT.0)CALL MZDROP(IXCONS,JN,' ')
      DO 20 IVO = 1,NVOLUM
         JVO = LQ(JVOLUM-IVO)
         IQ(JVO) = IBCLR(IQ(JVO),0)
         NIN   = Q(JVO+3)
         IF(NIN.GT.0) THEN
            ISEARC = Q(JVO+1)
            IF(ISEARC.LT.0) THEN
               Q(JVO+1)=-Q(LQ(LQ(JVO-NIN-1))+1)
            ENDIF
            NINL = IQ(JVO-2)
            DO 10 I = NIN+1,NINL
               JN=LQ(JVO-I)
               IF(JN.GT.0)CALL MZDROP(IXCONS,JN,' ')
   10       CONTINUE
         ENDIF
   20 CONTINUE
*
      CALL GGCLOS
  999 END
