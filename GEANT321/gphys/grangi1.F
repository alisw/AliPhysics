*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:33  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_NEVER)
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GRANGI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates the stopping range , in cm.                        *
C.    *     (i.e. the maximum step due to the continuous energy loss)  *
C.    *      The stopping range is the integral of the inverse of the  *
C.    *      DE/DX table.                                              *
C.    *     There are 4 tables : electron,positron,muon,proton         *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Author      M.Maire   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
*
      DO 50 IEKBIN=1,NEK1
*
* *** Electrons
*
         JLOSS = LQ(JMA-1)
         JRANG = LQ(JMA-15)
         IF (Z.LT.1.) THEN
            Q(JRANG+IEKBIN) = BIG
         ELSE IF ((IEKBIN.EQ.1).OR.(Q(JLOSS+IEKBIN-1).LE.0.)) THEN
            Q(JRANG+IEKBIN) = 0.
         ELSE
            DEBIN=ELOW(IEKBIN) - ELOW(IEKBIN-1)
            FRACT= Q(JLOSS+IEKBIN)/Q(JLOSS+IEKBIN-1) - 1.
            IF (ABS(FRACT).GE.0.1) THEN
               DRANG=DEBIN/(Q(JLOSS+IEKBIN)-Q(JLOSS+IEKBIN-1))
               DRANG=DRANG*LOG(1.+FRACT)
            ELSE
               DRANG=DEBIN/Q(JLOSS+IEKBIN-1)
               DRANG= DRANG*(1. - 0.5*FRACT + 0.3333*FRACT*FRACT)
            ENDIF
            Q(JRANG+IEKBIN)=Q(JRANG+IEKBIN-1) + DRANG
         ENDIF
*
* *** Positons
*
         JLOSS = LQ(JMA-1) + NEK1
         JRANG = LQ(JMA-15) + NEK1
         IF (Z.LT.1.) THEN
            Q(JRANG+IEKBIN) = BIG
         ELSE IF ((IEKBIN.EQ.1).OR.(Q(JLOSS+IEKBIN-1).LE.0.)) THEN
            Q(JRANG+IEKBIN) = 0.
         ELSE
            DEBIN=ELOW(IEKBIN) - ELOW(IEKBIN-1)
            FRACT= Q(JLOSS+IEKBIN)/Q(JLOSS+IEKBIN-1) - 1.
            IF (ABS(FRACT).GE.0.1) THEN
               DRANG=DEBIN/(Q(JLOSS+IEKBIN)-Q(JLOSS+IEKBIN-1))
               DRANG=DRANG*LOG(1.+FRACT)
            ELSE
               DRANG=DEBIN/Q(JLOSS+IEKBIN-1)
               DRANG= DRANG*(1. - 0.5*FRACT + 0.3333*FRACT*FRACT)
            ENDIF
            Q(JRANG+IEKBIN)=Q(JRANG+IEKBIN-1) + DRANG
         ENDIF
*
* *** Muons
*
         JLOSS = LQ(JMA-2)
         JRANG = LQ(JMA-16)
         IF (Z.LT.1.) THEN
            Q(JRANG+IEKBIN) = BIG
         ELSE IF ((IEKBIN.EQ.1).OR.(Q(JLOSS+IEKBIN-1).LE.0.)) THEN
            Q(JRANG+IEKBIN) = 0.
         ELSE
            DEBIN=ELOW(IEKBIN) - ELOW(IEKBIN-1)
            FRACT= Q(JLOSS+IEKBIN)/Q(JLOSS+IEKBIN-1) - 1.
            IF (ABS(FRACT).GE.0.1) THEN
               DRANG=DEBIN/(Q(JLOSS+IEKBIN)-Q(JLOSS+IEKBIN-1))
               DRANG=DRANG*LOG(1.+FRACT)
            ELSE
               DRANG=DEBIN/Q(JLOSS+IEKBIN-1)
               DRANG= DRANG*(1. - 0.5*FRACT + 0.3333*FRACT*FRACT)
            ENDIF
            Q(JRANG+IEKBIN)=Q(JRANG+IEKBIN-1) + DRANG
         ENDIF
*
* *** Protons
*
         JLOSS = LQ(JMA-3)
         JRANG = LQ(JMA-16) + NEK1
         IF (Z.LT.1.) THEN
            Q(JRANG+IEKBIN) = BIG
         ELSE IF ((IEKBIN.EQ.1).OR.(Q(JLOSS+IEKBIN-1).LE.0.)) THEN
            Q(JRANG+IEKBIN) = 0.
         ELSE
            DEBIN=ELOW(IEKBIN) - ELOW(IEKBIN-1)
            FRACT= Q(JLOSS+IEKBIN)/Q(JLOSS+IEKBIN-1) - 1.
            IF (ABS(FRACT).GE.0.1) THEN
               DRANG=DEBIN/(Q(JLOSS+IEKBIN)-Q(JLOSS+IEKBIN-1))
               DRANG=DRANG*LOG(1.+FRACT)
            ELSE
               DRANG=DEBIN/Q(JLOSS+IEKBIN-1)
               DRANG= DRANG*(1. - 0.5*FRACT + 0.3333*FRACT*FRACT)
            ENDIF
            Q(JRANG+IEKBIN)=Q(JRANG+IEKBIN-1) + DRANG
         ENDIF
   50 CONTINUE
*
      END
#endif
