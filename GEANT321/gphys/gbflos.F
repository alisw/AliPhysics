*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      FUNCTION GBFLOS(T,C)
C.    ******************************************************************
C.    *                                                                *
C.    *  It calculates a weight factor for the positron continuous     *
C.    *       bremsstrahlung energy loss.                              *
C.    *                                                                *
C.    *  Input:T kinetic energy (GeV)                                  *
C.    *        C cut for bremsstrahlung (BCUTE)                        *
C.    *        ETA ratio of positron to electron brems loss            *
C.    *  Output: weight factor                                         *
C.    *                                                                *
C.    *  ===>   positron loss = FLOS * electron loss   <===            *
C.    *                                                                *
C.    *  Author   L. Urban                                             *
C.    *                                                                *
C.    ******************************************************************
C.
      ETA = GCBREM(T)
      IF(ETA.LT.0.0001) THEN
         F=1.E-10
      ELSEIF(ETA.GT.0.9999) THEN
         F=1.
      ELSE
C
         E0=C/T
         IF(E0.GT.1.) E0=1.
C
         IF(E0.LT.1.E-8) THEN
            F=1.
         ELSE
            ALFA=1./ETA
            E1=1.-E0
            F=ETA*(1.-E1**ALFA)/E0
         ENDIF
C
      ENDIF
*
      GBFLOS=F
      END
