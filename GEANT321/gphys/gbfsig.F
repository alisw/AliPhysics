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
      FUNCTION GBFSIG(T,C)
C.    ******************************************************************
C.    *                                                                *
C.    *  It calculates a weight factor for the positron discrete       *
C.    *       (hard) bremsstrahlung cross section                      *
C.    *                                                                *
C.    *  Input:T kinetic energy (GeV)                                  *
C.    *        C cut for bremsstrahlung (BCUTE)                        *
C.    *        ETA ratio of positron to electron brems cross section   *
C.    *  Output: weight factor                                         *
C.    *                                                                *
C.    *  ===>  positron sigma = Fsig * electron sigma  <===            *
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
C
         IF(E0.GE.1.) THEN
            F=1.E-10
         ELSE
            ALFA=1./ETA
            E1=1.-E0
            F=ETA*E1**ALFA/E1
         ENDIF
      ENDIF
C
      GBFSIG=F
      END
