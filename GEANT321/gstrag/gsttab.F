*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:39  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GSTTAB(GAMMA,NT,EN,FN)
C   ********************************************************************
C   *                                                                  *
C   *  dN^2/dEdx table for a gamma factor                              *
C   *                                                                  *
C   *  ==> Called by : GSTINI                                          *
C   *      Authors   : P. Nevski                                       *
C   *                                                                  *
C   ********************************************************************
#include "geant321/gcstra.inc"
      EXTERNAL GSTDN
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION DGAMM, DGAM2, DBET2
#endif
 
      DIMENSION EN(1),FN(1)
      PARAMETER (RENER=0.05)
*
*     ------------------------------------------------------------------
*
C   ----  final tables for dN/dx integral vs Energy for a gamma   ------
      DGAMM = GAMMA
      DGAM2 = DGAMM**2
      GAMMA2= DGAM2
      DBET2 = 1-1/DGAM2
      BETA2 = DBET2
      S=0
      NT=0
      ES=0
      EI=EMAX
C   Loop over I=NTAB,1,-1
      DO 10 I=NTAB,1,-1
         EO=EI
         EI=EMM(I)
         S=S+GXGINT(GSTDN,EI,EO,EPPS)
C      Skip Unless ABS(ES-EI)>RENER \ I==1
         IF ((ABS(ES-EI).LE.RENER).AND.(I.NE.1))GO TO 10
         ES=EI
         NT=NT+1
         EN(NT)=EI
         FN(NT)=S
   10 CONTINUE
      DNDX=S
      END
