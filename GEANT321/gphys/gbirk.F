*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GBIRK(EDEP)
*
************************************************************************
*                                                                      *
* *** apply BIRK's saturation law to energy deposition ***             *
* *** only organic scintillators implemented in this version MODEL=1   *
*                                                                      *
* Note : the material is assumed ideal, which means that impurities    *
*        and aging effects are not taken into account                  *
*                                                                      *
* algorithm : edep = destep / (1. + RKB*dedx + C*(dedx)**2)            *
*                                                                      *
* the values of RKB and C can be entered via :                         *
*                                                                      *
* call gstpar(imate,'BIRK1',value) to set the model (must be 1 or 2)   *
* call gstpar(imate,'BIRK2',value) to set RKB                          *
* call gstpar(imate,'BIRK3',value) to set C                            *
*                                                                      *
* the basic units of the coefficient are g/(Mev*cm**2)                 *
* because the de/dx is obtained in Mev/cm                              *
*                                                                      *
* exp. values from NIM 80 (1970) 239-244 :                             *
*                                                                      *
* RKB = 0.013  g/mev*cm**2  and  C = 9.6e-6  g**2/(Mev**2)(cm**4)      *
*                                                                      *
*    ==>Called by : <USER>,GUSTEP                                      *
*       Author    R.Brun , N.Van Eijndhoven    *********               *
*                                                                      *
************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcmate.inc"
#include "geant321/gctrak.inc"
C.
C.    ------------------------------------------------------------------
C.
      EDEP=DESTEP
C
C --- no saturation law for neutral particles ---
      IF (ABS(CHARGE) .LE. 1.E-10) GO TO 99
C
C --- get the values for the BIRK coefficients stored via GSTPAR ---
      JTMN=LQ(JTM)
      IF (JTMN .LE. 0)GO TO 99
      MODEL=Q(JTMN+27)
      IF(MODEL.GT.2)GO TO 99
      RKB=Q(JTMN+28)/DENS
      C=Q(JTMN+29)/DENS**2
C
C --- correction for particles with more than 1 charge unit ---
C --- based on alpha particle data (only apply for MODEL=1) ---
      IF(MODEL.EQ.1)THEN
         IF (ABS(CHARGE) .GE. 2.) RKB=RKB*7.2/12.6
      ENDIF
C
      DEDXCM=1000.*Q(JLOSS+IEKBIN)
C
      EDEP=DESTEP/(1.+RKB*DEDXCM+C*DEDXCM**2)
C
  99  RETURN
      END
 
