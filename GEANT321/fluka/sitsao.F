*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:05  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE SITSAO.FOR
*COPY SITSAO
*
*=== sitsao ===========================================================*
*
      FUNCTION SITSAO ( E, IZ, A )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Original from A. Fasso`, CERN-TIS                                *
*                                                                      *
*     Created on 20 september 1991 by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 20-sep-91     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
*----------------------------------------------------------------------*
*   Calculates proton-nucleus inelastic cross-section by the formula   *
*   of Letaw, Silberberg and Tsao, Astrophys. J. Suppl. 51, 271 (1983) *
*----------------------------------------------------------------------*
*   E = Energy in MeV, IZ = Atomic number, A = Mass number (At. weight)*
*----------------------------------------------------------------------*
*                                                formula (2)
      SITSAO = 45.D+00 * A**0.7D+00 * ( 1.D+00 + 0.016D+00 * SIN
     &       ( 5.3D+00 - 2.63D+00 * LOG(A)) )
*                                                formula (5)
      IF(E .LT. 2000.D+00) SITSAO = SITSAO * ( 1.D+00 - 0.62D+00
     &                            / EXP ( E / 200.D+00 )
     &                            * SIN ( 10.9D+00 / E**0.28D+00 ) )
      IF (IZ .EQ. 2) THEN
         SITSAO=SITSAO*0.8D+00
      ELSE IF (IZ .EQ. 4) THEN
         SITSAO=SITSAO*(1.D+00 + 0.75D+00/EXP(E/75.D+00))
      END IF
      RETURN
      END
