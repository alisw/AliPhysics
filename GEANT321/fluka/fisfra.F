*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE FISFRA.FOR
*COPY FISFRA
*
*=== fisfra ===========================================================*
*
      SUBROUTINE FISFRA ( JA, JZ, U, EREC, UMO, GAMCM, ETACM )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  FISsion FRAgments emission:                                         *
*                                                                      *
*                                                                      *
*      Subroutine to pick post fission parameters for nucleus          *
*      JZ,JA excited to U and recoiling with Erec.                     *
*                                                                      *
*  Created  on  08 april 1993  by  A.Ferrari & P.Sala, INFN - Milan    *
*                                                                      *
*  Last change  on  14-apr-93  by  Alfredo Ferrari, INFN - Milan       *
*                                                                      *
*  This routine is just the FLUKA equivalent of the routines FISSED and*
*  FISDIS of LAHET, taken thanks to R.E.Prael                          *
*                                                                      *
*  Input variables:                                                    *
*     JA = Mass number of the fissioning nucleus                       *
*     JZ = Atomic number of the fissioning nucleus                     *
*     U  = Excitation energy (MeV) of the fissioning nucleus           *
*     Erec = Recoil kinetic energy (MeV) of the residual nucleus       *
*            The recoil direction is given by Coslf0 (i)               *
*     Umo  = invariant mass of the fissioning system (mass+excitation) *
*     Gamcm= CMS Lorentz boost "gamma"                                 *
*     Etacm= CMS Lorentz boost "eta"                                   *
*                                                                      *
*----------------------------------------------------------------------*
*
      STOP 'FISFRA'
*      RETURN
*=== Fisfra ===========================================================*
      END
