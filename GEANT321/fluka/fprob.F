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
*$ CREATE FPROB.FOR
*COPY FPROB
*
*=== fprob ============================================================*
*
      FUNCTION FPROB (Z,A,U)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  Fission PROBability:                                                *
*                                                                      *
*  Created  on  8 april 1993   by     Alfredo Ferrari & Paola Sala     *
*                                                INFN - Milan          *
*                                                                      *
*  Last change  on  09-apr-93  by  Alfredo Ferrari, INFN - Milan       *
*                                                                      *
*  This function is just a translation into Fluka of the same function *
*  of LAHET kindly provided by R.E. Prael                              *
*                                                                      *
*     Input variables:                                                 *
*                       Z = present nucleus atomic number              *
*                       A =    "       "     mass    "                 *
*                       U =    "       "    excitation energy (MeV)    *
*                                                                      *
*----------------------------------------------------------------------*
*
      FPROB = 0
      STOP 'FPROB'
*      RETURN
*=== End of function Fprob ============================================*
      END
