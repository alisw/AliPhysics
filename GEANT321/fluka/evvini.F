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
*$ CREATE EVVINI.FOR
*COPY EVVINI
*
*=== evvini ===========================================================*
*
      SUBROUTINE EVVINI (WHAT,SDUM)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Version march 1990   by          Alfredo Ferrari, INFN-Milan     *
*                                                                      *
*     Last change  on  20-apr-92   by  Alfredo Ferrari, INFN-Milan     *
*     This subroutine initialises the quark model subroutines and the  *
*     evaporation module (Evap5 from Hetc-Kfa).                        *
*                                                                      *
*     input variables:                                                 *
*        what has currently no meaning                                 *
*        sdum = material name of a heavy ion beam                      *
*                                                                      *
*     output variables: none                                           *
*                                                                      *
*     Common blocks:                                                   *
*                                                                      *
*     /print/ printing control for the nucevt package                  *
*     /inpdat/ parameters for the nucevt package                       *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/inpdat2.inc"
#include "geant321/parevt.inc"
      CHARACTER*8 SDUM
      DIMENSION WHAT(6)
      COMMON/FKPRIN/IPRI,INIT
      IF (.NOT.LEVPRT) LDEEXG=.FALSE.
      CALL DATAR3
      CALL HADDEN
      CALL RCHANV
      CALL CHANWT
      ONEDUM = ONEONE
      ZERDUM = ZERZER
      NUD = NUDISV ( ONEDUM, 1, ZERDUM, ZERDUM, DPOWER, ZERDUM )
      CALL BERTTP
      CALL INCINI
      IPRI=0
      INIT=0
      LTBAMJ=0
      ASBAMJ = 0.5D+00
      B8BAMJ = 0.4D+00
      RETURN
      END
