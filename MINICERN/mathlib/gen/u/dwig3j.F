*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DWIG3J(A,B,C,X,Y,Z)
#include "gen/imp128.inc"
      CHARACTER*6 NAME
      DWIG3J=0
      NAME  ='DWIG3J'
      CALL MTLPRT(NAME,'U111',
     +'not available on this machine - see documentation')
      RETURN
 
      ENTRY DCLEBG(A,B,C,X,Y,Z)
      DCLEBG=0
      NAME  ='DCLEBG'
      CALL MTLPRT(NAME,'U111',
     +'not available on this machine - see documentation')
      RETURN
 
      ENTRY DWIG6J(A,B,C,X,Y,Z)
      DWIG6J=0
      NAME  ='DWIG6J'
      CALL MTLPRT(NAME,'U111',
     +'not available on this machine - see documentation')
      RETURN
 
      ENTRY DRACAW(A,B,C,X,Y,Z)
      DRACAW=0
      NAME  ='DRACAW'
      CALL MTLPRT(NAME,'U111',
     +'not available on this machine - see documentation')
      RETURN
 
      ENTRY DJAHNU(A,B,C,X,Y,Z)
      DJAHNU=0
      NAME  ='DJAHNU'
      CALL MTLPRT(NAME,'U111',
     +'not available on this machine - see documentation')
      RETURN
 
      ENTRY DWIG9J(A,B,C,P,Q,R,X,Y,Z)
      DWIG9J=0
      NAME  ='DWIG9J'
      CALL MTLPRT(NAME,'U111',
     +'not available on this machine - see documentation')
      RETURN
      END
#endif
