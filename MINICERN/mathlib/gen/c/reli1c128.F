*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:11  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DELI1C(AKP)
C
#include "gen/imp128.inc"
C
      DELI1C=0
      CALL MTLPRT('DELI1C','C347.0','Not available on this machine')
      RETURN

      ENTRY DELI2C(AKP,A,B)
      DELI2C=0
      CALL MTLPRT('DELI2C','C347.0','Not available on this machine')
      RETURN

      ENTRY DELI3C(AKP,AK2,P)
      DELI3C=0
      CALL MTLPRT('DELI3C','C347.0','Not available on this machine')
      RETURN

      ENTRY DELIGC(AKP,P,A,B)
      DELIGC=0
      CALL MTLPRT('DELIGC','C347.0','Not available on this machine')
      RETURN

      ENTRY DELIKC(AKP,P,A,B)
      ENTRY DELLIK(AKP,P,A,B)
      DELIKC=0
      DELLIK=0
      CALL MTLPRT('DELIKC','C347.0','Not available on this machine')
      RETURN

      ENTRY DELIEC(AKP,P,A,B)
      ENTRY DELLIE(AKP,P,A,B)
      DELIEC=0
      DELLIE=0
      CALL MTLPRT('DELIEC','C347.0','Not available on this machine')
      RETURN
      END
#endif
