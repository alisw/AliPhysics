*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:28  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/tmpro.F"
#else
      SUBROUTINE TMPRO (TEXT)
C
C CERN PROGLIB#         TMPRO           .VERSION KERNFOR  4.37  930715
C ORIG. 30/05/91, JZ
C
C     Fortran interface : print a prompt string
C
      CHARACTER    TEXT*(*)

      LGTEXT = LEN (TEXT)
#if defined(CERNLIB_QMVAX)
      CALL TMPROI (%REF(TEXT),LGTEXT)
#endif
#if !defined(CERNLIB_QMVAX)
      CALL TMPROI (TEXT,LGTEXT)
#endif
      RETURN
      END
#endif
