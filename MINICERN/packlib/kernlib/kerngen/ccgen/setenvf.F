*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:26  mclareni
* Kernlib
*
*
*FCA :          19/03/99  17:54:30  by  Federico Carminati
*               include the correct pilot.h
#include "../kerngen/pilot.h"
      INTEGER FUNCTION SETENVF (NAME, TEXT)
C
C CERN PROGLIB# Z265    SETENVF         .VERSION KERNFOR  4.36  930602
C ORIG. 22/02/91, JZ
C Fortran interface routine to setenv (BSD) or putenv (SYS V)
C
      CHARACTER*(*)  NAME, TEXT
      INTEGER        SETENI

      LGNM = LNBLNK (NAME)
      LGTX = LNBLNK (TEXT)
      LGTX = MAX (1, LGTX)
      IF (LGNM.EQ.0)               GO TO 91

      SETENVF = SETENI (NAME,TEXT,LGNM,LGTX)
      RETURN

   91 SETENVF = -1
      END
