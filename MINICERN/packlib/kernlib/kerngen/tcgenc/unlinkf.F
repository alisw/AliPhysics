*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/unlinkf.F"
#else
      SUBROUTINE UNLINKF (CHNAME)
C
C CERN PROGLIB# Z265    UNLINKF         .VERSION KERNFOR  4.34  930114
C ORIG. 20/10/92, JZ
C
      CHARACTER    CHNAME*(*)
      LOGICAL      THERE

C--           Does the file exist?

      INQUIRE (FILE=CHNAME,EXIST=THERE)
      IF (.NOT.THERE)              RETURN

C--           Find a free LUN

      DO  17 LUN=8,49
      INQUIRE (UNIT=LUN,OPENED=THERE)
      IF (.NOT.THERE)              GO TO 21
   17 CONTINUE
      RETURN

C--           Delete the file

   21 OPEN (LUN,FILE=CHNAME,STATUS='OLD'
     +,        FORM='FORMATTED', IOSTAT=ISTAT)
      IF (ISTAT.EQ.0)              GO TO 49

      OPEN (LUN,FILE=CHNAME,STATUS='OLD'
     +,        FORM='UNFORMATTED', IOSTAT=ISTAT)
      IF (ISTAT.EQ.0)              GO TO 49
      CLOSE (LUN)
      RETURN

   49 CLOSE (LUN, STATUS='DELETE')
      RETURN
      END
#endif
