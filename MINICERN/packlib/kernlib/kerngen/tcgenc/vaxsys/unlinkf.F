*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:41  mclareni
* Kernlib
*
*
      SUBROUTINE UNLINKF (CHNAME)
C
C CERN PROGLIB# Z265    UNLINKF         .VERSION KERNVAX  2.45  940228
C ORIG. 20/10/92, JZ, mod. J. Shiers, 15/02/94
C
      CHARACTER    CHNAME*(*)


      ISTAT = LIB$DELETE_FILE(CHNAME)
      RETURN
      END
