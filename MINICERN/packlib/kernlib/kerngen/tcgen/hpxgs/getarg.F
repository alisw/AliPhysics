*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:57  mclareni
* Kernlib
*
*
      SUBROUTINE GETARG (JARG, CHARG)
C
C CERN PROGLIB#         GETARG          .VERSION KERNHPX  1.00  910201
C ORIG.  13/03/91  JZ
C
C-    interface to HPX routine IGETARG

      CHARACTER    CHARG*(*)


      NCHMX = LEN (CHARG)
      NCH   = IGETARG (JARG, CHARG, NCHMX)

      IF (NCH.LT.NCHMX)  CHARG(NCH+1:) = ' '
      RETURN
      END
