*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:57  mclareni
* Kernlib
*
*
      SUBROUTINE FFSIZE (LUN,NWREC,JREC,ISTAT)
C
C CERN PROGLIB#         FFSIZE          .VERSION KERNHPX  1.01  920304
C ORIG.  01/02/91  JZ
C
C-    Get the size of a disk file connected to Fortran unit LUN


      LUNPTR = FNUM (LUN)
      IF (LUNPTR.EQ.0)             GO TO 91

      CALL CFSIZE (LUNPTR,0,NWREC,JREC,ISTAT)
      RETURN

   91 ISTAT = -7
      JREC  = 0
      RETURN
      END
