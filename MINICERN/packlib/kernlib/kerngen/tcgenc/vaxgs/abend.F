*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:39  mclareni
* Kernlib
*
*
      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNVAX  2.37  910718
C ORIG. 10/06/91 J.C.Marin
C
      INTEGER IVAL

      IVAL = '080F8FF4'X
      CALL LIB$STOP (IVAL)
      END
