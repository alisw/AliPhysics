*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:38  mclareni
* Kernlib
*
*
      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C

      CALL EXIT (127)
      STOP
      END
