*
* $Id$
*
* $Log$
* Revision 1.2.4.1  2002/11/26 16:50:55  hristov
* Merging NewIO with v3-09-04
*
* Revision 1.2  2002/10/14 14:57:11  hristov
* Merging the VirtualMC branch to the main development branch (HEAD)
*
* Revision 1.1.2.1  2002/07/11 17:14:49  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
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
