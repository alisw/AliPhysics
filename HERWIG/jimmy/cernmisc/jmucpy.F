*
* $Id$
*
* $Log$
* Revision 1.1  2005/07/29 16:15:46  jmb
* Include the various CERNLIB functions jimmy needs
*
* Revision 1.1.1.1  1996/02/15 17:50:19  mclareni
* Kernlib
*
*
      SUBROUTINE JMUCPY (A,B,N)
C
C CERN PROGLIB# V301    UCOPY           .VERSION KERNFOR  4.40  940929
C ORIG. 01/01/65 JZ
C
      REAL*8  A(*),B(*)
C
C--                NO OVERLAP  OR  BEGINNING OF A ON END OF B
      IF (N.EQ.0) RETURN
         DO 21 I=1,N
   21 B(I)=A(I)
      RETURN
      END
