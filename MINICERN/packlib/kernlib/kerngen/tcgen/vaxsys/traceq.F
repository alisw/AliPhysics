*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:59  mclareni
* Kernlib
*
*
      SUBROUTINE TRACEQ (LUN,LEV)
C
C CERN PROGLIB# N105    TRACEQ          .VERSION KERNVAX  1.09  810226
C ORIG. 26/02/81, J. BETTELS, CERN
C
      INTEGER*4 SYS$CRELOG
      CHARACTER ERRFIL*10

      ERRFIL = 'FOR000.DAT'
      J  = 6
      ND = LUN
      DO 9  I=1,3
         NN = ND/10
         NR = ND - 10*NN
         IF (NR.NE.0)  ERRFIL(J:J)=CHAR(NR+ICHAR('0'))
         ND = NN
         J  = J - 1
    9 CONTINUE
      II = SYS$CRELOG (%VAL (2),'SYS$ERROR',ERRFIL,)
      CALL LIB$SIGNAL (%VAL (0))
      RETURN
      END
