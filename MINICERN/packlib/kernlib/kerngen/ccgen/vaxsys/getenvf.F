*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
* Kernlib
*
*
      SUBROUTINE GETENVF (CHSYM,CHTEXT)

C CERN PROGLIB# Z265    GETENVF         .VERSION KERNVAX  2.36  910306
C ORIG. 31/07/90, JZ
C MOD.  11/01/91, JZ, change name to GETENVF

C-    Simulate UNIX system call

      COMMON /SLATE/ NTEXT, ISTAT, ISLATE(38)
      CHARACTER    CHSYM*(*), CHTEXT*(*)

      NALL = 0
      NSYM = LNBLNK (CHSYM)
      IF (NSYM.LE.0)  THEN
          ISTAT = -7
          NALL  =  0
          GO TO 27
        ENDIF

      IF (NSYM.NE.4)               GO TO 21
      IF (CHSYM(1:4).NE.'HOME')    GO TO 21
      ISTAT = LIB$SYS_TRNLOG ('SYS$LOGIN',NALL,CHTEXT,,,%VAL(0))
      GO TO 24

   21 ISTAT = LIB$GET_SYMBOL (CHSYM(1:NSYM),CHTEXT,NALL)
   24 IF (ISTAT)  THEN
          ISTAT = 0
        ELSE
          NALL = 0
        ENDIF

   27 CHTEXT(NALL+1:) = ' '
      NTEXT = NALL
      RETURN
      END
