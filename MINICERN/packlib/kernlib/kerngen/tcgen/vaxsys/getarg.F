*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:59  mclareni
* Kernlib
*
*
      INTEGER FUNCTION GETARG (JARG,TEXT)

C CERN PROGLIB# Z265    GETARG          .VERSION KERNVAX  2.45  940228
C ORIG. 16/02/94 MM+JS+JZ

C-    Simulate UNIX system calls GETARG, IARGC

      INCLUDE '($JPIDEF)'

      CHARACTER       CHWK*1024
      COMMON /CHSLAT/ CHWK
      COMMON /SLATE/ NDSLAT, NESLAT, ISLATE(38)
      CHARACTER    TEXT*(*)
      INTEGER*2    LGR

      CHARACTER    COMD*256
      INTEGER      NARGS, NCOMD
      SAVE         NARGS, NCOMD, COMD
      DATA NARGS  / -1 /

C----              Doing GETARG

      GETARG = 0
      MODE = 0
      IF (NARGS.LT.0)              GO TO 71
   22 TEXT = ' '
      IF (JARG.EQ.0)               GO TO 41
      IF (JARG.GT.NARGS)           RETURN
      IF (JARG.LE.0)               RETURN

      JNX = 1
      JL  = 0
   24 JL  = JL + 1
      JGO = JNX
      JF  = ICNEXT (COMD,JGO,NCOMD)
      IF (COMD(JF:JF).EQ.'"')  THEN
          JF  = JF + 1
          JNX = ICFIND ('"', COMD,JF,NCOMD)
          N   = JNX - JF
          JNX = JNX + 1
        ELSE
          JNX = NESLAT
          N   = NDSLAT
        ENDIF
      IF (JL.LT.JARG)              GO TO 24
      IF (N.GT.0)  TEXT(1:N) = COMD(JF:JF+N-1)
      RETURN

C--           doing GETARG (0, TEXT), get image-name as for Unix

   41 NDSLAT = LIB$GETJPI (JPI$_IMAGNAME,,,,CHWK,%REF(LENSTR))
      IF (.NOT.NDSLAT)             RETURN
      ISTAT = LIB$TRIM_FILESPEC (CHWK(1:LENSTR),TEXT,LEN(TEXT),)
      NDSLAT = 0
      RETURN

C----              Doing IARGC

      ENTRY  IARGC()

      MODE = 1
      IF (NARGS.LT.0)              GO TO 71
   68 IARGC = NARGS
      RETURN

C----              Load the calling command

   71 COMD = ' '
      CALL LIB$GET_FOREIGN (COMD,,LGR)

      NARGS = 0
      IF (LGR.EQ.0)                GO TO 77

      NCOMD = LNBLNK (COMD(1:LGR))
      IF (NCOMD.EQ.0)              GO TO 77

      JGO   = 1
   74 NARGS = NARGS + 1
      JF  = ICNEXT (COMD,JGO,NCOMD)
      IF (COMD(JF:JF).EQ.'"')  THEN
          JF  = JF + 1
          JGO = ICFIND ('"', COMD,JF,NCOMD) + 1
        ELSE
          JGO = NESLAT
        ENDIF
      IF (JGO.LE.NCOMD)            GO TO 74

   77 IF (MODE.EQ.0)               GO TO 22
      GO TO 68
      END
