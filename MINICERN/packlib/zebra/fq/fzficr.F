*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZFICR (MODE,CHOPT)

C-    Crack (MODE=0) or pack the options for FZFILE

#include "zebra/fzcf.inc"
#include "fzficc.inc"
C--------------    End CDE                             --------------
      CHARACTER    CHOPT*(*)
      PARAMETER   (NOPT=23, NOPCH=14)
      DIMENSION    IOPCH(NOPCH)
      EQUIVALENCE (IOPCH(1),IOPTVF(5))
      CHARACTER    OPTSTR*(NOPT), CHX(NOPCH-1)*10, CHWK*1

      DATA  OPTSTR / '0123MCADT*KLYFXNSUIORPQ' /
C-                    _:.=+=.:_1_:.=+=.:_2_:.

C-                        _:.=+=.:_1_:.=
C-                        MCADT*KLY
      DATA CHX    /     'C-        '
     +,                 'A--       '
     +,                 'D-+-      '
     +,                 'T-- -     '
     +,                 '*-----    '
     +,                 'K---??+   '
     +,                 'L---+++-  '
     +,                 'Y----++-- '
     +,                 'F--iddd---'
     +,                 'Xiiii++iii'
     +,                 'N+++++d+++'
     +,                 'Si++++++++'
     +,                 'U++-------' /
C-                        MCADT*KLY


C--                Blow the option string

      IF (MODE.NE.0)               GO TO 41
      CALL UOPTC (CHOPT,OPTSTR,IOPTVF)
      IFLERR = 0

C--                Check compatibility

      DO  27  JP=1,9
      IF (IOPCH(JP).EQ.0)         GO TO 27

      DO  24  JS=JP+1,NOPCH
      CHWK = CHX(JS-1)(JP+1:JP+1)
      IF (IOPCH(JS).EQ.0)  THEN
          IF (CHWK.EQ.'i')  IOPCH(JS) = -1
          GO TO 24
        ENDIF

      IF (CHWK.EQ.'-')  THEN
          IOPCH(JS) = 0
          IFLERR = IFLERR + 1
        ENDIF
   24 CONTINUE
   27 CONTINUE

      IF (IOPTI+IOPTO.EQ.0)  IOPTI=1

      NEOF = -1
      IF (IOPTVF(4) .NE.0)  NEOF=3
      IF (IOPTVF(3) .NE.0)  NEOF=2
      IF (IOPTVF(2) .NE.0)  NEOF=1
      IF (IOPTVF(1) .NE.0)  NEOF=0
      IOPTVF(1) = NEOF

      IF (IOPTM+IOPTC+IOPTA+IOPTD+IOPTT+IOPTAS.EQ.0)  IOPTAS=-1
      RETURN

C----              Construct text string of options used

   41 CHOPT  = ' '
      N = 0
      DO  44  J=5,NOPT
      IF (IOPTVF(J).EQ.1)  THEN
          N = N + 1
          CHOPT(N:N) = OPTSTR(J:J)
        ENDIF
   44 IOPTVF(J) = IABS(IOPTVF(J))

      NEOF = IOPTVF(1)
      IF (NEOF.GE.0)  THEN
          N = N + 2
          CALL CSETDI (NEOF,CHOPT,N,N)
        ENDIF
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
