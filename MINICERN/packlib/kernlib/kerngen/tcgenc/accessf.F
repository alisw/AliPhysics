*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      LOGICAL FUNCTION ACCESSF (FNAME, MODE)
C
C CERN PROGLIB# Z265    ACCESSF         .VERSION KERNFOR  4.39  940228
C ORIG. 17/02/94, JZ
C
      CHARACTER    FNAME*(*)
      LOGICAL      THERE

C--           Does the file exist?

      INQUIRE (FILE=FNAME,EXIST=THERE)
      IF (.NOT.THERE)              GO TO 39

      MRW = JBYT (MODE,1,2)
      IF (MRW.LT.2)                GO TO 49

C--           Find a free LUN

      DO  17 LUN=8,49
      INQUIRE (UNIT=LUN,OPENED=THERE)
      IF (.NOT.THERE)              GO TO 21
   17 CONTINUE
      GO TO 39

C--           Try to open the file

C-        R_OK    4   test for read permission
C-        W_OK    2   test for write permission
C-        X_OK    1   test for execute (search) permission
C-        F_OK    0   test for presence of file


   21 IF (MRW.LT.4)                GO TO 31
      OPEN (LUN,FILE=FNAME,STATUS='OLD'
     +,        FORM='FORMATTED', IOSTAT=ISTAT)
      IF (ISTAT.EQ.0)              GO TO 31

      OPEN (LUN,FILE=FNAME,STATUS='OLD'
     +,        FORM='UNFORMATTED', IOSTAT=ISTAT)
      IF (ISTAT.NE.0)              GO TO 38

   31 IF (MRW.EQ.4)                GO TO 48
      OPEN (LUN,FILE=FNAME,STATUS='UNKNOWN'
     +,        FORM='FORMATTED', IOSTAT=ISTAT)
      IF (ISTAT.EQ.0)              GO TO 48

      OPEN (LUN,FILE=FNAME,STATUS='UNKNOWN'
     +,        FORM='UNFORMATTED', IOSTAT=ISTAT)
      IF (ISTAT.EQ.0)              GO TO 48

   38 CLOSE (LUN)
   39 ACCESSF = .FALSE.
      RETURN

   48 CLOSE (LUN)
   49 ACCESSF = .TRUE.
      RETURN
      END
