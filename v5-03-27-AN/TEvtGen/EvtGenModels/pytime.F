 
C*********************************************************************
 
C...PYTIME
C...Finds current date and time.
C...Since this task is not standardized in Fortran 77, the routine
C...is dummy, to be replaced by the user. Examples are given for
C...the Fortran 90 routine DEC Fortran 77, GNU Linux g77, and what to 
C...do if you do not have access to suitable routines.
 
      SUBROUTINE PYTIME(IDATI)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C      INTEGER PYK,PYCHGE,PYCOMP
C      CHARACTER*8 ATIME
C...Local array.
      INTEGER IDATI(6),IDTEMP(3)
 
C...Example 0: if you do not have suitable routines.
      DO 100 J=1,6
      IDATI(J)=0
  100 CONTINUE
 
C...Example 1: Fortran 90 routine.
C      INTEGER IVAL(8)
C      CALL DATE_AND_TIME(VALUES=IVAL)
C      IDATI(1)=IVAL(1)
C      IDATI(2)=IVAL(2)
C      IDATI(3)=IVAL(3)
C      IDATI(4)=IVAL(5)
C      IDATI(5)=IVAL(6)
C      IDATI(6)=IVAL(7)
 
C...Example 2: DEC Fortran 77. AIX.
C      CALL IDATE(IMON,IDAY,IYEAR)
C      IDATI(1)=IYEAR
C      IDATI(2)=IMON
C      IDATI(3)=IDAY
C      CALL ITIME(IHOUR,IMIN,ISEC)
C      IDATI(4)=IHOUR
C      IDATI(5)=IMIN
C      IDATI(6)=ISEC
 
C...Example 3: DEC Fortran, IRIX, IRIX64.
C      CALL IDATE(IMON,IDAY,IYEAR)
C      IDATI(1)=IYEAR
C      IDATI(2)=IMON
C      IDATI(3)=IDAY
C      CALL TIME(ATIME)
C      IHOUR=0
C      IMIN=0
C      ISEC=0
C      READ(ATIME(1:2),'(I2)') IHOUR
C      READ(ATIME(4:5),'(I2)') IMIN
C      READ(ATIME(7:8),'(I2)') ISEC
C      IDATI(4)=IHOUR
C      IDATI(5)=IMIN
C      IDATI(6)=ISEC
 
C...Example 4: GNU LINUX libU77, SunOS.
      CALL IDATE(IDTEMP)
      IDATI(1)=IDTEMP(3)
      IDATI(2)=IDTEMP(2)
      IDATI(3)=IDTEMP(1)
      CALL ITIME(IDTEMP)
      IDATI(4)=IDTEMP(1)
      IDATI(5)=IDTEMP(2)
      IDATI(6)=IDTEMP(3)
 
C...Common code to ensure right century.
      IDATI(1)=2000+MOD(IDATI(1),100)
 
      RETURN
      END
