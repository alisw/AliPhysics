*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:40  mclareni
* Kernlib
*
*
      SUBROUTINE DATIMH(ND,NT)
C
C CERN PROGLIB# Z007    DATIMH          .VERSION KERNMPW  1.00  891208
C
      INTEGER NT(2), ND(2)
      CHARACTER*8 CTIME
      CHARACTER*6 CDATE
C
      CALL TIME(CTIME)
      CALL IDATE(IMONTH,IDAY,IYEAR)
      WRITE(CDATE,'(3I2.2)') IDAY, IMONTH, IYEAR
      NT(1) = 0
      NT(2) = 0
      ND(1) = 0
      ND(2) = 0
C
      NT(2) = ICHAR(CTIME(8:8))
      NT(2) = IOR(ISHFT(ICHAR(CTIME(7:7)),8),NT(2))
      NT(2) = IOR(ISHFT(ICHAR('.'),16),NT(2))
      NT(2) = IOR(ISHFT(ICHAR(CTIME(5:5)),24),NT(2))
      NT(1) = ICHAR(CTIME(4:4))
      NT(1) = IOR(ISHFT(ICHAR('.'),8),NT(1))
      NT(1) = IOR(ISHFT(ICHAR(CTIME(2:2)),16),NT(1))
      NT(1) = IOR(ISHFT(ICHAR(CTIME(1:1)),24),NT(1))
C
      ND(2) = ICHAR(CDATE(6:6))
      ND(2) = IOR(ISHFT(ICHAR(CDATE(5:5)),8),ND(2))
      ND(2) = IOR(ISHFT(ICHAR('/'),16),ND(2))
      ND(2) = IOR(ISHFT(ICHAR(CDATE(4:4)),24),ND(2))
      ND(1) = ICHAR(CDATE(3:3))
      ND(1) = IOR(ISHFT(ICHAR('/'),8),ND(1))
      ND(1) = IOR(ISHFT(ICHAR(CDATE(2:2)),16),ND(1))
      ND(1) = IOR(ISHFT(ICHAR(CDATE(1:1)),24),ND(1))
      END
