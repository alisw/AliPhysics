*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:40  mclareni
* Kernlib
*
*
      SUBROUTINE DATIMH (JDATE,JTIME)
C
C CERN PROGLIB# Z007    DATIMH          .VERSION KERNVAX  1.15  820721
C ORIG. 20/07/82 JZ
C
      DIMENSION JDATE(2), JTIME(2)

      COMMON /SLATE/ RES(10), JTM, JY,JM,JD, M(26)

      DATA IDOT /4H.   /

      CALL IDATE (JM,JD,JY)
      CALL UBLOW (8H00/00/00,M,8)
      CALL USET  (JD,M,1,2)
      CALL USET  (JM,M,4,5)
      CALL USET  (JY,M,7,8)
      CALL UBUNCH (M,JDATE,8)

      CALL TIME  (JTM)
      CALL UBLOW (JTM,M,8)
      M(3) = IDOT
      M(6) = IDOT
      CALL UBUNCH (M,JTIME,8)
      RETURN
      END
