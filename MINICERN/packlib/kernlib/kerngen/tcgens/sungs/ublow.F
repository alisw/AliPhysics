*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:22  mclareni
* Kernlib
*
*
      SUBROUTINE UBLOW (MS,MT,NCHP)
C
C CERN PROGLIB# M409    UBLOW           .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C

      DIMENSION    MS(99), MT(99), NCHP(9)
C-  - PARAMETER    (IBLAN1 = X'00202020')
      PARAMETER    (IBLAN1 =    2105376 )
C-  - PARAMETER    (MASK1  = X'FF000000')
      PARAMETER    (MASK1  =  -16777216 )


      NCH = NCHP(1)
      IF   (NCH)             91, 29, 11
   11 NWS    = RSHIFT (NCH,2)
      NTRAIL = AND (NCH,3)
      JT     = 0
      IF (NWS.EQ.0)          GO TO 26

C--                Unpack the initial complete words

      DO 24 JS=1,NWS
      MWD      = MS(JS)
      MT(JT+1) = OR (IBLAN1,AND(MASK1,MWD))
      MT(JT+2) = OR (IBLAN1,AND(MASK1,LSHIFT(MWD, 8)))
      MT(JT+3) = OR (IBLAN1,AND(MASK1,LSHIFT(MWD,16)))
      MT(JT+4) = OR (IBLAN1,          LSHIFT(MWD,24) )
   24 JT = JT + 4

      IF (NTRAIL.EQ.0)       RETURN

C--                Unpack the trailing word

   26 MWD = MS(NWS+1)

      DO 28 JS=1,NTRAIL
      MT(JT+1) = OR (IBLAN1,AND(MASK1,MWD))
      MWD = LSHIFT (MWD,8)
   28 JT = JT + 1
   29 RETURN

   91 CALL ABEND
      END
