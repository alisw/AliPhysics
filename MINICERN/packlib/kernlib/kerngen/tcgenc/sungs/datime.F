*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:39  mclareni
* Kernlib
*
*
      SUBROUTINE DATIME (ND,NT)
C
C CERN PROGLIB# Z007    DATIME          .VERSION KERNSUN  1.02  891025
C ORIG. 22/09/88, JZ
C MOD.  27/07/89  add /SLATE/, JZ
C
C-    CALL DATIME (ND,NT)   RETURNS INTEGER DATE   ND = YYMMDD
C-                                  INTEGER TIME   NT =   HHMM
C-    returns also ISL  (1) y  (2) m  (3) d  (4) h  (5) m  (6) s
C
      COMMON /SLATE/ ISL(34), MM(6)


      CALL IDATE (MM)
      CALL ITIME (MM(4))

      DO 24  J=1,6
   24 ISL(J) = MM(J)
      ISL(1) = MM(3)
      ISL(3) = MM(1)

      MM(3) = MOD (MM(3),100)

      ND = (MM(3)*100+MM(2))*100 + MM(1)
      NT =            MM(4) *100 + MM(5)
      RETURN
      END
