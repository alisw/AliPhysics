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
C CERN PROGLIB# Z007    DATIME          .VERSION KERNVAX  2.32  890926
C ORIG. 9/11/79 JZ
C MOD.  27/07/89  add /SLATE/, JZ
C

      COMMON /SLATE/ JY,JM,JD, JHR,JMN,JSC, ITM(4),M(30)

C-        1 2 3 4 5 6 7 8
C-        H H : M M : S S

      CALL TIME  (ITM)
      CALL UBLOW (ITM,M,8)

      DO 27 J=1,8
   27 M(J) = (M(J).AND.255) - 48
      JHR  = M(1)*10 + M(2)
      JMN  = M(4)*10 + M(5)
      JSC  = M(7)*10 + M(8)
      NT   = 100*JHR + JMN

      CALL IDATE (JM,JD,JY)
      ND = (JY*100+JM)*100+JD

      JY = JY + 1900
      RETURN
      END
