*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHSDEC(JSEL)
c	=======================

c	To select allowed decay modes

*KEEP,LUDAT3.
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      SAVE /LUDAT3/
*KEND.

      IF (JSEL.EQ.0) RETURN

c	Define e- e+ for rho

      MDME(673,2) = 0
      KFDP(673,1) = 11
      KFDP(673,2) = -11
      KFDP(673,3) = 0
      KFDP(673,4) = 0
      KFDP(673,5) = 0

c	Define e- e+ for omega

      MDME(674,1) = 0
      MDME(676,1) = 0
      MDME(677,1) = 0
      MDME(678,1) = 0
      MDME(675,2) = 0
      KFDP(675,1) = 11
      KFDP(675,2) = -11
      KFDP(675,3) = 0
      KFDP(675,4) = 0
      KFDP(675,5) = 0

c	Define e- e+ for phi

      MDME(679,1) = 0
      MDME(680,1) = 0
      MDME(682,1) = 0
      MDME(683,1) = 0
      MDME(684,1) = 0
      MDME(685,1) = 0
      MDME(686,1) = 0
      MDME(681,2) = 0
      KFDP(681,1) = 11
      KFDP(681,2) = -11
      KFDP(681,3) = 0
      KFDP(681,4) = 0
      KFDP(681,5) = 0

c	Select e- e+ for J/psi

      MDME(688,1) = 0		
      MDME(689,1) = 0		

c	GO TO 1700

c	Select gamma gamma for pi0 and eta

      MDME(639,1) = 0

      MDME(641,1) = 0
      MDME(642,1) = 0
      MDME(643,1) = 0
      MDME(644,1) = 0
      MDME(645,1) = 0

      GO TO 2000

c	Select Dalitz for pi0 and eta

1500	CONTINUE

      MDME(638,1) = 0		
      MDME(639,1) = 1

      MDME(640,1) = 0
      MDME(641,1) = 0
      MDME(642,1) = 0
      MDME(643,1) = 0
      MDME(644,1) = 1
      MDME(645,1) = 0

1700	CONTINUE

c	Select full BR for pi0

      MDME(638,1) = 1		
      MDME(639,1) = 1



2000	CONTINUE
      RETURN
      END
