*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHWDIS
c	=================

c	Disable weak decays

*KEEP,LUDAT3.
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      SAVE /LUDAT3/
*KEND.

      MDCY(LUCOMP(310),1) = 0		! Disable K0s decay.

      RETURN
      END
