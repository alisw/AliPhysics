*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
c
c	Federico Antinori Productions is proud to present...
c
c                 S H A K E R
c	
c	Central Rapidity Phase Space Cocktail Event Generator
c
c	=====================================================

c	version 0/05

c	09.12.91, FA: 0/03	Pre-release
c	28.02.92, FA: 0/04	SH2BOD and SHMTSC


      SUBROUTINE SH2BOD(AMP,AMDA,AMDB,PP,PDA,PDB)
c	===========================================

c	Two-body decay of a parent of mass AMP and lab 4-momentum PP to
c	a daughter of mass AMDA and one of mass AMDB.
c	PDA and PDB contain the lab 4-momenta of the daughters

      DIMENSION PP(4),PDA(4),PDB(4)
      DIMENSION PACM(4),PBCM(4)

      IF ((AMDA+AMDB).GT.AMP) RETURN
      PACM(4) = (AMP**2+AMDA**2-AMDB**2)/2./AMP
      PBCM(4) = (AMP**2+AMDB**2-AMDA**2)/2./AMP
      PCM = SQRT(PACM(4)**2-AMDA**2)
      CT = 2.*RLU(0.)-1.
      ST = SQRT(1.-CT**2)
      PHI = 2.*3.14159*RLU(0.)
      PACM(1) = PCM*ST*COS(PHI)
      PACM(2) = PCM*ST*SIN(PHI)
      PACM(4) = PCM*CT
      PBCM(1) = -PACM(1)
      PBCM(2) = -PACM(2)
      PBCM(3) = -PACM(3)
      CALL LORENB(AMP,PP,PACM,PDA)
      CALL LORENB(AMP,PP,PBCM,PDB)

      RETURN
      END
