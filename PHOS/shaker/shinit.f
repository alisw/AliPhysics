*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHINIT
c	=================

c	Program Initialization
c	eta/pi0, rho/pi0, omega/pi0, phi/pi0 from asymptotic value
c       & mT-scaling
c	[E.g.: V.Hedberg, LUNFD6/(NFFL-7073)/1987;
c       M.Bourquin and J.M.Gaillard, Nucl. Phys. B114 (1976),334]
c	p/pi and K/pi from Tevatron
c	[T.Alexopoulos et al.: Phys. Rev. Lett. 64 (1990), 991]
c	psi and Drell-Yan are rates, independent of dN/dy, for one unit
c       of rapidity, for central events.
c       Cross section data from J. Schukraft, priv.comm.
c	Drell-Yan rate is for M > 1 GeV
c	psi and Drell-Yan rates include BR --> e+ e-

*KEEP,SHRUNP.
      COMMON /SHRUNP/ VMAJ,IMIN,NRUN,NEVTOT
*KEEP,SHPHYP.
      COMMON /SHPHYP/ JWEI,NDNDY,YLIM,PTLIM,JWEAK,JPI0,JETA,JPIC,JPRO,
     +                  JKAC,JKA0,JRHO,JOME,JPHI,JPSI,JDRY
*KEEP,SHRATS.
      COMMON /SHRATS/ RETAPI,RPROPI,RKACPI,RRHOPI,ROMEPI,RPHIPI
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEEP,SHPSDY.
      COMMON /SHPSDY/ RNPSI,RNDRY
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.

c	LUND initialization

      MSTU(11) = 3			! output file
      MSTU(4)  = 200000		! /LUJETS/ size

c	SHAKER initialization

      VMAJ = 0.0	! Major version number
      IMIN = 5	! Minor version number

c	particle ratios initialization

      RETAPI = 0.17	! eta over pi0 ratio
      RPROPI = 0.074	! p over pi ratio
      RKACPI = 0.112	! K over pi ratio
      RRHOPI = 0.15	! rho over pi0 ratio
      ROMEPI = 0.14	! omega over pi0 ratio
      RPHIPI = 0.016	! phi over pi0 ratio

      RNPSI  = 0.0018 ! J/psi rate * BR
      RNDRY  = 0.0005 ! Drell-Yan rate * BR for M > 1 GeV



      CALL SHRATE	
      IF (JWEAK.EQ.0) CALL SHWDIS
      CALL SHRNDV
      CALL SHSDEC(1)	! Select decays?

      RETURN
      END
