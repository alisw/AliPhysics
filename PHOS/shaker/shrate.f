*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHRATE
c	=================

c	Compute rates for various particles.
c	!!! eta --> gamma gamma only !!!

*KEEP,SHPHYP.
      COMMON /SHPHYP/ JWEI,NDNDY,YLIM,PTLIM,JWEAK,JPI0,JETA,JPIC,JPRO,
     +                  JKAC,JKA0,JRHO,JOME,JPHI,JPSI,JDRY
*KEEP,SHRATS.
      COMMON /SHRATS/ RETAPI,RPROPI,RKACPI,RRHOPI,ROMEPI,RPHIPI
*KEEP,SHPRAT.
      COMMON /SHPRAT/ PI0R,ETAR,RHOR,OMER,PHIR,PSIR,DRYR
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEEP,SHPSDY.
      COMMON /SHPSDY/ RNPSI,RNDRY
*KEND.

      CHA = NDNDY*2*YLIM
      NCHA = CHA
      NPIC = JPIC*NCHA/(1+JPRO*RPROPI+JKAC*RKACPI+JKA0*RKACPI)
      NPRO = JPRO*RPROPI*NPIC
      NKAC = JKAC*RKACPI*NPIC
      NKA0 = JKA0*RKACPI*NPIC

      IF (JPIC.EQ.0) NPI0=JPI0*NCHA/2
      IF (JPIC.EQ.1) NPI0=JPI0*NPIC/2
      NETA=JETA*RETAPI*0.389*NPI0  !eta/pi0 * npi0 * BR[eta --> gamma gamma]

c	Initialize /SHPRAT/ rates...

c	... to half of charged multiplicity for pi0

      PI0R = FLOAT(NCHA)/2.

c	... to pi0 * eta/pi0 * BR --> gamma gamma for eta

      ETAR = PI0R*RETAPI*.389

c	... to pi0 * vect/pi0 * BR --> e+e- for vector mesons

      RHOR = PI0R*RRHOPI*4.44E-5	
      OMER = PI0R*ROMEPI*7.07E-5
      PHIR = PI0R*RPHIPI*3.11E-4

c	... to calculated rate for J/psi and Drell-Yan

      PSIR = RNPSI
      DRYR = RNDRY


      RETURN
      END
