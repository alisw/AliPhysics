*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:37  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_ASHO)
*CMZ :  3.21/02 29/03/94  15.41.24  by  S.Giani
*-- Author :
      SUBROUTINE GASHO(P,XXMASS,GSTEP,DE)
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcasho.inc"
C-----------------------------------------------------------------------
C      ASHO model for energy loss straggling
C      The main subprograms are:
C      - MIXMGO - prepares the parameters of a mixture
C          or compound used (called from GIASHO);
C      - GOSCMG - prepares the parameters of the oscillators;
C      - GDIFMG - calculates the differential spectrum of the
C          energy loss;
C      The main initial data bank is in COMMON block GCASHO contai-
C      ning the parameters of a number of elements.
C-----------------------------------------------------------------------
C      Update history: 16-8-93
C
C-----------------------------------------------------------------------
C      The main initial common block and its components.
C      NELM - number of basic elements;
C      ZELM - atomic numbers of the elements;
C      AELM - atomic weights of the elements;
C      DELM - densities of the elements (g/cm3);
C      E0ELM - ionization potentials of the elements (keV);
C      NSELM - numbers of the electron shells in the atoms;
C      ZSELM - numbers of electrons in each of the shells;
C      ESELM - binding energies of the shell electrons (keV).
C       This information may be included in JMATE by GSMATE or GMATE.
C      *)   The density values for gases are given at
C           0 C 1 atm. Those for solids are not precise and
C           should be substituted for particular
C           samples used.
C      **)  The current number of the elements used is equal to 16.
C      ***) To change the energy scale from keV to GeV, as it is
C           used in GEANT, it is needed to add 1.E-6 to E0ELM, ESELM,
C           PLASM, TRNSMA, 6. (in RESMGO), SL (in REAMGO), SL (in
C           SLDMGO) and extract 1.E-6 in DE (in DELMGO).
C-----------------------------------------------------------------------
C      GCASHO contains also the kinematical parame-
C      ters as well as some thermodynamical ones:
C      PLIN (P/mc), PLOG (log10(PLIN)), BE2 (velocity/c squared),
C      the path length STEP of the particle,
C      and the medium plasma energy
C      PLASM. TRNSMA is the maximum transferable energy (now
C      it is set to correspond to TCUT in GEANT !).
C      In GEANT the initial parameters are P, E, XMASS, STEP
C      and DENSIT.
C-----------------------------------------------------------------------
C      In addition, it contains the parameters of the oscillators:
C      - BOSC is the array of "ksi/I",
C      - AOSC is the array of "lnA",
C      - NOSC is the number of the oscillators,
C      - EOSC is the array of the oscillator energies,
C      - IOSC is the upper array in integers,
C      - ZOSC is the array of the oscillator weights,
C      - EMEAN is the mean energy loss.
C-----------------------------------------------------------------------
C      Finally, GCASHO contains the parameters of the energy loss
C      distribution:
C      - CMGO is the differential or integral distribution coeff.,
C      - NMGO is the number of bins,
C      - NMGOMA is the maximum of possible NMGO,
C      - EMGO is the energy corresponding to one bin,
C      - EMGOMI is the energy corresponding to 0th bin.
C-----------------------------------------------------------------------
      PARAMETER (DNMGOM=2000)
C-----------------------------------------------------------------------
      JASHO = LQ(JMA-20)
      NSMED = Q(JASHO+1)
      ZMED  = Q(JASHO+2)
      AMED  = Q(JASHO+3)
      ALFA  = Q(JASHO+4)
      E0MED = Q(JASHO+5)
      DO 10 KMED=1,NSMED
         ZSMED(KMED) = Q(JASHO+5+KMED)
         ESMED(KMED) = Q(JASHO+5+NSMED+KMED)
   10 CONTINUE
      PLASM = 0.028817*SQRT((ZMED/AMED)*DENS)
      STEP  = GSTEP
      NMGOMA = DNMGOM
C-----------------------------------------------------------------------
C      In GEANT the initial kinematic parameters : P, E, XMASS.
C      The following gets the kinematic parameters of the particle:
C      PLIN (P/mc), PLOG (lg(PLIN)), BE2 (velocity/c squared) and
C      its path length STEP.
C-----------------------------------------------------------------------
      PLIN = P/XXMASS
      PLOG = LOG10(PLIN)
      BE2 = PLIN**2/(1.+PLIN**2)
C-----------------------------------------------------------------------
C      Now we set TRNSMA to correspond to TCUT in GEANT.
C-----------------------------------------------------------------------
      TRNSMA = DCUTE*1E6
C-----------------------------------------------------------------------
C      Here is the calculation of the parameters of the oscillators
C-----------------------------------------------------------------------
      CALL GOSCMG
C-----------------------------------------------------------------------
C      Here is the calculation of the diffrential distribution
C-----------------------------------------------------------------------
      CALL GDIFMG(DE)
      DE = 1E-6*DE
C-----------------------------------------------------------------------
      END
#endif
