*
* $Id$
*
* $Log$
* Revision 1.5  1998/02/09 15:59:47  japost
*   Fixed a problem on AIX 4 xlf, caused by max(double,float).
*
* Revision 1.4  1998/02/06 16:46:57  japost
* Fix a wrong parenthesis.
*
* Revision 1.3  1998/02/06 16:22:24  japost
*   Protected a square root from a negative argument.
*   This root was added there in previous changes, and not deleted from its
* old position. In its old position it was protected from being negative, but in
* its new position it was not.
*
*   Deleted the same square root from its old position, as it was redundant.
*
* Revision 1.2  1996/03/13 12:03:24  ravndal
* Tranverse momentum conservation
*
* Revision 1.1.1.1  1995/10/24 10:21:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 21/02/95  11.53.59  by  S.Giani
*-- Author :
#if defined(CERNLIB_HPUX)
$OPTIMIZE OFF
#endif
      SUBROUTINE GPAIRG
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Simulates e+e- pair production by photons.                    *
C.    *                                                                *
C.    *  The secondary electron energies are sampled using the         *
C.    *  Coulomb corrected BETHE-HEITLER cross-sections.For this the   *
C.    *   modified version of the random number techniques of          *
C.    *   BUTCHER and MESSEL (NUCL.PHYS,20(1960),15) are employed.     *
C.    *                                                                *
C.    *  NOTE :                                                        *
C.    *  (1) Effects due to the breakdown of the BORN approximation at *
C.    *      low energies are ignored.                                 *
C.    *  (2) The differential cross-section implicitly takes account   *
C.    *      of pair production in both nuclear and atomic electron    *
C.    *      fields. However, triplet production is not generated.     *
C.    *                                                                *
C.    *    ==>Called by : GTGAMA                                       *
C.    *       Authors    G.Patrick, L.Urban  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcking.inc"
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
      DIMENSION NTYPEL(2)
      DIMENSION RNDM(2)
      LOGICAL ROTATE
      PARAMETER (ONE=1,ONETHR=ONE/3,EMAS2=2*EMASS)
C.
C.    ------------------------------------------------------------------
C.
C             If not enough energy : no pair production
C
      EGAM   = VECT(7)
      IF (EGAM.LT.EMAS2) GO TO 999
C
      KCASE  = NAMEC(6)
      IF(IPAIR.NE.1) THEN
         ISTOP  = 2
         NGKINE = 0
         DESTEP = DESTEP + EGAM
         VECT(7)= 0.
         GEKIN  = 0.
         GETOT  = 0.
         GO TO 999
      ENDIF
C
C             For low energy photons approximate the electron energy by
C             sampling from a uniform distribution in the interval
C             EMASS -> EGAM/2.
C
      IF (EGAM.LE.2.1E - 03)THEN
         CALL GRNDM(RNDM,1)
         EEL1   = EMASS + (RNDM(1)*(0.5*EGAM - EMASS))
         X=EEL1/EGAM
         GO TO 20
      ENDIF
C
      Z3=Q(JPROB+2)
      F=8.*Q(JPROB+3)
      IF(EGAM.GT.0.05) F=F+8.*Q(JPROB+4)
      X0=EMASS/EGAM
      DX=0.5-X0
      DMIN=544.*X0/Z3
      DMIN2=DMIN*DMIN
      IF(DMIN.LE.1.)THEN
         F10=42.392-7.796*DMIN+1.961*DMIN2-F
         F20=41.405-5.828*DMIN+0.8945*DMIN2-F
      ELSE
         F10=42.24-8.368*LOG(DMIN+0.952)-F
         F20=F10
      ENDIF
C
C             Calculate limit for screening variable,DELTA, to ensure
C             that screening rejection functions always remain
C             positive.
C
      DMAX=EXP((42.24-F)/8.368)-0.952
C
C             Differential cross-section factors which form
C             the coefficients of the screening functions.
C
      DSIG1=DX*DX*F10/3.
      DSIG2=0.5*F20
      BPAR   = DSIG1 / (DSIG1 + DSIG2)
C
C             Decide which screening rejection function to use and
C             sample the electron/photon fractional energy BR.
C
   10 CALL GRNDM(RNDM,2)
      IF(RNDM(1).LT.BPAR)THEN
         X=0.5-DX*RNDM(2)**ONETHR
         IREJ=1
      ELSE
         X=X0+DX*RNDM(2)
         IREJ   = 2
      ENDIF
C
C             Calculate DELTA ensuring positivity.
C
      D=0.25*DMIN/(X*(1.-X))
      IF(D.GE.DMAX) GOTO 10
      D2=D*D
C
C             Calculate F1 and F2 functions using approximations.
C             F10 and F20 are the F1 and F2 functions calculated for the
C             DELTA=DELTA minimum.
C
      IF(D.LE.1.)THEN
         F1=42.392-7.796*D+1.961*D2-F
         F2=41.405-5.828*D+0.8945*D2-F
      ELSE
         F1=42.24-8.368*LOG(D+0.952)-F
         F2=F1
      ENDIF
      IF(IREJ.NE.2)THEN
         SCREJ=F1/F10
      ELSE
         SCREJ=F2/F20
      ENDIF
C
C             Accept or reject on basis of random variate.
C
      CALL GRNDM(RNDM,1)
      IF(RNDM(1).GT.SCREJ) GOTO 10
      EEL1=X*EGAM
C
C             Successful sampling of first electron energy.
C
C             Select charges randomly.
C
   20 NTYPEL(1) = 2
      CALL GRNDM(RNDM,2)
      IF (RNDM(1).GT.0.5) NTYPEL(1) = 3
      NTYPEL(2) = 5 - NTYPEL(1)
C
C             Generate electron decay angles with respect to a Z-axis
C             defined along the parent photon.
C             PHI is generated isotropically and THETA is assigned
C             a universal angular distribution
C
      EMASS1 = EMASS
      THETA  = GBTETH(EEL1, EMASS1, X)*EMASS/EEL1
      SINTH  = SIN(THETA)
      COSTH  = COS(THETA)
      PHI    = TWOPI*RNDM(2)
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
C
C             Rotate tracks into GEANT system
C
      CALL GFANG(VECT(4),COSAL,SINAL,COSBT,SINBT,ROTATE)
C
C            Polar co-ordinates to momentum components.
C
      NGKINE = 0
      TEL1 = EEL1 - EMASS
      PEL1 = SQRT(MAX((EEL1+REAL(EMASS))*TEL1,0.))
      IF(TEL1.GT.CUTELE) THEN
         NGKINE = NGKINE + 1
         GKIN(1,NGKINE) = PEL1 * SINTH * COSPHI
         GKIN(2,NGKINE) = PEL1 * SINTH * SINPHI
         GKIN(3,NGKINE) = PEL1 * COSTH
         GKIN(4,NGKINE) = EEL1
         GKIN(5,NGKINE) = NTYPEL(1)
         TOFD(NGKINE)=0.
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
         IF(ROTATE)
     +   CALL GDROT(GKIN(1,NGKINE),COSAL,SINAL,COSBT,SINBT)
      ELSE
         DESTEP = DESTEP + TEL1
         IF(NTYPEL(1).EQ.2) CALL GANNI2
      ENDIF
C
C             Momentum vector of second electron. Recoil momentum of
C             target nucleus/electron ignored.
C
      EEL2=EGAM-EEL1
      TEL2=EEL2-EMASS
      IF(TEL2.GT.CUTELE) THEN
         PEL2 = SQRT((EEL2+EMASS)*TEL2)
         NGKINE = NGKINE + 1
         SINTH=SINTH*PEL1/PEL2
         COSTH=SQRT(MAX(0.,1.-SINTH**2))
         GKIN(1,NGKINE)=-PEL2*SINTH*COSPHI
         GKIN(2,NGKINE)=-PEL2*SINTH*SINPHI
         GKIN(3,NGKINE)=PEL2*COSTH
         GKIN(4,NGKINE)=EEL2
         GKIN(5,NGKINE) = NTYPEL(2)
         TOFD(NGKINE)=0.
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
         IF(ROTATE)
     +   CALL GDROT(GKIN(1,NGKINE),COSAL,SINAL,COSBT,SINBT)
      ELSE
         DESTEP = DESTEP + TEL2
         IF(NTYPEL(2).EQ.2) CALL GANNI2
      ENDIF
      ISTOP = 1
      IF(NGKINE.EQ.0) ISTOP = 2
 999  END
#if defined(CERNLIB_HPUX)
$OPTIMIZE ON
#endif
