*
* $Id$
*
* $Log$
* Revision 1.3  1997/06/20 18:32:55  japost
* A summary of the problem buf fix and its original report:
*
* ------------------------------------------------------------------------------
* From: Shawn McKee - (313) 764-4395 <mckee@pooh.physics.lsa.umich.edu>
* Subject: RE: Problem with FLUKA interface
* Date: Thu, 12 Jun 97 18:09:31 METDST
*
*
* I have seen this error with FLUKA.  I traced it to
* the following statement:
*
*      IF (RNDEVT .GE. SINE/FSIG) THEN
*
* Apparently FSIG (for this event) was not initialized
* properly (it is 0.0).  My fix was to insert a check
* on FSIG in the FLUFIN routine source.
*
* ------------------------------------------------------------------------------
* From: pniemine@estwm0.wm.estec.esa.nl
* Subject: Problem with FLUKA interface - cern.heplib #5667
* Date: Thu, 12 Jun 97 17:21:14 METDST
*
* In trying to run a GEANT program with FLUKA interface in our VMS-AXP (Alpha)
* cluster, I get the following error message (originating from FLUFIN) once
* a hadronic interaction is to be simulated:
*
* arithmetic trap, floating/decimal divide by zero
* ------------------------------------------------------------------------------
*
* Revision 1.2  1996/02/27 14:12:57  ravndal
* Correction for overstopped antiprotons
*
* Revision 1.1.1.1  1995/10/24 10:19:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 19/05/94  13.35.12  by  S.Ravndal
*FCA :          05/01/99  09:47:25  by  Federico Carminati
*               Added Xi- to the list of particles that are forced
*               to decay
*-- Author :
      SUBROUTINE FLUFIN
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcphys.inc"
#include "geant321/gctrak.inc"
#include "geant321/gsecti.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcunit.inc"
#include "geant321/dimpar.inc"
 
#if !defined(CERNLIB_SINGLE)
#include "geant321/finuct.inc"
#endif
#include "geant321/finuc.inc"
      REAL RNDM(1)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION AOCMBM, AMSS , ZTAR, RHO , ZLIN, ZLEL, ZLRAD,
     +ZUL
#endif
      COMMON / FKMAPA / AOCMBM (MXXMDF), AMSS (MXXMDF), ZTAR   (MXXMDF),
     +                  RHO    (MXXMDF), ZLIN (MXXMDF), ZLEL   (MXXMDF),
     +                  ZLRAD  (MXXMDF), ZUL  (MXXMDF), MEDIUM (MXXRGN),
     +                  MULFLG (MXXMDF),IFCOMP(MXXMDF), MSSNUM (MXXMDF),
     +                  NREGS, NMATF, MTBSNM
#if !defined(CERNLIB_SINGLE)
#include "geant321/part2t.inc"
#endif
#include "geant321/part2.inc"
#if !defined(CERNLIB_SINGLE)
#include "geant321/comcont.inc"
#endif
#include "geant321/comcon.inc"
#if !defined(CERNLIB_SINGLE)
#include "geant321/fheavyt.inc"
#endif
#include "geant321/fheavy.inc"
#include "geant321/paprop.inc"
#if !defined(CERNLIB_SINGLE)
#include "geant321/papropt.inc"
#endif
#include "geant321/gfkdis.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION POO,EKE,TXI,TYI,TZI,AMM,WE,ONE,PGEANT,DMOD
#endif
      PARAMETER (ONE=1)
      DIMENSION IGTOFL(49),IFLTOG(39),IHVTOG(6),ZSAMP(50)
      DATA IGTOFL / 0, 0, 0, 0, 0, 0,23,13,14,12, 15,16, 8, 1, 2,19, 0,
     +17,21,22, 20, 34, 36, 38, 9,18, 31, 32, 33, 35, 37, 39, 17*0/
 
      DATA IFLTOG /14,15, 3, 2, 4, 4, 1,13,25, 5, 6,10, 8, 9,11,12,18,
     +26,16,21, 19,20, 7, 7*0, 27, 28, 29, 22, 30, 23, 31, 24, 32/
      DATA IHVTOG /13,14,45,46,49,47/
*
      NP = 0
      NPHEAV = 0
*
*    Stopped particles:
*    o also xi- are forced to decay
*    o Neutral particles are sent to GHSTOP
*    o pi+ and K+/K- are forced to decay
*    o pi-, antiprotons and antineutrons are sent to FLUKA
*      for annihilation (not here but later in this routine)
      IF ((IGF.EQ.1).OR.
     +      (GEKIN.EQ.0..AND.ITRTYP.EQ.3.AND.IPART.NE.25)) THEN
         CALL GHEISH
         IGF = 0
         GOTO 999
      ELSE IF (GEKIN.EQ.0..AND.
     +   (IPART.EQ.8.OR.IPART.EQ.12.OR.IPART.EQ.11.OR.IPART.EQ.23)) THEN
         CALL GDECAY
         NMEC=NMEC+1
         LMEC(NMEC)=5
         ISTOP=1
         GOTO 999
      ENDIF
*
      IF (IFINIT(5) .EQ. 0) CALL FLINIT
      INT=0
      IJ=IGTOFL(IPART)
      IF(IJ.EQ.0) GOTO 110
      NMEC = NMEC + 1
      EKE = GEKIN
      TXI = VECT(4)
      TYI = VECT(5)
      TZI = VECT(6)
      DMOD = ONE/SQRT(TXI**2+TYI**2+TZI**2)
      TXI = TXI*DMOD
      TYI = TYI*DMOD
      TZI = TZI*DMOD
      WE  = 1.
      JMA = LQ(JMATE-NMAT)
      NCOMP = ABS (Q(JMA+11))
      AMM = Q(JMA+6)
      JMIXT = LQ(JMA-5)
 
*    Antiprotons, antineutrons and pi- are sent to
*    eventv for annihilation
      IF (GEKIN.EQ.0..AND.
     +         (IPART.EQ.15.OR.IPART.EQ.9.OR.IPART.EQ.25)) THEN
         IF(NCOMP.LE.1) THEN
            AMSS(1) = Q(JMA+6)
            ZTAR(1) = Q(JMA+7)
            MSSNUM(1) = 0
            RHO(1) = Q(JMA+8)
         ELSE
            ZSAMP(1) = 0.
            DO 10 I=1,NCOMP
               ZSAMP(I+1) = ZSAMP(I) + Q(JMIXT+NCOMP+I)
   10       CONTINUE
            CALL GRNDM(RNDM,1)
            ZCONT=ZSAMP(NCOMP+1)*RNDM(1)
            DO 20 I=1,NCOMP
               IF(ZCONT.LE.ZSAMP(I+1)) GO TO 30
   20       CONTINUE
            I = NCOMP
   30       CONTINUE
            AMSS(1)   = Q(JMIXT+I)
            MSSNUM(1) = 0
            ZTAR(1)   = Q(JMIXT+NCOMP+I)
            RHO(1)    = Q(JMIXT+2*NCOMP+I)*DENS
         END IF
         EKE = 1E-9
         POO=SQRT(EKE*(EKE+2*AM(IJ)))
         CALL EVENTV(IJ,POO,EKE,TXI,TYI,TZI,WE,1)
         NMEC=NMEC+1
         LMEC(NMEC)=17
         GOTO 80
      ELSE IF (GEKIN.LE.CUTHAD .AND. ITRTYP.EQ.4) THEN
         DESTEP = DESTEP + GEKIN
         GEKIN  = 0.
         GETOT  = AMASS
         VECT(7) = 0.
         ISTOP = 1
         LMEC(NMEC)=30
         GO TO 110
      ENDIF
*
      CALL GRNDM(RNDM,1)
      RNDEVT=RNDM(1)

      IF (FSIG.LE.0.) THEN
         SRAT = 1.E7
      ELSE
         SRAT = SINE/FSIG
      ENDIF
      IF ( RNDEVT .GE. SRAT) THEN
         IF (GEKIN .GT. 0.02) THEN
            POO=SQRT(EKE*(EKE+2*AM(IJ)))
         ELSE
            GO TO 110
         END IF
         INT=1
         LMEC(NMEC)=13
         IF(NCOMP.LE.1) THEN
            CALL NUCREL(IJ,POO,EKE,TXI,TYI,TZI,AMM,WE)
         ELSE
            CALL GRNDM(RNDM,1)
            RCONT=ELXNOR*RNDM(1)
            DO 40  I=1,NCOMP
               IF(RCONT.LE.CABELX(I)) GO TO 50
   40       CONTINUE
            I=NCOMP
   50       CONTINUE
            CALL NUCREL(IJ,POO,EKE,TXI,TYI,TZI,ONE*Q(JMIXT+I),WE)
         END IF
      ELSE
         LMEC(NMEC)=20
         IF (IHADR.EQ.2) THEN
            ISTOP = 2
            DESTEP = DESTEP + GETOT
            GO TO 110
         ENDIF
         IF (GEKIN .GT. 0.02) THEN
            POO=SQRT(EKE*(EKE+2*AM(IJ)))
         ELSE
            IF ((IJ.EQ.2 .OR. IJ.EQ.9 .OR. IJ.EQ.14 .OR. IJ.EQ.16)
     +            .AND. GEKIN .GT. 0.0) THEN
               POO=SQRT(EKE*(EKE+2*AM(IJ)))
            ELSE
               NMEC=NMEC-1
               GO TO 110
            END IF
         END IF
         INT=2
         IF(NCOMP.LE.1) THEN
            AMSS(1) = Q(JMA+6)
            ZTAR(1) = Q(JMA+7)
            MSSNUM(1) = 0
            RHO(1) = Q(JMA+8)
         ELSE
            CALL GRNDM(RNDM,1)
            RCONT=ANXNOR*RNDM(1)
            DO 60  I=1,NCOMP
               IF(RCONT.LE.CABINX(I)) GO TO 70
   60       CONTINUE
            I=NCOMP
   70       CONTINUE
            AMSS(1)   = Q(JMIXT+I)
            MSSNUM(1) = 0
            ZTAR(1)   = Q(JMIXT+NCOMP+I)
            RHO(1)    = Q(JMIXT+2*NCOMP+I)*DENS
         END IF
         CALL EVENTV(IJ,POO,EKE,TXI,TYI,TZI,WE,1)
      END IF
*
   80 IF(NP.EQ.1.AND.NPHEAV.EQ.0.AND.KPART(1).EQ.IJ) THEN
         VECT(4)=CXR(1)
         VECT(5)=CYR(1)
         VECT(6)=CZR(1)
         VECT(7)=SQRT(TKI(1)*(TKI(1)+2*AMASS))
         GETOT=TKI(1)+AMASS
         GEKIN=TKI(1)
      ELSE
         ISTOP=1
         NSTAK1 = MIN(NP,MXGKIN-NGKINE)
         IF(NP.GT.NSTAK1) THEN
            WRITE(CHMAIL,10000) NP-NSTAK1
            CALL GMAIL(0,0)
         ENDIF
         DO 90  K=1,NSTAK1
            NGKINE = NGKINE + 1
            IF (KPART(K) .EQ. 24 .OR. KPART(K) .EQ. 25) THEN
               KPART(K) = 19
               CALL GRNDM(RNDM,1)
               IF (RNDM(1) .GT. 0.5) KPART(K) = 12
            END IF
            IGEPAR = IFLTOG(KPART(K))
            JPA = LQ(JPART-IGEPAR)
            AGEMAS = Q(JPA+7)
            PGEANT = SQRT(TKI(K)*(TKI(K)+2*AGEMAS))
            GKIN(1,NGKINE)=CXR(K)*PGEANT
            GKIN(2,NGKINE)=CYR(K)*PGEANT
            GKIN(3,NGKINE)=CZR(K)*PGEANT
            GKIN(4,NGKINE)=TKI(K)+AGEMAS
            GKIN(5,NGKINE)=IGEPAR
            TOFD(NGKINE)=0.0
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
   90    CONTINUE
*
         NSTAK2 = MIN(NPHEAV,MXGKIN-NGKINE)
         IF(NPHEAV.GT.NSTAK2) THEN
            WRITE(CHMAIL,10100) NPHEAV-NSTAK2
            CALL GMAIL(0,0)
         ENDIF
         DO 100 K=1,NSTAK2
            NGKINE = NGKINE + 1
            IGEPAR = IHVTOG(KHEAVY(K))
            JPA = LQ(JPART-IGEPAR)
            AGEMAS = Q(JPA+7)
            PGEANT = SQRT(TKHEAV(K)*(TKHEAV(K)+2*AGEMAS))
            GKIN(1,NGKINE)=CXHEAV(K)*PGEANT
            GKIN(2,NGKINE)=CYHEAV(K)*PGEANT
            GKIN(3,NGKINE)=CZHEAV(K)*PGEANT
            GKIN(4,NGKINE)=TKHEAV(K)+AGEMAS
            GKIN(5,NGKINE)=IGEPAR
            TOFD(NGKINE)=0.0
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
  100    CONTINUE
*
         KCASE=NAMEC(12)
      END IF
  110 CONTINUE
      ZINTHA = GARNDM(DUMMY)
      SLHADR = SLENG
      STEPHA = 1.0E10
10000 FORMAT(' **** FLUFIN: Stack overflow, ',I6,' particles lost')
10100 FORMAT(' **** FLUFIN: Stack overflow, ',I6,
     +' heavy particles lost')
  999 END
