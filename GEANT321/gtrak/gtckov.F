*PATCH,GEXAM1.
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*DECK, GTCKOV
*CMZ :          12/09/95  11.07.14  by  S.Ravndal
*CMZ :  3.21/04 13/12/94  15.17.13  by  S.Giani
*-- Author :
      SUBROUTINE GTCKOV
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   This routine is called to follow the Cherenkov photons       *
C.    *   created during the tracking of charged particles and         *
C.    *   simulate the relevant processes along the way, until either  *
C.    *   the photon is absorbed or exits the detector. Processes      *
C.    *   currently simulated are absorption in-flight, and reflection *
C.    *   /transmission/absorption at a medium boundary. There are two *
C.    *   boundary types: dielectric-metal and dielectric-dielectric.  *
C.    *   For each of these there is a continuum of reflectivity       *
C.    *   and of surface quality from mirror finish to matte. The      *
C.    *   surface model is contained in routine GHSURF.                *
C.    *                                                                *
C.    *   ==>Called by : GTRACK                                        *
C.    *      Authors    F.Carminati, R.Jones ************              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcvol1.inc"
#include "geant321/gcunit.inc"

#if !defined(CERNLIB_OLD)
#include "geant321/gcvdma.inc"
#endif

*
* ** The following common is in GTMEDI. LSAMVL is set to true if
* ** we did not change volume yet
*
      COMMON/GCCHAN/LSAMVL
      LOGICAL LSAMVL
*
      DIMENSION R(3),D(3),U(3),QQ(3),vin(3)
      DIMENSION VBOU(3)
      LOGICAL LOLDTR
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
      PARAMETER (MXPUSH=10)
      SAVE RIN1,EFFIC
C.
C.    ------------------------------------------------------------------
*
* *** Update local pointers if medium has changed
*
      LOLDTR=.FALSE.
      IF(IUPD.EQ.0)THEN
         IUPD  = 1
         JTCKOV = LQ(JTM-3)
         IF(JTCKOV.EQ.0) THEN
*
* *** This Cerenkov photon has crossed into a black medium.
* *** Just absorb it.
            IPROC = 101
            STEP  = 0.
            SLABS = 0.
            ISTOP = 2
            DESTEP = 0.
            GOTO 110
         ENDIF
         NPCKOV = Q(JTCKOV+1)
         JABSCO = LQ(JTCKOV-1)
         JEFFIC = LQ(JTCKOV-2)
         JINDEX = LQ(JTCKOV-3)
         JPOLAR = LQ(JSTAK-1)
      ENDIF
      IF(SLENG.LE.0.) THEN
*
* *** Calculate GEKRAT for the particle
         IF(VECT(7).GE.Q(JTCKOV+NPCKOV+1)) THEN
            GEKRAT=1.
            IEKBIN=NPCKOV-1
         ELSEIF(VECT(7).LT.Q(JTCKOV+2)) THEN
*
* *** Particle below energy threshold ?  Short circuit
* *** This should never happen because the photons are generated
* *** only above threshold
*
*            GEKIN = 0.
*            GETOT = 0.
*            VECT(7)= 0.
*            ISTOP = 2
*            NMEC = 1
*            LMEC(1)= 30
*            GO TO 110

            GEKRAT=0.
            IEKBIN=1
         ELSE
            JMIN = 1
            JMAX = NPCKOV
   10       JMED = (JMIN+JMAX)/2
            IF(Q(JTCKOV+JMED+1).LT.VECT(7)) THEN
               JMIN = JMED
            ELSE
               JMAX = JMED
            ENDIF
            IF(JMAX-JMIN.GT.1) GO TO 10
            IEKBIN = JMIN
            GEKRAT = (VECT(7) - Q(JTCKOV+IEKBIN+1))/
     +      (Q(JTCKOV+IEKBIN+2)-Q(JTCKOV+IEKBIN+1))
         ENDIF
         GEKRT1=1.-GEKRAT
         RIN1=Q(JINDEX+IEKBIN)*GEKRT1+Q(JINDEX+IEKBIN+1)*GEKRAT
         EFFIC=Q(JEFFIC+IEKBIN)*GEKRT1+Q(JEFFIC+IEKBIN+1)*GEKRAT
         STEPLA=Q(JABSCO+IEKBIN)*GEKRT1+Q(JABSCO+IEKBIN+1)*GEKRAT
      ENDIF
*
* *** Compute current step size
*
      IPROC  = 103
      STEP   = STEMAX
*
*  **   Step limitation due to in flight absorbtion ?
*
      IF (ILABS.GT.0) THEN
         SLABS  = STEPLA*ZINTLA
         IF (SLABS.LT.STEP) THEN
            STEP  = SLABS
            IPROC = 101
         ENDIF
      ENDIF
*
      IF (STEP.LT.0.) STEP = 0.
*
*  **   Step limitation due to geometry ?
*
      STEPT=0.
      IF (STEP.GE.SAFETY) THEN
         CALL GTNEXT
         IF (IGNEXT.NE.0) THEN
*
* **    We are going to cross a boundary, so we need to simulate
* **    boundary effects and to know what is on the other side.
* **    For the moment save the current vector in the geometry tree.
*
#if !defined(CERNLIB_OLD)
         if(mycoun.gt.1.and.nfmany.gt.0)then
           nlevel=manyle(nfmany)
           do 99 i=1,nlevel
             names(i)=manyna(nfmany,i)
             number(i)=manynu(nfmany,i)
 99        continue
           call glvolu(nlevel,names,number,ier)
           if(ier.ne.0)print *,'Fatal error in GLVOLU'
           ingoto=0
         endif
#endif
            NLEVL1 = NLEVEL
            DO 20   I=1,NLEVEL
               NAMES1(I) = NAMES(I)
               NUMBR1(I) = NUMBER(I)
               LVOLU1(I) = LVOLUM(I)
   20       CONTINUE
*
* *** This is different from the other tracking routines.
* *** We get to the boundary and then we just jump over it
* *** So, linear transport till we are very near the boundary
*
#if !defined(CERNLIB_IBM)
            STEP = MAX(SNEXT-PREC,0.)
#endif
#if defined(CERNLIB_IBM)
            STEP = MAX(SNEXT-2.*PREC,0.)
#endif
C
C In case of Cherenkovs beeing near a corner, particle holding is
C prevented by a bigger step size. The value of VBOU is only used to
C calculate the correct surface normal (by GLISUR and GGPERP). The
C tracking of the Cherenkov is still using STEP
C
            IF(SNEXT.GT.0.) THEN
               DO 25 I=1,3
                  VBOU(I)=VECT(I)+SNEXT*VECT(I+3)
   25          CONTINUE
            END IF
C
            IF(STEP.GT.0.) THEN
               DO 30 I=1,3
                  VECT(I)=VECT(I)+STEP*VECT(I+3)
   30          CONTINUE
            ENDIF
            STEPT=STEP
#if !defined(CERNLIB_IBM)
            STEP  = SNEXT - STEP + PREC
#endif
#if defined(CERNLIB_IBM)
            STEP  = SNEXT - STEP + 2.*PREC
#endif
            IPROC = 0
            INWVOL= 2
            NMEC =  1
            LMEC(1)=1
         ENDIF
*
*        Update SAFETY in stack companions, if any
*        This may well not work.
         IF (IQ(JSTAK+3).NE.0) THEN
            DO 40 IST = IQ(JSTAK+3),IQ(JSTAK+1)
               JST = JSTAK +3 +(IST-1)*NWSTAK
               Q(JST+11) = SAFETY
   40       CONTINUE
            IQ(JSTAK+3) = 0
         ENDIF
*
      ELSE
         IQ(JSTAK+3) = 0
      ENDIF
*
* *** Linear transport
*
      VIN(1) = VECT(1)
      VIN(2) = VECT(2)
      VIN(3) = VECT(3)
      IF (INWVOL.EQ.2) THEN
         NBPUSH = 0
   50    DO 60 I = 1,3
            VECTMP  = VECT(I) +STEP*VECT(I+3)
            IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
               IF(VECT(I+3).NE.0.) THEN
                  VECTMP = VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*
     +            EPSMAC
                  IF(NMEC.GT.0) THEN
                     IF(LMEC(NMEC).EQ.104) NMEC=NMEC-1
                  ENDIF
                  NMEC=NMEC+1
                  LMEC(NMEC)=104
#if defined(CERNLIB_DEBUG)
                  WRITE(CHMAIL, 10000)
                  CALL GMAIL(0,0)
                  WRITE(CHMAIL, 10100) GEKIN, NUMED, STEP, SNEXT
                  CALL GMAIL(0,0)
10000 FORMAT(' Boundary correction in GTCKOV: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 FORMAT(31X,E10.3,1X,I10,1X,E10.3,1X,E10.3,1X)
#endif
               ENDIF
            ENDIF
            VOUT(I) = VECTMP
   60    CONTINUE
         NUMED1=NUMED
         CALL GTMEDI(VOUT,NUMED2)
         LOLDTR=.TRUE.
         IF(NUMED2.LE.0) THEN
            VECT(1) = VOUT(1)
            VECT(2) = VOUT(2)
            VECT(3) = VOUT(3)
            GO TO 110
         ENDIF
         JVO=LQ(JVOLUM-LVOLUM(NLEVEL))
         IF(LSAMVL.AND.Q(JVO+2).NE.12.) THEN
*
* *** In spite of our efforts we have not crossed the boundary
* *** we increase the step size and try again
*
            NBPUSH = NBPUSH + 1
            IF (NBPUSH.LE.MXPUSH) THEN
              STEP = STEP + NBPUSH*PREC
              GOTO 50
            ELSE
              INWVOL = 0
            ENDIF

         ENDIF
         IF(NUMED1.EQ.NUMED2) THEN
*
* *** If we are in the same medium, nothing needs to be done!
*
            VECT(1)=VOUT(1)
            VECT(2)=VOUT(2)
            VECT(3)=VOUT(3)
            IPROC=0
         ELSE
            JTM2 = LQ(JTMED-NUMED2)
            IF(IQ(JTM2-2).GE.3) THEN
               JTCKV2 = LQ(JTM2-3)
            ELSE
               JTCKV2 = 0
            ENDIF
            IF(JTCKV2.GT.0) THEN
               NPCKV2 = Q(JTCKV2+1)
               JABSC2 = LQ(JTCKV2-1)
               JEFFI2 = LQ(JTCKV2-2)
               JINDX2 = LQ(JTCKV2-3)
               IF(VECT(7).GE.Q(JTCKV2+NPCKV2+1)) THEN
                  GEKRT2=1.
                  IEKBI2=NPCKV2-1
               ELSEIF(VECT(7).LT.Q(JTCKV2+2)) THEN
                  GEKRT2=0.
                  IEKBI2=1
               ELSE
                  JMIN = 1
                  JMAX = NPCKV2
   64             JMED = (JMIN+JMAX)/2
                  IF(Q(JTCKV2+JMED+1).LT.VECT(7)) THEN
                     JMIN = JMED
                  ELSE
                     JMAX = JMED
                  ENDIF
                  IF(JMAX-JMIN.GT.1) GO TO 64
                  IEKBI2 = JMIN
                  GEKRT2 = (VECT(7) - Q(JTCKV2+IEKBI2+1))/
     +            (Q(JTCKV2+IEKBI2+2)-Q(JTCKV2+IEKBI2+1))
               ENDIF
               GEKRT1=1.-GEKRT2
               ABSCO2=Q(JABSC2+IEKBI2)*GEKRT1+Q(JABSC2+IEKBI2+1)*GEKRT2
               EFFIC2=Q(JEFFI2+IEKBI2)*GEKRT1+Q(JEFFI2+IEKBI2+1)*GEKRT2
               IF(JINDX2.GT.0) THEN
                 RIN2=Q(JINDX2+IEKBI2)*GEKRT1+Q(JINDX2+IEKBI2+1)*GEKRT2
               ELSE
                 RIN2=0.
               ENDIF
               IPROC = 102
            ELSE
               ISTOP=2
               DESTEP=0
               NMEC = NMEC+1
               LMEC(NMEC)=30
               VECT(1)=VOUT(1)
               VECT(2)=VOUT(2)
               VECT(3)=VOUT(3)
               GOTO 110
            ENDIF
         ENDIF
      ELSE
         DO 70 I = 1,3
            VECT(I)  = VECT(I) +STEP*VECT(I+3)
   70    CONTINUE
      ENDIF
*
      STEP = STEPT + STEP
      SLENG = SLENG +STEP
*
* *** Update time of flight
*
      TOFG = TOFG +STEP*RIN1/CLIGHT
*
* *** Update interaction probabilities
*
      IF (ILABS.GT.0)    ZINTLA = ZINTLA -STEP/STEPLA
*
      IF (IPROC.EQ.0) GO TO 110
      NMEC = NMEC+1
      LMEC(NMEC) = IPROC
*
*  ** Absorbtion in flight ?
*
      IF (IPROC.EQ.101) THEN
         ISTOP=2
         CALL GRNDM(RNDM,1)
         IF(RNDM.LT.EFFIC) THEN
*
* ***  Destep =/= 0 means that the photon has been detected
*
            DESTEP=VECT(7)
         ELSE
            DESTEP=0.
         ENDIF
*
*  ** Boundary action?
*
      ELSE IF (IPROC.EQ.102) THEN
         IF(JINDX2.EQ.0) THEN
*
* *** Case dielectric -> metal
*
            CALL GRNDM(RNDM,1)
            IF(RNDM.LT.ABSCO2) THEN
*
* *** Photon is absorbed in the next medium
*
               NMEC=NMEC+1
               LMEC(NMEC)=101
               ISTOP = 2
               CALL GRNDM(RNDM,1)
               IF(RNDM.LT.EFFIC2) THEN
*
* ***  Destep =/= 0 means that the photon has been detected
*
                  DESTEP=VECT(7)
               ELSE
                  DESTEP = 0.
               END IF
               VECT(1) = VOUT(1)
               VECT(2) = VOUT(2)
               VECT(3) = VOUT(3)
               GOTO 110
            ELSE
*
* *** Photon is reflected (no polarization effects)
*
               CALL GLISUR(VECT,VBOU,NUMED1,NUMED2,U,PDOTU,IERR)
               IF (IERR.NE.0) THEN
                  WRITE(CHMAIL,10200) IERR
                  CALL GMAIL(0,0)
                  GO TO 110
               END IF
*
* ** Restore old volume tree, the photon does not cross the boundary
*
*              CALL GLVOLU(NLEVL1,NAMES1,NUMBR1,IERR)
               CALL GTMEDI(VIN,N)
               LOLDTR=.FALSE.
*
*
* *** Reflect, but make sure we are in the old volume before
*
               NBPUSH = 0
   80          NBPUSH = NBPUSH+1
               IF(NBPUSH.GT.MXPUSH) THEN
                  WRITE(CHMAIL,10300) NTMULT,NSTEP
                  CALL GMAIL(0,0)
                  ISTOP=1
                  GOTO 110
               ELSE
                  CALL GINVOL(VECT,ISAME)
                  IF(ISAME.EQ.0) THEN
                     PRECN = NBPUSH*PREC
                     VECT(1) = VECT(1) - PRECN*VECT(4)
                     VECT(2) = VECT(2) - PRECN*VECT(5)
                     VECT(3) = VECT(3) - PRECN*VECT(6)
                     GO TO 80
                  ENDIF
               ENDIF
*
               NMEC=NMEC+1
               LMEC(NMEC)=106
               EDOTU = POLAR(1)*U(1) + POLAR(2)*U(2) + POLAR(3)*U(3)
               VECT(4) = +VECT(4) - 2*PDOTU*U(1)
               VECT(5) = +VECT(5) - 2*PDOTU*U(2)
               VECT(6) = +VECT(6) - 2*PDOTU*U(3)
               POLAR(1) = -POLAR(1) + 2*EDOTU*U(1)
               POLAR(2) = -POLAR(2) + 2*EDOTU*U(2)
               POLAR(3) = -POLAR(3) + 2*EDOTU*U(3)
               INWVOL  = 0
            ENDIF
         ELSE
*
*  case dielectric-dielectric:
*
            CALL GLISUR(VECT,VBOU,NUMED1,NUMED2,U,PDOTU,IERR)
            IF (IERR.NE.0) THEN
               WRITE(CHMAIL,10200) IERR
               CALL GMAIL(0,0)
               GO TO 110
            END IF
            EDOTU = POLAR(1)*U(1) + POLAR(2)*U(2) + POLAR(3)*U(3)
            COST1 = -PDOTU
            IF (ABS(COST1).LT.1.) THEN
               SINT1 = SQRT((1-COST1)*(1+COST1))
               SINT2 = SINT1*RIN1/RIN2
            ELSE
               SINT1 = 0.0
               SINT2 = 0.0
            END IF
            IF (SINT2.GE.1) THEN
*
* ***  Simulate total internal reflection
*
               NMEC=NMEC+1
               LMEC(NMEC)=106
*
* ** Restore old volume tree, the photon does not cross the boundary
*
*              CALL GLVOLU(NLEVL1,NAMES1,NUMBR1,IERR)
               CALL GTMEDI(VIN,N)
               LOLDTR=.FALSE.
*
* *** Reflect, but make sure we are in the old volume before
*
               NBPUSH = 0
   90          NBPUSH = NBPUSH+1
               IF(NBPUSH.GT.MXPUSH) THEN
                  WRITE(CHMAIL,10300) NTMULT,NSTEP
                  CALL GMAIL(0,0)
                  ISTOP=1
                  GOTO 110
               ELSE
                  CALL GINVOL(VECT,ISAME)
                  IF(ISAME.EQ.0) THEN
                     PRECN = NBPUSH*PREC
                     VECT(1) = VECT(1) - PRECN*VECT(4)
                     VECT(2) = VECT(2) - PRECN*VECT(5)
                     VECT(3) = VECT(3) - PRECN*VECT(6)
                     GO TO 90
                  ENDIF
               ENDIF
*
               VECT(4) = +VECT(4) - 2*PDOTU*U(1)
               VECT(5) = +VECT(5) - 2*PDOTU*U(2)
               VECT(6) = +VECT(6) - 2*PDOTU*U(3)
               POLAR(1) = -POLAR(1) + 2*EDOTU*U(1)
               POLAR(2) = -POLAR(2) + 2*EDOTU*U(2)
               POLAR(3) = -POLAR(3) + 2*EDOTU*U(3)
               INWVOL = 0
            ELSE
*
* ***  Calculate amplitude for transmission (Q = P x U)
*
               COST2 = SIGN(1.0,COST1)*SQRT((1-SINT2)*(1+SINT2))
               IF (SINT1.GT.1.E-4) THEN
                  QQ(1) = VECT(5)*U(3) - VECT(6)*U(2)
                  QQ(2) = VECT(6)*U(1) - VECT(4)*U(3)
                  QQ(3) = VECT(4)*U(2) - VECT(5)*U(1)
                  EPERP1 = (POLAR(1)*QQ(1) + POLAR(2)*QQ(2) + POLAR(3)*
     +            QQ(3))
                  ENORM  = SQRT(EPERP1**2+EDOTU**2)
                  EPERP1 = EPERP1/ENORM
                  EPARL1 = EDOTU/ENORM
               ELSE
                  QQ(1) = POLAR(1)
                  QQ(2) = POLAR(2)
                  QQ(3) = POLAR(3)
*
*     Here we follow Jackson's conventions and we set the parallel
*     component = 1 in case of a ray perpendicular to the surface
                  EPERP1 = 0.
                  EPARL1 = 1.
               END IF
               IF(COST1.NE.0.) THEN
                  S1 = RIN1*COST1
                  EPERP2 = 2*S1*EPERP1/(RIN1*COST1+RIN2*COST2)
                  EPARL2 = 2*S1*EPARL1/(RIN2*COST1+RIN1*COST2)
                  E2 = EPERP2**2 + EPARL2**2
                  S2 = RIN2*COST2*E2
                  TCOEF = S2/S1
               ELSE
                  TCOEF = 0.
               ENDIF
               CALL GRNDM(RNDM,1)
               IF (RNDM.GT.TCOEF) THEN
*
* *** Simulate reflection
*
                  NMEC=NMEC+1
                  LMEC(NMEC)=106
*
* *** Restore old volume tree, the photon does not cross the boundary
*
*                 CALL GLVOLU(NLEVL1,NAMES1,NUMBR1,IERR)
                  CALL GTMEDI(VIN,N)
                  LOLDTR=.FALSE.
*
* *** Reflect, but make sure we are in the old volume before
*
                  NBPUSH = 0
  100             NBPUSH = NBPUSH+1
                  IF(NBPUSH.GT.MXPUSH) THEN
                     WRITE(CHMAIL,10300) NTMULT,NSTEP
                     CALL GMAIL(0,0)
                     ISTOP=1
                     GOTO 110
                  ELSE
                     CALL GINVOL(VECT,ISAME)
                     IF(ISAME.EQ.0) THEN
                        PRECN = NBPUSH*PREC
                        VECT(1) = VECT(1) - PRECN*VECT(4)
                        VECT(2) = VECT(2) - PRECN*VECT(5)
                        VECT(3) = VECT(3) - PRECN*VECT(6)
                        GO TO 100
                     ENDIF
                  ENDIF
*
                  VECT(4) = VECT(4) - 2*PDOTU*U(1)
                  VECT(5) = VECT(5) - 2*PDOTU*U(2)
                  VECT(6) = VECT(6) - 2*PDOTU*U(3)
                  IF(SINT1.GT.1E-4) THEN
                     EPARL2 = RIN2*EPARL2/RIN1-EPARL1
                     EPERP2 = EPERP2-EPERP1
                     E2 = EPERP2**2 + EPARL2**2
                     R(1) = U(1) + PDOTU*VECT(4)
                     R(2) = U(2) + PDOTU*VECT(5)
                     R(3) = U(3) + PDOTU*VECT(6)
                     EABS = SQRT(E2)*SINT1
                     CPARL = EPARL2/EABS
                     CPERP = EPERP2/EABS
                     POLAR(1) = CPARL*R(1) - CPERP*QQ(1)
                     POLAR(2) = CPARL*R(2) - CPERP*QQ(2)
                     POLAR(3) = CPARL*R(3) - CPERP*QQ(3)
                  ELSEIF(RIN2.GT.RIN1) THEN
*
* *** Case of ray perpendicular to the surface. No change or
* *** an inversion of phase.
                     POLAR(1) = -POLAR(1)
                     POLAR(2) = -POLAR(2)
                     POLAR(3) = -POLAR(3)
                  ENDIF
                  INWVOL = 0
               ELSE
*
* *** Simulate transmission/refraction
*
                  NMEC=NMEC+1
                  LMEC(NMEC)=107
                  VECT(1) = VOUT(1)
                  VECT(2) = VOUT(2)
                  VECT(3) = VOUT(3)
                  IF(SINT1.GT.1E-4) THEN
                     ALPHA = COST1-COST2*(RIN2/RIN1)
                     D(1) = VECT(4) + ALPHA*U(1)
                     D(2) = VECT(5) + ALPHA*U(2)
                     D(3) = VECT(6) + ALPHA*U(3)
                     DABS = SQRT(D(1)**2+D(2)**2+D(3)**2)
                     VECT(4) = D(1)/DABS
                     VECT(5) = D(2)/DABS
                     VECT(6) = D(3)/DABS
                     PDOTU = -COST2
                     R(1) = U(1) - PDOTU*VECT(4)
                     R(2) = U(2) - PDOTU*VECT(5)
                     R(3) = U(3) - PDOTU*VECT(6)
                     EABS = SQRT(E2)
                     CPARL = EPARL2/(EABS*SINT2)
                     CPERP = EPERP2/(EABS*SINT1)
                     POLAR(1) = CPARL*R(1) + CPERP*QQ(1)
                     POLAR(2) = CPARL*R(2) + CPERP*QQ(2)
                     POLAR(3) = CPARL*R(3) + CPERP*QQ(3)
                  ENDIF
                  GEKRAT = GEKRT2
                  IEKBIN = IEKBI2
                  STEPLA = ABSCO2
                  EFFIC = EFFIC2
                  RIN1 = RIN2
               END IF
            END IF
         END IF
      ENDIF
*                                                             END GTCKOV
  110 IF(LOLDTR) CALL GLVOLU(NLEVL1,NAMES1,NUMBR1,IERR)
10200  FORMAT(' **** GTCKOV: error from GLISUR = ',I10)
10300  FORMAT(' **** GTCKOV: unable to reflect at NTMULT = ',
     +        I8,' step No. ',I8,' photon abandoned!')
      END
