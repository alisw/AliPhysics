*
* $Id$
*
* $Log$
* Revision 1.2  1997/04/29 09:54:58  japost
* Re-dimensioning and filling common arrays to fix bug/error.
*
* More details can be found in this error report, which caused the changes:
*
* Date: Tue, 15 Apr 1997 19:03:49 +0200
* Simone's forward:
* ----------
*
* > From: Dieter Heck <heck@ik3.fzk.de>
* > To: SGIANI@cernvm.cern.ch
* > Subject: GHEISHA_ERROR
* > Date: Tuesday, April 08, 1997 2:09 PM
* >
* >
* > Dear Simone Giani,
* >
* > again I have detected an error within one of the GHEISHA routines.
* > It is the routine GHETUN. This routine contains two arrays keeping the
* > strangenes numbers (SNUM) and baryon numbers (BNUM) of the particles to
* > be treated. Both arrays are dimensioned with 32, while in GHEISHA also
* > OMAGA(-) and ANTI-OMEAG(-) particles are admitted with particle codes 33
* > and 34. In case of treatment of omega particles, the array boundaries
* > are exceeded, which gives an error stop when running with the bounds
* > check option. To correct this error, the arrays must be dimensionend
* > with 34, and the DATA statements have to be filled with the
* > corresponding values:
* >
* >       DIMENSION RNDM(4),SNUM(34),BNUM(34),REDDEC(7)
* >       DATA SNUM/9*0.,1.,0.,0.,-1.,4*0.,-1.,1.,-1.,-1.,-1.,1.,1.,1.,
* >      $          -2.,-2.,2.,2.,3*0.,-3.,3./
* >       DATA BNUM/13*0.,1.,-1.,1.,-1.,1.,-1.,1.,1.,1.,-1.,-1.,-1.,
* >      $          1.,1.,-1.,-1.,2.,3.,4.,1.,-1./
* >
* > Hopefully this correction is helpful for the GHEISHA users.
* >
* > Kind regards,
* >
* >   Dieter Heck
* >
* > =========================================================================
* > From:   Dieter Heck,   Institut fuer Kernphysik 3,
* > Forschungszentrum Karlsruhe,  Postfach 3640,  D-76021 Karlsruhe,  Germany
* > Tel: (49) 7247-82-3777   Fax: (49) 7247-82-4075   e-mail: heck@ik3.fzk.de
* > =========================================================================
*
* Revision 1.1.1.1  1995/10/24 10:21:12  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.40  by  S.Giani
*-- Author :
      SUBROUTINE GHETUN(NT)
C**
C** TUNING OF THE HIGH ENERGY COLLISION MODEL:
C**
C** 1. AVOID THAT PI0 IS LEADING PARTICLE.
C** 2. SOME FINE-TUNING FOR THE NUMBER OF PRODUCED PROTONS AND
C**    NEUTRONS.
C** 3. INTRODUCE A FLAVOUR DEPENDENT CORRECTION FOR SINGLE PARTICLE
C**    SPECTRA.
C** 4. FINE-TUNING OF LEADING PARTICLE SPECTRA AND MOMENTUM
C**    CONSERVATION.
C**
#include "geant321/mxgkgh.inc"
#include "geant321/s_consts.inc"
#include "geant321/s_curpar.inc"
#include "geant321/s_result.inc"
#include "geant321/s_prntfl.inc"
#include "geant321/s_blank.inc"
C
C
      DIMENSION RNDM(4),SNUM(34),BNUM(34),REDDEC(7)
      DATA SNUM/9*0.,1.,0.,0.,-1.,4*0.,-1.,1.,-1.,-1.,-1.,1.,1.,1.,
     $          -2.,-2.,2.,2.,3*0.,-3.,3./
      DATA BNUM/13*0.,1.,-1.,1.,-1.,1.,-1.,1.,1.,1.,-1.,-1.,-1.,
     $          1.,1.,-1.,-1.,2.,3.,4.,1.,-1./
 

#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION AHMF,BHMF
#endif
C**
      MX=MXGKPV-10
      MX1=MX+1
      MX2=MX+2
      MX3=MX+3
      MX4=MX+4
      MX5=MX+5
      MX6=MX+6
      MX7=MX+7
      MX8=MX+8
      MX9=MX+9
      NT1=NT
      IF(NT1.GT.MXGKPV-10) NT1=MXGKPV-10
      NT=NT1
C
      CALL GRNDM(RNDM,1)
      IF(EK.LT.(25.+RNDM(1)*75.)) GOTO 15
C
C**  IF PI0 IS THE HIGHEST MOMENTUM PARTICLE, INTERCHANGE IT WITH A
C**  CHARGED PION.
C
      CALL GRNDM(RNDM,4)
      REDEN = -0.7 + 0.29*LOG10(EK)
      REDAT = 1. - 0.4000*LOG10(ATNO2)
      PMAX  = -200.
      PMAPIP= -200.
      PMAPI0= -200.
      PMAPIM= -200.
      IPMAX = 0
      MAXPIP= 0
      MAXPI0= 0
      MAXPIM= 0
      IF(RNDM(1).GT.(ATNO2/100.-0.28).AND.ABS(NCH).GT.0.5) THEN
         DO 46 I=1,NT1
            IPHMF=IFIX(PV(8,I)+0.1)
            CALL LENGTX(I,PPP)
            IF(PPP.GT.PMAX) THEN
               PMAX=PPP
               IPMAX=I
            ENDIF
            IF(IPHMF.EQ.7) THEN
               IF(PPP.GT.PMAPIP) THEN
                  PMAPIP=PPP
                  MAXPIP=I
               ENDIF
            ENDIF
            IF(IPHMF.EQ.8) THEN
               IF(PPP.GT.PMAPI0) THEN
                  PMAPI0=PPP
                  MAXPI0=I
               ENDIF
            ENDIF
            IF(IPHMF.EQ.9) THEN
               IF(PPP.GT.PMAPIM) THEN
                  PMAPIM=PPP
                  MAXPIM=I
               ENDIF
            ENDIF
   46    CONTINUE
      ENDIF
C**
C**   SOME ADDITIONAL TUNING OF THE NUMBER OF GREY TRACK PARTICLES
C**
      IF(NT1.GT.2) THEN
      DO 47 I=3,NT1
         IPHMF=IFIX(PV(8,I)+0.1)
         IF(IPHMF.EQ.14.OR.IPHMF.EQ.16.OR.IPHMF.GE.30) THEN
            CALL LENGTX(I,PPP)
            IF(PPP.LT.1.5) THEN
               IF(RNDM(2).LT.REDEN.OR.RNDM(3).LT.REDAT) THEN
                  PV(1,I) = 0.
                  PV(2,I) = 0.
                  PV(3,I) = 0.
                  PV(4,I) = ABS(PV(5,I))
               ENDIF
            ENDIF
         ENDIF
   47 CONTINUE
      ENDIF
C**
      IF(MAXPI0.EQ.0)    GOTO 10
      IF(PMAPI0.LT.PMAX) GOTO 10
      IF(RNDM(4).LT.PMAPI0/P) THEN
      IF(NCH.GT.0.5.AND.MAXPIP.NE.0) THEN
         DO 49 J=5,10
           PV(J,MX1)=PV(J,MAXPI0)
           PV(J,MAXPI0)=PV(J,MAXPIP)
           PV(J,MAXPIP)=PV(J,MX1)
   49    CONTINUE
      ENDIF
      IF(NCH.LT.-0.5.AND.MAXPIM.NE.0) THEN
         DO 56 J=5,10
           PV(J,MX1)=PV(J,MAXPI0)
           PV(J,MAXPI0)=PV(J,MAXPIM)
           PV(J,MAXPIM)=PV(J,MX1)
   56    CONTINUE
      ENDIF
      ENDIF
C
   10 CONTINUE
C**
C** CHECK TOTAL BARYON- NUMBER AND
C** SKIP ZERO MOMENTUM PARTICLES
C**
      BNTOT=-BNUM(IPART)-ATNO2
      DO 57 I=1,NT1
         IPHMF=IFIX(PV(8,I)+0.1)
         BNTOT=BNTOT+BNUM(IPHMF)
   57 CONTINUE
      BNTOT=1.+BNTOT/ATNO2
      IF(ATNO2.LT.1.5) BNTOT=0.
      CALL GRNDM(RNDM,1)
      IF(ATNO2.GT.(75.+RNDM(1)*25.)) BNTOT=0.
C**
      II=0
      DO 12 I=1,NT1
         CALL LENGTX(I,PPP)
         IF(PPP.GT.1.E-6) THEN
            IPHMF=IFIX(PV(8,I)+0.1)
            IF(BNTOT.GT.0.3) THEN
            IF(IPHMF.EQ.14.OR.IPHMF.EQ.16.OR.IPHMF.GE.30) THEN
               CALL GRNDM(RNDM,1)
               IF(RNDM(1).LT.0.5.AND.PPP.LT.1.2) GOTO 12
            ENDIF
            ENDIF
            II=II+1
            DO 11 J=1,10
               PV(J,II)=PV(J,I)
   11       CONTINUE
         ENDIF
   12 CONTINUE
      NT1=II
      NT=NT1
C**
C**   EXACT MOMENTUM CONSERVATION AND SOME CORRECTIONS FOR SINGLE
C**   PARTICLE SPECTRA AT HIGH ENERGIES ONLY
C
   15 PV(1,MX1) = P*PX
      PV(2,MX1) = P*PY
      PV(3,MX1) = P*PZ
      PV(4,MX1) = EN
      PV(5,MX1) = ABS(AMAS)
      PV(6,MX1) = NCH
      PV(1,MX2) = 0.
      PV(2,MX2) = 0.
      PV(3,MX2) = 0.
      PV(4,MX2) = MP
      PV(5,MX2) = MP
      PV(6,MX2) = 0.
C
      IF(NPRT(4)) THEN
         WRITE(NEWBCD,2000)
         WRITE(NEWBCD,2001) MX1,(PV(J,MX1),J=1,6)
         WRITE(NEWBCD,2001) MX2,(PV(J,MX2),J=1,6)
      ENDIF
C
      DO 58 J=1,10
         PV(J,MX9) = 0.
   58 CONTINUE
      CALL ADD(MX1,MX2,MX3)
      CALL LOR(MX1,MX3,MX4)
      CALL LOR(MX2,MX3,MX5)
      LEDPAR=0
      REDPAR=0.
      GESPAR=0.
      SNUM1=SNUM(IPART)
      IF(IPART.EQ.11.OR.IPART.EQ.12) THEN
        CALL GRNDM(RNDM,1)
        SNUM1=1.
        IF(RNDM(1).LT.0.5) SNUM1=-1.
      ENDIF
      DO 20 I=1,NT1
         IPHMF=IFIX(PV(8,I)+0.1)
         IF(IPHMF.LE.6.OR.IPHMF.GT.32) GOTO 19
         CALL LENGTX(I,PPP)
         IF(PPP.LT.1.E-3) GOTO 19
         CALL LOR(I,MX3,MX6)
         CALL ANG(MX4,MX6,COST,TETA)
         SNUM2=SNUM(IPHMF)
         IF(IPHMF.EQ.11.OR.IPHMF.EQ.12) THEN
            CALL GRNDM(RNDM,1)
            SNUM2=1.
            IF(RNDM(1).LT.0.5) SNUM2=-1.
         ENDIF
         IF(COST.GT.0.) THEN
            HFMAS=ABS(AMAS)
            REDDEC(1)=ABS(HFMAS    -ABS(PV(5,I)))
            REDDEC(2)=ABS(NCH-PV(6,I))
            REDDEC(3)=ABS(SNUM1      -SNUM2)
            REDDEC(4)=ABS(BNUM(IPART)-BNUM(IPHMF))
         ELSE
            HFMAS=MP
            REDDEC(1)=ABS(HFMAS     -ABS(PV(5,I)))
            REDDEC(2)=ABS(ZNO2/ATNO2-PV(6,I))
            REDDEC(3)=ABS(SNUM2)
            REDDEC(4)=ABS(1.-BNUM(IPHMF))
         ENDIF
         REDDEC(6)=REDDEC(1)+REDDEC(2)+REDDEC(3)+REDDEC(4)
         SBQWGT=REDDEC(6)
         IF(HFMAS.LT.0.2) THEN
            SBQWGT=(SBQWGT-2.5)*0.10
            IF(IPHMF.EQ.8) SBQWGT=0.15
         ELSE IF (HFMAS.LT.0.6) THEN
            SBQWGT=(SBQWGT-3.0)*0.10
         ELSE
            SBQWGT=(SBQWGT-2.0)*0.10
            IF(IPHMF.EQ.8) SBQWGT=0.15
         ENDIF
         CALL LENGTX(MX6,PPP)
         IF(SBQWGT.GT.0. .AND. PPP.GT.1.E-6) THEN
         PLHMF=PPP*COST
         PTHMF=PPP*SQRT(1.-COST*COST)
         PLHMF=PLHMF*(1.-SBQWGT)
         CALL CROSS3(MX4,MX6,MX8)
         CALL CROSS3(MX8,MX4,MX8)
         CALL LENGTX(MX4,PPP)
         PV(1,MX7)=PV(1,MX4)*PLHMF/PPP
         PV(2,MX7)=PV(2,MX4)*PLHMF/PPP
         PV(3,MX7)=PV(3,MX4)*PLHMF/PPP
         CALL LENGTX(MX8,PPP)
         PV(1,MX8)=PV(1,MX8)*PTHMF/PPP
         PV(2,MX8)=PV(2,MX8)*PTHMF/PPP
         PV(3,MX8)=PV(3,MX8)*PTHMF/PPP
         CALL ADD3(MX7,MX8,MX6)
         CALL LENGTX(MX6,PPP)
         AHMF=PPP
         BHMF=PV(5,I)
         PV(4,MX6)=DSQRT(AHMF**2+BHMF**2)
C
         IF(NPRT(4))
     $      WRITE(NEWBCD,3001) I,(PV(J,I),J=1,8),SBQWGT
C
         CALL LOR(MX6,MX5,I)
C
         IF(NPRT(4))
     $      WRITE(NEWBCD,3001) I,(PV(J,I),J=1,8),SBQWGT
         ENDIF
C
         IF(IPHMF.EQ.8) GOTO 19
         CALL SUB3(MXGKPV,I,MX8)
         CALL LENGTX(MX8,PPP)
         REDDEC(5) = PPP/P
         REDDEC(7)=REDDEC(5)*2./3. + REDDEC(6)/12.
         REDDEC(7) =  1.-REDDEC(7)
         IF(REDDEC(7) .LT.  0.) REDDEC(7) =  0.
         GESPAR=GESPAR+REDDEC(7)
         IF(REDDEC(6).LE.3.75) THEN
          IF(REDDEC(7) .GT.  REDPAR) THEN
             LEDPAR=I
             REDPAR=REDDEC(7)
          ENDIF
         ENDIF
         IF(NPRT(4))
     $      WRITE(NEWBCD,2001) I,(PV(J,I),J=1,6),PV(8,I),REDDEC
C
  19     CALL ADD3(MX9,I,MX9)
C
  20  CONTINUE
      IF(NPRT(4))
     $   WRITE(NEWBCD,1001) LEDPAR,REDPAR,GESPAR
C**
C** APPLY CORRECTION ON LEADING PARTICLE
C**
      CALL GHEPEC(LEDPAR)
C**
      RETURN
 1001 FORMAT(1H0,'SEARCH FOR LEADING PARTICLE, WEIGHT, TOTAL WEIGHT ',
     $ I5,3X,2F10.4)
 2000 FORMAT(1H ,'MOMENTUM CONSERVATION AT HIGH ENERGIES: ')
 2001 FORMAT(1H ,I3,2X,7F8.3/1H ,5X,7F8.3)
 3001 FORMAT(1H ,I3,2X,5F8.3,F5.1,F8.3,F5.1,F8.3)
      END
