*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GDRACK
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       RAY-TRACING                                              *
C.    *       Computation and tracking of all the light rays.          *
C.    *                                                                *
C.    *       Author: S.Giani                                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcflag.inc"
********************************************************************************
#include "geant321/gcdraw.inc"
#include "geant321/gconst.inc"
#include "geant321/gcmutr.inc"
#include "geant321/pawc.inc"
********************************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcparm.inc"
#include "geant321/gcsets.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
********************************************************************************
#include "geant321/gcfdim.inc"
#include "geant321/gcpixe.inc"
#include "geant321/gcrayt.inc"
#include "geant321/gcvdma.inc"
 
      CHARACTER*4 NAME
      DIMENSION VVVMIN(80),VVVMAX(80)
      SAVE VVVMIN,VVVMAX,ZIMPRE
********************************************************************************
      COMMON/GCCHAN/LSAMVL
      LOGICAL LSAMVL
*
      DIMENSION CUTS(10),MECA(5,13)
      EQUIVALENCE (CUTS(1),CUTGAM),(MECA(1,1),IPAIR)
      SAVE PRECOR
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
C.
C.    ------------------------------------------------------------------
         FINVIS=1./ISFILL
         LLL=1
         ZLN=0.
         IF(IMAP.EQ.1)ZIMPRE=LIMPRE
         IF(IMAP.EQ.2)THEN
           ZRATIO=LIMPRE/ZIMPRE
         ENDIF
         DO 21 IIIJ=1,LIMPRE
           IF(IMAP.EQ.1)THEN
            VVVMIN(IIIJ)=10000.
            VVVMAX(IIIJ)=0.
           ELSEIF(IMAP.EQ.2)THEN
            ZLN=ZLN+1
            IF((ZLN-ZRATIO).GT.0.001)THEN
              ZLN=ZLN-ZRATIO
              LLL=LLL+1
            ENDIF
           ENDIF
           VVV=FINVIS
           IYYY=1
           IF(IIIJ.NE.1)THEN
             UUU=UUU+FINVIS
             IXXX=IXXX+1
           ENDIF
          DO 23 JJJI=1,LIMPRE
           IF(IMAP.EQ.2)THEN
            IF(LLL.EQ.1)THEN
             VVMA=VVVMAX(LLL)
             VVMI=VVVMIN(LLL)
            ELSE
             VVMA=MAX(VVVMAX(LLL-1),VVVMAX(LLL))
             VVMI=MIN(VVVMIN(LLL-1),VVVMIN(LLL))
            ENDIF
            IF(((VVV-(VVMA+ZNMAP1)).GT.0.001).OR.
     +        ((VVV-(VVMI-ZNMAP1)).LT.-0.001))GOTO 22
           ENDIF
            IF(IIIJ.EQ.1.AND.JJJI.EQ.1)THEN
               MYTRME=NUMED
               ime=0
*               print *,vect(1),vect(2),vect(3),'vertex from gtrack'
*               print *,vect(4),vect(5),vect(6),'impulse from gtrack'
               GOTO 9
            ENDIF
                                 XPINTS=ZROTS(1,4)+ZROTS(1,1)*
     +                           UUU+ZROTS(1,2)*VVV+ZROTS(1,3)*
     +                           ZUV
                                 YPINTS=ZROTS(2,4)+ZROTS(2,1)*
     +                           UUU+ZROTS(2,2)*VVV+ZROTS(2,3)*
     +                           ZUV
                                 ZPINTS=ZROTS(3,4)+ZROTS(3,1)*
     +                           UUU+ZROTS(3,2)*VVV+ZROTS(3,3)*
     +                           ZUV
            JON=0
            ISSEEN=0
            IME=0
            SLENG=0.
            NUMED=MYTRME
            NSTEP=0
            INFROM=0
            nlevel=1
       CALL UCTOH('PERS',IPERS,4,4)
       XCOSXS=(SIN(GTHETA*DEGRAD))*(COS(GPHI*DEGRAD))
       YCOSYS=(SIN(GTHETA*DEGRAD))*(SIN(GPHI*DEGRAD))
       ZCOSZS=COS(GTHETA*DEGRAD)
       VDX=XCOSXS
       VDY=YCOSYS
       VDZ=ZCOSZS
       VECT(1)=XPINTS
       VECT(2)=YPINTS
       VECT(3)=ZPINTS
       IF(IPERS.EQ.IPRJ)THEN
        CONMOD=1./SQRT(((XPINTS-FPINTX)**2)+((YPINTS-FPINTY)**2)+
     +                 ((ZPINTS-FPINTZ)**2))
        XCOSXS=-(XPINTS-FPINTX)*CONMOD
        YCOSYS=-(YPINTS-FPINTY)*CONMOD
        ZCOSZS=-(ZPINTS-FPINTZ)*CONMOD
       ENDIF
       VECT(4) = -XCOSXS
       VECT(5) = -YCOSYS
       VECT(6) = -ZCOSZS
       IF(VECT(1).LE.XCUT.OR.VECT(2).LE.YCUT.OR.VECT(3).LE.ZCUT)THEN
         CALL GTMEDI (VECT, NUMED)
       ENDIF
 9    CONTINUE
 
      ISTOP = 0
      EPSCUR = EPSMAC
      NSTOUT = 0
      INWOLD = 0
      LSAMVL = .FALSE.
      NUMOLD=0
*
* *** Check validity of tracking medium and material parameters
*
   10 IF (NUMED.NE.NUMOLD) THEN
         NUMOLD = NUMED
         IUPD   = 0
         JTM    = LQ(JTMED- NUMED)
         EPSIL    = Q(JTM + 13)
         PRECOR   = MIN(0.1*EPSIL, 0.0010)
*
         NMAT     = Q(JTM + 6)
         JMA   = LQ(JMATE-NMAT)
         DENS = Q(JMA +8)
      ENDIF
        IF(ISCOLO.EQ.-10.OR.IMYSE.EQ.0)ISCOLO=MOD(NUMED,6)+2
        IF(ISLSTY.EQ.-10.OR.IMYSE.EQ.0)ISLSTY=4
        IF(ISSEEN.EQ.-10.OR.IMYSE.EQ.0)THEN
         IF(IME.EQ.1)THEN
         IF(DENS.LT.0.00130)THEN
          ISSEEN=0
         ELSE
          ISSEEN=1
         ENDIF
         ENDIF
        ENDIF
        IF(IME.EQ.0)IME=1
*
      IF(LSAMVL) THEN
*
*       If now the particle is entering in the same volume where
*       it was exiting from last step, and if it has done this for
*       more than 5 times, we decrease the precision of tracking
         NSTOUT=NSTOUT+1
         IF(MOD(NSTOUT,5).EQ.0) THEN
            EPSCUR=NSTOUT*EPSMAC
*            WRITE(CHMAIL,10000)ITRA,ISTAK,NTMULT,NAPART
*10000          FORMAT(' *** GTRACK *** Boundary loop: track ',
*     +         I4,' stack ',I4,' NTMULT ',I5,1X,5A4)
*            CALL GMAIL(1,0)
*            WRITE(CHMAIL,10100) EPSCUR
*10100          FORMAT('                Precision now set to ',G10.3)
*            CALL GMAIL(0,1)
         ENDIF
      ELSE
         NSTOUT = 0
         EPSCUR = EPSMAC
      ENDIF
*
      INWVOL = 1
*
* *** Compute SET and DET number if volume is sensitive
*
      IF (JSET.GT.0) CALL GFINDS
*
*    Clear step dependent variables
*
   80 NMEC   = 0
      STEP   = 0.
      DESTEL = 0.
      DESTEP = 0.
      NGKINE = 0
      IGNEXT = 0
      INWOLD = INWVOL
      PREC   = MAX(PRECOR,MAX(ABS(VECT(1)),ABS(VECT(2)),
     +                        ABS(VECT(3)),SLENG)*EPSCUR)
*
*     Give control to user at entrance of volume (INWVOL=1)
*
      IF (INWVOL.EQ.1) THEN
         CALL GDSTEP
         IF (ISTOP.NE.0) GO TO 22
         INWVOL = 0
      ENDIF
*
* *** Propagate particle up to next volume boundary or end of track
*
      INGOTO = 0
      NLEVIN = NLEVEL
      IF (IPARAM.NE.0) THEN
         IF (GEKIN.LE.PACUTS(ITRTYP)) THEN
            NMEC = NMEC+1
            LMEC(NMEC) = 26
            ISTOP = 2
#if !defined(CERNLIB_USRJMP)
            CALL GUPARA
#endif
#if defined(CERNLIB_USRJMP)
            CALL JUMPT0(JUPARA)
#endif
            GO TO 90
         ENDIF
      ENDIF
 
         CALL GDNINO
 
      STLOSS=STEP
*
*     Check for possible endless loop
*
   90 NSTEP = NSTEP +1
 
      IF (NSTEP.GT.MAXNST) THEN
         IF (ISTOP.EQ.0) THEN
            ISTOP = 99
            NMEC  = NMEC +1
            LMEC(NMEC) = 30
*            WRITE(CHMAIL,10200) MAXNST
*            CALL GMAIL(1,0)
*            WRITE(CHMAIL,10300)ITRA,ISTAK,NTMULT,(NAPART(I),I=1,5),
*     +      TOFG*1.E9
*            CALL GMAIL(0,1)
*10200       FORMAT(' *** GTRACK *** More than ',I6,
*     +             ' steps, tracking abandoned!')
*10300       FORMAT('                Track ',I4,' stack ',I4,' NTMULT ',
*     +             I5,1X,5A4,1X,'Time of flight ',F10.3,' ns')
         ENDIF
      ENDIF
*
* *** Give control to user at end of each tracking step
*
      SAFETY = SAFETY -STEP
      CALL GDSTEP
*
      IF (ISTOP.NE.0) GO TO 22
*
*      Renormalize direction cosines
*
      CMOD = 1./SQRT(VECT(4)*VECT(4)+VECT(5)*VECT(5)+VECT(6)*VECT(6))
      VECT(4) = VECT(4)*CMOD
      VECT(5) = VECT(5)*CMOD
      VECT(6) = VECT(6)*CMOD
*
      IF (INWVOL.EQ.0) GO TO 80
*
      IF (NJTMAX.GT.0) THEN
         CALL GSTRAC
         IF (NLEVIN.EQ.0) GO TO 100
         GO TO 22
      ELSE
         IF (NLEVIN.GE.NLEVEL) THEN
            INFROM = 0
         ELSE
            IF (NLEVIN.EQ.0) GO TO 100
            INFROM = LINDEX(NLEVIN+1)
         ENDIF
         IF (NLEVIN.NE.NLEVEL) INGOTO = 0
         NLEVEL = NLEVIN
*
         CALL GTMEDI (VECT, NUMED)
         IF (NUMED.NE.0) THEN
            SAFETY = 0.
            CALL UHTOC(NAMES(NLEVEL),4,NAME,4)
*            print *,NAME
           IF(IMYSE.EQ.1)THEN
            CALL GFIND(NAME,'SEEN',ISSEEN)
            CALL GFIND(NAME,'COLO',ISCOLO)
            CALL GFIND(NAME,'LSTY',ISLSTY)
           ENDIF
*            print *,isseen,iscolo,islsty
            GO TO 10
         ELSE
            ISSEEN=1
         ENDIF
      ENDIF
*
*     Track outside setup, give control to user (INWVOL=3)
*
  100 INWVOL = 3
      ISTOP  = 1
      ISET   = 0
      IDET   = 0
      NMEC   = 0
      STEP   = 0.
      DESTEL = 0.
      DESTEP = 0.
      NGKINE = 0
      NLCUR  = NLEVEL
      NLEVEL = 1
      CALL GDSTEP
      NLEVEL = NLCUR
 22        CONTINUE
          IF(IMAP.EQ.1)THEN
           IF(JON.EQ.1)THEN
             IF(VVV.LT.VVVMIN(IIIJ))THEN
               VVVMIN(IIIJ)=VVV
             ELSEIF(VVV.GT.VVVMAX(IIIJ))THEN
               VVVMAX(IIIJ)=VVV
             ENDIF
           ENDIF
          ENDIF
           VVV=VVV+FINVIS
           IYYY=IYYY+1
   23     CONTINUE
   21    CONTINUE
*                                                             END GTRACK
  999 END
