*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDRPRT(IPART,IMATE,STEP,MPOINT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Landau fluctuations and/or delta-rays generation .            *
C.    *  print the value of some relevant parameters (see phys 332) :  *
C.    *  Tmax, zeta, kapa and the regime for fluctuations              *
C.    *                                                                *
C.    *             Input parameters                                   *
C.    *  IPART   Geant particle number                                 *
C.    *  IMATE   Geant material number                                 *
C.    *  STEP    cm                                                    *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author     M.Maire    *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcmulo.inc"
*
      PARAMETER (MXPOIN=100)
      DIMENSION TCUT(5)
      DIMENSION TKINV(MXPOIN),DEDXMV(MXPOIN),DNTOTV(MXPOIN)
      CHARACTER REGIME*14,KUNIT*4,NAPART*8,NAMATE*8
*
      PARAMETER (DGEV=0.153536 E-3)
*
* *** incident particle
      JPA = LQ(JPART-IPART)
      CALL UHTOC (IQ(JPA+1),4,NAPART,8)
      AMASS = Q(JPA+7)
      CHARG = Q(JPA+8)
      EAMAS = EMASS/AMASS
*
* *** material
      JMA = LQ(JMATE-IMATE)
      IF(JMA.LE.0.OR.IMATE.GT.IQ(JMATE-2)) THEN
         WRITE(CHMAIL,10000) IMATE
         GOTO 999
      ENDIF
      CALL UHTOC (IQ(JMA+1),4,NAMATE,8)
      AMT = Q(JMA+6)
      ZMT = Q(JMA+7)
      RO  = Q(JMA+8)
      JPR = LQ(JMA-4)
*
* *** mean ionization potential (GeV)
*     POTI=16E-9*ZMT**0.9
      POTI=Q(JPR+9)
*
      CONS = DGEV*CHARG*CHARG*STEP*RO*ZMT/AMT
*
      WRITE (CHMAIL,10200) NAPART,NAMATE,STEP
      CALL GMAIL(1,0)
      WRITE (CHMAIL,10100)
      CALL GMAIL(0,1)
      WRITE (CHMAIL,10300)
      CALL GMAIL(0,1)
*
* *** loop on kinetic energy
      NPOINT = MIN(MPOINT,MXPOIN)
      TKINI = ELOW(1)
      TKEND = ELOW(NEK1)
      STEPL = (TKEND/TKINI)**(1./(NPOINT-1))
      TKINV(1) = TKINI
      DO 10 JTK=2,NPOINT
         TKINV(JTK) = TKINV(JTK-1)*STEPL
   10 CONTINUE
      CALL GFTMAT(IMATE,IPART,'LOSS',NPOINT,TKINV,DEDXMV,TCUT,IXST)
      CALL GFTMAT(IMATE,IPART,'DRAY',NPOINT,TKINV,DNTOTV,TCUT,IXST)
      DO 20  I=1,NPOINT
         TKIN  = TKINV(I)
         DEDXM = DEDXMV(I)
         DNTOT = DNTOTV(I)
         TKM = TKIN/AMASS
         GAMA = TKM + 1.
         GAM2 = GAMA**2
         BTA2 = TKM*(TKM+2.)/GAM2
*        maximum transferable energy (gev)
         TMAX = (2*EMASS*GAM2*BTA2)/(1.+2*GAMA*EAMAS+EAMAS**2)
         IF(IPART.EQ.3) TMAX=0.5*TMAX
*
* ***    mean de/dx on material
         DEGEV = STEP*DEDXM*1.E-3
         DEPOT = DEGEV/POTI
*
* ***    total number of delta ray production
         DNTOT = DNTOT*STEP
*
* ***    low energy tranfert
         ZETA = CONS/BTA2
         DNLOW = ZETA/POTI
*
* ***    high energy transfert
         DNHIG = ZETA/TMAX
*
* ***    regime
         DNMIN = MIN(DNLOW,DEPOT)
         IF(DNMIN.GE.50.) THEN
            IF (DNHIG.GE.10. ) THEN
               REGIME = 'Gauss'
            ELSEIF (DNHIG.GE.0.01) THEN
               REGIME = 'Vavilov'
            ELSE
               REGIME = 'Landau'
            ENDIF
         ELSE
            IFREE = 1
            IF(DNMIN.GE.0.01) THEN
               REGIME = 'Urban/'
               IFREE  = 7
            ENDIF
            IF(DNMIN.GE.1..AND.DNMIN.LT.30.) THEN
               REGIME(IFREE:) = 'ASHO/'
               IFREE  = 12
            ENDIF
            REGIME(IFREE:) = 'PAI'
         ENDIF
*
* ***    print
         CALL GEVKEV (TKIN,TKUN,KUNIT)
         WRITE (CHMAIL,10400) TKUN,KUNIT,DNTOT,DEPOT,DNLOW,DNHIG,
     +   REGIME
         CALL GMAIL(0,0)
   20 CONTINUE
*
10000 FORMAT(' **** GDRPRT: Material N ',I3,' has not been defined')
10100 FORMAT (9X,'For an explication of the ''regimes'' field ',
     +           'see the GEANT manual')
10200 FORMAT (5X,'*** Energy loss straggling for ',A8,' in ',A8,
     +           ', step ',F8.4,' cm')
10300 FORMAT (4X,'   Ekin ','     N d-rays','      dE/I    ',
     +    '    xi/I    ','   xi/Emax ','    regime')
10400 FORMAT (4X,F5.1,A4,4E12.3,4X,A14)
*
  999 END
