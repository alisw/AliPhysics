*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :          12/06/95  15.03.22  by  S.Ravndal
*-- Author :
      SUBROUTINE GXPHYS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Physics parameters control commands                       *
C.    *                                                                *
C.    *       Author:    R.Brun      **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmulo.inc"
      DIMENSION UCUTS(10),ULCUTS(10)
      EQUIVALENCE(UCUTS(1),CUTGAM)
      DIMENSION MECA(5,13)
      EQUIVALENCE (MECA(1,1),IPAIR)
      CHARACTER*6 CUTNAM(10)
      CHARACTER*4 CEN(10)
      CHARACTER*32 CHPATL
      CHARACTER*(*) CHNUMB
      PARAMETER (CHNUMB='1234567890')
      DATA CUTNAM/'CUTGAM','CUTELE','CUTNEU','CUTHAD','CUTMUO',
     +            'BCUTE' ,'BCUTM' ,'DCUTE' ,'DCUTM' ,'PPCUTM'/
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
      IF(CHPATL.EQ.'ANNI')THEN
         CALL KUGETI(IANNI)
*
      ELSEIF(CHPATL.EQ.'AUTO')THEN
         CALL KUGETI(IGAUTO)
*
      ELSEIF(CHPATL.EQ.'BREM')THEN
         CALL KUGETI(IBREM)
*
      ELSEIF(CHPATL.EQ.'CKOV')THEN
         CALL KUGETI(ICKOV)
*
      ELSEIF(CHPATL.EQ.'COMP')THEN
         CALL KUGETI(ICOMP)
*
      ELSEIF(CHPATL.EQ.'DCAY')THEN
         CALL KUGETI(IDCAY)
*
      ELSEIF(CHPATL.EQ.'DRAY')THEN
         CALL KUGETI(IDRAY)
*
      ELSEIF(CHPATL.EQ.'ERAN')THEN
         CALL KUGETR(EKMIN)
         CALL KUGETR(EKMAX)
         CALL KUGETI(NEKBIN)
         NEKBIN=MIN(NEKBIN,199)
*
      ELSEIF(CHPATL.EQ.'HADR')THEN
         CALL KUGETI(IHADR)
*
      ELSEIF(CHPATL.EQ.'LABS')THEN
         CALL KUGETI(ILABS)
*
      ELSEIF(CHPATL.EQ.'LOSS')THEN
         CALL KUGETI(ILOSS)
         IF(ILOSS.EQ.2.OR.ILOSS.EQ.0)THEN
            IDRAY=0
         ELSE
            IDRAY=1
         ENDIF
*
      ELSEIF(CHPATL.EQ.'MULS')THEN
         CALL KUGETI(IMULS)
*
      ELSEIF(CHPATL.EQ.'MUNU')THEN
         CALL KUGETI(IMUNU)
*
      ELSEIF(CHPATL.EQ.'PAIR')THEN
         CALL KUGETI(IPAIR)
*
      ELSEIF(CHPATL.EQ.'PFIS')THEN
         CALL KUGETI(IPFIS)
*
      ELSEIF(CHPATL.EQ.'PHOT')THEN
         CALL KUGETI(IPHOT)
*
      ELSEIF(CHPATL.EQ.'RAYL')THEN
         CALL KUGETI(IRAYL)
*
      ELSEIF(CHPATL.EQ.'STRA')THEN
         CALL KUGETI(ISTRA)
*
      ELSEIF(CHPATL.EQ.'SYNC')THEN
         CALL KUGETI(ISYNC)
*
      ELSEIF(CHPATL.EQ.'CUTS')THEN
         IF(NPAR.LE.0)THEN
            WRITE(LOUT,10000)
10000       FORMAT(/,' Current PHYSICS parameters:',/)
            DO 10 I=1,10
               CALL GEVKEV(UCUTS(I),ULCUTS(I),CEN(I))
               WRITE(LOUT,10100)CUTNAM(I),ULCUTS(I),CEN(I)
10100          FORMAT(5X,A,' = ',F7.2,1X,A)
   10       CONTINUE
            GO TO 999
         ENDIF
         CALL KUGETR(CUTGAM)
         CALL KUGETR(CUTELE)
         CALL KUGETR(CUTHAD)
         CALL KUGETR(CUTNEU)
         CALL KUGETR(CUTMUO)
         CALL KUGETR(BCUTE)
         CALL KUGETR(BCUTM)
         CALL KUGETR(DCUTE)
         CALL KUGETR(DCUTM)
         CALL KUGETR(PPCUTM)
         CALL KUGETR(TOFMAX)
         CALL KUGETR(GCUTS(1))
         IF(BCUTE.LE.0.)BCUTE=CUTGAM
         IF(BCUTM.LE.0.)BCUTM=CUTGAM
         IF(DCUTE.LE.0.)DCUTE=CUTELE
         IF(DCUTM.LE.0.)DCUTM=CUTELE
         IF(PPCUTM.LT.4.*EMASS)PPCUTM=4.*EMASS
*
      ELSEIF(CHPATL.EQ.'DRPRT')THEN
         CALL KUGETI(IPART)
         CALL KUGETI(IMATE)
         CALL KUGETR(STEP)
         CALL KUGETI(NPOINT)
         CALL GDRPRT(IPART,IMATE,STEP,NPOINT)
*
      ELSEIF(CHPATL.EQ.'PHYSI')THEN
         IF(JTMED.GT.0)THEN
            DO 30 I=1,IQ(JTMED-2)
               JTM=LQ(JTMED-I)
               IF(JTM.LE.0)GO TO 30
               IF(IQ(JTM-2).EQ.0)THEN
                  CALL MZPUSH(IXCONS,JTM,10,0,'I')
                  GO TO 30
               ENDIF
               DO 20 J=1,10
                  JTMI=LQ(JTM-J)
                  IF(JTMI.GT.0)THEN
                     CALL MZDROP(IXCONS,JTMI,' ')
                  ENDIF
   20          CONTINUE
   30       CONTINUE
            CALL UCOPY(CUTGAM,Q(JTMED+1),10)
            DO 40 I=1,13
               Q(JTMED+10+I)=MECA(1,I)
   40       CONTINUE
         ENDIF
         IF(JMATE.LE.0)GO TO 999
         DO 60 I=1,IQ(JMATE-2)
            JMA=LQ(JMATE-I)
            IF(JMA.LE.0)GO TO 60
            DO 50 J=1,IQ(JMA-2)
               IF(J.EQ.4.OR.J.EQ.5)GO TO 60
               JM=LQ(JMA-J)
               IF(JM.LE.0)GO TO 50
               CALL MZDROP(IXCONS,JM,'L')
   50       CONTINUE
   60    CONTINUE
         CALL MZGARB (IXCONS, 0)
         CALL GPHYSI
      ENDIF
*
  999 END
