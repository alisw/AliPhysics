*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.43  by  S.Giani
*-- Author :
*$ CREATE HADEVV.FOR
*COPY HADEVV
*
*=== hadevv ===========================================================*
*
      SUBROUTINE HADEVV ( NHAD, KPROJ, KTARG, PPROJ, EPROJ, UMO )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*    Modified version of Hadevt created by Alfredo Ferrari, INFN-Milan *
*                                                                      *
*    Last change  on  20-jun-93  by  Alfredo Ferrari, INFN - MIlan     *
*                                                                      *
*  Hadevt90: kinematics completed revised by A. Ferrari, before it was *
*            always wrong every time the second jet to be sampled was  *
*            a "parjet". A few other bugs corrected: maybe others are  *
*            still in!!!                                               *
*----------------------------------------------------------------------*
*
C
C     GENERATE HADRON PRODC
C     GENERATE HADRON PRODUCTION EVENT IN  KPROJ - KTARG  COLLISION
C     WITH  LAB PROJECTILE MOMENTUM  PPROJ
C     INCLUDING MESON MESON AND MESON ANTIBARYON COLLISIONS
C
C********************************************************************
C
#include "geant321/auxpar.inc"
#include "geant321/balanc.inc"
#include "geant321/cmsres.inc"
#include "geant321/finpar.inc"
#include "geant321/hadpar.inc"
#include "geant321/inpdat2.inc"
#include "geant321/part.inc"
#include "geant321/qquark.inc"
      COMMON /FKINVT/PNUC(3),INUCVT
      COMMON /FKPRIN/ IPRI, INIT
      REAL RNDM(2)
      LOGICAL LISSU, LQTARG, LQPROJ
*
      SAVE UNON,UNOM,UNOMS
      DATA UNON/2.0D0/
      DATA UNOM/1.5D0/
      DATA UNOMS/0.5D0/
*
C
C*******************************************************************"
C
C     KINEMATICS
C
C********************************************************************
C
 8899 CONTINUE
* Ijproj = paprop numbering
      IJPROJ = KPTOIP (KPROJ)
      IJTARG = KPTOIP (KTARG)
      AMPROJ = AM(KPROJ)
      AMTAR = AM(KTARG)
*  The usual gamma and sqrt[beta**2/(1-beta**2)]=eta=gamma*beta factors
*  for the CMS system
*
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM  = PPROJ / UMO
C
*or      IF(IPRI.EQ.1) WRITE(LUNOUT,101)KPROJ,KTARG,PPROJ,AMPROJ,AMTAR,
*or     &EPROJ,UMO,GAMCM,BGCM
*or 101  FORMAT(2I5,10F11.5)
C
C********************************************************************
C
C     SELECTION OF  QUARK - DIQUARK - CHAINS
C
C********************************************************************
C
*
*  Ibproj = baryonic charge of the projectile
*
      IBPROJ = IBAR(KPROJ)
*
*  Ibtarg = baryonic charge of the target nucleon
*
      IBTARG = IBAR(KTARG)
*
*  Ipq1,ipq2,ipq3 = quarks of the projectile
*
      IQP1 = MQUARK(1,IJPROJ)
      IQP2 = MQUARK(2,IJPROJ)
      IQP3 = MQUARK(3,IJPROJ)
*
*  Iqt1,iqt2,iqt3 = quarks of the projectile
*
      IQT1 = MQUARK(1,IJTARG)
      IQT2 = MQUARK(2,IJTARG)
      IQT3 = MQUARK(3,IJTARG)
*or      IF (IPRI .EQ. 1)
*or     &WRITE(LUNOUT,102)IBPROJ, IQP1,IQP2,IQP3,IQT1,IQT2,
*or     &IQT3
*or 102  FORMAT(12I10)
      IF (IBPROJ) 11, 12, 13
  11  CONTINUE
C
C********************************************************************
C
C     SELECTION OF CHAINS
C     ANTIBARYON - BARYON   COLLISION
C
C********************************************************************
C
*  +-------------------------------------------------------------------*
*  | The incoming projectile is an antibaryon!!!
*  |
         CALL GRNDM(RNDM,1)
         ISAM1 = 1.D0 + 3.D0*RNDM(1)
         GO TO (111,112,113),ISAM1
 111     CONTINUE
            IBF = IQP1
            IFF1 = IQP2
            IFF2 = IQP3
            GO TO 114
 112     CONTINUE
            IBF = IQP2
            IFF1 = IQP1
            IFF2 = IQP3
            GO TO 114
 113     CONTINUE
            IBF = IQP3
            IFF1 = IQP1
            IFF2 = IQP2
 114     CONTINUE
         CALL GRNDM(RNDM,1)
         ISAM2 = 1.D0 + 3.D0*RNDM(1)
         GO TO (115,116,117),ISAM2
 115     CONTINUE
            IBB = IQT1
            IFB1 = IQT2
            IFB2 = IQT3
            GO TO 118
 116     CONTINUE
            IBB = IQT2
            IFB1 = IQT1
            IFB2 = IQT3
            GO TO 118
 117     CONTINUE
            IBB = IQT3
            IFB1 = IQT1
            IFB2 = IQT2
 118     CONTINUE
         GO TO 14
*  | Quark selection for incoming antibaryon has been completed
*  +-->-->-->-->-->-->-->-->-->--> go to 14 continue
 
  12  CONTINUE
*  +-------------------------------------------------------------------*
*  | The incoming projectile is a meson!!!
*  |
         IF (IBTARG)712,812,912
*  |  +----------------------------------------------------------------*
*  |  | The target nucleon is a baryon (meson projectile)
*  |  |
 912        CONTINUE
C
C********************************************************************
C
C     SELECTION OF CHAINS
C     MESON - BARYON  COLLISION
C
C********************************************************************
C
            CALL GRNDM(RNDM,1)
            ISAM3 = 1.D0 + 2.D0*RNDM(1)
            GO TO (121,122),ISAM3
 121        CONTINUE
               IFF = IQP1
               IBF = IQP2
               GO TO 123
 122        CONTINUE
               IFF = IQP2
               IBF = IQP1
 123        CONTINUE
            CALL GRNDM(RNDM,1)
            ISAM4 = 1.D0 + 3.D0*RNDM(1)
            GO TO (124,125,126),ISAM4
 124        CONTINUE
               GO TO (1241,1242),ISAM3
1241           CONTINUE
                  IBB = IQT1
                  IFB1 = IQT2
                  IFB2 = IQT3
                  GO TO 127
1242           CONTINUE
                  IBB1 = IQT2
                  IBB2 = IQT3
                  IFB = IQT1
                  GO TO 127
 125        CONTINUE
               GO TO (1251,1252),ISAM3
1251           CONTINUE
                  IBB = IQT2
                  IFB1 = IQT1
                  IFB2 = IQT3
                  GO TO 127
1252           CONTINUE
                  IBB1 = IQT1
                  IBB2 = IQT3
                  IFB = IQT2
                  GO TO 127
 126        CONTINUE
               GO TO (1261,1262),ISAM3
1261           CONTINUE
                  IBB = IQT3
                  IFB1 = IQT1
                  IFB2 = IQT2
                  GO TO 127
1262           CONTINUE
                  IBB1 = IQT1
                  IBB2 = IQT2
                  IFB = IQT3
 127        CONTINUE
            GO TO 14
*  |  | Quark selection for incoming meson and baryon target completed
*  |  +-->-->-->-->-->-->-->-->-->--> go to 114 continue
 
*  |  +----------------------------------------------------------------*
*  |  | The target nucleon is a meson (meson projectile)
*  |  |
 812        CONTINUE
C===============================================================
C
C     SELECTION OF CHAINS
C        MESON MESON COLLISIONS
C
C================================================================
            CALL GRNDM(RNDM,1)
            ISAM3 = 1.D0 + 2.D0*RNDM(1)
            GO TO (1218,1228),ISAM3
1218        CONTINUE
               IFF = IQP1
               IBF = IQP2
               IBB = IQT1
               IFB = IQT2
               GO TO 1238
1228        CONTINUE
               IFF = IQP2
               IBF = IQP1
               IBB = IQT2
               IFB = IQT1
1238        CONTINUE
            GO TO 14
*  |  | Quark selection for incoming meson and meson target completed
*  |  +-->-->-->-->-->-->-->-->-->--> go to 14 continue
 
*  |  +----------------------------------------------------------------*
*  |  | The target nucleon is an anti-baryon (meson projectile)
*  |  |
 712        CONTINUE
C=================================================================
C
C     SELECTION OF CHAINS
C      MESON ANTIBARYON COLLISIONS
C
C==================================================================
            ISAM3 = 2
            IFF = IQP1
            IBF = IQP2
            CALL GRNDM(RNDM,1)
            ISAM4 = 1.D0 + 3.D0*RNDM(1)
            GO TO (1247,1257,1267),ISAM4
1247        CONTINUE
               GO TO (12417,12427),ISAM3
12417          CONTINUE
                  IBB = IQT1
                  IFB1 = IQT2
                  IFB2 = IQT3
                  GO TO 1277
12427          CONTINUE
                  IBB1 = IQT2
                  IBB2 = IQT3
                  IFB = IQT1
                  GO TO 1277
1257        CONTINUE
               GO TO (12517,12527),ISAM3
12517          CONTINUE
                  IBB = IQT2
                  IFB1 = IQT1
                  IFB2 = IQT3
                  GO TO 1277
12527          CONTINUE
                  IBB1 = IQT1
                  IBB2 = IQT3
                  IFB = IQT2
                  GO TO 1277
1267        CONTINUE
               GO TO (12617,12627),ISAM3
12617          CONTINUE
                  IBB = IQT3
                  IFB1 = IQT1
                  IFB2 = IQT2
                  GO TO 1277
12627          CONTINUE
                  IBB1 = IQT1
                  IBB2 = IQT2
                  IFB = IQT3
1277        CONTINUE
            GO TO 14
*  |  | Quark selection for incoming meson and a-baryon target completed
*  |  +-->-->-->-->-->-->-->-->-->--> go to 14 continue
*  |
*  |                                     end meson projectile
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |  The incoming projectile is a baryon!!!
*  |
  13  CONTINUE
C
C********************************************************************
C
C     SELECTION OF CHAINS
C     BARYON - BARYON   COLLISION
C
C********************************************************************
C
         CALL GRNDM(RNDM,1)
         ISAM5 = 1.D0 + 3.D0*RNDM(1)
         GO TO (131,132,133),ISAM5
 131     CONTINUE
            IBF = IQP1
            IFF1 = IQP2
            IFF2 = IQP3
            GO TO 134
 132     CONTINUE
            IBF = IQP2
            IFF1 = IQP1
            IFF2 = IQP3
            GO TO 134
 133     CONTINUE
            IBF = IQP3
            IFF1 = IQP1
            IFF2 = IQP2
 134     CONTINUE
         CALL GRNDM(RNDM,1)
         ISAM6 = 1.D0 + 3.D0*RNDM(1)
         GO TO (135,136,137),ISAM6
 135     CONTINUE
            IFB = IQT1
            IBB1 = IQT2
            IBB2 = IQT3
            GO TO 138
 136     CONTINUE
            IFB = IQT2
            IBB1 = IQT1
            IBB2 = IQT3
            GO TO 138
 137     CONTINUE
            IFB = IQT3
            IBB1 = IQT1
            IBB2 = IQT2
 138     CONTINUE
*  |  | Quark selection for incoming baryon and baryon target completed
*  +  |-->-->-->-->-->-->-->-->-->--> go to 14 continue
  14  CONTINUE
*  | Quark selection completed
*  +-------------------------------------------------------------------*
 
*or      IF (IPRI.EQ.1) WRITE(LUNOUT,102)IFF,IBF,IFF1,IFF2,IFB1,IFB2,
*or     &IFB,IBB,IBB1,IBB2
C
C********************************************************************
C
C*** SAMPLING MOMENTUM FRACTIONS OF QUARKS AND DIQUARKS
C
C********************************************************************
C
      IXPXT = 0
  25  CONTINUE
      IXPXT = IXPXT + 1
*  +-------------------------------------------------------------------*
*  | Selection of xp and xt from beta distribution:
*  |   xp and xt are then used to select the fraction of momentum and
*  |   energy for each jet, according to the following relations:
*  |   (where we assume to use xp, xt for the jet n. 1)
*  |
*  |      xxp = 1 - xp
*  |      xxt = 1 - xt
*  |      Ech1 = umo * (xp  + xt ) / 2
*  |      Ech2 = umo * (xxp + xxt) / 2
*  |      Pch1 = umo * (xp  - xt ) / 2
*  |      Pch2 = umo * (xxp - xxt) / 2
*  |      Amch1 = umo * sqrt (xp  * xt )
*  |      Amch2 = umo * sqrt (xxp * xxt)
*  |
      IF (IBPROJ) 21,22,23
  21  CONTINUE
*  |  Note for antibaryon projectile xp and xt are sampled from the
*  |  same distribution, ===> no difference in exchanging them!
         UNO = UNON
         XP  = BETARN(HLFHLF,UNO)
         XXP = 1.D0 - XP
         UNO = UNON
         XT  = BETARN(HLFHLF,UNO)
         XXT = 1.D0 - XT
         GO TO 24
  23  CONTINUE
*  |  Note for baryon projectile xp and xt are sampled from the
*  |  same distribution, ===> no difference in exchanging them!
         UNO = UNON
         XP  = BETARN(HLFHLF,UNO)
         XXP = 1.D0 - XP
         UNO = UNON
         XT  = BETARN(HLFHLF,UNO)
         XXT = 1.D0 - XT
         GO TO 24
  22  CONTINUE
*  |  Note for meson projectile xp and xt are not sampled from the
*  |  same distribution, ===> difference in exchanging them!
         UNO = UNOM
         IF (IFF.EQ.3 .OR. IFF.EQ.-3) UNO = UNOMS
         XP  = BETARN(HLFHLF,UNO)
         XXP = 1.D0 - XP
         IF (IBTARG .EQ. 0) GO TO 2288
            UNO = UNON
2288     CONTINUE
         XT  = BETARN(HLFHLF,UNO)
         XXT = 1.D0 - XT
  24  CONTINUE
*  |  Now xp and xt have been selected from the appropriate beta
*  |  distributions, xxp and xxt are the complements to one of xp and xt
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |  From here to 1124 it is likely to be obsolete (inucvt is now used
*  |  nowhere in the code, the same for pnuc
*  |
*or      RNDMVV=RNDM(V)
*or      IF (INUCVT.EQ.0) GO TO 1124
*or      IF (RNDMVV.LT.PNUC(1)) GO TO 1124
*or      XT=2.D0*BETARN(0.5D0,UNO+6.D0)
*or      XXT=1.D0-XT
*or      IF (XXT.LE.0.D0) XXT=RNDM(V)
*or      IF (RNDMVV.LT.PNUC(2)) GO TO 1124
*or      XT=3.D0*BETARN(0.5D0,UNO+12.D0)
*or      XXT=1.D0-XT
*or      IF (XXT.LE.0.D0) XXT=RNDM(V)
*or 1124 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*or      IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XP,XT,XXP,XXT
*or 103  FORMAT (10F10.5)
C
C********************************************************************
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C
C********************************************************************
C
****
****===================================================================*
*  |  Now selecting the kinematical parameters for the two jets:
*  |
*  |      amch1 = invariant mass of the 1st jet
*  |      ech1  = total energy of the 1st jet in CMS
*  |      pch1  = total momentum of the 1st jet in CMS (with sign)
*  |
*  |      amch2 = invariant mass of the 2nd jet
*  |      ech2  = total energy of the 2nd jet in CMS
*  |      pch2  = total momentum of the 2nd jet in CMS (with sign)
*  |
*  |  The following relations must be fulfilled:
*  |
*  |      ech1 + ech2 = umo (energy in CMS = inv. mass of the system)
*  |      ech1 = sqrt (pch1**2 + amch1**2)
*  |      ech2 = sqrt (pch2**2 + amch2**2)
*  |      pch1 + pch2 = 0
*  |
****===================================================================*
****
      IF (IBPROJ) 31,32,33
*  +-------------------------------------------------------------------*
*  |   antibaryon projectile!!! (note that only antibaryon projectile
*  |   baryon target is allowed, antibaryon projectile meson target is
*  |   considered as meson projectile antibaryon target)
*  |
  31  CONTINUE
C
C********************************************************************
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C*** ANTINUCLEON-NUCLEUON
C*** LONG ANTIDIQUARK - DIQUARK CHAIN
C
C********************************************************************
C
*  |  Ibb, ifb1, ifb2 contain the quark numbers of the target, ibf,
*  |  iff1, iff2 the quark numbers of the projectile
*  |  iff...= forward chain, forward quark (diquark)
*  |  ifb...= forward chain, backward quark (diquark)
*  |  ibf...= backward chain, forward quark
*  |  ibb...= backward chain, backward quark
*  |  By definition all i..f.. come from the projectile and all
*  |  i..b.. from the target
*  |  Of course, since we are treating an antibaryon projectile and
*  |  a baryon target all i..f.. are antiquark and all i..b.. are
*  |  quark
*  |  Of course the two following cards are equivalent to
*  |     IIFF1 = IABS (IFF1) ...
         IIFF1 = -IFF1
         IIFF2 = -IFF2
         IF (IIFF1 .EQ. IFB1) GO TO 3111
         IF (IIFF1 .EQ. IFB2) GO TO 3112
         IF (IIFF2 .EQ. IFB1) GO TO 3113
         IF (IIFF2 .EQ. IFB2) GO TO 3114
*  |    Get the index and the mass of the pseudoscalar meson corre-
*  |    sponding to the lowest energy for chain 2 ("b")
         IIBF = IABS(IBF)
         IBPS = IMPS(IIBF,IBB)
         AMBPS = AM(IBPS)
*  |  *****************************************************************
*  |    New version: of course, as it is explained below it is not    *
*  |    correct to comment the "go to 3115", however it is not correct*
*  |    also to use it in its original form, we need to compute a de- *
*  |    tailed threshold, also for the Amff value: we can believe that*
*  |    the lowest threshold is given by the two scalar mesons resul- *
*  |    ting from the combinations of the 4 quarks (iff1,iff2,ibf1,   *
*  |    ibf2).  Remember that the Imps(i,j) array gives the index of  *
*  |    the pseudoscalar meson with antiquark -i and quark j, the same*
*  |    apply to the Imve array but for vector mesons                 *
*  |    But, since bamjev it is likely to produce at least one baryon *
*  |    and one antibaryon when called with Iopt=5 (since it hadroni- *
*  |    zes a chain with a diquark and an anti-diquark at the extremi-*
*  |    ties) a more realistic threshold could be to check for the    *
*  |    masses of the baryon-antibaryon combinations corresponding to *
*  |    a uubar or a ddbar sea pair added to the original diquarks    *
*  |  *****************************************************************
*  |  Selection of the mass threshold from the two pseudoscalar mesons
*        IMPS11 = IMPS(IIFF1,IFB1)
*        IMPS21 = IMPS(IIFF2,IFB2)
*        IMPS12 = IMPS(IIFF1,IFB2)
*        IMPS22 = IMPS(IIFF2,IFB1)
*  |  Amff is selected as the maximum of the two possible meson configu-
*  |  rations, to be sure that no problem will result during frag-
*  |  mentation
*        AMFF = MAX ( AM (IMPS11) + AM (IMPS21), AM (IMPS12) +
*    &                AM (IMPS22) )
*  |  Of course at least two mesons must be produced
*  |  First check that the total invariant mass is enough (it must
*  |  be larger than the two meson masses and the pseudoscalar
*  |  meson mass of the second chain)
*        AMFF = AMFF
*  |  Selection of the mass threshold from the two baryon configura-
*  |  tions
         CALL BKLASS (-1, IFF1, IFF2, IA1F8, IA1F10 )
         CALL BKLASS ( 1, IFB1, IFB2,  I1F8,  I1F10 )
         CALL BKLASS (-2, IFF1, IFF2, IA2F8, IA2F10 )
         CALL BKLASS ( 2, IFB1, IFB2,  I2F8,  I2F10 )
         AMFF = MIN ( AM (IA1F8) + AM (I1F8), AM (IA2F8) + AM (I2F8) )
     &        + 0.3D+00
*  |  +----------------------------------------------------------------*
*  |  |  New treatment: check the mass threshold
         IF ( AMFF + AMBPS .LT. UMO ) THEN
            XSQ   = SQRT(XXP*XXT)
            AMCH1 = UMO*XSQ
            NNCH1 = 0
            AMCH2 = UMO*UMO*XP*XT
*  |  |  +------------------------------------------------------------*
*  |  |  |  Check if we have enough energy for the "f" jet
            IF ( AMCH1 .LE. AMFF .OR. AMCH2 .LE. AMBPS * AMBPS ) THEN
               IXPXT = IXPXT + 1
               IF ( IXPXT .LT. 5 ) GO TO 25
*  |  |  |   if amch1 < amfps xp and xt are resampled, but if we are
*  |  |  |   resampling too often force amch1 to be above
*  |  |  |   the minimum required, anyway the kinematic region
*  |  |  |   allowed for xp, xt seems to be marginal
*  |  |  *-->-->-->-->-->-->-->-->-->--> xp, xt resampling
               XSQ1 = ( AMFF  / UMO )**2
               XSQ2 = ( AMBPS / UMO )**2
               XXXMN = XSQ1
               XXYMX = 0.5D+00 * ( 1.D+00 + XSQ1 - XSQ2 )
               XXYMN = SQRT ( 1.D+00 - XSQ1 / ( XXYMX * XXYMX ) )
               XXXMN = MAX  ( XSQ1, XXYMX * ( ONEONE - XXYMN ) )
               XXXMX = MIN  ( 1.D+00, XXYMX * ( 1.D+00 + XXYMN ) )
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
 3161          CONTINUE
                  CALL GRNDM(RNDM,1)
                  XXP = XXXMN + ( XXXMX - XXXMN ) * RNDM (1)
                  XP  = 1.D+00 - XXP
                  XXYMN = XSQ1 / XXP
                  XXYMX = 1.D+00 - XSQ2 / XP
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( XXYMN .GT. XXYMX ) THEN
                     XXXMN = XXP
                     GO TO 3161
*  |  |  |  |-<|--<--<--<--< no allowed interval for xxt, resample
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               CALL GRNDM(RNDM,1)
               XXT = XXYMN + ( XXYMX - XXYMN ) * RNDM (1)
               XT  = 1.D+00 - XXT
               XSQ   = SQRT(XXP*XXT)
               AMCH1 = UMO*XSQ
               NNCH1 = 0
            END IF
*  |  |  |
*  |  |  +------------------------------------------------------------*
            IOPBAM = 5
            GO TO 3116
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |           We are in troubles: the selected quark combinations
*  |  |           for the two chains are unphysical since the invariant
*  |  |           mass is too low to produce the required three mesons
*  |  |           First try to change the quark combinations:
         ELSE
            CALL GRNDM(RNDM,1)
            IRNDM = 1.D+00 + RNDM (1)
            LQPROJ = .FALSE.
            LQTARG = .FALSE.
            GO TO (3171,3181) IRNDM
*  |  |  +-------------------------------------------------------------*
*  |  |  |        Try to change one of the projectile quarks in the
*  |  |  |        first chain
 3171       CONTINUE
               LQPROJ = .TRUE.
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |     The third antibaryon quark can combine with the
*  |  |  |  |     first or/and the second quark of the target diquark
               IF ( -IBF .EQ. IFB1 .OR. -IBF .EQ. IFB2 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  Make a random choiche of the quark to be substituted
                  CALL GRNDM(RNDM,1)
                  IF ( RNDM (1) .LT. 0.5D+00 ) THEN
                     IBF0 = IBF
                     IBF  = IFF1
                     IFF1 = IBF0
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     IBF0 = IBF
                     IBF  = IFF2
                     IFF2 = IBF0
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  GO TO 31
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  No possibility to solve the situation changing one of
*  |  |  |  |  the projectile quarks inside the "f" chain, try with
*  |  |  |  |  the projectile quarks (if not yet tried)
               ELSE
                  IF ( .NOT. LQTARG ) GO TO 3181
                  GO TO 3191
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |        Try to change one of the target quarks in the
*  |  |  |        first chain (for uud and udd targets this is usually
*  |  |  |        useless, but it has been included for the sake of
*  |  |  |        completness )
 3181       CONTINUE
               LQTARG = .TRUE.
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |     The third target quark can combine with the first
*  |  |  |  |     or/and the second quark of the projectile diquark
               IF ( IBB .EQ. -IFF1 .OR. IBB .EQ. -IFF2 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  Make a random choiche of the quark to be substituted
                  CALL GRNDM(RNDM,1)
                  IF ( RNDM (1) .LT. 0.5D+00 ) THEN
                     IBB0 = IBB
                     IBB  = IFB1
                     IFB1 = IBB0
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     IBB0 = IBB
                     IBB  = IFB2
                     IFB2 = IBB0
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  GO TO 31
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  No possibility to solve the situation changing one of
*  |  |  |  |  the target quarks inside the "f" chain, try with
*  |  |  |  |  the projectile quarks (if not yet tried)
               ELSE
                  IF ( .NOT. LQPROJ ) GO TO 3171
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 3191       CONTINUE
*  |  |  If we are here we cannot perform an interaction conserving
*  |  |  all additive quantum numbers: ignore one (typically it is
*  |  |  strangeness) and go on
            WRITE (LUNERR,*)
     &      ' *** Hadevv, impossible interaction, kp,kt, Umo',
     &        KPROJ,KTARG,UMO
            WRITE (LUNOUT,*)
     &      ' *** Hadevv, impossible interaction, kp,kt, Umo',
     &        KPROJ,KTARG,UMO
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  Now we want to get the indexes of the pseudo-scalar and vector
*  |  mesons which can be created from the first (forward) chain.
*  |  Of course if one of the ff quark is the antiquark of one of
*  |  fb quarks then the mesons are defined completely by the remaining
*  |  two quarks: we have anyway to deal also with the situation
*  |  were no ff quark is the antiquark of a fb quark (the goto 3115...)
*  |  which clearly requires at least two mesons!!! So no Amfps and
*  |  Amfv can be defined in this situation, and the "go to 3115" was
*  |  really important to assure that we are not producing single
*  |  particle jets which are not possible, we have only to substitute
*  |  the 1.5 threshold with the proper threshold for the two
*  |  mesons (at least the two scalar mesons): this is also true for
*  |  the Amff value which of course should be equal to the sum
*  |  of the two scalar mesons resulting from the combinations
*  |  of the 4 quarks (see above)
3111     CONTINUE
            IIFF2 = IABS(IFF2)
            IFPS  = IMPS(IIFF2,IFB2)
            IAIFF = IFF2
            IFB   = IFB2
            IFPS2 = IMPS(IIFF1,IFB1)
            IFV   = IMVE(IIFF2,IFB2)
            GO TO 3117
3112     CONTINUE
            IIFF2 = IABS(IFF2)
            IFPS  = IMPS(IIFF2,IFB1)
            IAIFF = IFF2
            IFB   = IFB1
            IFPS2 = IMPS(IIFF1,IFB2)
            IFV   = IMVE(IIFF2,IFB1)
            GO TO 3117
3113     CONTINUE
            IIFF1 = IABS(IFF1)
            IFPS  = IMPS(IIFF1,IFB2)
            IAIFF = IFF1
            IFB   = IFB2
            IFPS2 = IMPS(IIFF2,IFB1)
            IFV   = IMVE(IIFF1,IFB2)
            GO TO 3117
3114     CONTINUE
            IIFF1 = IABS(IIFF1)
            IFPS  = IMPS(IIFF1,IFB1)
            IAIFF = IFF1
            IFB   = IFB1
            IFPS2 = IMPS(IIFF2,IFB2)
            IFV   = IMVE(IIFF1,IFB1)
3117     CONTINUE
*  Amfps, amfv are the masses of the pseudoscalar and vector mesons
*  corresponding to the two unpaired quarks of the 1 (f) chain
         AMFPS = AM(IFPS)
         AMFV  = AM(IFV )
         AMFPS2= AM(IFPS2)
         NNCH1 = 0
C     ATTENTION THIS MIGHT LEAD TO TOO LOW ANNIHILATION MULTIPLICITIES
C        AMFF = AMFV+0.3D0
*        AMFF = 2.3D0
         AMFF  = MAX ( AMFV + 0.3D+00, AMFV + AMFPS2 + 0.1D+00 )
*  | This expression for Amff0 (threshold for a complete diquark-diquark
*  | jet) is not correct: this jet will fragment at least into two
*  | hadrons corresponding to the lowest energy hadrons with respecti-
*  | vely the first and the second diquark. Of corse one will be an
*  | antibaryon. It must be checked that bamjet produces a baryon-
*  | antibaryon pair any time it is called with iopt = 5 (corresponding
*  | to the fragmentation of a complete diquark-antidiquark jet)
*  | It is also questionable if it is correct to call bamjet with
*  | iopt=5 whenever energetically possible: this clearly implies
*  | that no annihilation takes place (an antibaryon will emerge
*  | from the interaction). This is a question to be addressed
*  | to Hannes!!
*        AMFF0 = 2.D+00 * AMFF
*  | A tentative way could be to check for the masses of the
*  | baryon-antibaryon combinations corresponding to a uubar
*  | or a ddbar sea pair
         CALL BKLASS (-1, IFF1, IFF2, IA1F8, IA1F10 )
         CALL BKLASS ( 1, IFB1, IFB2,  I1F8,  I1F10 )
         CALL BKLASS (-2, IFF1, IFF2, IA2F8, IA2F10 )
         CALL BKLASS ( 2, IFB1, IFB2,  I2F8,  I2F10 )
         AMFF0 = MIN ( AM (IA1F8) + AM (I1F8), AM (IA2F8) + AM (I2F8) )
     &         + 0.3D+00
*  |
*  |  here xxt and xxp are used for the first jet
*  |
         XSQ  = SQRT(XXP*XXT)
         AMCH1 = UMO*XSQ
         AAPS  = IFPS
         AAV  = IFV
*or         IF (IPRI.EQ.1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2,AAPS,AAV,
*or     &                                   AMFPS,AMFV
 
*  |  +----------------------------------------------------------------*
*  |  |
         IF (AMCH1 .LT. AMFPS) GO TO 25
*  |  |   if amch1 < amfps xp and xt are resampled
*  |  *-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
*  |  Kinematical parameters xxp, xxt for the 1st jet
         IF (AMCH1 .GT. AMFV ) GO TO 3151
C*** PRODUCE AMFPS
            AMCH1 = AMFPS
            NNCH1 = -1
            XSQ = AMFPS/UMO
*  |  |  recalculating xxp, xp
            XXP = XSQ**2/XXT
            XP  = 1.D0 - XXP
            GO TO 3153
3151     CONTINUE
         IF (AMCH1 .GT. AMFF) GO TO 3153
C*** PRODUCE AMFV
            AMCH1 = AMFV
            NNCH1 = 1
            XSQ = AMFV/UMO
*  |  |  recalculating xxp, xp
            XXP = XSQ**2/XXT
            XP  = 1.D0 - XXP
            GO TO 3153
3153     CONTINUE
         IF (AMCH1 .LE. AMFF0) THEN
            IOPBAM = 3
         ELSE
            IOPBAM = 5
         END IF
         GO TO 3116
3116     CONTINUE
C
C********************************************************************
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C     ANTINUCLEON NUCLEON
C*** SHORT ANTIQUARK - QUARK CHAIN
C
C********************************************************************
C
         IIBF = IABS(IBF)
         IBPS = IMPS(IIBF,IBB)
         IBV  = IMVE(IIBF,IBB)
         AMBPS = AM(IBPS)
         AMBV  = AM(IBV)
         NNCH2 = 0
         AMBB = AMBV + 0.3D0
         AAPS = IBPS
         AAV  = IBV
*  |  +----------------------------------------------------------------*1
*  |  |
*  |  | Now commented, is useless!!!
*         IF (XP .LE. 0.D0 .OR. XT .LE. 0.D0) GO TO 25
*  |  | resample xp and xt if one is negative or zero
*  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
*  |  Now xp and xt are used for the second jet
         XXSQ  = SQRT(XP*XT)
         AMCH2 = UMO*XXSQ
*or         IF (IPRI.EQ.1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &                  ,AAPS,AAV,AMBPS,AMBV
 
*  |  +----------------------------------------------------------------*
*  |  |
         IF (AMCH2 .LT. AMBPS) THEN
            IXPXT = IXPXT + 1
            GO TO 25
*  |  |   if amch1 < amfps xp and xt are resampled
*  |  *-->-->-->-->-->-->-->-->-->--> xp, xt resampling
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 
         IF (AMCH2 .GT. AMBV ) GO TO 3121
C*** PRODUCE AMBPS
* For Prof. Ranft: here there was a "large" mistake, amch2 = ambps
* was missing
            AMCH2 = AMBPS
            NNCH2 = -1
            XXSQ  = AMBPS/UMO
*or            IF (INUCVT .EQ. 1) GO TO 3123
            GO TO 31236
3121     CONTINUE
         IF (AMCH2 .GT. AMBB) GO TO 3123
C*** PRODUCE AMBV
            AMCH2 = AMBV
            NNCH2 = 1
            XXSQ  = AMBV/UMO
*or            IF (INUCVT .EQ. 1) GO TO 3123
            GO TO 31236
*  |  +----------------------------------------------------------------*
*  |  |  Here adjusting kinematics!!!!!
*  |  |  Now, chain 2 is a single particle jet, so we have to reset the
*  |  |  kinematical parameters xp, xxp, xt and xxt
*  |  |
31236       CONTINUE
            XXSQ2  = XXSQ  * XXSQ
*  |  |  +-------------------------------------------------------------*
*  |  |  |  If also chain 1 is a parjet (nnch1 .ne. 0) then the paramet.
*  |  |  |  to be recomputed are xp and xt, from alpha and beta, in
*  |  |  |  such a way to conserve the original momentum direction
*  |  |  |
            IF (NNCH1 .NE. 0) THEN
               XSQ2 = XSQ * XSQ
               HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 * (XXSQ2 - 2.D0)
     &                 - 2.D0 * XSQ2 * XXSQ2
               DDIFF = SQRT (HELP + 1.D0)
               SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XXSQ2)
               ALPHA = (SSUM + DDIFF) * 0.5D0
               BETA  = (SSUM - DDIFF) * 0.5D0
               IF (XP .GT. XT) THEN
                  XP  = ALPHA
                  XT  = BETA
               ELSE
                  XT  = ALPHA
                  XP  = BETA
               END IF
               XXP = 1.D0 - XP
               XXT = 1.D0 - XT
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            ELSE
*  |  |  +-------------------------------------------------------------*
*  |  |  |  If chain 1 is not a parjet (nnch1 .eq. 0) then the paramet.
*  |  |  |  xp, xt are to be recomputed in such a way to conserve the
*  |  |  |  original momentum direction and modulus
*  |  |  |
               DDIFF = XP  - XT
               SSUM  = SQRT (4.D0 * XXSQ2 + DDIFF**2)
               XP  = (SSUM + DDIFF) * 0.5D0
               XT  = (SSUM - DDIFF) * 0.5D0
               XXP = 1.D0 - XP
               XXT = 1.D0 - XT
               XSQ   = SQRT (XXP * XXT)
               AMCH1 = XSQ * UMO
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |                             end kinematics correction
*  |  +----------------------------------------------------------------*
 
3123     CONTINUE
C
         PCH1 = UMO*(XXP - XXT)*.5D0
         ECH1 = UMO*(XXP + XXT)*.5D0
         GAMCH1 = ECH1/AMCH1
         BGCH1  = PCH1/AMCH1
         PCH2 = UMO*(XP - XT)*.5D0
         ECH2 = UMO*(XP + XT)*.5D0
         GAMCH2 = ECH2/AMCH2
         BGCH2  = PCH2/AMCH2
         GO TO 34
*  |              end kin. sel. for antibaryon projectile
*  +-->-->-->-->-->-->-->-->-->--> go to 34
 
*  +-------------------------------------------------------------------*
*  |   baryon projectile!!!
*  |
  33  CONTINUE
C
C********************************************************************
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C*** NUCLEON - NUCLEON
C*** FORWARD DIQUARK - QUARK CHAIN
C
C********************************************************************
C
         GO TO 332
*  |
*  +-->-->-->-->-->-->-->-->-->--> jump # 1 to 332
 
*  +--<--<--<--<--<--<--<--<--<--< here from jump # 3
*  |
 
3310     CONTINUE
         CALL BKLASS (IFB,IFF1,IFF2,IF8,IF10)
         AMF8  = AM(IF8)
         AMF10 = AM(IF10)
         AMFF  = AMF10 + 0.3D0
         NNCH1 = 0
*  |  here xxp, xt are used for the jet # 1
         XSQ = SQRT(ABS(XXP*XT))
         AMCH1 = UMO*XSQ
         AA8   = IF8
         AA10  = IF10
*or         IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &                  ,AA8,AA10,AMF8,AMF10
 
*  |  +----------------------------------------------------------------*
*  |  | This check added by A. Ferrari, to avoid negative x(x)p or
*  |  | x(x)t !!!!!?????? Maybe also "go to 33366" if we want to create
*  |  | the jet anyway
* I (A. Ferrari) think this check was missing and is
* actually needed, else we can get negative energies
         IF (AMCH1 .LT. AMF8) GO TO 25
*  |  |   if amch1 < amf8 xp and xt are resampled
*  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
         IF (AMCH1 .GT. AMF10) GO TO 331
C*** PRODUCE AMF8
            AMCH1 = AMF8
            NNCH1 = -1
            XSQ = AMF8/UMO
            GO TO 33366
 331     CONTINUE
         IF (AMCH1 .GT. AMFF) GO TO 333
C*** PRODUCE AMF10
            AMCH1 = AMF10
            NNCH1 = 1
            XSQ = AMF10/UMO
            GO TO 33366
*  |  +----------------------------------------------------------------*
*  |  |  Here adjusting kinematics!!!!!
*  |  |  Now, chain 1 is a single particle jet, so we have to reset the
*  |  |  kinematical parameters xp, xxp, xt and xxt
*  |  |
33366       CONTINUE
            XSQ2  = XSQ  * XSQ
*  |  |  +-------------------------------------------------------------*
*  |  |  |  If also chain 2 is a parjet (nnch2 .ne. 0) then the paramet.
*  |  |  |  to be recomputed are xxp and xt, from alpha and beta, in
*  |  |  |  such a way to conserve the original momentum direction
*  |  |  |
            IF (NNCH2 .NE. 0) THEN
               XXSQ2 = XXSQ * XXSQ
               HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 * (XXSQ2 - 2.D0)
     &                 - 2.D0 * XSQ2 * XXSQ2
               DDIFF = SQRT (HELP + 1.D0)
               SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XSQ2)
               ALPHA = (SSUM + DDIFF) * 0.5D0
               BETA  = (SSUM - DDIFF) * 0.5D0
               IF (XXP .GT. XT) THEN
                  XXP = ALPHA
                  XT  = BETA
               ELSE
                  XT  = ALPHA
                  XXP = BETA
               END IF
               XP  = 1.D0 - XXP
               XXT = 1.D0 - XT
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            ELSE
*  |  |  +-------------------------------------------------------------*
*  |  |  |  If chain 2 is not a parjet (nnch2 .eq. 0) then the paramet.
*  |  |  |  xxp,xt have to be recomputed in such a way to conserve the
*  |  |  |  original momentum direction and modulus
*  |  |  |
               DDIFF = XXP - XT
               SSUM  = SQRT (4.D0 * XSQ2 + DDIFF**2)
               XXP = (SSUM + DDIFF) * 0.5D0
               XT  = (SSUM - DDIFF) * 0.5D0
               XP  = 1.D0 - XXP
               XXT = 1.D0 - XT
               XXSQ  = SQRT (XP * XXT)
               AMCH2 = XXSQ * UMO
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |                             end kinematics correction
*  |  +----------------------------------------------------------------*
 
 333     CONTINUE
*  |  we are using xp and xxt for chain 2
         PCH1 = UMO*(XXP - XT)*.5D0
         ECH1 = UMO*(XXP + XT)*.5D0
         GAMCH1 = ECH1/AMCH1
         BGCH1  = PCH1/AMCH1
         PCH2 = UMO*(XP - XXT)*.5D0
         ECH2 = UMO*(XP + XXT)*.5D0
         GAMCH2 = ECH2/AMCH2
         BGCH2  = PCH2/AMCH2
         GO TO 34
*  |              end kin. sel. for baryon projectile
*  +-->-->-->-->-->-->-->-->-->--> go to 34
 
C
C********************************************************************
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C     NUCLEON NUCLEON
C*** BACKWARD QUARK - DIQUARK CHAIN
C
C********************************************************************
C
*  +--<--<--<--<--<--<--<--<--<--< here from the previous jump # 1
*  |
 332     CONTINUE
*  |  Starting from chain # 2!!!
         CALL BKLASS (IBF,IBB1,IBB2,IB8,IBIO)
         NNCH2 = 0
         AMB8  = AM(IB8)
         AMB10 = AM(IBIO)
         AMBB  = AMB10 + 0.3D0
*  |  here xp, xxt are used for the jet # 2
         XXSQ  = SQRT(XP*XXT)
         AMCH2 = UMO*XXSQ
         AAB8  = IB8
         AAB10 = IBIO
*or         IF (IPRI.EQ.1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &                  ,AAB8,AAB10,AMB8,AMB10
 
*  |  +----------------------------------------------------------------*
*  |  | This check added by A. Ferrari, to avoid negative x(x)p or
*  |  | x(x)t !!!!!??????  Maybe also "go to 335" if we want to create
*  |  | the jet anyway
* I (A. Ferrari) think this check was missing and is
* actually needed, else we can get negative energies
         IF (AMCH2 .LT. AMB8) GO TO 25
*  |  |   if amch2 < amb8 xp and xt are resampled
*  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
         IF (AMCH2 .GT. AMB10) GO TO 334
C*** PRODUCE AMB8
            AMCH2 = AMB8
            NNCH2 = -1
            XXSQ = AMB8/UMO
            XP   = XXSQ**2/XXT
            XXP  = 1.D0 - XP
            GO TO 335
 334     CONTINUE
         IF (AMCH2 .GT. AMBB) GO TO 335
C*** PRODUCE AMB10
            AMCH2 = AMB10
            NNCH2 = 1
            XXSQ = AMB10/UMO
            XP   = XXSQ**2/XXT
            XXP  = 1.D0 - XP
 
C     PCH1=UMO*(XXP-XT)*.5D0
C     ECH1=UMO*(XXP+XT)*.5D0
C     GAMCH1=ECH1/AMCH1
C     BGCH1=PCH1/AMCH1
C     GO TO 335
 
 335     CONTINUE
         GO TO 3310
*  |
*  +-->-->-->-->-->-->-->-->-->--> jump # 3 to 3310
 
*  +-------------------------------------------------------------------*
*  |   meson projectile!!!
*  |
  32  CONTINUE
C
C********************************************************************
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C*** MESON NUCLEON
C
C********************************************************************
C
         IF (IBTARG)3277,3288,3299
*  |  +----------------------------------------------------------------*
*  |  |  meson projectile, baryon target!!!!
*  |  |
3299     CONTINUE
            GO TO (321,325),ISAM3
*  |  |  +-------------------------------------------------------------*
*  |  |  |  meson projectile, baryon target, isam3 = 1
*  |  |  |
 321        CONTINUE
C*** MESON NUCLEON Q(XXP)-QQ(XXT)+AQ(XP-Q(XT)
               GO TO 322
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 4 to 322
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 5
*  |  |  |
3259           CONTINUE
C=================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C      MESON NUCLEON
C*** FIRST LONG Q(XXP)-QQ(XXT) CHAIN
C
C===================================================================
               CALL BKLASS (IFF,IFB1,IFB2,IF8,IFIO)
               AMF8  = AM(IF8)
               AMF10 = AM(IFIO)
               NNCH1 = 0
               AMFF  = AMF10 + 0.3D0
*  |  here xxp, xxt are used for the jet # 1
               XSQ   = SQRT(XXP*XXT)
               AMCH1 = UMO*XSQ
               AA8   = IF8
               AA10  = IFIO
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AA8,AA10,AMF8,AMF10
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH1 .LT. AMF8) GO TO 25
*  |  |  |  | if amch1 < amf8 xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH1 .GT. AMF10) GO TO 3211
C*** PRODUCE AMF8
                  AMCH1 = AMF8
                  NNCH1 = -1
                  XSQ = AMF8/UMO
                  GO TO 32136
3211           CONTINUE
               IF (AMCH1 .GT. AMFF) GO TO 3213
C*** PRODUCE AMF10
                  AMCH1 = AMF10
                  NNCH1 = 1
                  XSQ = AMF10/UMO
                  GO TO 32136
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Here adjusting kinematics!!!!!
*  |  |  |  |  Now, chain 1 is a single particle jet, so we have to
*  |  |  |  |  reset the kinematical parameters xp, xxp, xt and xxt
*  |  |  |  |
32136             CONTINUE
                  XSQ2  = XSQ  * XSQ
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If also chain 2 is a parjet (nnch2 .ne. 0) then the
*  |  |  |  |  |  param. to be recomputed are xxp and xxt, from alpha
*  |  |  |  |  |  and beta, in such a way to conserve the original
*  |  |  |  |  |  momentum direction
*  |  |  |  |  |
                  IF (NNCH2 .NE. 0) THEN
                     XXSQ2 = XXSQ * XXSQ
                     HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 *
     &                      (XXSQ2 - 2.D0) - 2.D0 * XSQ2 * XXSQ2
                     DDIFF = SQRT (HELP + 1.D0)
                     SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XSQ2)
                     ALPHA = (SSUM + DDIFF) * 0.5D0
                     BETA  = (SSUM - DDIFF) * 0.5D0
                     IF (XXP .GT. XXT) THEN
                         XXP = ALPHA
                         XXT = BETA
                     ELSE
                         XXT = ALPHA
                         XXP = BETA
                     END IF
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  ELSE
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If chain 2 is not a parjet (nnch2 .eq. 0) then the
*  |  |  |  |  |  paramet. xxp,xxt have to be recomputed in such a way
*  |  |  |  |  |  to conserve the original momentum direction and
*  |  |  |  |  |  modulus
*  |  |  |  |  |
                     DDIFF = XXP - XXT
                     SSUM  = SQRT (4.D0 * XSQ2 + DDIFF**2)
                     XXP = (SSUM + DDIFF) * 0.5D0
                     XXT  = (SSUM - DDIFF) * 0.5D0
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
                     XXSQ  = SQRT (XP * XT)
                     AMCH2 = XXSQ * UMO
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |                       end kinematics correction
*  |  |  |  +----------------------------------------------------------*
 
3213           CONTINUE
               PCH1 = UMO*(XXP - XXT)*.5D0
               ECH1 = UMO*(XXP + XXT)*.5D0
               GAMCH1 = ECH1/AMCH1
               BGCH1  = PCH1/AMCH1
               PCH2 = UMO*(XP - XT)*.5D0
               ECH2 = UMO*(XP + XT)*.5D0
               GAMCH2 = ECH2/AMCH2
               BGCH2  = PCH2/AMCH2
               GO TO 34
*  |  |  | end kin. sel. for meson proj. (baryon target), isam3 = 1
*  |  |  +-->-->-->-->-->-->-->-->-->--> go to 34
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 4
*  |  |  |
 322           CONTINUE
C===============================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C       MESON NUCLEON
C*** SHORT AQ(XP)-Q(XT) CHAIN
C
C================================================================
               IIBF = IABS(IBF)
               IBPS = IMPS(IIBF,IBB)
               IBV  = IMVE(IIBF,IBB)
               AMBPS = AM(IBPS)
               AMBV  = AM(IBV)
               NNCH2 = 0
               AMBB  = AMBV + 0.3D0
*  |  here xp, xt are used for the jet # 2
               XXSQ  = SQRT(XP*XT)
               AMCH2 = UMO*XXSQ
               AAPS  = IBPS
               AAV   = IBV
*or               IF (IPRI.EQ.1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AAPS,AAV,AMBPS,AMBV
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH2 .LT. AMBPS) GO TO 25
*  |  |  |  | if amch2 < ambps xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH2 .GT. AMBV) GO TO 3221
C*** PRODUCE AMBPS
                  AMCH2 = AMBPS
                  NNCH2 = -1
                  XXSQ  = AMBPS/UMO
                  XT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 3223
                  XXT = 1.D0 - XT
                  GO TO 3223
3221           CONTINUE
               IF (AMCH2 .GT. AMBB) GO TO 3223
C*** PRODUCE AMBV
                  AMCH2 = AMBV
                  NNCH2 = 1
                  XXSQ  = AMBV/UMO
                  XT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 3223
                  XXT = 1.D0 - XT
 
C     PCH1=UMO*(XXP-XXT)*.5D0
C     ECH1=UMO*(XXP+XXT)*.5D0
C     GAMCH1=ECH1/AMCH1
C     BGCH1=PCH1/AMCH1
C     GO TO 3223
 
3223           CONTINUE
               GO TO 3259
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 5 to 3259
 
*  |  |  +-------------------------------------------------------------*
*  |  |  |  meson projectile, baryon target, isam3 = 2
*  |  |  |
 325        CONTINUE
C=================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C*** MESON NUCLEUS
C*** FORWARD ANTIQUARK-DIQUARK CHAIN
C
C=================================================================
               GO TO 326
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 6 to 326
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 7
*  |  |  |
3250           CONTINUE
C*** MESON NUCLEON FORWARD AQ(XXP)-Q(XT) AND BACKWARD CHAINS Q(XP)-QQ(
C*** XXT) CHAINS
C=====================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C     MESON NUCLEON
C*** FORWARD AQ(XXP)-Q(XT) CHAIN
C
C====================================================================
               IIFF = IABS(IFF)
               IFPS = IMPS(IIFF,IFB)
               IFV  = IMVE(IIFF,IFB)
               AMFPS = AM(IFPS)
               AMFV  = AM(IFV)
               NNCH1 = 0
               AMFF = AMFV + 0.3D0
*  |  |  |  here xxp, xt are used for the jet # 1
               XSQ  = SQRT(XXP*XT)
               AMCH1 = UMO*XSQ
               AAPS  = IFPS
               AAV   = IFV
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AAPS,AAV,AMFPS,AMFV
 
*  |  |  +-------------------------------------------------------------*
*  |  |  |
               IF (AMCH1 .LT. AMFPS) GO TO 25
*  |  |  | if amch1 < amfps xp and xt are resampled
*  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH1 .GT. AMFV) GO TO 3251
C*** PRODUCE AMFPS
                  AMCH1 = AMFPS
                  NNCH1 = -1
                  XSQ = AMFPS/UMO
                  GO TO 32536
3251           CONTINUE
               IF (AMCH1.GT.AMFF) GO TO 3253
C*** PRODUCE AMFV
                  AMCH1 = AMFV
                  NNCH1 = 1
                  XSQ = AMFV/UMO
                  GO TO 32536
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Here adjusting kinematics!!!!!
*  |  |  |  |  Now, chain 1 is a single particle jet, so we have to
*  |  |  |  |  reset the kinematical parameters xp, xxp, xt and xxt
*  |  |  |  |
32536             CONTINUE
                  XSQ2  = XSQ  * XSQ
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If also chain 2 is a parjet (nnch2 .ne. 0) then the
*  |  |  |  |  |  param. to be recomputed are xxp and xt, from alpha
*  |  |  |  |  |  and beta, in such a way to conserve the original
*  |  |  |  |  |  momentum direction
*  |  |  |  |  |
                  IF (NNCH2 .NE. 0) THEN
                     XXSQ2 = XXSQ * XXSQ
                     HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 *
     &                      (XXSQ2 - 2.D0) - 2.D0 * XSQ2 * XXSQ2
                     DDIFF = SQRT (HELP + 1.D0)
                     SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XSQ2)
                     ALPHA = (SSUM + DDIFF) * 0.5D0
                     BETA  = (SSUM - DDIFF) * 0.5D0
                     IF (XXP .GT. XT) THEN
                         XXP = ALPHA
                         XT  = BETA
                     ELSE
                         XT  = ALPHA
                         XXP = BETA
                     END IF
                     XP  = 1.D0 - XXP
                     XXT = 1.D0 - XT
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  ELSE
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If chain 2 is not a parjet (nnch2 .eq. 0) then the
*  |  |  |  |  |  paramet. xxp,xt have to be recomputed in such a way
*  |  |  |  |  |  to conserve the original momentum direction and
*  |  |  |  |  |  modulus
*  |  |  |  |  |
                     DDIFF = XXP - XT
                     SSUM  = SQRT (4.D0 * XSQ2 + DDIFF**2)
                     XXP = (SSUM + DDIFF) * 0.5D0
                     XT  = (SSUM - DDIFF) * 0.5D0
                     XP  = 1.D0 - XXP
                     XXT = 1.D0 - XT
                     XXSQ  = SQRT (XP * XXT)
                     AMCH2 = XXSQ * UMO
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |                       end kinematics correction
*  |  |  |  +----------------------------------------------------------*
 
3253           CONTINUE
               PCH1 = UMO*(XXP - XT)*.5D0
               ECH1 = UMO*(XXP + XT)*.5D0
               GAMCH1 = ECH1/AMCH1
               BGCH1  = PCH1/AMCH1
               PCH2 = UMO*(XP - XXT)*.5D0
               ECH2 = UMO*(XP + XXT)*.5D0
               GAMCH2 = ECH2/AMCH2
               BGCH2  = PCH2/AMCH2
               GO TO 34
*  |  |  | end kin. sel. for meson proj. (baryon target), isam3 = 2
*  |  |  +-->-->-->-->-->-->-->-->-->--> go to 34
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 6
*  |  |  |
 326           CONTINUE
C*** BACKWARD Q(XP)-QQ(XXT) CHAIN
               CALL BKLASS (IBF,IBB1,IBB2,IB8,IBIO)
               AMB8  = AM(IB8)
               AMB10 = AM(IBIO)
               NNCH2 = 0
               AMBB  = AMB10 + 0.3D0
C===================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C      MESON NUCLEON
C*** BACKWARD QUARK -DIQUARK CHAIN
C
C====================================================================
*  |  |  |  here xp, xxt are used for the jet # 2
               XXSQ  = SQRT(XP*XXT)
               AMCH2 = UMO*XXSQ
               AA8   = IB8
               AA10  = IBIO
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XXSQ,AMCH1,AMCH2
*or     &         ,AA8,AA10,AMB8,AMB10
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH2 .LT. AMB8 ) GO TO 25
*  |  |  |  |  if amch2 < amb8 xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH2 .GT. AMB10) GO TO 3261
C*** PRODUCE AMB8
                  AMCH2 = AMB8
                  NNCH2 = -1
                  XXSQ = AMB8/UMO
                  XXT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 3263
                  XT  = 1.D0 - XXT
                  GO TO 3263
3261           CONTINUE
               IF (AMCH2 .GT. AMBB) GO TO 3263
C*** PRODUCE AMB10
                  AMCH2 = AMB10
                  NNCH2 = 1
                  XXSQ = AMB10/UMO
                  XXT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 3263
                  XT  = 1.D0 - XXT
 
C     PCH1=UMO*(XXP-XT)*.5D0
C     ECH1=UMO*(XXP+XT)*.5D0
C     GAMCH1=ECH1/AMCH1
C     BGCH1=PCH1/AMCH1
 
                  GO TO 3263
3263           CONTINUE
               GO TO 3250
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 7 to 3250
*  |  +----------------------------------------------------------------*
*  |
*  |  +----------------------------------------------------------------*
*  |  |  Meson projectile, meson target!!!
*  |  |
3288     CONTINUE
 
C================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C           MESON MESON
C
C================================================================
            CALL GRNDM(RNDM,1)
            IF (RNDM(1) .LE. 0.5D0) GO TO 93288
               XT  = XXT
               XXT = 1.D0 - XT
93288       CONTINUE
            GO TO (3218,3258),ISAM3
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Meson projectile, meson target, isam3 = 1
*  |  |  |
3218        CONTINUE
C==================================================================
C*** MESON MESON Q(XXP)-AQ(XXT)+AQ(XP)-Q(XT)
C=================================================================
               GO TO 3228
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 8 to 3228
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 9
*  |  |  |
32598          CONTINUE
C=================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C      MESON MESON
C*** FIRST LONG Q(XXP)-AQ(XXT) CHAIN
C
C===================================================================
               IIFB = IABS(IFB)
               IFPS = IMPS(IIFB,IFF)
               IFV  = IMVE(IIFB,IFF)
               AMFPS = AM(IFPS)
*  |  |  |  Of course AMPV seems to be a mistyping for AMFV
*              AMPV  = AM(IFV)
               AMFV  = AM(IFV)
               NNCH1 = 0
               AMFF = AMFV + 0.3D0
*  |  |  |  here we are using xxp, xxt for the jet # 1
               XSQ  = SQRT(XXP*XXT)
               AMCH1 = UMO*XSQ
               AAPS  = IFPS
               AAV   = IFV
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AAPS,AAV,AMFPS,AMFV
 
*  |  |  +-------------------------------------------------------------*
*  |  |  |
               IF (AMCH1 .LT. AMFPS) GO TO 25
*  |  |  | if amch2 < amfps xp and xt are resampled
*  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH1 .GT. AMFV ) GO TO 32118
C*** PRODUCE AMFPS
                  AMCH1 = AMFPS
                  NNCH1 = -1
                  XSQ = AMFPS/UMO
                  GO TO 32133
32118          CONTINUE
               IF (AMCH1 .GT. AMFF) GO TO 32138
C*** PRODUCE AMFV
                  AMCH1 = AMFV
                  NNCH1 = 1
                  XSQ = AMFV/UMO
                  GO TO 32133
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Here adjusting kinematics!!!!!
*  |  |  |  |  Now, chain 1 is a single particle jet, so we have to
*  |  |  |  |  reset the kinematical parameters xp, xxp, xt and xxt
*  |  |  |  |
32133             CONTINUE
                  XSQ2  = XSQ  * XSQ
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If also chain 2 is a parjet (nnch2 .ne. 0) then the
*  |  |  |  |  |  param. to be recomputed are xxp and xxt, from alpha
*  |  |  |  |  |  and beta, in such a way to conserve the original
*  |  |  |  |  |  momentum direction
*  |  |  |  |  |
                  IF (NNCH2 .NE. 0) THEN
                     XXSQ2 = XXSQ * XXSQ
                     HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 *
     &                      (XXSQ2 - 2.D0) - 2.D0 * XSQ2 * XXSQ2
                     DDIFF = SQRT (HELP + 1.D0)
                     SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XSQ2)
                     ALPHA = (SSUM + DDIFF) * 0.5D0
                     BETA  = (SSUM - DDIFF) * 0.5D0
                     IF (XXP .GT. XXT) THEN
                         XXP = ALPHA
                         XXT = BETA
                     ELSE
                         XXT = ALPHA
                         XXP = BETA
                     END IF
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  ELSE
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If chain 2 is not a parjet (nnch2 .eq. 0) then the
*  |  |  |  |  |  paramet. xxp,xxt have to be recomputed in such a way
*  |  |  |  |  |  to conserve the original momentum direction and
*  |  |  |  |  |  modulus
*  |  |  |  |  |
                     DDIFF = XXP - XXT
                     SSUM  = SQRT (4.D0 * XSQ2 + DDIFF**2)
                     XXP = (SSUM + DDIFF) * 0.5D0
                     XXT = (SSUM - DDIFF) * 0.5D0
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
                     XXSQ  = SQRT (XP * XT)
                     AMCH2 = XXSQ * UMO
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |                       end kinematics correction
*  |  |  |  +----------------------------------------------------------*
 
32138          CONTINUE
               PCH1 = UMO*(XXP - XXT)*.5D0
               ECH1 = UMO*(XXP + XXT)*.5D0
               GAMCH1 = ECH1/AMCH1
               BGCH1  = PCH1/AMCH1
               PCH2 = UMO*(XP - XT)*.5D0
               ECH2 = UMO*(XP + XT)*.5D0
               GAMCH2 = ECH2/AMCH2
               BGCH2  = PCH2/AMCH2
               GO TO 348
*  |  |  | end kin. sel. for meson proj. (meson target), isam3 = 1
*  |  |  +-->-->-->-->-->-->-->-->-->--> go to 348
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 8
*  |  |  |
3228           CONTINUE
C===============================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C       MESON MESON
C*** SHORT AQ(XP)-Q(XT) CHAIN
C
C================================================================
               IIBF = IABS(IBF)
               IBPS = IMPS(IIBF,IBB)
               IBV  = IMVE(IIBF,IBB)
               AMBPS = AM(IBPS)
               AMBV  = AM(IBV)
               NNCH2 = 0
               AMBB  = AMBV + 0.3D0
*  |  |  |  here we are using xp, xt for the jet # 2
               XXSQ  = SQRT(XP*XT)
               AMCH2 = UMO*XXSQ
               AAPS = IBPS
               AAV  = IBV
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XXSQ,AMCH1,AMCH2
*or     &         ,AAPS,AAV,AMBPS,AMBV
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH2 .LT. AMBPS) GO TO 25
*  |  |  |  | if amch2 < ambps xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH2 .GT. AMBV ) GO TO 32218
C*** PRODUCE AMBPS
                  AMCH2 = AMBPS
                  NNCH2 = -1
                  XXSQ  = AMBPS/UMO
                  XT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 32238
                  XXT = 1.D0 - XT
                  GO TO 32238
32218          CONTINUE
               IF (AMCH2 .GT. AMBB) GO TO 32238
C*** PRODUCE AMBV
                  AMCH2 = AMBV
                  NNCH2 = 1
                  XXSQ  = AMBV/UMO
                  XT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 32238
                  XXT = 1.D0 - XT
 
C     PCH1=UMO*(XXP-XXT)*.5D0
C     ECH1=UMO*(XXP+XXT)*.5D0
C     GAMCH1=ECH1/AMCH1
C     BGCH1=PCH1/AMCH1
C     GO TO 32238
 
32238          CONTINUE
               GO TO 32598
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 9 to 32598
 
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Meson projectile, meson target, isam3 = 2
*  |  |  |
3258        CONTINUE
C=================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C*** MESON MESON
C*** FORWARD ANTIQUARK-DIQUARK CHAIN
C
C=================================================================
               GO TO 3268
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 10 to 3268
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 11
*  |  |  |
32508          CONTINUE
C===================================================================
C*** MESON MESON FORWARD AQ(XXP)-Q(XT) AND BACKWARD CHAINS Q(XP)-AQ(
C*** XXT) CHAINS
C====================================================================
C=====================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C     MESON MESON
C*** FORWARD AQ(XXP)-Q(XT) CHAIN
C
C====================================================================
               IIFF = IABS(IFF)
               IFPS = IMPS(IIFF,IFB)
               IFV  = IMVE(IIFF,IFB)
               AMFPS = AM(IFPS)
               AMFV  = AM(IFV)
               NNCH1 = 0
               AMFF  = AMFV + 0.3D0
*  |  |  |  here we are using xxp, xt for the jet # 1
               XSQ   = SQRT(XXP*XT)
               AMCH1 = UMO*XSQ
               AAPS = IFPS
               AAV  = IFV
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AAPS,AAV,AMFPS,AMFV
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH1 .LT. AMFPS) GO TO 25
*  |  |  |  | if amch1 < amfps xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH1 .GT. AMFV) GO TO 32518
C*** PRODUCE AMFPS
                  AMCH1 = AMFPS
                  NNCH1 = -1
                  XSQ = AMFPS/UMO
                  GO TO 32535
32518          CONTINUE
               IF (AMCH1 .GT. AMFF) GO TO 32538
C*** PRODUCE AMFV
                  AMCH1 = AMFV
                  NNCH1 = 1
                  XSQ = AMFV/UMO
                  GO TO 32535
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Here adjusting kinematics!!!!!
*  |  |  |  |  Now, chain 1 is a single particle jet, so we have to
*  |  |  |  |  reset the kinematical parameters xp, xxp, xt and xxt
*  |  |  |  |
32535             CONTINUE
                  XSQ2  = XSQ  * XSQ
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If also chain 2 is a parjet (nnch2 .ne. 0) then the
*  |  |  |  |  |  param. to be recomputed are xxp and xt, from alpha
*  |  |  |  |  |  and beta, in such a way to conserve the original
*  |  |  |  |  |  momentum direction
*  |  |  |  |  |
                  IF (NNCH2 .NE. 0) THEN
                     XXSQ2 = XXSQ * XXSQ
                     HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 *
     &                      (XXSQ2 - 2.D0) - 2.D0 * XSQ2 * XXSQ2
                     DDIFF = SQRT (HELP + 1.D0)
                     SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XSQ2)
                     ALPHA = (SSUM + DDIFF) * 0.5D0
                     BETA  = (SSUM - DDIFF) * 0.5D0
                     IF (XXP .GT.  XT) THEN
                         XXP = ALPHA
                         XT  = BETA
                     ELSE
                         XT  = ALPHA
                         XXP = BETA
                     END IF
                     XP  = 1.D0 - XXP
                     XXT = 1.D0 - XT
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  ELSE
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If chain 2 is not a parjet (nnch2 .eq. 0) then the
*  |  |  |  |  |  paramet. xxp,xt have to be recomputed in such a way
*  |  |  |  |  |  to conserve the original momentum direction and
*  |  |  |  |  |  modulus
*  |  |  |  |  |
                     DDIFF = XXP - XT
                     SSUM  = SQRT (4.D0 * XSQ2 + DDIFF**2)
                     XXP = (SSUM + DDIFF) * 0.5D0
                     XT  = (SSUM - DDIFF) * 0.5D0
                     XP  = 1.D0 - XXP
                     XXT = 1.D0 - XT
                     XXSQ  = SQRT (XP * XXT)
                     AMCH2 = XXSQ * UMO
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |                       end kinematics correction
*  |  |  |  +----------------------------------------------------------*
 
32538          CONTINUE
               PCH1 = UMO*(XXP - XT)*.5D0
               ECH1 = UMO*(XXP + XT)*.5D0
               GAMCH1 = ECH1/AMCH1
               BGCH1  = PCH1/AMCH1
               PCH2 = UMO*(XP - XXT)*.5D0
               ECH2 = UMO*(XP + XXT)*.5D0
               GAMCH2 = ECH2/AMCH2
               BGCH2  = PCH2/AMCH2
               GO TO 348
*  |  |  | end kin. sel. meson proj. (meson target), isam3 = 2
*  |  |  +-->-->-->-->-->-->-->-->-->--> go to 348
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 10
*  |  |  |
3268           CONTINUE
C*** BACKWARD Q(XP)-AQ(XXT) CHAIN
               IIBB = IABS(IBB)
               IBPS = IMPS(IIBB,IBF)
               IBV  = IMVE(IIBB,IBF)
               AMBPS = AM(IBPS)
               AMBV  = AM(IBV)
               NNCH2 = 0
               AMBB  = AMBV + 0.3D0
C===================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C      MESON MESON
C*** BACKWARD QUARK -ANTIQUARK CHAIN
C
C====================================================================
*  |  |  |  here we are using xp, xxt for the jet # 2
               XXSQ  = SQRT(XP*XXT)
               AMCH2 = UMO*XXSQ
               AAPS  = IBPS
               AAV   = IBV
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AAPS,AAV,AMBPS,AMBV
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH2 .LT. AMBPS) GO TO 25
*  |  |  |  | if amch2 < ambps xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH2 .GT. AMBV) GO TO 32618
C*** PRODUCE AMBPS
                  AMCH2 = AMBPS
                  NNCH2 = -1
                  XXSQ = AMBPS/UMO
                  XXT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 32638
                  XT = 1.D0 - XXT
                  GO TO 32638
32618          CONTINUE
               IF (AMCH2 .GT. AMBB) GO TO 32638
C*** PRODUCE AMBV
                  AMCH2 = AMBV
                  NNCH2 = 1
                  XXSQ = AMBV /UMO
                  XXT  = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 32638
                  XT  = 1.D0 - XXT
 
C     PCH1=UMO*(XXP-XT)*.5D0
C     ECH1=UMO*(XXP+XT)*.5D0
C     GAMCH1=ECH1/AMCH1
C     BGCH1=PCH1/AMCH1
 
                  GO TO 32638
32638          CONTINUE
               GO TO 32508
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 11 to 32508
*  |  +----------------------------------------------------------------*
 
 348  CONTINUE
*or      IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XP,XT,XXP,XXT
*or      IF (IPRI .EQ. 1) WRITE(LUNOUT,103)AMCH1,AMCH2,PCH1,PCH2,ECH1,
*or     &ECH2,GAMCH1,GAMCH2,BGCH1,BGCH2
      GO TO 34
*  |  | end kin. sel. meson proj. (meson target)
*  |  +-->-->-->-->-->-->-->-->-->--> go to 34
 
*  |  +----------------------------------------------------------------*
*  |  | meson projectile, antibaryon target!!!
*  |  |
3277     CONTINUE
C=================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C       MESON ANTIBARYON
C
C===============================================================
C====================================================================
C====================================================================
C     TO BE CORRECTED
C
C====================================================================
C`===================================================================
            GO TO (3217,3257),ISAM3
*  |  |  +-------------------------------------------------------------*
*  |  |  |  meson projectile, antibaryon target, isam = 2
*  |  |  |
3257        CONTINUE
C*** MESON NUCLEON AQ(XXP)-AQAQ(XXT)+Q(XP)-AQ(XT)
               GO TO 3227
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 12 to 3227
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 13
*  |  |  |
32597 CONTINUE
C=================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C      MESON ANTINUCLEON
C*** FIRST LONG AQ(XXP)-AQAQ(XXT) CHAIN
C
C===================================================================
               CALL BKLASS (IBF,IBB1,IBB2,IF8,IFIO)
               AMF8  = AM(IF8)
               AMF10 = AM(IFIO)
               NNCH1 = 0
               AMFF  = AMF10 + 0.3D0
*  |  |  |  here we are using xxp, xxt for the jet # 1
               XSQ   = SQRT(XXP*XXT)
               AMCH1 = UMO*XSQ
               AA8   = IF8
               AA10  = IFIO
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AA8,AA10,AMF8,AMF10
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH1 .LT. AMF8) GO TO 25
*  |  |  |  | if amch1 < amf8 xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH1 .GT. AMF10) GO TO 32117
C*** PRODUCE AMF8
                  AMCH1 = AMF8
                  NNCH1 = -1
                  XSQ = AMF8/UMO
                  GO TO 32135
32117          CONTINUE
               IF (AMCH1 .GT. AMFF) GO TO 32137
C*** PRODUCE AMF10
                  AMCH1 = AMF10
                  NNCH1 = 1
                  XSQ = AMF10/UMO
                  GO TO 32135
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Here adjusting kinematics!!!!!
*  |  |  |  |  Now, chain 1 is a single particle jet, so we have to
*  |  |  |  |  reset the kinematical parameters xp, xxp, xt and xxt
*  |  |  |  |
32135             CONTINUE
                  XSQ2  = XSQ  * XSQ
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If also chain 2 is a parjet (nnch2 .ne. 0) then the
*  |  |  |  |  |  param. to be recomputed are xxp and xxt, from alpha
*  |  |  |  |  |  and beta, in such a way to conserve the original
*  |  |  |  |  |  momentum direction
*  |  |  |  |  |
                  IF (NNCH2 .NE. 0) THEN
                     XXSQ2 = XXSQ * XXSQ
                     HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 *
     &                      (XXSQ2 - 2.D0) - 2.D0 * XSQ2 * XXSQ2
                     DDIFF = SQRT (HELP + 1.D0)
                     SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XSQ2)
                     ALPHA = (SSUM + DDIFF) * 0.5D0
                     BETA  = (SSUM - DDIFF) * 0.5D0
                     IF (XXP .GT. XXT) THEN
                         XXP = ALPHA
                         XXT = BETA
                     ELSE
                         XXT = ALPHA
                         XXP = BETA
                     END IF
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  ELSE
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If chain 2 is not a parjet (nnch2 .eq. 0) then the
*  |  |  |  |  |  paramet. xxp,xxt have to be recomputed in such a way
*  |  |  |  |  |  to conserve the original momentum direction and
*  |  |  |  |  |  modulus
*  |  |  |  |  |
                     DDIFF = XXP - XXT
                     SSUM  = SQRT (4.D0 * XSQ2 + DDIFF**2)
                     XXP = (SSUM + DDIFF) * 0.5D0
                     XXT = (SSUM - DDIFF) * 0.5D0
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
                     XXSQ  = SQRT (XP * XT)
                     AMCH2 = XXSQ * UMO
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |                       end kinematics correction
*  |  |  |  +----------------------------------------------------------*
 
32137          CONTINUE
               PCH1 = UMO*(XXP - XXT)*.5D0
               ECH1 = UMO*(XXP + XXT)*.5D0
               GAMCH1 = ECH1/AMCH1
               BGCH1  = PCH1/AMCH1
               PCH2 = UMO*(XP - XT)*.5D0
               ECH2 = UMO*(XP + XT)*.5D0
               GAMCH2 = ECH2/AMCH2
               BGCH2  = PCH2/AMCH2
               GO TO 34
*  |  |  | end kin. sel. meson proj. (abaryon target), isam3 = 2
*  |  |  +-->-->-->-->-->-->-->-->-->--> go to 34
 
*  |  |  +--<--<--<--<--<--<--<--<--<--< here from jump # 12
*  |  |  |
3227           CONTINUE
C===============================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C       MESON ANTINUCLEON
C*** SHORT AQ(XP)-Q(XT) CHAIN
C
C================================================================
              IIFB = IABS(IFB)
              IBPS = IMPS(IIFB,IFF)
              IBV  = IMVE(IIFB,IFF)
              AMBPS = AM(IBPS)
              AMBV  = AM(IBV)
              NNCH2 = 0
              AMBB = AMBV + 0.3D0
*  |  |  | here we are using xp,xt for jet # 2
              XXSQ  = SQRT(XP*XT)
              AMCH2 = UMO*XXSQ
              AAPS  = IBPS
              AAV   = IBV
*or              IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &        ,AAPS,AAV,AMBPS,AMBV
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH2 .LT. AMBPS) GO TO 25
*  |  |  |  | if amch2 < ambps xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH2 .GT. AMBV) GO TO 32217
C*** PRODUCE AMBPS
                  AMCH2 = AMBPS
                  NNCH2 = -1
                  XXSQ = AMBPS/UMO
                  XT   = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 32237
                  XXT = 1.D0 - XT
                  GO TO 32237
32217          CONTINUE
               IF (AMCH2 .GT. AMBB) GO TO 32237
C*** PRODUCE AMBV
                  AMCH2 = AMBV
                  NNCH2 = 1
                  XXSQ = AMBV/UMO
                  XT   = XXSQ**2/XP
*or                  IF (INUCVT .EQ. 1) GO TO 32237
                  XXT = 1.D0 - XT
 
C     PCH1=UMO*(XXP-XXT)*.5D0
C     ECH1=UMO*(XXP+XXT)*.5D0
C     GAMCH1=ECH1/AMCH1
C     BGCH1=PCH1/AMCH1
C     GO TO 32237
 
32237          CONTINUE
               GO TO 32597
*  |  |  |
*  |  |  +-->-->-->-->-->-->-->-->-->--> jump # 13 to 32597
 
*  |  |  +-------------------------------------------------------------*
*  |  |  |  meson projectile, antibaryon target, isam = 1
*  |  |  |
3217           CONTINUE
C=================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C*** MESON ANTINUCLEUS
C*** BACKWARD    QUARK-  ANTIQUARK CHAIN
C
C=================================================================
32507          CONTINUE
C=====================================================================
C
C     MESON ANTINUCLEON
C*** BACKWARD Q(XP)-AQ(XT) CHAIN
C
C====================================================================
               IIBB = IABS(IBB)
               IFPS = IMPS(IIBB,IBF)
               IFV  = IMVE(IIBB,IBF)
               AMFPS = AM(IFPS)
               AMFV  = AM(IFV)
               NNCH1 = 0
               AMFF = AMFV + 0.3D0
*  |  |  | here we are using xp,xt for jet # 1
               XSQ  = SQRT(XP*XT)
               AMCH1 = UMO*XSQ
               AAPS  = IFPS
               AAV   = IFV
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AAPS,AAV,AMFPS,AMFV
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH1 .LT. AMFPS) GO TO 25
*  |  |  |  | if amch1 < amfps xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH1 .GT. AMFV) GO TO 32517
C*** PRODUCE AMFPS
                  AMCH1 = AMFPS
                  NNCH1 = -1
                  XSQ = AMFPS/UMO
                  XP  = XSQ**2/XT
                  XXP = 1.D0 - XP
                  GO TO 32537
32517          CONTINUE
               IF (AMCH1 .GT. AMFF) GO TO 32537
C*** PRODUCE AMFV
                  AMCH1 = AMFV
                  NNCH1 = 1
                  XSQ = AMFV/UMO
                  XP  = XSQ**2/XT
                  XXP = 1.D0 - XP
                  GO TO 32537
32537          CONTINUE
               GO TO 3267
3267           CONTINUE
C***FORWARD AQ(XXP)-AQAQ(XXT) CHAIN
               CALL BKLASS(IFF,IFB1,IFB2,IB8,IBIO)
               AMB8  = AM(IB8)
               AMB10 = AM(IBIO)
               NNCH2 = 0
               AMBB  = AMB10 + 0.3D0
C===================================================================
C
C*** KINEMATICAL PARAMETERS OF BOTH CHAINS IN CMS
C               MESON ANTINUCLEON
C*** FORWARD ANTIQUARK -ANTIDIQUARK CHAIN
C
C====================================================================
*  |  |  | here we are using xxp,xXt for jet # 2
               XXSQ = SQRT(XXP*XXT)
               AMCH2 = UMO*XXSQ
               AA8   = IB8
               AA10  = IBIO
*or               IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XSQ,AMCH1,AMCH2
*or     &         ,AA8,AA10,AMB8,AMB10
 
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF (AMCH2 .LT. AMB8 ) GO TO 25
*  |  |  |  | if amch2 < amb8 xp and xt are resampled
*  |  |  |  +-->-->-->-->-->-->-->-->-->--> xp, xt resampling
 
               IF (AMCH2 .GT. AMB10) GO TO 32617
C*** PRODUCE AMB8
                  AMCH2 = AMB8
                  NNCH2 = -1
                  XXSQ = AMB8/UMO
*or                  IF (INUCVT .EQ. 1) GO TO 32637
                  GO TO 32636
32617          CONTINUE
               IF (AMCH2 .GT. AMBB) GO TO 32637
C*** PRODUCE AMB10
                  AMCH2 = AMB10
                  NNCH2 = 1
                  XXSQ = AMB10/UMO
*or                  IF (INUCVT .EQ. 1) GO TO 32637
 
C     PCH1=UMO*(XXP-XT)*.5D0
C     ECH1=UMO*(XXP+XT)*.5D0
C     GAMCH1=ECH1/AMCH1
C     BGCH1=PCH1/AMCH1
* Here there was a "large" mistake in the old Hadevt!!!
                  GO TO 32636
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Here adjusting kinematics!!!!!
*  |  |  |  |  Now, chain 2 is a single particle jet, so we have to
*  |  |  |  |  reset the kinematical parameters xp, xxp, xt and xxt
*  |  |  |  |
32636             CONTINUE
                  XXSQ2 = XXSQ  * XXSQ
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If also chain 1 is a parjet (nnch1 .ne. 0) then the
*  |  |  |  |  |  param. to be recomputed are xxp and xxt, from alpha
*  |  |  |  |  |  and beta, in such a way to conserve the original
*  |  |  |  |  |  momentum direction
*  |  |  |  |  |
                  IF (NNCH1 .NE. 0) THEN
                     XSQ2  = XSQ * XSQ
                     HELP  = XSQ2 * (XSQ2 - 2.D0) + XXSQ2 *
     &                      (XXSQ2 - 2.D0) - 2.D0 * XSQ2 * XXSQ2
                     DDIFF = SQRT (HELP + 1.D0)
                     SSUM  = SQRT (HELP + 1.D0 + 4.D0 * XXSQ2)
                     ALPHA = (SSUM + DDIFF) * 0.5D0
                     BETA  = (SSUM - DDIFF) * 0.5D0
                     IF (XXP .GT.  XXT) THEN
                         XXP = ALPHA
                         XXT = BETA
                     ELSE
                         XXT = ALPHA
                         XXP = BETA
                     END IF
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  ELSE
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  If chain 1 is not a parjet (nnch1 .eq. 0) then the
*  |  |  |  |  |  paramet. xxp,xxt have to be recomputed in such a way
*  |  |  |  |  |  to conserve the original momentum direction and
*  |  |  |  |  |  modulus
*  |  |  |  |  |
                     DDIFF = XXP - XXT
                     SSUM  = SQRT (4.D0 * XXSQ2 + DDIFF**2)
                     XXP = (SSUM + DDIFF) * 0.5D0
                     XXT = (SSUM - DDIFF) * 0.5D0
                     XP  = 1.D0 - XXP
                     XT  = 1.D0 - XXT
                     XSQ   = SQRT (XP * XT)
                     AMCH1 = XSQ * UMO
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |                       end kinematics correction
*  |  |  |  +----------------------------------------------------------*
 
32637          CONTINUE
               PCH1 = UMO*(XP - XT)*.5D0
               ECH1 = UMO*(XP + XT)*.5D0
               GAMCH1 = ECH1/AMCH1
               BGCH1  = PCH1/AMCH1
               PCH2 = UMO*(XXP - XXT)*.5D0
               ECH2 = UMO*(XXP + XXT)*.5D0
               GAMCH2 = ECH2/AMCH2
               BGCH2  = PCH2/AMCH2
               GO TO 34
*  |  |  | end kin. sel. meson proj. (abaryon target), isam3 = 1
*  |  |  +-->-->-->-->-->-->-->-->-->--> go to 34
*  |  +----------------------------------------------------------------*
  34  CONTINUE
*  |            end of kinematical selections!!!!!!!!!
*  +-------------------------------------------------------------------*
*or      IF (IPRI .EQ. 1) WRITE(LUNOUT,103)XP,XT,XXP,XXT
*or      IF (IPRI .EQ. 1) WRITE(LUNOUT,103)AMCH1,AMCH2,PCH1,PCH2,
*or     &ECH1,ECH2,GAMCH1,GAMCH2,BGCH1,BGCH2
C
C********************************************************************
C
C*** MC SAMPLING OF FORWARD CHAIN
C
C********************************************************************
C
      IF (IBPROJ) 41,42,43
C==================================================================
C         FORWARD CHAIN OF ANTIBARYON BARYON
C==================================================================
  41  CONTINUE
         IF (NNCH1) 4111,4112,4113
4111     CONTINUE
            ICH1 = IFPS
            GO TO 4114
4113     CONTINUE
            ICH1 = IFV
            GO TO 4114
4112     CONTINUE
            IF (IOPBAM .EQ. 5) THEN
               IAIFF1 = IABS(IFF1) + 6
               IAIFF2 = IABS(IFF2) + 6
*or               IF (IPRI .EQ. 1)WRITE(LUNOUT,991)IFB1,IFB2,IAIFF1
*or 991           FORMAT (' BAMJEV 4112',5I5)
               CALL BAMJEV(IHAD,IFB1,IFB2,IAIFF1,IAIFF2,AMCH1,5)
            ELSE
               IAIFF = IABS(IAIFF) + 6
               CALL BAMJEV(IHAD,IFB,IAIFF,IFB,IFB,AMCH1,IOPBAM)
            END IF
            GO TO 4115
4114  CONTINUE
*or      IF (IPRI .EQ. 1) WRITE(LUNOUT,992)ICH1
*or  992 FORMAT (' PARJET 4114 ',5I5)
            CALL PARJET(IHAD,ICH1)
4115     CONTINUE
C
C
C
C     CALL DECAY(IHAD,2)
         GO TO 44
42    CONTINUE
C=====================================================================
C*** FORWARD CHAIN OF MESON NUCLEON
C======================================================================
         IF (IBTARG) 427,428,429
 429     CONTINUE
            GO TO (421,422),ISAM3
 421        CONTINUE
               IF (NNCH1) 4211,4212,4213
4211           CONTINUE
                  ICH1 = IF8
                  GO TO 4214
4213           CONTINUE
                  ICH1 = IFIO
                  GO TO 4214
4212           CONTINUE
                  CALL BAMJEV(IHAD,IFF,IFB1,IFB2,IFF,AMCH1,4)
                  GO TO 4215
4214           CONTINUE
                  CALL PARJET(IHAD,ICH1)
4215           CONTINUE
C     CALL DECAY(IHAD,2)
               GO TO 44
422         CONTINUE
C*** IFA - IFB AQ - Q EXISTIERT NICHT
               IF (NNCH1) 4221,4222,4223
4221           CONTINUE
                  ICH1 = IFPS
                  GO TO 4224
4223           CONTINUE
                  ICH1 = IFV
                  GO TO 4224
4222           CONTINUE
                  IAIFF = IABS(IFF) + 6
                  CALL BAMJEV(IHAD,IFB,IAIFF,IFB,IFB,AMCH1,3)
                  GO TO 4225
4224           CONTINUE
                  CALL PARJET(IHAD,ICH1)
4225           CONTINUE
C     CALL DECAY(IHAD,2)
               GO TO 44
C===================================================================
C
C   FORWARD CHAIN OF MESON MESON
C
C====================================================================
C=====================================================================
 428     CONTINUE
            GO TO (4218,4228),ISAM3
4218        CONTINUE
               IF (NNCH1) 42118,42128,42138
42118          CONTINUE
                  ICH1 = IFPS
                  GO TO 42148
42138          CONTINUE
                  ICH1 = IFV
                  GO TO 42148
42128          CONTINUE
                  IAIFB = IABS(IFB) + 6
                  CALL BAMJEV(IHAD,IFF,IAIFB,IFB2,IFF,AMCH1,3)
                  GO TO 42158
42148          CONTINUE
                  CALL PARJET(IHAD,ICH1)
42158          CONTINUE
C     CALL DECAY(IHAD,2)
               GO TO 44
4228        CONTINUE
C*** IFA - IFB AQ - Q EXISTIERT NICHT
               IF (NNCH1) 42218,42228,42238
42218          CONTINUE
                  ICH1 = IFPS
                  GO TO 42248
42238          CONTINUE
                  ICH1 = IFV
                  GO TO 42248
42228          CONTINUE
                  IAIFF = IABS(IFF) + 6
                  CALL BAMJEV(IHAD,IFB,IAIFF,IFB,IFB,AMCH1,3)
                  GO TO 42258
42248          CONTINUE
                  CALL PARJET(IHAD,ICH1)
42258          CONTINUE
C     CALL DECAY(IHAD,2)
               GO TO 44
C================================================================
C
C      FORWARD CHAIN OF MESON ANTIBARYON
C
C==================================================================
 427     CONTINUE
C=================================================================
C    TO BE CORRECTED
C=================================================================
            GO TO (4217,4227),ISAM3
4227        CONTINUE
               IF (NNCH1) 42117,42127,42137
42117          CONTINUE
                  ICH1 = IF8
                  GO TO 42147
42137          CONTINUE
                  ICH1 = IFIO
                  GO TO 42147
42127          CONTINUE
                  IAIBF = IABS(IBF) + 6
                  IAIBB1= IABS(IBB1) + 6
                  IAIBB2= IABS(IBB2) + 6
*or                  IF (IPRI.EQ.1) WRITE(LUNOUT,994)IAIBF,IAIBB1,IAIBB2
*or  994 FORMAT(' BAMJEV 43147',5I5)
                  CALL BAMJEV(IHAD,IAIBF,IAIBB1,IAIBB2,IBB,AMCH1,4)
                  GO TO 42157
42147          CONTINUE
*or                  IF (IPRI.EQ.1) WRITE(LUNOUT,993)ICH1
*or 993              FORMAT (' PARJET 42147',5I5)
                  CALL PARJET(IHAD,ICH1)
42157          CONTINUE
C     CALL DECAY(IHAD,2)
               GO TO 44
4217        CONTINUE
C*** IFA - IFB AQ - Q EXISTIERT NICHT
               IF (NNCH2) 42217,42227,42237
42217          CONTINUE
                  ICH1 = IB8
                  GO TO 42247
42237          CONTINUE
                  ICH1 = IBIO
                  GO TO 42247
42227          CONTINUE
                  IAIFF = IABS(IFF) + 6
                  IAIFB1= IABS(IFB1) + 6
                  IAIFB2= IABS(IFB2) + 6
*or                  IF (IPRI.EQ.1) WRITE(LUNOUT,995)IAIFF,IAIFB1,IAIFB2
*or 995              FORMAT (' BAMJEV 42227',5I5)
                  CALL BAMJEV(IHAD,IAIFF,IAIFB1,IAIFB2,IFB,AMCH2,4)
                  GO TO 42257
42247          CONTINUE
*or                  IF (IPRI.EQ.1) WRITE(LUNOUT,996)ICH1
*or 996              FORMAT ('PARJET 42247',5I5)
                  CALL PARJET(IHAD,ICH1)
42257          CONTINUE
C     CALL DECAY(IHAD,2)
               GO TO 44
  43  CONTINUE
C==================================================================
C*** FORWARD CHAIN OF NUCLEON NUCLEON
C===================================================================
         IF (NNCH1) 431,432,433
 431     CONTINUE
            ICH1 = IF8
            GO TO 434
 433     CONTINUE
            ICH1 = IF10
            GO TO 434
 432     CONTINUE
            CALL BAMJEV(IHAD,IFB,IFF1,IFF2,IFB,AMCH1,4)
            GO TO 435
 434     CONTINUE
            CALL PARJET(IHAD,ICH1)
 435     CONTINUE
C     CALL DECAY(IHAD,2)
  44  CONTINUE
C*** TURN CHAINS AROUND IF NECESSARY
      IF (IBPROJ) 51,52,53
  51  CONTINUE
         GO TO 55
  52  CONTINUE
C*** MESON NUCLEON
         GO TO (521,522),ISAM3
 521     CONTINUE
            GO TO 54
 522     CONTINUE
C*** TURN JET
            GO TO 55
  53  CONTINUE
C*** NUCLEON-NUCLEON
C*** TURN JET
         GO TO 55
  55  CONTINUE
C*** TURN JET AROUND
         DO 56 I=1,IHAD
            PZF(I) = -PZF(I)
  56     CONTINUE
  54  CONTINUE
      IIIHAD = IHAD
C*** AND  INT. CHAIN TRANSVERSE MOMENTA
      B3SAVE = B3BAMJ
* This is consistent with b3bamj = 10 for an initial value of 6
      B3BAMJ = 1.666666666666667D+00 * B3BAMJ
*     B3BAMJ = 5.D+00 / ( LOG10 ( 1.D+00 + ( UMO / E00 )**1.5D+00 )
*    &       + 1.D+00 )
      CALL GRNDM(RNDM,2)
      ES  = -2.D0/(B3BAMJ**2)*LOG(RNDM(1)*RNDM(2))
      B3BAMJ = B3SAVE
*     HPS = SQRT(ES*ES+2.D0*ES*AMCH1)
      HPS = SQRT(ES*ES+2.D0*ES*AM(1))
      CALL SFECFE(SFE,CFE)
*
*  tentative guess
*
      PTXCH1 = HPS * CFE
      PTYCH1 = HPS * SFE
*  +-------------------------------------------------------------------*
*  |                  Loop to establish the transverse momentum
      GO TO 6171
6170  CONTINUE
         PTXCH1 = 0.75D0 * PTXCH1
         PTYCH1 = 0.75D0 * PTYCH1
6171     CONTINUE
         IHAD = IIIHAD
*  |         The following two cards provide momentum conservation for
*  |         x and y components
         PTXCH2 = -PTXCH1
         PTYCH2 = -PTYCH1
         BGCH1X = PTXCH1/AMCH1
         BGCH1Y = PTYCH1/AMCH1
         ACH1 = BGCH1**2-(PTXCH1**2+PTYCH1**2)/AMCH1**2
         IF (ACH1.LE.0.D0) GO TO 6170
*  |
*  +--<--<--<--<--< if Pt is too large loop again on the forward jet
         BGCH2X = PTXCH2/AMCH2
         BGCH2Y = PTYCH2/AMCH2
         ACH2   = BGCH2**2-(PTXCH2**2+PTYCH2**2)/AMCH2**2
      IF (ACH2.LE.0.D0) GO TO 6170
*  |
*  +--<--<--<--<--<--<--<--<--<--< if Pt is too large loop again
      BGCH1Z = SQRT(ACH1)
      BGCH1Z = SIGN(BGCH1Z,BGCH1)
      BGCH2Z = SQRT(ACH2)
      BGCH2Z = SIGN(BGCH2Z,BGCH2)
 
      CALL LORTRA(IHAD,1,GAMCH1,BGCH1X,BGCH1Y,BGCH1Z)
C==============================================================
C
C*** TRANSFORM FORWARD JET INTO CMS
C
C================================================================
      IHAD = IIIHAD
      NAUX = IHAD
C===============================================================
C
C*** SAMPLING OF BACKWARD CHAIN
C
C===============================================================
      IF (IBPROJ) 61,62,63
  61  CONTINUE
C================================================================
C      BACKWARD CHAIN OF ANTINUCLEON NUCLEON
C=================================================================
         IF (NNCH2) 6111,6112,6113
6111     CONTINUE
            ICH2 = IBPS
            GO TO 6114
6113     CONTINUE
            ICH2 = IBV
            GO TO 6114
6112     CONTINUE
            IAIBF = IABS(IBF) + 6
            CALL BAMJEV(IHAD,IBB,IAIBF,IBB,IBB,AMCH2,3)
            GO TO 6115
6114     CONTINUE
            CALL PARJET(IHAD,ICH2)
6115     CONTINUE
C
C
C     CALL DECAY(IHAD,2)
         GO TO 64
  62  CONTINUE
C================================================================
C*** BACKWARD CHAIN OF MESON - BARYON
C==================================================================
         IF (IBTARG) 627,628,629
 629     CONTINUE
            GO TO (621,622),ISAM3
 621        CONTINUE
               IF (NNCH2) 6211,6212,6213
6211           CONTINUE
                 ICH2 = IBPS
                 GO TO 6214
6213          CONTINUE
                 ICH2 = IBV
                 GO TO 6214
6212          CONTINUE
                 IAIBF = IABS(IBF) + 6
                 CALL BAMJEV(IHAD,IBB,IAIBF,IBB,IBB,AMCH2,3)
                 GO TO 6215
6214          CONTINUE
                 CALL PARJET(IHAD,ICH2)
6215          CONTINUE
C     CALL DECAY(IHAD,2)
              GO TO 64
622        CONTINUE
              IF (NNCH2) 6221,6222,6223
6221          CONTINUE
                 ICH2 = IB8
                 GO TO 6224
6223          CONTINUE
                 ICH2 = IBIO
                 GO TO 6224
6222          CONTINUE
                 CALL BAMJEV(IHAD,IBF,IBB1,IBB2,IBF,AMCH2,4)
                 GO TO 6225
6224          CONTINUE
                 CALL PARJET(IHAD,ICH2)
6225             CONTINUE
C     CALL DECAY(IHAD,2)
                 GO TO 64
C==================================================================
C
C     BACKWARD CHAIN OF MESON MESON
C
C===================================================================
C        TO BE CORRECTED
C===================================================================
 628    CONTINUE
           GO TO(6218,6228),ISAM3
6218       CONTINUE
              IF (NNCH2) 62118,62128,62138
62118         CONTINUE
                 ICH2 = IBPS
                 GO TO 62148
62138         CONTINUE
                 ICH2 = IBV
                 GO TO 62148
62128         CONTINUE
                 IAIBF = IABS(IBF) + 6
                 CALL BAMJEV(IHAD,IBB,IAIBF,IBB,IBB,AMCH2,3)
                 GO TO 62158
62148         CONTINUE
                 CALL PARJET(IHAD,ICH2)
62158         CONTINUE
C     CALL DECAY(IHAD,2)
              GO TO 64
6228       CONTINUE
              IF (NNCH2) 62218,62228,62238
62218         CONTINUE
                 ICH2 = IBPS
                 GO TO 62248
62238         CONTINUE
                 ICH2 = IBV
                 GO TO 62248
62228         CONTINUE
                 IAIBB = IABS(IBB) + 6
                 CALL BAMJEV(IHAD,IBF,IAIBB,IBB2,IBF,AMCH2,3)
                 GO TO 62258
62248         CONTINUE
                 CALL PARJET(IHAD,ICH2)
62258 CONTINUE
C     CALL DECAY(IHAD,2)
              GO TO 64
C================================================================
C
C       BACKWARD CHAIN OF MESON ANTIBARYON
C=================================================================
C     TO BE CORRECTED
C=================================================================
 627     CONTINUE
            GO TO(6217,6227),ISAM3
6227        CONTINUE
               IF (NNCH2) 62117,62127,62137
62117          CONTINUE
                 ICH2 = IBPS
                 GO TO 62147
62137         CONTINUE
                 ICH2 = IBV
                 GO TO 62147
62127         CONTINUE
                 IAIFB = IABS(IFB) + 6
*or                 IF (IPRI.EQ.1) WRITE(LUNOUT,997)IFF,IAIFB
*or 997             FORMAT (' BAMJEV 62127',5I5)
                 CALL BAMJEV(IHAD,IFF,IAIFB,IBB,IBB,AMCH2,3)
                 GO TO 62157
62147         CONTINUE
*or                 IF (IPRI.EQ.1) WRITE(LUNOUT,998)ICH2
*or 998             FORMAT ('PARJET 62147',5I5)
                 CALL PARJET(IHAD,ICH2)
62157         CONTINUE
C     CALL DECAY(IHAD,2)
              GO TO 64
6217       CONTINUE
              IF (NNCH1) 62217,62227,62237
62217         CONTINUE
                 ICH2 = IFPS
                 GO TO 62247
62237         CONTINUE
                 ICH2 = IFV
                 GO TO 62247
62227         CONTINUE
                 IAIBB = IABS(IBB) + 6
*or                 IF (IPR1.EQ.1) WRITE(LUNOUT,9911)IBF,IAIBB
*or9911             FORMAT (' BAMJEV 62227',5I5)
                 CALL BAMJEV(IHAD,IBF,IAIBB,IAIFB2,IBF,AMCH1,3)
                 GO TO 62257
62247         CONTINUE
*or                 IF (IPRI .EQ. 1) WRITE(LUNOUT,9912)ICH2
*or9912             FORMAT ('PARJET 62247',5I5)
                 CALL PARJET(IHAD,ICH2)
62257         CONTINUE
C     CALL DECAY(IHAD,2)
              GO TO 64
63    CONTINUE
C==================================================================
C*** BACKWARD CHAIN OF BARYON BARYON
C==================================================================
         IF (NNCH2) 631,632,633
631      CONTINUE
            ICH2 = IB8
            GO TO 634
633      CONTINUE
            ICH2 = IBIO
            GO TO 634
632      CONTINUE
            CALL BAMJEV(IHAD,IBF,IBB1,IBB2,IBF,AMCH2,4)
            GO TO 635
 634     CONTINUE
            CALL PARJET(IHAD,ICH2)
 635     CONTINUE
C     CALL DECAY(IHAD,2)
*
*  We arrive here after jet creation: created particles are in
*  /finpar/ common (there are ihad particles)
*
  64  CONTINUE
C*** TURN CHAIN AROUND IF NECESSARY
      IF (IBPROJ) 71,72,73
  71  CONTINUE
         GO TO 75
  72  CONTINUE
         GO TO (721,722),ISAM3
 721     CONTINUE
C*** TURN JET
            GO TO 75
 722     CONTINUE
            GO TO 74
  73  CONTINUE
C*** KEEP JET
         GO TO 74
  75  CONTINUE
C*** TURN JET AROUND
         DO 76 I=1,IHAD
            PZF(I) = -PZF(I)
  76     CONTINUE
  74  CONTINUE
C================================================================
C
C*** TRANSFORM BACKWARD JET INTO CMS
C
C=================================================================
      NAUX   = NAUX+1
 
 
      CALL LORTRA(IHAD,NAUX,GAMCH2,BGCH2X,BGCH2Y,BGCH2Z)
      NAUX = IHAD + NAUX - 1
      DO 181 I=1,NAUX
         PXR(I)  = PXA(I)
         PYR(I)  = PYA(I)
         PZR(I)  = PZA(I)
         AMR(I)  = AMA(I)
         ICHR(I) = ICHA(I)
         ANR(I)  = ANA(I)
         IBARR(I)= IBARA(I)
         NRER(I) = NREA(I)
         HER(I)  = HEPA(I)
 181  CONTINUE
 
      NRES=NAUX
      CALL DECAUX(NAUX,3)
 
 
*or      IF (IPRI.EQ.1) WRITE(LUNOUT,85)(I,NREA(I),ICHA(I),IBARA(I),
*or     &ANA(I)PXA(I),PYA(I),PZA(I),HEPA(I),AMA(I),I=1,NAUX)
 
      EVZ = 0.D0
      PVX = 0.D0
      PVY = 0.D0
      PVZ = 0.D0
      ICCU = 0
      IBBU = 0
      ISSU = 0
      LISSU = .TRUE.
 
C*** TRANSFORM INTO LABSYSTEM
*  +-------------------------------------------------------------------*
*  |   particles from /auxpar/ common  are transformed back in the lab
*  |   system (which is actually the system of the target nucleon with
*  |   the projectile along the z-axis)
*  |   and put in /hadpar/ common
*  |
*  | The transformation is:
*  |     Elab  =  Ecms * gamma + ETAzlab * Pzcms
*  |     Pzlab = Pzcms * gamma + ETAzlab * Ecms
*  |                note ETAzlab = -ETAzcms!!!!
*  |
      DO 81 I=1,NAUX
         HEPH(I) = GAMCM*HEPA(I) + BGCM*PZA(I)
         PZH(I)  = GAMCM*PZA(I)  + BGCM*HEPA(I)
         PXH(I)  = PXA(I)
         PYH(I)  = PYA(I)
         AMH(I)  = AMA(I)
         ICHH(I) = ICHA(I)
         ANH(I)  = ANA(I)
         IBARH(I)= IBARA(I)
         NREH(I) = NREA(I)
         EVZ = EVZ + HEPH(I)
         PVX = PVX + PXH(I)
         PVY = PVY + PYH(I)
         PVZ = PVZ + PZH(I)
         ICCU = ICCU + ICHH(I)
         IBBU = IBBU + IBARH(I)
         IJNREH = KPTOIP ( NREH (I) )
         IF (IJNREH .LE. 0 .OR. IJNREH .GT. 39) THEN
            WRITE (LUNOUT,*)' Hadevt: Ijnreh = 0, > 39 after decay!!',
     &      IJNREH,NREH(I),I,HEPH(I)
            WRITE (LUNERR,*)' Hadevt: Ijnreh = 0, > 39 after decay!!',
     &      IJNREH,NREH(I),I,HEPH(I)
            LISSU = .FALSE.
         ELSE
            DO 8011 J=1,3
               ISSU = ISSU + IQSCHR (MQUARK(J,IJNREH))
8011        CONTINUE
         END IF
  81  CONTINUE
*  |
*  +-------------------------------------------------------------------*
 
      NHAD = NAUX
      ICHTOT = ICH(KPROJ) + ICH(KTARG)
      IBTOT  = IBPROJ + IBTARG
      ISTOT  = 0
      DO 8111 J=1,3
         ISTOT = ISTOT + IQSCHR(MQUARK(J,IJPROJ))
     &         + IQSCHR(MQUARK(J,IJTARG))
8111  CONTINUE
*  +-------------------------------------------------------------------*
*  |
      IF (ICCU .NE. ICHTOT) THEN
*  | write an error message and then resample!!!
         WRITE(LUNOUT,*)' Hadevt: charge conservation failure: ',
     &                  ' ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT',
     &                  ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT
         WRITE(LUNERR,*)' Hadevt: charge conservation failure: ',
     &                  ' ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT',
     &                  ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT
         LRESMP = .TRUE.
      ELSE IF (IBBU .NE. IBTOT) THEN
*  | write an error message and then resample!!!
         WRITE(LUNOUT,*)' Hadevt: baryon conservation failure: ',
     &                  ' ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT',
     &                  ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT
         WRITE(LUNERR,*)' Hadevt: baryon conservation failure: ',
     &                  ' ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT',
     &                  ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT
         LRESMP = .TRUE.
      ELSE IF (ISSU .NE. ISTOT .AND. LISSU) THEN
*  | write an error message and then resample!!!
         WRITE(LUNOUT,*)' Hadevt: strangeness conservation failure: ',
     &                  ' ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT',
     &                  ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT
         WRITE(LUNERR,*)' Hadevt: strangeness conservation failure: ',
     &                  ' ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT',
     &                  ICCU,ICHTOT,IBBU,IBTOT,ISSU,ISTOT
         LRESMP = .TRUE.
      ELSE
      END IF
*  |
*  +-------------------------------------------------------------------*
      EPSEPS = MAX ( 10.D+00*ANGLGB, 1.D-12 )
*  +-------------------------------------------------------------------*
*  |
      IF (ABS(EVZ-AMTAR-EPROJ)/(AMTAR+EPROJ) .GT. EPSEPS) THEN
      ELSE IF (ABS(PVX)/PPROJ .GT. EPSEPS) THEN
      ELSE IF (ABS(PVY)/PPROJ .GT. EPSEPS) THEN
      ELSE IF (ABS(PVZ-PPROJ)/PPROJ .GT. EPSEPS) THEN
      ELSE
         GO TO 90
      END IF
*  |
*  +-------------------------------------------------------------------*
 
C
C********************************************************************
C
C*** PRINT AND TEST ENERGY CONSERVATION
C
C********************************************************************
C
*or      PVZ = 0.D0
*or      EVZ = 0.D0
*or      PVX = 0.D0
*or      PVY = 0.D0
*or      ICCU = 0
*or      IBBU = 0
*or      DO 82 I=1,NHAD
*or         PVX = PVX + PXH(I)
*or         PVY = PVY + PYH(I)
*or         PVZ = PVZ + PZH(I)
*or         EVZ = EVZ + HEPH(I)
*or         ICCU = ICCU + ICHH(I)
*or         IBBU = IBBU + IBARH(I)
*or  82  CONTINUE
*or      IF (IBTOT  .NE. IBBU) GO TO 9999
*or      IF (ICHTOT .NE. ICCU) GO TO 9999
*or      IF (ABS(PVX).GE.0.01D0) GO TO 9999
*or      IF (ABS(PVY).GE.0.01D0) GO TO 9999
*or      IF ((PVZ.GT.1.02D0*PPROJ).OR.(PVZ.LT.0.98D0*PPROJ)) GO TO 9999
*or      IF (IPRI.NE.1) GO TO 90
9999  CONTINUE
*
*  If a failure occured the event is resampled!!!
*
      GO TO 8899
*or      IF (IPRI.EQ.0) GO TO 8899
*or      WRITE(LUNOUT,83)NHAD,KPROJ,KTARG,PPROJ,EPROJ,PVX,PVY,PVZ,EVZ,
*or     &ICCU,IBBU,NHAD
*or  83  FORMAT (3I5,6F12.6,3I5)
*or      DO 84 I=1,NHAD
*or         WRITE(LUNOUT,85)I,NREH(I),ICHH(I),IBARH(I),ANH(I),PXH(I),
*or     &   PYH(I),PZH(I),HEPH(I),AMH(I)
*or  85     FORMAT (4I5,A8,5F12.6)
*or  84  CONTINUE
*
*  If a failure occured the event is resampled!!!
*
*or      IF (IPRI.EQ.0)GO TO 8899
 
  90  CONTINUE
      RETURN
      END
