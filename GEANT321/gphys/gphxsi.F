*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:31  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GPHXSI
*.
*.    ******************************************************************
*.    *  Creates PHXS bank containing x-section constants              *
*.    *                                                                *
*.    *    ==>CALLED BY : GPHYSI                                       *
*.    *       AUTHOR : J. CHWASTOWSKI                                  *
*.    *                                                                *
*.    ******************************************************************
*.
*      The photoelectric effect x-section bank
*       1                 - Number of elements neded to create medium = NZ
*       2      <-> NZ+1   -  Z of elements
*       NZ+1   <-> 2*NZ+1 -  0
*       2+2*NZ <-> 3*NZ+1 -  weight of x-section constants
*       3*NZ+2 <-> top    -  X-section constants
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcpmxz.inc"
#include "geant321/gcphxs.inc"
#include "geant321/gc10ev.inc"
      CHARACTER*1 CHSHEL, CHGROU
      DIMENSION EUP(MAXINT),TMP(MAXPOW,MAXINT),ESHL(24)
      PARAMETER (NXSB=73,NXSBF=70)
      Z = Q(JMA+7)
      IPHXSP = 0
* Check Z range of validity for Sandia parametrization
      IF(Z.GE.1.AND.Z.LE.MAXELZ) THEN
*
* Find number of elements neded to create current medium
c
* Is this medium a mixture ?
         NMIX = Q(JMA+11)
         IF(NMIX.GT.1) THEN
            NZ = 0
            JMIXT=LQ(JMA-5)
            DO 10 I = 1,NMIX
               ZCUR = Q(JMIXT+NMIX+I)
               IF(ZCUR.NE.INT(ZCUR)) THEN
*
* When Z is non integer you need 2 elements
                  NZ=NZ+2
               ELSE
                  NZ=NZ+1
               ENDIF
   10       CONTINUE
            CALL GWORK(3*NZ+1)
            WS(1) = NZ
*
* Calculate weigths
*
            K = 1
            DO 20 I = 1,NMIX
               ZCUR = Q(JMIXT+NMIX+I)
               HZCUR = INT(ZCUR)
               WS(1+K) = HZCUR
               IF(ZCUR.NE.HZCUR) THEN
                  WS(1+2*NZ+K) = (HZCUR+1.-ZCUR)*Q(JMIXT+2*NMIX+I)
                  K = K+1
                  WS(1+K) = HZCUR+1
                  WS(1+2*NZ+K) = (ZCUR-HZCUR)*Q(JMIXT+2*NMIX+I)
               ELSE
                  WS(1+2*NZ+K) = Q(JMIXT+2*NMIX+I)
               ENDIF
               K = K+1
   20       CONTINUE
*
* Do Z values repeat ?
*
            K = NZ
            DO 40 I = 1,NZ-1
               Z1 = WS(1+I)
               IF(Z1.GT.0.0) THEN
                  DO 30 J = I+1,NZ
                     Z2 = WS(1+J)
                     IF(Z1.EQ.Z2) THEN
*                        WS(1+NZ+I) = WS(1+NZ+I)+WS(1+NZ+J)
                        WS(1+2*NZ+I) = WS(1+2*NZ+I)+WS(1+2*NZ+J)
                        WS(1+J) = -WS(1+J)
                        K = K-1
                     ENDIF
   30             CONTINUE
               ENDIF
   40       CONTINUE
* Now you can book a pht. eff. x-sec. constant bank and hang it as
* a first struc. link from JPHOT.
* From this bank you hang the banks for each separate Z !!!
            NW = NZ*NXSB+2
            IF(K.EQ.1) NW = NXSB+2
            CALL MZBOOK(IXCONS,JPHXS,JPHOT,-1,'PHXS',NZ,NZ,NW,3,0)
* Fill JPHXS bank and calculates the weights
            Q(JPHXS+1) = K
            NZN = K
            K = 1
            DO 50 I = 1,NZ
               IF(WS(1+I).GT.0.0) THEN
                  Q(JPHXS+1+K) = WS(1+I)
*                  Q(JPHXS+1+K+NZN) = WS(1+I+NZ)
                  Q(JPHXS+1+K+2*NZN) = WS(1+I+2*NZ)
                  K = K+1
               ENDIF
   50       CONTINUE
         ELSE
* Current medium consists of one "element"
            NZ = 1
            ZCUR = Z
            HZCUR = INT(ZCUR)
            IF(MOD(ZCUR,HZCUR).EQ.0) THEN
* Now you can book a pht. eff. x-sec. constant bank and hang it as
* a first struc. link from JPHOT.
               NW = NXSB+2
               CALL MZBOOK(IXCONS,JPHXS,JPHOT,-1,'PHXS',1,1,NW,3,0)
               Q(JPHXS+1) = 1.
               Q(JPHXS+2) = ZCUR
               Q(JPHXS+4) = 1.
            ELSE
* Somebody is cheating. We need two elements
               NZ = NZ+1
               NW = 2*NXSB+2
               CALL MZBOOK(IXCONS,JPHXS,JPHOT,-1,'PHXS',2,2,NW,3,0)
               Q(JPHXS+1) = 2.
               Q(JPHXS+2) = HZCUR
               Q(JPHXS+3) = HZCUR+1
               Q(JPHXS+6) = (HZCUR+1.-ZCUR)
               Q(JPHXS+7) = (ZCUR-HZCUR)
            ENDIF
         ENDIF
* We passed the bank booking phase and we can start real work
         NZ = Q(JPHXS+1)
* Create a temporary working space
         IBINS = 5*NZ+NZ*(MAXPOW+1)*(MAXINT+1)
         CALL GWORK(IBINS)
         DO 60 I = 1,IBINS
            WS(I) = 0.0
   60    CONTINUE
         KS = 5*NZ
         L = 0
         DO 170 N = 1,NZ
            JPHXS = LQ(JPHOT-1)
            ZCUR = Q(JPHXS+1+N)
            IZ = ZCUR
            DO 70 I = 1,24
               ESHL(I) = 0.0
   70       CONTINUE
            CALL GFSHLS(ZCUR,ESHL,NSHL)
            DO 80 I = 1,NSHL
               ESHL(I) = ESHL(I)*1.E-3
   80       CONTINUE
* Use Sandia data.
* Find out the interval upper limit.
            IMAX = 0
            DO 90 I = 1,MAXINT
               CHSHEL=CRNGUP(I,IZ)(1:1)
               CHGROU=CRNGUP(I,IZ)(2:2)
               IF(CHSHEL.EQ.'I') THEN
                  EUP(I) = 55.E11
                  IMAX = I
                  GO TO 100
               ELSEIF(CHSHEL.EQ.'K') THEN
                  EUP(I) = ESHL(1)
               ELSEIF(CHSHEL.EQ.'L') THEN
                  IF(CHGROU.EQ.'1') THEN
                     EUP(I) = ESHL(2)
                  ELSEIF(CHGROU.EQ.'2') THEN
                     EUP(I) = ESHL(3)
                  ELSEIF(CHGROU.EQ.'3') THEN
                     EUP(I) = ESHL(4)
                  ELSE
                     CHMAIL = ' GPHXSI Error: Inconsistent data for L '
     +               //'shell.  Maximum coded L3. Found shell '
     +               //'name: '//CRNGUP(I,IZ)
                     CALL GMAIL(0,0)
                  ENDIF
               ELSEIF(CHSHEL.EQ.'M') THEN
                  IF(CHGROU.EQ.'1') THEN
                     EUP(I) = ESHL(5)
                  ELSEIF(CHGROU.EQ.'2') THEN
                     EUP(I) = ESHL(6)
                  ELSEIF(CHGROU.EQ.'3') THEN
                     EUP(I) = ESHL(7)
                  ELSEIF(CHGROU.EQ.'4') THEN
                     EUP(I) = ESHL(8)
                  ELSEIF(CHGROU.EQ.'5') THEN
                     EUP(I) = ESHL(9)
                  ELSE
                     CHMAIL = ' GPHXSI Error: Inconsistent data for M '
     +               //'shell.  Maximum coded M5. Found shell '
     +               //'name: '//CRNGUP(I,IZ)
                     CALL GMAIL(0,0)
                  ENDIF
               ELSEIF(CHSHEL.EQ.'N') THEN
                  IF(CHGROU.EQ.'1') THEN
                     EUP(I) = ESHL(10)
                  ELSEIF(CHGROU.EQ.'2') THEN
                     EUP(I) = ESHL(11)
                  ELSEIF(CHGROU.EQ.'3') THEN
                     EUP(I) = ESHL(12)
                  ELSEIF(CHGROU.EQ.'4') THEN
                     EUP(I) = ESHL(13)
                  ELSEIF(CHGROU.EQ.'5') THEN
                     EUP(I) = ESHL(14)
                  ELSEIF(CHGROU.EQ.'6') THEN
                     EUP(I) = ESHL(15)
                  ELSEIF(CHGROU.EQ.'7') THEN
                     EUP(I) = ESHL(16)
                  ELSE
                     CHMAIL = ' GPHXSI Error: Inconsistent data for N '
     +               //'shell. Maximum coded N7. Found shell name:'
     +               //' '//CRNGUP(I,IZ)
                     CALL GMAIL(0,0)
                  ENDIF
               ELSEIF(CHSHEL.EQ.'O') THEN
                  IF(CHGROU.EQ.'1') THEN
                     EUP(I) = ESHL(17)
                  ELSEIF(CHGROU.EQ.'2') THEN
                     EUP(I) = ESHL(18)
                  ELSEIF(CHGROU.EQ.'3') THEN
                     EUP(I) = ESHL(19)
                  ELSEIF(CHGROU.EQ.'4') THEN
                     EUP(I) = ESHL(20)
                  ELSEIF(CHGROU.EQ.'5') THEN
                     EUP(I) = ESHL(21)
                  ELSE
                     CHMAIL = ' GPHXSI Error: Inconsistent data for O '
     +               //'shell.  Maximum code O5. Found shell name:'
     +               //' '//CRNGUP(I,IZ)
                     CALL GMAIL(0,0)
                  ENDIF
               ELSEIF(CHSHEL.EQ.'P') THEN
                  IF(CHGROU.EQ.'1') THEN
                     EUP(I) = ESHL(22)
                  ELSEIF(CHGROU.EQ.'2') THEN
                     EUP(I) = ESHL(23)
                  ELSEIF(CHGROU.EQ.'3') THEN
                     EUP(I) = ESHL(24)
                  ELSE
                     CHMAIL = ' GPHXSI Error: Inconsistent data for P '
     +               //'shell.  Maximum coded P3. Found shell '
     +               //'name: '//CRNGUP(I,IZ)
                     CALL GMAIL(0,0)
                  ENDIF
               ELSE
                  READ(CRNGUP(I,IZ),'(F6.0)') EUP(I)
               ENDIF
               IF(EUP(I).EQ.0.0) THEN
                  WRITE(CHMAIL,'(A44,I5)') ' GPHXSI Error: Upper limit '
     +            //'= 0. Interval:',I
                  CALL GMAIL(0,0)
                  STOP 14
               ENDIF
   90       CONTINUE
  100       CONTINUE
            DO 120 I = 1,IMAX
               DO 110 J = 1,MAXPOW
                  TMP(J,I) = COFS(J,I,IZ)*Q(JPHXS+1+2*NZ+N)
  110          CONTINUE
  120       CONTINUE
* Copy upper limits and coofficients to a work bank
            K = KS+1
            WS(K) = GPOMIN(IZ)
            WS(K+1) = 0.0
            WS(K+2) = 0.0
            WS(K+3) = 0.0
            WS(K+4) = 0.0
            K = K+5
            DO 140 I = 1,IMAX
               WS(K) = EUP(I)
               K = K+1
               DO 130 J = 1,MAXPOW
                  WS(K) = TMP(J,I)
                  K = K+1
  130          CONTINUE
  140       CONTINUE
            KS = KS+(MAXINT+1)*(MAXPOW+1)
            IWS(NZ+N) = IMAX
            IWS(2*NZ+N) = 0
            IWS(3*NZ+N) = 0
*            WS(N) = Q(JPHXS+1+2*NZ+N)
* Create element x-section & final state bank
            CALL MZBOOK(IXCONS,JPHFN,JPHXS,-N,'PHFN',0,0,IMAX*5+1,3,0)
            Q(JPHFN+1) = IMAX
* Update pointer
            JPHFN = JPHFN+1
            KFN = 1
* Copy energy & x-section parameters
            DO 160 I = 1,IMAX
               Q(JPHFN+KFN) = EUP(I)
               KFN = KFN+1
               DO 150 J = 1,MAXPOW
c                  Q(JPHFN+KFN) = TMP(J,I)*Q(JPHXS+1+2*NZ+N)
                  Q(JPHFN+KFN) = TMP(J,I)
                  KFN = KFN+1
  150          CONTINUE
  160       CONTINUE
* Get shells decay parameters
            CALL GFSHDC(N,ZCUR)
  170    CONTINUE
* Now find intervals and calculate the coofs for each
         K = 0
         JPHXS = LQ(JPHOT-1)
         IF(NZ.LT.2) THEN
* Simple element so life is easy and nice
            Q(JPHXS+5) = IMAX+1
            JPHXS6 = JPHXS+6
            Q(JPHXS6) = GPOMIN(IZ)
            Q(JPHXS6+1) = 0.0
            Q(JPHXS6+2) = 0.0
            Q(JPHXS6+3) = 0.0
            Q(JPHXS6+4) = 0.0
            JPHXS6 = JPHXS6+5
            DO 190 I = 1,IMAX
               Q(JPHXS6+K) = EUP(I)
               K = K+1
               DO 180 J = 1,MAXPOW
                  Q(JPHXS6+K) = TMP(J,I)
                  K = K+1
  180          CONTINUE
  190       CONTINUE
            NPSH = NW-5*(IMAX+1)-5
         ELSE
* More elements. It will not be so easy
* The following code is difficult and probably there are better solutions
* but I could think only about this one.
            IPHXSP = JPHXS+2+3*NZ
            IPIMAX = JPHXS+2+3*NZ
            IOFFST = (MAXPOW+1)*(MAXINT+1)
            DO 250 II = 1,NZ*(MAXINT+1)
* Find the interval for which the upper limit is the smallest
               AMINV = 1.E20
               DO 200 I = 1,NZ
                  K = IWS(2*NZ+I)
                  IPOINT = 5*NZ+1+(I-1)*IOFFST+K
                  IEUP = 4*NZ+I
                  WS(IEUP) = WS(IPOINT)
                  IF(WS(IEUP).LT.AMINV) AMINV = WS(IEUP)
  200          CONTINUE
               L = 0
               DO 210 I = 1,NZ
                  IF(WS(4*NZ+I).LE.AMINV) L = L+1
  210          CONTINUE
               IF(L.LT.1) THEN
                  CHMAIL = ' GPHXSI Error: Zero intervals found.'
                  CALL GMAIL(0,0)
                  STOP 16
               ENDIF
* Copy to JPHXS bank
               Q(IPIMAX) = Q(IPIMAX)+1.
               Q(IPHXSP+1) = AMINV
               IPHXSP = IPHXSP+1
               DO 230 J = 1,MAXPOW
                  QS = 0.0
                  DO 220 I = 1,NZ
                     K = IWS(2*NZ+I)
                     IPOINT = 5*NZ+1+(I-1)*IOFFST+K+J
c                     QS = QS+WS(I)*WS(IPOINT)
                     QS = QS+WS(IPOINT)
  220             CONTINUE
                  Q(IPHXSP+1) = QS
                  IPHXSP = IPHXSP+1
  230          CONTINUE
               IF(L.EQ.NZ.AND.AMINV.EQ.55.E11) GO TO 260
* Update local pointers
               DO 240 I = 1,NZ
                  IEUP = 4*NZ+I
                  IF(WS(IEUP).LE.AMINV) THEN
                     INZ2 = 2*NZ+I
                     IF(IWS(3*NZ+I).LT.IWS(I+NZ)) THEN
                        IWS(INZ2) = IWS(INZ2)+5
                        IWS(3*NZ+I) = IWS(3*NZ+I)+1
                     ENDIF
                  ENDIF
  240          CONTINUE
  250       CONTINUE
* It is THE END of th x-secs. part, however you may not believe it.
  260       CONTINUE
            NIT=Q(IPIMAX)
            IIT = -1
  261       IIT = IIT+1
*
* It just may happen that some of the energy limits are the same
               IENE=IIT*5+1
               IF(ABS(Q(IPIMAX+IENE)-Q(IPIMAX+IENE+5)).LT.
     +                                Q(IPIMAX+IENE)*5E-4) THEN
                  DO 262 II=1,(NIT-IIT-1)*5
                     Q(IPIMAX+IENE+II)=Q(IPIMAX+IENE+II+5)
  262             CONTINUE
                  NIT=NIT-1
                  IIT=IIT-1
               ENDIF
            IF(IIT.LT.NIT-1) GO TO 261
            Q(IPIMAX)=NIT
            NPSH = NW-Q(IPIMAX)*5-3*NZ-2
         ENDIF
* Return unused locations in JPHXSI bank to the system
         IF(NPSH.GT.0) THEN
            JPHXS = LQ(JPHOT-1)
            CALL MZPUSH(IXCONS,JPHXS,0,-NPSH,'R')
         ENDIF
      ELSEIF(Z.GT.100.) THEN
* Just in case we got called
         CHMAIL = ' GPHXSI Error: Z > 100. No Sandia parameters. '
         CALL GMAIL(0,0)
      ENDIF
      END
