*CMZ :  2.03/01 20/08/98  17.46.44  by  Federico Carminati
*CMZ :  2.00/01 20/04/98  10.51.01  by  Unknown
*CMZ :  2.00/00 16/04/98  12.54.34  by  Unknown
*CMZ :  1.05/00 08/06/94  09.19.45  by  Nick van Eijndhoven (RUU/CERN)
*-- Author :    Nick van Eijndhoven (RUU/CERN)   08/06/94
      SUBROUTINE SXPART
C
C *** DEFINITION OF ADDITIONAL PARTICLES AND DECAY MODES ***
C *** NVE 08-JUL-1994 SAP UTRECHT ***
C
C CALLED BY : UGINIT
C ORIGIN    : NICK VAN EIJNDHOVEN
C
      DIMENSION BRATIO(6),MODE(6)
C
C --- Define additional particles
      CALL GSPART(33,'OMEGA(783)',3,0.782, 0.,7.836E-23,0.,0)
      CALL GSPART(34,'PHI(1020)' ,3,1.019, 0.,1.486E-22,0.,0)
      CALL GSPART(35,'D +'       ,4,1.870, 1.,1.066E-12,0.,0)
      CALL GSPART(36,'D -'       ,4,1.870,-1.,1.066E-12,0.,0)
      CALL GSPART(37,'D 0'       ,3,1.865, 0.,4.2E-13  ,0.,0)
      CALL GSPART(38,'ANTI D 0'  ,3,1.865, 0.,4.2E-13  ,0.,0)
      CALL GSPART(42,'RHO +'     ,4,0.768, 1.,4.353E-24,0.,0)
      CALL GSPART(43,'RHO -'     ,4,0.768,-1.,4.353E-24,0.,0)
      CALL GSPART(44,'RHO 0'     ,3,0.768, 0.,4.353E-24,0.,0)
C
C --- Define additional decay modes ---
C --- omega(783) ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=33
      BRATIO(1)=89.
      BRATIO(2)=8.5
      BRATIO(3)=2.5
      MODE(1)=70809
      MODE(2)=107
      MODE(3)=908
      CALL GSDK(IPA,BRATIO,MODE)
C --- phi(1020) ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=34
      BRATIO(1)=49.
      BRATIO(2)=34.4
      BRATIO(3)=12.9
      BRATIO(4)=2.4
      BRATIO(5)=1.3
      MODE(1)=1112
      MODE(2)=1610
      MODE(3)=4407
      MODE(4)=90807
      MODE(5)=1701
      CALL GSDK(IPA,BRATIO,MODE)
C --- D+ ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=35
      BRATIO(1)=25.
      BRATIO(2)=25.
      BRATIO(3)=25.
      BRATIO(4)=25.
      MODE(1)=80809
      MODE(2)=120808
      MODE(3)=111208
      MODE(4)=110809
      CALL GSDK(IPA,BRATIO,MODE)
C --- D- ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=36
      BRATIO(1)=25.
      BRATIO(2)=25.
      BRATIO(3)=25.
      BRATIO(4)=25.
      MODE(1)=90908
      MODE(2)=110909
      MODE(3)=121109
      MODE(4)=120908
      CALL GSDK(IPA,BRATIO,MODE)
C --- D0 ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=37
      BRATIO(1)=33.
      BRATIO(2)=33.
      BRATIO(3)=33.
      MODE(1)=809
      MODE(2)=1208
      MODE(3)=1112
      CALL GSDK(IPA,BRATIO,MODE)
C --- Anti D0 ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=38
      BRATIO(1)=33.
      BRATIO(2)=33.
      BRATIO(3)=33.
      MODE(1)=809
      MODE(2)=1109
      MODE(3)=1112
      CALL GSDK(IPA,BRATIO,MODE)
C --- rho+ ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=42
      BRATIO(1)=100.
      MODE(1)=807
      CALL GSDK(IPA,BRATIO,MODE)
C --- rho- ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=43
      BRATIO(1)=100.
      MODE(1)=907
      CALL GSDK(IPA,BRATIO,MODE)
C --- rho0 ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=44
      BRATIO(1)=100.
      MODE(1)=707
      CALL GSDK(IPA,BRATIO,MODE)
C --- jpsi ---
      DO KZ=1,6
         BRATIO(KZ)=0
         MODE(KZ)=0
      ENDDO
      IPA=113
      BRATIO(1)=50.
      BRATIO(2)=50.
      MODE(1)=506
      MODE(2)=605
      CALL GSDK(IPA,BRATIO,MODE)
C --- upsilon ---
      IPA=114
      CALL GSDK(IPA,BRATIO,MODE)
C --- phi ---
      IPA=115
      CALL GSDK(IPA,BRATIO,MODE)
C
      RETURN
      END
*CMZ :          30/12/98  16.23.44  by  Rene Brun
*-- Author :    Rene Brun   30/12/98
      subroutine setclip(name,xmin,xmax,ymin,ymax,zmin,zmax)
*
*       define the clip box parameters for volume name
*
#undef CERNLIB_GEANT321_GCUNIT_INC
#include "geant321/gcunit.inc"
#undef CERNLIB_GEANT321_GCMUTR_INC
#include "geant321/gcmutr.inc"
#undef CERNLIB_GEANT321_GCGOBJ_INC
#include "geant321/gcgobj.inc"
*KEND.
      character *(*) name
      real xmin,xmax,ymin,ymax,zmin,zmax
*
         IHOLE=0
         NCVOLS=NCVOLS+1
         IF(NCVOLS.EQ.MULTRA)THEN
            WRITE(CHMAIL, 10000)
10000      FORMAT(' *** GXDRAW ***:',
     +            ' No more space to store MCVOL information.')
            CALL GMAIL(0,0)
            WRITE(CHMAIL, 10100)
10100      FORMAT(' *** GXDRAW ***: Please reset MCVOL')
            CALL GMAIL(0,0)
            return
         ENDIF
         IF(XMIN.GE.XMAX.OR.YMIN.GE.YMAX.OR.ZMIN.GE.ZMAX)THEN
            WRITE(CHMAIL,10200)
10200       FORMAT(' Wrong Box limits. Check values ')
            CALL GMAIL(0,0)
            return
         ENDIF
****SG
         GNNVV(NCVOLS)=NAME
         GNASH(NCVOLS)='BOX'
         GXMIN(NCVOLS)=XMIN
         GXMAX(NCVOLS)=XMAX
         GYMIN(NCVOLS)=YMIN
         GYMAX(NCVOLS)=YMAX
         GZMIN(NCVOLS)=ZMIN
         GZMAX(NCVOLS)=ZMAX
         IF(GXMIN(NCVOLS).GT.-99999.)IHOLE=1
* Resetting Mcvol mode
         IF(GNNVV(NCVOLS).EQ.'.')THEN
            IHOLE=0
            DO 10 JJ=1,NCVOLS
               GNNVV(JJ)=' '
               GXMIN(JJ)=-100000
               GXMAX(JJ)=-99999
               GYMIN(JJ)=-100000
               GYMAX(JJ)=-99999
               GZMIN(JJ)=-100000
               GZMAX(JJ)=-99999
   10       CONTINUE
            NCVOLS=0
         ENDIF
      end

      subroutine setbomb(boom)
*
*      set the bomb factor
*
#undef CERNLIB_GEANT321_GCMUTR_INC
#include "geant321/gcmutr.inc"
*KEND.
*
      real boom
*
      gboom = boom
      end
