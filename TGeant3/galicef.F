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
