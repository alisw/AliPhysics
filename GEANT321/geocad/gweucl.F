*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GWEUCL (LUN,FILNAM,TOPVOL,NUMBER,NLEVEL)
*
*
*     ******************************************************************
*     *                                                                *
*     *  Write out the geometry of the detector in EUCLID file format  *
*     *                                                                *
*     *       filnam : will be with the extension .euclid              *
*     *       topvol : volume name of the starting node                *
*     *       number : copy number of topvol (relevant for gsposp)     *
*     *       nlevel : number of  levels in the tree structure         *
*     *                to be written out, starting from topvol         *
*     *                                                                *
*     *       Author : M. Maire                                        *
*     *                                                                *
*     ******************************************************************
*
*
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
*
      CHARACTER*(*) FILNAM
      CHARACTER*80  FILEXT
      CHARACTER    CARD*80
      CHARACTER*4  TOPVOL
      CHARACTER*20 NATMED, NAMATE
      CHARACTER*4  NAME, MOTHER, SHAPE(16), KONLY
*
      DIMENSION PAR(100), ATT(20)
*
      DATA SHAPE/'BOX ','TRD1','TRD2','TRAP','TUBE','TUBS','CONE',
     +           'CONS','SPHE','PARA','PGON','PCON','ELTU','HYPE',
     +           'GTRA','CTUB'/
*
*
* *** The output filnam name will be with extension '.euclid'
      IF(INDEX(FILNAM,'.').EQ.0) THEN
         IT=LNBLNK(FILNAM)
      ELSE
         IT=INDEX(FILNAM,'.')-1
      ENDIF
#if !defined(CERNLIB_IBM)
      FILEXT=FILNAM(1:IT)//'.euclid'
#endif
#if defined(CERNLIB_IBM)
      FILEXT='/'//FILNAM(1:MIN(IT,8))//' EUCLID A1'
      CALL CLTOU(FILEXT)
#endif
*
      OPEN (UNIT=LUN,FILE=FILEXT,STATUS='UNKNOWN',FORM='FORMATTED')
*
* *** Initialisation of the working space
      IADVOL = NVOLUM
      IADTMD = IADVOL + NVOLUM
      IADROT = IADTMD + NTMED
      IF(JROTM.GT.0) THEN
         NROTM  = IQ(JROTM-2)
      ELSE
         NROTM = 0
      ENDIF
      NWTOT  = IADROT + NROTM
      CALL GWORK (NWTOT)
      CALL VZERO (IWS(1),NWTOT)
      MLEVEL = NLEVEL
      IF (NLEVEL.LE.0) MLEVEL = 20
*
* *** find the top volume and put it in the stak
      NUMBR = NUMBER
      IF (NUMBER.LE.0) NUMBR = 1
      CALL GFPARA (TOPVOL,NUMBR,1,NPAR,NATT,PAR,ATT)
      IF (NPAR.LE.0) THEN
         WRITE (CHMAIL,11100) TOPVOL,NUMBR
         CALL GMAIL (0,0)
         RETURN
      ENDIF
*
*     authorized shape ?
      CALL GLOOK (TOPVOL,IQ(JVOLUM+1),NVOLUM,IVO)
      JVO = LQ(JVOLUM - IVO)
      ISH =  Q(JVO + 2)
      IF (ISH.GT.12) THEN
         WRITE (CHMAIL,11100) TOPVOL,NUMBR
         CALL GMAIL (0,0)
         RETURN
      ENDIF
*
      LEVEL  = 1
      NVSTAK = 1
      IWS(NVSTAK)     = IVO
      IWS(IADVOL+IVO) = LEVEL
      IVSTAK = 0
*
* *** Flag all volumes and fill the stak
*
   10 CONTINUE
*
*     pick the next volume in stak
      IVSTAK = IVSTAK + 1
      IVO   = ABS(IWS(IVSTAK))
      JVO   = LQ(JVOLUM - IVO)
*
*     flag the tracking medium
      NUMED =  Q(JVO + 4)
      IWS(IADTMD + NUMED) = 1
*
*     get the daughters ...
      LEVEL = IWS(IADVOL+IVO)
      IF (LEVEL.LT.MLEVEL) THEN
         LEVEL = LEVEL + 1
         NIN = Q(JVO + 3)
*
*        from division ...
         IF (NIN.LT.0) THEN
            JDIV = LQ(JVO  - 1)
            IVIN =  Q(JDIV + 2)
            NVSTAK = NVSTAK + 1
            IWS(NVSTAK)      = -IVIN
            IWS(IADVOL+IVIN) =  LEVEL
*
*        from position ...
         ELSE IF (NIN.GT.0) THEN
            DO 20 IN=1,NIN
               JIN  = LQ(JVO - IN)
               IVIN =  Q(JIN + 2 )
               JVIN = LQ(JVOLUM - IVIN)
               ISH  =  Q(JVIN + 2)
*              authorized shape ?
               IF (ISH.LE.12) THEN
*                 not yet flagged ?
                  IF (IWS(IADVOL+IVIN).EQ.0) THEN
                     NVSTAK = NVSTAK + 1
                     IWS(NVSTAK)      = IVIN
                     IWS(IADVOL+IVIN) = LEVEL
                  ENDIF
*                 flag the rotation matrix
                  IROT =  Q(JIN + 4 )
                  IF (IROT.GT.0) IWS(IADROT+IROT) = 1
               ENDIF
   20       CONTINUE
         ENDIF
      ENDIF
*
*     next volume in stak ?
      IF (IVSTAK.LT.NVSTAK) GO TO 10
*
* *** Write down the tracking medium definition
*
      CARD = '!       Tracking medium'
      WRITE (LUN,10000) CARD
*
      DO 30 ITM = 1,NTMED
         IF (IWS(IADTMD+ITM).GT.0) THEN
            JTM  = LQ(JTMED-ITM)
            CALL UHTOC (IQ(JTM+1),4,NATMED,20)
            IMAT =  Q(JTM+6)
            JMA  = LQ(JMATE-IMAT)
            IF(JMA.LE.0) THEN
               NAMATE = ' '
               WRITE(CHMAIL,11300) ITM, NATMED(1:LNBLNK(NATMED))
               CALL GMAIL(1,1)
            ELSE
               CALL UHTOC (IQ(JMA+1),4,NAMATE,20)
            ENDIF
            CARD = ' '
            WRITE (CARD,10100) ITM,NATMED,IMAT,NAMATE
            WRITE (LUN,'(A)') CARD
         ENDIF
   30 CONTINUE
*
* *** Write down the rotation matrix
*
      CARD = '!       Reperes'
      WRITE (LUN,10000) CARD
*
      DO 40 IRM = 1,NROTM
         IF (IWS(IADROT+IRM).GT.0) THEN
            JRM  = LQ(JROTM-IRM)
            CARD = ' '
            WRITE (CARD,10200) IRM,(Q(JRM+K),K=11,16)
            WRITE (LUN,'(A)') CARD
         ENDIF
   40 CONTINUE
*
* *** Write down the volume definition
*
      CARD = '!       Volumes'
      WRITE (LUN,10000) CARD
*
      DO 50 IVSTAK = 1,NVSTAK
         IVO = IWS(IVSTAK)
         IF (IVO.GT.0) THEN
            CALL UHTOC (IQ(JVOLUM+IVO),4,NAME,4)
            JVO  = LQ(JVOLUM-IVO)
            ISH   = Q(JVO+2)
            NMED  = Q(JVO+4)
            IF (IVSTAK.GT.1) NPAR  = Q(JVO+5)
            CARD = ' '
            IF (NPAR.GT.0) THEN
               IF (IVSTAK.GT.1) CALL UCOPY (Q(JVO+7),PAR(1),NPAR)
               CALL GCKPAR (ISH,NPAR,PAR)
               WRITE (CARD,10300) NAME,SHAPE(ISH),NMED,NPAR
               WRITE (LUN,'(A)') CARD
               WRITE (LUN,10400) (PAR(K),K=1,NPAR)
            ELSE
               WRITE (CARD,10300) NAME,SHAPE(ISH),NMED,NPAR
               WRITE (LUN,'(A)') CARD
            ENDIF
         ENDIF
   50 CONTINUE
*
* *** Write down the division of volumes
*
      CARD = '!       Divisions'
      WRITE (LUN,10000) CARD
*
      DO 60 IVSTAK = 1,NVSTAK
         IVO = ABS(IWS(IVSTAK))
         JVO  = LQ(JVOLUM-IVO)
         ISH  =  Q(JVO+2)
         NIN  =  Q(JVO+3)
*        this volume is divided ...
         IF (NIN.LT.0) THEN
            JDIV = LQ(JVO-1)
            IAXE =  Q(JDIV+1)
            IVIN =  Q(JDIV+2)
            NDIV =  Q(JDIV+3)
            C0   =  Q(JDIV+4)
            STEP =  Q(JDIV+5)
            JVIN = LQ(JVOLUM-IVIN)
            NMED =  Q(JVIN+4)
            CALL UHTOC (IQ(JVOLUM+IVO ),4,MOTHER,4)
            CALL UHTOC (IQ(JVOLUM+IVIN),4,NAME  ,4)
            CARD = ' '
            IF ((STEP.LE.0.).OR.(ISH.GE.11)) THEN
*              volume with negative parameter or gsposp or PGON ...
               WRITE (CARD,10500) NAME,MOTHER,NDIV,IAXE
            ELSEIF ((NDIV.LE.0).OR.(ISH.EQ.10)) THEN
*              volume with negative parameter or gsposp or PARA ...
               NDVMX = ABS(NDIV)
               WRITE (CARD,10600) NAME,MOTHER,STEP,IAXE,NMED,NDVMX
            ELSE
*              normal volume : all kind of division are equivalent
               WRITE (CARD,10700) NAME,MOTHER,STEP,IAXE,C0,NMED,NDIV
            ENDIF
            WRITE (LUN,'(A)') CARD
         ENDIF
   60 CONTINUE
*
* *** Write down the the positionnement of volumes
*
      card = '!       Positionnements'
      WRITE (LUN,10000) CARD
*
      DO 80 IVSTAK = 1,NVSTAK
         IVO = ABS(IWS(IVSTAK))
         CALL UHTOC (IQ(JVOLUM+IVO ),4,MOTHER,4)
         JVO  = LQ(JVOLUM-IVO)
         NIN  =  Q(JVO+3)
*        this volume has daughters ...
         IF (NIN.GT.0) THEN
            DO 70 IN=1,NIN
               JIN  = LQ(JVO-IN)
               IVIN =  Q(JIN +2)
               NUMB =  Q(JIN +3)
               IROT =  Q(JIN +4)
               X    =  Q(JIN +5)
               Y    =  Q(JIN +6)
               Z    =  Q(JIN +7)
               KONLY = 'ONLY'
               IF (Q(JIN+8).NE.1.) KONLY = 'MANY'
               CALL UHTOC (IQ(JVOLUM+IVIN),4,NAME  ,4)
               JVIN = LQ(JVOLUM-IVIN)
               ISH  =  Q(JVIN+2)
               CARD = ' '
*              gspos or gsposp ?
               NDATA = IQ(JIN-1)
               IF (NDATA.EQ.8) THEN
                  WRITE (CARD,10800) NAME,NUMB,MOTHER,X,Y,Z,IROT,KONLY
                  WRITE (LUN,'(A)') CARD
               ELSE
                  NPAR =  Q(JIN+9)
                  CALL UCOPY (Q(JIN+10),PAR(1),NPAR)
                  CALL GCKPAR (ISH,NPAR,PAR)
                  WRITE (CARD,10900) NAME,NUMB,MOTHER,X,Y,Z,IROT,KONLY,
     +            NPAR
                  WRITE (LUN,'(A)') CARD
                  WRITE (LUN,10400) (PAR(K),K=1,NPAR)
               ENDIF
   70       CONTINUE
         ENDIF
   80 CONTINUE
*
      WRITE (LUN,11000)
      CLOSE (LUN)
*
      WRITE (CHMAIL,11200) FILEXT(1:IT+9)
      CALL GMAIL (1,1)
*
10000 FORMAT (1H!,/,A,/,1H!)
*
10100 FORMAT ('TMED',2(1X,I3,1X,1H',A20,1H'))
10200 FORMAT ('ROTM',1X,I3,6(1X,F8.3))
10300 FORMAT ('VOLU',2(1X,1H',A4,1H'),2(1X,I3))
10400 FORMAT (      (5X,6(1X,F11.5)))
10500 FORMAT ('DIVN',2(1X,1H',A4,1H'),2(1X,I3))
10600 FORMAT ('DIVT',2(1X,1H',A4,1H'),1X,F11.5,3(1X,I3))
10700 FORMAT ('DVT2',2(1X,1H',A4,1H'),1X,F11.5,1X,I3,1X,F11.5,2(1X,I3))
10800 FORMAT ('POSI',1X,1H',A4,1H',1X,I3,1X,1H',A4,1H',3(1X,F11.5),1X,I3
     &              ,1X,1H',A4,1H')
10900 FORMAT ('POSP',1X,1H',A4,1H',1X,I3,1X,1H',A4,1H',3(1X,F11.5),1X,I3
     &              ,1X,1H',A4,1H',1X,I3)
11000 FORMAT ('END')
*
11100 FORMAT (' *** GWEUCL *** top volume : ',A4,' number :',I3,
     &        ' can not be a valid root')
11200 FORMAT (' *** GWEUCL *** file: ',A,' is now written out')
11300 FORMAT (' *** GWEUCL *** material not defined for tracking ',
     +        'medium ',I5,' ',A)
*
      END
