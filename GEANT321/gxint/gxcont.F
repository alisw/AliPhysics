*CMZ :          05/02/96  17.20.51  by  S.Ravndal
*CMZ :  3.21/02 07/07/94  19.47.38  by  S.Giani
*-- Author :
      SUBROUTINE GXCONT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      General control commands                                  *
C.    *                                                                *
C.    *       Authors:   R.Brun      **********                        *
C.    *                  P.Zanarini  **********                        *
C.    *                  S.Giani     **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/pawc.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gclist.inc"
#include "geant321/gcstak.inc"
#include "geant321/gcdraw.inc"
      CHARACTER*20 CHNAME
      CHARACTER*4  CHMECA,CHDEB,CHSET,CHDET,CHMECS(20)
      CHARACTER*32 CHPATL
      CHARACTER*32 CHPAR,NAMATE,NAPAR
      CHARACTER*32 UBUF,WMAT,AA,ZZ
      CHARACTER*32 BRATIO,MODE
C.
C.    ------------------------------------------------------------------
C.
      IWKSTY = IGIWTY(1)
      CALL KUPATL(CHPATL,NPAR)
*
      IF (CHPATL.EQ.'KINE') THEN
         CALL KUGETI(IKINE)
         NPK=NPAR-1
         DO 10 I=1,NPK
            CALL KUGETR(PKINE(I))
   10    CONTINUE
*
      ELSEIF (CHPATL.EQ.'SORD') THEN
         CALL KUGETI(ISTORD)
*
*     ELSEIF (CHPATL.EQ.'PATR') THEN
*        CALL KUGETI(NJTMAX)
*        CALL KUGETI(NJTMIN)
*        WRITE(CHMAIL,'('' *** GXCONT: Parallel tracking disabled, '',
*    +                  ''see documentation'')')
*        CALL GMAIL(1,1)
*
      ELSEIF (CHPATL.EQ.'RUNG') THEN
         CALL KUGETI(IDRUN)
         CALL KUGETI(IDEVT)
*
      ELSEIF (CHPATL.EQ.'GTIME') THEN
         CALL KUGETR(TIMINT)
         CALL KUGETR(TIMEND)
         CALL KUGETI(ITIME)
*
      ELSEIF (CHPATL.EQ.'TRACK') THEN
         IF(JXYZ.NE.0)CALL MZDROP(IXDIV,JXYZ,' ')
         IF(JHITS.NE.0)CALL MZDROP(IXDIV,JHITS,' ')
         CALL MZGARB(IXDIV,0)
         CALL GTREVE
*
      ELSEIF (CHPATL.EQ.'TRIGGER') THEN
         CALL KUGETI(N)

         IF( NEVENT .EQ. 10000000 .OR. NEVENT .NE. IEVENT) THEN
            NEVENT = N
         ELSE
            NEVENT = NEVENT + N
         ENDIF

 20      CONTINUE

         IF( IEVENT .LT. NEVENT ) THEN

            CALL GTRIGC
            IEVENT=IEVENT+1
            CALL GTRIGI
            CALL GTRIG
            CALL GRNDMQ(IQ(JRUNG+19),IQ(JRUNG+20),0,'G')
            IQ(JRUNG+29)=IEVENT

            GOTO 20

         ENDIF

*         DO 20 I=1,N
*            CALL GTRIGC
*            IEVENT=IEVENT+1
*            CALL GTRIGI
*            CALL GTRIG
*            CALL GRNDMQ(IQ(JRUNG+19),IQ(JRUNG+20),0,'G')
*            IQ(JRUNG+29)=IEVENT
*   20    CONTINUE
*
      ELSEIF (CHPATL.EQ.'SWITCH') THEN
         CALL KUGETI(I)
         IF(I.GE.1.AND.I.LE.10) THEN
            IF(NPAR.LT.2)THEN
               WRITE(CHMAIL,'('' ISWIT('',I2,'')= '',I3)')I,ISWIT(I)
               CALL GMAIL(0,0)
            ENDIF
            CALL KUGETI(ISWIT(I))
         ENDIF
*
      ELSEIF (CHPATL.EQ.'PRINT') THEN
         CALL KUGETC(CHNAME,NCH)
         N=0
         CALL KUGETI(N)
         CALL GPRINT(CHNAME,N)
*
      ELSEIF (CHPATL.EQ.'DEBUG') THEN
         CALL KUGETC(CHDEB,NCH)
         IF(CHDEB.EQ.'ON')THEN
            IDEBUG=1
            IDEMIN=1
            IDEMAX=1000000
            ITIME=1
         ELSE
            IDEBUG=0
            IDEMIN=0
            IDEMAX=0
         ENDIF
*
      ELSEIF (CHPATL.EQ.'OUTPUT_LP') THEN
         CALL KUGETI(LOUT)
*
      ELSEIF (CHPATL.EQ.'PHITS') THEN
         CALL KUGETC(CHSET,NCH)
         CALL KUGETC(CHDET,NCH)
         CALL KUGETI(NUMHIT)
         IPKHIT=NUMHIT
         CALL GPHITS(CHSET,CHDET)
         IPKHIT=0
*
      ELSEIF (CHPATL.EQ.'PDIGI') THEN
         CALL KUGETC(CHSET,NCH)
         CALL KUGETC(CHDET,NCH)
         CALL GPDIGI(CHSET,CHDET)
*
      ELSEIF (CHPATL.EQ.'MZLOGL') THEN
         IXSTOR=0
         LEVEL=0
         CALL KUGETI(LEVEL)
         CALL MZLOGL(IXSTOR,LEVEL)
*
      ELSEIF (CHPATL.EQ.'PMATE') THEN
         CALL KUGETI(NUMB)
         IF(IWKSTY.GE.1.AND.IWKSTY.LE.10) THEN
            CALL GPMATX(NUMB)
         ELSE
            CALL GPMATE(NUMB)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'SMATE') THEN
         CALL KUGETI(IMAT)
         CALL KUGETC(NAMATE,NCH)
         CALL KUGETR(A)
         CALL KUGETR(Z)
         CALL KUGETR(DENS)
         CALL KUGETR(RADL)
         CALL KUGETR(ABSL)
         CALL KUGETV(UBUF,LLOW,LHIGH)
         CALL KUGETI(NWBUF)
         CALL GSMATE(IMAT,NAMATE,A,Z,DENS,RADL,ABSL,QQ(LLOW),NWBUF)
*
      ELSEIF (CHPATL.EQ.'SMIXT') THEN
         CALL KUGETI(IMAT)
         CALL KUGETC(NAMATE,NCH)
         CALL KUGETV(AA,LLOW,LHIGH)
         CALL KUGETV(ZZ,LLO1,LHIG1)
         CALL KUGETR(DENS)
         CALL KUGETI(NLMAT)
         CALL KUGETV(WMAT,LLO2,LHIG2)
         CALL GSMIXT(IMAT,NAMATE,QQ(LLOW),QQ(LLO1),DENS,NLMAT,QQ(LLO2))
*
      ELSEIF (CHPATL.EQ.'PRMAT') THEN
         CALL KUGETI(IMATE)
         CALL KUGETI(IPART)
         CALL KUGETC(CHMECA,NCH)
         CALL GPRMAT(IMATE,IPART,CHMECA,NEKBIN,ELOW)
*
      ELSEIF (CHPATL.EQ.'DRMAT') THEN
         CALL KUGETI(IMATE)
         CALL KUGETI(IPART)
         IF(NPAR.GT.2) THEN
            NMEC = NPAR-2
            DO 11 KMEC=1,NMEC
               CALL KUGETC(CHMECS(KMEC),NCH)
  11        CONTINUE
         ELSE
            NMEC = 1
            CHMECS(1)='ALL'
         ENDIF
         CALL GDRMAT(IMATE,IPART,CHMECS,NMEC)
*
      ELSEIF (CHPATL.EQ.'PLMAT') THEN
         CALL KUGETI(IMATE)
         CALL KUGETI(IPART)
         CALL KUGETC(CHMECA,NCH)
         IDM=0
         CALL KUGETI(IDM)
         CALL GPLMAT(IMATE,IPART,CHMECA,NEKBIN,ELOW,IDM)
*
      ELSEIF (CHPATL.EQ.'SPART') THEN
         CALL KUGETI(IPART)
         CALL KUGETC(NAPAR,NCH)
         CALL KUGETI(ITRTYP)
         CALL KUGETR(AMASS)
         CALL KUGETR(CHARGE)
         CALL KUGETR(TLIFE)
         CALL KUGETV(UBUF,LLOW,LHIGH)
         CALL KUGETI(NWBUF)
         CALL KUGETV(BRATIO,LLO1,LHI1)
         CALL KUGETV(MODE,LLO2,LHI2)
         CALL GSPART(IPART,NAPAR,ITRTYP,AMASS,CHARGE,TLIFE,
     +               QQ(LLOW),NWBUF)
         CALL GSDK(IPART,QQ(LLO1),IQQ(LLO2))
*
      ELSEIF (CHPATL.EQ.'PPART') THEN
         CALL KUGETI(NUMB)
         IF(IWKSTY.GE.1.AND.IWKSTY.LE.10) THEN
            CALL GPPARX(NUMB)
         ELSE
            CALL GPPART(NUMB)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'PRKINE') THEN
         CALL KUGETI(NUMB)
         IF(IWKSTY.GE.1.AND.IWKSTY.LE.10) THEN
            CALL GPKINX(NUMB)
         ELSE
            CALL GPKINE(NUMB)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'STPAR') THEN
         CALL KUGETI(ITMED)
         CALL KUGETC(CHPAR,NCH)
         CALL KUGETR(PARVAL)
         CALL GSTPAR(ITMED,CHPAR,PARVAL)
*
      ELSEIF (CHPATL.EQ.'RNDM') THEN
        IF (NPAR.EQ.0) THEN
          CALL GRNDMQ(IS1,IS2,0,'G')
          WRITE(6, 10000) IS1, IS2
        ELSEIF (NPAR.EQ.2) THEN
          CALL KUGETI(NRNDM(1))
          CALL KUGETI(NRNDM(2))
          IF ( NRNDM(2) .NE. 0 ) THEN
            CALL GRNDMQ(NRNDM(1),NRNDM(2),1,'S')
          ELSE
            ISEQ=NRNDM(1)
            CALL GRNDMQ(NRNDM(1),NRNDM(2),ISEQ,'Q')
            CALL GRNDMQ(NRNDM(1),NRNDM(2),ISEQ,'S')
          ENDIF
        ELSE
          WRITE(6, 10100)

10100 FORMAT(' Please give 2 seeds to initialize the random ',
     +       'number generator')

10000 FORMAT('  (E154) RANDOM NUMBER = ',2I12)
         ENDIF
      ENDIF
*
   30 END
