*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
*
      SUBROUTINE GXOBJ(BRCLAS,BRNAME,BRPATH,OBNAME,OBCLAS,STEXT,LTEXT)
*
***********************************************************************
*                                                                     *
*      This routine scans the Geant data structures to return each    *
*       time the next object found. It is used in the new user in-    *
*       terface based on KUIP-MOTIF.                                  *
*                                                                     *
*      Author :    S. Giani     ********************                  *
*                                                                     *
***********************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcdraw.inc"
*
      CHARACTER*(*) BRCLAS,BRNAME,BRPATH,OBNAME,OBCLAS,STEXT,LTEXT
      CHARACTER*4 NAMV,NAMS,IUSET,IUDET
      CHARACTER*20 NATMED,NAPART,NAMM
      DIMENSION UBUF(1),VBUF(10),WBUF(10)
      SAVE IDPOIN,ISET
*
      IF(OBNAME.EQ.' ')IDPOIN=0
*
      IF(BRCLAS.EQ.'VOLU')THEN
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NVOLUM+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            CALL GFVOLU(IDPOIN,NAMV,NAMS)
            OBNAME=NAMV
            OBCLAS=NAMS
            WRITE(STEXT,'(I12)') IDPOIN
            LTEXT=OBCLAS
         ENDIF
      ELSEIF(BRCLAS.EQ.'MATE')THEN
   10    CONTINUE
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NMATE+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            CALL GFMATE(IDPOIN,NAMM,A,Z,DENS,RADL,ABSL,UBUF,NWBUF)
            JMA = LQ(JMATE-IDPOIN)
            IFF = Q(JMA+11)
            IF(A.GE.0.)THEN
               OBNAME=NAMM
               IF(IFF.LE.1)THEN
                  OBCLAS='Elem'
               ELSE
                  OBCLAS='Mixt'
               ENDIF
               WRITE(STEXT,'(I12)') IDPOIN
               LTEXT=OBCLAS
            ELSE
               GOTO 10
            ENDIF
         ENDIF
      ELSEIF(BRCLAS.EQ.'TMED')THEN
   20    CONTINUE
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NTMED+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            CALL GFTMED(IDPOIN,NATMED,NMAT,ISVOL,IFIELD,FIELDM, TMAXFD,
     +      STEMAX,DEEMAX,EPSIL,STMIN,UBUF,NWBUF)
            IF(NMAT.NE.0)THEN
               OBNAME=NATMED
               OBCLAS='Med'
               WRITE(STEXT,'(I12)') IDPOIN
               LTEXT=OBCLAS
            ELSE
               GOTO 20
            ENDIF
         ENDIF
      ELSEIF(BRCLAS.EQ.'PART')THEN
   30    CONTINUE
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NPART+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            CALL GFPART(IDPOIN,NAPART,ITRTYP,AMASS,CHARGE,TLIFE,
     +      VBUF,NWBUF)
            IF(ITRTYP.NE.0)THEN
               OBNAME=NAPART
               OBCLAS='Part'
               WRITE(STEXT,'(I12)') IDPOIN
               LTEXT=OBCLAS
            ELSE
               GOTO 30
            ENDIF
         ENDIF
      ELSEIF(BRCLAS.EQ.'KINE')THEN
   40    CONTINUE
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NTRACK+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            CALL GFKINE(IDPOIN,VERT,PVERT,IPART,NVERT,WBUF,NWBUF)
            IF(IPART.NE.0)THEN
               WRITE(OBNAME,'(I12)') IDPOIN
               OBCLAS='Kine'
               WRITE(STEXT,'(I12)') IDPOIN
               LTEXT=OBCLAS
            ELSE
               GOTO 40
            ENDIF
         ENDIF
      ELSEIF(BRCLAS.EQ.'HITS')THEN
        IF(BRPATH.EQ.' ')THEN
         IF(JSET.NE.0)THEN
           NSET=IQ(JSET-1)
         ELSE
           NSET=0
         ENDIF
   45    CONTINUE
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NSET+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            CALL UHTOC(IQ(JSET+IDPOIN),4,IUSET,4)
            CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
            IF(ISET.NE.0)THEN
               OBNAME=IUSET
               OBCLAS='Hitset'
               WRITE(STEXT,'(I12)') IDPOIN
               LTEXT=OBCLAS
            ELSE
               GOTO 45
            ENDIF
         ENDIF
        ELSE
         IF(JSET.NE.0)THEN
           NSET=IQ(JSET-1)
           IUSET=BRPATH(2:5)
           CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
           NDET=IQ(LQ(JSET-ISET)-1)
         ELSE
           NDET=0
         ENDIF
   46    CONTINUE
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NDET+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            CALL UHTOC(IQ(LQ(JSET-ISET)+IDPOIN),4,IUDET,4)
            CALL GLOOK(IUDET,IQ(LQ(JSET-ISET)+1),NDET,IDET)
            IF(IDET.NE.0)THEN
               OBNAME=IUDET
               OBCLAS='Hitdet'
               WRITE(STEXT,'(I12)') IDPOIN
               LTEXT=OBCLAS
            ELSE
               GOTO 46
            ENDIF
         ENDIF
        ENDIF
      ELSEIF(BRCLAS.EQ.'ROTM')THEN
   50    CONTINUE
         IDPOIN=IDPOIN+1
         IF((IDPOIN.EQ.(IQ(JROTM-2)+1)).OR.JROTM.LE.0)THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            JR=LQ(JROTM-IDPOIN)
            IF(JR.GT.0)THEN
               WRITE(OBNAME,'(I12)') IDPOIN
               OBCLAS='Rmatr'
               WRITE(STEXT,'(I12)') IDPOIN
               LTEXT=OBCLAS
            ELSE
               GOTO 50
            ENDIF
         ENDIF
      ELSEIF(BRCLAS.EQ.'VIEW')THEN
         IDPOIN=IDPOIN+1
         IF(IDPOIN.EQ.(NKVIEW+1))THEN
            IDPOIN=0
            OBNAME=' '
         ELSE
            WRITE(OBNAME,'(I12)') IDPOIN
            OBCLAS='VB'
            WRITE(STEXT,'(I12)') IDPOIN
            LTEXT=OBCLAS
         ENDIF
      ELSE
         PRINT*,'BRCLAS NOT VALID'
      ENDIF
*
      END
