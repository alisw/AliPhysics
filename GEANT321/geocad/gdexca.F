*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 13/09/94  18.53.20  by  S.Ravndal
*-- Author :
      SUBROUTINE GDEXCA(NAME,NBINS)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Based on GDRAW, calculates parameters of each volume        *
C.    *    Areas marked JV + NH                                        *
C.    *                                                                *
C.    *    Called by GTXSET                                            *
C.    *                                                                *
C.    *       Authors : R.Brun, A.McPherson, P.Zanarini,   *********   *
C.    *                 J.Salt, S.Giani , J. Vuoskoski,  N. Hoimyr     *
C.    ******************************************************************
C.
C JV
#include "geant321/gcsetf.inc"
C
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gconst.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcdlin.inc"
#include "geant321/gcmutr.inc"
*JS
#include "geant321/gcgobj.inc"
C+SEQ,CGHPAR.
#include "geant321/gchiln.inc"
#include "geant321/gcspee.inc"
*JS
*
*
C   this by jv
      DIMENSION PARMJV(9), POSJV(3)
C         if volume is divided jdvinf(level) is 1
      DIMENSION JDVINF(0:15)
      CHARACTER*4 JVVOLU,JVVOLD, NAME,NAMEE2
      CHARACTER*10 VOLNAM
      DIMENSION X(3),ATT(10)
      DIMENSION LVOLS(15),LINDX(15),LNAMES(15)
      DIMENSION GPAR(50,15)
*     INTEGER START, OLDOLD, PASS
C
      IF(JCADNT.EQ.0) THEN
         CALL MZBOOK(IXSTOR,JCADNT,JCADNT,1,'CADI',1,1,0,2,-1)
         CALL MZBOOK(IXSTOR,JBUF1,
     +               JCADNT,-1,'CAD1',0,0,NVOLUM,2,-1)
      ENDIF
      DO 10 JV=1,NVOLUM
         IQ(JBUF1+JV)=0
   10 CONTINUE
      JDVINF(0)=0
      JLEVEL=0
      MYSEEN=1
      JVVOLU='----'
C
C            Set IOBJ to VOLUME
C
      IOBJ=1
C
C             Save /GCVOLU/ if necessary
C
      IFCVOL=0
      IF (NLEVEL.NE.0) THEN
         CALL GSCVOL
         IFCVOL=1
      ENDIF
      IF (NLEVEL.LT.0) NLEVEL=IABS(NLEVEL)
C
C             Start of general code
C
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF(IVO.LE.0)GO TO 210
C
C             Theta, phi and psi angles are normalized in [0-360[ range
C
*
      JVO=LQ(JVOLUM-IVO)
C
C             Initialize JIN to switch correctly CALL GFPARA/GFIPAR
C
      JIN=0
C
      LEVSEE=1000
C
      IF (IDRNUM.NE.0) GO TO 30
C
C             Initialize for new geometry structure
C
      IF (JGPAR.EQ.0) CALL GMEDIN
      CALL GLMOTH(NAME,1,NLEV,LVOLS,LINDX)
      DO 20 J=1, NLEV
         LNAMES(J)=IQ(JVOLUM+LVOLS(J))
   20 CONTINUE
      NLEV=NLEV+1
      CALL UCTOH(NAME,LNAMES(NLEV),4,4)
      LINDX(NLEV)=1
      CALL GLVOLU(NLEV, LNAMES, LINDX, IER)
C
      NLVTOP=NLEVEL
C
   30 CONTINUE
C
      NLMIN=NLEVEL
      NLMAX=NLEVEL
C
      IF (IDRNUM.NE.0) GO TO 70
C
      CALL GFPARA(NAME,1,1,NPAR,NATT,GPAR(1,NLEVEL),ATT)
C
      IF (NPAR.LE.0) GO TO 220
C
      DO 60  LLL=1,NLEVEL
         DO 50 I=1,3
            GTRAN(I,LLL)=0.0
            X(I)=0.0
            DO 40 J=1,3
               K=(I-1)*3+J
               GRMAT(K,LLL)=0.0
   40       CONTINUE
            K=I*4-3
            GRMAT(K,LLL)=1.0
   50    CONTINUE
         GRMAT(10,LLL)=0.0
   60 CONTINUE
C
C             Ready for general case code
C
   70 CONTINUE
*SG
*    Taking volume name and shape from Zebra Structure
*
      IVOLNA=IQ(JVOLUM+IVO)
      ISHAPE=Q(JVO+2)
*SG
C
      IF (IDRNUM.NE.0) GO TO 80
C
      IF (NLEVEL.EQ.NLVTOP) GO TO 90
C
   80 CONTINUE
C
      IF (IDRNUM.NE.0.AND.JIN.EQ.0) THEN
         CALL UHTOC(NAMES(NLEVEL),4,NAMEE2,4)
         CALL GFPARA(NAMEE2,NUMBER(NLEVEL),1,NPAR,
     +            NATT,GPAR(1,NLEVEL),ATT)
      ELSE
         NPAR=Q(JVO+5)
         NATT=Q(JVO+6)
         JATT=JVO+7+NPAR
         CALL UCOPY(Q(JATT),ATT,NATT)
      ENDIF
C
   90 CONTINUE
C
      WORK=ATT(1)
      SEEN=ATT(2)
 
*
      LINSTY=ATT(3)
      LINWID=ATT(4)
      LINCOL=ATT(5)
      LINFIL=ATT(6)
      CALL MVBITS(LINCOL,0,4,LINATT,7)
      CALL MVBITS(LINWID,0,3,LINATT,11)
      CALL MVBITS(LINSTY,0,3,LINATT,14)
      CALL MVBITS(LINFIL,0,3,LINATT,17)
*
***SG
*
*    New logic scanning the geometrical tree:
*    A volume can set bounds OR be compared with bounds;
*    this can happen only IF a relationship mother-daughters exists.
*
*    Optimization for Hidden Volume Removal:
*    POS and DIV cases are considered at the same time.
*
*     IF(START.EQ.1)THEN
*
*        IF(NLEVEL.GT.OLDOLD)THEN
*           IF(PASS.NE.0)THEN
*              S1=0
*              S2=0
*              S3=0
*              SS1=0
*              SS2=0
*              SS3=0
*              SRAGMX=0
*              SRAGMN=0
*              PASS=0
*              IF(SEEN.EQ.0.OR.SEEN.EQ.-1)PASS=1
*              OLDOLD=NLEVEL
*           ENDIF
*
*        ELSE IF(NLEVEL.LE.OLDOLD)THEN
*           S1=0
*           S2=0
*           S3=0
*           SS1=0
*           SS2=0
*           SS3=0
*           SRAGMX=0
*           SRAGMN=0
*           PASS=0
*           IF(SEEN.EQ.0.OR.SEEN.EQ.-1)PASS=1
*           OLDOLD=NLEVEL
*        ENDIF
*     ENDIF
*
*     IF(OLDOLD.EQ.0.AND.(SEEN.EQ.1.OR.SEEN.EQ.-2))THEN
*        START=1
*        PASS=0
*        OLDOLD=NLEVEL
*     ENDIF
C
C             WORK attribute enabled ?
C
      IF(WORK.LE.0.)GO TO 200
C
C             SEEN attribute processing
C
      IF (SEEN.LT.50.) GO TO 100
      ISEENL=SEEN/10.+0.5
      SEEN=ISEENL-10
  100 CONTINUE
C
      IF(NLEVEL.LE.LEVSEE)LEVSEE=1000
      IF(SEEN.EQ.-1.)GO TO 200
      IF (NLEVEL.GT.LEVSEE) GO TO 200
      IF(SEEN.EQ.0.)GO TO 150
      IF (SEEN.EQ.-2.) LEVSEE=NLEVEL
C
*              Standard Mode: Output to SET
*
C-----------------------JV----Mod NH---------------------------
C        Get positioning variables
C
      IF(NLEVEL.LT.JLEVEL)THEN
         CALL GSATT(JVVOLU,'SEEN',MYSEEN)
         Q(JBUF1+MYIVO)=0
      ENDIF
      JVO=LQ(JVOLUM-IVO)
      NIN=Q(JVO+3)
      IF(NIN.LT.0) THEN
         JDVINF(NLEVEL)=1
      ELSE
         JDVINF(NLEVEL)=0
      ENDIF
C
C            if division
      IF (JDVINF(NLEVEL-1).EQ.1) THEN
C
         IF (IQ(JBUF1+IVO).LT.NBINS) THEN
            DO 110 JV=1, 9
               PARMJV(JV)=GRMAT(JV,NLEVEL)
  110       CONTINUE
C
            DO 120 JV=1, 3
               POSJV(JV)=GTRAN(JV,NLEVEL)
  120       CONTINUE
C
C        Appends new name VOLNAM to each volume, with index.
C
            IQ(JBUF1+IVO)=IQ(JBUF1+IVO)+1
            WRITE(VOLNAM(1:5),10200)IVOLNA
            WRITE(VOLNAM(6:10),'(I4.0)')IQ(JBUF1+IVO)
C
C        Call SHAPE to SET routines
C
C        Updates SET block sequence number:
            N1=N1+1
            CALL GETSHP(ISHAPE,GPAR(1,NLEVEL))
C
C        Position the volumes
C
            N1=N1+1
            CALL GPOSI(PARMJV,POSJV,VOLNAM,LINCOL)
C
         ELSE
            JVVOLD=JVVOLU
            CALL UHTOC(IVOLNA,4,JVVOLU,4)
            IF(JVVOLD.NE.JVVOLU)MYSEEN=ATT(2)
            CALL GSATT(JVVOLU,'SEEN',-1)
            JLEVEL=NLEVEL
            MYIVO=IVO
            GOTO 200
         ENDIF
C
C       normal volumes
      ELSE
         DO 130 JV=1, 9
            PARMJV(JV)=GRMAT(JV,NLEVEL)
  130    CONTINUE
C
         DO 140 JV=1, 3
            POSJV(JV)=GTRAN(JV,NLEVEL)
  140    CONTINUE
C
C        Appends new name VOLNAM to each volume, with index.
C
         IQ(JBUF1+IVO)=IQ(JBUF1+IVO)+1
         WRITE(VOLNAM(1:5),10200)IVOLNA
         WRITE(VOLNAM(6:10),'(I4.0)')IQ(JBUF1+IVO)
C
C        Call SHAPE to SET routines
C
C        Updates SET block sequence number:
         N1=N1+1
         CALL GETSHP(ISHAPE,GPAR(1,NLEVEL))
C
C        Position the volumes
C
         N1=N1+1
         CALL GPOSI(PARMJV,POSJV,VOLNAM,LINCOL)
      ENDIF
C------------------------------------------------------------------------
C         Output of material list
C
      IF (IQ(JBUF1+IVO).EQ.1) THEN
         NTRMED=Q(JVO+4)
         CALL GPTSET (IVOLNA, NTRMED)
      ENDIF
C------------------------------END JV------------------------------------
C
*JS
*
*     Logic has been modified  >>>>>
*
*
*JS
      IF(SEEN.EQ.-2.)GO TO 200
C
  150 CONTINUE
C
***   IF (IDRNUM.NE.0) GO TO 999
C
C             Skip User shapes (not yet implemented)
C
**      ISEARC=Q(JVO+1)
C
C             Now go down the tree
C
      NIN=Q(JVO+3)
      IF(NIN.EQ.0) GO TO 200
      IF(NIN.LT.0) GO TO 170
C
C             Contents placed by GSPOS
C
      IN=0
      IF(NLMAX.GT.NLEVEL) IN=LINDEX(NLEVEL+1)
      IN=IN+1
      IF(IN.GT.NIN.AND.NLEVEL.EQ.NLMIN) GO TO 230
*
      IF(IN.GT.NIN) GO TO 190
*
      CALL GMEPOS(JVO,IN,X,0)
*
      NPAR=IQ(JGPAR+NLEVEL)
      DO 160 I=1,NPAR
         GPAR(I,NLEVEL)=Q(LQ(JGPAR-NLEVEL)+I)
  160 CONTINUE
*
      IVO=LVOLUM(NLEVEL)
      JVO=LQ(JVOLUM-IVO)
      NLMAX=NLEVEL
      GO TO 70
C
  170 CONTINUE
C
C             Contents by division
C
      IN=0
      IF(NLMAX.GT.NLEVEL) IN=LINDEX(NLEVEL+1)
      IN=IN+1
      CALL GMEDIV(JVO,IN,X,0)
*
      IF (IN.EQ.0) GO TO 190
*
      NPAR=IQ(JGPAR+NLEVEL)
      DO 180 I=1,NPAR
         GPAR(I,NLEVEL)=Q(LQ(JGPAR-NLEVEL)+I)
  180 CONTINUE
*
      IF (IN.EQ.0) GO TO 190
*
      IVO=LVOLUM(NLEVEL)
      JVO=LQ(JVOLUM-IVO)
      NLMAX=NLEVEL
      GO TO 70
C
  190 CONTINUE
      NLMAX=NLEVEL
  200 CONTINUE
      NLEVEL=NLEVEL-1
      IF(NLEVEL.LT.NLMIN) GO TO 230
      IVO=LVOLUM(NLEVEL)
      JVO=LQ(JVOLUM-IVO)
      GO TO 150
C
  210 WRITE(CHMAIL,10000)NAME
      CALL GMAIL(0,0)
      GO TO 230
C
  220 CONTINUE
C
C             TOP OF THE TREE HAS PARAMETERS SET BY GSPOSP.
C             BUT GDRAW DOES NOT HAVE ACCESS TO THE IN BANK
C             WHICH PLACED IT IN ITS MOTHER.
C
      WRITE(CHMAIL,10100) NAME
      CALL GMAIL(0,0)
C
  230 CONTINUE
*
***SG
*
*JS
      IF(KCGST.EQ.-9)THEN
         KSHIFT=0
         IF(JCG.NE.0)CALL MZDROP(IXSTOR,JCG,' ')
         IF(JCGOBJ.NE.0)CALL MZDROP(IXSTOR,JCGOBJ,' ')
         CALL GDCGRS
         IF(JCGCOL.NE.0)CALL MZDROP(IXSTOR,JCGCOL,' ')
         LARECG(1)=0
         CALL MZGARB(IXSTOR+1,0)
         NCLAS2=0
         NCLAS3=0
      ENDIF
      ICUT=0
      IF (IFCVOL.EQ.1) THEN
         CALL GFCVOL
      ELSE
         NLEVEL=0
      ENDIF
C
C             Restore permanent value of color and return
C
      CALL GDCOL(0)
      IOBJ=0
C
10000 FORMAT(' *** GDEXCA *** : Volume ',A4,' does not exist')
10100 FORMAT(' *** GDEXCA *** : Top of tree ',A4,' parameters defined',
     +       '  by GSPOSP - info not available to GDEXCA.')
10200 FORMAT(A4,'_')
      END
