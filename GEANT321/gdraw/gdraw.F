*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDRAW(NAME,UTHET,UPHI,UPSI,UU0,UV0,SU,SV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine draws the object called NAME, with its      *
C.    *       contents, at the screen point (UU0,UV0), with the        *
C.    *       screen factors SU and SV acting on the U and V           *
C.    *       dimensions respectively;                                 *
C.    *       the object is rotated by an angle UTHET along Y-axis     *
C.    *       and UPHI along Z-axis and the resulting 2-D picture      *
C.    *       is also rotated by an angle UPSI along the line of       *
C.    *       projection (i.e. the normal to the 2-D view plane).      *
C.    *                                                                *
C.    *       If IDRNUM<>0 then /GCVOLU/ is already filled by GLVOLU   *
C.    *       and a special case is handled (GDRAW called by GDRVOL).  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDPRTR, GDRAWC, GDRAWX,     *
C.    *                   GDSPEC                                       *
C.    *       Authors : R.Brun, A.McPherson, P.Zanarini,   *********   *
C.    *                 J.Salt, S.Giani                                *
C.    ******************************************************************
C.
***SG**************************************************************************
*                                                                             *
*      Most important improvements in this new version :                      *
*                                                                             *
*      Problems with the number of faces are solved, so that it's             *
*         possible to use all the memory available; moreover it's             *
*         now possible to use HIDE ON on CRAY-like machines;                  *
*                                                                             *
*      Problems with number of volumes now only depends upon the              *
*         size of Zebra store : a message will tell the number of             *
*         words you need more; problems in iterated drawings have             *
*         been solved too.                                                    *
*                                                                             *
*      A NEW SIZE EVALUATION is performed separately for Hide Structure and   *
*         Wire Structure with a resolution of 1 word ; this is the new logic: *
*         create immediately HIDE and WIRE structures and perform the         *
*         drawing WHILE evaluating the memory used; if memory booked in the   *
*         zebra store is not enough, then go on evaluating the number of      *
*         words needed and print it.                                          *
*                                                                             *
*      Multi-colour view of the different parts of a detector is now          *
*         available in Hidden Line Removal; a new bank is created for this.   *
*         Enjoy clipping now !!                                               *
*         Different line styles and width work as well !!                     *
*                                                                             *
*      Zooming is now possible in Hidden Line Removal too; Dspec works        *
*         even when Cvol is on and Seen attribute setting has been            *
*         optimized.                                                          *
*                                                                             *
*      Speed in drawing divided volumes can be increased by a factor about    *
*         linear with the number of volumes (a factor 30 for 900 tubes)!!!    *
*         At the same time, the number of words used can be decreased by a    *
*         factor 50!!! For example, in Gexam1 is possible to draw 22500 tubes *
*         using much less than 800000 words.                                  *
*                                                                             *
*      HIDDEN FACE REMOVAL algorithm has been implemented; it allows to       *
*         increase speed and decrease memory used by the same factors as      *
*         above for any kind of drawings!!! For example, it's now possible    *
*         to draw the complete L3 geometry using less than 3 Mwords (before   *
*         we needed 1 Gigabyte !!!)                                           *
*                                                                             *
*      The new command CVOL has been created: it allows to clip EACH          *
*         VOLUME in the detector by a different SHAPE; moreover, it's possible*
*         to clip twice each volume. You can clip each volume by the following*
*         shapes: BOX , TUBE, CONE, SPHE !                                    *
*                                                                             *
*      The new command SHIFT has been created: it allows to translate each    *
*         volume in the detector into a more visible place; for each volume,  *
*         the last shift you asked is performed.                              *
*                                                                             *
*      The option 'one' has been implemented in the shift command to allow    *
*         the shifting of a single copy for each volume; the new command      *
*         BOMB has been created to allow 'exploded' view of detectors.        *
*                                                                             *
*      A new bank is booked to create CG objects: it's now possible to draw   *
*         in HIDE ON the following shapes as well: PCON, TUBS, SPHE, TRAP     *
*         and Pyramids as TRAP having 4 edges in the same point. Moreover,    *
*         the new shape CTUB has been created even in hide on.                *
*                                                                             *
*      The number of words to draw TUBS, CONS, PCON, PGON can be decreased    *
*         even by a factor 10; moreover, several problems about edge visi-    *
*         bility have been solved.                                            *
*                                                                             *
*      A REFLECTION algorithm has been implemented for hidden line removal    *
*         as well; finally, a new logic scanning the geometrical tree has     *
*         been created to simulate, without alterate, the date structure.     *
*                                                                             *
*      A new SURFACE SHADING algorithm has been written to fill faces with    *
*         solid colours with varying intensity according to the light         *
*         inclination. Please see details in the documentation of the         *
*         routines which are concerned.                                       *
*                                                                             *
***SG**************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gconst.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcmutr.inc"
#include "geant321/gcgobj.inc"
#if defined(CERNLIB_CG)
#include "geant321/cghpar.inc"
#endif
#include "geant321/gchiln.inc"
#include "geant321/gcspee.inc"
*
*    - The following common to be used by GXPICK
*
      COMMON /QUEST/ IQUEST(100)
      COMMON/GCVHLP/NVLAST
      COMMON/SP3D/ISPFLA
      COMMON/INIFIR/NFIRST
*
      CHARACTER*4 NAME,NAMEE2
      DIMENSION X(3),ATT(10)
      DIMENSION LVOLS(15),LINDX(15),LNAMES(15)
      DIMENSION GPAR(50,15)
#if defined(CERNLIB_CG)
      DIMENSION V(3),T(4,3)
#endif
      SAVE IFIRST, LFIRST
      DATA IFIRST,LFIRST/2*0/
*       Save info for GXPICK in case is needed
      CALL UCTOH(NAME, NVLAST, 4, 4)
***      CALL IGSET('SYNC',1.)
*
*       Hidden flag 'ON" (Default)
*
      CALL UCTOH('ON  ',IFLH,4,4)
*
***SG
*
      IF(NFIRST.EQ.0)THEN
        CALL GDCOTA
        NFIRST=1
      ENDIF
#if defined(CERNLIB_CG)
      IF(IHIDEN.EQ.IFLH)THEN
         IF(LFIRST.EQ.0) THEN
*      Link area for the banks
            CALL MZLINT(IXSTOR,'/GCHILN/',LARECG,ICLIP1,ICLIP2)
            LFIRST = 1
         ELSE
            IF(JCGCOL.NE.0) CALL MZDROP(IXSTOR, JCGCOL, ' ')
            IF(JCGOBJ.NE.0) CALL MZDROP(IXSTOR, JCGOBJ, ' ')
            IF(JCOUNT.NE.0) CALL MZDROP(IXSTOR, JCOUNT, ' ')
            IF(JCLIPS.NE.0) CALL MZDROP(IXSTOR, JCLIPS, ' ')
            IF(IMPOIN.NE.0) CALL MZDROP(IXSTOR, IMPOIN, ' ')
            IF(IMCOUN.NE.0) CALL MZDROP(IXSTOR, IMCOUN, ' ')
            IF(JSIX.NE.0) CALL MZDROP(IXSTOR, JSIX, ' ')
            IF(JSIY.NE.0) CALL MZDROP(IXSTOR, JSIY, ' ')
            IF(JSIZ.NE.0) CALL MZDROP(IXSTOR, JSIZ, ' ')
            IF(JPXC.NE.0) CALL MZDROP(IXSTOR, JPXC, ' ')
            IF(JPYC.NE.0) CALL MZDROP(IXSTOR, JPYC, ' ')
            IF(JPZC.NE.0) CALL MZDROP(IXSTOR, JPZC, ' ')
         ENDIF
         LARECG(1)=1
*
*                     Initialization
*
*  NWHS1: n. of words for Hide Structure
*  NWWS1: n. of words for Wire Structure
*  NWFLAG: Indicates if the size of CG bank is precise
*          =0 , it is
*          =-9, it isn't
*  IPAS: This flag indicates if the Size Evaluation has been performed
*        =0 , it does not
*        =1 , it does
*  NOBJ:  Counter of CG objects
*  NUVO:  Counter of CG volumes
*  II: Counter for volumes' line attributes
*  KGG: Flag for booking line attributes bank
*  LSTEP: Number of CG objects forming each volume
*  IFACST: Flag indicating final status of Hide Structure
*          =0 , it's ok
*          <0 , internal error
*          >0 , total number of words needed for Hide Structure
*  NCLAS2: Total number of volumes
*  S1...SS3: Min and Max values of volume scope
*  SRAGMX,SRAGMN: Max values of volume scope along R
*  NFILT= n. of words for HIDE+totalWIRE structures+CG+Line
*  NTCUR= n. of words for HIDE+instantWIRE structures+CG+Line
*  KSTART: Flag for Hidden Volume Removal
*  IOLDOL: Nlevel of last volume setting bounds for scope
*
*         NWHS1=0
*         NWWS1=0
*         NWFLAG=0
         IPER=0
         IPEOLD=0
         NOBJ=0
         NUVO=0
         IPAS=0
         II=0
         KGG=0
         LSTEP=1
         IFACST=0
*         IWORK=0
         NCLAS1=0
         NCLAS2=0
         NCLAS3=0
         IIIIII=0
*    Initialization of Hidden Volume and Face Removal
         S1=0
         S2=0
         S3=0
         SS1=0
         SS2=0
         SS3=0
         SRAGMX=0
         SRAGMN=0
         RAINT1=0
         RAINT2=0
         ISCOP=0
         KSTART=0
         IOLDOL=0
*    Initialisation for Shift
         NIET=0
         IOLDSU=0
         PORGX=0
         PORGY=0
         PORGZ=0
         DO 10  J=1,15
            POX(J)=0
            POY(J)=0
            POZ(J)=0
   10    CONTINUE
         DO 20  J=1,100
            IVECVO(J)=0
   20    CONTINUE
         IVOOLD=0
*         IMENO=0
         IPRELE=0
*    Resetting
         IHPOIN=0
         IWPOIN=0
         ICLIP1=0
         ICLIP2=0
         IVOLNA=0
         LPASS=0
         NWHS=0
         MFLA=0
         MVENLE=0
         MVECOL=0
         LFLA=0
         LVENLE=0
         LVEWID=0
         LFFLA=0
         LFENLE=0
         LVEFIL=0
         IXCG=IXSTOR+1
         JMEMT1=0
         JMEMT3=0
         JMEMT2=0
*    Resetting
         IF(JCG.NE.0)THEN
            CALL MZDROP(IXSTOR, JCG, ' ')
            CALL MZGARB(IXSTOR+1,0)
         ENDIF
*    Booking bank to create CG objects
         CALL MZNEED(IXCG,30000,'G')
         CALL MZBOOK(IXCG,JCGOBJ,JCGOBJ,1,'CGOB',0,0,30000,3,-1)
         CALL MZNEED(IXCG,33000,'G')
         CALL MZBOOK(IXCG,JCLIPS,JCLIPS,1,'CGCLIP',0,0,33000,3,-1)
         ICLIP1=JCLIPS+1
         ICLIP2=JCLIPS+16500
         JMEMT1=IQUEST(11)*.013
         IF(JMEMT1.LT.10000)JMEMT1=10000
         CALL MZNEED(IXCG,JMEMT1,'G')
         CALL MZBOOK(IXCG,JCOUNT,JCOUNT,1,'CGCONT',0,0,JMEMT1,2,-1)
         IQ(JCOUNT+1)=1
         IQ(JCOUNT+2)=4000
         IQ(JCOUNT+3)=8000
         LLEP=ABS(LEP)
         IF(LLEP.NE.1)THEN
            CALL MZNEED(IXCG,54000,'G')
            CALL MZBOOK(IXCG,JSIX,JSIX,1,'XGEN',0,0,9000,3,-1)
            CALL MZBOOK(IXCG,JSIY,JSIY,1,'YGEN',0,0,9000,3,-1)
            CALL MZBOOK(IXCG,JSIZ,JSIZ,1,'ZGEN',0,0,9000,3,-1)
            CALL MZBOOK(IXCG,JPXC,JPXC,1,'XPAR',0,0,9000,3,-1)
            CALL MZBOOK(IXCG,JPYC,JPYC,1,'YPAR',0,0,9000,3,-1)
            CALL MZBOOK(IXCG,JPZC,JPZC,1,'ZPAR',0,0,9000,3,-1)
            JMEMT3=JMEMT1
            CALL MZNEED(IXCG,JMEMT3,'G')
            CALL MZBOOK(IXCG,IMCOUN,IMCOUN,1,'SHCONT',0,0,JMEMT3,2,-1)
            IQ(IMCOUN+1)=1
            IQ(IMCOUN+2)=4000
            IQ(IMCOUN+3)=8000
            JMEMT2=IQUEST(11)*.1
            CALL MZNEED(IXCG,JMEMT2,'G')
            CALL MZBOOK(IXCG,IMPOIN,IMPOIN,1,'SHAFAC',0,0,JMEMT2,2,-1)
         ENDIF
*    Resetting
         CALL GDCGRS
         NFILT=0
         NTCUR=0
         NTNEX=0
         NAIN=0
         ITSTCU=0
      ENDIF
#endif
*
***SG
C
C            Set IOBJ to VOLUME
C
      IOBJ=1
C
      IF (IFIRST.NE.0) GO TO 40
C
      IFIRST=1
      DPHI=PI/20.
      PHI=0.
C
      DO 30 I=1,40
         GSIN(I)=SIN(PHI)
         GCOS(I)=COS(PHI)
         PHI=PHI+DPHI
   30 CONTINUE
C
      GSIN(41)=GSIN(1)
      GCOS(41)=GCOS(1)
C
   40 CONTINUE
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
C             If in cut-mode then open the GDRAWV line buffer
C             else reset ICUT that could have been set by GDRAWC/GDRAWX
C
      IF (ICUTFL.EQ.1) THEN
         CALL GDRAWV(0.,0.,-1)
      ELSE
         ICUT=0
      ENDIF
C
C             Start of general code
C
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF(IVO.LE.0)GO TO 280
C
C             Theta, phi and psi angles are normalized in [0-360[ range
C
 
      GTHETA=MOD(ABS(UTHET),360.)
      GPHI=MOD(ABS(UPHI),360.)
      GPSI=MOD(ABS(UPSI),360.)
      GU0=UU0
      GV0=UV0
      GSCU=SU
      GSCV=SV
      IMOD=0
*
#if defined(CERNLIB_CG)
*              Set Transformation Matrix T for CG Package
*
      IF(IHIDEN.EQ.IFLH)THEN
         V(1)=GTHETA
         V(2)=GPHI
         V(3)=GPSI
         CALL GDCGVW(V,T)
         CALL CGTSET(NTRCG,T,IREP)
         IF(IREP.EQ.-1)THEN
            WRITE(CHMAIL,10200)
            CALL GMAIL(0,0)
         ENDIF
         IF(IREP.EQ.-2)THEN
            WRITE(CHMAIL,10300)
            CALL GMAIL(0,0)
         ENDIF
*
*    Obtaining the IMOD flag for setting the run mode
*
         IF(IPAS.EQ.0)THEN
            IMOD=0
            IF(ICUT.NE.0)IMOD=1
            IF(IHOLE.EQ.1)IMOD=2
         ELSE
            IMOD=3
            IF(ICUT.NE.0)IMOD=4
            IF(IHOLE.EQ.1)IMOD=5
            IF(ICUT.NE.0.OR.IHOLE.EQ.1)THEN
               IFCG=4
               ILCG=3
            ENDIF
         ENDIF
      ENDIF
#endif
*JS
*  77  CONTINUE
      SINPSI=SIN(GPSI*DEGRAD)
      COSPSI=COS(GPSI*DEGRAD)
      GU0=UU0
      GV0=UV0
      NGVIEW=0
      JVO=LQ(JVOLUM-IVO)
C
C             Initialize JIN to switch correctly CALL GFPARA/GFIPAR
C
      JIN=0
C
      LEVSEE=1000
C
      IF (IDRNUM.NE.0) GO TO 70
C
C             Initialize for new geometry structure
C
      IF (JGPAR.EQ.0) CALL GMEDIN
      CALL GLMOTH(NAME,1,NLEV,LVOLS,LINDX)
      DO 50 J=1, NLEV
         LNAMES(J)=IQ(JVOLUM+LVOLS(J))
   50 CONTINUE
      NLEV=NLEV+1
      CALL UCTOH(NAME,LNAMES(NLEV),4,4)
      LINDX(NLEV)=1
      DO 60   KLEV=2,NLEV
         JVOF = LQ(JVOLUM-LVOLS(KLEV-1))
         NIN  = Q(JVOF+3)
         IF(NIN.GT.0) THEN
            JIN = LQ(JVOF-LINDX(KLEV))
            ICOPY = Q(JIN+3)
         ELSE
            ICOPY = 1
         ENDIF
         LINDX(KLEV) = ICOPY
   60 CONTINUE
      CALL GLVOLU(NLEV, LNAMES, LINDX, IER)
C
      NLVTOP=NLEVEL
C
   70 CONTINUE
C
      NLMIN=NLEVEL
      NLMAX=NLEVEL
C
      IF (IDRNUM.NE.0) GO TO 110
C
      CALL GFPARA(NAME,1,1,NPAR,NATT,GPAR(1,NLEVEL),ATT)
C
      IF (NPAR.LE.0) GO TO 290
C
      DO 100 LLL=1,NLEVEL
         DO 90 I=1,3
            GTRAN(I,LLL)=0.0
            X(I)=0.0
            DO 80 J=1,3
               K=(I-1)*3+J
               GRMAT(K,LLL)=0.0
   80       CONTINUE
            K=I*4-3
            GRMAT(K,LLL)=1.0
   90    CONTINUE
         GRMAT(10,LLL)=0.0
  100 CONTINUE
C
C             Ready for general case code
C
  110 CONTINUE
*SG
*    Taking volume name and shape from Zebra Structure
*
      IMENO=IVOLNA
      IVOLNA=IQ(JVOLUM+IVO)
      ISHAPE=Q(JVO+2)
*SG
      GSCU=SU
      GSCV=SV
C
      IF (IDRNUM.NE.0) GO TO 120
C
      IF (NLEVEL.EQ.NLVTOP) GO TO 130
C
  120 CONTINUE
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
  130 CONTINUE
C
      WORK=ATT(1)
      SEEN=ATT(2)
*
      LINSTY=ATT(3)
      LINWID=ATT(4)
      LINCOL=ATT(5)
      LINFIL=ATT(6)
      IF(LLEP.EQ.1)LINFIL=0
*SG
*    New logic setting the line attributes
*
      IF(NLEVEL.EQ.1)THEN
***               CALL GDCOTA
         IF(LINFIL.LT.0)THEN
            LINFIL=ABS(LINFIL)
            CALL ISFACI(LINFIL)
            CALL ISFAIS(1)
            CALL IGBOX(0.,20.,20.,0.)
            CALL ISFAIS(0)
            LINFIL=2
         ENDIF
      ENDIF
      IF(LINCOL.LT.2)THEN
         IF(MFLA.EQ.1.AND.NLEVEL.GT.MVENLE)LINCOL=MVECOL
         IF(NLEVEL.LE.MVENLE)MFLA=0
         IF(LINCOL.LT.0)THEN
            MVECOL=ABS(LINCOL)
            LINCOL=MVECOL
            MVENLE=NLEVEL
            MFLA=1
         ENDIF
      ENDIF
      IF(LINWID.LT.2)THEN
         IF(LFLA.EQ.1.AND.NLEVEL.GT.LVENLE)LINWID=LVEWID
         IF(NLEVEL.LE.LVENLE)LFLA=0
         IF(LINWID.LT.0)THEN
            LVEWID=ABS(LINWID)
            LINWID=LVEWID
            LVENLE=NLEVEL
            LFLA=1
         ENDIF
      ENDIF
      IF(LINFIL.LT.1)THEN
         IF(LFFLA.EQ.1.AND.NLEVEL.GT.LFENLE)LINFIL=LVEFIL
         IF(NLEVEL.LE.LFENLE)LFFLA=0
         IF(LINFIL.LT.0)THEN
            LVEFIL=ABS(LINFIL)
            LINFIL=LVEFIL
            LFENLE=NLEVEL
            LFFLA=1
         ENDIF
      ENDIF
*SG
      CALL MVBITS(LINCOL,0,8,LINATT,16)
      CALL MVBITS(LINWID,0,3,LINATT,7)
      CALL MVBITS(LINSTY,0,3,LINATT,10)
      CALL MVBITS(LINFIL,0,3,LINATT,13)
*
***SG
*
*    New logic scanning the geometrical tree:
*    A volume can set bounds OR be compared with bounds;
*    this can happen only IF a relationship mother-daughters exists.
*
*    Optimization for Hidden Volume and Face Removal:
*    POS and DIV cases are considered at the same time.
*
      IF(IPAS.EQ.0)GOTO 170
      IPORLI=0
      ISUBLI=0
      IF(KSTART.EQ.1)THEN
*
         IF(NLEVEL.GT.IOLDOL)THEN
            IF(LPASS.EQ.0)THEN
               ISUBLI=1
            ELSE
               IPORLI=1
               S1=0
               S2=0
               S3=0
               SS1=0
               SS2=0
               SS3=0
               SRAGMX=0
               SRAGMN=0
               RAINT1=0
               RAINT2=0
               ISCOP=0
               LPASS=0
               IF(SEEN.EQ.0.OR.SEEN.EQ.-1)LPASS=1
               IOLDOL=NLEVEL
            ENDIF
*
         ELSE IF(NLEVEL.LE.IOLDOL)THEN
            IPORLI=1
            S1=0
            S2=0
            S3=0
            SS1=0
            SS2=0
            SS3=0
            SRAGMX=0
            SRAGMN=0
            RAINT1=0
            RAINT2=0
            ISCOP=0
            LPASS=0
            IF(SEEN.EQ.0.OR.SEEN.EQ.-1)LPASS=1
            IOLDOL=NLEVEL
         ENDIF
*
         IF(NLEVEL.LE.IPRELE)THEN
            IF(IVOLNA.NE.IMENO)THEN
               NIET=2
               IF(SEEN.EQ.0.OR.IPORLI.EQ.1.OR.SEEN.EQ.-1)THEN
                  DO 140 I=1,15
                     IF((NLEVEL-I).GE.1)THEN
                        IF(POX(NLEVEL-I).NE.0.OR.POY(NLEVEL-I)
     +                  .NE.0.OR. POZ(NLEVEL-I).NE.0)THEN
                           PO1=POX(NLEVEL-I)
                           PO2=POY(NLEVEL-I)
                           PO3=POZ(NLEVEL-I)
                           GOTO 150
                        ENDIF
                     ENDIF
  140             CONTINUE
                  POX(NLEVEL)=0.
                  POY(NLEVEL)=0.
                  POZ(NLEVEL)=0.
                  PORGX=0.
                  PORGY=0.
                  PORGZ=0.
                  GOTO 160
  150             CONTINUE
                  POX(NLEVEL)=PO1
                  POY(NLEVEL)=PO2
                  POZ(NLEVEL)=PO3
               ENDIF
  160          CONTINUE
            ENDIF
         ENDIF
         IPRELE=NLEVEL
      ENDIF
*
      IF(IOLDOL.EQ.0.AND.(SEEN.EQ.1.OR.SEEN.EQ.-2
     +   .OR.SEEN.EQ.9))THEN
         KSTART=1
         IPORLI=1
         LPASS=0
         IOLDOL=NLEVEL
         IPRELE=NLEVEL
      ENDIF
  170 CONTINUE
*
***SG
*
C
C             WORK attribute enabled ?
C
      IF(WORK.LE.0.)GO TO 270
C
C             SEEN attribute processing
C
      IF (SEEN.LT.50.) GO TO 180
      ISEENL=SEEN/10.+0.5
      SEEN=ISEENL-10
  180 CONTINUE
      IF(NLEVEL.LE.LEVSEE)LEVSEE=1000
      IF(SEEN.EQ.-1.)GO TO 270
      IF (NLEVEL.GT.LEVSEE) GO TO 270
      IF(SEEN.EQ.0.)GO TO 220
      IF (SEEN.EQ.-2.) LEVSEE=NLEVEL
*JS
*
*     Logic has been modified  >>>>>
*
C
C       For the Standard Mode:  Draw the shape
C       For CG Mode : Make a CG-Object for each shape
C
      IF(IHIDEN.EQ.IFLH)THEN
*
*              CG Mode:
*
*
*  Case : divided without clipping
*
****SG
*      Optimization for setting seen attributes
*
         NIN=Q(JVO+3)
         IF(SEEN.EQ.9.AND.NIN.NE.0)THEN
            LPASS=1
            ITSTCU=0
            ICGP=0
            IIIIII=1
            GOTO 190
         ELSE
            IIIIII=0
         ENDIF
         IF(NIN.LT.0.AND.((ISHAPE.GT.1.AND.ISHAPE.LT.5).OR.
     +   (ISHAPE.GT.9.AND.ISHAPE.LT.13)))THEN
            IF(IPORLI.EQ.1)ISCOP=1
         ENDIF
         IF(SEEN.EQ.1.AND.NIN.LT.0)GOTO 190
         IF(SEEN.NE.-2.AND.NIN.LT.0.AND.(IMOD.EQ.0.OR.IMOD.EQ.3))
     +    GOTO 240
  190    CONTINUE
*
#if defined(CERNLIB_CG)
*
*      Creating, clipping and counting CG objects
*      Inserting the visible ones in Hide and Wire Structures
*
         IF(IIIIII.EQ.1)GOTO 200
         IMSE=IMOD
         IF(ISHAPE.EQ.11)LSTEP=GPAR(4,NLEVEL)-1
         IF(ISHAPE.EQ.12)LSTEP=GPAR(3,NLEVEL)-1
         CALL GDCGOB(IMSE,ISHAPE,GPAR(1,NLEVEL),NOBJ,NWWS,IVOLNA,
     +LSTEP)
  200    CONTINUE
         NOBJ=NOBJ+LSTEP
         NUVO=NUVO+1
         LSTEP=1
         IF(IPAS.NE.0)THEN
            IF(NCLAS2.GT.1000)THEN
               IPER=(100*NOBJ)/NCLAS2
               IF(IPER.EQ.10.OR.IPER.EQ.20.OR.IPER.EQ.30.OR.IPER.EQ.
     +         40.OR.IPER.EQ.50.OR.IPER.EQ.60.OR.IPER.EQ.70.OR.IPER
     +         .EQ.80.OR.IPER.EQ.90)THEN
                  IF(IPER.NE.IPEOLD)THEN
                     WRITE(CHMAIL,11800)IPER
                     CALL GMAIL(0,0)
                     IPEOLD=IPER
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
*
*            Setting line attributes volume by volume
*
         IF(IHIDEN.EQ.IFLH) THEN
            IF(IPAS.NE.0.AND.KGG.EQ.1) THEN
**           IF(ITSTCU.NE.0.AND.IVFUN.NE.0) THEN
               IF((ITSTCU.NE.0).OR. ((ISHAPE.EQ.11.OR.ISHAPE.EQ.12)
     +         .AND.(ICGP.NE.0)))THEN
                  IF(ISHAPE.EQ.11.OR.ISHAPE.EQ.12)LSTEP=ICGP
                  LL=II+1
                  II=II+LSTEP
                  DO 210 KHH=LL,II
                     IQ(JCGCOL+KHH)=LINATT
  210             CONTINUE
                  LSTEP=1
               ENDIF
            ENDIF
         ENDIF
*
*
*      Logic has been modified again :
*      do the size evaluation while creating Hide Structure
*      do the same for Wire Structure
*
*      If number of words booked for Hide Structure or for
*      Wire Structure is not sufficient, then evaluate the
*      the right number of words needed and send a mail.
*
         IF(KCGST.EQ.-9)THEN
*
*      Ten words more for safety
*
            NWHS1 = NCLAS1+ 10
            CALL CGHEVA(Q(IHPOIN),HISI)
            IF(HISI.GT.NWHS1)NWHS1=HISI
            IWORH = NWHS1 - NWHS
            IWORH1= IWORH * 1.666666
            WRITE(CHMAIL,11500)IWORH1
            CALL GMAIL(0,0)
            GOTO 320
         ELSE IF(KCGST.EQ.-10)THEN
            NWWS1 = NCLAS3+ 10
            IWORW = NWWS1 - NWWS
            IWORW1= IWORW * 2.5
            WRITE(CHMAIL,11400)IWORW1
            CALL GMAIL(0,0)
            GOTO 320
         ELSEIF(KCGST.EQ.-4.OR.KCGST.EQ.-1.OR.KCGST.EQ.-2) THEN
*      Exiting without having made evaluation of size
            GOTO 320
         ENDIF
****SG
*    Case: Volume placed by GSPOS, not clipped and 'closed'
*
*XX
*         IF(NIN.GT.0.AND.(IMOD.EQ.0.OR.IMOD.EQ.3).
*     +   AND.KSHIFT.EQ.0.AND.GBOOM.EQ.0.)  GOTO 150
*XX
#endif
      ELSE
*
*              Standard Mode: Draw the shape
*
         CALL IGPID(1,'Pick',IVO,' ')
         CALL IGPID(2,'Pick',IQ(JVOLUM+IVO),' ')
         IF (ICUTFL.EQ.1) THEN
            CALL GDRWSC(ISHAPE,GPAR(1,NLEVEL))
         ELSE
            CALL GDRAWS(ISHAPE,GPAR(1,NLEVEL))
         ENDIF
*
      ENDIF
*JS
      JVO=LQ(JVOLUM-IVO)
C
      IF(SEEN.EQ.-2.)GO TO 270
C
  220 CONTINUE
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
      IF(NIN.EQ.0) GO TO 270
      IF(NIN.LT.0) GO TO 240
C
C             Contents placed by GSPOS
C
      IN=0
      IF(NLMAX.GT.NLEVEL) IN=LINDEX(NLEVEL+1)
      IN=IN+1
      IF(IN.GT.NIN.AND.NLEVEL.EQ.NLMIN) GO TO 300
*
      IF(IN.GT.NIN) GO TO 260
*
      CALL GMEPOS(JVO,IN,X,0)
      JIN = LQ(JVO-IN)
*
      NPAR=IQ(JGPAR+NLEVEL)
      DO 230 I=1,NPAR
         GPAR(I,NLEVEL)=Q(LQ(JGPAR-NLEVEL)+I)
  230 CONTINUE
*
      IVO=LVOLUM(NLEVEL)
      JVO=LQ(JVOLUM-IVO)
      NLMAX=NLEVEL
      GO TO 110
C
  240 CONTINUE
C
C             Contents by division
C
      IN=0
      IF(NLMAX.GT.NLEVEL) IN=LINDEX(NLEVEL+1)
      IN=IN+1
      CALL GMEDIV(JVO,IN,X,0)
      JIN = LQ(JVO-IN)
*
      IF (IN.EQ.0) GO TO 260
*
      NPAR=IQ(JGPAR+NLEVEL)
      DO 250 I=1,NPAR
         GPAR(I,NLEVEL)=Q(LQ(JGPAR-NLEVEL)+I)
  250 CONTINUE
*
      IF (IN.EQ.0) GO TO 260
*
      IVO=LVOLUM(NLEVEL)
      JVO=LQ(JVOLUM-IVO)
      NLMAX=NLEVEL
      GO TO 110
C
  260 CONTINUE
      NLMAX=NLEVEL
  270 CONTINUE
      NLEVEL=NLEVEL-1
      IF(NLEVEL.LT.NLMIN) GO TO 300
      IVO=LVOLUM(NLEVEL)
      JVO=LQ(JVOLUM-IVO)
      GO TO 220
C
  280 WRITE(CHMAIL,10000)NAME
      CALL GMAIL(0,0)
      GO TO 300
C
  290 CONTINUE
C
C             TOP OF THE TREE HAS PARAMETERS SET BY GSPOSP.
C             BUT GDRAW DOES NOT HAVE ACCESS TO THE IN BANK
C             WHICH PLACED IT IN ITS MOTHER.
C
      WRITE(CHMAIL,10100) NAME
      CALL GMAIL(0,0)
C
  300 CONTINUE
*
***SG
*
#if defined(CERNLIB_CG)
*
*  CG Mode
*
      IF(IHIDEN.EQ.IFLH) THEN
*
*    In CG Mode the program flow has two scanning of the geometrical tree:
*    the first one is to count the number of volumes (IPAS=0);
*    the second one is to compute volumes' visibility with Hidden Volume
*    and Face Removal, inserting them in the Hide and Wire structures if
*    the size of Zebra store is sufficient (IPAS=1).
*
         IF(IPAS.EQ.0)THEN
*
*    Creating a bank for setting line attributes volume by volume.
*    The number of words needed is just equal to the total number
*      of visible volumes.
*
            NCLAS2=NOBJ
            IF(KGG.EQ.0)THEN
               CALL MZNEED(IXCG,NCLAS2+10,'G')
*
*    Take everything is left but leave 100,000 words just in case
*
               MEMO=IQUEST(11)-(IQUEST(11)*.11)
               IF(IQUEST(11).LE.0)THEN
                  WRITE(CHMAIL,11300)NCLAS2+10
                  CALL GMAIL(0,0)
                  GOTO 320
               ENDIF
               CALL MZBOOK(IXCG,JCGCOL,JCGCOL, 1,'LINE',0,0,NCLAS2+10,
     +         2,-1)
            ENDIF
            KGG=KGG+1
            NOBJ=0
            NUVO=0
            IPAS=1
*
            NFILT=63000+NCLAS2+10
*
*    Use max Zebra store for Hide and Wire structures
*
            NWHS=0.6*MEMO
            NWWS=0.4*MEMO
            CALL MZBOOK(IXCG,JCG,JCG,1,'CG',0,0,MEMO,3,-1)
            CALL GDCGRS
*
            IHS=1
            IHPOIN=JCG+1
*
*    Creating the Hidden Structure
*
            IF(NWHS.LE.LHHEAD)NWHS=LHHEAD+1
***SG
            CALL CGHCRE(NTRCG,0,DUMMY,DUMMY,NWHS,Q(IHPOIN))
            NTCUR=NWHS
            IOLDCU=NTCUR
            NFILT=NFILT+NTCUR
            GOTO 40
         ENDIF
         IF(NOBJ.EQ.0)GOTO 320
****SG
*    Closing the Hidden Structure
*
         IHPOIN=JCG+1
*
*     Last size evaluation for Hide Structure
*
         CALL CGHEND(Q(IHPOIN),IFACST,RSHD)
         IF(IFACST.GT.0)THEN
            NWHS1 = IFACST+ 10
            IWORH = NWHS1 - NWHS
            IWORH1= IWORH * 1.666666
            WRITE(CHMAIL,11500)IWORH1
            CALL GMAIL(0,0)
            GOTO 320
         ENDIF
*
* IFACST shouldn't be negative now
*
         IF(IFACST.LT.0)THEN
            WRITE(CHMAIL,10500)
            CALL GMAIL(0,0)
            GOTO 320
         ENDIF
*
*     Setting the right colours
*               and
*     Drawing the CG Objects
*
         IF(ILCG.LT.IFCG)THEN
            WRITE(CHMAIL,10400)
            CALL GMAIL(0,0)
         ELSE
            IF(LEP.GE.0.AND.ISPFLA.NE.1)THEN
               WRITE(CHMAIL,11900)
               CALL GMAIL(0,0)
            ENDIF
***            call write_dxf_sect_entity( 1 )
            DO 310 K=IFCG,ILCG
               IF(IDVIEW.NE.0)THEN
                  IF(ILCG.GT.1000)THEN
                     IPER=(100*K)/ILCG
                     IF(IPER.EQ.10.OR.IPER.EQ.20.OR.IPER.EQ.30.OR.IPER.
     +               EQ. 40.OR.IPER.EQ.50.OR.IPER.EQ.60.OR.IPER.EQ.70.O
     +               R.IPER .EQ.80.OR.IPER.EQ.90)THEN
                        IF(IPER.NE.IPEOLD)THEN
                           WRITE(CHMAIL,12000)IPER
                           CALL GMAIL(0,0)
                           IPEOLD=IPER
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
               LINATT=IQ(JCGCOL+K-IFCG+1)
               IWPOIN=JCG+IQ(JCOUNT+K)
               IF(ISPFLA.EQ.1)ISPFLA=2
               IF(LLEP.NE.1)THEN
                  MMPOIN=IMPOIN+IQ(IMCOUN+K)
                  CALL GD16V((IQ(JCOUNT+K)-1),IQ(MMPOIN))
               ELSE
                  CALL GD16V((IQ(JCOUNT+K)-1),0)
               ENDIF
               IF(ISPFLA.EQ.2)ISPFLA=1
  310       CONTINUE
         ENDIF
***        call write_dxf_sect_entity( 2 )
*
*    Printing statistics
*
         IF(LEP.GE.0.AND.ISPFLA.NE.1)THEN
            NFILT=NFILT-NWHS+RSHD
            INFILT=NFILT
            IRSHD=RSHD
            MEMO1=MEMO+63000+NCLAS2+10+JMEMT2+JMEMT1+JMEMT3+54000
            RSWR =NFILT-RSHD-63000-NCLAS2-10
            IRSWR=RSWR
            ICGOB=63000
            RATIO=RSHD/RSWR
            JJIIKK=0
            JPARGE=0
            IF(LLEP.NE.1)THEN
               JJIIKK=ILCG-IFCG+1
               JPARGE=54000
            ENDIF
            WRITE(CHMAIL,10600)INFILT+NTNEX+ILCG-IFCG+1+JJIIKK+JPARGE
            CALL GMAIL(0,0)
            WRITE(CHMAIL,10700)MEMO1
            CALL GMAIL(0,0)
            WRITE(CHMAIL,10900)NCLAS2+10
            CALL GMAIL(0,0)
            WRITE(CHMAIL,10800)ICGOB
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11000)IRSHD
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11100)IRSWR+ILCG-IFCG+1
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11200)NTNEX+JJIIKK+JPARGE
            CALL GMAIL(0,0)
*         WRITE(CHMAIL,10799)RATIO
*         CALL GMAIL(0,0)
            WRITE(CHMAIL,11600)ILCG-IFCG+1
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11700)NCLAS2
            CALL GMAIL(0,0)
         ENDIF
*
      ENDIF
*
*    Dropping + resetting parameters
*
  320 CONTINUE
      IF(IHIDEN.EQ.IFLH)THEN
         ICUT=0
         IF(JCG.NE.0)CALL MZDROP(IXSTOR,JCG,' ')
         IF(JCGOBJ.NE.0)CALL MZDROP(IXSTOR,JCGOBJ,' ')
         CALL GDCGRS
         IF(JCGCOL.NE.0)CALL MZDROP(IXSTOR,JCGCOL,' ')
         IF(JCOUNT.NE.0)CALL MZDROP(IXSTOR,JCOUNT,' ')
         IF(JCLIPS.NE.0)CALL MZDROP(IXSTOR,JCLIPS,' ')
         IF(IMPOIN.NE.0)CALL MZDROP(IXSTOR,IMPOIN,' ')
         IF(IMCOUN.NE.0)CALL MZDROP(IXSTOR,IMCOUN,' ')
         IF(JSIX.NE.0) CALL MZDROP(IXSTOR, JSIX, ' ')
         IF(JSIY.NE.0) CALL MZDROP(IXSTOR, JSIY, ' ')
         IF(JSIZ.NE.0) CALL MZDROP(IXSTOR, JSIZ, ' ')
         IF(JPXC.NE.0) CALL MZDROP(IXSTOR, JPXC, ' ')
         IF(JPYC.NE.0) CALL MZDROP(IXSTOR, JPYC, ' ')
         IF(JPZC.NE.0) CALL MZDROP(IXSTOR, JPZC, ' ')
         LARECG(1)=0
         CALL MZGARB(IXSTOR+1,0)
*         NWHS1=0
*         NWFLAG=0
         NCLAS1=0
         NCLAS2=0
         NCLAS3=0
      ENDIF
*
****SG
*
#endif
      IF (IFCVOL.EQ.1) THEN
         CALL GFCVOL
      ELSE
         NLEVEL=0
      ENDIF
C
C             If in cut-mode then close the GDRAWV line buffer
C
      IF (ICUTFL.EQ.1) CALL GDRAWV(0.,0.,0)
C
C             Restore permanent value of color and return
C
*      CALL GDCOL(0)
      IOBJ=0
C
10000 FORMAT(' *** GDRAW *** : Volume ',A4,' does not exist')
10100 FORMAT(' *** GDRAW *** : Top of tree ',A4,' parameters defined',
     +       '  by GSPOSP - info not available to GDRAW.')
10200 FORMAT(' *** GDRAW *** : Illegal Transformation Matrix',
     +       ' Number NTRCG ')
10300 FORMAT(' *** GDRAW *** :  >>> Det (T) = 0  ')
10400 FORMAT(' *** GDRAW *** :  Warning! Volume is destroyed.')
*SG
10500 FORMAT(' *** GDRAW *** :  Internal error, please report to',
     +       ' GEANT support team')
10600 FORMAT(' *** GDRAW *** : Total memory used    =',I10,' words.')
10700 FORMAT(' *** GDRAW *** : Total memory booked  =',I10,' words.')
10800 FORMAT(' *** GDRAW *** : Memory used for CGOB =',I10,' words.')
10900 FORMAT(' *** GDRAW *** : Memory used for LINE =',I10,' words.')
11000 FORMAT(' *** GDRAW *** : Memory used for HIDE =',I10,' words.')
11100 FORMAT(' *** GDRAW *** : Memory used for WIRE =',I10,' words.')
11200 FORMAT(' *** GDRAW *** : Memory used for SHAD =',I10,' words.')
11300 FORMAT(' *** GDRAW *** : Memory needed for the LINE attributes',
     +       ' =',I10,' words.')
*10799 FORMAT(' *** GDRAW *** : HIDE/WIRE=',F4.2,'.')
11400 FORMAT(' *** GDRAW *** : Please, increase size of Zebra store',
     +       ' by',I10,' words to create WIRE structure.')
11500 FORMAT(' *** GDRAW *** : Please, increase size of Zebra store',
     +       ' by',I10,' words to create HIDE structure.')
11600 FORMAT(' *** GDRAW *** : Visible volumes      =',I10,'.')
11700 FORMAT(' *** GDRAW *** : Total   volumes      =',I10,'.')
11800 FORMAT(' *** GDRAW *** :',I2,'% of volumes analysed.')
11900 FORMAT(' *** GDRAW *** : Now the drawing is starting !')
12000 FORMAT(' *** GDRAW *** :',I2,'% of volumes drawn.')
*
*SG
      END
