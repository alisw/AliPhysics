*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GRIN(CHOBJT,IDVERS,CHOPT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to read GEANT object(s) fromin the RZ file       *
C.    *         at the Current Working Directory (See RZCDIR)          *
C.    *       The data structures from disk are read in memory         *
C.    *           (VOLU,ROTM,TMED,MATE,SETS,PART,SCAN)                 *
C.    *                                                                *
C.    *       CHOBJ  The type of object to be read:                    *
C.    *              MATE read JMATE structure                         *
C.    *              TMED read JTMED structure                         *
C.    *              VOLU read JVOLUM structure                        *
C.    *              ROTM read JROTM structure                         *
C.    *              SETS read JSET  structure                         *
C.    *              PART read JPART structure                         *
C.    *              SCAN read LSCAN structure                         *
C.    *              INIT read all initialisation structures           *
C.    *                                                                *
C.    *       IDVERS is a positive integer which specifies the version *
C.    *           number of the object(s).                             *
C.    *                                                                *
C.    *       CHOPT List of options (none for the time being)          *
C.    *                                                                *
C.    *                                                                *
C.    *    The RZ data base has been created via GRFILE/GROUT          *
C.    *                                                                *
C.    *                                                                *
C.    *      Example.                                                  *
C.    *                                                                *
C.    *      CALL GRFILE(1,'Geometry.dat',' ')                         *
C.    *      CALL GRIN ('VOLU',1,' ')                                  *
C.    *      CALL GRIN ('MATE',1,' ')                                  *
C.    *      CALL GRIN ('TMED',1,' ')                                  *
C.    *      CALL GRIN ('ROTM',1,' ')                                  *
C.    *      CALL GRIN ('PART',1,' ')                                  *
C.    *      CALL GRIN ('SCAN',1,' ')                                  *
C.    *      CALL GRIN ('SETS',1,' ')                                  *
C.    *                                                                *
C.    *      The same result can be achieved by:                       *
C.    *      CALL GRFILE(1,'Geometry.dat','I')                         *
C.    *                                                                *
C.    *    ==>Called by : <USER> ,GRFILE                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcnum.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctime.inc"
*      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
*     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
*     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
      COMMON/QUEST/IQUEST(100)
      PARAMETER (NLINIT=9,NLKINE=2,NLTRIG=6,NMKEY=22)
      DIMENSION JNAMES(20),KEYS(2),LINIT(NLINIT),LKINE(NLKINE)
      DIMENSION LTRIG(NLTRIG),IXD(NMKEY)
      DIMENSION LINK(NMKEY),IVERSI(NMKEY),LDIV(2),IRESUL(NMKEY)
      DIMENSION IDOLD(8), IDNEW(8), VEROLD(8), VERNEW(8)
      EQUIVALENCE (JNAMES(1),JDIGI)
      CHARACTER*4 CHOBJ,CKEY,KNAMES(NMKEY)
      CHARACTER*(*) CHOPT,CHOBJT
      DATA KNAMES/'DIGI','DRAW','HEAD','HITS','KINE','MATE','PART',
     +     'ROTM','RUNG','SETS','STAK','STAT','TMED','NULL','VERT',
     +     'VOLU','JXYZ','NULL','NULL','NULL','SCAN','NULL'/
      DATA IXD/2,1,2,2,2,8*1,2,2,1,2,3*0,1,0/
      DATA LINIT/2,6,7,8,9,10,13,16,21/
      DATA LKINE/5,15/
      DATA LTRIG/1,3,4,5,15,17/
      DATA IDNEW / 8*0 /
      DATA VERNEW / 8*0. /
C.
C.    ------------------------------------------------------------------
C.
      IQUEST(1)=0
      LDIV(1)  =IXCONS
      LDIV(2)  =IXDIV
      KVOL=JVOLUM
      CHOBJ=CHOBJT
*
*     Save old JRUNG dates and versions
      IF(JRUNG.GT.0) THEN
         DO 10 J=1,8
            IDOLD(J) = IQ(JRUNG+10+J)
            VEROLD(J) = Q(JRUNG+20+J)
   10    CONTINUE
      ENDIF
*
      IOPTI=INDEX(CHOPT,'i')+INDEX(CHOPT,'I')
      IOPTT=INDEX(CHOPT,'t')+INDEX(CHOPT,'T')
      IOPTK=INDEX(CHOPT,'k')+INDEX(CHOPT,'K')
      IOPTQ=INDEX(CHOPT,'q')+INDEX(CHOPT,'Q')
*
      IF(CHOBJ.EQ.'INIT') THEN
         CHOBJ='*'
         IOPTI=1
         IOPTT=0
         IOPTK=0
         IF(JDRAW.NE.0)CALL MZDROP(IXSTOR,JDRAW,'L')
         IF(JSTAK.NE.0)CALL MZDROP(IXSTOR,JSTAK,'L')
         IF(JGPAR.NE.0)CALL MZDROP(IXSTOR,JGPAR,'L')
      ELSEIF(CHOBJ.EQ.'TRIG') THEN
         CHOBJ='*'
         IOPTI=0
         IOPTT=1
         IOPTK=0
      ELSEIF(CHOBJ.EQ.'KINE') THEN
         CHOBJ='*'
         IOPTI=0
         IOPTT=0
         IOPTK=1
      ENDIF
*
      IF(CHOBJ.EQ.'*') THEN
         IF(IOPTI.NE.0) THEN
            DO 20 J=1, NLINIT
               LINK(J)=LINIT(J)
   20       CONTINUE
            NLINK=NLINIT
         ELSEIF(IOPTT.NE.0) THEN
            DO 30 J=1, NLTRIG
               LINK(J)=LTRIG(J)
   30       CONTINUE
            NLINK=NLTRIG
         ELSEIF(IOPTK.NE.0) THEN
            DO 40 J=1, NLKINE
               LINK(J)=LKINE(J)
   40       CONTINUE
            NLINK=NLKINE
         ENDIF
      ELSE
         NLINK=0
         DO 100 J=1, NMKEY
            IF(CHOBJ.EQ.KNAMES(J)) THEN
               IF(IOPTI.NE.0) THEN
                  DO 50 L=1, NLINIT
                     IF(LINIT(L).EQ.J) GOTO 80
   50             CONTINUE
                  GOTO 90
               ELSEIF(IOPTT.NE.0) THEN
                  DO 60 L=1, NLTRIG
                     IF(LTRIG(L).EQ.J) GOTO 80
   60             CONTINUE
                  GOTO 90
               ELSEIF(IOPTK.NE.0) THEN
                  DO 70 L=1, NLKINE
                     IF(LKINE(L).EQ.J) GOTO 80
   70             CONTINUE
                  GOTO 90
               ENDIF
   80          NLINK=1
               LINK(1)=J
               GOTO 110
*
   90          WRITE(CHMAIL,10000) CHOBJ, CHOPT
               CALL GMAIL(0,0)
               GOTO 999
*
            ENDIF
  100    CONTINUE
      ENDIF
*
  110 IF(NLINK.EQ.0) THEN
         WRITE(CHMAIL,10100) CHOBJ
         CALL GMAIL(0,0)
         GOTO 999
      ENDIF
*
      DO 120 J=1, NLINK
         IVERSI(J)=0
         IRESUL(J)=0
  120 CONTINUE
*
      IKEY=0
  130 CONTINUE
      IKEY=IKEY+1
      CALL RZINK(IKEY,9999,'S')
      IF(IQUEST(1).NE.0) THEN
         IQUEST(1)=0
         GOTO 150
      ENDIF
      CALL UHTOC(IQUEST(21),4,CKEY,4)
      DO 140 I=1,NLINK
         NKEY=ABS(LINK(I))
         IF(CKEY.EQ.KNAMES(NKEY))THEN
            IF(IDVERS.NE.0.AND.IDVERS.NE.IQUEST(22)) THEN
               LINK(I)=-ABS(LINK(I))
               GOTO 130
            ENDIF
            KEYS(1)=IQUEST(21)
            KEYS(2)=IQUEST(22)
            IDIV=LDIV(IXD(NKEY))
            IF(NKEY.LE.20)THEN
               IF(JNAMES(NKEY).NE.0)THEN
                  CALL MZDROP(IDIV,JNAMES(NKEY),'L')
                  JNAMES(NKEY)=0
               ENDIF
               CALL RZIN(IDIV,JNAMES(NKEY),1,KEYS,9999,' ')
            ELSE
               NKL=NKEY-20
               IF(ISLINK(NKL).NE.0)THEN
                  CALL MZDROP(IDIV,ISLINK(NKL),'L')
                  ISLINK(NKL)=0
               ENDIF
               CALL RZIN(IDIV,ISLINK(NKL),1,KEYS,9999,' ')
            ENDIF
            IF(IQUEST(1).EQ.0) THEN
               IVERSI(I)=IQUEST(22)
               IRESUL(I)=1
            ENDIF
         ENDIF
  140 CONTINUE
      GOTO 130
*
  150 NIN=0
      DO 160 I=1,NLINK
         IF(IRESUL(I).EQ.1) THEN
            IF(IOPTQ.EQ.0) THEN
            WRITE(CHMAIL,10200) KNAMES(LINK(I)), IVERSI(I)
            CALL GMAIL(0,0)
            ENDIF
            NIN=NIN+1
         ELSEIF(LINK(I).GT.0) THEN
            IF(IOPTQ.EQ.0) THEN
            WRITE(CHMAIL,10300) KNAMES(LINK(I))
            CALL GMAIL(0,0)
            ENDIF
         ELSEIF(LINK(I).LT.0) THEN
            IF(IOPTQ.EQ.0) THEN
            WRITE(CHMAIL,10400) KNAMES(-LINK(I)), IDVERS
            CALL GMAIL(0,0)
            ENDIF
         ENDIF
  160 CONTINUE
*
      IF(NIN.EQ.0) THEN
         WRITE(CHMAIL,10500)
         CALL GMAIL(0,0)
         GOTO 999
      ENDIF
*
      IF(KVOL.NE.JVOLUM)THEN
         NVOLUM=IQ(JVOLUM-1)
         CALL MZGARB(IXCONS,0)
         CALL GGDVLP
         CALL GGNLEV
      ENDIF
*
      IF(JVOLUM.GT.0) THEN
         NLEVEL=0
         NVOLUM=0
         DO 170 J=1, IQ(JVOLUM-2)
            IF(LQ(JVOLUM-J).EQ.0) GOTO 180
            NVOLUM=NVOLUM+1
  170    CONTINUE
  180    CONTINUE
      ENDIF
*
      IF(JTMED.NE.0 )THEN
         CALL UCOPY(Q(JTMED+1),CUTGAM,10)
         NTMED=IQ(JTMED-2)
      ENDIF
*
      IF(JPART.NE.0 ) NPART  = IQ(JPART-2)
      IF(JVERTX.NE.0) NVERTX = IQ(JVERTX-2)
      IF(JMATE.NE.0 ) NMATE  = IQ(JMATE-2)
      IF(JROTM.NE.0 ) NROTM  = IQ(JROTM-2)
      IF(JDRAW.GT.0 ) THEN
         NKVIEW = IQ(JDRAW-2)
      ELSE
         NKVIEW = 0
C
C             Book JDRAW structure for view banks
C
         CALL MZBOOK(IXCONS,JDRAW,JDRAW,1,'DRAW',0,0,0,3,0)
      ENDIF
 
C
      IF(JHEAD.GT.0)THEN
         IDRUN=IQ(JHEAD+1)
         IDEVT=IQ(JHEAD+2)
      ENDIF
      IF(JRUNG.GT.0) THEN
*
*             Here we deal with version numbers If JRUNG has been read in,
*             then save the version numbers of the new JRUNG and restore
*             the current version number for KINE, HITS and DIGI
         DO 210 J=1, NLINK
            IF(IRESUL(J).EQ.1) THEN
               NKEY = ABS(LINK(J))
               IF(KNAMES(NKEY).EQ.'RUNG') THEN
                  DO 190 I=1,8
                     IDNEW(I) = IQ(JRUNG+10+I)
                     VERNEW(I) = Q(JRUNG+20+I)
  190             CONTINUE
*
*             And we put back the old version numbers because,
*             in principle, KINE, HITS and DIGI have not be read in
                  DO 200 I=3,8
                     IQ(JRUNG+10+I) = IDOLD(I)
                     Q(JRUNG+20+I) = VEROLD(I)
  200             CONTINUE
               ENDIF
            ENDIF
  210    CONTINUE
*
*            And here we do it again for KINE, HITS and DIGI
         DO 220 J=1, NLINK
            IF(IRESUL(J).EQ.1) THEN
               NKEY = ABS(LINK(J))
               IF(KNAMES(NKEY).EQ.'KINE') THEN
                  IF(IDNEW(3).GT.0) THEN
                     IQ(JRUNG+13) = IDNEW(3)
                     IQ(JRUNG+14) = IDNEW(4)
                     Q(JRUNG+23) = VERNEW(3)
                     Q(JRUNG+24) = VERNEW(4)
                  ENDIF
               ELSEIF(KNAMES(NKEY).EQ.'HITS') THEN
                  IF(IDNEW(5).GT.0) THEN
                     IQ(JRUNG+15) = IDNEW(5)
                     IQ(JRUNG+16) = IDNEW(6)
                     Q(JRUNG+25) = VERNEW(5)
                     Q(JRUNG+26) = VERNEW(6)
                  ENDIF
               ELSEIF(KNAMES(NKEY).EQ.'DIGI') THEN
                  IF(IDNEW(7).GT.0) THEN
                     IQ(JRUNG+17) = IDNEW(7)
                     IQ(JRUNG+18) = IDNEW(8)
                     Q(JRUNG+27) = VERNEW(7)
                     Q(JRUNG+28) = VERNEW(8)
                  ENDIF
               ELSEIF(KNAMES(NKEY).EQ.'MATE'.OR.
     +                KNAMES(NKEY).EQ.'TMED') THEN
                  IF(VERNEW(1).NE.0) THEN
*                We know which version number we are reading
                     IF(VERNEW(1).LT.GVERSN) THEN
                        WRITE(CHMAIL,10600)
     +                        KNAMES(NKEY),VERNEW(1),GVERSN
                        CALL GMAIL(0,0)
                        WRITE(CHMAIL,10700)
                        CALL GMAIL(0,0)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  220    CONTINUE
      ENDIF
*
10000 FORMAT(' *** GRIN *** Data structure ',A4,' not read in ',
     +       'phase ',A)
10100 FORMAT(' *** GRIN *** Unknown key ',A4)
10200 FORMAT(' *** GRIN *** Data structure ',A4,' version ',I10,
     +       ' successfully read in ')
10300 FORMAT(' *** GRIN *** Data structure ',A4,' was not found')
10400 FORMAT(' *** GRIN *** Data structure ',A4,' version ',I10,
     +       ' was not found')
10500 FORMAT(' *** GRIN *** Nothing found to read !')
10600 FORMAT(' *** GRIN *** ',A4,' data structure ',
     +       'version ',F6.4,' current version is ',F6.4)
10700 FORMAT('              Please call subroutine GPHYSI before ',
     +       'tracking')
  999 END
