*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 18/01/95  19.56.13  by  S.Giani
*-- Author :
      SUBROUTINE GLVOLU (NLEV, LNAM, LNUM, IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Loads the common block GCVOLU for the volume at lebel NLEV   *
C.    *   as described by the lists of names (LNAM) and numbers (LNUM) *
C.    *                                                                *
C.    *   The routine is optimized and does not re-compute the part of *
C.    *   history already available in GCVOLU.                         *
C.    *                                                                *
C.    *   IER returns non zero in case of fatal error                  *
C.    *                                                                *
C.    *   Called by : 'User', GDRVOL                                   *
C.    *   Authors   : S.Banerjee, F.Bruyant, A.McPherson               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
      PARAMETER  (NLVMAX=15)
      INTEGER    LNUM(*), LNAM(*), IDTYP(3,12)
      DIMENSION  LVOLS(NLVMAX), LINDX(NLVMAX)
      REAL       XC(3)
      CHARACTER*4 KNAME
      SAVE IDTYP
C.
      DATA IDTYP / 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 2, 3, 1,
     +             2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 4, 3, 1, 1, 1,
     +             2, 3, 1, 2, 3, 1/
C.
C.    ------------------------------------------------------------------
*
      IER = 0
      NLEVL=NLEV
      IF (NLEVL.LE.0.OR.NLEVL.GT.NLVMAX) GO TO 910
      IF (JGPAR.EQ.0) THEN
         CALL MZBOOK (IXCONS, JGPAR, JGPAR, 1, 'GPAR', NLVMAX, 0,
     +                NLVMAX, 2, 0)
      ENDIF
      IF (NLEVEL.EQ.0)                 GO TO 20
*
* *** Scan tree from top to bottom to
*       check if some levels are already filled
*
      NLMX = MIN (NLEVL, NLEVEL)
      NLEVEL = 0
      DO 10 I = 1,NLMX
         IF (LNAM(I).NE.NAMES(I))  GO TO 15
         IF (LNUM(I).NE.NUMBER(I)) GO TO 15
         NLEVEL = NLEVEL +1
   10 CONTINUE
      IF (NLEVL.GT.NLEVEL) GO TO 95
      nlevel=0
   15 IF (NLEVEL.NE.0)    GO TO 95
*
* *** Special case, first volume
*
   20 IF (JVOLUM.EQ.0)    GO TO 920
      IF (IQ(JVOLUM+1).EQ.LNAM(1)) THEN
         IVO = 1
      ELSE
#if defined(CERNLIB_DEBUG)
         WRITE (CHMAIL, 7000) LNAM(1)
         CALL GMAIL (0, 0)
#endif
         IF (IQ(JVOLUM-1).LE.1) GO TO 920
         DO 25 IV=2,IQ(JVOLUM-1)
            IF (IQ(JVOLUM+IV).EQ.LNAM(1)) THEN
               IVO = IV
               GO TO 30
            ENDIF
   25    CONTINUE
         WRITE (CHMAIL, 8000) LNAM(1)
         CALL GMAIL (0, 0)
         GO TO 999
      ENDIF
   30 NLEVEL = 1
      JVO = LQ(JVOLUM-IVO)
      LVOLUM(NLEVEL) = IVO
      NAMES(NLEVEL)  = IQ(JVOLUM+IVO)
      NUMBER(NLEVEL) = LNUM(1)
      GONLY(NLEVEL)  = 1.
      IF (LQ(JVO).EQ.0) THEN
         NLDEV(1) = NLVMAX
      ELSE
         NLDEV(1) = 1
      ENDIF
      IF (IVO.EQ.1) THEN
         LINDEX(NLEVEL) = 1
         LINMX (NLEVEL) = 1
         NLDM = 0
         IQ(JGPAR+NLEVEL) = Q(JVO+5)
         LQ(JGPAR-NLEVEL) = JVO + 6
      ELSE
         CALL UHTOC(NAMES,4,KNAME,4)
         CALL GLMOTH (KNAME, NUMBER, NLDM, LVOLS, LINDX)
         IF (NLDM.GT.0) THEN
            JVOM = LQ(JVOLUM-LVOLS(NLDM))
            NIN  = Q(JVOM+3)
            IF (NIN.LT.0) THEN
               LINDEX(NLEVEL) = LNUM(1)
            ELSE
               LINMX(NLEVEL)  = NIN
               DO 70 IN = 1, NIN
                  JIN  = LQ(JVOM-IN)
                  IF (IFIX(Q(JIN+2)).NE.LVOLUM(1)) GO TO 70
                  IF (IFIX(Q(JIN+3)).NE.LNUM(1))   GO TO 70
                  LINDEX(NLEVEL) = IN
                  GO TO 75
   70          CONTINUE
               GO TO 920
            ENDIF
   75       JPAR = LQ(LQ(JVOLUM-LVOLS(1)))
            IF (NLDM.GT.1) THEN
               DO 76 ILEV = 2, NLDM
                  IF (IQ(JPAR+1).EQ.0) THEN
                     JPAR = LQ(JPAR-LINDX(ILEV))
                     IF (JPAR.EQ.0) GO TO 77
                  ELSE IF (IQ(JPAR-3).GT.1) THEN
                     JPAR = LQ(JPAR-LINDX(ILEV))
                  ELSE
                     JPAR = LQ(JPAR-1)
                  ENDIF
   76          CONTINUE
            ENDIF
            IF (NIN.GT.0) THEN
               JPAR = LQ(JPAR-IN)
               IF (JPAR.EQ.0) GO TO 77
            ELSE
               NDIV = IQ(JPAR+1)
               LINMX(NLEVEL) = NDIV
               IF (LINDEX(1).GT.NDIV) THEN
                  NL1  = 1
                  NAME = IQ(JVOLUM+LVOLS(NLDM))
                  GO TO 950
               ENDIF
               IF (IQ(JPAR-3).GT.1) THEN
                  IF (LINDEX(1).GT.0) THEN
                     JPAR = LQ(JPAR-LINDEX(1))
                  ELSE
                     JPAR = LQ(JPAR-1)
                  ENDIF
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
            ENDIF
            IQ(JGPAR+NLEVEL) = IQ(JPAR+5)
            LQ(JGPAR-NLEVEL) = JPAR + 5
            GO TO 78
   77       NPAR = Q(JVO+5)
            IF (NPAR.EQ.0.AND.NIN.GT.0) THEN
               IQ(JGPAR+NLEVEL) = Q(JIN+9)
               LQ(JGPAR-NLEVEL) = JIN+9
            ELSE
               IQ(JGPAR+NLEVEL) = NPAR
               LQ(JGPAR-NLEVEL) = JVO + 6
            ENDIF
         ELSE
            LINDEX(NLEVEL) = 1
            LINMX(NLEVEL)  = 1
            IQ(JGPAR+NLEVEL) = Q(JVO+5)
            LQ(JGPAR-NLEVEL) = JVO + 6
         ENDIF
      ENDIF
   78 CONTINUE
*
      DO 90 I = 1,3
         GTRAN(I,1) = 0.
         DO 80 J = 1,3
            K = (I-1)*3 +J
            GRMAT(K,1) = 0.
   80    CONTINUE
         K = I*4 -3
         GRMAT(K,1) = 1.
   90 CONTINUE
      GRMAT(10,1) = 0.
      IF (NLEVL.GT.1) THEN
         GO TO 100
      ELSE
         GO TO 990
      ENDIF
*
* *** Check if there are volumes up in the tree where development
*           structure exists
*
   95 IF (LVOLUM(1).EQ.1.OR.NLDEV(1).EQ.1) THEN
         NLDM = 0
      ELSE
         CALL UHTOC(NAMES,4,KNAME,4)
         CALL GLMOTH (KNAME, NUMBER, NLDM, LVOLS, LINDX)
      ENDIF
*
*  ** Next level
*
  100 CONTINUE
      IVO = LVOLUM(NLEVEL)
      JVO = LQ(JVOLUM-IVO)
      NLD = NLDEV(NLEVEL)
      NIN = Q(JVO+3)
      IF (NIN.EQ.0) GO TO 930
      NL1 = NLEVEL +1
*
      IF (NIN.GT.0) THEN
*
*  *     Content obtained by positioning
*
         DO 110 IN=1,NIN
            JIN=LQ(JVO-IN)
            IVOT=Q(JIN+2)
            IF (IQ(JVOLUM+IVOT).NE.LNAM(NL1)) GO TO 110
            INUM = Q(JIN+3)
            IF (INUM.EQ.LNUM(NL1)) GO TO 115
  110    CONTINUE
         GO TO 940
  115    IF (NLEVEL.GE.NLD) THEN
*           (case with JVOLUM structure locally developed)
            JPAR = LQ(LQ(JVOLUM-LVOLUM(NLD)))
            DO 120 ILEV = NLD, NLEVEL
               IF (IQ(JPAR+1).EQ.0) THEN
                  IF (ILEV.EQ.NLEVEL) THEN
                     JPAR = LQ(JPAR-IN)
                  ELSE
                     JPAR = LQ(JPAR-LINDEX(ILEV+1))
                  ENDIF
                  IF (JPAR.EQ.0) GO TO 125
               ELSE IF (IQ(JPAR-3).GT.1) THEN
                  JPAR = LQ(JPAR-LINDEX(ILEV+1))
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
  120       CONTINUE
            JPAR = JPAR + 5
            NPAR = IQ(JPAR)
            GO TO 130
         ELSE IF (NLDM.GT.0) THEN
            JPAR = LQ(LQ(JVOLUM-LVOLS(1)))
            IF (NLDM.GT.1) THEN
               DO 121 ILEV = 2, NLDM
                  IF (IQ(JPAR+1).EQ.0) THEN
                     JPAR = LQ(JPAR-LINDX(ILEV))
                     IF (JPAR.EQ.0)   GO TO 125
                  ELSE IF (IQ(JPAR-3).GT.1) THEN
                     JPAR = LQ(JPAR-LINDX(ILEV))
                  ELSE
                     JPAR = LQ(JPAR-1)
                  ENDIF
  121          CONTINUE
            ENDIF
            DO 122 ILEV = 1, NL1
               IF (IQ(JPAR+1).EQ.0) THEN
                  IF (ILEV.EQ.NL1) THEN
                     JPAR = LQ(JPAR-IN)
                  ELSE
                     JPAR = LQ(JPAR-LINDEX(ILEV))
                  ENDIF
                  IF (JPAR.EQ.0) GO TO 125
               ELSE IF (IQ(JPAR-3).GT.1) THEN
                  JPAR = LQ(JPAR-LINDEX(ILEV))
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
  122       CONTINUE
            JPAR = JPAR + 5
            NPAR = IQ(JPAR)
            GO TO 130
         ENDIF
*        (normal case)
  125    JVOT = LQ(JVOLUM-IVOT)
         NPAR = Q(JVOT+5)
         IF (NPAR.EQ.0) THEN
            JPAR = JIN + 9
            NPAR = Q(JPAR)
         ELSE
            JPAR = JVOT + 6
         ENDIF
*
  130    IROTT = Q(JIN+4)
         NINSK = NIN
         GONLY(NL1) = Q(JIN+8)
         CALL GTRMUL (GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), Q(JIN+5), IROTT
     +,               GTRAN(1,NL1), GRMAT(1,NL1))
*
      ELSE
*
*  *     This section for divided objects
*
         JDIV = LQ(JVO-1)
         IVOT = Q(JDIV+2)
         IF (LNAM(NL1).NE.IQ(JVOLUM+IVOT)) GO TO 960
         JVOT = LQ(JVOLUM-IVOT)
         IF (NLEVEL.GT.NLD) THEN
*           (case with JVOLUM structure locally developed)
            JPAR = LQ(LQ(JVOLUM-LVOLUM(NLD)))
            DO 135 ILEV = NLD, NLEVEL-1
               IF (IQ(JPAR+1).EQ.0) THEN
                  JPAR = LQ(JPAR-LINDEX(ILEV+1))
                  IF (JPAR.EQ.0) GO TO 140
               ELSE IF (IQ(JPAR-3).GT.1) THEN
                  JPAR = LQ(JPAR-LINDEX(ILEV+1))
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
               IF (ILEV.EQ.NLEVEL-1) THEN
                  NDIV = IQ(JPAR+1)
                  ORIG =  Q(JPAR+2)
                  STEP =  Q(JPAR+3)
               ENDIF
  135       CONTINUE
            GO TO 145
         ELSE IF (NLD.EQ.NLEVEL) THEN
            JPAR = LQ(LQ(JVOLUM-LVOLUM(NLD)))
         ELSE IF (NLDM.GT.0) THEN
            JPAR = LQ(LQ(JVOLUM-LVOLS(1)))
            IF (NLDM.GT.1) THEN
               DO 136 ILEV = 2, NLDM
                  IF (IQ(JPAR+1).EQ.0) THEN
                     JPAR = LQ(JPAR-LINDX(ILEV))
                     IF (JPAR.EQ.0) GO TO 140
                  ELSE IF (IQ(JPAR-3).GT.1) THEN
                     JPAR = LQ(JPAR-LINDX(ILEV))
                  ELSE
                     JPAR = LQ(JPAR-1)
                  ENDIF
  136          CONTINUE
            ENDIF
            DO 137 ILEV = 1, NLEVEL
               IF (IQ(JPAR+1).EQ.0) THEN
                  JPAR = LQ(JPAR-LINDEX(ILEV))
                  IF (JPAR.EQ.0) GO TO 140
               ELSE IF (IQ(JPAR-3).GT.1) THEN
                  JPAR = LQ(JPAR-LINDEX(ILEV))
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
               IF (ILEV.EQ.NLEVEL) THEN
                  NDIV = IQ(JPAR+1)
                  ORIG =  Q(JPAR+2)
                  STEP =  Q(JPAR+3)
               ENDIF
  137       CONTINUE
            GO TO 145
         ELSE
            JPAR = 0
         ENDIF
*        (normal case)
  140    NDIV = Q(JDIV+3)
         ORIG = Q(JDIV+4)
         STEP = Q(JDIV+5)
  145    IN   = LNUM(NL1)
         IF (IN.LT.1.OR.IN.GT.NDIV) THEN
            NAME = NAMES(NLEVEL)
            GO TO 950
         ENDIF
*
         IF (JPAR.NE.0) THEN
            IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-IN)
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
            JPAR = JPAR + 5
            NPAR = IQ(JPAR)
         ELSE
            NPAR = Q(JVOT+5)
            JPAR = JVOT + 6
         ENDIF
         GONLY(NL1) = GONLY(NLEVEL)
*
         IAXIS = Q(JDIV+1)
         ISH   = Q(JVO+2)
         IDT   = IDTYP(IAXIS,ISH)
         NINSK = NDIV
*
         IF (IDT.EQ.1) THEN
            DO 151 I = 1, 3
  151       XC(I) = 0.
            XC(IAXIS) = ORIG + (IN - 0.5) * STEP
            IF (ISH.EQ.4.OR.(ISH.EQ.10.AND.IAXIS.NE.1)) THEN
               CALL GCENT (IAXIS, XC)
            ENDIF
            IF (GRMAT(10,NLEVEL).EQ.0.0) THEN
               DO 152 I = 1, 3
  152          GTRAN(I,NL1) = GTRAN(I,NLEVEL)+XC(I)
               DO 153 I = 1, 10
  153          GRMAT(I,NL1) = GRMAT(I,NLEVEL)
            ELSE
               CALL GTRMUL (GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), XC, 0,
     +                      GTRAN(1,NL1), GRMAT(1,NL1))
            ENDIF
*
         ELSE IF (IDT.EQ.3.OR.IDT.EQ.4) THEN
            IF (IDT.EQ.3) THEN
               PH0  = DEGRAD * (ORIG + (IN - 0.5) * STEP)
               CPHR = COS (PH0)
               SPHR = SIN (PH0)
            ELSE
               PH0  = 0.0
               CPHR = 1.0
               SPHR = 0.0
            ENDIF
            DO 154 I = 1, 3
               GTRAN(I  ,NL1) = GTRAN(I  ,NLEVEL)
               GRMAT(I  ,NL1) = GRMAT(I  ,NLEVEL)*CPHR
     +                        + GRMAT(I+3,NLEVEL)*SPHR
               GRMAT(I+3,NL1) = GRMAT(I+3,NLEVEL)*CPHR
     +                        - GRMAT(I  ,NLEVEL)*SPHR
               GRMAT(I+6,NL1) = GRMAT(I+6,NLEVEL)
  154       CONTINUE
            IF (PH0.EQ.0.0.AND.GRMAT(10,NLEVEL).EQ.0.0) THEN
               GRMAT(10,NL1) = 0.0
            ELSE
               GRMAT(10,NL1) = 1.0
            ENDIF
*
         ELSE
            DO 155 I = 1, 3
  155       GTRAN(I,NL1) = GTRAN(I,NLEVEL)
            DO 156 I = 1, 10
  156       GRMAT(I,NL1) = GRMAT(I,NLEVEL)
         ENDIF
      ENDIF
*
  200 LINDEX(NL1) = IN
      LVOLUM(NL1) = IVOT
      NAMES(NL1)  = LNAM(NL1)
      NUMBER(NL1) = LNUM(NL1)
      LINMX(NL1)  = NINSK
      IF (LQ(LQ(JVOLUM-IVOT)).EQ.0) THEN
         NLDEV(NL1) = NLD
      ELSE
         NLDEV(NL1) = NL1
      ENDIF
      IQ(JGPAR+NL1) = NPAR
      LQ(JGPAR-NL1) = JPAR
      NLEVEL = NL1
      IF (NLEVEL.EQ.NLEVL) GO TO 990
      GO TO 100
*
* *** Error messages
*
  910 IER = 1
      WRITE (CHMAIL, 1000) NLEV
      CALL GMAIL (0, 0)
      GO TO 999
*
  920 IER = 2
      WRITE (CHMAIL, 2000) LNAM(1)
      CALL GMAIL (0, 0)
      GO TO 999
*
  930 IER = 3
      WRITE (CHMAIL, 3000) NLEVEL,NLEV,NAMES(NLEVEL)
      CALL GMAIL (0, 0)
      GO TO 999
*
  940 IER = 4
      WRITE (CHMAIL, 4000) LNAM(NL1),NL1,NAMES(NLEVEL)
      CALL GMAIL (0, 0)
      GO TO 999
*
  950 IER = 5
      WRITE (CHMAIL, 5000) NL1,LNUM(NL1),NAME,NDIV
      CALL GMAIL (0, 0)
      GO TO 999
*
  960 IER = 6
      WRITE (CHMAIL, 6000) NL1,LNAM(NL1),IQ(JVOLUM+IVOT)
      CALL GMAIL (0, 0)
      GO TO 999
*
  990 CONTINUE
#if defined(CERNLIB_DEBUG)
      WRITE (CHMAIL, 1001) NLEVEL
      CALL GMAIL (0, 0)
      DO 991 I1=1,NLEVEL,5
         I2 = I1 + 4
         IF (I2.GT.NLEVEL) I2 = NLEVEL
         WRITE (CHMAIL, 1003) (NAMES(I),NUMBER(I),LVOLUM(I),LINDEX(I),
     +   I=I1,I2)
         CALL GMAIL (0, 0)
  991 CONTINUE
      DO 992 I= 1,NLEVEL
         WRITE (CHMAIL, 1002) (GTRAN(J,I),J=1,3),(GRMAT(J,I),J=1,10)
         CALL GMAIL (0, 0)
  992 CONTINUE
 1001 FORMAT (' GLVOLU : NLEVEL =',I3)
 1003 FORMAT (5(1X,A4,3I3))
 1002 FORMAT (1X,13F9.4)
#endif
*
 1000 FORMAT (' GLVOLU : called with useless Level # ',I5)
 2000 FORMAT (' GLVOLU : Volume ',A4,' not top of tree, or no tree')
 3000 FORMAT (' GLVOLU : at Level ',I3,' of ',I3,' there are no',
     *        ' contents for Volume ',A4)
 4000 FORMAT (' GLVOLU : Volume ',A4,' for Level ',I3,
     *        ' does not exist in Volume ',A4)
 5000 FORMAT (' GLVOLU : at Level ',I3,' asked for #',I3,
     *        ' in divided Volume ',A4,' which has ',I3,' divisions.')
 6000 FORMAT (' GLVOLU : at Level ',I3,' user name ',A4,
     *        ' not equal to name ',A4,' of division.')
#if defined(CERNLIB_DEBUG)
 7000 FORMAT (' GLVOLU : Warning, ',A4,' not top of tree',
     *        ' you should reset NLEVEL to 0 before tracking !')
#endif
 8000 FORMAT (' GLVOLU : Volume ',A4,' Level 1 does not exist')
*                                                             END GLVOLU
  999 END
