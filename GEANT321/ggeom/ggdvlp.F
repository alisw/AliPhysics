*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GGDVLP
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GGDVLP                                                *
C.    *                                                                *
C.    *   DeVeLoPs locally the JVOLUM structure so as not to have to   *
C.    *   recompute variable parameters and division specifications    *
C.    *   at tracking time.                                            *
C.    *                                                                *
C.    *   Called by : GGCLOS                                           *
C.    *   Authors   : S.Banerjee, F.Bruyant                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
C.
      PARAMETER (NLVMAX=15, NPAMAX=50)
C.
      COMMON /GCDVLP/ LREF(2), JVO, LJVOM(NLVMAX), LLVODV(NLVMAX), LVODV
C.
      INTEGER     ILINK(NLVMAX), LVAR(NPAMAX), NLINK(NLVMAX)
      INTEGER     NUMB(NLVMAX+1), NBGN(NLVMAX+1), IVOL(5200)
      EQUIVALENCE (IVOL(1), WS(1))
      REAL        PAR(NPAMAX)
      LOGICAL     BTEST
C.
C.    ------------------------------------------------------------------
*
      CALL MZLINT (IXSTOR, '/GCDVLP/', LREF, JVO, LVODV)
*
      CALL MZFORM ('VODV', '1I 2F 2I -F', IOVODV)
*
* *** Create the volume tree in memory
*
      IF (NVOLUM.LT.1) GO TO 990
      NUMB(1) = 1
      NBGN(1) = 1
      IVOL(1) = 1
      NLVT = 1
*
   10 NLEV = NLVT + 1
      NUMB(NLEV) = 0
      NBGN(NLEV) = NBGN(NLVT) + NUMB(NLVT)
      DO 30 I = 1, NUMB(NLVT)
         IVO  = IVOL(NBGN(NLVT) + I - 1)
         JVO  = LQ(JVOLUM-IVO)
         NIN  = Q(JVO+3)
         IF (NIN.GT.0) THEN
            DO 20 IN = 1, NIN
               JIN  = LQ(JVO-IN)
               IVOT = Q(JIN+2)
               DO 15 I1 = 1, NUMB(NLEV)
                  IVOS = IVOL(NBGN(NLEV)+I1-1)
                  IF (IVOS.EQ.IVOT) GO TO 20
   15          CONTINUE
               NUMB(NLEV) = NUMB(NLEV) + 1
               I1         = NBGN(NLEV) + NUMB(NLEV) - 1
               IVOL(I1)   = IVOT
   20       CONTINUE
         ELSE IF (NIN.LT.0) THEN
            JDIV = LQ(JVO-1)
            IVOT = Q(JDIV+2)
            DO 25 I1 = 1, NUMB(NLEV)
               IVOS = IVOL(NBGN(NLEV)+I1-1)
               IF (IVOS.EQ.IVOT) GO TO 30
   25       CONTINUE
            NUMB(NLEV) = NUMB(NLEV) + 1
            I1         = NBGN(NLEV) + NUMB(NLEV) - 1
            IVOL(I1)   = IVOT
         ENDIF
   30 CONTINUE
*
      IF (NUMB(NLEV).GT.0) THEN
         NLVT = NLEV
         GO TO 10
      ENDIF
*
      NLEV  = 1
      NLEVS = 1
*
* *** Loop over volumes in a given level
*
  110 IF (NLEVS.GT.NLVT) GO TO 990
      IF (NUMB(NLEVS).GT.0) THEN
         NUMB(NLEVS) = NUMB(NLEVS) - 1
         IVOM  = IVOL(NBGN(NLEVS)+NUMB(NLEVS))
      ELSE
         NLEVS = NLEVS + 1
         GO TO 110
      ENDIF
*
* *** Check if current volume should be locally developed
*
      LJVOM(1) = LQ(JVOLUM-IVOM)
      IF(LQ(LJVOM(1)).GT.0) THEN
         CALL MZDROP(IXCONS,LQ(LJVOM(1)),'L')
      ENDIF
      NIN  = Q(LJVOM(1)+3)
#if defined(CERNLIB_DEBUGG)
      WRITE (CHMAIL,2000) IQ(JVOLUM+IVOM),NIN
      CALL GMAIL (0, 0)
 2000 FORMAT (' GGDVLP : Volume ',A4,'  NIN = ',I3)
#endif
*
*  **  Skip it if not a possible source of local development
*
      IF (NIN.EQ.0)                  GO TO 110
      IF (LQ(LJVOM(1)).NE.0)         GO TO 110
      IF (BTEST(IQ(LJVOM(1)),1))     GO TO 110
      CALL GGVCHK (LJVOM(1), 0, NVAR, LVAR)
      IF (NVAR.NE.0) GO TO 110
*
*  **  Otherwise, analyze contents
*
      IOK = 0
      IF (NIN.LT.0) GO TO 120
*
*   *   Current volume has contents defined by position
*
      DO 119 IN = 1,NIN
         CALL GGVCHK (LJVOM(1), IN, NVAR, LVAR)
         IF (NVAR.NE.0) THEN
*
*         In case a content is found with variable parameters,
*          initialize development
*
            NLINK(1) = NIN
            ILINK(1) = IN
            GO TO 200
         ENDIF
  119 CONTINUE
      GO TO 110
*
*   *   Current volume is divided
*
  120 CALL GGVCHK (LJVOM(1), 1, NVAR, LVAR)
      IF (NVAR.EQ.0) GO TO 110
*
*        If cells have variable sizes, initialize development
*
      IOK  = 2
      JDIV = LQ(LJVOM(1)-1)
      NLINK(1) = Q(JDIV+3)
      ILINK(1) = 1
*
  200 CALL MZBOOK (IXCONS, LLVODV(1), LJVOM(1), 0, 'VODV',
     +                         NLINK(1), NLINK(1), 1, 2, 1)
#if defined(CERNLIB_DEBUGG)
      WRITE (CHMAIL, 2002) ILINK(1), NLINK(1), NVAR
      CALL GMAIL (0, 0)
 2002 FORMAT (' GGDVLP :         I, N, NVAR = ',3I5)
#endif
*
* *** Complete development for current content at current level
*
      IF (IOK.EQ.0) GO TO 250
      IOK = 0
      IQ(LLVODV(1)+1) = NLINK(1)
      GO TO 260
*
  210 NIN = Q(LJVOM(NLEV)+3)
      IF (NIN.LT.0) THEN
         CALL GGVCHK (LJVOM(NLEV), 1, NVAR, LVAR)
      ELSE
         CALL GGVCHK (LJVOM(NLEV), ILINK(NLEV), NVAR, LVAR)
      ENDIF
      IF (NVAR.EQ.0) GO TO 290
#if defined(CERNLIB_DEBUGG)
      WRITE (CHMAIL, 2004) NLEV, ILINK(NLEV), NVAR
      CALL GMAIL (0, 0)
 2004 FORMAT (' GGDVLP :              LEVEL,I,NVAR,NIN = ',4I5)
#endif
      IF (NIN.LT.0)  GO TO 260
*
*  **  Compute actual parameters for current content
*
*   *   Case with contents obtained by position
*
  250 CALL GGPPAR (LJVOM(NLEV), ILINK(NLEV), NVAR,LVAR, LLVODV(NLEV),
     +                                                     NPAR, PAR)
      JIN = LQ(LJVOM(NLEV)-ILINK(NLEV))
      GO TO 270
*
*   *   Case with contents obtained by division
*
  260 CALL GGDPAR (LJVOM(NLEV), ILINK(NLEV), NVAR,LVAR, LLVODV(NLEV),
     +                                                     NPAR, PAR)
      JIN = LQ(LJVOM(NLEV)-1)
*
*  **  Initialize next level down
*
  270 IVO = Q(JIN+2)
      JVO = LQ(JVOLUM-IVO)
      NIN = Q(JVO+3)
      IF (NIN.GE.0) THEN
         IF (NIN.NE.0) IOK = 1
         NL = NIN
      ELSE
         IOK = 2
*
*   *   Current content is divided, compute division specifications
*
         CALL GGDSPE (JVO, NPAR, PAR, NL, NDIV, ORIG, STEP)
      ENDIF
*
      IQ(JVO) = IBSET(IQ(JVO),1)
      CALL MZBOOK (IXCONS, LVODV, LLVODV(NLEV), -ILINK(NLEV), 'VODV',
     +                  NL, NL, NPAR+5, IOVODV, 3)
      IF (IOK.EQ.2) THEN
         IQ(LVODV+1) = NDIV
          Q(LVODV+2) = ORIG
          Q(LVODV+3) = STEP
#if defined(CERNLIB_DEBUGG)
          WRITE (CHMAIL, 2006) NDIV, ORIG, STEP
          CALL GMAIL (0, 0)
 2006     FORMAT (' GGDVLP : After GGDSPE, NDIV ORIG STEP = ',I4,2F10.4)
#endif
      ENDIF
      IQ(LVODV+4) = IVO
      IQ(LVODV+5) = NPAR
      CALL UCOPY (PAR, Q(LVODV+6), NPAR)
#if defined(CERNLIB_DEBUGG)
      WRITE (CHMAIL, 2008)
      CALL GMAIL (0, 0)
 2008 FORMAT (' GGDVLP : Store parameters into development structure')
#endif
*
  290 IF (IOK.EQ.0) THEN
         IF (ILINK(NLEV).EQ.NLINK(NLEV)) THEN
*
*          Go one level up
*
            IF (NLEV.EQ.1) GO TO 110
            NLEV = NLEV -1
            GO TO 290
         ENDIF
*
*       Analyze next content
*
         ILINK(NLEV) = ILINK(NLEV) +1
         GO TO 210
      ENDIF
*
*      A new level has been initialized, start analyzing it
*
      NLEV = NLEV +1
      LJVOM(NLEV)  = JVO
      LLVODV(NLEV) = LVODV
      ILINK(NLEV)  = 1
      NLINK(NLEV)  = NL
      IOK = 0
      GO TO 210
*
  990 DO 991 IVO = 1, NVOLUM
         JVO = LQ(JVOLUM-IVO)
         IF (BTEST(IQ(JVO),1)) THEN
            IF (LQ(JVO).GT.0) THEN
               CALL MZDROP (IXCONS, LQ(JVO), ' ')
               WRITE (CHMAIL, 1001) IQ(JVOLUM+IVO)
            ENDIF
         ENDIF
  991 CONTINUE
      LREF(1) = 0
*
 1001 FORMAT (' GGDVLP : Unnecessary development at volume  ',A4)
*                                                             END GGDVLP
      END
