*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GGDETV (ISET, IDET)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Routine - to compute and store the list of volumes which    *
C.    *    permit to identify uniquely any detector volume specified   *
C.    *    by  the set number ISET, the detector number IDET and the   *
C.    *    corresponding list of volume copy numbers                   *
C.    *            - to compute and store the physical path(s) through *
C.    *    the JVOLUM data structure down to the given detector volume *
C.    *                                                                *
C.    *    ==>Called by : GHCLOS                                       *
C.    *         Author  F.Bruyant  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
C.
      PARAMETER (NLVMAX=15,NSKMAX=20,NVMAX=20)
      INTEGER  IVOSK(NSKMAX,NLVMAX-1), LIMUL(NLVMAX), LINAM(NLVMAX)
     +,        LIST(2), NSK(NLVMAX-1)
      EQUIVALENCE (LINAM(1),WS(1)), (LIMUL(1),WS(NLVMAX+1)), (IVOSK(1,1)
     +,   WS(2*NLVMAX+1)), (NSK(1),WS((NSKMAX+2)*(NLVMAX-1)+3))
     +,   (LIST(1),WS((NSKMAX+3)*(NLVMAX-1)+3))
C.
C.         -------------------------------------------------------------
C.
      JS = LQ(JSET-ISET)
      JD = LQ(JS-IDET)
C
C     Check that JD bank has been created by GSDETV (not GSDET)
C      or has not been already processed.
C
      IF (IQ(JD+9).NE.-1) GO TO 999
      IQ(JD+9) = -2
      IHDET = IQ(JS+IDET)
      IF (IDEBUG.NE.0) THEN
         WRITE (CHMAIL, 1001) IHDET
         CALL GMAIL (0,0)
 1001    FORMAT (' GGDETV : Detector ',A4)
      ENDIF
C
C     Check that current detector is not an alias
C
      IALI = IQ(JD+10)
      IF (IALI.NE.0) GO TO 200
      NSOL = 0
      NV = 0
 
      NLIST = 0
      CALL VZERO (NSK, NLVMAX-1)
      NLEV = 1
      LINAM(1) = IHDET
      MULT1 = 1
   10 IVOS = IUCOMP (LINAM(NLEV), IQ(JVOLUM+1), NVOLUM)
C
C     Search for detector parents up to top of tree
C
   20 IF (IVOS.EQ.1) GO TO 60
C
      DO 40 IVO=1,NVOLUM
         IF (IVO.EQ.IVOS) GO TO 40
         JVO = LQ(JVOLUM-IVO)
         NIN = Q(JVO+3)
         IF (NIN.EQ.0) GO TO 40
         IF (NSOL.GT.0) THEN
C           Skip mother banks already found
            IF (IUCOMP (IVO, IVOSK(1,NLEV), NSK(NLEV)) .NE. 0) GO TO 40
         ENDIF
C
         IF (NIN.LT.0) THEN
C           Division case
            JDIV = LQ(JVO-1)
            IF (IFIX(Q(JDIV+2)).NE.IVOS) GO TO 40
            MULTI = ABS(Q(JDIV+3))
            IF (MULTI.EQ.0)  MULTI = 255
         ELSE
C           Position case
            MULTI = 0
            DO 30 IN=1,NIN
               JIN = LQ(JVO-IN)
               IF (IFIX(Q(JIN+2)).NE.IVOS) GO TO 30
               MULTI = MAX(MULTI, IFIX(Q(JIN+3)))
   30       CONTINUE
            IF (MULTI.EQ.0) GO TO 40
         ENDIF
C
C     New level found
C
         LIMUL(NLEV) = MULTI
         IF (NLEV.EQ.NLVMAX) GO TO 920
         IF (NSK(NLEV).EQ.NSKMAX) GO TO 930
         NSK(NLEV) = NSK(NLEV) +1
         IVOSK(NSK(NLEV),NLEV) = IVO
         NLEV = NLEV +1
         LINAM(NLEV) = IQ(JVOLUM+IVO)
         IVOS = IVO
         GO TO 20
C
   40 CONTINUE
C
C     No more path found at current level
C
      IF (NSK(NLEV).EQ.0) GO TO 910
      IF (NSK(NLEV).GT.1.OR.LIMUL(NLEV+1).GT.1) THEN
         DO 50 N = 1,NSK(NLEV)
            IVO = IVOSK(N,NLEV)
            NANEW = IQ(JVOLUM+IVO)
            IPJD = JD +10
            IF (NV.GT.0) THEN
               DO 49 I = 1,NV
                  IF (NANEW.EQ.IQ(IPJD+1)) GO TO 50
                  IPJD = IPJD +2
   49          CONTINUE
            ENDIF
            IF (NV.EQ.NVMAX) GO TO 940
            NV = NV +1
            IQ(IPJD+1) = NANEW
   50    CONTINUE
      ENDIF
      GO TO 90
C
C     Store current solution
C
   60 NSOL = NSOL +1
      LIMUL(NLEV) = 0
      IF (LIMUL(1).GT.MULT1) MULT1 = LIMUL(1)
#if defined(CERNLIB_DEBUGG)
      IF (IDEBUG.NE.0) THEN
         WRITE (CHMAIL, 1002) NSOL, NLEV
         CALL GMAIL (0,0)
         WRITE (CHMAIL, 1012) (LINAM(I),LIMUL(I),I=1,NLEV)
         CALL GMAIL (0,0)
 1002    FORMAT (' GGDETV DEBUG : NSOL,NLEV ',2I3)
 1012    FORMAT (15(1X,A4,I3))
      ENDIF
#endif
C
      DO 80 N = NLEV,1,-1
         LIST(NLIST+1) = LINAM(N)
         LIST(NLIST+2) = LIMUL(N)
         IF (N.EQ.NLEV)  LIST(NLIST+2) = NLEV
         NLIST = NLIST +2
   80 CONTINUE
      IF (NLEV.LT.3) GO TO 100
      NLEV = NLEV -1
C
   90 NSK(NLEV) = 0
      NLEV = NLEV -1
      IF (NLEV.GT.0) GO TO 10
C
  100 IF (MULT1.GT.1) THEN
        NV = NV +1
        IQ(JD+9+2*NV) = LINAM(1)
      ENDIF
C
C     Perform final operations on JD bank
C
      NW = 0
      IF (NV.EQ.0) GO TO 150
C
C     Compute maximum multiplicities
C
      DO 120 N = 1,NLIST,2
         IPJD = JD +10
         DO 110 I = 1,NV
            IF (IQ(IPJD+1).EQ.LIST(N))
     +          IQ(IPJD+2)=MAX(IQ(IPJD+2),LIST(N+1))
            IPJD = IPJD +2
  110    CONTINUE
  120 CONTINUE
C
      IF (IDEBUG.NE.0) THEN
         I2 = 0
  125    I1 = I2 + 1
         I2 = I1 + 14
         IF (I2.GT.NV) I2 = NV
         WRITE (CHMAIL, 1003) (IQ(JD+10+I),I=2*I1-1,2*I2)
         CALL GMAIL (0, 0)
         IF (I2.LT.NV)  GO TO 125
 1003    FORMAT (10X,15(1X,A4,I3))
      ENDIF
C
C     Compute corresponding bit numbers for packing
C
      IPJD = JD +10
      K = 32
      DO 140 N = 1,NV
         NBITS = 0
  130    NBITS = NBITS +1
         IF (IQ(IPJD+2).GT.2**NBITS-1) GO TO 130
         IF (NBITS.GE.32)  NBITS = 0
         IQ(IPJD+2) = NBITS
         IPJD = IPJD +2
         IF (NBITS.EQ.0) THEN
            K = 32
            NW = NW +1
         ELSE
            K = K +NBITS
            IF (K.LE.32) GO TO 140
            K = NBITS
            NW = NW +1
         ENDIF
  140 CONTINUE
C
  150 IQ(JD+1) = NW
      IQ(JD+2) = NV
      IQ(JD+9) = NSOL
C
      NDATA = 10 +2*NV +NLIST
      ND = IQ(JD-1)
      CALL MZPUSH (IXCONS, JD, 0, NDATA-ND, 'I')
      CALL UCOPY (LIST, IQ(JD+2*NV+11), NLIST)
#if defined(CERNLIB_DEBUGG)
      IF (IDEBUG.NE.0) THEN
         ND1=MIN(10,NDATA)
         WRITE (CHMAIL, 1004) NDATA,(IQ(JD+I),I=1,ND1)
         CALL GMAIL (0,0)
         DO 160 II=ND1+1,NDATA,10
            ND2=MIN(II+9,NDATA)
            WRITE (CHMAIL, 1005) (IQ(JD+I),I=II,ND2)
            CALL GMAIL (0,0)
 160     CONTINUE
 1004    FORMAT (' GGDETV DEBUG : NDATA ',I3,'  JD --> ',10I4)
 1005    FORMAT (10(1X,A4,I4))
      ENDIF
#endif
      GO TO 999
C
C     Current detector IDET is an alias
C
  200 CONTINUE
      IF (IDEBUG.NE.0) THEN
         IHALI = IQ(JS+IALI)
         WRITE (CHMAIL, 1006) IHALI
         CALL GMAIL (0,0)
 1006    FORMAT ('       Alias of detector ',A4)
      ENDIF
C
      IDM = IQ(JD+10)
      JDM = LQ(JS-IDM)
      NDM = IQ(JDM-1)
      ND = IQ(JD-1)
      CALL MZPUSH (IXCONS, JD, 0, NDM-ND, 'I')
      NWHI = IQ(JD+7)
      NWDI = IQ(JD+8)
      JS = LQ(JSET-ISET)
      JDM  = LQ(JS-IDM)
      CALL UCOPY (IQ(JDM+1), IQ(JD+1), NDM)
      IQ(JD+7) = NWHI
      IQ(JD+8) = NWDI
      IQ(JD+10) = IDM
      GO TO 999
C
C     Errors
C
  910 WRITE (CHMAIL, 1000) LINAM(NLEV)
      CALL GMAIL (0,0)
 1000 FORMAT (' GGDETV : Hanging volume ',A4)
      GO TO 990
  920 CHMAIL=' GGDETV : Parameter NLVMAX too small'
      CALL GMAIL (0,0)
      GO TO 990
  930 CHMAIL=' GGDETV : Parameter NSKMAX too small'
      CALL GMAIL (0,0)
      GO TO 990
  940 CHMAIL=' GGDETV : NVMAX (= size of NUMBV) too small'
      CALL GMAIL (0,0)
  990 IEOTRI = 1
C
  999 RETURN
      END
