*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_CG)
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCGWV(IWOFFS)
************************************************************************
*                                                                      *
*     Name: GDCGWV    extracted from CG routine CGVWIR                 *
*     Author: E. Chernyaev                       Date:    09.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Vizualisation of WIRE-object                           *
*                                                                      *
*     References: CGVHED                                               *
*                                                                      *
*     Input: Q(JCG+IWOFFS+*) - WIRE-object                             *
*               Q(JCG+*) - HIDE-structure                              *
*                                                                      *
*     Output: none                                                     *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/cggpar.inc"
#include "geant321/cghpar.inc"
#include "geant321/cgctra.inc"
#include "geant321/cgcedg.inc"
#include "geant321/gcdraw.inc"
*
       REAL P1(3),P2(3),AX(2),AY(2)
*      SAVE FIRST,JPOLD
*      LOGICAL FIRST
*      DATA FIRST/.TRUE./
*SG
*      DIMENSION TESTX(2),TESTY(2)
*SG
*      IF(FIRST) THEN
*         JPOLD=LQ(LQ(JDRAW-3)-5)
*         FIRST=.FALSE.
*      ENDIF
      IF (Q(JCG+IWOFFS+KCGSIZ) .LE. 0.)   THEN
         WRITE(CHMAIL,10000)
         GOTO 50
      ENDIF
      IF (Q(JCG+KHSIZE) .LE. LHHEAD)  THEN
         WRITE(CHMAIL,10100)
         GOTO 50
      ENDIF
*SG
*  Icolar is equal to 1 even if hide is off, but linbuf can be = 2 only
*  if fill area option is switched on AND hide is on.
*
*      LINBUF=0
*      DO 111 IJ=1,2
*        TESTX(IJ)=0
*        TESTY(IJ)=0
* 111  CONTINUE
*      LINFIL=IBITS(LINATT,17,3)
*      IF(LINFIL.GT.0)THEN
*        CALL GDRAWV(TESTX,TESTY,-1)
*      ENDIF
*SG
      NT = Q(JCG+KHNT)
      IF (NT.LE.0 .OR. NT.GT.NTMAX)     GOTO 50
      NEDGE  = Q(JCG+IWOFFS+KCGNF)
*      IF(NEDGE.EQ.0)GOTO 49
      JXYZ1   = Q(JCG+KHJXYZ)
      JFA    = Q(JCG+KHJFA)
      JPFA   = Q(JCG+KHJPFA)
      JDFA   = Q(JCG+KHJDFA)
      JTRE   = Q(JCG+KHJTRE)
      JSTA   = Q(JCG+KHJSTA)
      JALE   = Q(JCG+KHJALE)
      JARI   = Q(JCG+KHJARI)
      NFACE  = Q(JCG+KHNFAC)
      DO 40  NE=1,NEDGE
         J = LCGHEA + (NE-1)*LCGEDG
*         IEA = Q(JCG+IWOFFS+J+KCGAE)
         P1(1) = Q(JCG+IWOFFS+J+KCGX1)
         P1(2) = Q(JCG+IWOFFS+J+KCGY1)
         P1(3) = Q(JCG+IWOFFS+J+KCGZ1)
         P2(1) = Q(JCG+IWOFFS+J+KCGX2)
         P2(2) = Q(JCG+IWOFFS+J+KCGY2)
         P2(3) = Q(JCG+IWOFFS+J+KCGZ2)
         CALL CGVEDG(NT,P1,P2,IVIS)
         IF (IVIS .LT. 0) GOTO 40
         IF (NFACE .EQ. 0) GOTO 30
         CALL CGVHED( Q(JCG+JXYZ1),Q(JCG+JFA),Q(JCG+JPFA),
     +   Q(JCG+JDFA),Q(JCG+JTRE),Q(JCG+JSTA),
     +   Q(JCG+JALE),Q(JCG+JARI))
*           D R A W   E D G E
         IF (NPART) 40 ,30 ,10
   10    DO 20  I=1,NPART+1
            IF (I .EQ. 1) T1 = 0.
            IF (I .NE. 1) T1 = TEND(I-1)
            IF (I .NE. NPART+1) T2 = TSTRT(I)
            IF (I .EQ. NPART+1) T2 = 1.
            IF (T2-T1 .LT. TDEL) GOTO 20
            AX(1) = XA + T1*XDELT
            AY(1) = YA + T1*YDELT
            AX(2) = XA + T2*XDELT
            AY(2) = YA + T2*YDELT
*            JPNOW=LQ(LQ(JDRAW-3)-5)
*            IF(JPNOW.NE.JPOLD) THEN
*               PRINT*,'JPOLD, JPNEW before 1st GDRAWV ',JPOLD,JPNEW
*               Q(JCG+100000000)=1.
*            ENDIF
            CALL GDRAWV(AX,AY,2)
*            IF(JPNOW.NE.JPOLD) THEN
*               PRINT*,'JPOLD, JPNEW after 1st GDRAWV ',JPOLD,JPNEW
*               Q(JCG+100000000)=1.
*            ENDIF
   20    CONTINUE
         GOTO 40
*           L I N E   I S   F U L L Y   V I S I B L E
   30    CONTINUE
         AX(1) = AA(1)
         AY(1) = AA(2)
         AX(2) = BB(1)
         AY(2) = BB(2)
*            IF(JPNOW.NE.JPOLD) THEN
*               PRINT*,'JPOLD, JPNEW before 2nd GDRAWV ',JPOLD,JPNEW
*               Q(JCG+100000000)=1.
*            ENDIF
         CALL GDRAWV(AX,AY,2)
*            IF(JPNOW.NE.JPOLD) THEN
*               PRINT*,'JPOLD, JPNEW after 2nd GDRAWV ',JPOLD,JPNEW
*               Q(JCG+100000000)=1.
*            ENDIF
   40 CONTINUE
*SG
*      IF(LINFIL.EQ.1)CALL GDRAWV(TESTX,TESTY,0)
*   49 CONTINUE
*SG
   50 RETURN
*
10000    FORMAT('  First word of WIRE less or equal 0 ')
10100    FORMAT(' N. of words of Hidden Structure less or equal 18 ')
*
      END
#endif
