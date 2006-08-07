* $Log$
* Revision 1.3  2005/12/02 14:50:54  whalley
* Changes for new CTEQ code/AB sets
*
* Revision 1.2  2005/10/18 10:34:52  whalley
* small changes from cdf/d0 comments - renamubg conflicting name
* adding pftopdg, etc.
*
* Revision 1.1.1.1  2005/05/06 14:54:44  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
* Mathlib gen
*
*
      FUNCTION DDILOG_LHA(X)
      implicit real*8 (a-h,o-z)
      DIMENSION C(0:19)
      PARAMETER (Z1 = 1d0, HF = Z1/2d0)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1d0) THEN
       H=PI6
      ELSEIF(X .EQ. -1d0) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2d0) THEN
        Y=-1/(1d0+T)
        S=1d0
        A=-PI3+HF*(LOG(-T)**2-LOG(1d0+1d0/T)**2)
       ELSEIF(T .LT. -1d0) THEN
        Y=-1d0-T
        S=-1d0
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1d0+1d0/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1d0+T)/T
        S=1d0
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1d0+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1d0+T)
        S=-1d0
        A=HF*LOG(1d0+T)**2
       ELSE IF(T .LE. 1d0) THEN
        Y=T
        S=1d0
        A=0d0
       ELSE
        Y=1d0/T
        S=-1d0
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DDILOG_LHA=H
      RETURN
      END
*
* $Id$
*
* $Log$
* Revision 1.3  2005/12/02 14:50:54  whalley
* Changes for new CTEQ code/AB sets
*
* Revision 1.2  2005/10/18 10:34:52  whalley
* small changes from cdf/d0 comments - renamubg conflicting name
* adding pftopdg, etc.
*
* Revision 1.1.1.1  2005/05/06 14:54:44  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
      FUNCTION DGAUSS_LHA(F,A,B,EPS)

      implicit real*8 (a-h,o-z)
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DGAUSS_LHA')
      DIMENSION W(12),X(12)
      PARAMETER (Z1 = 1d0, HF = Z1/2d0, CST = 5d0*Z1/1000d0)
      DATA X
     1        /0.96028 98564 97536 23168 35608 68569 47D0,
     2         0.79666 64774 13626 73959 15539 36475 83D0,
     3         0.52553 24099 16328 98581 77390 49189 25D0,
     4         0.18343 46424 95649 80493 94761 42360 18D0,
     5         0.98940 09349 91649 93259 61541 73450 33D0,
     6         0.94457 50230 73232 57607 79884 15534 61D0,
     7         0.86563 12023 87831 74388 04678 97712 39D0,
     8         0.75540 44083 55003 03389 51011 94847 44D0,
     9         0.61787 62444 02643 74844 66717 64048 79D0,
     A         0.45801 67776 57227 38634 24194 42983 58D0,
     B         0.28160 35507 79258 91323 04605 01460 50D0,
     C         0.95012 50983 76374 40185 31933 54249 58D-1/

      DATA W
     1        /0.10122 85362 90376 25915 25313 54309 96D0,
     2         0.22238 10344 53374 47054 43559 94426 24D0,
     3         0.31370 66458 77887 28733 79622 01986 60D0,
     4         0.36268 37833 78361 98296 51504 49277 20D0,
     5         0.27152 45941 17540 94851 78057 24560 18D-1,
     6         0.62253 52393 86478 92862 84383 69943 78D-1,
     7         0.95158 51168 24927 84809 92510 76022 46D-1,
     8         0.12462 89712 55533 87205 24762 82192 02D0,
     9         0.14959 59888 16576 73208 15017 30547 48D0,
     A         0.16915 65193 95002 53818 93120 79030 36D0,
     B         0.18260 34150 44923 58886 67636 67969 22D0,
     C         0.18945 06104 55068 49628 53967 23208 28D0/

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
       H=0
       write(*,*) NAME,'D103.1','TOO HIGH ACCURACY REQUIRED'
       GO TO 99
      END IF
   99 DGAUSS_LHA=H
      RETURN
      END

      SUBROUTINE VZERO (A,N)
C
C CERN PROGLIB# F121    VZERO           .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71, modif. 24/05/87 to set integer zero
C                 modif. 25/05/94 to depend on QINTZERO
C
      DIMENSION A(*)
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= 0d0
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.3  2005/12/02 14:50:54  whalley
* Changes for new CTEQ code/AB sets
*
* Revision 1.2  2005/10/18 10:34:52  whalley
* small changes from cdf/d0 comments - renamubg conflicting name
* adding pftopdg, etc.
*
* Revision 1.1.1.1  2005/05/06 14:54:44  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
      FUNCTION LENOCC_LHA (CHV)
C
C CERN PROGLIB# M507    LENOCC          .VERSION KERNFOR  4.21  890323
C ORIG. March 85, A.Petrilli, re-write 21/02/89, JZ
C
C-    Find last non-blank character in CHV

      CHARACTER    CHV*(*)

      N = LEN(CHV)

      DO 17  JJ= N,1,-1
      IF (CHV(JJ:JJ).NE.' ') GO TO 99
   17 CONTINUE
      JJ = 0

   99 LENOCC_LHA = JJ
      RETURN
      END
*
* $Id$
*
* $Log$
* Revision 1.3  2005/12/02 14:50:54  whalley
* Changes for new CTEQ code/AB sets
*
* Revision 1.2  2005/10/18 10:34:52  whalley
* small changes from cdf/d0 comments - renamubg conflicting name
* adding pftopdg, etc.
*
* Revision 1.1.1.1  2005/05/06 14:54:44  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
      SUBROUTINE CLTOU_LHA (CHV)
C
C CERN PROGLIB# M432    CLTOU           .VERSION KERNFOR  4.21  890323
C ORIG. 11/02/86 A. PETRILLI
C NEW    9/02/89 JZ, for speed
C
C-    Convert character string CHV from lower to upper case.

      CHARACTER    CHV*(*)
      DO 19  JJ=1,LEN(CHV)
          J = ICHAR(CHV(JJ:JJ))
          IF (J.LT.97)       GO TO 19
          IF (J.GE.123)      GO TO 19
          CHV(JJ:JJ) = CHAR(J-32)
   19 CONTINUE
      END
*
      SUBROUTINE TIMEX_LHA (T)
C
C CERN PROGLIB# Z007    TIMEX   DUMMY   .VERSION KERNFOR  4.05  821202
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      T = 9.
      RETURN
      END
*
* $Id$
*
* $Log$
* Revision 1.3  2005/12/02 14:50:54  whalley
* Changes for new CTEQ code/AB sets
*
* Revision 1.2  2005/10/18 10:34:52  whalley
* small changes from cdf/d0 comments - renamubg conflicting name
* adding pftopdg, etc.
*
* Revision 1.1.1.1  2005/05/06 14:54:44  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
      SUBROUTINE FLPSOR_LHA(A,N)
C
C CERN PROGLIB# M103    FLPSOR          .VERSION KERNFOR  3.15  820113
C ORIG. 29/04/78
C
C   SORT THE ONE-DIMENSIONAL FLOATING POINT ARRAY A(1),...,A(N) BY
C   INCREASING VALUES
C
C-    PROGRAM  M103  TAKEN FROM CERN PROGRAM LIBRARY,  29-APR-78
C
      DIMENSION A(N)
      COMMON /SLATE/ LT(20),RT(20)
      INTEGER R,RT
C
      LEVEL=1
      LT(1)=1
      RT(1)=N
   10 L=LT(LEVEL)
      R=RT(LEVEL)
      LEVEL=LEVEL-1
   20 IF(R.GT.L) GO TO 200
      IF(LEVEL) 50,50,10
C
C   SUBDIVIDE THE INTERVAL L,R
C     L : LOWER LIMIT OF THE INTERVAL (INPUT)
C     R : UPPER LIMIT OF THE INTERVAL (INPUT)
C     J : UPPER LIMIT OF LOWER SUB-INTERVAL (OUTPUT)
C     I : LOWER LIMIT OF UPPER SUB-INTERVAL (OUTPUT)
C
  200 I=L
      J=R
      M=(L+R)/2
      X=A(M)
  220 IF(A(I).GE.X) GO TO 230
      I=I+1
      GO TO 220
  230 IF(A(J).LE.X) GO TO 231
      J=J-1
      GO TO 230
C
  231 IF(I.GT.J) GO TO 232
      W=A(I)
      A(I)=A(J)
      A(J)=W
      I=I+1
      J=J-1
      IF(I.LE.J) GO TO 220
C
  232 LEVEL=LEVEL+1
      IF((R-I).GE.(J-L)) GO TO 30
      LT(LEVEL)=L
      RT(LEVEL)=J
      L=I
      GO TO 20
   30 LT(LEVEL)=I
      RT(LEVEL)=R
      R=J
      GO TO 20
   50 RETURN
      END

      SUBROUTINE DATIMH_LHA (ND,NT)
C
C CERN PROGLIB# Z007    DATIMH  DUMMY   .VERSION KERNFOR  4.03  821008
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      DIMENSION ND(9), NT(9)
*      DIMENSION M(8)

      do i=1,9
         ND(i)=0
         NT(i)=0
      enddo
*      CALL UBLOW (8H29/09/79,M,8)
*      CALL UBUNCH           (M,ND,8)
*      CALL UBLOW (8H12.00.00,M,8)
*      CALL UBUNCH           (M,NT,8)
      RETURN
      END
c**************************************************************************
c these next added from CERNLIB to allow some photon and pion sets to work
c nothing to do with QCDNUM mrw 9/12/2004
c*************************************************************************8
c*
* $Id$
*
* $Log$
* Revision 1.3  2005/12/02 14:50:54  whalley
* Changes for new CTEQ code/AB sets
*
* Revision 1.2  2005/10/18 10:34:52  whalley
* small changes from cdf/d0 comments - renamubg conflicting name
* adding pftopdg, etc.
*
* Revision 1.1.1.1  2005/05/06 14:54:44  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.1  1996/02/15 17:48:17  mclareni
* Kernlib
*
*
      DOUBLE PRECISION FUNCTION DGAMMA_LHA(X)
      LOGICAL MFLAG,RFLAG
      REAL SX
      DOUBLE PRECISION X,U,F,ZERO,ONE,THREE,FOUR,PI
      DOUBLE PRECISION C(0:24),H,ALFA,B0,B1,B2
      DATA ZERO /0.0D0/, ONE /1.0D0/, THREE /3.0D0/, FOUR /4.0D0/
c#if defined(CERNLIB_NUMHIPRE)
c      DATA NC /24/
c      DATA PI    /3.14159 26535 89793 23846 26433 83D0/
c      DATA C( 0) /3.65738 77250 83382 43849 88068 39D0/
c      DATA C( 1) /1.95754 34566 61268 26928 33742 26D0/
c      DATA C( 2) / .33829 71138 26160 38915 58510 73D0/
c      DATA C( 3) / .04208 95127 65575 49198 51083 97D0/
c      DATA C( 4) / .00428 76504 82129 08770 04289 08D0/
c      DATA C( 5) / .00036 52121 69294 61767 02198 22D0/
c      DATA C( 6) / .00002 74006 42226 42200 27170 66D0/
c      DATA C( 7) / .00000 18124 02333 65124 44603 05D0/
c      DATA C( 8) / .00000 01096 57758 65997 06993 06D0/
c      DATA C( 9) / .00000 00059 87184 04552 00046 95D0/
c      DATA C(10) / .00000 00003 07690 80535 24777 71D0/
c      DATA C(11) / .00000 00000 14317 93029 61915 76D0/
c      DATA C(12) / .00000 00000 00651 08773 34803 70D0/
c      DATA C(13) / .00000 00000 00025 95849 89822 28D0/
c      DATA C(14) / .00000 00000 00001 10789 38922 59D0/
c      DATA C(15) / .00000 00000 00000 03547 43620 17D0/
c      DATA C(16) / .00000 00000 00000 00168 86075 04D0/
c      DATA C(17) / .00000 00000 00000 00002 73543 58D0/
c      DATA C(18) / .00000 00000 00000 00000 30297 74D0/
c      DATA C(19) /-.00000 00000 00000 00000 00571 22D0/
c      DATA C(20) / .00000 00000 00000 00000 00090 77D0/
c      DATA C(21) /-.00000 00000 00000 00000 00005 05D0/
c      DATA C(22) / .00000 00000 00000 00000 00000 41D0/
c      DATA C(23) /-.00000 00000 00000 00000 00000 03D0/
c      DATA C(24) / .00000 00000 00000 00000 00000 01D0/
c#endif
c#if defined(CERNLIB_NUMLOPRE)
      DATA NC /15/
      DATA PI    /3.14159 26535 89793 24D0/
      DATA C( 0) /3.65738 77250 83382 44D0/
      DATA C( 1) /1.95754 34566 61268 27D0/
      DATA C( 2) / .33829 71138 26160 39D0/
      DATA C( 3) / .04208 95127 65575 49D0/
      DATA C( 4) / .00428 76504 82129 09D0/
      DATA C( 5) / .00036 52121 69294 62D0/
      DATA C( 6) / .00002 74006 42226 42D0/
      DATA C( 7) / .00000 18124 02333 65D0/
      DATA C( 8) / .00000 01096 57758 66D0/
      DATA C( 9) / .00000 00059 87184 05D0/
      DATA C(10) / .00000 00003 07690 81D0/
      DATA C(11) / .00000 00000 14317 93D0/
      DATA C(12) / .00000 00000 00651 09D0/
      DATA C(13) / .00000 00000 00025 96D0/
      DATA C(14) / .00000 00000 00001 11D0/
      DATA C(15) / .00000 00000 00000 04D0/
c#endif
      U=X
      IF(X .LE. ZERO) THEN
       IF(X .EQ. INT(X)) THEN
        CALL KERMTR_LHA('C305.1',LGFILE,MFLAG,RFLAG)
        IF(MFLAG) THEN
         SX=X
         IF(LGFILE .EQ. 0) THEN
          WRITE(*,100) SX
         ELSE
          WRITE(LGFILE,100) SX
         END IF
        END IF
        IF(.NOT.RFLAG) CALL ABEND_LHA
        DGAMMA_LHA=ZERO
        RETURN
       ELSE
        U=ONE-U
       END IF
      END IF
      F=ONE
      IF(U .LT. THREE) THEN
       DO 1 I = 1,INT(FOUR-U)
       F=F/U
    1  U=U+ONE
      ELSE
       DO 2 I = 1,INT(U-THREE)
       U=U-ONE
    2  F=F*U
      END IF
      U=U-THREE
      H=U+U-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 3 I = NC,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    3 B1=B0
      U=F*(B0-H*B2)
      IF(X .LT. ZERO) U=PI/(SIN(PI*X)*U)
      DGAMMA_LHA=U
      RETURN
  100 FORMAT(1X,'DGAMMA ... ARGUMENT IS NON-POSITIVE INTEGER = ',E15.1)
      END
c========================================================================
	  SUBROUTINE KERSET_LHA(ERCODE,LGFILE,LIMITM,LIMITR)
	  PARAMETER(KOUNTE  =  28)
	  CHARACTER*6	      ERCODE,	CODE(KOUNTE)
	  LOGICAL	      MFLAG,	RFLAG
	  INTEGER	      KNTM(KOUNTE),	  KNTR(KOUNTE)
	  DATA      LOGF      /  0  /
	  DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 100, 100 /
	  DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 100, 100 /
	  DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 100, 100 /
	  DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 100, 100 /
	  DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 100, 100 /
	  DATA      CODE(6), KNTM(6), KNTR(6)  / 'C205.3', 100, 100 /
	  DATA      CODE(7), KNTM(7), KNTR(7)  / 'C305.1', 100, 100 /
	  DATA      CODE(8), KNTM(8), KNTR(8)  / 'C308.1', 100, 100 /
	  DATA      CODE(9), KNTM(9), KNTR(9)  / 'C312.1', 100, 100 /
	  DATA      CODE(10),KNTM(10),KNTR(10) / 'C313.1', 100, 100 /
	  DATA      CODE(11),KNTM(11),KNTR(11) / 'C336.1', 100, 100 /
	  DATA      CODE(12),KNTM(12),KNTR(12) / 'C337.1', 100, 100 /
	  DATA      CODE(13),KNTM(13),KNTR(13) / 'C341.1', 100, 100 /
	  DATA      CODE(14),KNTM(14),KNTR(14) / 'D103.1', 100, 100 /
	  DATA      CODE(15),KNTM(15),KNTR(15) / 'D106.1', 100, 100 /
	  DATA      CODE(16),KNTM(16),KNTR(16) / 'D209.1', 100, 100 /
	  DATA      CODE(17),KNTM(17),KNTR(17) / 'D509.1', 100, 100 /
	  DATA      CODE(18),KNTM(18),KNTR(18) / 'E100.1', 100, 100 /
	  DATA      CODE(19),KNTM(19),KNTR(19) / 'E104.1', 100, 100 /
	  DATA      CODE(20),KNTM(20),KNTR(20) / 'E105.1', 100, 100 /
	  DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.1', 100, 100 /
	  DATA      CODE(22),KNTM(22),KNTR(22) / 'E208.2', 100, 100 /
	  DATA      CODE(23),KNTM(23),KNTR(23) / 'F010.1', 100,   0 /
	  DATA      CODE(24),KNTM(24),KNTR(24) / 'F011.1', 100,   0 /
	  DATA      CODE(25),KNTM(25),KNTR(25) / 'F012.1', 100,   0 /
	  DATA      CODE(26),KNTM(26),KNTR(26) / 'F406.1', 100,   0 /
	  DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.1', 100, 100 /
	  DATA      CODE(28),KNTM(28),KNTR(28) / 'G100.2', 100, 100 /
	  LOGF  =  LGFILE
	  IF(ERCODE .EQ. ' ')  THEN
	     L  =  0
	  ELSE
	     DO 10  L = 1, 6
		IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
  10		CONTINUE
  12	     CONTINUE
	  ENDIF
	  DO 14     I  =  1, KOUNTE
	     IF(L .EQ. 0)  GOTO 13
	     IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
  13	     KNTM(I)  =  LIMITM
	     KNTR(I)  =  LIMITR
  14	     CONTINUE
	  RETURN
	  ENTRY KERMTR_LHA(ERCODE,LOG,MFLAG,RFLAG)
	  LOG  =  LOGF
	  DO 20     I  =  1, KOUNTE
	     IF(ERCODE .EQ. CODE(I))  GOTO 21
  20	     CONTINUE
	  WRITE(*,1000)  ERCODE
	  CALL ABEND_LHA
	  RETURN
  21	  RFLAG  =  KNTR(I) .GE. 1
	  IF(RFLAG  .AND.  (KNTR(I) .LT. 100))  KNTR(I)  =  KNTR(I) - 1
	  MFLAG  =  KNTM(I) .GE. 1
	  IF(MFLAG  .AND.  (KNTM(I) .LT. 100))  KNTM(I)  =  KNTM(I) - 1
	  IF(.NOT. RFLAG)  THEN
	     IF(LOGF .LT. 1)  THEN
		WRITE(*,1001)  CODE(I)
	     ELSE
		WRITE(LOGF,1001)  CODE(I)
	     ENDIF
	  ENDIF
	  IF(MFLAG .AND. RFLAG)  THEN
	     IF(LOGF .LT. 1)  THEN
		WRITE(*,1002)  CODE(I)
	     ELSE
		WRITE(LOGF,1002)  CODE(I)
	     ENDIF
	  ENDIF
	  RETURN
1000	  FORMAT(' KERNLIB LIBRARY ERROR. ' /
     +  	 ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',
     +  	 ' ERROR MONITOR. RUN ABORTED.')
1001	  FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     +  	 'CONDITION ',A6)
1002	  FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
	  END
c========================================================================
      SUBROUTINE ABEND_LHA
C
C CERN PROGLIB# Z035	ABEND		.VERSION KERNVAX  1.10  811126

      STOP '*** ABEND ***'
      END
