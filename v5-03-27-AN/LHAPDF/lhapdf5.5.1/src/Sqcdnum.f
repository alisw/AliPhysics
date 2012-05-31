! -*- F90 -*-


! $Log$                                                                 
! Revision 1.3  2005/12/02 14:50:54  whalley                            
! Changes for new CTEQ code/AB sets                                     
!                                                                       
! Revision 1.2  2005/10/18 10:34:52  whalley                            
! small changes from cdf/d0 comments - renamubg conflicting name        
! adding pftopdg, etc.                                                  
!                                                                       
! Revision 1.1.1.1  2005/05/06 14:54:44  whalley                        
! Initial CVS import of the LHAPDF code and data sets                   
!                                                                       
! Revision 1.1.1.1  1996/04/01 15:02:05  mclareni                       
! Mathlib gen                                                           
!                                                                       
!                                                                       
      FUNCTION DDILOG_LHA(X) 
      implicit real*8 (a-h,o-z) 
      DIMENSION C(0:19) 
      PARAMETER (Z1 = 1d0, HF = Z1/2d0) 
      PARAMETER (PI = 3.14159265358979324D0) 
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12) 
                                                                        
      DATA C( 0) / 0.42996693560813697D0/ 
      DATA C( 1) / 0.40975987533077105D0/ 
      DATA C( 2) /-0.01858843665014592D0/ 
      DATA C( 3) / 0.00145751084062268D0/ 
      DATA C( 4) /-0.00014304184442340D0/ 
      DATA C( 5) / 0.00001588415541880D0/ 
      DATA C( 6) /-0.00000190784959387D0/ 
      DATA C( 7) / 0.00000024195180854D0/ 
      DATA C( 8) /-0.00000003193341274D0/ 
      DATA C( 9) / 0.00000000434545063D0/ 
      DATA C(10) /-0.00000000060578480D0/ 
      DATA C(11) / 0.00000000008612098D0/ 
      DATA C(12) /-0.00000000001244332D0/ 
      DATA C(13) / 0.00000000000182256D0/ 
      DATA C(14) /-0.00000000000027007D0/ 
      DATA C(15) / 0.00000000000004042D0/ 
      DATA C(16) /-0.00000000000000610D0/ 
      DATA C(17) / 0.00000000000000093D0/ 
      DATA C(18) /-0.00000000000000014D0/ 
      DATA C(19) /+0.00000000000000002D0/ 
                                                                        
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
!                                                                       
! $Id: Sqcdnum.f 365 2008-09-02 09:12:20Z buckley $                     
!                                                                       
! $Log$                                                                 
! Revision 1.3  2005/12/02 14:50:54  whalley                            
! Changes for new CTEQ code/AB sets                                     
!                                                                       
! Revision 1.2  2005/10/18 10:34:52  whalley                            
! small changes from cdf/d0 comments - renamubg conflicting name        
! adding pftopdg, etc.                                                  
!                                                                       
! Revision 1.1.1.1  2005/05/06 14:54:44  whalley                        
! Initial CVS import of the LHAPDF code and data sets                   
!                                                                       
! Revision 1.1.1.1  1996/04/01 15:02:13  mclareni                       
! Mathlib gen                                                           
!                                                                       
!                                                                       
      FUNCTION DGAUSS_LHA(F,A,B,EPS) 
                                                                        
      implicit real*8 (a-h,o-z) 
      CHARACTER NAME*(*) 
      PARAMETER (NAME = 'DGAUSS_LHA') 
      DIMENSION W(12),X(12) 
      PARAMETER (Z1 = 1d0, HF = Z1/2d0, CST = 5d0*Z1/1000d0) 
      DATA X                                                            &
     &        /0.96028985649753623168356086856947D0,              &
     &         0.79666647741362673959155393647583D0,              &
     &         0.52553240991632898581773904918925D0,              &
     &         0.18343464249564980493947614236018D0,              &
     &         0.98940093499164993259615417345033D0,              &
     &         0.94457502307323257607798841553461D0,              &
     &         0.86563120238783174388046789771239D0,              &
     &         0.75540440835500303389510119484744D0,              &
     &         0.61787624440264374844667176404879D0,              &
     &         0.45801677765722738634241944298358D0,              &
     &         0.28160355077925891323046050146050D0,              &
     &         0.95012509837637440185319335424958D-1/             
                                                                  
      DATA W                                                      &
     &        /0.10122853629037625915253135430996D0,              &
     &         0.22238103445337447054435599442624D0,              &
     &         0.31370664587788728733796220198660D0,              &
     &         0.36268378337836198296515044927720D0,              &
     &         0.27152459411754094851780572456018D-1,             &
     &         0.62253523938647892862843836994378D-1,             &
     &         0.95158511682492784809925107602246D-1,             &
     &         0.12462897125553387205247628219202D0,              &
     &         0.14959598881657673208150173054748D0,              &
     &         0.16915651939500253818931207903036D0,              &
     &         0.18260341504492358886676366796922D0,              &
     &         0.18945061045506849628539672320828D0/              
                                                                        
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
                                                                        
      SUBROUTINE VZERO_LHA (A,N) 
!                                                                       
! CERN PROGLIB# F121    VZERO           .VERSION KERNFOR  4.40  940929  
! ORIG. 01/07/71, modif. 24/05/87 to set integer zero                   
!                 modif. 25/05/94 to depend on QINTZERO                 
!                                                                       
      DIMENSION A(*) 
      IF (N.LE.0)  RETURN 
      DO 9 I= 1,N 
    9 A(I)= 0d0 
      RETURN 
      END                                           
                                                                        
!                                                                       
! $Id: Sqcdnum.f 365 2008-09-02 09:12:20Z buckley $                     
!                                                                       
! $Log$                                                                 
! Revision 1.3  2005/12/02 14:50:54  whalley                            
! Changes for new CTEQ code/AB sets                                     
!                                                                       
! Revision 1.2  2005/10/18 10:34:52  whalley                            
! small changes from cdf/d0 comments - renamubg conflicting name        
! adding pftopdg, etc.                                                  
!                                                                       
! Revision 1.1.1.1  2005/05/06 14:54:44  whalley                        
! Initial CVS import of the LHAPDF code and data sets                   
!                                                                       
! Revision 1.1.1.1  1996/02/15 17:49:49  mclareni                       
! Kernlib                                                               
!                                                                       
!                                                                       
      FUNCTION LENOCC_LHA (CHV) 
!                                                                       
! CERN PROGLIB# M507    LENOCC          .VERSION KERNFOR  4.21  890323  
! ORIG. March 85, A.Petrilli, re-write 21/02/89, JZ                     
!                                                                       
!-    Find last non-blank character in CHV                              
                                                                        
      CHARACTER    CHV*(*) 
                                                                        
      N = LEN(CHV) 
                                                                        
      DO 17  JJ= N,1,-1 
      IF (CHV(JJ:JJ).NE.' ') GO TO 99 
   17 END DO 
      JJ = 0 
                                                                        
   99 LENOCC_LHA = JJ 
      RETURN 
      END                                           
!                                                                       
! $Id: Sqcdnum.f 365 2008-09-02 09:12:20Z buckley $                     
!                                                                       
! $Log$                                                                 
! Revision 1.3  2005/12/02 14:50:54  whalley                            
! Changes for new CTEQ code/AB sets                                     
!                                                                       
! Revision 1.2  2005/10/18 10:34:52  whalley                            
! small changes from cdf/d0 comments - renamubg conflicting name        
! adding pftopdg, etc.                                                  
!                                                                       
! Revision 1.1.1.1  2005/05/06 14:54:44  whalley                        
! Initial CVS import of the LHAPDF code and data sets                   
!                                                                       
! Revision 1.1.1.1  1996/02/15 17:49:43  mclareni                       
! Kernlib                                                               
!                                                                       
!                                                                       
      SUBROUTINE CLTOU_LHA (CHV) 
!                                                                       
! CERN PROGLIB# M432    CLTOU           .VERSION KERNFOR  4.21  890323  
! ORIG. 11/02/86 A. PETRILLI                                            
! NEW    9/02/89 JZ, for speed                                          
!                                                                       
!-    Convert character string CHV from lower to upper case.            
                                                                        
      CHARACTER    CHV*(*) 
      DO 19  JJ=1,LEN(CHV) 
          J = ICHAR(CHV(JJ:JJ)) 
          IF (J.LT.97) EXIT
          IF (J.GE.123) EXIT
          CHV(JJ:JJ) = CHAR(J-32) 
   19 END DO 
      END                                           
!                                                                       
      SUBROUTINE TIMEX_LHA (T) 
!                                                                       
! CERN PROGLIB# Z007    TIMEX   DUMMY   .VERSION KERNFOR  4.05  821202  
!                                                                       
!-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE     
                                                                        
      T = 9. 
      RETURN 
      END                                           
!                                                                       
! $Id: Sqcdnum.f 365 2008-09-02 09:12:20Z buckley $                     
!                                                                       
! $Log$                                                                 
! Revision 1.3  2005/12/02 14:50:54  whalley                            
! Changes for new CTEQ code/AB sets                                     
!                                                                       
! Revision 1.2  2005/10/18 10:34:52  whalley                            
! small changes from cdf/d0 comments - renamubg conflicting name        
! adding pftopdg, etc.                                                  
!                                                                       
! Revision 1.1.1.1  2005/05/06 14:54:44  whalley                        
! Initial CVS import of the LHAPDF code and data sets                   
!                                                                       
! Revision 1.1.1.1  1996/02/15 17:49:44  mclareni                       
! Kernlib                                                               
!                                                                       
!                                                                       
      SUBROUTINE FLPSOR_LHA(A,N) 
!                                                                       
! CERN PROGLIB# M103    FLPSOR          .VERSION KERNFOR  3.15  820113  
! ORIG. 29/04/78                                                        
!                                                                       
!   SORT THE ONE-DIMENSIONAL FLOATING POINT ARRAY A(1),...,A(N) BY      
!   INCREASING VALUES                                                   
!                                                                       
!-    PROGRAM  M103  TAKEN FROM CERN PROGRAM LIBRARY,  29-APR-78        
!                                                                       
      DIMENSION A(N) 
      COMMON /SLATE/ LT(20),RT(20) 
      INTEGER R,RT 
!                                                                       
      LEVEL=1 
      LT(1)=1 
      RT(1)=N 
   10 L=LT(LEVEL) 
      R=RT(LEVEL) 
      LEVEL=LEVEL-1 
   20 IF(R.GT.L) GO TO 200 
      IF(LEVEL.LE.0) THEN 
        GO TO 50 
      ELSE 
        GO TO 10 
      ENDIF 
!                                                                       
!   SUBDIVIDE THE INTERVAL L,R                                          
!     L : LOWER LIMIT OF THE INTERVAL (INPUT)                           
!     R : UPPER LIMIT OF THE INTERVAL (INPUT)                           
!     J : UPPER LIMIT OF LOWER SUB-INTERVAL (OUTPUT)                    
!     I : LOWER LIMIT OF UPPER SUB-INTERVAL (OUTPUT)                    
!                                                                       
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
!                                                                       
  231 IF(I.GT.J) GO TO 232 
      W=A(I) 
      A(I)=A(J) 
      A(J)=W 
      I=I+1 
      J=J-1 
      IF(I.LE.J) GO TO 220 
!                                                                       
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
!                                                                       
! CERN PROGLIB# Z007    DATIMH  DUMMY   .VERSION KERNFOR  4.03  821008  
!                                                                       
!-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE     
                                                                        
      DIMENSION ND(9), NT(9) 
!      DIMENSION M(8)                                                   
                                                                        
      do i=1,9 
         ND(i)=0 
         NT(i)=0 
      enddo 
!      CALL UBLOW (8H29/09/79,M,8)                                      
!      CALL UBUNCH           (M,ND,8)                                   
!      CALL UBLOW (8H12.00.00,M,8)                                      
!      CALL UBUNCH           (M,NT,8)                                   
      RETURN 
      END                                           
!***********************************************************************
! these next added from CERNLIB to allow some photon and pion sets to wo
! nothing to do with QCDNUM mrw 9/12/2004                               
!***********************************************************************
!*                                                                      
! $Id: Sqcdnum.f 365 2008-09-02 09:12:20Z buckley $                     
!                                                                       
! $Log$                                                                 
! Revision 1.3  2005/12/02 14:50:54  whalley                            
! Changes for new CTEQ code/AB sets                                     
!                                                                       
! Revision 1.2  2005/10/18 10:34:52  whalley                            
! small changes from cdf/d0 comments - renamubg conflicting name        
! adding pftopdg, etc.                                                  
!                                                                       
! Revision 1.1.1.1  2005/05/06 14:54:44  whalley                        
! Initial CVS import of the LHAPDF code and data sets                   
!                                                                       
! Revision 1.1.1.1  1996/02/15 17:48:17  mclareni                       
! Kernlib                                                               
!                                                                       
!                                                                       
      DOUBLE PRECISION FUNCTION DGAMMA_LHA(X) 
      LOGICAL MFLAG,RFLAG 
      REAL SX 
      DOUBLE PRECISION X,U,F,ZERO,ONE,THREE,FOUR,PI 
      DOUBLE PRECISION C(0:24),H,ALFA,B0,B1,B2 
      DATA ZERO /0.0D0/, ONE /1.0D0/, THREE /3.0D0/, FOUR /4.0D0/ 
!#if defined(CERNLIB_NUMHIPRE)                                          
!      DATA NC /24/                                                     
!      DATA PI    /3.14159 26535 89793 23846 26433 83D0/                
!      DATA C( 0) /3.65738 77250 83382 43849 88068 39D0/                
!      DATA C( 1) /1.95754 34566 61268 26928 33742 26D0/                
!      DATA C( 2) / .33829 71138 26160 38915 58510 73D0/                
!      DATA C( 3) / .04208 95127 65575 49198 51083 97D0/                
!      DATA C( 4) / .00428 76504 82129 08770 04289 08D0/                
!      DATA C( 5) / .00036 52121 69294 61767 02198 22D0/                
!      DATA C( 6) / .00002 74006 42226 42200 27170 66D0/                
!      DATA C( 7) / .00000 18124 02333 65124 44603 05D0/                
!      DATA C( 8) / .00000 01096 57758 65997 06993 06D0/                
!      DATA C( 9) / .00000 00059 87184 04552 00046 95D0/                
!      DATA C(10) / .00000 00003 07690 80535 24777 71D0/                
!      DATA C(11) / .00000 00000 14317 93029 61915 76D0/                
!      DATA C(12) / .00000 00000 00651 08773 34803 70D0/                
!      DATA C(13) / .00000 00000 00025 95849 89822 28D0/                
!      DATA C(14) / .00000 00000 00001 10789 38922 59D0/                
!      DATA C(15) / .00000 00000 00000 03547 43620 17D0/                
!      DATA C(16) / .00000 00000 00000 00168 86075 04D0/                
!      DATA C(17) / .00000 00000 00000 00002 73543 58D0/                
!      DATA C(18) / .00000 00000 00000 00000 30297 74D0/                
!      DATA C(19) /-.00000 00000 00000 00000 00571 22D0/                
!      DATA C(20) / .00000 00000 00000 00000 00090 77D0/                
!      DATA C(21) /-.00000 00000 00000 00000 00005 05D0/                
!      DATA C(22) / .00000 00000 00000 00000 00000 41D0/                
!      DATA C(23) /-.00000 00000 00000 00000 00000 03D0/                
!      DATA C(24) / .00000 00000 00000 00000 00000 01D0/                
!#endif                                                                 
!#if defined(CERNLIB_NUMLOPRE)                                          
      DATA NC /15/ 
      DATA PI    /3.14159265358979324D0/ 
      DATA C( 0) /3.65738772508338244D0/ 
      DATA C( 1) /1.95754345666126827D0/ 
      DATA C( 2) / .33829711382616039D0/ 
      DATA C( 3) / .04208951276557549D0/ 
      DATA C( 4) / .00428765048212909D0/ 
      DATA C( 5) / .00036521216929462D0/ 
      DATA C( 6) / .00002740064222642D0/ 
      DATA C( 7) / .00000181240233365D0/ 
      DATA C( 8) / .00000010965775866D0/ 
      DATA C( 9) / .00000000598718405D0/ 
      DATA C(10) / .00000000030769081D0/ 
      DATA C(11) / .00000000001431793D0/ 
      DATA C(12) / .00000000000065109D0/ 
      DATA C(13) / .00000000000002596D0/ 
      DATA C(14) / .00000000000000111D0/ 
      DATA C(15) / .00000000000000004D0/ 
!#endif                                                                 
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
!=======================================================================
      SUBROUTINE KERSET_LHA(ERCODE,LGFILE,LIMITM,LIMITR) 
      PARAMETER(KOUNTE  =  28) 
      CHARACTER*6          ERCODE,    CODE(KOUNTE) 
      LOGICAL          MFLAG,    RFLAG 
      INTEGER          KNTM(KOUNTE),      KNTR(KOUNTE) 
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
   10       CONTINUE 
   12        CONTINUE 
      ENDIF 
      DO 14     I  =  1, KOUNTE 
         IF(L .EQ. 0)  GOTO 13 
         IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14 
   13        KNTM(I)  =  LIMITM 
         KNTR(I)  =  LIMITR 
   14        CONTINUE 
      RETURN 
      ENTRY KERMTR_LHA(ERCODE,LOG,MFLAG,RFLAG) 
      LOG  =  LOGF 
      DO 20     I  =  1, KOUNTE 
         IF(ERCODE .EQ. CODE(I))  GOTO 21 
   20        CONTINUE 
      WRITE(*,1000)  ERCODE 
      CALL ABEND_LHA 
      RETURN 
   21     RFLAG  =  KNTR(I) .GE. 1 
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
 1000     FORMAT(' KERNLIB LIBRARY ERROR. ' /                           &
     &       ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',             &
     &       ' ERROR MONITOR. RUN ABORTED.')                            
 1001     FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',       &
     &       'CONDITION ',A6)                                           
 1002     FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6) 
      END                                           
!=======================================================================
      SUBROUTINE ABEND_LHA 
!                                                                       
! CERN PROGLIB# Z035    ABEND        .VERSION KERNVAX  1.10  811126     
                                                                        
      STOP '*** ABEND ***' 
      END                                           
