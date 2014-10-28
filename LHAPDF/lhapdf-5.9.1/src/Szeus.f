!-----------------------------------------------------------------------
      SUBROUTINE DVZERO (A,N) 
!                                                                       
! CERN PROGLIB# F121    VZERO           .VERSION KERNFOR  4.40  940929  
! ORIG. 01/07/71, modif. 24/05/87 to set integer zero                   
!                 modif. 25/05/94 to depend on QINTZERO                 
!                                                                       
      implicit real*8 (a-h,o-z) 
      DIMENSION A(*) 
      IF (N.LE.0)  RETURN 
      DO 9 I= 1,N 
    9 A(I)= 0d0 
      RETURN 
      END                                           
