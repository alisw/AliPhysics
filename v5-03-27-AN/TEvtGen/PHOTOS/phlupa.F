      SUBROUTINE PHLUPA(IPOINT)
C.----------------------------------------------------------------------
C.
C.    PHLUPA:   debugging tool
C.
C.    Purpose:  NONE, eventually may printout content of the 
C.              /PHOEVT/ common
C.
C.    Input Parameters:   Common /PHOEVT/ and /PHNUM/ 
C.                        latter may have number of the event. 
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was                          Created at:  30/05/93
C.                                                Last Update: 10/08/93
C.
C.----------------------------------------------------------------------
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      COMMON /PHNUM/ IEV
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      DIMENSION SUM(5)
      IF (IPOINT.LT.3000) RETURN
      IOUT=56
      IF (IEV.LT.1000) THEN
      DO I=1,5
        SUM(I)=0.0D0
      ENDDO 
      WRITE(PHLUN,*) 'EVENT NR=',IEV,
     $            'WE ARE TESTING /PHOEVT/ at IPOINT=',IPOINT
      WRITE(PHLUN,10)
      I=1
      WRITE(PHLUN,20) IDPHO(I),PPHO(1,I),PPHO(2,I),PPHO(3,I),
     $                      PPHO(4,I),PPHO(5,I),JDAPHO(1,I),JDAPHO(2,I)
      I=2
      WRITE(PHLUN,20) IDPHO(I),PPHO(1,I),PPHO(2,I),PPHO(3,I),
     $                      PPHO(4,I),PPHO(5,I),JDAPHO(1,I),JDAPHO(2,I)
      WRITE(PHLUN,*) ' '
      DO I=3,NPHO
      WRITE(PHLUN,20) IDPHO(I),PPHO(1,I),PPHO(2,I),PPHO(3,I),
     $                      PPHO(4,I),PPHO(5,I),JMOPHO(1,I),JMOPHO(2,I)
        DO J=1,4
          SUM(J)=SUM(J)+PPHO(J,I)
        ENDDO
      ENDDO
      SUM(5)=SQRT(ABS(SUM(4)**2-SUM(1)**2-SUM(2)**2-SUM(3)**2))
      WRITE(PHLUN,30) SUM
 10   FORMAT(1X,'  ID      ','p_x      ','p_y      ','p_z      ',
     $                   'E        ','m        ',
     $                   'ID-MO_DA1','ID-MO DA2' )
 20   FORMAT(1X,I4,5(F9.3),2I9)
 30   FORMAT(1X,' SUM',5(F9.3))
      ENDIF
      END
