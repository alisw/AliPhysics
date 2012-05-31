      SUBROUTINE PHOINF
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays general INFo
C.
C.    Purpose:  Print PHOTOS info
C.
C.    Input Parameters:   PHOLUN
C.
C.    Output Parameters:  PHOVN1, PHOVN2
C.
C.    Author(s):  B. van Eijk                     Created at:  12/04/90
C.                                                Last Update: 02/10/93
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IV1,IV2,IV3
      INTEGER PHOVN1,PHOVN2
      COMMON/PHOVER/PHOVN1,PHOVN2
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      LOGICAL INTERF,ISEC,IFTOP
      REAL*8 FINT,FSEC
      COMMON /PHOKEY/ FSEC,FINT,INTERF,ISEC,IFTOP
      REAL*8 ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
C--
C--   PHOTOS version number and release date
      PHOVN1=202
      PHOVN2=150399
C--
C--   Print info
      WRITE(PHLUN,9000)
      WRITE(PHLUN,9020)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9030)
      IV1=PHOVN1/100
      IV2=PHOVN1-IV1*100
      WRITE(PHLUN,9040) IV1,IV2
      IV1=PHOVN2/10000
      IV2=(PHOVN2-IV1*10000)/100
      IV3=PHOVN2-IV1*10000-IV2*100
      WRITE(PHLUN,9050) IV1,IV2,IV3
      WRITE(PHLUN,9030)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9060) 
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9070)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9020)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9064) INTERF,ISEC,IFTOP,ALPHA,XPHCUT
      WRITE(PHLUN,9010)
      IF (INTERF) WRITE(PHLUN,9061)
      IF (ISEC)   WRITE(PHLUN,9062)
      IF (IFTOP)  WRITE(PHLUN,9063)
      WRITE(PHLUN,9080)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9020)
      RETURN
 9000 FORMAT(1H1)
 9010 FORMAT(1H ,'*',T81,'*')
 9020 FORMAT(1H ,80('*'))
 9030 FORMAT(1H ,'*',26X,26('='),T81,'*')
 9040 FORMAT(1H ,'*',28X,'PHOTOS, Version: ',I2,'.',I2,T81,'*')
 9050 FORMAT(1H ,'*',28X,'Released at:  ',I2,'/',I2,'/',I2,T81,'*')
 9060 FORMAT(1H ,'*',18X,'PHOTOS QED Corrections in Particle Decays',
     &T81,'*')
 9061 FORMAT(1H ,'*',18X,'option with interference is active       ',
     &T81,'*')
 9062 FORMAT(1H ,'*',18X,'option with double photons is active     ',
     &T81,'*')
 9063 FORMAT(1H ,'*',18X,'emision in t tbar production is active   ',
     &T81,'*')
 9064 FORMAT(1H ,'*',18X,'Internal input parameters:',T81,'*'
     &,/,    1H ,'*',T81,'*'
     &,/,    1H ,'*',18X,'INTERF=',L2,'  ISEC=',L2,'  IFTOP=',L2,T81,'*'
     &,/,    1H ,'*',18X,'ALPHA_QED=',F8.5,'   XPHCUT=',F8.5,T81,'*')
 9070 FORMAT(1H ,'*',9X,'Monte Carlo Program - by E. Barberio, B. van Ei
     &jk and Z. Was',T81,'*',/,
     &      1H ,'*',9X,'From version 2.0 on - by E.B. and Z.W.',T81,'*')
 9080 FORMAT( 1H ,'*',9X,' ',T81,'*',/,
     &  1H ,'*',9X,
     & ' WARNING (1): /HEPEVT/ is not anymore the standard common block'
     & ,T81,'*',/,
     &  1H ,'*',9X,' ',T81,'*',/,
     &  1H ,'*',9X,
     & ' PHOTOS expects /HEPEVT/ to have REAL*4 variables. To change to'
     & ,T81,'*',/,  1H ,'*',9X,
     & ' REAL*8 modify its declaration in subr. PHOTOS_GET PHOTOS_SET:'
     & ,T81,'*',/,  1H ,'*',9X,
     & '      REAL*8  d_h_phep,  d_h_vhep'
     & ,T81,'*',/,  1H ,'*',9X,
     & ' WARNING (2): check dims. of /hepevt/ /phoqed/ /ph_hepevt/.'
     & ,T81,'*',/,  1H ,'*',9X,
     & ' HERE:                     d_h_nmxhep=4000  and  NMXHEP=10000'
     & ,T81,'*')
      END
