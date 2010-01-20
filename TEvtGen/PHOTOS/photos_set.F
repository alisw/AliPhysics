

      SUBROUTINE PHOTOS_SET
C.----------------------------------------------------------------------
C.
C.    Setter for PHOTOS:   
C.
C.    Purpose:  Copies /PH_HEPEVT/ into /HEPEVT/
C.              
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was                          Created at:  21/07/98
C.                                                Last Update: 21/07/98
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  d_h_nmxhep         ! maximum number of particles
      PARAMETER ( d_h_nmxhep=4000)
      REAL*8  d_h_phep,  d_h_vhep ! to be real*4 or *8  depending on host
      INTEGER d_h_nevhep,d_h_nhep,d_h_isthep,d_h_idhep,d_h_jmohep,
     $        d_h_jdahep
      COMMON /hepevt/
     $      d_h_nevhep,               ! serial number
     $      d_h_nhep,                 ! number of particles
     $      d_h_isthep(d_h_nmxhep),   ! status code
     $      d_h_idhep(d_h_nmxhep),    ! particle ident KF
     $      d_h_jmohep(2,d_h_nmxhep), ! parent particles
     $      d_h_jdahep(2,d_h_nmxhep), ! childreen particles
     $      d_h_phep(5,d_h_nmxhep),   ! four-momentum, mass [GeV]
     $      d_h_vhep(4,d_h_nmxhep)    ! vertex [mm]
* ----------------------------------------------------------------------
!      LOGICAL d_h_qedrad
!      COMMON /phoqed/ 
!     $     d_h_qedrad(d_h_nmxhep)    ! Photos flag
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL QEDRAD
      COMMON/PHOQED/QEDRAD(NMXHEP)
      INTEGER K,L

         d_h_nevhep=  nevhep             ! serial number
         d_h_nhep  =  nhep               ! number of particles
         DO K=1,nhep
           d_h_isthep(k)    =isthep(k)   ! status code
           d_h_idhep(k)     =idhep(k)    ! particle ident KF
           d_h_jmohep(1,k)  =jmohep(1,k) ! parent particles
           d_h_jdahep(1,k)  =jdahep(1,k) ! childreen particles
           d_h_jmohep(2,k)  =jmohep(2,k) ! parent particles
           d_h_jdahep(2,k)  =jdahep(2,k) ! childreen particles
           DO l=1,4
           d_h_phep(l,k)    =phep(l,k)   ! four-momentum, mass [GeV]
           d_h_vhep(l,k)    =vhep(l,k)   ! vertex [mm]
           ENDDO
           d_h_phep(5,k)    =phep(5,k)   ! four-momentum, mass [GeV]
         ENDDO
      END
