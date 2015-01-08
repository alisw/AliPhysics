      BLOCK DATA DUMMY
C     Initialize the HEPEVT COMMON block
C     This helps to resolve the hepevt_ symbol when we link the library

#include "stdhep.inc"

      DATA nevhep/0/                   ! serial number
      DATA nhep/0/                     ! number of particles
      DATA (isthep(I),I=1,nmxhep)/nmxhep*0/           ! status code
      DATA (idhep(I),I=1,nmxhep)/nmxhep*0/            ! particle ident KF
      DATA (jmohep(1,I),I=1,nmxhep)/nmxhep*0/         ! parent particles
      DATA (jmohep(2,I),I=1,nmxhep)/nmxhep*0/         ! parent particles
      DATA (jdahep(1,I),I=1,nmxhep)/nmxhep*0/         ! childreen particles
      DATA (jdahep(2,I),I=1,nmxhep)/nmxhep*0/         ! childreen particles
      DATA (phep(1,I),I=1,nmxhep)/nmxhep*0/           ! four-momentum, mass [GeV]
      DATA (phep(2,I),I=1,nmxhep)/nmxhep*0/           ! four-momentum, mass [GeV]
      DATA (phep(3,I),I=1,nmxhep)/nmxhep*0/           ! four-momentum, mass [GeV]
      DATA (phep(4,I),I=1,nmxhep)/nmxhep*0/           ! four-momentum, mass [GeV]
      DATA (phep(5,I),I=1,nmxhep)/nmxhep*0/           ! four-momentum, mass [GeV]
      DATA (vhep(1,I),I=1,nmxhep)/nmxhep*0/           ! vertex [mm]
      DATA (vhep(2,I),I=1,nmxhep)/nmxhep*0/           ! vertex [mm]
      DATA (vhep(3,I),I=1,nmxhep)/nmxhep*0/           ! vertex [mm]
      DATA (vhep(4,I),I=1,nmxhep)/nmxhep*0/           ! vertex [mm]

      END
