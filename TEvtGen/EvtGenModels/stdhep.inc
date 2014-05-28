C*******************************************************
C  ---------------------------------------------------------------------
C     The stdhep event common
C  ---------------------------------------------------------------------
C This is the "official" event common of the HEP community for MC events.
C  NEVHEP      = event number
C  NHEP        = number of particles ("lines")
C  ISTHEP(N)   = stability code; 0 if stable, 1 if decays.
C  IDHEP(I)    = particle ID code; see PDG '92, III.73
C  JMOHEP(1,I) = line number of mother particle
C  JMOHEP(2,I) = line number of last mother particle, if formed by string
C  JDAHEP(1,I) = line number of 1st daughter particle
C  JDAHEP(2,I) = line number of last daughter particle - MUST BE CONSECUTIVE!
C  PHEP(5,I)   = px,py,pz,E,m of particle
C  VHEP(4,I)   = x,y,z,L of production vertex

      INTEGER  nmxhep         ! maximum number of particles
      PARAMETER ( nmxhep=10000)
      REAL*8  phep,  vhep ! to be real*4 or *8  depending on host
      INTEGER nevhep,nhep,isthep,idhep,jmohep,
     $        jdahep
      COMMON /hepevt/
     $      nevhep,               ! serial number
     $      nhep,                 ! number of particles
     $      isthep(nmxhep),   ! status code
     $      idhep(nmxhep),    ! particle ident KF
     $      jmohep(2,nmxhep), ! parent particles
     $      jdahep(2,nmxhep), ! childreen particles
     $      phep(5,nmxhep),   ! four-momentum, mass [GeV]
     $      vhep(4,nmxhep)    ! vertex [mm]

c$$$      INTEGER NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
c$$$      REAL    PHEP,VHEP
c$$$      PARAMETER (NMXHEP=2000)
c$$$      COMMON/STDHEP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c$$$     *JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
c$$$      SAVE /STDHEP/
