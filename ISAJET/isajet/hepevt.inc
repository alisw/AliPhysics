#ifndef CERNLIB_ISAJET_HEPEVT_INC
#define CERNLIB_ISAJET_HEPEVT_INC
*
*
* hepevt.inc
*
      INTEGER NMXHEP
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     $JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      REAL PHEP,VHEP
      SAVE /HEPEVT/
C...      NEVHEP       - event number
C...      NHEP         - number of entries in this event
C...      ISTHEP(..)   - status code
C...      IDHEP(..)    - particle ID, P.D.G. standard
C...      JMOHEP(1,..) - position of mother particle in list
C...      JMOHEP(2,..) - position of second mother particle in list
C...      JDAHEP(1,..) - position of first daughter in list
C...      JDAHEP(2,..) - position of last daughter in list
C...      PHEP(1,..)   - x momentum in GeV/c
C...      PHEP(2,..)   - y momentum in GeV/c
C...      PHEP(3,..)   - z momentum in GeV/c
C...      PHEP(4,..)   - energy in GeV
C...      PHEP(5,..)   - mass in GeV/c**2
C...      VHEP(1,..)   - x vertex position in mm
C...      VHEP(2,..)   - y vertex position in mm
C...      VHEP(3,..)   - z vertex position in mm
C...      VHEP(4,..)   - production time in mm/c

#endif
