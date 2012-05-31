      SUBROUTINE PHOMAK(IPPAR,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOMAK:   PHOtos MAKe
C.
C.    Purpose:  Single or double bremstrahlung radiative corrections  
C.              are generated in  the decay of the IPPAR-th particle in 
C.              the  HEP common /PH_HEPEVT/. Example of the use of 
C.              general tools.
C.
C.    Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
C.                                /PH_HEPEVT/ and the common itself
C.
C.    Output Parameters:  Common  /PH_HEPEVT/, either  with  or  without
C.                                particles added.
C.
C.    Author(s):  Z. Was,                         Created at:  26/05/93
C.                                                Last Update: 
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION DATA
      REAL*8 PHORAN
      INTEGER IP,IPPAR,NCHARG
      INTEGER WTDUM,IDUM,NHEP0
      INTEGER NCHARB,NEUDAU
      REAL*8 RN,WT,PHINT
      LOGICAL BOOST
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL INTERF,ISEC,IFTOP
      REAL*8 FINT,FSEC
      COMMON /PHOKEY/ FSEC,FINT,INTERF,ISEC,IFTOP
C--
      IP=IPPAR
      IDUM=1
      NCHARG=0
C--
        CALL PHOIN(IP,BOOST,NHEP0)
        CALL PHOCHK(JDAHEP(1,IP))
        WT=0.0D0
        CALL PHOPRE(1,WT,NEUDAU,NCHARB)
        IF (WT.EQ.0.0D0) RETURN
        RN=PHORAN(WTDUM)
C PHODO is caling PHORAN, thus change of series if it is moved before if
        CALL PHODO(1,NCHARB,NEUDAU)
        IF (INTERF) WT=WT*PHINT(IDUM)/FINT
        DATA=WT 
        IF (WT.GT.1.0D0) CALL PHOERR(3,'WT_INT',DATA)
      IF (RN.LE.WT) THEN 
        CALL PHOOUT(IP,BOOST,NHEP0)
      ENDIF
      RETURN
      END
