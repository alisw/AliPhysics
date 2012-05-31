      SUBROUTINE PHOIN(IP,BOOST,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOIN:   PHOtos INput
C.
C.    Purpose:  copies IP branch of the common /PH_HEPEVT/ into /PHOEVT/
C.              moves branch into its CMS system.
C.
C.    Input Parameters:       IP:  pointer of particle starting branch
C.                                 to be copied
C.                        BOOST:   Flag whether boost to CMS was or was 
C     .                            not performed.
C.
C.    Output Parameters:  Commons: /PHOEVT/, /PHOCMS/
C.
C.    Author(s):  Z. Was                          Created at:  24/05/93
C.                                                Last Update: 16/11/93
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER IP,IP2,I,FIRST,LAST,LL,NA
      LOGICAL BOOST
      INTEGER J,NHEP0
      DOUBLE PRECISION BET(3),GAM,PB
      COMMON /PHOCMS/ BET,GAM
      LOGICAL INTERF,ISEC,IFTOP
      REAL*8 FINT,FSEC
      COMMON /PHOKEY/ FSEC,FINT,INTERF,ISEC,IFTOP
C--
C let's calculate size of the little common entry
        FIRST=JDAHEP(1,IP)
        LAST =JDAHEP(2,IP)
        NPHO=3+LAST-FIRST+NHEP-NHEP0
        NEVPHO=NPHO
C let's take in decaying particle
           IDPHO(1)=IDHEP(IP)
           JDAPHO(1,1)=3
           JDAPHO(2,1)=3+LAST-FIRST
           DO I=1,5
             PPHO(I,1)=PHEP(I,IP)
           ENDDO
C let's take in eventual second mother
         IP2=JMOHEP(2,JDAHEP(1,IP))
         IF((IP2.NE.0).AND.(IP2.NE.IP)) THEN 
           IDPHO(2)=IDHEP(IP2)
           JDAPHO(1,2)=3
           JDAPHO(2,2)=3+LAST-FIRST
           DO I=1,5
             PPHO(I,2)=PHEP(I,IP2)
           ENDDO
         ELSE
           IDPHO(2)=0
           DO I=1,5
             PPHO(I,2)=0.0D0
           ENDDO
         ENDIF
C let's take in daughters
        DO LL=0,LAST-FIRST
           IDPHO(3+LL)=IDHEP(FIRST+LL)
           JMOPHO(1,3+LL)=JMOHEP(1,FIRST+LL)
           IF (JMOHEP(1,FIRST+LL).EQ.IP) JMOPHO(1,3+LL)=1
           DO I=1,5
             PPHO(I,3+LL)=PHEP(I,FIRST+LL)
           ENDDO
        ENDDO
        IF (NHEP.GT.NHEP0) THEN
C let's take in illegitimate daughters
        NA=3+LAST-FIRST 
        DO LL=1,NHEP-NHEP0
           IDPHO(NA+LL)=IDHEP(NHEP0+LL)
           JMOPHO(1,NA+LL)=JMOHEP(1,NHEP0+LL)
           IF (JMOHEP(1,NHEP0+LL).EQ.IP) JMOPHO(1,NA+LL)=1
           DO I=1,5
             PPHO(I,NA+LL)=PHEP(I,NHEP0+LL)
           ENDDO
        ENDDO
C--        there is NHEP-NHEP0 daugters more.
           JDAPHO(2,1)=3+LAST-FIRST+NHEP-NHEP0
        ENDIF
        CALL PHLUPA(1)
        CALL PHCORK(0)
C special case of t tbar production process
        IF(IFTOP) CALL PHOTWO(0)
        BOOST=.FALSE.
C--   Check whether parent is in its rest frame...
      IF (     (ABS(PPHO(4,1)-PPHO(5,1)).GT.PPHO(5,1)*1.D-8)
     $    .AND.(PPHO(5,1).NE.0))                            THEN
        BOOST=.TRUE.
C--
C--   Boost daughter particles to rest frame of parent...
C--   Resultant neutral system already calculated in rest frame !
        DO 10 J=1,3
   10   BET(J)=-PPHO(J,1)/PPHO(5,1)
        GAM=PPHO(4,1)/PPHO(5,1)
        DO 30 I=JDAPHO(1,1),JDAPHO(2,1)
          PB=BET(1)*PPHO(1,I)+BET(2)*PPHO(2,I)+BET(3)*PPHO(3,I)
          DO 20 J=1,3
   20     PPHO(J,I)=PPHO(J,I)+BET(J)*(PPHO(4,I)+PB/(GAM+1.D0))
   30   PPHO(4,I)=GAM*PPHO(4,I)+PB
C--    Finally boost mother as well
          I=1   
          PB=BET(1)*PPHO(1,I)+BET(2)*PPHO(2,I)+BET(3)*PPHO(3,I)
          DO J=1,3
            PPHO(J,I)=PPHO(J,I)+BET(J)*(PPHO(4,I)+PB/(GAM+1.D0))
          ENDDO
          PPHO(4,I)=GAM*PPHO(4,I)+PB
      ENDIF
C special case of t tbar production process
        IF(IFTOP) CALL PHOTWO(1)
      CALL PHLUPA(2)
      END 
