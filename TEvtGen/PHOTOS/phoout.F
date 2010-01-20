      SUBROUTINE PHOOUT(IP,BOOST,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOOUT:   PHOtos OUTput
C.
C.    Purpose:  copies back IP branch of the common /PH_HEPEVT/ from 
C.              /PHOEVT/ moves branch back from its CMS system.
C.
C.    Input Parameters:       IP:  pointer of particle starting branch
C.                                 to be given back.
C.                        BOOST:   Flag whether boost to CMS was or was 
C     .                            not performed.
C.
C.    Output Parameters:  Common /PHOEVT/, 
C.
C.    Author(s):  Z. Was                          Created at:  24/05/93
C.                                                Last Update:
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
      INTEGER IP,LL,FIRST,LAST,I
      LOGICAL BOOST
      INTEGER NN,J,K,NHEP0,NA
      DOUBLE PRECISION BET(3),GAM,PB
      COMMON /PHOCMS/ BET,GAM
      IF(NPHO.EQ.NEVPHO) RETURN
C--   When parent was not in its rest-frame, boost back...
      CALL PHLUPA(10)
      IF (BOOST) THEN
        DO 110 J=JDAPHO(1,1),JDAPHO(2,1)
          PB=-BET(1)*PPHO(1,J)-BET(2)*PPHO(2,J)-BET(3)*PPHO(3,J)
          DO 100 K=1,3
  100     PPHO(K,J)=PPHO(K,J)-BET(K)*(PPHO(4,J)+PB/(GAM+1.D0))
  110   PPHO(4,J)=GAM*PPHO(4,J)+PB
C--   ...boost photon, or whatever else has shown up
        DO NN=NEVPHO+1,NPHO
          PB=-BET(1)*PPHO(1,NN)-BET(2)*PPHO(2,NN)-BET(3)*PPHO(3,NN)
          DO 120 K=1,3
  120     PPHO(K,NN)=PPHO(K,NN)-BET(K)*(PPHO(4,NN)+PB/(GAM+1.D0))
          PPHO(4,NN)=GAM*PPHO(4,NN)+PB
        ENDDO
      ENDIF
        FIRST=JDAHEP(1,IP)
        LAST =JDAHEP(2,IP)
C let's take in original daughters
        DO LL=0,LAST-FIRST
         IDHEP(FIRST+LL) = IDPHO(3+LL)
           DO I=1,5
             PHEP(I,FIRST+LL) = PPHO(I,3+LL)
           ENDDO
        ENDDO
C let's take newcomers to the end of HEPEVT.
        NA=3+LAST-FIRST
        DO LL=1,NPHO-NA
         IDHEP(NHEP0+LL) = IDPHO(NA+LL)
         ISTHEP(NHEP0+LL)=ISTPHO(NA+LL)
         JMOHEP(1,NHEP0+LL)=IP
         JMOHEP(2,NHEP0+LL)=JMOHEP(2,JDAHEP(1,IP))
         JDAHEP(1,NHEP0+LL)=0
         JDAHEP(2,NHEP0+LL)=0
           DO I=1,5
             PHEP(I,NHEP0+LL) = PPHO(I,NA+LL)
           ENDDO
        ENDDO
        NHEP=NHEP+NPHO-NEVPHO
        CALL PHLUPA(20)
      END 
