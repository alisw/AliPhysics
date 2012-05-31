      FUNCTION PHINT1(IDUM)
C.----------------------------------------------------------------------
C.
C.    PHINT:   PHotos INTerference (Old version kept for tests only.
C.
C.    Purpose:  Calculates interference between emission of photons from
C.              different possible chaged daughters stored in
C.              the  HEP common /PHOEVT/.  
C.
C.    Input Parameter:    commons /PHOEVT/ /PHOMOM/ /PHOPHS/
C.    
C.
C.    Output Parameters:  
C.                        
C.
C.    Author(s):  Z. Was,                         Created at:  10/08/93
C.                                                Last Update: 15/03/99
C.
C.----------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 PHINT,phint1
      REAL*8 PHOCHA
      INTEGER IDUM
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL*8 PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL*8 XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL*8 MPASQR,XX,BETA
      LOGICAL IFINT
      INTEGER K,IDENT 
C
      DO  K=JDAPHO(2,1),JDAPHO(1,1),-1
         IF(IDPHO(K).NE.22) THEN
           IDENT=K
           GOTO 20
         ENDIF
      ENDDO
 20   CONTINUE
C check if there is a photon
      IFINT= NPHO.GT.IDENT
C check if it is two body + gammas reaction
      IFINT= IFINT.AND.(IDENT-JDAPHO(1,1)).EQ.1
C check if two body was particle antiparticle
      IFINT= IFINT.AND.IDPHO(JDAPHO(1,1)).EQ.-IDPHO(IDENT)
C check if particles were charged
      IFINT= IFINT.AND.PHOCHA(IDPHO(IDENT)).NE.0
C calculates interference weight contribution
      IF(IFINT) THEN
        MPASQR = PPHO(5,1)**2
        XX=4.D0*MCHSQR/MPASQR*(1.D0-XPHOTO)/(1.D0-XPHOTO+(MCHSQR-MNESQR)
     &     /MPASQR)**2
         BETA=SQRT(1.D0-XX)
         PHINT  = 2D0/(1D0+COSTHG**2*BETA**2)
      ELSE
         PHINT  = 1D0
      ENDIF
       phint1=1
      END
