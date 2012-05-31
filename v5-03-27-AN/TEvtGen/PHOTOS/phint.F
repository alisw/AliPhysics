
      FUNCTION PHINT(IDUM)
C.----------------------------------------------------------------------
C.
C.    PHINT:   PHotos INTerference
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
C.                                                Last Update: 
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
      REAL*8 MPASQR,XX,BETA,PQ1(4),PQ2(4),PPHOT(4)
      REAL*8 SS,PP2,PP,E1,E2,Q1,Q2,COSTHE
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
C check if two body was particle antiparticle (we improve on it !
C      IFINT= IFINT.AND.IDPHO(JDAPHO(1,1)).EQ.-IDPHO(IDENT)
C check if particles were charged
      IFINT= IFINT.AND.abs(PHOCHA(IDPHO(IDENT))).GT.0.01D0
C check if they have both charge
      IFINT= IFINT.AND.
     $       abs(PHOCHA(IDPHO(JDAPHO(1,1)))).gt.0.01D0
C calculates interference weight contribution
      IF(IFINT) THEN
        MPASQR = PPHO(5,1)**2
        XX=4.D0*MCHSQR/MPASQR*(1.-XPHOTO)/(1.-XPHOTO+(MCHSQR-MNESQR)/
     &     MPASQR)**2
         BETA=SQRT(1.D0-XX)
         PHINT  = 2D0/(1D0+COSTHG**2*BETA**2)
         SS =MPASQR*(1.D0-XPHOTO)
         PP2=((SS-MCHSQR-MNESQR)**2-4*MCHSQR*MNESQR)/SS/4
         PP =SQRT(PP2)
         E1 =SQRT(PP2+MCHSQR)
         E2 =SQRT(PP2+MNESQR)
         PHINT= (E1+E2)**2/((E2+COSTHG*PP)**2+(E1-COSTHG*PP)**2)
C
      q1=PHOCHA(IDPHO(JDAPHO(1,1)))
      q2=PHOCHA(IDPHO(IDENT))
      do k=1,4
       pq1(k)=ppho(k,JDAPHO(1,1))
       pq2(k)=ppho(k,JDAPHO(1,1)+1)
       pphot(k)=ppho(k,npho)
      enddo
       costhe=(pphot(1)*pq1(1)+pphot(2)*pq1(2)+pphot(3)*pq1(3))
       costhe=costhe/sqrt(pq1(1)**2+pq1(2)**2+pq1(3)**2)
       costhe=costhe/sqrt(pphot(1)**2+pphot(2)**2+pphot(3)**2)
C
! --- this IF checks whether JDAPHO(1,1) was MCH or MNE. 
! --- COSTHG angle (and in-generation variables) may be better choice 
! --- than costhe. note that in the formulae below amplitudes were 
! --- multiplied by (E2+COSTHG*PP)*(E1-COSTHG*PP). 
        IF (costhg*costhe.GT.0) then

         PHINT= (q1*(E2+COSTHG*PP)-q2*(E1-COSTHG*PP))**2
     &         /(q1**2*(E2+COSTHG*PP)**2+q2**2*(E1-COSTHG*PP)**2)
        ELSE

         PHINT= (q1*(E1-COSTHG*PP)-q2*(E2+COSTHG*PP))**2
     &         /(q1**2*(E1-COSTHG*PP)**2+q2**2*(E2+COSTHG*PP)**2)
        ENDIF
      ELSE
         PHINT  = 1D0
      ENDIF
         phint1=1
      END
