*
*===program crint======================================================*
*
C      OPTIONS/ EXTEND_SOURCE
C      SUBROUTINE CRINT
      SUBROUTINE DT_PRODUCEEVENT(ENERGY_SL, NPARTICLES)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL ENERGY_SL
      INTEGER INIT
      REAL ne,etest,prob,slump
      SAVE

* Call the init sub routine in the first event
      DATA INIT /0/

      PARAMETER (NMXHKK=200000)

      COMMON /DTIONT/ LINP,LOUT,LDAT

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

*     event flag
      COMMON /DTEVNO/ NEVENT, ICASCA

      IF(INIT.EQ.0) THEN
         OPEN (UNIT = 50, file = "my.input")    
	 LINP = 50
         CALL DT_DTUINI(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)
*        Init called, make sure it's not called again
         INIT = 1
      ENDIF
*-----------------------------------------------------------------------
*     generation of one event
      NEVENT = 1
      KKMAT = -1

*   If an energy-range has been defined with the ENERGY input-card the
*   laboratory energy ELAB can be set to any value within that range,..
C        ELAB = DT_RNDM(EPN)*(EPN-0.5D7)+0.5D7

*   ..otherwise it has to coincide with EPN.
C        ELAB = EPN

      ELAB = ENERGY_SL

*   sampling of one event

*     TEST

      CALL DT_KKINC(NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,ELAB,KKMAT,IREJ)

      IF (IREJ.NE.0) RETURN

c     Return the number of particles produced
      
c     Fill the particle info 
      CALL DT_GETPARTICLES(NPARTICLES)

      END


      SUBROUTINE DT_GETPARTICLES(NPARTICLES)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER pid,qch,q_sum,Ntpc,Nfinal,NACCEPT,IPART,RES
      DOUBLE PRECISION yrap,pt,mass,mt,etot
      DOUBLE PRECISION pt_cut_tpc
      PARAMETER(pt_cut_tpc=0.050)

      SAVE
*
* COMMON /DTEVT1/ :
*                   NHKK         number of entries in common block
*                   NEVHKK       number of the event
*                   ISTHKK(i)    status code for entry i
*                   IDHKK(i)     identifier for the entry
*                                (for particles: identifier according
*                                 to the PDG numbering scheme)
*                   JMOHKK(1,i)  pointer to the entry of the first mother
*                                of entry i
*                   JMOHKK(2,i)  pointer to the entry of the second mother
*                                of entry i
*                   JDAHKK(1,i)  pointer to the entry of the first daughter
*                                of entry i
*                   JDAHKK(2,i)  pointer to the entry of the second daughter
*                                of entry i
*                   PHKK(1..3,i) 3-momentum
*                   PHKK(4,i)    energy
*                   PHKK(5,i)    mass
*
* event history

      PARAMETER (NMXHKK=200000)

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHIST(2,NMXHKK)

      DOUBLE PRECISION SLPX, SLPY, SLPZ, SLE, SLM
      INTEGER SLPID, SLCHARGE
      COMMON /DPMJETPARTICLE/ SLPX(NMXHKK), SLPY(NMXHKK), SLPZ(NMXHKK),
     &       SLE(NMXHKK), SLM(NMXHKK), SLPID(NMXHKK), SLCHARGE(NMXHKK)


C     >> Set Counter to Zero

      Nfinal=0
      
      DO 42 I=1, NHKK
c      I = IPART

CC       >> Remove all non-final-state particles
        IF(.not.(ISTHKK(I).eq.1.or.ISTHKK(I).eq.-1.or.
     $ISTHKK(I).eq.1001)) GOTO 42

C	>> Find Particle Charge, qch
        IF((ABS(ISTHKK(I)).eq.1).and.(IDHKK(I).ne.80000))THEN
C         >> final state ptcles except nuclei

          qch=IPHO_CHR3(IDHKK(I),1)/3
        ELSEIF(IDHKK(I).eq.80000)THEN
C         >> final state nuclei
          qch=IDXRES(I)
        ELSE
C         >> not a final state particle, qch not interesting
          qch=-999
        ENDIF

	Nfinal = Nfinal + 1
	SLPX(Nfinal) = PHKK(1,I)
        SLPY(Nfinal) = PHKK(2,I)
        SLPZ(Nfinal) = PHKK(3,I)
        SLE(Nfinal) = PHKK(4,I)
        SLM(Nfinal) = PHKK(5,I)
        SLPID(Nfinal) = IDHKK(I)
        SLCHARGE(Nfinal) = qch

 42     CONTINUE
        NPARTICLES = Nfinal
  
      END

      SUBROUTINE DT_USRHIS(MODE)
c Dummy to make the linker happy
      END

