c$$$      PROGRAM MAIN
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$C
c$$$      parameter        (npawc=250000)
c$$$      COMMON/PAWC/ HMEM(npawc)
c$$$      COMMON/LUN / Lunst
c$$$      character*80 file_name
c$$$      data Lunst,LunF,LunD/ 10, 20, 21 /  ! Lun_File, Lun_Data
c$$$C      
c$$$      COMMON /pp2init/ Iproc1,Iproc2, sqrtS, aMmin, aMmax, Weight, ProcXsec(20)
c$$$      COMMON /pp2evnt/ ievent,Iproc,Npart,iParticleNumber,pParticle(5,10),iParticleCode(10),iParticleStatus(10) 
c$$$      data Iproc1,Iproc2,sqrtS,aMmin,aMmax / 2, 2, 7000., 0.270,  3.0 /
c$$$      data ProcXsec /      
c$$$C
c$$$C Iproc=   1     2        3     4       5               !   Process numbers
c$$$C       pi+pi-  f0(500)  rho   f0(980) f2(1270)         !   Process mnemonic
c$$$     +   3.0,    1.5,    1.0,   0.3,    1.0,     15*0.0 / ! Cross sections of the processes
c$$$C...........................................................................................
c$$$      CALL HLIMIT(npawc)
c$$$      CALL HROPEN(Lunst,'ions','pp2pxp.hbook','N',1024,IOSTAT)
c$$$      CALL BOOKIN
c$$$C      
c$$$      Nevnt = 100000
c$$$      Weight= 1.      !  All events will have Weight= 1.  
c$$$C-    Weight= 2.      !  The weighted events 
c$$$C   
c$$$      do iev=1,Nevnt
c$$$      call pp2pxp_gen
c$$$      if (1000*(iev/1000).eq.iev) write(*,*) 'iev,Iproc,Weight=',iev,Iproc,Weight
c$$$      enddo
c$$$C      
c$$$  100 call hrout(0,icycle,' ')
c$$$      call hrend('ions')
c$$$      close (Lunst)
c$$$C      
c$$$  200 format(a80)
c$$$      STOP'OK'
c$$$      END
c$$$C
c$$$      function pp2rnd(I)
c$$$      real rndm,pp2rnd
c$$$      pp2rnd=rndm(-1.)
c$$$      return
c$$$      end
C
      SUBROUTINE pp2pxp_gen      
C                                                                      S.Sadovsky 24.09.2013
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  ............................................................................................................   
      COMMON /pp2init/ Iproc1,Iproc2, sqrtS, aMmin, aMmax, Weight, ProcXsec(20), F2Polarization(8)
C     f2 polarization: |D0|^2, |D-|^2, |D+|^2, |D--|^2, |D++|^2, phase(D-,D0), phase(D--,D0), phase(D++,D+)
      COMMON /pp2evnt/ ievent,Iproc,Npart,iParticleNumber,pParticle(5,10),iParticleCode(10),iParticleStatus(10) 
C 
C     Iproc1 - first process type in the loop (continuum production of pipi, resonance production of pipi, ...)
C     Iproc2 - last process type in the process loop
C     Iproc  - process type of the current event
C     Npart  - number of particles in the final state of central system
C     sqrtS  - center of mass energy
C     aMmin  - minimal mass of centreally prodused system
C     aMmax  - maximal mass of centrally produced system
C     Weight - event Weight (if at initialisation stage Weight=1. then all generated events will have Weight=1.)   
C
C     ievent               -  the current event number
C     iParticleNumber      -  number of particles in the stack
C     pParticle(5,10)      -  array of particle 5-momentums (Px, Py, Pz, E, m)
C     iParticleCode(10)    -  array of particle PDG codes
C     iParticleStatus(10)  -  array of particle status codes
C
C     Status codes:
C     iParticleStatus =  1 - final particle
C                     = 21 - first Pomeron  (PDG code = 29) 
C                     = 22 - second Pomeron (PDG code = 29) 
C                     = 99 - centrally produced system (PDG code = 99) 
C
C     Stack structure:
C     1. Outgoing proton 1 with status code=1 (PDG code = 2212)
C     2. Outgoing proton 2 with status code=1 (PDG code = 2212)
C     3. Pomeron 1 with status code=21        (PDG code =  29)
C     4. Pomeron 2 with status code=22        (PDG code =  29)
C     5. Centrally produced system with status code=99 (PDG code = 99)
C     6. Final state particles with status code=1 
C
C     PDG codes of the final state particles:
C         e-  =   11
C         mu- =   13
C         pi+ =  211
C         pi0 =  111
C         K+  =  321
C         K0  =  311
C         K0s =  310
C         K0L =  130
C         p   = 2212
C  ................................................................................   
C
      dimension P(6,20), PX(6), P1(6), P2(6), Q1(6), Q2(6), POM1(6),POM2(6)
      save  PX, P1, P2 
C
      real Efmas , P_x, P_y, P_z, P_e, P_t, cos_GJ, phi_TY, Dm1, Dm2
      real WeighF,cosT,phi,eta,Y,phiD(10)     
      data Amp / 0.9382720  /            ! Proton mass and LHC CM energy
      data Npart    / 2 /
C
      logical Lstart
      data    Lstart / .false. /
C
      if(.not.Lstart ) then
                       ievent= 0		      
              P1(1) = 0.
	      P1(2) = 0.
	      P1(3) = dsqrt(0.25*sqrtS**2 - Amp**2)
	      P1(4) = 0.5*sqrtS
	      P1(5) = Amp
C	      
	      P2(1) = 0.
	      P2(2) = 0.
	      P2(3) =-dsqrt(0.25*sqrtS**2 - Amp**2)
	      P2(4) = 0.5*sqrtS
	      P2(5) = Amp
C
C --- Initial stack of the event p1p2 -> q1q2 PP => q1q2 PX  ---
C	
      iParticleCode  (1) = 2212    !  q1 - secondary proton  
      iParticleStatus(1) = 1
C
      iParticleCode  (2) = 2212    !  q2 - secondary proton
      iParticleStatus(2) = 1
C      	
      iParticleCode  (3) = 29     !  P  - Pomeron
      iParticleStatus(3) = 21
C
      iParticleCode  (4) = 29     !  P  - Pomeron
      iParticleStatus(4) = 22
C
      iParticleCode  (5) = 99      !  X  - centrally produces system
      iParticleStatus(5) = 99
C
      Istack_0 = 5
C      
      Lstart=.true.     
      endif 	     
C                  --- The main event cycle ---
C
      ievent= ievent+1  
      call pp2pxp(Iproc,P1,P2,Q1,Q2,PX,Npart,P,WeighD)
C  
C            --- Filling the rest of event stack ---
C
C     Stack structure:
C     1. Outgoing proton 1 with status code=1 (PDG code = 2212)
C     2. Outgoing proton 2 with status code=1 (PDG code = 2212)
C     3. Pomeron 1 with status code=21        (PDG code =  29)
C     4. Pomeron 2 with status code=22        (PDG code =  29)
C     5. Centrally produced system with status code=99 (PDG code = 99)
C     6. Final state particles with status code=1 (pi1)
C     7. Final state particles with status code=1 (pi2)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

c$$$      do j=1,4
c$$$      pParticle(j,1) = P1(j)-Q1(j)       !  Outgoing proton 1
c$$$      pParticle(j,2) = P2(j)-Q2(j)       !  Outgoing proton 1
c$$$      pParticle(j,3) = Q1(j)  !  Pomeron 1 
c$$$      pParticle(j,4) = Q2(j)  !  Pomeron 2 
c$$$      pParticle(j,5) = PX(j)        !  Centrally produced system
c$$$      enddo           
      do j=1,4
      pParticle(j,1) = Q1(j)        !  Outgoing proton 1
      pParticle(j,2) = Q2(j)        !  Outgoing proton 1
      pParticle(j,3) = P1(j)-Q1(j)  !  Pomeron 1 
      pParticle(j,4) = P2(j)-Q2(j)  !  Pomeron 2 
      pParticle(j,5) = PX(j)        !  Centrally produced system
      POM1(j)=P1(j)-Q1(j)
      POM2(j)=P2(j)-Q2(j)

      enddo           
C      
      pParticle(5,1) = Amp          !  Mass of Outgoing proton 1
      pParticle(5,2) = Amp          !  Mass of Outgoing proton 1
      pParticle(5,5) = PX(5)        !  Mass of Centrally produced system
c      POM1(5)=dsqrt(POM1(4)**2-POM1(3)**2-POM1(2)**2-POM1(1)**2)
c      POM1(5)=dsqrt(POM2(4)**2-POM2(3)**2-POM2(2)**2-POM2(1)**2)

C            
      do ipi=1,2
      do j=1,5
      pParticle(j,Istack_0+ipi) = P(j,ipi)
      enddo
      enddo
c$$$      write(*,*) pParticle(1,Istack_0+1), pParticle(1,Istack_0+2)
c$$$      write(*,*) pParticle(2,Istack_0+1), pParticle(2,Istack_0+2)
c$$$      write(*,*) pParticle(3,Istack_0+1), pParticle(3,Istack_0+2)
c$$$      write(*,*) pParticle(4,Istack_0+1), pParticle(4,Istack_0+2)
c$$$      write(*,*) pParticle(5,Istack_0+1), pParticle(5,Istack_0+2)

C
      iParticleCode  (6) = 211      !  pi+
      iParticleStatus(6) = 1
      
      iParticleCode  (7) =-211      !  pi-
      iParticleStatus(7) = 1
C      
      iParticleNumber    = Istack_0 + Npart
C     .....................................
C
C
C --- Histogramming: Centrally produces system X (PX - the 4-momentum of system) ---
C
c$$$      Efmas = PX(5)
c$$$      call hf1(100, Efmas, 1.)
c$$$      do j =1,Npart                 !  Loop over the particles
c$$$         do i =1,4
c$$$         Dm1 =  P(i,j)
c$$$         call hf1(100+i,Dm1,1.)   
c$$$         enddo
c$$$C
c$$$         call KinPar(P(1,j),cosT,phi,eta,Y)
c$$$         call hf1(111, cosT, 1.)
c$$$         call hf1(112, phi , 1.)
c$$$         call hf1(113, eta , 1.)
c$$$         call hf1(114, Y   , 1.)
c$$$         call hf2(115, eta , phi, 1.)
c$$$         phiD(j)=phi
c$$$      enddo
c$$$      call hf1(116, phiD(1)-phiD(2), 1.)
c$$$C      
c$$$      WeighF= WeighD
c$$$      P_x   = PX(1)
c$$$      P_y   = PX(2)
c$$$      P_z   = PX(3)
c$$$      P_e   = PX(4)
c$$$      P_t   = dsqrt(PX(1)**2-PX(2)**2)
c$$$C      
c$$$      call hf1(200, Efmas,WeighF )
c$$$      call hf1(201, P_x  ,WeighF )
c$$$      call hf1(202, P_y  ,WeighF )
c$$$      call hf1(203, P_z  ,WeighF )
c$$$      call hf1(204, P_e  ,WeighF )
c$$$      call hf1(205, P_t  ,WeighF )
c$$$C      
c$$$      call KinPar(PX(1),cosT,phi,eta,Y)
c$$$      call hf1(211, cosT, WeighF )
c$$$      call hf1(212, phi , WeighF )
c$$$      call hf1(213, eta , WeighF )
c$$$      call hf1(214, Y   , WeighF )
c$$$C
c$$$C --- PWA angular disributions ---
c$$$C
c$$$      ipos = 1
c$$$      if (P(6,2).gt.0.) ipos = 2  
c$$$C                  
c$$$      CALL PP2PXPi(P1,PX,P(1,ipos),P(1,3-ipos),COSGJ,PHITY)   ! pi+p-
c$$$      cos_GJ = COSGJ
c$$$      phi_TY = PHITY
c$$$C      
c$$$      if (abs(Y).lt.1.2) then   
c$$$                         call hf1(221,cos_GJ,1.)
c$$$                         call hf1(222,phi_TY,1.)
c$$$			 call hf2(223,cos_GJ,phi_TY,1.)
c$$$          endif 
C      
C --- Like charge particles ---
C
C-    if (P(6,1)/Dabs(P(6,1)).eq.P(6,2)/Dabs(P(6,2))) then
C-        call hf1(1200,Efmas,1.)
C-        call hf1(1201,P_x,1.)
C-        call hf1(1202,P_y,1.)
C-        call hf1(1203,P_z,1.)
C-        call hf1(1204,P_e,1.)
C-        call hf1(1205,P_t,1.)
C-        call hf1(1211,cos_GJ,1.)
C-        call hf1(1212,phi_TY,1.)
C-    endif
C      
C --- UnLike charge particles ---
C
C      if (P(6,1)/Dabs(P(6,1)).eq.-P(6,2)/Dabs(P(6,2))) then
C         call hf1(2200,Efmas, 1.)
C         call hf1(2201, P_x , 1.)
C         call hf1(2202, P_y , 1.)
C         call hf1(2203, P_z , 1.)
C         call hf1(2204, P_e , 1.)
C         call hf1(2205, P_t , 1.)
C	  call hf2(2207,Efmas ,P_t, 1.)
C	  call hf1(2211, cosT, 1.)
C         call hf1(2212, phi , 1.)
C         call hf1(2213, eta , 1.)
C         call hf1(2214, Y   , 1.)
C         call hf1(2221,cos_GJ,1.)
C         call hf1(2222,phi_TY,1.)  
C	 call hf2(2223,cos_GJ,phi_TY,1.)	 
C         call hf2(2224,Efmas ,cos_GJ,1.)
C	 call hf2(2225,Efmas ,phi_TY,1.)   
C      endif
C      
C --- UnLike pion pairs ---
C
C      if (P(6,1).eq.-P(6,2).and.Dabs(P(6,1)).eq.2) then
C         call hf1(4200,Efmas, 1.)
C         call hf1(4201, P_x , 1.)
C         call hf1(4202, P_y , 1.)
C         call hf1(4203, P_z , 1.)
C         call hf1(4204, P_e , 1.)
C         call hf1(4205, P_t , 1.)
C	  call hf1(4206, phiD(1)-phiD(2),1.)
C	  call hf2(4207, Efmas ,P_t, 1.)
C
C
C	  call hf1(4211, cosT, 1.)
C         call hf1(4212, phi , 1.)
C         call hf1(4213, eta , 1.)
C         call hf1(4214, Y   , 1.)
C	 
C         call hf1(4221,cos_GJ,1.)
C         call hf1(4222,phi_TY,1.)
C	  call hf2(4223,cos_GJ,phi_TY,1.)	 
C         call hf2(4224,Efmas ,cos_GJ,1.)
C	  call hf2(4225,Efmas ,phi_TY,1.)
C      endif               
C
      return
      end
C
      subroutine pp2pxp(Iproc,P1,P2,Q1,Q2,PX,Npart,P,WeighD)     
      implicit double precision (A-H,O-Z)   
C
      COMMON /pp2init/ Iproc1, Iproc2, sqrtS, aMmin, aMmax, Weight, ProcXsec(20), F2Polarization(8)
C     f2 polarization: |D0|^2, |D-|^2, |D+|^2, |D--|^2, |D++|^2, phase(D-,D0), phase(D--,D0), phase(D++,D+)
C     Iproc - process number
      dimension P1(6), P2(6), Q1(6), Q2(6), PX(6), P(6,20), Xsec(20), Prob(20)
      dimension POM1(6),POM2(6), PION1(6), PION2(6)
      real Xm, Weighf, X(10)
      real*8 Dzero, Dplusminus,cosgj,phity, UserPolarizationD
      data Nsec / 20 /                                  !  Maximum process number 
C
C     Resonanse parameters:
      dimension Jspn(10),ResM(10),ResW(10),WtMax(11) 
      data Nres,Jspn,ResM,ResW /   8,                   !  Number of resonances, Max = 10 
C
C     f0(500)  rho     f0(980)  f2(1270)
     +   0,     1,      0,       2,      2,      2,      2,      2,      0 ,    0,
     +   0.446, 0.7755, 0.990,   1.2751,  1.2751, 1.2751, 1.2751, 1.2751,    2*0.0,  
     +   0.350, 0.1491, 0.075,   0.1851,  0.1851, 0.1851, 0.1851, 0.1851,    2*0.0       /
C
      save Prob
C
      logical Lstart, LMinBias, Weight_1! Weight_1 means Weight=1. for all events
      data    Lstart, LMinBias, Weight_1, ampic  /.false., .false., .false., 0.13956755 /
C
      if (.not.Lstart) then
      Weight = 1.
      write(*,*) 'Generated processes numbers from to:',Iproc1, Iproc2
      write(*,*) '---------------------------------------------------------------------'
C      
      if (Iproc1.eq.Iproc2) then
                            Iproc=Iproc1
          else
	  LMinBias =.true.
          Nproc = Iproc2-Iproc1+1
	  if (Iproc2.gt.Nsec) stop 'Error in process number' 
C
C   ---   Process Probability calcelation ---
C
          SumP = 0.
          do j = Iproc1,Iproc2
          SumP = SumP + ProcXsec(j)
          enddo
C
          Sumj = 0. 
          do j = Iproc1,Iproc2
	  Sumj   = Sumj + ProcXsec(j)
          Prob(j)= Sumj / SumP
	  write(*,*) 'Channel probability',j,Prob(j)
          enddo
	  write(*,*)
      endif  
C      
C --- WtMax calculation --- 
C
      if(Weight.eq.1.) Weight_1 =.true.

      if(.not.Weight_1) then
                   write(*,*) ' Start generation of the Weighted events'
		   write(*,*) '========================================='
		   write(*,*)
	 else
         Ipro = 1
         Xmd  = 0.600 
         WtMax(Ipro) = 3.4*WtMax2pi(Xmd)
C
         do Ires= 1,Nres  
         Ipro  =  1+Ires 
         Jsp   = Jspn(Ires)
         aMres = ResM(Ires)
         aWres = ResW(Ires)
         Xmd   = aMres
         WtMax(Ipro) = 2.*3.4*WtMax2pi(Xmd)*BW_fun(aMres,aWres,Xmd,Ampic,Jsp)  
         enddo
C      
         write(*,*) ' All generated events have the Weight =1.'
	 write(*,*) '................................................' 
C	 
         do j=1,Nres+1
         write(*,*) 'Iproc,WtMax(j)=',j,WtMax(j)
         enddo
	 write(*,*) '================================================'
	 write(*,*) 
      endif     
C         
      Lstart = .true.
      endif
C      
C --- Choose the Process number ---
C
   10 if (LMinBias) then
               CALL PP2RANLUX(X(1),1) 
               do j= Iproc1,Iproc2
               if(X(1).lt.Prob(j)) go to 20
               enddo
   20 Iproc = j 
      endif 
C   
   30 call gener_cPPX( P1, P2, Q1, Q2, PX, WeighD)   
      Xmd = PX(5)
C
      Xm     = Xmd
      Weighf = WeighD
c$$$      call hf1(10, Xm, 1.) 
c$$$      call hf1(20, Xm, Weighf)
C
      if (Xmd.lt.aMmin.or.Xmd.gt.aMmax) go to 30
C      
c$$$      call hf1(27, Weighf, 1.) 
c$$$      call hf1(28, Xm, Weighf) 
C
      if (Iproc.gt.1) go to 200   ! Cicle over the processes
C
C     ---  pi+/pi- Continuum (Phase Space) ---
C
  100 Npart  = 2
      P(5,1) = ampic
      P(5,2) = ampic	   
      if ( PX(5) .lt. P(5,1)+P(5,2) ) go to 30  
C
      Wpipi  = WPP2pi(Xmd)             
      WeighD = WeighD*Wpipi 
C
      if (Weight_1) then
                    CALL PP2RANLUX(X(1),1)
		    if(X(1)*WtMax(Iproc).gt.WeighD) go to 30
C
		    if(WtMax(Iproc).lt.WeighD) then 
                    write(*,*) 'WtMax, Wt=',WtMax(Iproc),WeighD
c		    stop'Err. WeighD'
		    endif
      WeighD = 1.
      endif				    
      go to 1000
  
C     ---  pi+/pi- resonances  ---
C  
  200 Npart  = 2
      P(5,1) = ampic
      P(5,2) = ampic
      if ( PX(5) .lt. P(5,1)+P(5,2) ) go to 10  
            
      Wpipi  = WPP2pi(Xmd)
      WeighD = WeighD*Wpipi         
C      
      Ires   = Iproc-1
      Jsp    = Jspn(Ires)
      aMres  = ResM(Ires)
      aWres  = ResW(Ires)
C
      WeighR = BW_fun(aMres, aWres, Xmd, Ampic, Jsp)
      WeighD = WeighD*WeighR
c      write(*,*) 'DRgen: I started to generate resonance'
C
      if (Weight_1) then
                    CALL PP2RANLUX(X(1),1)
		    if(X(1)*WtMax(Iproc).gt.WeighD) go to 30
C 
                    if(WtMax(Iproc).lt.WeighD) then 
                    write(*,*) 'WtMax, Wt=',WtMax(Iproc),WeighD
c		    stop'Err. WeighD'
		    endif
      WeighD = 1.
      endif				    
C      
 1000 Weight = WeighD
 1001 call decays(PX,P(1,1),P(1,2))
C 
      if(Ires.ge.5)then         !polarized f2
         do j=1,4
            POM1(j)=P1(j)-Q1(j)
            POM2(j)=P2(j)-Q2(j)
            PION1(j)=P(j,1)
            PION2(j)=P(j,2)
         enddo
c         POM1(5) = dsqrt(POM1(4)**2-POM1(3)**2-POM1(2)**2-POM1(1)**2)
c         POM2(5) = dsqrt(POM2(4)**2-POM2(3)**2-POM2(2)**2-POM2(1)**2)
         PION1(5)= ampic
         PION2(5)= ampic
         if(Ires.eq.5)then      !f2(1270) D0 wave in pom-pom system
            call PP2PXPi(POM1,PX,PION1,COSGJ,PHITY)
            WeiD=Dzero(COSGJ,PHITY)
         endif
         if(Ires.eq.6)then      !f2(1270) D+D- wave in pom-pom system
            call PP2PXPi(POM1,PX,PION1,COSGJ,PHITY)
            WeiD=Dplusminus(COSGJ,PHITY)
         endif
         if(Ires.eq.7)then      !f2(1270) D0 wave in prot-prot system
            call PP2PXPi(P1,PX,PION1,COSGJ,PHITY)
            WeiD=Dzero(COSGJ,PHITY)
         endif
         if(Ires.eq.8)then      !f2(1270) D+D- wave in prot-prot system
            call PP2PXPi(P1,PX,PION1,COSGJ,PHITY)
            WeiD=Dplusminus(COSGJ,PHITY)
         endif
         if(Ires.eq.9)then      !f2(1270) user-setted polarization in pom-pom system
            call PP2PXPi(POM1,PX,PION1,COSGJ,PHITY)
            WeiD=UserPolarizationD(COSGJ,PHITY, F2Polarization)
         endif
         if(Ires.eq.10)then      !f2(1270) user-setted polarization in prot-prot system
            call PP2PXPi(P1,PX,PION1,COSGJ,PHITY)
            WeiD=UserPolarizationD(COSGJ,PHITY, F2Polarization)
         endif
         CALL PP2RANLUX(X(1),1)
         if(X(1).gt.WeiD) go to 1001
C     
         if(1.0.lt.WeiD) then 
            write(*,*) 'Wt=',WeiD, 'wmax=1.'
c     evd		    stop'Err. WeighD'
         endif
      endif
      Weight = WeighD
c
      return
      end  
C
      function WtMax2pi(X)
      implicit double precision (A-H,O-Z)
C      
      dimension p(10)
C-    data p /  93274., -3.0002, 3359.8, -7733.2, 7694.7,
      data p /  9.3274, -3.0002, 3359.8, -7733.2, 7694.7,
     +         -3115.6, -80.988, 510.56, -159.83, 16.601 /  
C      
      WtMax2pi= X*X*P(1)*exp(p(2)*X) * (1.+ p(3)*X + p(4)*X**2 +
     +  p(5)*X**3 + p(6)*X**4 + p(7)*X**5 + p(8)*X**6 + p(9)*X**7 + p(10)*X**8) 
C
      return
      end  
C
      function WtMax2pi_old(X)      !   M.Albrow
      implicit double precision (A-H,O-Z)
C      
      dimension par(12)
      data par /  -0.19764E+02, -3.3263,   0.26724,  -0.011915,
     +		 -10.339,       -13.424,    11.925,    -7.4001,
     +		  0.090513,    0.082573, -0.033527, -0.0041702/
C
      WtMax2pi= par(1)*x*exp(par(2)*X+par(3)*X**2+par(4)*X**3)*
     +(0.1+par(5)*X+par(6)*X**2+par(7)*X**3+par(8)*X**4+par(9)*X**5+
     + par(10)*X**6+par(11)*X**7+par(12)*X**8) 
C
C-    WtMax2pi=p(1)       !  Weight = Const.
      return
      end         
C
C
      function WPP2pi(Xmd)
C
C     Diff. cross section of PP --> pi+pi- 
C
      implicit double precision (A-H,O-Z)
      data ampic / 0.13956755 / 
C
      if (Xmd.lt.2.*ampic) then
                           WPP2pi = 0.
          else		   
          p_pi   = dsqrt( dabs(Xmd*Xmd/4. - ampic**2) )
          WPP2pi = p_pi/Xmd * dexp( -2.*p_pi**2 )
      endif
      return
      end             
C      
      function BW_fun(Rmas, Rwid, Xm, Ampi, Jspn)
      implicit double precision (A-H,O-Z) 
C     The classical relativistic BW function for spin 0:  S.Sadovsky 17.08.2011
C                 
           if (Jspn.eq.2) then 
                          BW_fun = BW_f2_gams(Rmas, Rwid, Xm, Ampi)
      else if (Jspn.eq.1) then 
	                  BW_fun = BW_omg_gams(Rmas, Rwid, Xm, Ampi)
                          else         
			  BW_fun = BW_fun0    (Rmas, Rwid, Xm, Ampi)
                          endif
      return
      end 
C
      function BW_fun0(Rmas, Rwid, Xm, Ampi)
      implicit double precision (A-H,O-Z) 
C     The classical relativistic BW function for spin 0:  S.Sadovsky 17.08.2011
C                 
      BW_fun0 = Xm*Rmas*Rwid**2/((Xm*Xm-Rmas**2)**2+(Rmas*Rwid)**2)
      return
      end   
C
      function BW_omg_gams(Rmas, Rwid, Xm, Ampi)
      implicit double precision (A-H,O-Z) 
C
C     The GAMS relativistic BW function for spin 1 part.:  S.Sadovsky 17.08.2011 
C     Eur.Phys.J. A3, 361-371 (1998)
C     Input:   Xm   - 2pi mass in GeV
C              Rmas - resonanse mass  in GeV
C              Rwid - resonanse width in GeV
C
C     data r0  /  197.326 /  !  r0=1 Fm (in 1/MeV)     - hadron radius in Blatt-Weiskopf factor 
      data r0  /0.197326  /  !  r0=1 Fm (in 1/GeV) 
C                        					      
      ppi0  = sqrt(abs((Rmas/2.)**2-Ampi**2))
      ppiX  = sqrt(abs((Xm/2.)**2  -Ampi**2))
C     
      D2pi0 = 1.+(ppi0/r0)**2	                   ! Blatte-Weiskopf factor for omega
      D2piX = 1.+(ppiX/r0)**2
C
      Gamma = Rwid*(ppiX/ppi0)**3*D2pi0/D2piX                  
      BW    = Rmas**2*Gamma**2/((Xm*Xm-Rmas**2)**2+(Rmas*Gamma)**2)
      BW_omg_gams = Bw*Xm**2/ppiX
      return
      end      
C 
      function BW_f2_gams(Rmas, Rwid, Xm, Ampi)
      implicit double precision (A-H,O-Z) 
C      
C     The GAMS relativistic BW function for spin 2 part.:  S.Sadovsky 17.08.2011 
C     Eur.Phys.J. A3, 361-371 (1998)
C     Input:   Xm   - 2pi mass in MeV
C              Rmas - resonanse mass  in GeV
C              Rwid - resonanse width in GeV
C
C     data r0  /  197.326 /  !  r0=1 Fm (in 1/MeV)
      data r0  /0.197326  /  !  r0=1 Fm (in 1/GeV)
C                        					      
      ppi0  = sqrt(abs((Rmas/2.)**2-Ampi**2))
      ppiX  = sqrt(abs((Xm/2.)**2  -Ampi**2))
C     
      D2pi0 = 9.+3.*(ppi0/r0)**2+(ppi0/r0)**4
      D2piX = 9.+3.*(ppiX/r0)**2+(ppiX/r0)**4
C
      Gamma = Rwid*(ppiX/ppi0)**5*D2pi0/D2piX                  
      BW    = Rmas**2*Gamma**2/((Xm*Xm-Rmas**2)**2+(Rmas*Gamma)**2)
      BW_f2_gams = Bw*Xm**2/ppiX
      return
      end 
C      
      subroutine pp2pxp_old(P1,P2,PX,Npart,P,WeighD)
      implicit double precision (A-H,O-Z)      
      dimension P1(6), P2(6), PX(6), P(6,10)
      real Xm
C
      COMMON/RLRCOM3/ PLAB(4,30)
C
      call generate_all(WeighD)
      xm = dsqrt( plab(4,7)**2-plab(1,7)**2-plab(2,7)**2-plab(3,7)**2 )

C-    write(*,*) 'Xm=', Xm
C
c$$$      call hf1(10, Xm, 1.) 
      return
      end  

C-------------------------------------------------------------------------------------------------

      subroutine generate_all(Weight_all)

      real *8 weight_all
      real *8 weightd
      real *8 pbeam(4), pfast(4),  pslow(4), pres(4)

      real *8 xmd

      real *8 pbeam1(4), pbeam1_cm(4), pslow_cm(4)

c-----------  use rndh  (integration) or pp2ranlux (for weight=1 events) 
      common /mode_random/ imodrn

c---------- from ffread cards ----------------------
      common /reaction_type/ ityper
      common /fill_decay/ ihist_decay
      common /mass_limits/ xma, xmb
      common /mass_limits_1/ xma1, xmb1, xmstep1 ! course - binned  histograms (== PWA-bin)
      common /mass_limits_2/ xma2, xmb2, xmstep2 ! fine-binned histograms ( 4-MeV )

c------- temporary slots
      common /ampl_tmp/ name_ampl_tmp
      character *60 name_ampl_tmp
      common /bw_tmp/ cn_re_tmp, cn_im_tmp, parbw_tmp(20)
      common /begincoh_tmp/ nrank_tmp
      common /phase_space_only/  ips, iwaves

c-------- cross section structure integer and floating parameters, defined from cards
      parameter  (ncohmax=20,nbwmax=10, nparbwmax=15 )
      parameter  ( nampmax = 100 )
      common /cross_sect_struct_int/ ncoh, namp(ncohmax),  
     *              nbw( ncohmax,nampmax), iw(ncohmax,nampmax)
      common /cross_sect_struct_float/ 
     *  parbw(nparbwmax, ncohmax, nampmax, nbwmax), 
     *        cn(ncohmax, nampmax, nbwmax) 
      complex cn

      common /non_repeated_waves/ namp1, iw1(nampmax)

c-------- weight for n-meson phase space ---------
      common /gener_weight/ weight_gen

c----------- amplitudes, names,   quantum numbers from calc_ampl (physics-dependent)
      common /calc_ampl_amplits/ ampl(nampmax)
      complex  ampl
      common /calc_ampl_samplits/ sampl(nampmax)
      real    sampl

      common /calc_ampl_names/ name_calc_ampl(nampmax)
      character *60  name_calc_ampl
      
c------------------------------------ total spin -- spin projection -- reflectivity(if def.) --
      common /calc_ampl_quantum_numb/ J1(nampmax), JM(nampmax), IP(nampmax), IETA(nampmax)  
      common /calc_ampl_totals/ namp_tot

c------------------ stored during initialization  integrals and max weights
      parameter  ( nmasspnt = 3000 )
      common /samplits/  samp( nampmax, nmasspnt )
      common /weight_max/ weight_max(nmasspnt)
      common /weight_sum/ ssigma(nmasspnt)
      common /weight_sum_max/ ssigma_max

      COMMON/RLRCOM3/PLAB(4,30)
      real *8 PLAB

      common /decay_lab/  pdec_lab(4,10)
      real *8  pdec_lab

      common /decay_cms_momentum/ pdec(4,10) 
      real *8 pdec

      common /decay_cms_nmesons/ nmes 

      real *8 pdec1(4,10)
      real *8 anGJ(3,3)

      common /ampl2_contrib/ sampl2(nampmax)
      real sampl2

      common /ampl3_contrib/ sampl3(ncohmax, nampmax, nbwmax )
      real sampl3

      common /total_intensities/  sintens(nampmax)
      real sintens 

      common /total_density_matrix/  rho_big(nampmax, nampmax )
      complex rho_big

      complex rho1

      parameter(nbmax1 = 300)
      real val1(nbmax1,nampmax, nampmax ),val2(nbmax1,nampmax, nampmax ),val3(nbmax1, nampmax, nampmax),
     *     val4(nbmax1,nampmax, nampmax ),val9(nbmax1,nampmax, nampmax )
      real val0(nbmax1, nampmax ), valtot(nbmax1)

      parameter (nbmax2 = 2000)
      real val1a(nbmax2,nampmax, nampmax ),val2a(nbmax2, nampmax, nampmax), val3a(nbmax2, nampmax, nampmax ),
     *     val4a(nbmax2,nampmax, nampmax ),val9a(nbmax2, nampmax, nampmax )
      real val0a(nbmax2,nampmax), valtota(nbmax2)

      common /decay_channel/ ichannel

      common /inits/ init0, init1, init2, init3, init4, init5
      common /fast_mode/ ifast_mode

      common /nmc_for_histos_normalization/ nmc_norm
      common /nmc_for_fine_binned_waves/ nmc_fine
      common /nmc_for_fine_binned_production_model/ nmc_prod

c-    common /x_limits/  xma1, xmb1,  xmstep1, nmsteps1
      common /z_limits/ z1min, z1max, z2min, z2max, ztotmin, ztotmax, ntot_z

      common / sampl_integrals_tab / ama_tab, amb_tab, stepx_tab,  ps_tab( 3000, 100 )
      real ama_tab, amb_tab, stepx_tab, ps_tab

      common /normalize_amp/ inormal
      
      common /pbeam_f/ pbeamf, sigma_pbeam
      
      character *80 namtmp2
      logical Lbook,Lbook2 

      common /production_type_ffread/ idp, icp, ilp, idp1, ilp1
      data idp, icp, ilp, idp1, ilp1 / 0,   1,   0,    0,    0 /       ! Sdv+ 
      data init0,  Lbook, Lbook2     /  1, .false.,  .false.   /       ! Sdv+ 
      data inormal, init5 /  1,  1       /                             ! Sdv+ 
      data namp_tot, ichannel, iwaves    /  1, 2, 2 /                  ! Sdv+ 
C      
      data xma,  xmb           / 0.270, 1.270 /       ! Sdv+  common /mass_limits/ xma, xmb	      
      data xma1, xmb1, xmstep1 / 0.270, 1.270, 0.04 / ! Sdv+  common /mass_limits_1/ xma1, xmb1, xmstep1 ! course - binned  histograms (== PWA-bin)
      data xma2, xmb2, xmstep2 / 0.270, 1.270, 0.02 / ! Sdv+  common /mass_limits_2/ xma2, xmb2, xmstep2 ! fine-binned histograms ( 4-MeV 

      data pbeamf, sigma_pbea  / 300., 0.5 /

      ifast_mode = 0   
C-    write(*,*) 'ifast_mode=',ifast_mode
      if(init0 .ne. 0) then

c-------- choosing of production type ... -------------
 
          ityper = -1
          nr_up = 0

          if( idp .ne. 0) then 
            ityper = 1 
            nr_up = nr_up + 1
          endif

          if( icp .ne. 0) then 
            ityper = 2 
            nr_up = nr_up + 1
          endif

          if( ilp .ne. 0) then 
            ityper = 3 
            nr_up = nr_up + 1
          endif

          if( idp1 .ne. 0) then 
            ityper = 4 
            nr_up = nr_up + 1
          endif

          if( ilp1 .ne. 0) then 
            ityper = 5 
            nr_up = nr_up + 1
          endif

          write(0,*) ' ityper = ', ityper

          if( ityper .eq. -1) then
             write(0,*) ' production mechanism: DP, CP, LP or DP1 is not chosen in cards, stop '
             stop
          endif

          if( nr_up .ne. 1) then
             write(0,*) 'More than one production mechanism: DP, CP, LP or DP1 in cards,  stop '
             stop
          endif
C-------------------------------------------------------------------------------------------------
CSdv+-
       write(*,*) ' inormal= ',inormal

       if(inormal .ne. 0) then

          imodrn = 1

          xma_sav = xma
          xmb_sav = xmb

          ama_tab = xma - 0.010
          amb_tab = xmb + 0.010
          stepx_tab = 0.010
          nst = ( amb_tab - ama_tab)/stepx_tab + 0.5

          write(0,*) ' ama_tab  = ', ama_tab,' amb_tab= ', amb_tab, ' nst = ', nst
          write(0,*) ' namp_tot = ', namp_tot
          write(0,*) ' ichannel = ', ichannel
          write(0,*) ' ityper   = ', ityper

          xma = ama_tab 
          xmb = amb_tab 
          write(0,*) ' xma = ', xma, ' xmb= ', xmb

          do k1 = 1, nst
          do i  = 1, namp_tot
            ps_tab(k1,i) = 0.
          enddo
          enddo

CSdv-     nmc_normal = 300000
          nmc_normal = 30

          do k = 1, nmc_normal     

          if( ityper .eq. 1) call gener_dp2( xmd,  weightd    ) ! xmd is OUTPUT here !!!
          if( ityper .eq. 2) call gener_cp2( xmd,  weightd    ) ! xmd is OUTPUT here !!!
CSdv-     if( ityper .eq. 3) call gener_lp2( xmd,  weightd, 0 ) ! xmd is OUTPUT here !!!
CSdv-     if( ityper .eq. 5) call gener_lp2( xmd,  weightd, 1 ) ! xmd is OUTPUT here !!!
          if( ityper .eq. 4) call gener_dp1( xmd,  weightd    ) ! xmd is OUTPUT here !!!
              
          call gener_calc_ampl_all(ichannel, 1, xmd)  ! iwaves = 1 always, here we have decay amplitudes       
c--------------------------------------------------------------------------------------------------------
          k1 = (xmd - ama_tab)/stepx_tab + 1.0

          write(*,*) ' xmd = ', xmd, ' ama_tab = ', ama_tab, ' stepx_tab = ' , stepx_tab, ' k1 = ', k1
              
              do i = 1, namp_tot
                ps_tab(k1,i) =  ps_tab(k1,i)  +  weightd *  weight_gen * sampl(i)
              enddo
              
            enddo

             do i  = 1, namp_tot
             do k1 = 1, nst
                ps_tab(k1,i) =  ps_tab(k1,i) / float(nmc_normal)*float(nst)
                  write(0,*) ' k1 = ', k1, ' i= ', i,' ps_tab = ', ps_tab(k1,i)
             enddo
             enddo

            xma = xma_sav ! return full-mass interval values !!!!!!!
            xmb = xmb_sav

          imodrn = 0
      inormal = 0
      endif   

c==========================================================================

          if( xma1 .eq. xmb1 ) then
              xma1 = xma
              xmb1 = xmb
          endif

          nmsteps1 =  (xmb1-xma1)/xmstep1+0.001

c========================================================
          if( xma2 .eq. xmb2 ) then
              xma2 = xma1
              xmb2 = xmb1
          endif

      if (.not.Lbook)  then
      write(*,*) 'call hbook1', xma1, xmb1, xmstep2

      nmsteps2 =  (xmb2-xma2)/xmstep2+0.001

C-      call hbook1( 9001,'M(X) (prod. kinem)x(dec.ph.sp.)x(waves) $',nmsteps1, xma1, xmb1, 0.)
C-      call hbook1( 9005,'M(X) (prod. kinem)x(dec.ph.sp.)         $',nmsteps1, xma1, xmb1, 0.)
C-
C-      call hbook1( 9003, 'M(X) (prod. kinem)                     $',nmsteps1, xma1, xmb1, 0.)
C-      call hbook1( 9020, 'M(X) generated  total                  $',nmsteps1, xma1, xmb1, 0.)
C-
C-      call hbook1( 9002, 'M(X) (dec.ph.sp.)x(waves)   fine bins  $',nmsteps2, xma2, xmb2, 0.)
C-      call hbook1( 9006, 'M(X) (dec.ph.sp.)     fine bins        $',nmsteps2, xma2, xmb2, 0.)
C-      call hbook1( 9009, 'M(X) (prod. kinem) fine bins           $',nmsteps2, xma2, xmb2, 0.)
C-
C-              call hbarx (9001)
C-              call hbarx (9005)
C-              call hbarx (9003)
C-              call hbarx (9020)
C-              call hbarx (9002)
C-              call hbarx (9006)
C-              call hbarx (9009)
C       
              if( ihist_decay .ne. 0) then
                  call ushist_react (ichannel)
              endif
       Lbook=.true.
       endif
CSdv-      
       write(*,*) ' ichannel, iwaves =', ichannel, iwaves
      
       if( ichannel .ne. 0 .and.  iwaves .ne. 0 ) then
c-----------------------------------------------------------
c      common / non_repeated_waves/ namp1, iw1(nampmax)
c      common /calc_ampl_names/ name_calc_ampl(nampmax)
        
       if (.not.Lbook2) then
          do i = 1, namp1
            do j = 1, namp1
              if(i .ne. j) then

              namtmp2 =  name_calc_ampl(iw1(i))(1:20)//' '//name_calc_ampl(iw1(j))(1:20)
              
c$$$              call hbook1( 10000+100*i+j, 'Real('//namtmp2(1:41)//')  $',nmsteps1 , xma1, xmb1,0.)
c$$$              call hbook1( 20000+100*i+j, 'Imag('//namtmp2(1:41)//')  $',nmsteps1 , xma1, xmb1,0.)
c$$$              call hbook1( 30000+100*i+j, 'Abs ('//namtmp2(1:41)//')  $',nmsteps1 , xma1, xmb1,0.)
c$$$              call hbook1( 40000+100*i+j, 'Phase('//namtmp2(1:41)//') $',nmsteps1 , xma1, xmb1,0.)
c$$$              call hbook1( 90000+100*i+j, 'Coher('//namtmp2(1:41)//') $',nmsteps1 , xma1, xmb1,0.)
c$$$              
c$$$              call hbook1( 310000+100*i+j, 'Real('//namtmp2(1:41)//')  $',nmsteps2 , xma2, xmb2,0.)
c$$$              call hbook1( 320000+100*i+j, 'Imag('//namtmp2(1:41)//')  $',nmsteps2 , xma2, xmb2,0.)
c$$$              call hbook1( 330000+100*i+j, 'Abs ('//namtmp2(1:41)//')  $',nmsteps2 , xma2, xmb2,0.)
c$$$              call hbook1( 340000+100*i+j, 'Phase('//namtmp2(1:41)//') $',nmsteps2 , xma2, xmb2,0.)
c$$$              call hbook1( 390000+100*i+j, 'Coher('//namtmp2(1:41)//') $',nmsteps2 , xma2, xmb2,0.)
                
              endif
            enddo
          enddo

          do i = 1, namp1

            namtmp2 =  name_calc_ampl(iw1(i))(1:20)

c$$$            call hbook1( 1000+i, namtmp2(1:40)//' $',nmsteps1 , xma1, xmb1,0.)
c$$$            call hbook1(      i, namtmp2(1:40)//' $',nmsteps1 , xma1, xmb1,0.)
c$$$            call hbook1( 2000+i, namtmp2(1:40)//' $',nmsteps2 , xma2, xmb2,0.)
c$$$
c$$$            call hbarx  (1000+i)
c$$$            call hbarx  (2000+i)

          enddo

      do ib = 1, ncoh
        do ia = 1, namp(ib)
          do ibw = 1, nbw(ib,ia)

             write(namtmp2,621)   name_calc_ampl(iw(ib,ia)), ibw
 621         format ( A40,' bw number ',I2)

c$$$c            call hbook1(9000000+10000*ib+100*ia+ibw , namtmp2(1:53)//' $',nmsteps1 , xma1, xmb1,0.)
c$$$
c$$$             call hbook1(1000000+10000*ib+100*ia+ibw , namtmp2(1:53)//' $',nmsteps1 , xma1, xmb1,0.)
c$$$             call hbook1(2000000+10000*ib+100*ia+ibw , namtmp2(1:53)//' $',nmsteps2 , xma2, xmb2,0.)
c$$$             call hbook1(3000000+10000*ib+100*ia+ibw , namtmp2(1:53)//' $',nmsteps2 , xma2, xmb2,0.)

           enddo
         enddo
       enddo
c                
c             call hbarx (0)

         Lbook2 = .true. 
       endif
    
c  density matrix in course binning --------------------------

CSdv+
              write(*,*) ' nmsteps1=', nmsteps1
 
              do ib1 = 1, nmsteps1
                 xmd = xma1 + xmstep1*(ib1-1) + 0.5*xmstep1

c------------- calculate and put to arrays ONLY big rho matrix, (economizing of time)

              call calc_sigma(xmd, weight_waves)

              do i = 1, namp1
                do j = 1, namp1
                    rho1= rho_big(iw1(i), iw1(j) )
c                  write(0,*) xmd, i,j, rho1
c                  if(i .ne. j) then
                    val1( ib1, i,j) = real(rho1)
                    val2( ib1, i,j) = aimag(rho1)
                    val3( ib1, i,j) = sqrt( real(rho1)**2 +
     *                                     aimag(rho1)**2 )
                  if(i .ne. j) then
                    val4( ib1, i,j) = atan2dmy(aimag(rho1),real(rho1) )
                    if( real(rho_big(iw1(i), iw1(i))) .gt. 0. .and.  
     *                  real(rho_big(iw1(j), iw1(j) )) .gt. 0. ) then
                      val9( ib1, i,j) =val3( ib1, i,j)/
     *  sqrt( real(rho_big(iw1(i), iw1(i))) *
     *        real(rho_big(iw1(j), iw1(j))) )
                    endif
                  endif
                enddo
              enddo

          enddo !imbin1 = 1, nmsteps1+1 

c----------------------------------
        endif                   !  iwaves.ne.0
c----------------------------------
      endif                     !  init0 .eq.0
  
  
c----------- only for fine bins  output   !!!!! -----------------
c   -------- expensive in CPU -----
CSdv+        
      write(*,*) 'init5, nmsteps2 =', init5, nmsteps2

      if( init5 .ne. 0 ) then  !   density matrix and intensities in "undistorted phase-space "  in fine bins 

c  density matrix in course binning --------------------------

              do ib1 = 1, nmsteps2
                xmd = xma2 + xmstep2*(ib1-1) + 0.5*xmstep2

c------------- calculate and put to arrays ONLY big rho matrix, (economizing of time)

              call calc_sigma(xmd, weight_waves)
	      
CSdv+
C-            write(*,*) ' xmd, weight_waves =', xmd, weight_waves  

              do i = 1, namp1
                do j = 1, namp1
                    rho1= rho_big(iw1(i), iw1(j) )
c                  write(0,*) xmd, i,j, rho1
c                  if(i .ne. j) then
                    val1a( ib1, i,j) = real(rho1)
                    val2a( ib1, i,j) = aimag(rho1)
                    val3a( ib1, i,j) = sqrt( real(rho1)**2 +
     *                                      aimag(rho1)**2 )
                  if(i .ne. j) then
                    val4a( ib1, i,j) = atan2dmy(aimag(rho1),real(rho1) )
                    if( real(rho_big(iw1(i), iw1(i))) .gt. 0. .and.  
     *                  real(rho_big(iw1(j), iw1(j) )) .gt. 0. ) then
                      val9a( ib1, i,j) =val3a( ib1, i,j)/
     *  sqrt( real(rho_big(iw1(i), iw1(i))) *
     *        real(rho_big(iw1(j), iw1(j))) )
                    endif
                  endif
                enddo
              enddo

          enddo ! imbin1 = 1, nmsteps1+1 

c--------------  use rndh for 3-4 particles phase-space -----------
               imodrn = 1
c------------------------------------------------------------------
               weight_gen = 1. ! if no channel defined

          if(  ichannel .ne. 0) then ! even if no waves defined...)

          do ib1 = 1, nmsteps2+1
c-------- here xmd  is fixed !!!!   -------------------------
            xmd = xma2 + xmstep2*(ib1-1) + 0.5*xmstep2
            
c------- integration of "undistorted cross-section" during init
c
c             nmc = 10000
c             nmc = 100
c             goto 888

            do imc = 1, nmc_fine

              call gener_calc_ampl_all(ichannel, iwaves,   xmd)

              if( iwaves.ne.0 ) call calc_sigma(xmd, weight_waves)
              
              xm = xmd

c----------------------------------------------
              if( iwaves.ne.0 ) then

              do i = 1, namp1
c                sintens(iw1(i)) = sintens(iw1(i)) + sampl2(iw1(i))
c$$$                 call hfill(2000+i, xm,  xtmp, sampl2(iw1(i))  )
c                if( i.eq. 1) write(0,*) ' xm = ', xm, ' sampl2 =', sampl2(iw1(i))

              enddo

cccccccccccccccccccccccccccccccc

      do ib = 1, ncoh
        do ia = 1, namp(ib)
          do ibw = 1, nbw(ib,ia)

            weight_tmp = sampl3(ib,ia,ibw)

c$$$            call hfill(2000000+10000*ib+100*ia+ibw, xm,  xtmp,   weight_tmp )
            
          enddo
        enddo
      enddo
      
cccccccccccccccccccccccccccccccc

c$$$               call hfill(9002, xm, xtmp, weight_waves  )

               endif
c--------------------------------------------------

c$$$               call hfill(9006, xm, xtmp, weight_gen  ) ! generator only

           enddo              ! nmc

              write(0,*) 'undistorted intensities for  xmd = ', xmd

              enddo ! mass bins

              write(0,*) 'undistorted intensities are prepared '

              write(0,*) 'start pure PROD kinematics at high mc-statistics '

c               
         k1 = 0
         k2 = 0
             
CSdv-    nmc_prod = 10 000 000
         nmc_prod =      1 000

         ntot_z = 0
         z1min  = 10000000.
         z1max  =-10000000.
         z2min  = 10000000.
         z2max  =-10000000.
         ztotmin= 10000000.
         ztotmax=-10000000.

         nquant = nmc_prod/100

         ifast_mode = 1

         write(0,*) ' ityper= ', ityper
         write(0,*) ' xma= ', xma, ' xmb= ', xmb

         do imc = 1, nmc_prod

         k1 =k1 + 1
         if( k1 .eq. nquant) then
            k1 = 0
            k2 = k2 + 1
          write(0,*)  imc,  ' events generated for fine PRODUCTION mass curve '
c         
          if( ityper .eq. 2) then
            write(0,*)  ' ntot_z= ', ntot_z
            write(0,*)  ' z1min = ', z1min, ' z2min = ',z2min 
            write(0,*)  ' z1max = ', z1max, ' z2max = ',z2max 
            write(0,*)  ' ztotmin = ', ztotmin, ' ztotmax = ', ztotmax
          endif
        endif
        
        if( ityper .eq. 1) call gener_dp2( xmd,  weightd     ) ! xmd is OUTPUT here !!!
        if( ityper .eq. 2) call gener_cp2( xmd,  weightd     ) ! xmd is OUTPUT here !!!
Csdv-   if( ityper .eq. 3) call gener_lp2( xmd,  weightd, 0  ) ! xmd is OUTPUT here !!!
CSdv-   if( ityper .eq. 5) call gener_lp2( xmd,  weightd, 1  ) ! xmd is OUTPUT here !!!
        if( ityper .eq. 4) call gener_dp1( xmd,  weightd     ) ! xmd is OUTPUT here !!!

        xm = xmd
               weight_tmp  = weightd

c              write(0,*) ' xm= ', xm, ' weight_tmp = ' , weight_tmp 

c              call hff1( 9009, n9009, xm,   weight_tmp  ) ! generator only
c$$$               call hfill(9009,  xm, xtmp,   weight_tmp  ) ! generator only

         enddo

         ifast_mode = 0

         write(0,*) 'end of  pure PROD kinematics at high mc-statistics '
c            

 888     continue

c=================================== norm of density - matrix ===========================

           if( iwaves.ne.0 ) then

             write(0,*) ' namp1 = ' , namp1

              do i = 1, namp1

c---------------- needed !! , production kinematics is taken into account ---

c$$$                CALL HOPERA (2000+i,'*',9009,  3000+i,1.,1.)

                write(0,*) ' before hunpak, ampl = ' ,i
c$$$                call hunpak(3000+i,  val0a( 1, i),' ',' ')
                write(0,*) ' after  hunpak, ampl = ' ,i
              enddo

CSdv+:   actually ends here   !!!!

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do ib = 1, ncoh
        do ia = 1, namp(ib)
          do ibw = 1, nbw(ib,ia)
c$$$                CALL HOPERA (2000000+10000*ib+100*ia+ibw,'*',9009, 3000000+10000*ib+100*ia+ibw  ,1.,1.)          
          enddo
        enddo
      enddo
      
cccccccccccccccccccccccccccccccc

            do ib1 = 1, nmsteps2

                xm1 = xma2 + xmstep2*(ib1-1) + 0.5*xmstep2 

c               write(0,*) ' ib1= ', ib1, ' xm1= ' , xm1

            do i = 1, namp1
c                  write(0,*) ' i= ', i 
c                 write(0,*) ' val1= ',val1( ib1, i,i)
c                 write(0,*)  ' val0= ',val0( ib1, i)
              do j = 1, namp1
                if(i .ne. j) then

                  denom = val1a( ib1, i,i) *  val1a( ib1, j,j)

                  if( denom .ne. 0.) then 

                  renorm =  sqrt(  val0a(ib1,i) * val0a(ib1,j) / denom ) 

                    val1a( ib1, i,j) = val1a( ib1, i,j)*renorm
                    val2a( ib1, i,j) = val2a( ib1, i,j)*renorm
                    val3a( ib1, i,j) = val3a( ib1, i,j)*renorm

c                   write(0,*) ' ib1= ', ib1, ' xm1= ' , xm1, ' val1= ',val1( ib1, i,j)
c                   write(0,*) ' val0= ',  val0(ib1,i), val0(ib1,j), ' renorm= ', renorm 

                    endif
                endif
              enddo
            enddo
            
            enddo ! nmsteps1

            write(0,*) ' before pak'
              do i = 1, namp1
                do j = 1, namp1
                  if(i .ne. j) then
c$$$                    call hpak(310000+100*i+j, val1a( 1, i,j) )
c$$$                    call hpak(320000+100*i+j, val2a( 1, i,j) )
c$$$                    call hpak(330000+100*i+j, val3a( 1, i,j) )
c$$$                    call hpak(340000+100*i+j, val4a( 1, i,j) )
c$$$                    call hpak(390000+100*i+j, val9a( 1, i,j) )
                  endif
                enddo
              enddo

              endif ! iwaves

c======================================================================

           endif ! ichannel. ne. 0
c----------------------------------
      endif ! init5 .ne. 0

c=========== end of initialisation phase ==============================
c
c== for each input CP kinematics : 

      if( init1 .ne. 0 ) then  ! usual generating mode, including over-all CP kinematics
c--------------  use pp2ranlux  for 3-4 particles phase-space -----------
               imodrn = 0
c              imodrn = 1
c=================  full cross-section  here ..., find maximal weight  =========================

c            do imc = 1, nmc

c       write(0,*) ' before gener_xp, ityper = ', ityper
   
        if( ityper .eq. 1) call gener_dp2( xmd,  weightd     ) ! xmd is OUTPUT here !!!
        if( ityper .eq. 2) call gener_cp2( xmd,  weightd     ) ! xmd is OUTPUT here !!!
CSdv-   if( ityper .eq. 3) call gener_lp2( xmd,  weightd, 0  ) ! xmd is OUTPUT here !!!
CSdv-   if( ityper .eq. 5) call gener_lp2( xmd,  weightd, 1  ) ! xmd is OUTPUT here !!!
        if( ityper .eq. 4) call gener_dp1( xmd,  weightd     ) ! xmd is OUTPUT here !!!
CSdv-
C        write(0,*) ' weightd = ', weightd, ' xmd= ', xmd

              weight_gen = 1 ! if no channel defined

cc            write(0,*) ' ichannel =', ichannel, ' iwaves= ', iwaves

              call gener_calc_ampl_all(ichannel, iwaves, xmd)

              if( iwaves.ne.0 ) call calc_sigma(xmd, weight_waves)

              write(*,*) ' weight_waves= ', weight_waves
              if( iwaves.ne.0 ) then 
                weight_all =  weightd  *  weight_waves
              else
                weight_all =  weightd  * weight_gen
              endif
c             write(0,*) ' weight_all= ', weight_all
      
c-----------------------------------------------
      endif ! init1 .ne. 0

c=======================================================================
      if(  init2 .ne. 0   ) then  !   fill sintensities just from MC 
c=========== fill histograms with "sintensities" ===================

        xm = xmd
        if( ichannel .ne. 0 .and. iwaves.ne.0 ) then
          do i = 1, namp1
c           sintens(iw1(i)) = sintens(iw1(i)) + sampl2(iw1(i))
            weight_tmp =  sampl2(iw1(i))*weightd
c$$$            call hfill(1000+i, xm,  xtmp,   weight_tmp )
          enddo

c--------------------------------
      do ib = 1, ncoh
        do ia = 1, namp(ib)
          do ibw = 1, nbw(ib,ia)

            weight_tmp =    sampl3(ib,ia,ibw)  *weightd

c$$$            call hfill(1000000+10000*ib+100*ia+ibw, xm,  xtmp,   weight_tmp )
            
          enddo
        enddo
      enddo
      
c--------------------------------
        endif

        weight_tmp = weight_all
c$$$        call hfill(9001, xm,  xtmp, weight_tmp  )

        weight_tmp =  weightd * weight_gen
c$$$        call hfill(9005, xm,  xtmp, weight_tmp  )
        
        weight_tmp = weightd
c$$$        call hfill(9003 , xm ,xdum , weight_tmp)
                    
c-      write(0,*) ' wei= ', weight
c       
c       ssigma(imbin) = ssigma(imbin) + weight
c       
c       if(  weight .gt. weight_max(imbin) ) 
c     *      weight_max(imbin) = weight
c                   
c-           enddo                    ! nmc

c      xmd = dsqrt( pres(4)**2 - pres(1)**2 - pres(2)**2 - pres(3)**2 )
c
c       imbin = (xmd - xma)/xmstep + 1
c
c       xxx = ssigma_max * pp2rnd(-1)
c
c----------------------------------
      endif ! init2 .ne. 0

c       if( xxx .le.  ssigma(imbin) ) then
c        ireject = 0
c       else
c        ireject = 1
c        return
c       endif

c----------- in case ireject .eq. 0 calculate decay, weighted by differential cross-section ---

c 3     continue

c       if( ichannel .eq.  3) then
c         call gener_2eta(xmd)
c         call calc_ampl_2eta
c       endif
cc--------------------------------------------------
c       if( ichannel .eq.  4) then
c         call gener_4pic(xmd)
c         call calc_ampl_4pic
c       endif
c
c       
c       call calc_sigma(xmd, weight)
c
c       xxx = weight_max(imbin) * pp2rnd(-1)
c
c       if( xxx .gt. weight  ) goto 3

c---------------------------------------------------------------------------
      if(  init3 .ne. 0) then  !   kinematics for outputting this event ...

c--------- finally, make appropriate rotations and boost into lab system here... ---


        if( ityper .eq. 1 .or. ityper .eq. 4) then
       do k = 1,4
         pbeam(k) =  PLAB(k,1)
         pslow(k) =  PLAB(k,4)
         pres(k) =  PLAB(k,7)
         pbeam1(k) = pbeam(k) 
       enddo
       endif

       if( ityper .eq. 2 .or. ityper .eq. 3 .or. ityper .eq. 5 ) then
       do k = 1,4
         pbeam(k) =  PLAB(k,1)
         pfast(k) =  PLAB(k,3)
         pslow(k) =  PLAB(k,4)
         pres (k) =  PLAB(k,7)
         pbeam1(k) = pbeam(k) - pfast(k) 
       enddo
       endif

       wmtarg = dsqrt( pslow(4)**2 - 
     *                 pslow(1)**2 - pslow(2)**2 - pslow(3)**2)

       CALL  BOOST2( xmd, pres(1), pbeam1(1),pbeam1_cm(1) )
       CALL  BOOST2( xmd, pres(1), pslow(1),pslow_cm(1) )
c      
       if( ityper .eq. 1 .or.  ityper .eq. 2 .or. ityper .eq. 4) then
         call GJrotM(pbeam1_cm(1), pslow_cm(1), anGJ(1,1))
       endif
       if( ityper .eq. 3 .or.  ityper .eq. 5 ) then
         call HSrotM(pbeam1_cm(1), pslow_cm(1), anGJ(1,1))
       endif    
       
c        nmes = 2
c        write(0,*) ' ityper = ', ityper
c        write(0,*) ' pdec (1) = ', pdec(1,1), pdec(2,1), pdec(3,1)
c        write(0,*) ' pdec (2) = ', pdec(1,2), pdec(2,2), pdec(3,2)
c        write(0,*) ' pdec (3) = ', pdec(1,3), pdec(2,3), pdec(3,3)
c        write(0,*) ' pdec (4) = ', pdec(1,4), pdec(2,4), pdec(3,4)
  
      do i=1,nmes ! if nmes = 0, do nothing ...

c-----   GJ-system-------CM-system-------	
         do j = 1,3
           pdec1(j,i) = 0.
           do k = 1,3
            pdec1(j,i) = pdec1(j,i) + anGJ(k,j)*pdec(k,i)
           enddo
         enddo
         pdec1(4,i) = pdec(4,i)
c------------------------------
         
         call BOOST1( xmd , pres(1) , pdec1(1,i) , pdec_lab(1,i ) ) ! LAB system
         
      enddo

c----------------------------------
      endif ! init3 .ne. 0

c=======================================================================
      if(  init4 .ne. 0) then  !   hopera and hrput ...

           if( ichannel .ne. 0 .and. iwaves.ne.0 ) then

             write(0,*) ' namp1 = ' , namp1

              do i = 1, namp1
                write(0,*) ' before hunpak, ampl = ' ,i
c               call hunpak(1000+i,  val_tmp(1) )
c$$$                call hunpak(1000+i,  val0( 1, i),' ',' ')
                write(0,*) ' after hunpak, ampl = ' ,i
              enddo

              do ib1 = 1, nmsteps1

                xm1 = xma1 + xmstep1*(ib1-1) + 0.5*xmstep1 

c               write(0,*) ' ib1= ', ib1, ' xm1= ' , xm1

              do i = 1, namp1
c                 write(0,*) ' i= ', i 
c                 write(0,*) ' val1= ',val1( ib1, i,i)
c                 write(0,*) ' val0= ',val0( ib1, i)
c
              do j = 1, namp1
                 if(i .ne. j) then

                  denom = val1( ib1, i,i) *  val1( ib1, j,j)

                  if( denom .ne. 0.) then 

                  renorm =  sqrt(  val0(ib1,i) * val0(ib1,j) / denom ) 

                  val1( ib1, i,j) = val1( ib1, i,j)*renorm
                  val2( ib1, i,j) = val2( ib1, i,j)*renorm
                  val3( ib1, i,j) = val3( ib1, i,j)*renorm

c                 write(0,*) ' ib1= ', ib1, ' xm1= ' , xm1, ' val1= ',val1( ib1, i,j)
c                 write(0,*) ' val0= ',  val0(ib1,i), val0(ib1,j), ' renorm= ', renorm 

                  endif
                endif
              enddo
            enddo
            
            enddo ! nmsteps1

                write(0,*) ' before pak'
              do i = 1, namp1
                do j = 1, namp1
                  if(i .ne. j) then
c$$$                  call hpak(10000+100*i+j, val1( 1, i,j) )
c$$$                  call hpak(20000+100*i+j, val2( 1, i,j) )
c$$$                  call hpak(30000+100*i+j, val3( 1, i,j) )
c$$$                  call hpak(40000+100*i+j, val4( 1, i,j) )
c$$$                  call hpak(90000+100*i+j, val9( 1, i,j) )
                  endif
                enddo
              enddo


c              do i = 1, namp1
c                call hpak(1000+i,  val0( 1, i))
c              enddo
c
c                call hpak(9001,  valtot(1) )

              write(0,*) ' before hopera'

         sum_tot1 = HSUM (9001)

         scale_1  = nmc_norm/sum_tot1

c--------------- 9001 is finally rescaled to nmc_norm  ------
c$$$         CALL HOPERA (9001,'+',9001, 9001, scale_1  ,0.)
	 
         sum_tot =   HSUM (9001)

c$$$         CALL HOPERA (9002,'*',9009,9004,1.,1.)

         sum_tot2 =   HSUM (9004)

         write(0,*) ' sum_tot2 = ', sum_tot2

         scale_2 = nmc_norm/sum_tot2*(xmstep1/xmstep2)


c--------------- 9001 is finally rescaled to nmc_norm  ------
c$$$         CALL HOPERA (9004,'+',9009, 9004, scale_2  ,0.)

         write(0,*) ' sum_tot = ', sum_tot
         write(0,*) '  scale_1  = ', scale_1

       do i = 1, namp1
c        CALL HOPERA (9002,'/',9001,8888,1.,1.)
c        CALL HOPERA (2000+i,'*',,2000+i,1.,1.)

c$$$         CALL HOPERA (3000+i,'+',9009,  3000+i, scale_2 ,0.)

c        CALL HOPERA (1000+i,'/',9001,999000+i,1.,1.)
c        CALL HOPERA (999000+i,'*',9020,8000+i,1.,1.)

c$$$         CALL HDELET (i)     ! for the mass-dep fit
c$$$
c$$$         CALL HOPERA (1000+i,'+',9020, i, scale_1  ,0.)
c$$$
c$$$         sum =  HSUM (i)

         write (0,651) i, name_calc_ampl(iw1(i)), sum,  sum/sum_tot*100. 
 651     format('Intens (',I2,' )', A40, ' sum = ',F14.3,'   frac= ', F9.6 ,' %'   )

       enddo

       do ib = 1, ncoh
         do ia = 1, namp(ib)
           do ibw = 1, nbw(ib,ia)
            
c        CALL HOPERA (900000+10000*ib+100*ia+ibw  ,'/',9001,9999999 ,1.,1.)
c        CALL HOPERA (9999999,'*',9020, 100000+10000*ib+100*ia+ibw  ,1.,1.) 

c$$$         CALL HOPERA ( 1000000+10000*ib+100*ia+ibw,'+',9020, 1000000+10000*ib+100*ia+ibw  , scale_1  ,0.)
c$$$         sum =  HSUM ( 1000000+10000*ib+100*ia+ibw )

         write (0,652) ib, ia,  name_calc_ampl(iw(ib,ia)), ibw,  sum,   sum/sum_tot*100. 
 652     format(' Nblock (',I2,'), Namp (',I2,' )', A40,' Nbw(' , I2, ') sum= ', F14.3, '   frac= ', F9.6 ,' %')

c        write (0,*) ' N_in_bl (',ib,') ia =  (',ia,') ibw(',ibw,') ' , sum,' frac= ', sum/sum_tot*100., ' % '           '

c$$$         CALL HOPERA ( 3000000+10000*ib+100*ia+ibw,'+',9009, 3000000+10000*ib+100*ia+ibw   , scale_2  ,0.)
    
           enddo
         enddo
       enddo
      
c  scale  non-diagonal elements to nmc_norm 

         do i = 1, namp1
           do j = 1, namp1
             if(i .ne. j) then

c$$$             CALL HOPERA (10000+100*i+j  ,'+', 9020, 10000+100*i+j , scale_1  ,0.)
c$$$             CALL HOPERA (20000+100*i+j  ,'+', 9020, 20000+100*i+j , scale_1  ,0.)
c$$$             CALL HOPERA (30000+100*i+j  ,'+', 9020, 30000+100*i+j , scale_1  ,0.)
c$$$
c$$$             CALL HOPERA (310000+100*i+j ,'+', 9009, 310000+100*i+j, scale_2  ,0.)
c$$$             CALL HOPERA (320000+100*i+j ,'+', 9009, 320000+100*i+j, scale_2  ,0.)
c$$$             CALL HOPERA (330000+100*i+j ,'+', 9009, 330000+100*i+j, scale_2  ,0.)

             endif
           enddo
         enddo
            
         endif ! iwaves.ne.0

         write(0,*)   'before hrput'
c$$$         CALL HRPUT(0,'test_cp.hbook4','N')
         write(0,*)   'after hrput'

c----------------------------------
       endif ! init4 .ne. 0

c------that's all -----------------------------
       
      return
      end

     
c------------ calculate differential cross-section weight -----------
c
      subroutine calc_sigma(xmd, weight_waves )
      real *8 xmd

c-------- cross section structure integer and floating parameters, defined from cards
      parameter  (ncohmax= 20, nbwmax=10, nparbwmax=15 )
      parameter  (nampmax=100 )
      common /cross_sect_struct_int/ ncoh, namp(ncohmax),
     *      nbw( ncohmax ,nampmax),  iw(ncohmax,nampmax)
      common /cross_sect_struct_float/
     *      parbw( nparbwmax, ncohmax, nampmax, nbwmax ),
     *      cn(ncohmax, nampmax, nbwmax) 
      complex  cn

      common / non_repeated_waves/ namp1, iw1(nampmax)

c-------- weight for n-meson phase space ---------
      common /gener_weight/ weight_gen

c----------- amplitudes, names,   quantum numbers from calc_ampl (physics-dependent)
      common /calc_ampl_amplits/ ampl(nampmax)
      complex  ampl

      common /calc_ampl_names/ name_calc_ampl(nampmax)
      character *60  name_calc_ampl
c------------------------------------ total spin --- spin projection -- reflectivity(if def.) --
      common /calc_ampl_quantum_numb/  J1 (nampmax),    JM (nampmax), IP(nampmax),  IETA(nampmax)  

      common /calc_ampl_totals/ namp_tot


c------------------ stored during initialization  integrals and max weights
      parameter  ( nmasspnt = 3000 )
      common /samplits/  samp( nampmax, nmasspnt )


c------ local variables , functions
      complex ccoh(ncohmax)
      complex ccoh_a(ncohmax,nampmax )
      complex ccoh_intens(ncohmax,nampmax )
      complex cbw, c_tot

      common /ampl2_contrib/ sampl2(nampmax)
      real sampl2

      common /ampl3_contrib/ sampl3(ncohmax, nampmax, nbwmax )
      real sampl3

c      common / total_intensities/  sintens(nampmax)
c      real sintens 

      common / total_density_matrix/  rho_big(nampmax, nampmax )
      complex rho_big

      common /decay_channel/ ichannel

      real *8 par_all(50)
      complex *16  bwuniv(10)
      integer it_bw
CSdv+       
      data namp1  /  2 /    ! common /cross_sect_struct_int/ ncoh, namp(ncohmax), ...
      data ncoh   /  2 /    ! common / non_repeated_waves/ namp1, iw1(nampmax)
C         
CSdv- weight_waves = 0.
      weight_waves = 1.
      return                !   Sdv+ !!!!!!!!!!!!!!!!!!!!

      do i = 1, namp1
        sampl2(iw1(i)) = 0.
      enddo

      do i = 1, namp1
        do j = 1, namp1
          rho_big(iw1(i), iw1(j) ) = (0.,0.)
        enddo
      enddo

CSdv+
      write(*,*) ' calc_sigma: namp1, ncoh = ', namp1, ncoh
      write(*,*) ' calc_sigma: namp =', namp
      write(*,*) ' calc_sigma: nbw  =', nbw

      do ib = 1, ncoh
        ccoh(ib) = cmplx(0.,0.)
        do ia = 1, namp(ib)
          ccoh_intens(ib,iw(ib,ia) ) = cmplx(0.,0.)
          ccoh_a(ib,iw(ib,ia) ) = cmplx(0.,0.) 
          do ibw = 1, nbw(ib,ia)

c           ccoh(ib) =  ccoh(ib) + cn(ib,ia,ibw) * ampl( iw(ib,ia) ) / sqrt( samp( iw(ib,ia), imbin ) ) *

             it_bw = parbw(1,ib,ia,ibw)
             do ip2 = 1, nparbwmax-1 ! 10 < 50 !!
               par_all(ip2) =  parbw(ip2+1, ib,ia,ibw)
             enddo

            call bwuniv_all( xmd, it_bw,  par_all(1),  bwuniv(1) )  ! new, standartized !!!

c             write(0,*) ' par_bw= ', (parbw(ll, ib, ia, ibw), ll=1, 15)
c             write(0,*) ' par_all= ', (par_all(ll), ll=1, 15)
c            write(0,*) ' xmd = ', xmd,' ibw=', ibw, ' it_bw=', it_bw,           'bw = ', bwuniv(1) 

             c_tot = cn(ib,ia,ibw) * ampl( iw(ib,ia) )/sqrt(abs(  psp_tab(xmd, iw(ib,ia)) )) * bwuniv(1) ! no K-matrix implemented yet...
c            c_tot = cn(ib,ia,ibw) * ampl( iw(ib,ia) ) * cbw( xmd, parbw(1,ib,ia,ibw) )

             sampl3(ib,ia,ibw) = ( real(c_tot)**2 + aimag(c_tot)**2 ) *  weight_gen

            ccoh(ib) =  ccoh(ib) +  c_tot
c     cn(ib,ia,ibw) * ampl( iw(ib,ia) ) * cbw( xmd, parbw(1,ib,ia,ibw) ) 
            
            ccoh_intens(ib,iw(ib,ia) ) =  ccoh_intens(ib,iw(ib,ia) ) + c_tot
 
c     * cn(ib,ia,ibw) * ampl( iw(ib,ia) ) *cbw( xmd, parbw(1,ib,ia,ibw))
            
            ccoh_a(ib,iw(ib,ia) ) =  ccoh_a(ib,iw(ib,ia) ) +
c     * cn(ib,ia,ibw) *  cbw( xmd, parbw(1,ib,ia,ibw) ) 
     * cn(ib,ia,ibw) *  bwuniv(1)
               
          enddo
          
          sampl2(iw(ib,ia)) = sampl2(iw(ib,ia)) +
     * conjg( ccoh_intens(ib,iw(ib,ia) ) )* ccoh_intens(ib,iw(ib,ia) ) *
     *  weight_gen


        enddo !ia = 1, namp(ib) 

        weight_waves = weight_waves +  conjg( ccoh(ib))*  ccoh(ib)

c-   fill  "big density matrix" here ...
        do ia1 = 1, namp(ib)
          do ia2 = 1, namp(ib)
            
            rho_big(iw(ib,ia1), iw(ib,ia2) ) =  
     *      rho_big(iw(ib,ia1), iw(ib,ia2) ) + 
     *  conjg( ccoh_a(ib,iw(ib,ia1) ))*  ccoh_a(ib,iw(ib,ia2))
            
          enddo
        enddo
c-    
      enddo
            
      weight_waves = weight_waves * weight_gen 
      
      return
      end
      
C-------------------------------------------------------------------------------------------------
      subroutine  bwuniv_all( xm, it, par_all,  bwuniv )
      integer it
      real *8 xm
      real *8 par_all(50)
      complex *16  bwuniv(10)

      complex*16  BWSIMP_D, BW1L_D,  BW2L_D,  BGEXP_D,  BGEXP_POW_D, BWA1_BOWLER_D

      !print *, "it = ", it

      if( it.eq. -1) then
        bwuniv(1) = 1.
      endif

      if( it.eq. 0) then
          bwuniv(1) = BWSIMP_D(xm , par_all(1) )
      endif

      if( it.eq. 1) then
          bwuniv(1) = BW1L_D(xm , par_all(1) )
      endif

      if( it.eq. 2) then
          bwuniv(1) = BW2L_D(xm , par_all(1) )
      endif

      if( it.eq. 3) then
          bwuniv(1) = BGEXP_D(xm , par_all(1) )
      endif

      if( it.eq. 4) then
          bwuniv(1) = BGEXP_POW_D(xm , par_all(1) )
      endif

      if( it.eq. 5) then
          bwuniv(1) = BWA1_BOWLER_D(xm , par_all(1) )
      endif

      if( it.eq. 6) then
          call BWA1_BOWLER_DECK_D(xm , par_all(1),  bwuniv(1) )
      endif

      if( it.eq. 7) then
          call BWGC_DECK_D(xm , par_all(1),  bwuniv(1) )
      endif

      if (it.eq.71) then
         call BWconformal(xm, par_all(1), bwuniv(1))
      end if

      if (it.eq.72) then
         call BWconformal_a2(xm, par_all(1), bwuniv(1))
      end if

      if (it.eq.73) then
         call BWgaus(xm, par_all(1), bwuniv(1))
      end if

      if (it.eq.74) then
         call BWconformalPoly(xm, par_all(1), bwuniv(1))
      end if

      if (it.eq.75) then
         call BWa2_a2prime_poly(xm, par_all(1), bwuniv(1))
      end if

      if (it.eq.76) then
         call BWa2_a2prime_poly_expo(xm, par_all(1), bwuniv(1))
      end if
      return
      end

C-----------------------------------------------
        COMPLEX *16  FUNCTION BWSIMP_D(AM, par_all)
        real *8  par_all(*) 
        real *8 am, am1,g1, s,a,b,c,den,bwre,bwim

	S = am*am
        am1 = par_all(1)
        g1  = par_all(2)

        A=AM1**2-S
        B=AM1*G1

c       C = AM1*G1

        C = dsqrt(AM1*G1)

        DEN = A*A + B*B
        BWRE= C*A/DEN
        BWIM= C*B/DEN

	BWSIMP_D =DCMPLX(BWRE,BWIM)

        RETURN
        END
	
C-------------------------------------------------------------------------------------------
c       complex FUNCTION BW1L(AM, AM0, G0, AMA, AMB, R, L) 
        complex *16 FUNCTION BW1L_D(AM, par_all) 
        real *8  par_all(50) 

        real *8  am,am0d, s, ama, amb,  r,  am0, g0, g

        real *8 a, b, c, den, bwre, bwim

        real *8 psl1_d

        ama = par_all(1)
        amb = par_all(2)
        r   = par_all(3)
        l   = par_all(4)
        am0 = par_all(5)
        g0  = par_all(6)

	S = AM*AM

        am0d = am0

        G = G0 * AM0 / AM * psl1_d ( am, ama, amb,  R, L ) / psl1_d ( am0d, ama, amb, R, L )

        A = AM0**2-S
        B = AM0*G
        C = DSQRT(AM0*G0)
	DEN = A*A + B*B
        BWRE= C*A/DEN
        BWIM= C*B/DEN
	
        BW1L_D=DCMPLX(BWRE,BWIM)

c       write(0,*) ' ama = ' , ama, ' amb = ', amb,' r=',r,' l=', l
c       write(0,*) ' am0, g0 = ', am0, g0
c       write(0,*) ' am= ', am, 'BW1L_D  = ',	BW1L_D

        RETURN
        END

C-----------------------------------------------------------------------------------------------------
        real *8   function   PSL1_D ( am,  am1, am2,  R, lorb  )

        real *8 am
        real *8 R
        real *8 am1,am2
	real *8 s,e
        real *8 p

c	real *8   fd1_d,fd2_d,fd3_d,fd4_d

        real *8 fdl_d

	ampor = am1 + am2
        S     = AM*AM

	if( am .gt. ampor ) then

	    E = ( S + am1**2 - am2**2 )/(2.*am)

	    P = dsqrt( dabs(E**2 - am1**2) ) 

            PSL1_D = p**(2*lorb+1)/FDL_D(p, r, lorb)

 	    else

	 PSL1_D = 0.

	endif

        RETURN
        END

c       complex FUNCTION BW2L(AM, AM0, G0, AM1A,AM1B,R1,L1, AM2A,AM2B,R2,L2,X2) 
        complex *16  FUNCTION BW2L_D( AM,  par_all) 
        real *8  par_all(50) 
          real *8  am,am0d, s, am1a, am1b, am2a, am2b, r1, r2, x2, am0, g0, g

          real *8 a, b, c, den, bwre, bwim

          real *8 psl1_d

          am1a = par_all(1)
          am1b = par_all(2)
          r1   = par_all(3)
          l1   = par_all(4)
          am2a = par_all(5)
          am2b = par_all(6)
          r2   = par_all(7)
          l2   = par_all(8)
          x2   = par_all(9)

          am0  = par_all(10)
          g0   = par_all(11)

	  S    = AM*AM
          am0d = am0

          G = G0 * AM0 / AM * ( 
     *   (1. - x2) * psl1_d ( am, am1a, am1b,  R1, L1 ) / psl1_d ( am0d, am1a, am1b, R1, L1 ) +
     *         x2  * psl1_d ( am, am2a, am2b , R2, L2 ) / psl1_d ( am0d, am2a, am2b, R2, L2 )
     *              )
        A=AM0**2-S
        B=AM0*G
        C=SQRT(AM0*G0)
	DEN = A*A + B*B
        BWRE=C*A/DEN
        BWIM=C*B/DEN
        BW2L_D=DCMPLX(BWRE,BWIM)
	
c         write(0,*) ' am1a = ' , am1a, ' am1b = ', am1b,' r1=',r1,' l1=', l1
c         write(0,*) ' am2a = ' , am2a, ' am2b = ', am2b,' r2=',r2,' l2=', l2
c         write(0,*) ' x2 = ', x2
c         write(0,*) ' am0, g0 = ', am0, g0
c         write(0,*) ' am= ', am, 'BW2L_D  = ',   BW2L_D

        RETURN
        END
C-------------------------------------------------------------

        real*8  function fdl_d( p, r, lorb)
        real*8  p, r, pr
	real*8  fd1_d,fd2_d,fd3_d,fd4_d

	pr = P*R

        if( lorb.eq.0) FDL_D  = 1.
        if( lorb.eq.1) FDL_D  = FD1_D(PR)
        if( lorb.eq.2) FDL_D  = FD2_D(PR)
        if( lorb.eq.3) FDL_D  = FD3_D(PR)
        if( lorb.eq.4) FDL_D  = FD4_D(PR)

        RETURN
        END
C------------------------------------
 	real *8 FUNCTION  FD1_D(X)
        real *8 X
	FD1_D = 1. + X*X
	RETURN
	END
C------------------------------------
	real *8 FUNCTION  FD2_D(X)
        real *8 X
	FD2_D = 9. + 3.*X*X + X*X*X*X
	RETURN
	END
C------------------------------------
	real *8 FUNCTION  FD3_D(X)
        real *8 X
	FD3_D = 225. + 45.*X*X + 6.*X*X*X*X + X*X*X*X*X*X
	RETURN
	END
C--------------------------------
        real *8 FUNCTION  FD4_D(X)
        real *8 X
        FD4_D = 11025. + 1575.*X*X + 135.*X*X*X*X + 10.*X*X*X*X*X*X + X*X*X*X*X*X*X*X
        RETURN
        END
C--------------------------------
C-------------------------------------------------------------------------------------------
c       complex FUNCTION BWA1_BOWLER(AM, AM0,G0) 
        complex  *16 FUNCTION BWA1_BOWLER_D(AM, par_all) 
        real *8  par_all(50) 

        real *8  am, am0, g0, s, g, a, b, c, den, bwre, bwim
        real *8  pav_rhopi_s_d

        am0 = par_all(1)
        g0  = par_all(2)
	S   = AM*AM
        A   = AM0**2-S

        G = G0 * pav_rhopi_s_d(am)/ pav_rhopi_s_d(am0)  * am0/am

        B = AM0 * G
        C = DSQRT(AM0*G0)
	DEN = A*A + B*B
        BWRE= C*A/DEN
        BWIM= C*B/DEN
        BWA1_BOWLER_D=DCMPLX(BWRE,BWIM)

        RETURN
        END

c---------------------------------------------------------------------

        subroutine BWA1_BOWLER_DECK_D(am , par_all,  bwuniv )
        real *8  par_all(50) 
        complex *16   bwuniv(10)

        real *8  am, am0, g0, s, g, a, b, c, den, bwre, bwim
        real *8  pav_rhopi_s_d
        complex *16 bwa1
        real *8 deck_simp
        parameter (ampi = 0.13957018)

        am0 = par_all(1)
        g0 =  par_all(2)

	S = AM*AM

        A=AM0**2-S

        G = G0 * pav_rhopi_s_d(am)/ pav_rhopi_s_d(am0)  * am0/am
C-      G = G0 * pav_rhopi_s_d(am)/ pav_rhopi_s_d(am0)

        B = AM0 * G
        C = DSQRT(AM0*G0)
	DEN  = A*A + B*B
        BWRE = C*A/DEN
        BWIM = C*B/DEN
        BWA1 = DCMPLX(BWRE,BWIM)
C------------------------------------------

        deck_simp = 1./( am**2 - ampi**2 )

        bwuniv(1) = bwa1
        bwuniv(2) = bwa1 * deck_simp * ( am0**2 - am**2 +  am0 * g/pav_rhopi_s_d(am) *   par_all(3) )
c       bwuniv(2) = bwa1 * deck_simp * am * ( am0**2 - am**2 +  am0 * g/pav_rhopi_s_d(am) *   par_all(3) )

        RETURN
        END

C---------------------------------------------

        subroutine BWGC_DECK_D(am , par_all,  bwuniv )
        real *8  par_all(50) 
        complex *16   bwuniv(10)

        real *8  am, am0, g0, s, g, a, b, c, den, bwre, bwim
        real *8  pav_rhopi_s_d
        complex *16 bwa1
        real *8 deck_simp
        parameter (ampi = 0.13957018)

        am0 = par_all(1)
        g0  = par_all(2)
	S   = AM*AM
        A   = AM0**2-S
        G   = G0

        B   = AM0 * G
        C   = DSQRT(AM0*G0)
	DEN = A*A + B*B
        BWRE= C*A/DEN
        BWIM= C*B/DEN
        BWA1=DCMPLX(BWRE,BWIM)
c-----------------------------

        deck_simp = 1./( am**2 - ampi**2   )

c------------   0.3 is arbitrary break-up momentum at resonance , anyway modified by  par_all(3)

        bwuniv(1) = bwa1
        bwuniv(2) = bwa1 * deck_simp * ( am0**2 - am**2 +  am0 * g / 0.3 *   par_all(3) )

        RETURN
        END

      subroutine BWconformal_a2(am, par, bwuniv)
      implicit none
      real*8 am, par(50)
      complex*16 bwuniv(10)

      ! Breit-Wigner deformed by a polynomial in the conformal variable
      ! omega given below

      complex*16 bw
      complex*16 omega
      complex*16  BWSIMP_D, BW1L_D,  BW2L_D,  BGEXP_D,  BGEXP_POW_D, BWA1_BOWLER_D

      real*8 am1a, am1b, r1, l1, am2a, am2b, r2, l2, x2
      real*8 parbw(11)
      real*8 linRe, linIm, quadRe, quadIm

      real*8 mPi, mEta, s0
      parameter(mPi = 0.139, mEta = 0.548, s0 = (mPi+mEta)**2)

      ! 9 parameters:
      !   as in BW2L_D
      ! 6 variables:
      !   BW mass
      !   BW width
      !   real part of linear term
      !   imag part of linear term
      !   real part of quadratic term
      !   imag part of quadratic term

      bw = BW2L_D(am, par)

      linRe = par(12)
      linIm = par(13)
      quadRe = par(14)
      quadIm = par(15)

      !print "(6f7.3)", par(1:6)

      omega = (am - complex(0d0,1d0)*sqrt(am*am - s0)) / (am + complex(0d0,1d0)*sqrt(am*am - s0))
      bwuniv(1) = bw * (1 + omega*cmplx(linRe,linIm,8) + omega**2*cmplx(quadRe,quadIm,8))

      !print "(2f7.3)", bwuniv(1)
      end

      SUBROUTINE BWGAUS(am, par, bwuniv)
      real*8 am, par(50)
      complex*16 bwuniv(10)

      real*8 mean, sigma

      mean = par(1)
      sigma = par(2)

      bwuniv(1) = sqrt(exp(-0.5*(am - mean)**2 / sigma**2))
      END

      subroutine BWconformalPoly(am, par, bwuniv)
      real*8 am, par(50)
      complex*16 bwuniv(10)

      ! A polynomial in the conformal variable omega given below
      ! 3 parameters:
      !   imag part of linear term
      !   real part of quadratic term
      !   imag part of quadratic term

      complex*16 bw
      complex*16 omega
      complex*16  BWSIMP_D, BW1L_D,  BW2L_D,  BGEXP_D,  BGEXP_POW_D, BWA1_BOWLER_D

      real*8 parbw(11)
      real*8 const, linRe, linIm, quadRe, quadIm

      real*8 mPi, mEta, s0
      parameter(mPi = 0.139, mEta = 0.548, s0 = (mPi+mEta)**2)

      const = par(1)
      linRe = par(2)
      linIm = par(3)
      quadRe = par(4)
      quadIm = par(5)

      !print "(6f7.3)", par(1:6)

      omega = (am - complex(0d0,1d0)*sqrt(am*am - s0)) / (am + complex(0d0,1d0)*sqrt(am*am - s0))
      ! the overall phase is tuned by the factor fit outside
      ! I can't think of a good way of not having the problem that multiplying everthing by a
      ! constant factor can be absorbed by the outside factor right now.  Let's hope it behaves.
      bwuniv(1) = const + omega*cmplx(linRe,linIm,8) + omega**2*cmplx(quadRe,quadIm,8)

      !print "(2f7.3)", bwuniv(1)
      end


      subroutine BWa2_a2prime_poly(am, par, bwuniv)
      implicit none
      real*8 am, par(50)
      complex*16 bwuniv(10)

      ! a2 BW, another BW, each deformed by poly in
      ! omega given below

      complex*16 bw_a2, bw_a2p
      complex*16 omega
      complex*16  BWSIMP_D, BW1L_D,  BW2L_D,  BGEXP_D,  BGEXP_POW_D, BWA1_BOWLER_D

      real*8 am1a, am1b, r1, l1, am2a, am2b, r2, l2, x2
      real*8 parbw(11)
      real*8 linRe, linIm, quadRe, quadIm, cstReP, cstImP, linReP, linImP, quadImP, quadReP

      real*8 mPi, mEta, s0
      parameter(mPi = 0.139, mEta = 0.548, s0 = (mPi+mEta)**2)

      ! 9 parameters:
      !   as in BW2L_D
      ! 14 variables:
      !   BW mass
      !   BW width
      !   real part of linear term
      !   imag part of linear term
      !   real part of quadratic term
      !   imag part of quadratic term
      !   BW' mass
      !   BW' width
      !   real part of constant
      !   imag part of constant
      !   real part of linear term
      !   imag part of linear term
      !   real part of quadratic term
      !   imag part of quadratic term

      bw_a2 = BW2L_D(am, par)
      linRe = par(12)
      linIm = par(13)
      quadRe = par(14)
      quadIm = par(15)

      bw_a2p = BWSIMP_D(am, par(16))
      cstReP = par(18)
      cstImP = par(19)
      linReP = par(20)
      linImP = par(21)
      quadReP = par(22)
      quadImP = par(23)

      !print "(6f7.3)", par(1:6)

      omega = (am - complex(0d0,1d0)*sqrt(am*am - s0)) / (am + complex(0d0,1d0)*sqrt(am*am - s0))
      bwuniv(1) = bw_a2 * (1 + omega*cmplx(linRe,linIm,8) + omega**2*cmplx(quadRe,quadIm,8))
     1          +  cmplx(cstReP,cstImP,8)*bw_a2p*(1 + omega*cmplx(linReP,linImP,8) + omega**2*cmplx(quadReP,quadImP,8))

      !print "(2f7.3)", bwuniv(1)
      end
C-------------
      subroutine BWconformal(am, par, bwuniv)
      real*8 am, par(50)
      complex*16 bwuniv(10)

      ! Breit-Wigner deformed by a polynomial in the conformal variable
      ! omega given below
      ! 6 parameters:
      !   BW mass
      !   BW width
      !   real part of linear term
      !   imag part of linear term
      !   real part of quadratic term
      !   imag part of quadratic term

      complex*16 bw
      complex*16 omega
      complex*16  BWSIMP_D, BW1L_D,  BW2L_D,  BGEXP_D,  BGEXP_POW_D, BWA1_BOWLER_D

      real*8 parbw(11)
      real*8 linRe, linIm, quadRe, quadIm

      real*8 mPi, mEta, s0
      ! could be smarter to start at eta'pi threshold, but I want to use the same parameters for pi-eta and pi-eta', ideally
      parameter(mPi = 0.139, mEta = 0.548, s0 = (mPi+mEta)**2)

      bw = BWSIMP_D(am, par)
      linRe = par(3)
      linIm = par(4)
      quadRe = par(5)
      quadIm = par(6)

      !print "(6f7.3)", par(1:6)

      omega = (am - complex(0d0,1d0)*sqrt(am*am - s0)) / (am + complex(0d0,1d0)*sqrt(am*am - s0))
      bwuniv(1) = bw * (1 + omega*cmplx(linRe,linIm,8) + omega**2*cmplx(quadRe,quadIm,8))

      !print "(2f7.3)", bwuniv(1)
      end
C-------------------------------------------------------------------------------------------

      real *8 function pav_rhopi_s_d(am)
      real *8 am, x1, x2, y1, y2, a, b
      real pav(401)

      data pav /
     *   0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     *   0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     *   0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     *   0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     *   0.0000E+00, 0.0000E+00, 0.1313E-09, 0.8940E-07, 0.6175E-06, 0.2020E-05, 0.4786E-05, 0.9458E-05, 0.1664E-04, 0.2701E-04,
     *   0.4131E-04, 0.6037E-04, 0.8510E-04, 0.1165E-03, 0.1558E-03, 0.2040E-03, 0.2628E-03, 0.3335E-03, 0.4179E-03, 0.5178E-03,
     *   0.6355E-03, 0.7732E-03, 0.9337E-03, 0.1120E-02, 0.1335E-02, 0.1583E-02, 0.1867E-02, 0.2194E-02, 0.2568E-02, 0.2995E-02,
     *   0.3483E-02, 0.4039E-02, 0.4673E-02, 0.5396E-02, 0.6220E-02, 0.7160E-02, 0.8233E-02, 0.9458E-02, 0.1086E-01, 0.1246E-01,
     *   0.1430E-01, 0.1641E-01, 0.1884E-01, 0.2163E-01, 0.2484E-01, 0.2853E-01, 0.3277E-01, 0.3759E-01, 0.4306E-01, 0.4917E-01,
     *   0.5591E-01, 0.6322E-01, 0.7100E-01, 0.7913E-01, 0.8752E-01, 0.9604E-01, 0.1046E+00, 0.1132E+00, 0.1218E+00, 0.1302E+00,
     *   0.1386E+00, 0.1469E+00, 0.1551E+00, 0.1631E+00, 0.1711E+00, 0.1790E+00, 0.1867E+00, 0.1944E+00, 0.2020E+00, 0.2095E+00,
     *   0.2169E+00, 0.2243E+00, 0.2315E+00, 0.2387E+00, 0.2458E+00, 0.2529E+00, 0.2599E+00, 0.2668E+00, 0.2737E+00, 0.2805E+00,
     *   0.2873E+00, 0.2940E+00, 0.3007E+00, 0.3073E+00, 0.3138E+00, 0.3204E+00, 0.3269E+00, 0.3333E+00, 0.3397E+00, 0.3461E+00,
     *   0.3525E+00, 0.3587E+00, 0.3650E+00, 0.3713E+00, 0.3775E+00, 0.3837E+00, 0.3898E+00, 0.3959E+00, 0.4020E+00, 0.4081E+00,
     *   0.4141E+00, 0.4201E+00, 0.4261E+00, 0.4320E+00, 0.4380E+00, 0.4439E+00, 0.4498E+00, 0.4556E+00, 0.4615E+00, 0.4673E+00,
     *   0.4731E+00, 0.4790E+00, 0.4847E+00, 0.4905E+00, 0.4962E+00, 0.5019E+00, 0.5076E+00, 0.5134E+00, 0.5189E+00, 0.5246E+00,
     *   0.5303E+00, 0.5359E+00, 0.5415E+00, 0.5471E+00, 0.5526E+00, 0.5582E+00, 0.5638E+00, 0.5693E+00, 0.5749E+00, 0.5804E+00,
     *   0.5858E+00, 0.5914E+00, 0.5968E+00, 0.6023E+00, 0.6077E+00, 0.6132E+00, 0.6186E+00, 0.6241E+00, 0.6294E+00, 0.6348E+00,
     *   0.6403E+00, 0.6456E+00, 0.6510E+00, 0.6563E+00, 0.6617E+00, 0.6671E+00, 0.6724E+00, 0.6777E+00, 0.6830E+00, 0.6882E+00,
     *   0.6936E+00, 0.6990E+00, 0.7041E+00, 0.7095E+00, 0.7149E+00, 0.7199E+00, 0.7252E+00, 0.7305E+00, 0.7356E+00, 0.7410E+00,
     *   0.7462E+00, 0.7514E+00, 0.7567E+00, 0.7619E+00, 0.7668E+00, 0.7723E+00, 0.7774E+00, 0.7826E+00, 0.7878E+00, 0.7930E+00,
     *   0.7982E+00, 0.8033E+00, 0.8084E+00, 0.8135E+00, 0.8188E+00, 0.8239E+00, 0.8291E+00, 0.8340E+00, 0.8393E+00, 0.8444E+00,
     *   0.8493E+00, 0.8547E+00, 0.8597E+00, 0.8649E+00, 0.8700E+00, 0.8750E+00, 0.8800E+00, 0.8851E+00, 0.8903E+00, 0.8953E+00,
     *   0.9005E+00, 0.9054E+00, 0.9105E+00, 0.9156E+00, 0.9205E+00, 0.9256E+00, 0.9308E+00, 0.9358E+00, 0.9408E+00, 0.9458E+00,
     *   0.9507E+00, 0.9560E+00, 0.9609E+00, 0.9659E+00, 0.9711E+00, 0.9760E+00, 0.9808E+00, 0.9860E+00, 0.9909E+00, 0.9960E+00,
     *   0.1001E+01, 0.1006E+01, 0.1011E+01, 0.1016E+01, 0.1021E+01, 0.1026E+01, 0.1031E+01, 0.1036E+01, 0.1041E+01, 0.1046E+01,
     *   0.1051E+01, 0.1056E+01, 0.1061E+01, 0.1066E+01, 0.1071E+01, 0.1076E+01, 0.1081E+01, 0.1085E+01, 0.1090E+01, 0.1096E+01,
     *   0.1100E+01, 0.1105E+01, 0.1110E+01, 0.1115E+01, 0.1120E+01, 0.1125E+01, 0.1130E+01, 0.1135E+01, 0.1140E+01, 0.1145E+01,
     *   0.1150E+01, 0.1154E+01, 0.1160E+01, 0.1164E+01, 0.1169E+01, 0.1174E+01, 0.1179E+01, 0.1184E+01, 0.1189E+01, 0.1194E+01,
     *   0.1199E+01, 0.1204E+01, 0.1208E+01, 0.1214E+01, 0.1218E+01, 0.1223E+01, 0.1228E+01, 0.1233E+01, 0.1238E+01, 0.1243E+01,
     *   0.1248E+01, 0.1253E+01, 0.1257E+01, 0.1262E+01, 0.1267E+01, 0.1272E+01, 0.1277E+01, 0.1282E+01, 0.1287E+01, 0.1292E+01,
     *   0.1296E+01, 0.1301E+01, 0.1306E+01, 0.1311E+01, 0.1316E+01, 0.1321E+01, 0.1326E+01, 0.1330E+01, 0.1336E+01, 0.1340E+01,
     *   0.1345E+01, 0.1350E+01, 0.1355E+01, 0.1359E+01, 0.1365E+01, 0.1369E+01, 0.1374E+01, 0.1379E+01, 0.1384E+01, 0.1389E+01,
     *   0.1394E+01, 0.1398E+01, 0.1404E+01, 0.1408E+01, 0.1412E+01, 0.1418E+01, 0.1422E+01, 0.1427E+01, 0.1432E+01, 0.1437E+01,
     *   0.1442E+01, 0.1447E+01, 0.1451E+01, 0.1457E+01, 0.1461E+01, 0.1466E+01, 0.1472E+01, 0.1475E+01, 0.1480E+01, 0.1486E+01,
     *   0.1490E+01, 0.1495E+01, 0.1500E+01, 0.1504E+01, 0.1510E+01, 0.1514E+01, 0.1518E+01, 0.1524E+01, 0.1529E+01, 0.1534E+01,
     *   0.1538E+01, 0.1542E+01, 0.1549E+01, 0.1552E+01, 0.1557E+01, 0.1562E+01, 0.1567E+01, 0.1573E+01, 0.1577E+01, 0.1581E+01,
     *   0.1586E+01, 0.1592E+01, 0.1595E+01, 0.1601E+01, 0.1605E+01, 0.1610E+01, 0.1616E+01, 0.1619E+01, 0.1625E+01, 0.1630E+01,
     *   0.1634E+01, 0.1639E+01, 0.1644E+01, 0.1648E+01, 0.1654E+01, 0.1658E+01, 0.1663E+01, 0.1668E+01, 0.1672E+01, 0.1678E+01,
     *   0.1682E+01, 0.1687E+01, 0.1692E+01, 0.1696E+01, 0.1701E+01, 0.1708E+01, 0.1710E+01, 0.1716E+01, 0.1721E+01, 0.1724E+01,
     *   0.1726E+01  /
         data xmin /0. /
         data xmax /4. /
         data step /0.01 /

         if( am .ge. xmin .and. am .le. xmax ) then
            
             nb = ( am - xmin )/step + 1 

             x1 = xmin + (nb-1)*step
             x2 = x1 + step	   

             y1 = pav( nb )
             y2 = pav( nb+1)

c            write(0,*) ' xm,k= ', xm,k, ' nb= ', nb, ' x1,x2,y1,y2= ',  x1,x2,y1,y2, ' xmass1nb=', xmass1(nb)

             a = (y2-y1)      /(x2-x1)
             b = (x2*y1-x1*y2)/(x2-x1)

c            write(0,*) ' am= ',am, ' nb= ',nb 
c            write(0,*) ' x12= ',x1,x2,' y12= ',y1,y2

             pav_rhopi_s_d = am*a + b 

            else
	    
            pav_rhopi_s_d = 0.

         endif

      return
      end
c------------------------------------------------------------------------
      subroutine BWa2_a2prime_poly_expo(am, par, bwuniv)
      implicit none
      real*8 am, par(50)
      complex*16 bwuniv(10)

      ! a2 BW, another BW, each deformed by poly in
      ! omega given below

      complex*16 bw_a2, bw_a2p, bg
      complex*16 omega
      complex*16  BWSIMP_D, BW1L_D,  BW2L_D,  BGEXP_D,  BGEXP_POW_D, BWA1_BOWLER_D

      real*8 am1a, am1b, r1, l1, am2a, am2b, r2, l2, x2
      real*8 parbw(11)
      real*8 linRe, linIm, quadRe, quadIm, cstReP, cstImP, linReP, linImP, quadImP, quadReP, cstReBG, cstImBG

      real*8 mPi, mEta, s0
      parameter(mPi = 0.139, mEta = 0.548, s0 = (mPi+mEta)**2)

      ! 9 parameters:
      !   as in BW2L_D
      ! 17 variables:
      !   BW mass
      !   BW width
      !   real part of linear term
      !   imag part of linear term
      !   real part of quadratic term
      !   imag part of quadratic term
      !   BW' mass
      !   BW' width
      !   real part of constant
      !   imag part of constant
      !   real part of linear term
      !   imag part of linear term
      !   real part of quadratic term
      !   imag part of quadratic term
      !   alpha (expo)
      !   real part
      !   imag part

      bw_a2 = BW2L_D(am, par)
      linRe = par(12)
      linIm = par(13)
      quadRe = par(14)
      quadIm = par(15)

      bw_a2p = BWSIMP_D(am, par(16))
      cstReP = par(18)
      cstImP = par(19)
      linReP = par(20)
      linImP = par(21)
      quadReP = par(22)
      quadImP = par(23)

      parbw(1) = par(16)
      parbw(2) = par(17)
      parbw(3) = par(24)
      bg = BGEXP_D(am, parbw)
      cstReBG = par(25)
      cstImBG = par(26)

      !print "(6f7.3)", par(1:6)

      omega = (am - complex(0d0,1d0)*sqrt(am*am - s0)) / (am + complex(0d0,1d0)*sqrt(am*am - s0))
      bwuniv(1) = bw_a2 * (1 + omega*cmplx(linRe,linIm,8) + omega**2*cmplx(quadRe,quadIm,8))
     1          +  cmplx(cstReP,cstImP,8)*bw_a2p*(1 + omega*cmplx(linReP,linImP,8) + omega**2*cmplx(quadReP,quadImP,8))
     2          +  cmplx(cstReBG,cstImBG,8)*bg

      !print "(2f7.3)", bwuniv(1)
      end
C------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine  gener_calc_ampl_all(ireact, iwaves,  xmd)
      real *8 xmd

c----------- amplitudes, names,   quantum numbers from calc_ampl (physics-dependent)
      parameter  ( nampmax = 100 ) ! 4pi-generator currently ...and eta-eta 

      common /calc_ampl_names/ name_calc_ampl(nampmax)
      character *60  name_calc_ampl

      common /calc_ampl_quantum_numb/  J1 (nampmax),  JM (nampmax),  IP(nampmax), IETA(nampmax)  

      common /calc_ampl_totals/ namp_tot

      common /calc_ampl_init/ ifq
      
      common /react_mass_limits/ xm_min_ir(10), xm_max_ir(20), step_ir(10)
      
      common /ps_all_ampl_react/  nreact_tot,  n_amp_calc(10), nam_amp_calc( 100, 10 ), nam_react(10)
      character *60  nam_amp_calc
      character *20 nam_react
      
      common /nmc_all_react/ nmc_ir(10)


      if( ireact .eq. 0) then
          write(0,*) '  gener_calc_ampl_all  initialization  ' 
      endif
      
      ireact1 = 0
      
c-------------------------------------------------------
      ireact1 = ireact1 + 1
      if( ireact .eq. ireact1 .or. ireact .eq. 0 ) then
        if( ireact .eq. 0) then
          nam_react(ireact1) = '2pi+2pi-'
          xm_min_ir(ireact1) = 0.8
          xm_max_ir(ireact1) = 2.5
          step_ir(ireact1)   = 0.010
          nmc_ir(ireact1)    = 10000
          ifq = 0
        endif
        
CSdv-        if(ireact .ge. 0) call gener_4pic(xmd)
CSdv-        if(iwaves .ne. 0) call calc_ampl_4pic
        
        if( ireact .eq. 0) then
          n_amp_calc(ireact1) = namp_tot
          write(0,*) ' namp_tot = ', namp_tot
          do ka = 1, namp_tot
            nam_amp_calc( ka, ireact1 ) = name_calc_ampl(ka)
          write(0,*) ' nam(' ,ka, ')  = ', name_calc_ampl(ka)
          enddo
        endif
      endif
      
c-------------------------------------------------------
          ireact1 = ireact1 + 1
      if( ireact .eq. ireact1 .or.  ireact .eq. 0) then
        if( ireact .eq. 0) then
          nam_react(ireact1)  = '3pic'
           xm_min_ir(ireact1) = 0.5
           xm_max_ir(ireact1) = 2.5
           step_ir(ireact1)   = 0.010
           nmc_ir(ireact1)    = 10000
           ifq = 0
        endif

CSdv-   if(ireact .ge. 0) call gener_3pic(xmd)
CSdv-   if(iwaves .ne. 0) call calc_ampl_3pic

        if( ireact .eq. 0) then
          n_amp_calc(ireact1) = namp_tot
          do ka = 1, namp_tot
            nam_amp_calc( ka, ireact1 ) = name_calc_ampl(ka)
          enddo
        endif
      endif

c-------------------------------------------------------
          ireact1 = ireact1 + 1
      if( ireact .eq. ireact1 .or.  ireact .eq. 0 ) then
        if( ireact .eq. 0) then
          nam_react(ireact1)  = '2eta'
           xm_min_ir(ireact1) = 1.100
           xm_max_ir(ireact1) = 2.500
           step_ir(ireact1)   = 0.001
           nmc_ir(ireact1)    = 1
           ifq = 0
        endif

        if(ireact .ge. 0) call gener_2eta(xmd)
        if(iwaves .ne. 0) call calc_ampl_2eta

        if( ireact .eq. 0) then
          n_amp_calc(ireact1) = namp_tot
          write(0,*) ' namp_tot = ', namp_tot
          do ka = 1, namp_tot
            nam_amp_calc( ka, ireact1 ) = name_calc_ampl(ka)
          write(0,*) ' nam(',ka , ')  = ', name_calc_ampl(ka)
          enddo
        endif
      endif

c-------------------------------------------------------
          ireact1 = ireact1 + 1
      if( ireact .eq. ireact1 .or.  ireact .eq. 0   ) then
        if( ireact .eq. 0) then
          nam_react(ireact1)  = 'etapic'
           xm_min_ir(ireact1) = 0.686
           xm_max_ir(ireact1) = 2.500
           step_ir(ireact1)   = 0.001
           nmc_ir(ireact1)    = 1
           ifq = 0
        endif

        if(ireact .ge. 0) call gener_etapic(xmd)
        if(iwaves .ne. 0) call calc_ampl_etapic

        if( ireact .eq. 0) then
          n_amp_calc(ireact1) = namp_tot
          do ka = 1, namp_tot
            nam_amp_calc( ka, ireact1 ) = name_calc_ampl(ka)
          enddo
        endif
      endif

c-------------------------------------------------------
          ireact1 = ireact1 + 1
      if( ireact .eq. ireact1 .or.  ireact .eq. 0   ) then
        if( ireact .eq. 0) then
          nam_react(ireact1)  = '5pic'
           xm_min_ir(ireact1) = 1.200
           xm_max_ir(ireact1) = 2.700
           step_ir(ireact1)   = 0.010
           nmc_ir(ireact1)    = 10000
           ifq = 1
        endif


CSdv-   if(ireact .ge. 0) call gener_5pic(xmd)
CSdv-   if(iwaves .ne. 0) call calc_ampl_5pic

        if( ireact .eq. 0) then

            write(0,*) '-=====================================-'
            write(0,*) ' after call calc_ampl_5pic namp_tot = ' ,  namp_tot
            write(0,*) '-=====================================-'

          n_amp_calc(ireact1) = namp_tot
          do ka = 1, namp_tot
            nam_amp_calc( ka, ireact1 ) = name_calc_ampl(ka)
             
            write(0,*) '-=====================================-'
            write(0,*) ' 5pic name ',ka,  name_calc_ampl(ka)
            write(0,*) '-=====================================-'

          enddo
        endif
      endif

      if(  ireact .eq. 0 ) then
        nreact_tot = ireact1
      endif

      if( ireact .eq. 0) then
          write(0,*) '  gener_calc_ampl_all  initialization  nreact_tot =  ' ,  nreact_tot  
      endif

      return
      end
C---------------------------------------------
      subroutine gener_dp1( xmd,  weightd )
      real *8 xmd, weightd
      real xx(10)

      real *8 pbeam, ambeam, am1, am2,  amtarg, amx, s, sqs, s1, s2 
      real *8 dams2min, dams2pl, sdams2pl, t1q, t2q
      real *8 traj_p1, traj_r1, traj_p, traj_r, FUNCBD

c     common / ts_xs/  t1, t2, x1, x2
c     real *8 t1, t2, x1, x2

      COMMON/RLRCOM3/PLAB(4,30)
      real *8 PLAB

      real *8 pcm(4), pcm1(4), pfastcm(4), pslowcm(4), pxcm(4), dxfeyn

      common /cm_vectors/  pa(4),pb(4),pa1(4),pb1(4),px(4)
      real *8 pa,pb,pa1,pb1,px

c     real *8	pdecprod(4,10)

      real *8 pbm(4),ax,ay

      common /vertex_simulated/ vertex(3)

c----------- c.m  of cent. produced system  

      real *8 pbeam1(4), pbeam2(4), pres(4), pbeam1_cm(4), pbeam2_cm(4), pslow_cm(4), pfast_cm(4)
      real *8 pbeam_cm(4), ptarg_cm(4)
      real *8 ygj1(4), ygj2(4)
      real *8 ctd

c----- cuts on CP variables, independent of the channel
      common /y_cent_cut/ yminc, ymaxc
      common /x_fast_slow_cuts/ x1min, x1max, x2min, x2max
      common /t_fast_slow_cuts/ t1min, t1max, t2min, t2max
      
      real *8 p_pi, amx2, wpipi
      real *8 q2a, q2b

c--------- begin of CP-generator  variables  commons, read from file with ascii cards...

      common /pbeam_f/ pbeamf,  sigma_pbeam
      
      common /itype_beam/  itypb1 

      common /mass_limits/ xma, xmb

      common /number_of_mc_events/ nmc1

      common /mode_reject/ imoderej

c--------- end of CP-generator  variables  commons ---------

      common /fast_geant_code/ icodegeant_fast 

c----------------------------------------------
       real *8   pb6(4),px6(4)

       COMMON /NUCLEONS/ Wmtarg,Wmrec
       DOUBLE  PRECISION  Wmtarg,Wmrec

C----------------------------------------
       common /mcbeaml/ pb1a(4),prec1(4)
       real*8 pb1a,prec1
c----------------------------------------
       real *8 ptot_mc(4), ptr_mc(4)
c----------------------------------------
       common / tinv_generated/ tinv

       real *8 pbeam_calc(4)

       real *4 t_calc, tprime_calc

c--------------- generated m_target, recoil  ----------------------------------
       common / target_recoil_from_cards/ amtarget, amrecoil

c---------------  m_target, recoil used in beam_calc --------------------------
       common / target_recoil_recon_from_cards/ amtarget_recon, amrecoil_recon

c===================================================
       common / tprime_generated/ t1_gener

c------------------------------------
       parameter (AMPIC = 0.13957018)

       real *8   ex, pxc, pbc, eb

       real *8 q2beam

       pbeam = pbeamf

c      write(0,*) ' itypb1 = ', itypb1

c      if( itypb1 .eq. 1) then
c    	 ambeam = 0.93827231
c    	 icodegeant_fast = 14
c      endif

       ambeam = AMPIC

       call  mcbeam(pb6,vx,vy,vz,ax,ay)

       vertex(1) = vx
       vertex(2) = vy
       vertex(3) = vz

       do k = 1,4
         plab(k,1) = pb6(k)
       enddo

       am_min = xma
       am_max = xmb

c------------------------------------------------------------------
       nmes  = 0 ! until gener is called and nmes is defined...

c      amtarg =0.93827231

       wmtarg = amtarget
       wmrec  = amrecoil

       plab(1,2) = 0.
       plab(2,2) = 0.
       plab(3,2) = 0.
       plab(4,2) = wmtarg

       do i = 1,4
          pcm(i) =  plab(i,1) +  plab(i,2) 
       enddo

       s = pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2
       sqs = dsqrt( s )

       CALL BOOST2(sqs, pcm(1), plab(1,1),  pa(1))
       CALL BOOST2(sqs, pcm(1), plab(1,2),  pb(1))

c--------------------------------------------------------

	CALL PP2RANLUX(xx(1), 1)
        xmd = xma + (xmb-xma) * xx(1)

c       write(0,*) ' xmd = ', xmd
c       write(0,*) ' pbeam = ',  pb6(1), pb6(2), pb6(3), pb6(4)
        
        q2beam = ambeam**2
        call generate_t(Tinv, xmd, Wmtarg, Wmrec, q2beam, pb6(1))
        
c       write(0,*) ' tinv = ', tinv

        call pxlab(pb6(1), Wmtarg, Wmrec, xmd, Tinv, px6(1))
        
        weightd = 1/xmd**1.7  ! for 190-GeV beam from pomeron trajectory effective behavior

        do k4 = 1,4
           PLAB(k4,7) =  px6(k4)
           PLAB(k4,4) = PLAB(k4,1) + PLAB(k4,2) -  PLAB(k4,7)
           pbeam1(k4) = PLAB(k4,1) 
           pbeam2(k4) = PLAB(k4,2) - PLAB(k4,4) 
           pres  (k4) = PLAB(k4,7) ! in case when decay is absent...
        enddo
        
c       xmd = dsqrt( pres(4)**2 - pres(1)**2-  pres(2)**2- pres(3)**2  )
        
        CALL BOOST2( xmd, pres(1), pbeam1(1),pbeam1_cm(1) )
        CALL BOOST2( xmd, pres(1), pbeam2(1),pbeam2_cm(1) )
        CALL BOOST2( xmd, pres(1), PLAB(1,4),pslow_cm(1) )
        CALL BOOST2( xmd, pres(1), PLAB(1,1),pbeam_cm(1) )
        CALL BOOST2( xmd, pres(1), PLAB(1,2),ptarg_cm(1) )

        call vec(pslow_cm(1),pbeam2_cm(1), ygj2(1)  )

c====================== check of t-generation , and beam_calc procedure 

c----------------- using beam_calc, assuming  different target, recoil  masses ..... -------

	if( amtarget_recon .ne. 0. .and. amrecoil_recon .ne. 0.) then
c	    write(0,*) ' tinv = ', tinv
c	    write(0,*) ' pb = ', pb(1), pb(2),pb(3),pb(4)
c	    write(0,*) ' px = ', px(1), px(2),px(3),px(4)

	write(0,*) ' amtarget_gener,  amrecoil_gener = ', amtarget,  amrecoil
	write(0,*) ' amtarget_recon,  amrecoil_recon = ', amtarget_recon,  amrecoil_recon
	write(0,*) '------------------------------------------------------------------------- '
		
c------------------- check and printout -------------------------------------
	call  beam3_calc(px6(1) , amtarget_recon,  amrecoil_recon,  pb6(1),  pbeam_calc(1),  t_calc, tprime_calc )

c	write(0,*) '=========================================================== '	

	tinv_recalc = ( pb6(4) - px6(4) ) **2 - ( pb6(1) - px6(1) ) **2 - ( pb6(2) - px6(2) ) **2 - ( pb6(3) - px6(3) ) **2

c       write(0,*) 'rexp: amx = ', amx, ' TMIN,TMAX = ', TMIN,TMAX

	S = AMPIC**2+Wmtarg**2+2.*PB6(4)*Wmtarg
	shs = dsqrt(s)
	ex  = (s+amx**2-Wmrec**2)/(2.*shs)
	pxc = dsqrt(ex**2-amx**2)
	eb  = (s+ampic**2-Wmtarg**2)/(2.*shs)
	pbc = dsqrt(eb**2-ampic**2)
C-----------------------
	tminPH = ampic**2+amx**2-2.*eb*ex-2.*pbc*pxc
	tmaxPH = ampic**2+amx**2-2.*eb*ex+2.*pbc*pxc
C----------------------------

c	write(0,*) ' tprime_calc = ', tprime_calc
c	write(0,*) ' tmin = ', abs(tmaxPH )
	
	t1_recalc = abs(tinv_recalc) - abs(tmaxPH )

c-	write(0,*) ' tprime_generated = ', t1, ' tprime_event = ', t1_recalc
c-	write(0,*) ' t_generated = ',  abs(tinv), ' t_event = ', abs(tinv_recalc)

	write(0,777) abs(tinv), abs(tinv_recalc),  abs(tinv_recalc)-abs(tinv)
 777	format ('         t_gener =  ', f15.9,'            t_event =  ', f15.9,'       t_event-t_gener        = ', f15.9) 

	write(0,776) t1_gener, t1_recalc,  t1_recalc- t1_gener
 776	format ('    tprime_gener =  ', f15.9,'       tprime_event =  ', f15.9,'  tprime_event-tprime_gener   = ', f15.9) 

	write(0,*) '------------------------------------------------------------------------- '	

	write(0,772)  abs(tinv) , t_calc, t_calc-abs(tinv)
 772	format ('    t_gener =  ', f15.9,'	 t_calc =  ', f15.9,'	     t_calc-t_gen = ', f15.9) 
c       
	de_gener =pb(4)-px(4)
	de_calc =pbeam_calc(4)-px(4)
	write(0,773)  de_gener,  de_calc, de_calc-de_gener 
 773	format ('ea-ec(gener)=  ', f15.9,' ea-ec(calc)   = ', f15.9,'  de(calc)-de(gener) = ', f15.9)
c       
c--------------------------------------

	write(0,*) '=========================================================== '	
	write(0,*) '=========================================================== '	

	endif
c------------------------------------------------------------------

        return
        end

c------------------------------------------------------------------
CSdv- Eta pi channel
C
      subroutine gener_etapic(xmd)
      real *8 xmd
      
c--------weight for n-meson phase space ---------
      common /gener_weight/ weight_gen
      
      common /decay_cms_momentum/ pdec(4,10) 
      real *8 pdec

      common /decay_cms_nmesons/ nmes 
      common /decay_geant_codes/ icodegeant(10) 

      PARAMETER (AMET = 0.54745)
      PARAMETER (AMPIC= 0.13956755)

      real X(2)

      real *8 P1,P2,E1,E2, AM12,am1,am2, PC1(4),PC2(4)
c---------------------------------------------
      common /mode_random/ imodrn
      
      AM1 = AMPIC
      AM2 = AMETA
      
      am12 = xmd

      nmes = 2 ! here the number of "physical particles" -- mesons is defined
      icodegeant(1) = 9
      icodegeant(2) = 17
C     
      E2=(AM12**2+AM2**2-AM1**2)/(2.*AM12)
      P2=DSQRT( DABS(E2**2-AM2**2) )
      
      E1=AM12-E2

      if( xmd .ge. AM1+AM2 ) then
        weight_gen = P2  !  d M , corrected !!!
      else
        weight_gen = 0.
        return
      endif

c      write(0,*) ' xmd = ', xmd, ' weight = ', weight_gen


          if( imodrn .eq. 0) then
            CALL PP2RANLUX(X(1),2)
          else
            call RNDH(X(1),2)
          endif

c     CALL PP2RANLUX(X(1),2)
      
      CALL SFO ( E2,P2,PC2(1), X(1), X(2) )
      
      DO I = 1,3
         PC1(I) = -PC2(I)
      ENDDO
      
         PC1(4) = E1
	 PC2(4) = E2		

	do i = 1,4
	  pdec(i,1) = PC1(i)
	  pdec(i,2) = PC2(i)
	enddo
c----------------------------------------

      return
      end

      subroutine calc_ampl_etapic

      common /decay_cms_momentum/ pdec(4,10) 
      real *8 pdec

c----------- amplitudes, names,   quantum numbers from calc_ampl (physics-dependent)
      parameter  ( nampmax = 100 )

      common /calc_ampl_amplits/ ampl(nampmax)
      complex  ampl

      common /calc_ampl_samplits/ sampl(nampmax)
      real  sampl

      common /calc_ampl_names/ name_calc_ampl(nampmax)
      character *60  name_calc_ampl

c------------------------------------ total spin --- spin projection -- reflectivity(if def.) --
      common /calc_ampl_quantum_numb/  J1 (nampmax),    JM (nampmax),  IP(nampmax), IETA(nampmax)  

      common /calc_ampl_init/ ifq

      common /calc_ampl_totals/ namp_tot

      real *8 pc12(4), pc1(4), pc2(4)
      real *8 samx, amx

      real *8 r1(4)

      real trr1(5), trrrr1(9)

      complex ampl_tmp( 20 )
      complex ci

      DO I = 1,4
        PC12(I)=pdec(I,1)+pdec(I,2)
        PC1(I)=pdec(I,1)
        PC2(I)=pdec(I,2)
      ENDDO
      
      CALL AMR(PC12(1),SAMX,AMX)
      
      do i=1,3
        r1(i)=PC1(i)
      enddo	
      
      call scl( r1(1),r1(1),rr1 )
      
      r1n = sqrt( rr1 )
      IF( r1N.EQ.0.) r1N = 0.0000001
      
      RINT = 4.94
      
      r1r = r1n*rint
      
      bp1r1 = sqrt( fd1(r1r) )
      bp2r1 = sqrt( fd2(r1r) )
      bp3r1 = sqrt( fd3(r1r) )
      bp4r1 = sqrt( fd4(r1r) )

	 kk = 0

C-----------FLAT------------
         kk = kk + 1

         ampl(kk) = 1.
         sampl(kk) = 1.

      if( ifq.eq.0) then
         name_calc_ampl(kk) = 'FLAT etapic'
         J1 (kk) = 0
         JM (kk) = 0  
      endif

C-----------0++ == FLAT ------------
        kk = kk + 1
        ampl(kk) = 1.
        sampl(kk) = 1.

      if( ifq.eq.0) then
        name_calc_ampl(kk) = '1-(0++)0- [c] [p]- s '
        J1 (kk) = 0
        JM (kk) = 0  
      endif

C-----------1-+ M=1  ------------

        kk  = kk + 1
        sampl(kk) = 0.

        do im = 1, 3
          ampl_tmp(im) = r1(im)/bp1r1
          sampl(kk) = sampl(kk)   +  ampl_tmp(im)*conjg(ampl_tmp(im))
        enddo
        sampl(kk) = sampl(kk) / 3.
        ampl (kk) = ampl_tmp(2)

      if( ifq.eq.0) then
        name_calc_ampl(kk) = '1-(1-+)1+ [c] [p]- p '
        J1 (kk) = 1
        JM (kk) = 1  
      endif

C-----------2++ M=1 ------------
	call w2ns( r1(1),r1(1),trr1(1))

        kk  = kk + 1
        sampl(kk) = 0.

        do im = 1, 5
          ampl_tmp(im) = trr1(im)/bp2r1
          sampl(kk) = sampl(kk)   +  ampl_tmp(im)*conjg(ampl_tmp(im))
        enddo
         sampl(kk) = sampl(kk) / 5.
         ampl (kk) = ampl_tmp(3)

      if( ifq.eq.0) then
        name_calc_ampl(kk) = '1-(2++)1+ [c] [p]- d '
        J1 (kk) = 2
        JM (kk) = 1
      endif
C-----------2++ M=1 ------------
	call  w4ns( r1(1),r1(1),r1(1),r1(1),trrrr1(1))

        kk = kk + 1

        sampl(kk) = 0.

        do im = 1, 9
          ampl_tmp(im) = trrrr1(im)/bp4r1
          sampl(kk) = sampl(kk)   +  ampl_tmp(im)*conjg(ampl_tmp(im))
        enddo
         sampl(kk) = sampl(kk) / 9.
         ampl(kk) =  ampl_tmp(3)

      if( ifq.eq.0) then
        name_calc_ampl(kk) = '1-(4++)1+ [c] [p]- g '
        J1 (kk) = 4
        JM (kk) = 1  
      endif

      if(ifq.eq.0) then 
        namp_tot = kk
        ifq = 1 
      endif
      
      return
      end

c--------------------------------------------------------------
      subroutine  ushist_etapic

        common /init_bin_hbook/ i_bin_hbook(300)

        character *16 namtail 
	character *17 cname,cname1 

c---------------------------------
        common /interval1/ xbins_pwa(100), nbins_pwa
c---------------------------------
c
        nbins_pwa = 8
c
        xbins_pwa(1) = 0.50
        xbins_pwa(2) = 1.00
        xbins_pwa(3) = 1.20
        xbins_pwa(4) = 1.40
        xbins_pwa(5) = 1.60
        xbins_pwa(6) = 1.90
        xbins_pwa(7) = 2.10
        xbins_pwa(8) = 2.50

        do i = 1, 300
           i_bin_hbook(i) = 0
        enddo

c------------ moments in eta- region ----

	do L = 0,8
	  if( L.ge.2) nnmax = 2
	  if( L .le. 2) nnmax = L

	 do M = 0,nnmax
cccccccccccccccccccccccccccccccccccccc

	Write( cname, 733 )  L,M
  733   format(' real H( ',I1,'  ',I1,') ' )

	Write( cname1, 734 )  L,M
  734   format(' imag H( ',I1,'  ',I1,') ' )

	id1 = 23100 + L*10+M

c$$$        call hbook1 (id1 ,cname//'for [c] [c]	$' , 40,   1.0, 3.0, 0.)
c$$$	call hbarx(  id1 )	   

	id1 = 23200 + L*10+M

c$$$        call hbook1( id1 ,cname1//'for [c] [c]   $' , 40,   1.0, 3.0 ,0.)
c$$$	call hbarx(  id1 )	   

        enddo
        enddo

        return
        end

c--------------------------------------------------------------
      subroutine  usfill_etapic(ves)
c     real *8 weight_all
      real ves
c-----------------------------------------------------
        common /interval1/ xbins_pwa(100), nbins_pwa
c-----------------------------------------------------
      common /decay_cms_momentum/ pdec(4,10) 
      real *8 pdec

	real ppmm(10,10), cty(10),sty(10)

      
      real *8 zgj(4)
      
      real *8 pc12(4), pc1(4), pc2(4)
      real *8 samx, amx

      real *8 r1(4)

C-------INV MASSES-----------
	    zgj(1) = 0.
	    zgj(2) = 0.
	    zgj(3) = 1.

      DO I = 1,4
        PC12(I)=pdec(I,1)+pdec(I,2)
        PC1(I)=pdec(I,1)
        PC2(I)=pdec(I,2)
      ENDDO
      
      CALL AMR(PC12(1),SAMX,AMX)
        xms = amx

      do i=1,3
        r1(i)=PC1(i)
      enddo	
      
      call scl( r1(1),r1(1),rr1 )
      
      r1n = sqrt( rr1 )
      IF( r1N.EQ.0.) r1N = 0.0000001

	  call scln( r1(1),zgj(1),var1r)
	     xx = r1(1)
	     yy = r1(2)

	  var2r = atan2(yy,xx)
	  beta  = acos( var1r )

	  cty(1) =  1.
	  cty(2) =-cos( var2r )
	  cty(3) = cos( 2.*var2r )
	  cty(4) =-cos( 3.*var2r )
	  cty(5) = cos( 4.*var2r )
	  cty(6) =-cos( 5.*var2r )
	  cty(7) = cos( 6.*var2r )
	  cty(8) =-cos( 7.*var2r )
	  cty(9) = cos( 8.*var2r )
 
	  sty(1) =   1.
	  sty(2) = sin( var2r )
	  sty(3) = sin( 2.*var2r )
	  sty(4) = sin( 3.*var2r )
	  sty(5) = sin( 4.*var2r )
	  sty(6) = sin( 5.*var2r )
	  sty(7) = sin( 6.*var2r )
	  sty(8) = sin( 7.*var2r )
	  sty(9) = sin( 8.*var2r )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	var1a = var1r

        call  ASLGF( 2, var1a, 0, 8, PPMM(1,1) ) 

        call  ASLGF( 2, var1a, 1, 8, PPMM(1,2) ) 

        call  ASLGF( 2, var1a, 2, 8, PPMM(1,3) ) 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       xxxx = PPMM(3,1 )*cty( 1)*ves/sqrt( 2.*2+1. )

	do L = 0,8
	  if( L.ge.2) nnmax = 2
	  if( L .le. 2) nnmax = L

	 do M = 0,nnmax
	 
	  xxxx =      PPMM(L+1,M+1 )*cty( M+1)*ves/sqrt( 2.*L+1. )
	  yyyy =      PPMM(L+1,M+1 )*sty( M+1)*ves/sqrt( 2.*L+1. )     

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	  XXXXMS = XMTOT + 0.0001
c	  
c	  if( xmid .ge. 0.02 ) xxxxms = xmid + 0.0001
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  x8 = 0.
	  
c======================================================

	  id1 = 23100 + L*10+M
c$$$	  call hfill( id1 , xms, x8,  xxxx )

	  id1 = 23200 + L*10+M
c$$$	  call hfill( id1 , xms, x8, yyyy )

c=====================================================
          enddo
          enddo

      return
      end
C
CSdv- End of the eta pi channel
C-----------------------------------------------------
    
      subroutine ushist_react (ichannel)
      
        if( ichannel .eq.  1) then
CSdv-          call ushist_4pic
        endif
        if( ichannel .eq.  2) then
CSdv-          call ushist_3pic
        endif
        if( ichannel .eq.  3) then
CSdv-          call ushist_2eta
        endif
        if( ichannel .eq.  4) then
            call ushist_etapic
        endif
      
      return
      end
      
C------------------------------------------------
      subroutine usfill_react(ichannel, weightf)

        if( ichannel .eq.  1) then
CSdv-          call usfill_4pic(weightf)
        endif
        if( ichannel .eq.  2) then
CSdv-          call usfill_3pic(weightf)
        endif
        if( ichannel .eq.  3) then
CSdv-          call usfill_2eta(weightf)
        endif
        if( ichannel .eq.  4) then
            call usfill_etapic(weightf)
        endif
      
      return
      end

C------------------------------------------------
      COMPLEX *16 FUNCTION BGEXP_D(AM, par_all)
      real *8  par_all(*) 
      real *8 am, am1, am2, alf, s, e, psq, ampor

           am1 = par_all(1)
           am2 = par_all(2)
           alf = par_all(3)
           ampor = am1 + am2

           BGEXP_D = (1.,0.)

           if( am.gt. ampor ) then

               S   = am*am
               E   = ( S + am1**2 - am2**2 )/(2.*am)
               PSQ = E**2- am1**2 
               BGEXP_D = dexp ( alf*PSQ )

           endif
	    
      return
      end
      
C------------------------------------------------
      subroutine gener_2eta(xmd)
      real *8 xmd
      
c--------weight for n-meson phase space ---------
      common /gener_weight/ weight_gen
      
      common /decay_cms_momentum/ pdec(4,10) 
      real *8 pdec

      common /decay_cms_nmesons/ nmes 
      common /decay_geant_codes/ icodegeant(10) 

      
      PARAMETER (AMET = 0.54745)

      real X(2)

      real *8 P1,P2,E1,E2, AM12,am1,am2, PC1(4),PC2(4)
      
c---------------------------------------------
      common /mode_random/ imodrn

      AM1  = AMET
      AM2  = AMET
      am12 = xmd

      nmes = 2 ! here the number of "physical particles" -- mesons is defined
      icodegeant(1) = 17
      icodegeant(2) = 17
C    
      E2=(AM12**2+AM2**2-AM1**2)/(2.*AM12)
      P2=DSQRT( DABS(E2**2-AM2**2) )
      
      E1=AM12-E2

      if( xmd .ge. AM1+AM2 ) then
        weight_gen = P2/AM12    ! better "compensated" 
      else
        weight_gen = 0.
        return
      endif

c     write(0,*) ' xmd = ', xmd, ' weight = ', weight_gen

      if( imodrn .eq. 0) then
        CALL PP2RANLUX(X(1),2)
      else
        call RNDH(X(1),2)
      endif

c     CALL PP2RANLUX(X(1),2)
      
      CALL SFO ( E2,P2,PC2(1), X(1), X(2) )
      
      DO I = 1,3
        PC1(I) = -PC2(I)
      ENDDO
      
        PC1(4) = E1
	PC2(4) = E2	       

	do i = 1,4
	  pdec(i,1) = PC1(i)
	  pdec(i,2) = PC2(i)
	enddo
c
      return
      end

c---------------------------------------------
      subroutine calc_ampl_2eta

      common /decay_cms_momentum/ pdec(4,10) 
      real *8 pdec

c----------- amplitudes, names,   quantum numbers from calc_ampl (physics-dependent)
      parameter  ( nampmax = 100 )

      common /calc_ampl_amplits/ ampl(nampmax)
      complex  ampl

      common /calc_ampl_samplits/ sampl(nampmax)
      real  sampl

      common /calc_ampl_names/ name_calc_ampl(nampmax)
      character *60  name_calc_ampl

c------------------------------------ total spin --- spin projection -- reflectivity(if def.) --
      common /calc_ampl_quantum_numb/  J1 (nampmax),    JM (nampmax),  IP(nampmax), IETA(nampmax)  

      common /calc_ampl_init/ ifq

      common /calc_ampl_totals/ namp_tot

      real *8 pc12(4), pc1(4), pc2(4)
      real *8 samx, amx

      real *8 r1(4)

      real trr1(5), trrrr1(9)

      complex ampl_tmp( 20 )
      complex ci

      DO I = 1,4
        PC12(I)=pdec(I,1)+pdec(I,2)
        PC1(I)=pdec(I,1)
        PC2(I)=pdec(I,2)
      ENDDO
      
      CALL AMR(PC12(1),SAMX,AMX)
 
      do i=1,3
        r1(i)=PC1(i)
      enddo	
      
      call scl( r1(1),r1(1),rr1 )
      
      r1n = sqrt( rr1 )
      IF( r1N.EQ.0.) r1N = 0.0000001
      
      RINT = 4.94
      r1r  = r1n*rint
      
      bp1r1 = sqrt( fd1(r1r) )
      bp2r1 = sqrt( fd2(r1r) )

      kk = 0

C-----------FLAT------------
      kk = kk + 1

      ampl(kk) = 1.
      sampl(kk) = 1.

      if( ifq.eq.0) then
        name_calc_ampl(kk) = 'FLAT 2eta'
        J1 (kk) = 0
        JM (kk) = 0  
      endif

C-----------0++ == FLAT ------------
        kk = kk + 1

        ampl(kk) = 1.
        sampl(kk) = 1.

      if( ifq.eq.0) then
        name_calc_ampl(kk) = '0++ eta eta'
        J1 (kk) = 0
        JM (kk) = 0  
      endif

C-----------2++ M=0 ------------
	call  w2ns( r1(1),r1(1),trr1(1))

        kk  = kk + 1
        sampl(kk) = 0.

        do im = 1, 5
          ampl_tmp(im) = trr1(im)/bp2r1
          sampl(kk) = sampl(kk)   +  ampl_tmp(im)*conjg(ampl_tmp(im))
        enddo
        sampl(kk) = sampl(kk) / 5.
        ampl (kk) = ampl_tmp(1)

      if( ifq.eq.0) then
        name_calc_ampl(kk) = '2++ M=0 eta eta'
        J1 (kk) = 2
        JM (kk) = 0  
      endif

      if(ifq.eq.0) then 
        namp_tot = kk
        ifq = 1 
      endif
      
      return
      end

c--------------------------------------------------------------
      subroutine  ushist_2eta

        common /init_bin_hbook/ i_bin_hbook(300)

        character *16 namtail 
	character *17 cname,cname1 

c---------------------------------
        common /interval1/ xbins_pwa(100), nbins_pwa
c---------------------------------
c
        nbins_pwa = 8
c
        xbins_pwa(1) = 0.50
        xbins_pwa(2) = 1.00
        xbins_pwa(3) = 1.20
        xbins_pwa(4) = 1.40
        xbins_pwa(5) = 1.60
        xbins_pwa(6) = 1.90
        xbins_pwa(7) = 2.10
        xbins_pwa(8) = 2.50

        do i = 1, 300
           i_bin_hbook(i) = 0
        enddo

c------------ moments in eta- region ----

	do L = 0,8
	  if( L.ge.2) nnmax = 2
	  if( L.le.2) nnmax = L

	  do M = 0,nnmax
cccccccccccccccccccccccccccccccccccccc

	  Write( cname, 733 )  L,M
  733     format(' real H( ',I1,'  ',I1,') ' )

	  Write( cname1, 734 )  L,M
  734     format(' imag H( ',I1,'  ',I1,') ' )

	  id1 = 23100 + L*10+M

c$$$          call hbook1 (id1 ,cname//'for [c] [c]   $' , 40,   1.0, 3.0, 0.)
c$$$	  call hbarx(  id1 )	     

	  id1 = 23200 + L*10+M

c$$$          call hbook1( id1 ,cname1//'for [c] [c]   $' , 40,   1.0, 3.0 ,0.)
c$$$	  call hbarx(  id1 )	     

          enddo
        enddo

       return
       end

c--------------------------------------------------------------
      subroutine  usfill_2eta(ves)
c     real *8 weight_all
      real ves
c-----------------------------------------------------
        common /interval1/ xbins_pwa(100), nbins_pwa
c-----------------------------------------------------
      common /decay_cms_momentum/ pdec(4,10) 
      real *8 pdec

      real ppmm(10,10), cty(10),sty(10)

      real *8 zgj(4)
      
      real *8 pc12(4), pc1(4), pc2(4)
      real *8 samx, amx

      real *8 r1(4)

C-------INV MASSES-----------
	    zgj(1) = 0.
	    zgj(2) = 0.
	    zgj(3) = 1.

      DO I = 1,4
        PC12(I)=pdec(I,1)+pdec(I,2)
        PC1(I)=pdec(I,1)
        PC2(I)=pdec(I,2)
      ENDDO
      
      CALL AMR(PC12(1),SAMX,AMX)
        xms = amx

      do i=1,3
        r1(i)=PC1(i)
      enddo	
      
      call scl( r1(1),r1(1),rr1 )
      
      r1n = sqrt( rr1 )
      IF( r1N.EQ.0.) r1N = 0.0000001

	  call scln( r1(1),zgj(1),var1r)
	     xx = r1(1)
	     yy = r1(2)

	  var2r = atan2(yy,xx)

	  beta = acos( var1r )

	  cty(1) =   1.
	  cty(2) =-cos( var2r )
	  cty(3) = cos( 2.*var2r )
	  cty(4) =-cos( 3.*var2r )
	  cty(5) = cos( 4.*var2r )
	  cty(6) =-cos( 5.*var2r )
	  cty(7) = cos( 6.*var2r )
	  cty(8) =-cos( 7.*var2r )
	  cty(9) = cos( 8.*var2r )

	  sty(1) =   1.
	  sty(2) = sin( var2r )
	  sty(3) = sin( 2.*var2r )
	  sty(4) = sin( 3.*var2r )
	  sty(5) = sin( 4.*var2r )
	  sty(6) = sin( 5.*var2r )
	  sty(7) = sin( 6.*var2r )
	  sty(8) = sin( 7.*var2r )
	  sty(9) = sin( 8.*var2r )

c---------------------------------------------------
	var1a = var1r

        call  ASLGF( 2, var1a, 0, 8, PPMM(1,1) ) 
        call  ASLGF( 2, var1a, 1, 8, PPMM(1,2) ) 
        call  ASLGF( 2, var1a, 2, 8, PPMM(1,3) ) 

c---------------------------------------------------

c       xxxx = PPMM(3,1 )*cty( 1)*ves/sqrt( 2.*2+1. )

	do L = 0,8
	  if( L.ge.2) nnmax = 2
	  if( L .le. 2) nnmax = L

	  do M = 0,nnmax

	  xxxx =      PPMM(L+1,M+1 )*cty( M+1)*ves/sqrt( 2.*L+1. )
	  yyyy =      PPMM(L+1,M+1 )*sty( M+1)*ves/sqrt( 2.*L+1. )     

c----------------------------------------------------------
c	  XXXXMS = XMTOT + 0.0001
c	  
c	  if( xmid .ge. 0.02 ) xxxxms = xmid + 0.0001
c----------------------------------------------------------

	  x8 = 0.
	  id1 = 23100 + L*10+M
c$$$	  call hfill( id1 , xms, x8,  xxxx )

	  id1 = 23200 + L*10+M
c$$$	  call hfill( id1 , xms, x8, yyyy )

c----------------------------------------------------------
          enddo
        enddo

      return
      end
      
c----------------------------------------------------------
        SUBROUTINE Beam3_calc(ptot4, amtarg,  amrec,  pbeam1,  pbeam4,  t, tprime )
        real *8 ptot4(4), pbeam1(4), pbeam4(4)
        real *8 pbeam1a

        real *8 tmp(5),test(5),dotd
        real *8 norm,xm2,a,b,c,cost2,c1,c2,Eb,pz,ps
        real *8 pb1
        real *8 xmcms,encms,pncms,enlmin
         
c       parameter ( amprot = 0.93827231 )
       	PARAMETER ( AMPIC  = 0.13956755 )
  
c       INCLUDE 'physcom.inc'
c       INCLUDE 'usdst2.inc'
c       ptot = Ptot1_2pi0(tmp)

        amprot = amtarg

        pbeam1a = dsqrt( pbeam1(1)**2 + pbeam1(2)**2 + pbeam1(3)**2 )
        axb = pbeam1(1)/pbeam1a
        ayb = pbeam1(2)/pbeam1a

        xm2 = ptot4(4)**2 - ptot4(1)**2 - ptot4(2)**2 - ptot4(3)**2

cc      write(0,*) ' xm2= ',xm2

cc	print*,'tmp1,2,3,4=',tmp(1),tmp(2),tmp(3),tmp(4)
cc	print*,'amp,amn=',amp,amn

        wm2 = xm2               ! output
c-
        norm=dsqrt(1.d0+dble(axb**2+ayb**2))
        cost2=( axb*ptot4(1) + ayb*ptot4(2) + ptot4(3) )/norm
        cost2=cost2**2
c-
c--------------------------------------------------
c       c1      =amn**2-amp**2-apic**2-xm**2+2.*ec*amp
c       a       =4.*((ec-amp)**2-ptot**2*cost2)
c       b       =4.*c1*(ec-amp)
c       c       =c1**2+4.*apic**2*ptot**2*cost2
c       eb      =(-dsqrt(b*b-4.*a*c)-b)/2./a
c       pb1     =dsqrt(eb*eb-apic**2)
c       t       =2.*(eb*ec-pb1*ptot*cost)-apic**2-xm**2
c                                       !Computation of Tmin
c       xmcms   =sqrt(2.*eb*amp+apic**2+amp**2)
c       encms   =(xmcms**2+amn**2-xmass**2)/2./xmcms
c       pncms   =sqrt(encms**2-amn**2)
c       enlmin  =((eb+amp)*encms-pb1*pncms)/xmcms
c       tmin    =2.*enlmin*amp-amp**2-amn**2
c       t       =t-tmin
c-------------------------------------------
c       amrec2 = amrec*amrec

        c1=ptot4(4)-amprot
        c2=2.*ptot4(4)*amprot - xm2 -ampic**2 + amrec**2 - amprot**2
c----
        a=4. * ( c1**2 - cost2 )
        b=4. * c1 * c2
        c=c2**2 + 4.*ampic**2*cost2
c----
        IF ( (b**2-4.*a*c).LT.0. ) THEN
              t=-1.
              RETURN
        ENDIF
	
        Eb = ( -b - dsqrt(b**2-4.*a*c) )/2./a

c       IF (eb.LE.30.OR.eb.GE.40) THEN
        IF (eb.LE. 1..OR.eb.GE.300.) THEN ! corrected for compass

            t=-1.
            RETURN
        END IF

c       write(0,*) 'beam_calc: eb= ',eb

        pbeam4(4)= Eb
        pb1=DSQRT(Eb**2-ampic**2)
        pz = pb1/norm
        pbeam4(1)=axb*pz
        pbeam4(2)=ayb*pz
        pbeam4(3)=pz

c-- t-calculation

        ps = (pbeam4(1)-ptot4(1))**2 + (pbeam4(2)-ptot4(2))**2 + (pbeam4(3)-ptot4(3))**2

        t=(pbeam4(4)-ptot4(4))**2 - ps ! was error
        t=-t
c                                        Computation of Tmin
c   another formula !!!
c
        xmcms   =dsqrt(2.*eb*amprot+ampic**2+amprot**2)
        encms   =(xmcms**2+amrec**2-xm2)/2./xmcms
        pncms   =dsqrt(encms**2-amrec**2)
c
c---------  this I don't understand... but this works...
c------- this is Lorentz transformation along beam direction ...
c
        enlmin  =((eb+amprot)*encms-pb1*pncms)/xmcms
        tmin    =2.*enlmin*amprot-amprot**2-amrec**2
c
c       write(0,*) ' t= ',t,' tmin= ',tmin

        tprime      =t-tmin
	
        RETURN
        END
C---------------------------------------------------------------------
      SUBROUTINE RNDH(X, N)
C
C Generate quasi-random (super-uniform) distribution in N dimensions.
C Scrambled Halton J.H. sequence. Source:
C Braaten E and Weller G 1979 J. Comp. Phys. 33 249-58
C "An improved Low-Discrepancy Sequence for Multidimentional
C  Quasi-Monte Carlo Integration"
C Note that 1-st coord with prime=2 coincides with RNDL 1-st coord.
C Modifications:
C 1. Coord 2,3 were changed i --> p-i
C 2. Coord 5-16 were randomized via RNDL (Sobol seq) coord 2-13.
C
      IMPLICIT NONE
      INTEGER N, NDIM, NPERM, NDIG, NCALL, ISEED, INIT, init9
      PARAMETER (NDIM = 16, NPERM = 381, NDIG = 32)
      INTEGER PRIM(NDIM),LOCP(NDIM),PERM(NPERM),ND(NDIM),ID(NDIG,NDIM)
      REAL X(N), OP(NDIM), T
      DATA PRIM /
     +  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53/
      DATA LOCP /
     +  1,  3,  6, 11, 18, 29, 42, 59, 78,101,130,161,198,239,282,329/
      DATA PERM /
     &  0, 1, 0, 1, 2, 0, 2, 4, 1, 3, 0, 4, 2, 6, 1, 5, 3, 0, 6, 8,
     &  3, 7, 2, 4, 9,10, 5, 1, 0, 7, 4,10,11, 5, 8, 2, 9, 3,12, 6,
     &  1, 0, 9,13, 5, 3,11,15, 7, 6,14,10, 2, 8,16,12, 4, 1, 0,10,
     &  5,14,12, 3,16, 7, 4,13, 8,17,11, 2,15, 6, 9,18, 1, 0,12,17,
     &  6,20, 9, 3,14,10,21,16, 5,13, 2, 7,18, 8,19,11, 1,15, 4,22,
     &  0,15, 8,22,11,25, 4,18,20, 6,27,13,23, 9,16, 2, 7,21,14,28,
     & 10,24, 3,17,12, 1,19, 5,26, 0,16,23, 8,12,27,19, 4,25,10, 2,
     & 17,21, 6,14,29,22, 7,15,30,26,11, 3,18,13,28,20, 5, 1, 9,24,
     &  0,19,10,28, 5,23,14,32,34,16,25, 7,30,12,21, 3,13,31, 4,22,
     & 17,35, 8,26,24, 6,33,15,20, 2,29,11,18,36, 9,27, 1, 0,21,31,
     & 11,26, 6,16,36, 3,23,33,13,28, 8,18,38,37,17, 7,27,12,32,22,
     &  2,39,19, 9,29,14,34,24, 4,40,20,10,30,15,35,25, 5, 1, 0,22,
     & 11,32,37,16,27, 6,14,35, 3,24,29, 8,40,19, 7,28,18,39,33,12,
     & 23, 2,20,41,10,31,25, 4,36,15, 1,21,13,34,38,17,30, 9,26,42,
     &  5, 0,24,35,12, 6,29,41,18,32, 9,21,44,26, 3,15,38,36,13, 2,
     & 25,42,19, 8,31,22,45,34,11,16,39,28, 5,10,33,46,23, 7,30,40,
     & 17,27, 4,14,37,20, 1,43, 0,27,14,40,33, 7,46,20,23,49,10,36,
     & 43,17,30, 4,31, 5,44,18,12,38,25,51,48,22,35, 9,15,41, 2,28,
     & 42,16,29, 3,21,47, 8,34,39,13,52,26, 6,32,19,45,24,50, 1,37,
     & 11/
C---
      DATA NCALL /0/, INIT /0/
      SAVE NCALL, INIT, PRIM, LOCP, PERM, ND, ID, OP
      INTEGER LOC, I, J, K, IT
C
C NCALL   number of the current point to be generated
C PRIM(I) is the I-th prime number, OP(I) = 1./PRIM(I)
C ID(*,I) is the presentation of NCALL in PRIM(I) - adic system,
C ND(I)   is the number of digits in this presentation.
C LOCP(I) location in PERM() of permutations of digits for PRIM(I)
C PERM()  table of permutations of digits for all used prime numbers
C
      IF (INIT.NE.0) GOTO 100
	 IF (N.GT.NDIM) STOP 'RNDC: N TOO LARGE'
	 DO 10 I = 1, NDIM
	 OP(I)   = 1./PRIM(I)
         ND(I)   = 1
   10	 ID(1,I) = 0
	 INIT    = 1
  100 CONTINUE
C
C Calculate NCALL+1 in P-adic system: N1 N2 N3 ... NN
C Then invert the order of digits and form I-th quasirandom number as
C X(I) = 0.NN ... N3 N2 N1 for ordinary Halton sequence,
C X(I) = 0.PERM(NN,I) ... PERM(N3,I) PERM(N2,I) PERM(N1,I)
C for scrambled sequence.
C

      NCALL = NCALL + 1

c      if( init9.eq.0) then
c      WRITE(*,710) NCALL
c      endif

      DO 200 I = 1, N
	DO 110 J = 1, ND(I)
	  ID(J,I) = ID(J,I) + 1
	  IF (ID(J,I).LT.PRIM(I)) GOTO 120
	  ID(J,I) = 0
  110   CONTINUE
	ND(I) = ND(I) + 1
	ID(ND(I),I) = 1
C
  120	T = 0.
	LOC = LOCP(I)
	DO 130 J = ND(I), 1, -1
  130	T = ( T + PERM(LOC+ID(J,I)) ) * OP(I)
	X(I) = T

c	  if( init9.eq.0) then
c	  WRITE(*,720) I,X(I),PRIM(I),(ID(K,I),K=ND(I),1,-1)
C	  WRITE(*,720) I,X(I),PRIM(I),(PERM(LOC+ID(K,I)),K=ND(I),1,-1)
c	  endif

  200 CONTINUE
c	   init9 = 1
      RETURN
c
c 710 FORMAT(/' NCALL =',I4)
  710 FORMAT(/' NCALL =',I6)
  720 FORMAT(' X(',I2,') =',F9.6,' P =',I3,' ID =',32I3)

      ENTRY RNDHI
C
C Reinitialize generator.
C
      NCALL = 0
      RETURN

      ENTRY RNDHIN (ISEED)
C
C Restart generator from the given point number ISEED.
C Next generated point will be point number ISEED+1.
C
c===================================================
      IF (INIT.NE.0) GOTO 190

	  DO 19 I = 1, NDIM
	  OP(I)   = 1./PRIM(I)
          ND(I)   = 1
   19	  ID(1,I) = 0
	  INIT    = 1
  190 CONTINUE

      write(0,*) 'rndh: after init'
c===================================================
      NCALL = ISEED
      DO 300 I = 1, NDIM
	IT = NCALL
	DO 310 J = 1, NDIG
	  ND(I)  = J
	  ID(J,I)= MOD(IT,PRIM(I))

	  write(0,*) 'i=,j=',i,j,'it=',it,'prim=',prim(I),'div=',id(j,i)

	  IT = IT/PRIM(I)
	  IF (IT.EQ.0) GOTO 300
  310   CONTINUE
        WRITE(0,*) 'RNDHIN: too many digits for coord #',I
	STOP 'RNDHIN: given number is too large'
	
  300   CONTINUE

	INIT = 1

c       print*,'rndhin:ncall=',ncall
c       print*,'rndhin:ndim=',ndim
c       print*,'rndhin:ndig=',ndig
c       print*,'nd(1) =',nd(1)
c       print*,'id(j,1) =',id(1,1),id(2,1),id(3,1),id(4,1),id(5,1)

      RETURN

      ENTRY RNDHOU (ISEED)
C
C Write out the state of the generator - number of last generated
C point ISEED.
C
      ISEED = NCALL
      RETURN
      END

C-------------------------------
        SUBROUTINE SFO(EJ,PJ,P,XX1,XX2)
        DOUBLE PRECISION P(4),EJ,PJ
	PARAMETER  PI = 3.141 592 653 589
        PH=2.*PI*XX2
        CT=2.*XX1-1.
        ST=SQRT(1.-CT**2)
        P(1)=PJ*ST*COS(PH)
        P(2)=PJ*ST*SIN(PH)
        P(3)=PJ*CT
        P(4)=EJ
        RETURN
        END
C------------------
        SUBROUTINE AMR(A,SAM,AM)
        DOUBLE PRECISION A(4),AM,SAM
        SAM=A(4)**2-A(1)**2-A(2)**2-A(3)**2
        AM=SQRT(SAM)
        RETURN
        END
C------------------
        FUNCTION  FD1(X)
	real X
        FD1 = 1. + X*X
        RETURN
        END
C------------------

        real *8 FUNCTION FD1d(X)
	real *8 X
        FD1d = 1. + X*X
        RETURN
        END

c------------------------------------
        real *8 FUNCTION FD2d(X)
	real *8 X
        FD2d = 9. + 3.*X*X + X*X*X*X
        RETURN
        END
C------------------	
        FUNCTION FD2(X)
	real X
        FD2 = 9. + 3.*X*X + X*X*X*X
        RETURN
        END	
C------------------------------------
        FUNCTION FD3(X)
	real  X
        FD3 = 225. + 45.*X*X + 6.*X*X*X*X + X*X*X*X*X*X
        RETURN
        END
C------------------------------------
        real *8 FUNCTION  FD3d(X)
	real *8 X
        FD3d = 225. + 45.*X*X + 6.*X*X*X*X + X*X*X*X*X*X
        RETURN
        END
C---------------------------------------------------------
        FUNCTION FD4(X)
	real X
        FD4  =  ( x**4 - 45.*x**2 + 105) **2 + 25.*x**2*(2.*x**2 - 21. )**2  
        RETURN
        END
C------------------------------------
        real *8 FUNCTION  FD4d(X)
	real *8 X
        FD4d =  ( x**4 - 45.*x**2 + 105) **2 + 25.*x**2*(2.*x**2 - 21. )**2  
        RETURN
        END
C
c------------------------------------
	subroutine GJrotm(ppi,pn,an)
	double precision ppi(4),pn(4),an(3,3),zn,yn
	zn=dsqrt(ppi(1)**2+ppi(2)**2+ppi(3)**2)
	do i=1,3
	 an(3,i)=ppi(i)/zn
	enddo
	an(2,1)=pn(2)*ppi(3)-pn(3)*ppi(2)
	an(2,2)=pn(3)*ppi(1)-pn(1)*ppi(3)
	an(2,3)=pn(1)*ppi(2)-pn(2)*ppi(1)
	yn=dsqrt(an(2,1)**2+an(2,2)**2+an(2,3)**2)
	do i=1,3
	 an(2,i)=an(2,i)/yn
	enddo
	an(1,1)=an(2,2)*an(3,3)-an(2,3)*an(3,2)
	an(1,2)=an(2,3)*an(3,1)-an(2,1)*an(3,3)
	an(1,3)=an(2,1)*an(3,2)-an(2,2)*an(3,1)
	return
	end

c--------------------------------------------------
          real function    psp_tab ( xm,  k   )
          real *8 xm 

          common / sampl_integrals_tab / ama_tab, amb_tab, stepx_tab,  ps_tab( 3000, 100 )
          real ama_tab, amb_tab, stepx_tab,  ps_tab
          common /normalize_amp/ inormal

          if(  inormal .eq. 0) then
            psp_tab  = 1.
            return
          endif
          
            if( xm.ge.ama_tab .and. xm.le.amb_tab ) then

                nb = ( xm - ama_tab - 0.5*stepx_tab  )/stepx_tab + 1 

                x1 = ama_tab+0.5*stepx_tab + (nb-1)*stepx_tab
                x2 = x1 + stepx_tab               

                y1 = ps_tab( nb,   k )
                y2 = ps_tab( nb+1, k )

                a  = (y2-y1)	 /(x2-x1)
                b  = (x2*y1-x1*y2)/(x2-x1)

                psp_tab = xm*a + b 

c              write(0,*)  ' xm = ', xm, ' nb = ', nb
c              write(0,*)  ' k = ', k, ' psp_tab = ', psp_tab

               if( psp_tab .eq. 0.)   psp_tab = 1. !  not to divide-by-zero

            else
 
            psp_tab = 1. !  not to divide-by-zero

            endif

          return
          end
c--------------------------------------
	subroutine scld(p1,p2,resd)
	double precision p1(4),p2(4),resd
	resd=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
	return
	end
cc-----------------------
	subroutine vnor(p1,res)
	double precision p1(4)
	res=sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
	return
	end
c------------------------------------
	subroutine vnord(p1,resd)
	double precision p1(4), resd
	resd=sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
	return
	end
c--------------------------
	subroutine  scln(p1,p2,res)
	double precision  p1(4),p2(4)
	call  scl(p1(1),p2(1),res1)
	call vnor(p1(1),an1)
	call vnor(p2(1),an2)
	res = res1/(an1*an2)
	return
	end
c------------------------------------------
	subroutine  sclnd(p1,p2,resd)
	double precision  p1(4),p2(4), resd, an1, an2, res1
	call  scld(p1(1),p2(1),res1)
	call vnord(p1(1),an1)
	call vnord(p2(1),an2)
	resd = res1/(an1*an2)
	return
	end

c----------------------------
	subroutine HSrotm(ppi,pn,an)
	double precision ppi(4),pn(4),an(3,3),zn,yn

	zn=dsqrt(pn(1)**2+pn(2)**2+pn(3)**2)

	do i=1,3
	 an(3,i)=-pn(i)/zn
	enddo

	an(2,1)=-(ppi(2)*pn(3)-ppi(3)*pn(2))
	an(2,2)=-(ppi(3)*pn(1)-ppi(1)*pn(3))
	an(2,3)=-(ppi(1)*pn(2)-ppi(2)*pn(1))

	yn=dsqrt(an(2,1)**2+an(2,2)**2+an(2,3)**2)

	do i=1,3
	 an(2,i)=an(2,i)/yn
	enddo

	an(1,1)=an(2,2)*an(3,3)-an(2,3)*an(3,2)
	an(1,2)=an(2,3)*an(3,1)-an(2,1)*an(3,3)
	an(1,3)=an(2,1)*an(3,2)-an(2,2)*an(3,1)

	return
	end

C------------------------
	subroutine cmpair(i,j,amij,pi)
	common /jack/ pgj(4,4)
	double precision pgj,pi(4),AMIJ
	do k=1,3
	PI(K)=0.5*(PGJ(K,I)-PGJ(K,J)+(PGJ(K,I)+PGJ(K,J))*
     *   (PGJ(4,J)-PGJ(4,I))/(AMIJ+PGJ(4,I)+PGJ(4,J)))
	enddo
	RETURN
	END

C------------------------------------------------
        COMPLEX *16  FUNCTION BGEXP_POW_D(AM, par_all)
        real *8  par_all(50) 
        real *8 am, am0, amn,  alf, bet

        am0 = par_all(1)
        
        alf = par_all(2)
        bet = par_all(3)
        
        
        amn = am0 + 0.5
        
        if( am.gt. am0 ) then
          BGEXP_POW_D =  ( ( am - am0 )/(amn - am0) )**alf * dexp( - bet*(am-amn) )
        else
          BGEXP_POW_D = (0.,0.)
        endif
        
c        write(0,*) ' am= ', am, 'BGEXP_POW_D = ',  BGEXP_POW_D 
        
        return
        end

c------------------------------------------------------------------------------------
	subroutine w4ns( a,b,c,d,res )
	real *8 a(4),b(4),c(4),d(4)
	real ab,ac,ad,bc,bd,cd, res(9)
c-------- order of projections: 0, 1x, 1y, 2x, 2y, 3x, 3y, 4x, 4y	

	call scl( a(1),b(1),ab )
	call scl( a(1),c(1),ac )
	call scl( a(1),d(1),ad )
	call scl( b(1),c(1),bc )
	call scl( b(1),d(1),bd )
	call scl( c(1),d(1),cd )


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 zw0n =
     -   ad*bc/35 + ac*bd/35 + ab*cd/35 - cd*a(3)*b(3)/7 - 
     -   bd*a(3)*c(3)/7 - ad*b(3)*c(3)/7 - bc*a(3)*d(3)/7 - 
     -   ac*b(3)*d(3)/7 - ab*c(3)*d(3)/7 + a(3)*b(3)*c(3)*d(3)

	zw0n = zw0n*24.*sqrt(35./8.)
	res(1) = zw0n

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	zw1n = 
     -  -cd*a(3)*b(1)/7 - cd*a(1)*b(3)/7 - bd*a(3)*c(1)/7 - 
     -   ad*b(3)*c(1)/7 - bd*a(1)*c(3)/7 - ad*b(1)*c(3)/7 - 
     -   bc*a(3)*d(1)/7 - ac*b(3)*d(1)/7 - ab*c(3)*d(1)/7 + 
     -   a(3)*b(3)*c(3)*d(1) - bc*a(1)*d(3)/7 - ac*b(1)*d(3)/7 - 
     -   ab*c(1)*d(3)/7 + a(3)*b(3)*c(1)*d(3) + a(3)*b(1)*c(3)*d(3) + 
     -   a(1)*b(3)*c(3)*d(3)
	   zw1n = zw1n*6.*sqrt(7.)

	res(2) = zw1n

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          zw1p = 
     -  -cd*a(3)*b(2)/7 - cd*a(2)*b(3)/7 - bd*a(3)*c(2)/7 - 
     -   ad*b(3)*c(2)/7 - bd*a(2)*c(3)/7 - ad*b(2)*c(3)/7 - 
     -   bc*a(3)*d(2)/7 - ac*b(3)*d(2)/7 - ab*c(3)*d(2)/7 + 
     -   a(3)*b(3)*c(3)*d(2) - bc*a(2)*d(3)/7 - ac*b(2)*d(3)/7 - 
     -   ab*c(2)*d(3)/7 + a(3)*b(3)*c(2)*d(3) + a(3)*b(2)*c(3)*d(3) + 
     -   a(2)*b(3)*c(3)*d(3)


	   zw1p = zw1p*6.*sqrt(7.)
	res(3) = zw1p

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	zw2n = 
     -  -cd*a(2)*b(1)/7 - cd*a(1)*b(2)/7 - bd*a(2)*c(1)/7 - 
     -   ad*b(2)*c(1)/7 - bd*a(1)*c(2)/7 - ad*b(1)*c(2)/7 - 
     -   bc*a(2)*d(1)/7 - ac*b(2)*d(1)/7 - ab*c(2)*d(1)/7 + 
     -   a(3)*b(3)*c(2)*d(1) + a(3)*b(2)*c(3)*d(1) + 
     -   a(2)*b(3)*c(3)*d(1) - bc*a(1)*d(2)/7 - ac*b(1)*d(2)/7 - 
     -   ab*c(1)*d(2)/7 + a(3)*b(3)*c(1)*d(2) + a(3)*b(1)*c(3)*d(2) + 
     -   a(1)*b(3)*c(3)*d(2) + a(3)*b(2)*c(1)*d(3) + 
     -   a(2)*b(3)*c(1)*d(3) + a(3)*b(1)*c(2)*d(3) + 
     -   a(1)*b(3)*c(2)*d(3) + a(2)*b(1)*c(3)*d(3) + a(1)*b(2)*c(3)*d(3)
	
	zw2n = zw2n*sqrt(7./2.)*4.

	res(4) = zw2n



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          zw2p=
     -  -cd*a(1)*b(1)/7 + cd*a(2)*b(2)/7 - bd*a(1)*c(1)/7 - 
     -   ad*b(1)*c(1)/7 + bd*a(2)*c(2)/7 + ad*b(2)*c(2)/7 - 
     -   bc*a(1)*d(1)/7 - ac*b(1)*d(1)/7 - ab*c(1)*d(1)/7 + 
     -   a(3)*b(3)*c(1)*d(1) + a(3)*b(1)*c(3)*d(1) + 
     -   a(1)*b(3)*c(3)*d(1) + bc*a(2)*d(2)/7 + ac*b(2)*d(2)/7 + 
     -   ab*c(2)*d(2)/7 - a(3)*b(3)*c(2)*d(2) - a(3)*b(2)*c(3)*d(2) - 
     -   a(2)*b(3)*c(3)*d(2) + a(3)*b(1)*c(1)*d(3) + 
     -   a(1)*b(3)*c(1)*d(3) - a(3)*b(2)*c(2)*d(3) - 
     -   a(2)*b(3)*c(2)*d(3) + a(1)*b(1)*c(3)*d(3) - a(2)*b(2)*c(3)*d(3)

	zw2p = zw2p*4.*sqrt(7./2.)

	res(5) = zw2p

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	zw3n = 

     -  -a(3)*b(2)*c(1)*d(1) - a(2)*b(3)*c(1)*d(1) - 
     -   a(3)*b(1)*c(2)*d(1) - a(1)*b(3)*c(2)*d(1) - 
     -   a(2)*b(1)*c(3)*d(1) - a(1)*b(2)*c(3)*d(1) - 
     -   a(3)*b(1)*c(1)*d(2) - a(1)*b(3)*c(1)*d(2) + 
     -   a(3)*b(2)*c(2)*d(2) + a(2)*b(3)*c(2)*d(2) - 
     -   a(1)*b(1)*c(3)*d(2) + a(2)*b(2)*c(3)*d(2) - 
     -   a(2)*b(1)*c(1)*d(3) - a(1)*b(2)*c(1)*d(3) - 
     -   a(1)*b(1)*c(2)*d(3) + a(2)*b(2)*c(2)*d(3)

	zw3n = 6.*zw3n

	res(6) = zw3n
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	zw3p = 
	
     -  a(3)*b(1)*c(1)*d(1) + a(1)*b(3)*c(1)*d(1) - 
     -   a(3)*b(2)*c(2)*d(1) - a(2)*b(3)*c(2)*d(1) + 
     -   a(1)*b(1)*c(3)*d(1) - a(2)*b(2)*c(3)*d(1) - 
     -   a(3)*b(2)*c(1)*d(2) - a(2)*b(3)*c(1)*d(2) - 
     -   a(3)*b(1)*c(2)*d(2) - a(1)*b(3)*c(2)*d(2) - 
     -   a(2)*b(1)*c(3)*d(2) - a(1)*b(2)*c(3)*d(2) + 
     -   a(1)*b(1)*c(1)*d(3) - a(2)*b(2)*c(1)*d(3) - 
     -   a(2)*b(1)*c(2)*d(3) - a(1)*b(2)*c(2)*d(3)

	zw3p = 6.*zw3p

	res(7) = zw3p
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	zw4n = 
	
     -  a(2)*b(1)*c(1)*d(1) + a(1)*b(2)*c(1)*d(1) + 
     -   a(1)*b(1)*c(2)*d(1) - a(2)*b(2)*c(2)*d(1) + 
     -   a(1)*b(1)*c(1)*d(2) - a(2)*b(2)*c(1)*d(2) - 
     -   a(2)*b(1)*c(2)*d(2) - a(1)*b(2)*c(2)*d(2)


	zw4n = zw4n*6.*sqrt(2.)

	res(8) = zw4n

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	zw4p = 
	
     -  a(1)*b(1)*c(1)*d(1) - a(2)*b(2)*c(1)*d(1) - 
     -   a(2)*b(1)*c(2)*d(1) - a(1)*b(2)*c(2)*d(1) - 
     -   a(2)*b(1)*c(1)*d(2) - a(1)*b(2)*c(1)*d(2) - 
     -   a(1)*b(1)*c(2)*d(2) + a(2)*b(2)*c(2)*d(2)

	zw4p = zw4p*24.*sqrt(1./8.)
	
	res(9) = zw4p

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 return
	 end


ccccccccccccccccccccccccccccccccccccccccc

	subroutine W3ns(a,b,c, res)
	double precision a(4),b(4),c(4)
	real  ab,ac,bc,res(7) 	
	
c-------- order of projections: 0, 1x, 1y, 2x, 2y, 3x, 3y

	call scl( a(1),b(1),ab )
	call scl( a(1),c(1),ac )
	call scl( b(1),c(1),bc )

	res(1) = a(3)*b(3)*c(3) - 1./5.*(c(3)*ab+b(3)*ac+a(3)*bc)

	res(1) = res(1)*5./2.*1.*1.

	res(2) =   a(1)*b(3)*c(3) + a(3)*b(1)*c(3) + 
     *             a(3)*b(3)*c(1)  - 
     *            1./5.*( c(1)*ab + b(1)*ac + a(1)*bc ) 

	res(2) = res(2)*5./2.*sqrt(1./12.)*sqrt(2.)



	res(3) =  a(2)*b(3)*c(3) + a(3)*b(2)*c(3) + 
     *             a(3)*b(3)*c(2)  - 
     *            1./5.*( c(2)*ab + b(2)*ac + a(2)*bc ) 

	res(3) = res(3)*5./2.*sqrt(1./12.)*sqrt(2.)


	res(4) = a(1)*b(1)*c(3) - a(2)*b(2)*c(3) + 
     *	         b(1)*c(1)*a(3) - b(2)*c(2)*a(3) + 
     *           a(1)*c(1)*b(3) - a(2)*c(2)*b(3)

c	res(4) = res(4)*15.*sqrt(1./120.)*sqrt(2.)
	res(4) = res(4)*5.*sqrt(1./120.)*sqrt(2.)


	res(5) =  a(1)*b(2)*c(3) + b(1)*a(2)*c(3) + 
     *	          b(1)*c(2)*a(3) + c(1)*b(2)*a(3) + 
     *            a(1)*c(2)*b(3) + c(1)*a(2)*b(3)

c	res(5) = res(5)*15.*sqrt(1./120.)*sqrt(2.)
	res(5) = res(5)*5.*sqrt(1./120.)*sqrt(2.)


	res(6) =  a(1)*b(1)*c(1) -
     *	          b(2)*c(2)*a(1) -
     *            a(2)*c(2)*b(1) - 
     *            a(2)*b(2)*c(1)
	res(6) = res(6)*15.*sqrt(1./720.)*sqrt(2.)

	res(7) = -a(2)*b(2)*c(2) +
     *	          b(1)*c(1)*a(2) +
     *            a(1)*c(1)*b(2) + 
     *            a(1)*b(1)*c(2)

	res(7) = res(7)*15.*sqrt(1./720.)*sqrt(2.)

	return
	end

C----------------------------------------------

	subroutine W2ns(a,b,res)!T2(AB)
	real *8 a(4),b(4)
	real ab, res(5) 	

c-------- order of projections: 0, 1x, 1y, 2x, 2y

ccc	call scl( a(1),b(1),ab )


cc	res(1) = a(3)*b(3)-1./3.*ab
	res(1) = a(3)*b(3)*2./3. - a(1)*b(1)*1./3. - a(2)*b(2)*1./3.
	res(1) = res(1)*3./2.*1.*1.

	res(2) = a(1)*b(3)+b(1)*a(3)
	res(2) = res(2)*3./2.*sqrt(1./6.)*sqrt(2.)

	res(3) = a(2)*b(3)+b(2)*a(3)
	res(3) = res(3)*3./2.*sqrt(1./6.)*sqrt(2.)

	res(4) = a(1)*b(1) - a(2)*b(2)
	res(4) = res(4)*3.*sqrt(1./24.)*sqrt(2.)


	res(5) = a(1)*b(2) + a(2)*b(1)
	res(5) = res(5)*3.*sqrt(1./24.)*sqrt(2.)


	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccc	

	subroutine W1ns(a,res)
	real *8 a(4)
	real res(3) 	

c-------- order of projections: 0, 1x, 1y
	res(1) = a(3)
	res(2) = a(1)
	res(3) = a(2)
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      VARIOUS CASES OF SPINS  SUMMING      c
cccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_000(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res
	
	call scl( a(1),b(1),res )

	return
	end
	
cccccccccccccccccccccccccccccccccccccccccc	

	subroutine JLS_101(c,a,b, res )
	real *8 a(4),b(4),c(4), q(4)
	real res(3)
	
	call vec( a(1),b(1), q(1) )
	 
	 res(1) = q(3)
	 res(2) = q(1)
	 res(3) = q(2)
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccc	

	subroutine JLS_202(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(5)
	
	call W2ns( a(1),b(1),res(1))

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_011(c,a,b, res )
	real *8 a(4),b(4),c(4), q(4)
	real res
	
	call vec( a(1),b(1), q(1) )
	call scl( c(1),q(1),res )

	return
	end	
		
cccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_110(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(3)

	call scl( a(1),b(1),ab )
	
	res(1) = ab*c(3)
	res(2) = ab*c(1)
	res(3) = ab*c(2)
		
	return
	end
	
cccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JLS_111(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(3)

	call scl( a(1),c(1),ac )
	call scl( b(1),c(1),bc )
	

	res(1) = a(3)*bc - b(3)*ac
	res(2) = a(1)*bc - b(1)*ac
	res(3) = a(2)*bc - b(2)*ac
	
	
	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JLS_112(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(3)

	call scl( a(1),c(1),ac )
	call scl( b(1),c(1),bc )
	call scl( a(1),b(1),ab )
	
	res(1) = a(3)*bc + b(3)*ac - 2./3.*c(3)*ab
	res(2) = a(1)*bc + b(1)*ac - 2./3.*c(1)*ab
	res(3) = a(2)*bc + b(2)*ac - 2./3.*c(2)*ab
	

	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_211(c,a,b, res )
	real *8 a(4),b(4),c(4), q(4)
	real res(5)

	call vec( a(1),b(1), q(1) )

	call W2ns( c(1),q(1),res(1))	
	
	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JLS_212(c,a,b, res )
	real *8 a(4),b(4),c(4), q1(4),q2(4)
	real res(5),res1(5),res2(5)

	call vec( a(1),c(1), q1(1) )
	call vec( b(1),c(1), q2(1) )

	call W2ns( a(1),q2(1),res1(1))	
	call W2ns( b(1),q1(1),res2(1))	
	
	do k = 1,5
	 res(k) = res1(k) + res2(k)  
	enddo
	
	return
	end
	

ccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_312(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(7)

	call W3ns( c(1),a(1),b(1),res(1))	
	
	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_022(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res

	call scl( a(1),b(1),ab )
	call scl( a(1),c(1),ac )
	call scl( b(1),c(1),bc )
	call scl( c(1),c(1),cc )   ! 06-sep-2000 D.R.

	res = ac*bc - 1./3.*ab*cc

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_121(c,a,b, res )
	real *8 a(4),b(4),c(4),q(4)
	real res(3),qc

        call vec( a(1),b(1), q(1) )
	
	call scl( q(1),c(1),qc )
	call scl( c(1),c(1),cc )

	res(1) = c(3)*qc - q(3)*cc/3.
	res(2) = c(1)*qc - q(1)*cc/3.
	res(3) = c(2)*qc - q(2)*cc/3.
	
	return
	end	

ccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_122(c,a,b, res )
	real *8 a(4),b(4),c(4),q1(4),q2(4)
	real res(3)

        call vec( a(1),c(1), q1(1) )
        call vec( b(1),c(1), q2(1) )
	
	call scl( a(1),c(1),ac )
	call scl( b(1),c(1),bc )

	res(1) = bc*q1(3) + ac*q2(3)	
	res(2) = bc*q1(1) + ac*q2(1)	
	res(3) = bc*q1(2) + ac*q2(2)	

	return
	end
	
cccccccccccccccccccccccccccccccccccccccccc	

	subroutine JLS_220(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(5),res1(5)

	call scl( a(1),b(1),ab )
	
	call W2ns( c(1),c(1),res1(1))

	   do k = 1,5
	    res(k) = res1(k)*ab	
	   enddo
	  
	return
	end
	

ccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JLS_221(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(5),res1(5),res2(5)

	call scl( a(1),c(1),ac )
	call scl( b(1),c(1),bc )
	
	call W2ns( a(1),c(1),res1(1))
	call W2ns( b(1),c(1),res2(1))

	do k = 1,5
	 res(k) = res1(k)*bc - res2(k)*ac    
	enddo
	  
	return
	end
	

ccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JLS_222(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(5),res1(5),res2(5),res3(5),res4(5)

	call scl( a(1),b(1),ab )
	call scl( a(1),c(1),ac )
	call scl( b(1),c(1),bc )
	call scl( c(1),c(1),cc )
	
	
	call W2ns( a(1),b(1),res1(1))
	call W2ns( a(1),c(1),res2(1))
	call W2ns( b(1),c(1),res3(1))
	call W2ns( c(1),c(1),res4(1))

	do k = 1,5
	res(k) = -2./3.*res1(k)*cc + res2(k)*bc + res3(k)*ac  -2./3.*res4(k)*ab  ! ERROR fixed !
	 
	enddo
	  
	return
	end
	

ccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JLS_321(c,a,b, res )
	real *8 a(4),b(4),c(4),q(4)
	real res(7)

	call vec( a(1),b(1),q(1) )
	
	call W3ns( c(1),c(1),q(1),res(1) )
	
	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JLS_322(c,a,b, res )
	real *8 a(4),b(4),c(4),q1(4),q2(4)
	real res(7),res1(7),res2(7)
	

	call vec( c(1),a(1),q1(1) )
	call vec( c(1),b(1),q2(1) )
	
	call W3ns( c(1),b(1),q1(1),res1(1) )
	call W3ns( c(1),a(1),q2(1),res2(1) )
	
	do k = 1,7
	  res(k) = res1(k) + res2(k)
	enddo  

	
	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JLS_422(c,a,b, res )
	real *8 a(4),b(4),c(4)
	real res(9)

	call W4ns( c(1),c(1),a(1),b(1), res(1) )
	
	return
	end
	
	
cccccccccccccccccccccccccccccccccccccccccccccccc
c     V+ P wavefunctions       c
cccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JL_10(c,a, res )
	real *8 a(4),c(4)
	real res(3)
	
	res(1) = a(3)
	res(2) = a(1)
	res(3) = a(2)
	
	return
	end	


cccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JL_01(c,a, res )
	real *8 a(4),c(4)
	real res

	call scl( a(1),c(1),res )
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JL_11(c,a, res )
	real *8 a(4),c(4), q(4)
	real res(3)
	
	call vec( a(1),c(1),q(1) )
	
	res(1) = q(3)
	res(2) = q(1)
	res(3) = q(2)
      
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine JL_21(c,a, res )
	real *8 a(4),c(4)
	real res(5)
	
	call W2ns( a(1),c(1), res(1) )      
	 
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JL_12(c,a, res )
	real *8 a(4),c(4)
	real res(3)
	
	call scl( a(1),c(1),ac )
	call scl( c(1),c(1),cc )
	
	res(1) = c(3)*ac - a(3)*cc/3.
	res(2) = c(1)*ac - a(1)*cc/3.
	res(3) = c(2)*ac - a(2)*cc/3.

c	res(1) = a(3)
c	res(2) = a(1)
c	res(3) = a(2)
	
	return
	end	


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JL_22(c,a, res )
	real *8 a(4),c(4),q(4)
	real res(5)
	
	call vec( a(1),c(1),q(1))
	
	call W2ns( c(1),q(1), res(1) )
	
	RETURN
	END
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JL_32(c,a, res )
	real *8 a(4),c(4)
	real res(7)

	call W3ns( c(1),c(1), a(1), res(1) )
	
	RETURN
	END
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine JL_33(c,a, res )
	real *8 a(4),c(4),q(4)
	real res(7)
	
	call vec( a(1),c(1),q(1))
	call W3ns( c(1),c(1), q(1), res(1) )
	
	RETURN
	END
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      subroutine gener_cPPX(P1, P2, Q1, Q2, PX, WeightD)
C                                                    
C                Sdv  --- 27.09.2013 ---
C
C     Correction of the gener_cp2 for ALICE experiment 
C     The PP --> X process            
C
C     Q1, Q2 - 5-momentums of the secondary protons         
C...............................................,,,,,......
      
      real *8 P1(6), P2(6), Q1(6), Q2(6), PX(6)
c      real *8 P1(6), P2(6), Q1(6), Q2(6)
      real *8 xmd, Weightd
      real *8 pbeam, ambeam, am1, am2, t1, t2, x1, x2, amtarg, amx, s, sqs, s1, s2 
      real *8 dams2min, dams2pl, sdams2pl, t1q, t2q
      real *8 traj_p1, traj_r1, traj_p, traj_r, FUNCBD

      COMMON/RLRCOM3/ PLAB(4,30)
      real *8 PLAB
      real *8 pcm(4), pcm1(4), pfastcm(4), pslowcm(4), pxcm(4), dxfeyn

      common /cm_vectors/  pa(4),pb(4),pa1(4),pb1(4),px0(4) !'px0' name used to avoid double naming
      real *8 pa,pb,pa1,pb1,px0

c----------- c.m  of cent. produced system  

      real *8 pbeam1(4), pbeam2(4), pres(4), pbeam1_cm(4), pbeam2_cm(4), pslow_cm(4), pfast_cm(4)
      real *8 pbeam_cm(4), ptarg_cm(4)
      real *8 ygj1(4), ygj2(4)
      real *8 ctd

c----- cuts on CP variables, independent of the channel
      common /y_cent_cut/ yminc, ymaxc
      common /x_fast_slow_cuts/ x1min, x1max, x2min, x2max
      common /t_fast_slow_cuts/ t1min, t1max, t2min, t2max

      real *8 p_pi, amx2, wpipi
      real *8 q2a, q2b
      real *8 Qmx1x2

      parameter (ampic =0.13956755)

c--------- begin of CP-generator  variables  commons, read from file with ascii cards...

      common /pbeam_f/ pbeamf,  sigma_pbeam     
      common /itype_beam/  itypb1 
      common /mass_limits/ xma, xmb
      common /number_of_mc_events/ nmc1
      common /mode_reject/ imoderej
      common /fill_fluxes/ ihist_fluxes

c--------- end of CP-generator  variables  commons ---------

      common /fast_geant_code/ icodegeant_fast 

      real *8 pbm(4),ax,ay
      
      real x1r, x2r

      common /vertex_simulated/ vertex(3)

      real X(10)

      common /fast_mode/ ifast_mode     
      common /init_fast/ init6 
      common /z_limits/ z1min, z1max, z2min, z2max, ztotmin, ztotmax, ntot_z
      
      data yminc, ymaxc               /-5.0, 5.0  /
      data x1min, x1max, x2min, x2max / 0.0, 1.0, 0.0, 1.0 /
      data t1min, t1max, t2min, t2max / 0.0, 5.0, 0.0, 5.0 /

c--------- copy from "cards" values 

C-    write(*,*) ' gener_cp2S:'

      if(  ifast_mode .ne. 0 .and. init6 .ne. 0)  goto 779

      init6 = 1
      icodegeant_fast = 14            ! proton_beam GEANT code

C=     write(0,*) ' itypb1 = ', itypb1, pbeam

       am_min = xma		      ! mass limits
       am_max = xmb		      !
c------------------------------------------------------------------
       nmes   = 0                     ! until gener is called and nmes is defined...

       ambeam = p1(5)                 ! Sdv: proton mass as beam mass 
       amtarg = p2(5)                 ! Sdv: proton mass as target mass 

       do i = 1,4
          pa (i) = p1(i)
	  pb (i) = p2(i)
          pcm(i) = p1(i) +  p2(i) 
       enddo

       s = pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2
       sqs = dsqrt( s )
C
C  --- calculation of the CM beam momenta ---
C

C-      write(*,*) ' pa= ', pa, dsqrt(pa(4)**2-pa(3)**2)
C-      write(*,*) ' pb= ', pb, dsqrt(pb(4)**2-pb(3)**2)
C
C--------------------------------------------------------

        zmin  =  0. ! defines max of abs(1-x), NO CUT here ...
        zmax  =  6.
        zmax1 =  6.

CSdv+
C-      z1min_set   = 1.90E-06     !  Sdv: we have here uncommented the above original statements
C-	z2min_set   = 1.90E-06     
C
        z1min_set   = 0.
	z2min_set   = 0.
c	
	z1max_set   = 10.00
	z2max_set   = 10.00
	
	ztotmin_set =  0.00        !  Sdv:  ztotmin_set and ztotmax_set are not clear at the moment
	ztotmax_set = 20.00
CSdv-	
 779    continue
  55    continue

  57    CALL PP2RANLUX(X(1),2)        !  random generator
	
        ntot_z = ntot_z  + 1
        z1     = z1min_set + (z1max_set - z1min_set) * X(1)
        z2     = z2min_set + (z2max_set - z2min_set) * X(2)
	ztot   = z1+z2
CSdv+
C-      write(*,*) ' gener_cp2: ztot, ztotmin_set, ztotmax_set =', ztot, ztotmin_set, ztotmax_set

        if( ztot .gt. ztotmax_set) goto 55
        if( ztot .lt. ztotmin_set) goto 55
c----------------------------------------
        x1 = 1. - exp(-z1)
        x2 = 1. - exp(-z2)
C	
	x1r= x1
	x2r= x2
c$$$	call hf1(11, x1r, 1.)
c$$$	call hf1(12, x2r, 1.)
c----------------------------------------
CSdv+
C-      write(*,*) ' gener_cp2: x1, x2 =', x1, x2, z1, z2
	
        b0a = 0.15
        b0b = 0.15

        call scatter_pomeron( pa(1), pa1(1), b0a, x1, ambeam)
        call scatter_pomeron( pb(1), pb1(1), b0b, x2, amtarg)

        do i  = 1,4
	q1(i) = pa1(i)
	q2(i) = pb1(i)
        px(i) = pa(i)+pb(i)-pa1(i)-pb1(i)
        px0(i)=px(i)
        enddo
	amx2  = px(4)**2-px(1)**2-px(2)**2-px(3)**2
	if (amx2.lt.0.0) go to 57 
	amx   = dsqrt(amx2)
	px(5) = amx
C
CSdv+	
	Qmx1x2= amx2/(s*(1.-x1)*(1.-x2))      !   Correct values 1-x1 and  1-x2  
C-      write(*,*) 'amx = ', amx, ' am_min, am_max = ',  am_min, am_max, Qmx1x2

CSdv-   if( amx.le.am_min .or. amx.ge.am_max ) goto 55

        if( z1 .lt. z1min ) z1min = z1
        if( z2 .lt. z2min ) z2min = z2

        if( z1 .gt. z1max ) z1max = z1
        if( z2 .gt. z2max ) z2max = z2

        if( ztot .lt. ztotmin ) ztotmin = ztot
        if( ztot .gt. ztotmax ) ztotmax = ztot

c-------- needed in "fast" mode !!! ---------------------------------
        xmd = dsqrt( amx2 )

        t1  = dams2min(   pa1(1), pa(1) )
        t2  = dams2min(   pb1(1), pb(1) )

        q2a = pa1(1)**2 + pa1(2)**2  
        q2b = pb1(1)**2 + pb1(2)**2  

        s1  = sdams2pl(   pa1(1), px0(1) )
        s2  = sdams2pl(   pb1(1), px0(1) )

        t1q = t1
        t2q = t2

        t1f = dabs(t1)
        t2f = dabs(t2)

C-      goto 765
C-----------------------

CSdv+   if( ihist_fluxes .ne. 0) then

c            weight1 =  traj_p1(x1, t1q) / exp( -b0a * q2a) 
c            weight2 =  traj_p1(x2, t2q) / exp( -b0b * q2b) 

C -will not work as
c 
c           weight1 =  traj_p(s, s2, t1q) / exp( -b0a * q2a) *(1-x1)*s/s2
c           weight2 =  traj_p(s, s1, t2q) / exp( -b0b * q2b) *(1-x2)*s/s1
C
c            weight  =  weight1 * weight2
C
c           write(0,*) ' x1, x2   ', x1, x2
c           write(0,*) ' s1, s2 = ', s1, s2, s*(1.-x2), s*(1.-x1)
c           write(0,*) 'Weights=',weight, weight1, weight2

c            xtmp1 = 1.-x1
c            xtmp2 = 1.-x2

C-	    call hfill(61 , xtmp1 ,xdum , weight1)
C-	    call hfill(62 , xtmp2 ,xdum , weight2)
C-	    call hfill(71 , xtmp1 ,xdum , weight1)
C-	    call hfill(72 , xtmp2 ,xdum , weight2)

c$$$            call hfill(63 , xtmp1 , t1f , weight1)
c$$$            call hfill(64 , xtmp2 , t2f , weight2)
CSdv-       endif
c-----------------------
 765    continue

CSdv+
C-      write(*,*) ' t1f, t1min, t1max =', t1f, t1min, t1max
C-      write(*,*) ' t2f, t2min, t2max =', t2f, t2min, t2max
C-      write(*,*) ' x1 , x1min, x1max =', x1 , x1min, x1max
C-      write(*,*) ' x2 , x2min, x2max =', x2 , x2min, x2max
C
        if( t1f .lt. t1min .or. t1f .gt. t1max ) goto 55
        if( t2f .lt. t2min .or. t2f .gt. t2max ) goto 55

        if( x1  .lt. x1min .or. x1  .gt. x1max ) goto 55
        if( x2  .lt. x2min .or. x2  .gt. x2max ) goto 55

        y_fast = 0.5 * dlog( (pa1(4) + pa1(3))/(pa1(4) - pa1(3)) )
        y_slow = 0.5 * dlog( (pb1(4) + pb1(3))/(pb1(4) - pb1(3)) )
        y_cent = 0.5 * dlog( ( px(4) +  px(3))/( px(4)  - px(3)) )
C
C-      write(*,*) 'y_cent, yminc, ymaxc  =', y_cent, yminc, ymaxc
C
        if( y_cent .lt. yminc .or. y_cent .gt. ymaxc) goto 55

C 1: -------- Pomeron-Pomeron, using (1-x)^al(t) formula 
c
c       weightd =  traj_p1(x1, t1q)*traj_p1(x2, t2q) / dexp( -b0a*q2a - b0b * q2b)
c
c
C 2:  ------- Pomeron-Reggeon, symmetric, using (1-x)^al(t) formula 
c
c       weightd = ( traj_p1(x1, t1q)*traj_r1(x2, t2q)  +  traj_r1(x1, t1q)*traj_p1(x2, t2q) )/ exp( -b0a*q2a - b0b * q2b)
c

C 3: -------- Pomeron-Pomeron, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
C
       weightd =  traj_p(s, s2, t1q)*traj_p(s, s1, t2q) * (1.-x1)*(1.-x2) *s**2/(s1*s2)/ exp( -b0a*q2a - b0b * q2b)

C 4: -------- Pomeron-Reggeon, symmetric, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
c                             weight, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
c
c      weightd = ( traj_p(s,s2,t1q)*traj_r(s,s1,t2q) + traj_r(s,s2,t1q)*traj_p(s,s1,t2q) )* (1.-x1)*(1.-x2)*s**2/(s1*s2)/exp( -b0a*q2a - b0b * q2b)
c
C 5: -------- Using FUNCB(x1) from  ComGeant instead of 1/(1-|x|)  and   exp(-8bt )
c    
c      weightd =  FUNCBD(x1) *  FUNCBD(x2) * (1.-x1)*(1.-x2) * dexp( 7.1*(t1+t2) ) / dexp( -b0a*q2a - b0b * q2b)
C
        weightf = weightd

c$$$        call hf1(19,weightf,1.)
c$$$	call hf1(21,  x1r , 1.)
c$$$	call hf1(22,  x2r , 1.)
c
       return
       end  
C---------------------------------------------------------

      subroutine gener_cp2S(P1, P2, PX, Weightd, Wpipi)
C                                                    
C                Sdv  --- 16.07.2012 ---
C
C     Correction of the gener_cp2 for ALICE experiment                              
C......................................................
      
       real *8 P1(6), P2(6), PX(6)
       real *8 xmd, Weightd
       real *8 pbeam, ambeam, am1, am2, t1, t2, x1, x2, amtarg, amx, s, sqs, s1, s2 
       real *8 dams2min, dams2pl, sdams2pl, t1q, t2q
       real *8 traj_p1, traj_r1, traj_p, traj_r, FUNCBD

       COMMON/RLRCOM3/ PLAB(4,30)
       real *8 PLAB
       real *8 pcm(4), pcm1(4), pfastcm(4), pslowcm(4), pxcm(4), dxfeyn

       common /cm_vectors/  pa(4),pb(4),pa1(4),pb1(4), px0(4)
       real *8 pa,pb,pa1,pb1,px0

c----------- c.m  of cent. produced system  

      real *8 pbeam1(4), pbeam2(4), pres(4), pbeam1_cm(4), pbeam2_cm(4), pslow_cm(4), pfast_cm(4)
      real *8 pbeam_cm(4), ptarg_cm(4)
      real *8 ygj1(4), ygj2(4)
      real *8 ctd

c----- cuts on CP variables, independent of the channel
      common /y_cent_cut/ yminc, ymaxc
      common /x_fast_slow_cuts/ x1min, x1max, x2min, x2max
      common /t_fast_slow_cuts/ t1min, t1max, t2min, t2max

      real *8 p_pi, amx2, wpipi
      real *8 q2a, q2b
      real *8 Qmx1x2

      parameter (ampic =0.13956755)

c--------- begin of CP-generator  variables  commons, read from file with ascii cards...

      common /pbeam_f/ pbeamf,  sigma_pbeam     
      common /itype_beam/  itypb1 
      common /mass_limits/ xma, xmb
      common /number_of_mc_events/ nmc1
      common /mode_reject/ imoderej
      common /fill_fluxes/ ihist_fluxes

c--------- end of CP-generator  variables  commons ---------

      common /fast_geant_code/ icodegeant_fast 

      real *8 pbm(4),ax,ay
      
      real x1r, x2r

      common /vertex_simulated/ vertex(3)

      real X(10)

      common /fast_mode/ ifast_mode     
      common /init_fast/ init6 
      common /z_limits/ z1min, z1max, z2min, z2max, ztotmin, ztotmax, ntot_z
      
      data yminc, ymaxc /  -5.0, 5.0  /
      data x1min, x1max, x2min, x2max / 0.0, 1.0, 0.0, 1.0 /
      data t1min, t1max, t2min, t2max / 0.0, 5.0, 0.0, 5.0 /

c--------- copy from "cards" values 

C-    write(*,*) ' gener_cp2S:'

      if(  ifast_mode .ne. 0 .and. init6 .ne. 0)  goto 779

      init6 = 1
      icodegeant_fast = 14            ! proton_beam GEANT code

C=     write(0,*) ' itypb1 = ', itypb1, pbeam

       am_min = xma		      ! mass limits
       am_max = xmb		      !
c------------------------------------------------------------------
       nmes   = 0                     ! until gener is called and nmes is defined...

       ambeam = p1(5)                 ! Sdv: proton mass as beam mass 
       amtarg = p2(5)                 ! Sdv: proton mass as target mass 

       do i = 1,4
          pa (i) = p1(i)
	  pb (i) = p2(i)
          pcm(i) = p1(i) +  p2(i) 
       enddo

       s = pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2
       sqs = dsqrt( s )
C
C  --- calculation of the CM beam momenta ---
C

C-      write(*,*) ' pa= ', pa, dsqrt(pa(4)**2-pa(3)**2)
C-      write(*,*) ' pb= ', pb, dsqrt(pb(4)**2-pb(3)**2)
C
C--------------------------------------------------------

        zmin  =  0. ! defines max of abs(1-x), NO CUT here ...
        zmax  =  6.
        zmax1 =  6.

CSdv+
C-      z1min_set   = 1.90E-06     !  Sdv: we have here uncommented the above original statements
C-	z2min_set   = 1.90E-06     
C
        z1min_set   = 0.
	z2min_set   = 0.
c	
	z1max_set   = 10.00
	z2max_set   = 10.00
	
	ztotmin_set =  0.00        !  Sdv:  ztotmin_set and ztotmax_set are not clear at the moment
	ztotmax_set = 20.00
CSdv-	
 779    continue
  55    continue

  57    CALL PP2RANLUX(X(1),2)        !  random generator
	
        ntot_z = ntot_z  + 1
        z1     = z1min_set + (z1max_set - z1min_set) * X(1)
        z2     = z2min_set + (z2max_set - z2min_set) * X(2)
	ztot   = z1+z2
CSdv+
C-      write(*,*) ' gener_cp2: ztot, ztotmin_set, ztotmax_set =', ztot, ztotmin_set, ztotmax_set

        if( ztot .gt. ztotmax_set) goto 55
        if( ztot .lt. ztotmin_set) goto 55
c----------------------------------------
        x1 = 1. - exp(-z1)
        x2 = 1. - exp(-z2)
C	
	x1r= x1
	x2r= x2
c$$$	call hf1(11, x1r, 1.)
c$$$	call hf1(12, x2r, 1.)
c----------------------------------------
CSdv+
C-      write(*,*) ' gener_cp2: x1, x2 =', x1, x2, z1, z2
	
        b0a = 8.
        b0b = 8.

        call scatter_pomeron( pa(1), pa1(1), b0a, x1, ambeam)
        call scatter_pomeron( pb(1), pb1(1), b0b, x2, amtarg)

        do i = 1,4
        px(i) = pa(i)+pb(i)-pa1(i)-pb1(i)
        px0(i) =px(i)
        enddo
	amx2  = px(4)**2-px(1)**2-px(2)**2-px(3)**2
	if (amx2.lt.0.0) go to 57 
	amx   = dsqrt(amx2)
	px(5) = amx
CSdv+	
	Qmx1x2= amx2/(s*(1.-x1)*(1.-x2))      !   Correct values 1-x1 and  1-x2  
C-      write(*,*) 'amx = ', amx, ' am_min, am_max = ',  am_min, am_max, Qmx1x2

CSdv-   if( amx.le.am_min .or. amx.ge.am_max ) goto 55

        if( z1 .lt. z1min ) z1min = z1
        if( z2 .lt. z2min ) z2min = z2

        if( z1 .gt. z1max ) z1max = z1
        if( z2 .gt. z2max ) z2max = z2

        if( ztot .lt. ztotmin ) ztotmin = ztot
        if( ztot .gt. ztotmax ) ztotmax = ztot

c-------- needed in "fast" mode !!! ---------------------------------
        xmd = dsqrt( amx2 )

        t1  = dams2min(   pa1(1), pa(1) )
        t2  = dams2min(   pb1(1), pb(1) )

        q2a = pa1(1)**2 + pa1(2)**2  
        q2b = pb1(1)**2 + pb1(2)**2  

        s1  = sdams2pl(   pa1(1), px0(1) )
        s2  = sdams2pl(   pb1(1), px0(1) )

        t1q = t1
        t2q = t2

        t1f = dabs(t1)
        t2f = dabs(t2)

C-      goto 765
C-----------------------

CSdv+   if( ihist_fluxes .ne. 0) then

            weight1 =  traj_p1(x1, t1q) / exp( -b0a * q2a) 
            weight2 =  traj_p1(x2, t2q) / exp( -b0b * q2b) 

C -will not work as
c 
c           weight1 =  traj_p(s, s2, t1q) / exp( -b0a * q2a) *(1-x1)*s/s2
c           weight2 =  traj_p(s, s1, t2q) / exp( -b0b * q2b) *(1-x2)*s/s1
C
            weight  =  weight1 * weight2
C
c           write(0,*) ' x1, x2   ', x1, x2
c           write(0,*) ' s1, s2 = ', s1, s2, s*(1.-x2), s*(1.-x1)
c           write(0,*) 'Weights=',weight, weight1, weight2

            xtmp1 = 1.-x1
            xtmp2 = 1.-x2

C-	    call hfill(61 , xtmp1 ,xdum , weight1)
C-	    call hfill(62 , xtmp2 ,xdum , weight2)
C-	    call hfill(71 , xtmp1 ,xdum , weight1)
C-	    call hfill(72 , xtmp2 ,xdum , weight2)

c$$$            call hfill(63 , xtmp1 , t1f , weight1)
c$$$            call hfill(64 , xtmp2 , t2f , weight2)
CSdv-       endif
c-----------------------
 765    continue

CSdv+
C-      write(*,*) ' t1f, t1min, t1max =', t1f, t1min, t1max
C-      write(*,*) ' t2f, t2min, t2max =', t2f, t2min, t2max
C-      write(*,*) ' x1 , x1min, x1max =', x1 , x1min, x1max
C-      write(*,*) ' x2 , x2min, x2max =', x2 , x2min, x2max
C
        if( t1f .lt. t1min .or. t1f .gt. t1max ) goto 55
        if( t2f .lt. t2min .or. t2f .gt. t2max ) goto 55

        if( x1  .lt. x1min .or. x1  .gt. x1max ) goto 55
        if( x2  .lt. x2min .or. x2  .gt. x2max ) goto 55

        y_fast = 0.5 * dlog( (pa1(4) + pa1(3))/(pa1(4) - pa1(3)) )
        y_slow = 0.5 * dlog( (pb1(4) + pb1(3))/(pb1(4) - pb1(3)) )
        y_cent = 0.5 * dlog( ( px(4) +  px(3))/( px(4)  - px(3)) )
C
C-      write(*,*) 'y_cent, yminc, ymaxc  =', y_cent, yminc, ymaxc
C
        if( y_cent .lt. yminc .or. y_cent .gt. ymaxc) goto 55

c--------------------------------------------------------------------

	p_pi = dsqrt( dabs(amx2/4. - ampic**2) )
	wpipi=  p_pi/dsqrt(amx2)* dexp( -2.*p_pi**2 ) 

c-------- Pomeron-Pomeron, using (1-x)^al(t) formula 
c
c        weightd =  traj_p1(x1, t1q)*traj_p1(x2, t2q) / dexp( -b0a*q2a - b0b * q2b)
c

c ------- Pomeron-Reggeon, symmetric, using (1-x)^al(t) formula 
c
c        weightd = ( traj_p1(x1, t1q)*traj_r1(x2, t2q)  +  traj_r1(x1, t1q)*traj_p1(x2, t2q) )/ exp( -b0a*q2a - b0b * q2b)
c

C-------- Pomeron-Pomeron, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
C
         weightd =  traj_p(s, s2, t1q)*traj_p(s, s1, t2q) * (1.-x1)*(1.-x2) *s**2/(s1*s2)/ exp( -b0a*q2a - b0b * q2b)

c-------- Pomeron-Reggeon, symmetric, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
c                             weight, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
c
c        weightd = ( traj_p(s,s2,t1q)*traj_r(s,s1,t2q) + traj_r(s,s2,t1q)*traj_p(s,s1,t2q) )* (1.-x1)*(1.-x2)*s**2/(s1*s2)/exp( -b0a*q2a - b0b * q2b)
c
c        write(0,*) ' weight = ', weight
c        weightd = wpipi
c
c        Using FUNCB(x1) from  ComGeant instead of 1/(1-|x|)  and   exp(-8bt )
c    
c        weightd =  FUNCBD(x1) *  FUNCBD(x2) * (1.-x1)*(1.-x2) * dexp( 7.1*(t1+t2) ) / dexp( -b0a*q2a - b0b * q2b)
c
c
c  using PP --> pi+pi- cross-section if needed ...
c
c       weightd =  weightd * wpipi

        weightf = weightd

c$$$        call hf1(19,weightf,1.)
c$$$	call hf1(21, x1r, 1.)
c$$$	call hf1(22, x2r, 1.)

        if( ifast_mode .ne. 0) return	

c$$$        call hfill(65 , xtmp1 , t1f  ,weightf)
c$$$        call hfill(66 , xtmp2 , t2f  ,weightf)

C-      CALL BOOST1(sqs, pcm(1), pa1(1), plab(1,3))
C-      CALL BOOST1(sqs, pcm(1), pb1(1), plab(1,4))
C-      CALL BOOST1(sqs, pcm(1), px(1),  plab(1,7))

c-----     introduce phi-dependence if needed ------------------------------------

        do k4 = 1,4
          pbeam1(k4) = PLAB(k4,1) - PLAB(k4,3) 
          pbeam2(k4) = PLAB(k4,2) - PLAB(k4,4) 

            pres(k4) = PLAB(k4,7)              ! in case when decay is absent...
        enddo
        
        xmd = dsqrt( pres(4)**2 - pres(1)**2-  pres(2)**2- pres(3)**2  )
        
        CALL  BOOST2( xmd, pres(1), pbeam1(1), pbeam1_cm(1))
        CALL  BOOST2( xmd, pres(1), pbeam2(1), pbeam2_cm(1))

        CALL  BOOST2( xmd, pres(1), PLAB(1,3), pfast_cm(1) )
        CALL  BOOST2( xmd, pres(1), PLAB(1,4), pslow_cm(1) )

        CALL  BOOST2( xmd, pres(1), PLAB(1,1), pbeam_cm(1) )
        CALL  BOOST2( xmd, pres(1), PLAB(1,2), ptarg_cm(1) )

        call vec(pfast_cm(1),pbeam1_cm(1), ygj1(1)  )
        call vec(pslow_cm(1),pbeam2_cm(1), ygj2(1)  )

        ct_12 =  ctd( ygj1(1), ygj2(1) )
c       write(0,*) ' ct_12=', ct_12

        th_12 =   acos(ct_12)

c-------------  usual dphi between pt-s --------------------------------

C-       xx = plab(1,3)
C-       yy = plab(2,3)
C
	 xx = pa1(1)
	 yy = pa1(2)	
         phipt1 = atan2dmy(yy,xx)
	 
C-       xx =  plab(1,4)
C-       yy =  plab(2,4)
C         	 
	 xx = pb1(1)
	 yy = pb1(2)	
         phipt2 = atan2dmy(yy,xx)

         dphi   = phipt2-phipt1
         if(dphi .lt. -180.) dphi = dphi+360.
         if(dphi .gt.  180.) dphi = dphi-360.
         dphi = abs(dphi)

c-------------  for eta, eta' - mesons... ------------------------------------
cc          weight = weight * sin(dphi*3.141592/180.)**2  * dabs(t1*t2)


c------- introduce  decay of X and unisotropis angular distributions
c        relative to 2 reggeon collision axis and phi_GJ, in case of M .ne. 0
c
c        call cp_reweight_simple ( plab(1,1), plab(1,3), plab(1,4), plab(1,7), plab(1,5), weight3 )
c
c        weight = weight * weight3
c
       return
       end
       
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C-----------------------------------------------------
CSdv  --- 05.06.2012 ---
C
      subroutine gener_cp2(xmd,weightd)
       real *8 xmd, weightd
       real *8 pbeam, ambeam, am1, am2, t1, t2, x1, x2, amtarg, amx, s, sqs, s1, s2 
       real *8 dams2min, dams2pl, sdams2pl, t1q, t2q
       real *8 traj_p1, traj_r1, traj_p, traj_r, FUNCBD

       COMMON/RLRCOM3/ PLAB(4,30)
       real *8 PLAB
       real *8 pcm(4), pcm1(4), pfastcm(4), pslowcm(4), pxcm(4), dxfeyn

       common /cm_vectors/  pa(4),pb(4),pa1(4),pb1(4),px(4)
       real *8 pa,pb,pa1,pb1,px

c----------- c.m  of cent. produced system  

      real *8 pbeam1(4), pbeam2(4), pres(4), pbeam1_cm(4), pbeam2_cm(4), pslow_cm(4), pfast_cm(4)
      real *8 pbeam_cm(4), ptarg_cm(4)
      real *8 ygj1(4), ygj2(4)
      real *8 ctd

c----- cuts on CP variables, independent of the channel
      common /y_cent_cut/ yminc, ymaxc
      common /x_fast_slow_cuts/ x1min, x1max, x2min, x2max
      common /t_fast_slow_cuts/ t1min, t1max, t2min, t2max

      real *8     p_pi, amx2, wpipi
      real *8 q2a, q2b

      parameter (ampic =0.13956755)

c--------- begin of CP-generator  variables  commons, read from file with ascii cards...

      common /pbeam_f/ pbeamf,  sigma_pbeam     
      common /itype_beam/  itypb1 
      common /mass_limits/ xma, xmb
      common /number_of_mc_events/ nmc1
      common /mode_reject/ imoderej
      common /fill_fluxes/ ihist_fluxes

c--------- end of CP-generator  variables  commons ---------

      common /fast_geant_code/ icodegeant_fast 

      real *8 pbm(4),ax,ay

      common /vertex_simulated/ vertex(3)

      real X(10)

      common /fast_mode/ ifast_mode     
      common /init_fast/ init6 
      common /z_limits/ z1min, z1max, z2min, z2max, ztotmin, ztotmax, ntot_z
      
      data yminc, ymaxc /  -3.0, 3.0  /
      data x1min, x1max, x2min, x2max / 0.0, 1.0, 0.0, 1.0 /
      data t1min, t1max, t2min, t2max / 0.0, 1.0, 0.0, 1.0 /

c--------- copy from "cards" values 

      write(*,*) ' gener_cp2:'

      if(  ifast_mode .ne. 0 .and. init6 .ne. 0)  goto 779

      init6 = 1

      pbeam = pbeamf

C=     write(0,*) ' itypb1 = ', itypb1, pbeam

       if( itypb1 .eq. 1) then	      ! proton_beam
           ambeam = 0.938272013    
           icodegeant_fast = 14
       endif
C
C-     if( itypb1 .eq. 2) then	      ! piminus_beam
C-         ambeam = 0.1395
C-         icodegeant_fast = 9	 
C-     endif
C
       am_min = xma		      ! mass limits
       am_max = xmb		      !

c------------------------------------------------------------------
       nmes   = 0                     ! until gener is called and nmes is defined...

       ambeam = 0.938272013           ! Sdv: proton mass as beam mass 
       amtarg = 0.938272013           ! Sdv: proton mass as target mass 

       call mcbeam(pbm,vx,vy,vz,ax,ay)

       vertex(1) = vx
       vertex(2) = vy
       vertex(3) = vz

       do k = 1,3
       plab(k,1) = pbm(k)
       enddo
       plab(4,1) = dsqrt( ambeam**2 + pbm(1)**2+pbm(2)**2+pbm(3)**2  )

       plab(1,2) = 0.
       plab(2,2) = 0.
       plab(3,2) = 0.
       plab(4,2) = amtarg

       do i = 1,4
       pcm(i)=  plab(i,1) +  plab(i,2) 
       enddo

       s = pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2
       sqs = dsqrt( s )
C
C  --- calculation of the CM beam momenta ---
C
       CALL BOOST2(sqs, pcm(1), plab(1,1), pa(1) )	!  pa in CM system
       CALL BOOST2(sqs, pcm(1), plab(1,2), pb(1) )	!  pb in CM system
       
       write(*,*) ' pa= ', pa, dsqrt(pa(4)**2-pa(3)**2)
       write(*,*) ' pb= ', pb, dsqrt(pb(4)**2-pb(3)**2)
C
C--------------------------------------------------------

        zmin  =  0. ! defines max of abs(1-x), NO CUT here ...
        zmax  =  6.
        zmax1 =  6.

C --------- pion beam, eta-eta mass-limits, compass beam
C
C-      z1min_set   = 1.00E-07
C-      z2min_set   = 1.90E-06
C-
C-      z1max_set   = 6.00
C-      z2max_set   = 6.00
C-
C-      ztotmin_set = 2.36
C-      ztotmax_set = 7.00
C.........................................................
CSdv+
        z1min_set   = 1.90E-06     !  Sdv: we have here uncommented the above original statements
	z2min_set   = 1.90E-06     
	
	z1max_set   = 6.00
	z2max_set   = 6.00
	
	ztotmin_set = 2.36         !  Sdv:  ztotmin_set and ztotmax_set are not clear at the moment
	ztotmax_set = 7.00
CSdv-	
 779    continue

  55    CALL PP2RANLUX(X(1),2)        !  random generator
	
        ntot_z = ntot_z  + 1
        z1     = z1min_set + (z1max_set - z1min_set) * X(1)
        z2     = z2min_set + (z2max_set - z2min_set) * X(2)
	ztot   = z1+z2
CSdv+
C-      write(*,*) ' gener_cp2: ztot, ztotmin_set, ztotmax_set =', ztot, ztotmin_set, ztotmax_set

        if( ztot .gt. ztotmax_set) goto 55
        if( ztot .lt. ztotmin_set) goto 55
c----------------------------------------
        x1 = 1. - exp(-z1)
        x2 = 1. - exp(-z2)
c----------------------------------------
CSdv+
C-      write(*,*) ' gener_cp2: x1, x2 =', x1, x2
	
        b0a = 8.
        b0b = 8.

        call scatter_pomeron( pa(1), pa1(1), b0a, x1, ambeam )
        call scatter_pomeron( pb(1), pb1(1), b0b, x2, amtarg )

        do i = 1,4
          px(i) = pa(i)+pb(i)-pa1(i)-pb1(i)
        enddo

        amx2 = px(4)**2-px(1)**2-px(2)**2-px(3)**2

C-      write(0,*) ' amx2 = ', amx2, ' am_min, am_max = ',  am_min, am_max 

        if( amx2 .le. am_min**2 .or. amx2 .ge. am_max**2 ) goto 55

        if( z1 .lt. z1min ) z1min = z1
        if( z2 .lt. z2min ) z2min = z2

        if( z1 .gt. z1max ) z1max = z1
        if( z2 .gt. z2max ) z2max = z2

        if( ztot .lt. ztotmin ) ztotmin = ztot
        if( ztot .gt. ztotmax ) ztotmax = ztot
CSdv+
C-      write(*,*) ' amx2 = ', amx2, ' am_min**2, am_max**2 = ',  am_min**2, am_max**2

c-------- needed in "fast" mode !!! ---------------------------------
        xmd = dsqrt( amx2 )

        t1  = dams2min(   pa1(1), pa(1) )
        t2  = dams2min(   pb1(1), pb(1) )

        q2a = pa1(1)**2 + pa1(2)**2  
        q2b = pb1(1)**2 + pb1(2)**2  

        s1  = sdams2pl(   pa1(1), px(1) )
        s2  = sdams2pl(   pb1(1), px(1) )

        t1q = t1
        t2q = t2

        t1f = dabs(t1)
        t2f = dabs(t2)

        goto 765
c------------------

        if( ihist_fluxes .ne. 0) then

            weight1 =  traj_p1(x1, t1q) / exp( -b0a * q2a) 
            weight2 =  traj_p1(x2, t2q) / exp( -b0b * q2b) 

c -will not work as
c 
c           weight1 =  traj_p(s, s2, t1q) / exp( -b0a * q2a) *(1-x1)*s/s2
c           weight2 =  traj_p(s, s1, t2q) / exp( -b0b * q2b) *(1-x2)*s/s1
c           weight  =  weight1 * weight2
c
c           write(0,*) ' x1, x2  ', x1, x2
c           write(0,*) ' s1, s2 = ', s1, s2, s*(1.-x2), s*(1.-x1)
c           write(0,*)   weight1, weight2

            xtmp1 = 1.-x1
            xtmp2 = 1.-x2

c$$$	    call hfill(61 , xtmp1 ,xdum , weight1)
c$$$	    call hfill(62 , xtmp2 ,xdum , weight2)
c$$$	    call hfill(71 , xtmp1 ,xdum , weight1)
c$$$	    call hfill(72 , xtmp2 ,xdum , weight2)
c$$$
c$$$            call hfill(63 , xtmp1 , t1f , weight1)
c$$$            call hfill(64 , xtmp2 , t2f , weight2)
        endif
c----------------------
 765    continue

CSdv+
C-      write(*,*) ' t1f, t1min, t1max =', t1f, t1min, t1max
C-      write(*,*) ' t2f, t2min, t2max =', t2f, t2min, t2max
C-      write(*,*) ' x1 , x1min, x1max =', x1 , x1min, x1max
C-      write(*,*) ' x2 , x2min, x2max =', x2 , x2min, x2max
C-
        if( t1f .lt. t1min .or. t1f .gt. t1max ) goto 55
        if( t2f .lt. t2min .or. t2f .gt. t2max ) goto 55

        if( x1  .lt. x1min .or. x1  .gt. x1max ) goto 55
        if( x2  .lt. x2min .or. x2  .gt. x2max ) goto 55

        y_fast = 0.5 * dlog( (pa1(4) + pa1(3))/(pa1(4) - pa1(3)) )
        y_slow = 0.5 * dlog( (pb1(4) + pb1(3))/(pb1(4) - pb1(3)) )
        y_cent = 0.5 * dlog( (px(4)  + px(3)) /(px(4)  - px(3)) )
CSdv+
C-      write(*,*) 'y_cent, yminc, ymaxc  =', y_cent, yminc, ymaxc

        if( y_cent .lt. yminc .or. y_cent .gt. ymaxc) goto 55

c--------------------------------------------------------------------

c	 p_pi = dsqrt( dabs(amx2/4. - ampic**2) )
c	 wpipi =  p_pi/ dsqrt(amx2)* dexp( -2.*p_pi**2 ) 

c-------- Pomeron-Pomeron, using (1-x)^al(t) formula 
c
c        weightd =  traj_p1(x1, t1q)*traj_p1(x2, t2q) / dexp( -b0a*q2a - b0b * q2b)
c

c ------- Pomeron-Reggeon, symmetric, using (1-x)^al(t) formula 
c
c        weightd = ( traj_p1(x1, t1q)*traj_r1(x2, t2q)  +  traj_r1(x1, t1q)*traj_p1(x2, t2q) )/ exp( -b0a*q2a - b0b * q2b)
c

c-------- Pomeron-Pomeron, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
c
         weightd =  traj_p(s, s2, t1q)*traj_p(s, s1, t2q) * (1.-x1)*(1.-x2) *s**2/(s1*s2)/ exp( -b0a*q2a - b0b * q2b)

c-------- Pomeron-Reggeon, symmetric, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
c                             weight, using invariant s1/2, s, t1/2 instead of x1/2, t1/2
c
c        weightd = ( traj_p(s,s2,t1q)*traj_r(s,s1,t2q) + traj_r(s,s2,t1q)*traj_p(s,s1,t2q) )* (1.-x1)*(1.-x2)*s**2/(s1*s2)/exp( -b0a*q2a - b0b * q2b)

c        write(0,*) ' weight = ', weight
c        weightd = wpipi

c        Using FUNCB(x1) from  ComGeant instead of 1/(1-|x|)  and   exp(-8bt )
c    

c        weightd =  FUNCBD(x1) *  FUNCBD(x2) * (1.-x1)*(1.-x2) * dexp( 7.1*(t1+t2) ) / dexp( -b0a*q2a - b0b * q2b)
c
c
c  using PP--> pi+pi- cross-section ...
c
c        weightd =  weightd * wpipi

        if( ifast_mode .ne. 0) return

        weightf = weightd

c$$$        call hfill(65 , xtmp1 , t1f  ,weightf)
c$$$        call hfill(66 , xtmp2 , t2f  ,weightf)

        CALL BOOST1(sqs, pcm(1), pa1(1), plab(1,3))
        CALL BOOST1(sqs, pcm(1), pb1(1), plab(1,4))
        CALL BOOST1(sqs, pcm(1), px(1),  plab(1,7))

c-----     introduce phi-dependence if needed ------------------------------------

        do k4 = 1,4
          pbeam1(k4) = PLAB(k4,1) - PLAB(k4,3) 
          pbeam2(k4) = PLAB(k4,2) - PLAB(k4,4) 

            pres(k4) = PLAB(k4,7)              ! in case when decay is absent...
        enddo
        
        xmd = dsqrt( pres(4)**2 - pres(1)**2-  pres(2)**2- pres(3)**2  )
        
        CALL  BOOST2( xmd, pres(1), pbeam1(1),pbeam1_cm(1) )
        CALL  BOOST2( xmd, pres(1), pbeam2(1),pbeam2_cm(1) )

        CALL  BOOST2( xmd, pres(1), PLAB(1,3),pfast_cm(1) )
        CALL  BOOST2( xmd, pres(1), PLAB(1,4),pslow_cm(1) )

        CALL  BOOST2( xmd, pres(1), PLAB(1,1),pbeam_cm(1) )
        CALL  BOOST2( xmd, pres(1), PLAB(1,2),ptarg_cm(1) )

        call vec(pfast_cm(1),pbeam1_cm(1), ygj1(1)  )
        call vec(pslow_cm(1),pbeam2_cm(1), ygj2(1)  )

        ct_12 =  ctd( ygj1(1), ygj2(1) )
c       write(0,*) ' ct_12=', ct_12

        th_12 =   acos(ct_12)

c-------------  usual dphi between pt-s --------------------------------

         yy =  plab(2,3)
         xx =  plab(1,3)
         phipt1 = atan2dmy(yy,xx)
         
         yy =  plab(2,4)
         xx =  plab(1,4)
         phipt2 = atan2dmy(yy,xx)

         dphi= phipt2-phipt1
         if(dphi .lt. -180.) dphi = dphi+360.
         if(dphi .gt.  180.) dphi = dphi-360.
         dphi = abs(dphi)

c-------------  for eta, eta' - mesons... ------------------------------------
cc          weight = weight * sin(dphi*3.141592/180.)**2  * dabs(t1*t2)


c------------- introduce  decay of X and unisotropis angular distributions
c              relative to 2 reggeon collision axis and phi_GJ, in case of M .ne. 0

c        goto 777
c
c        call cp_reweight_simple ( plab(1,1), plab(1,3), plab(1,4), plab(1,7), plab(1,5), weight3 )
c
c        weight = weight * weight3
c
c 777    continue
c
       return
       end
C-----------------------------------------------------------------------------------	
C
      subroutine gener_dp2( xmd,  weightd )
      real *8 xmd, weightd
      real *8 pbeam, ambeam, am1, am2, t1, t2, x1, x2, amtarg, amx, s, sqs, s1, s2 
      real *8 dams2min, dams2pl, sdams2pl, t1q, t2q
      real *8 traj_p1, traj_r1, traj_p, traj_r, FUNCBD
      real rx(10)
      common / ts_xs/  t1, t2, x1, x2
      COMMON/RLRCOM3/PLAB(4,30)
      real *8 PLAB

      real *8 pcm(4), pcm1(4), pfastcm(4), pslowcm(4), pxcm(4), dxfeyn

      common /cm_vectors/  pa(4),pb(4),pa1(4),pb1(4),px(4)
      real *8 pa,pb,pa1,pb1,px

c     real *8	pdecprod(4,10)

      real *8 pbm(4),ax,ay

      common /vertex_simulated/ vertex(3)

c----------- c.m  of cent. produced system  

      real *8 pbeam1(4), pbeam2(4), pres(4), pbeam1_cm(4), pbeam2_cm(4), pslow_cm(4), pfast_cm(4)
      real *8 pbeam_cm(4), ptarg_cm(4)
      real *8 ygj1(4), ygj2(4)
      real *8 ctd

c----- cuts on CP variables, independent of the channel
      common /y_cent_cut/ yminc, ymaxc
      common /x_fast_slow_cuts/ x1min, x1max, x2min, x2max
      common /t_fast_slow_cuts/ t1min, t1max, t2min, t2max

      real *8 p_pi, amx2, wpipi

      real *8 q2a, q2b

      parameter (ampic =0.13956755)

c--------- begin of CP-generator  variables  commons, read from file with ascii cards...

      common /pbeam_f/ pbeamf,  sigma_pbeam
      
      common /itype_beam/  itypb1 

      common /mass_limits/ xma, xmb

      common /number_of_mc_events/ nmc1

      common /mode_reject/ imoderej

c--------- end of CP-generator  variables  commons ---------

      common /fast_geant_code/ icodegeant_fast 

c--------- copy from "cards" values 

       pbeam = pbeamf

c      write(0,*) ' itypb1 = ', itypb1, pbeam

       if( itypb1 .eq. 1) then
           ambeam = 0.93827231
           icodegeant_fast = 14
       endif

       if( itypb1 .eq. 2) then
           ambeam = 0.1395
           icodegeant_fast = 8  ! piplus_beam
       endif

       call mcbeam(pbm,vx,vy,vz,ax,ay)

       vertex(1) = vx
       vertex(2) = vy
       vertex(3) = vz

       do k = 1,3
         plab(k,1) = pbm(k)
       enddo
       plab(4,1) = dsqrt( ambeam**2 + pbm(1)**2+pbm(2)**2+pbm(3)**2  )

       am_min = xma
       am_max = xmb
c------------------------------------------------------------------
       nmes = 0 ! until gener is called and nmes is defined...

       amtarg = 0.93827231

       plab(1,2) = 0.
       plab(2,2) = 0.
       plab(3,2) = 0.
       plab(4,2) =  amtarg

       do i = 1,4
         pcm(i) =  plab(i,1) +  plab(i,2) 
       enddo

       s = pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2
       sqs = dsqrt( s )

       CALL BOOST2(sqs, pcm(1), plab(1,1),  pa(1))
       CALL BOOST2(sqs, pcm(1), plab(1,2),  pb(1))

c--------------------------------------------------------

        zmin = 0. ! defines max of abs(1-x), NO CUT here ...
c
c       zmin =  log(1./( 1.-0.7 ) ) ! defines max of abs(1-x)
c
 55     continue

        call pp2ranlux( RX(1),1)
        z2 = zmin+ 20. * RX(1)

c-----------------------------
        x2 = 1. - exp(-z2)
c---------------------------

        b0b = 8.
        call scatter_pomeron( pb(1), pb1(1), b0b, x2,  amtarg )

        do i = 1,4
          px(i) = pa(i)+pb(i)-pb1(i)
        enddo

        t2 = dams2min(   pb1(1), pb(1)  )

        q2b = pb1(1)**2 +pb1(2)**2  

        s1 = px(4)**2 - px(1)**2-px(2)**2-px(3)**2

        t2q = t2

c       t1q = -q2a
c       t2q = -q2b

c       if( imoderej .eq. 0) then

        weight2 =  traj_p1(x2, t2q) / exp( -b0b * q2b) 

c  -will not work as 
c        weight1 =  traj_p(s, s2, t1q) / exp( -b0a * q2a) *(1-x1)*s/s2
c        weight2 =  traj_p(s, s1, t2q) / exp( -b0b * q2b) *(1-x2)*s/s1

c        weight = weight1 * weight2

         xtmp2 = 1.-x2

c        write(0,*) '  x2  ',  x2, xtmp2
c        write(0,*) '  t2q  ',  t2q
c        write(0,*) '  s1, s2 = ', s1, s2, s*(1.-x2), s*(1.-x1)

c        write(0,*)   ' weight2 = ', weight2

c$$$	 call hfill(62 , xtmp2 ,xdum ,weight2)
c$$$	 call hfill(72 , xtmp2 ,xdum ,weight2)

         t2f = dabs(t2)

c$$$         call hfill(64 , xtmp2 , t2f  ,weight2)

c        endif

        amx2 = px(4)**2-px(1)**2-px(2)**2-px(3)**2

        if( amx2 .le. am_min**2 .or. amx2 .ge. am_max**2 ) goto 55

c--------------------------------------------------------------------

c- actually, not needed ...
        p_pi = dsqrt( dabs(amx2/4. - ampic**2) )
        wpipi =  p_pi/ dsqrt(amx2)* dexp( -2.*p_pi**2 ) 

c-------- Pomeron, using (1-x)^al(t) formula 
c
c        weightd =  traj_p1(x2, t2q) / dexp( - b0b * q2b)
c


c ------- Reggeon,  using (1-x)^al(t) formula 
c        weightd = traj_r1(x2, t2q)/ exp( - b0b * q2b)
c

c-------- Pomeron, using invariant s1/2, s, t1/2 instead of x,1/2, t1/2
c
c        write(0,*) ' s= ', s, ' s1= ',s1, ' x2= ', x2, ' t2q= ',t2q,' q2b= ',q2b
c            ' 
         weightd =  traj_p(s, s1, t2q) * (1.-x2) *s/s1  / exp( - b0b * q2b)

c ------- Reggeon,  using invariant s1/2, s, t1/2 instead of x, 1/2, t1/2
c           weight, using invariant s1/2, s, t1/2 instead of x, 1/2, t1/2
c
c        weightd = traj_r(s, s1,  t2q) *(1.-x2) *s/s1  / exp(- b0b * q2b)

c        write(0,*) ' weight = ', weight
c          weightd = wpipi

c         Using FUNCB(x1) from  ComGeant instead of 1/(1-|x|) and exp(-8bt )
c
c        weightd = FUNCBD(x2) *(1.-x2) * dexp( 7.1 * t2 ) / dexp( - b0b * q2b)
c

c  using PP--> pi+pi- cross-section ...
c
c       weightd =  weightd * wpipi

        weightf = weightd

c$$$        call hfill(66 , xtmp2 , t2f  ,weightf)
c$$$        call hfill(67 , xtmp2 , t2f  ,weightf)

c       CALL BOOST1(sqs, pcm(1), pa1(1),  plab(1,3))
        CALL BOOST1(sqs, pcm(1), pb1(1),  plab(1,4))
        CALL BOOST1(sqs, pcm(1), px(1),   plab(1,7))

c-----     introduce phi-dependence if needed ------------------------------------

        do k4 = 1,4
          pbeam1(k4) = PLAB(k4,1) 
          pbeam2(k4) = PLAB(k4,2) - PLAB(k4,4) 

            pres(k4) = PLAB(k4,7) ! in case when decay is absent...

        enddo
        
        xmd = dsqrt( pres(4)**2 - pres(1)**2 - pres(2)**2 - pres(3)**2  )
        
        CALL BOOST2( xmd, pres(1), pbeam1(1),pbeam1_cm(1))
        CALL BOOST2( xmd, pres(1), pbeam2(1),pbeam2_cm(1))

        CALL BOOST2( xmd, pres(1), PLAB(1,4),pslow_cm(1) )

        CALL BOOST2( xmd, pres(1), PLAB(1,1),pbeam_cm(1) )
        CALL BOOST2( xmd, pres(1), PLAB(1,2),ptarg_cm(1) )

        call vec(pslow_cm(1),pbeam2_cm(1), ygj2(1)  )

       return
       end        
C--------------------------------------------------------------------------------------------
C
 	subroutine mcbeam(pb,vx,vy,vz,ax,ay)   ! Was mcbeam_had, but Sdv change it to mcbeam
c----------------------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c	Subroutine gives beam   pb = ( p(3) , E )
c	 and vertex for event
c---------------- needed for me -------------------
	double  precision  pb(4),eb,ax,ay,ampi
c--------------------------------------------------

      REAL*4 RSIGMA(3), Target_position, Target_thickness
      REAL*8 P_muon_beam, Angle_x, Angle_y

      common /pbeam_f/ pbeamf,  sigma_pbeam

      EXTERNAL PP2RANLUX

c--- not needed for me -----------------------
c      REAL*4 VERT(3)
c      REAL*8 P_beam(3), Azi_Ang(3)

* --- pion beam simulation : 
*
*    I took the data from Table in "Hadron beam" section of Lau Gatignon'
*    paper "The Modifications to the M2 Beam for COMPASS" :
*
*    RMS divergence at target : 0.3 x 0.6 mrad
*    Momentum Resolution  = 1.5 % !!! 
*    Beam spot size RMS at target :   2 x 2 mm
*    These numbers are mentioned as PRELIMINARY ! ( April 1999 )
*
*    Then, I mostly follow the procedure S.Prakhov did in compri.f 
*
*    instead of PP2RNDs :
*      call RNORML(RSIGMA,3)           ! CERNlib V120 gaussian-dist.random number
*      call RANMAR(RSIGMA,3)           ! CERNlib uniform random numbers
*
       ampi = 0.13956755

c      Beam

c     call RNORML(RSIGMA,3)            ! Random numbers are in RSIGMA
      CALL RNORMX(RSIGMA,3,PP2RANLUX)

c     P_muon_beam = 190. * (1. + 0.015 * RSIGMA(1) )          ! 1.5% momt.resolution
      P_muon_beam = pbeamf * (1. + sigma_pbeam * RSIGMA(1) )  ! 1.5% momt.resolution

c     write(0,*) 'sigma_pbeam =', sigma_pbeam, ' P_muon_beam = ' , P_muon_beam 
     
c     Angle_x = 0.0003 * RSIGMA(2)	! RMS = 0.3 mrad
c     Angle_y = 0.0006 * RSIGMA(3)	! RMS = 0.6 mrad

      Angle_x = 0.00
      Angle_y = 0.00

      pb(1) = P_muon_beam * DSIN(Angle_x)
      pb(2) = P_muon_beam * DSIN(Angle_y)
      pb(3) = P_muon_beam * DCOS(Angle_x) * DCOS(Angle_y)

*                                       ! here beam is along 3rd axis yet
*     call PXYZTP( P_beam, Azi_Ang )	! got Theta and Phi
*    			                ! Azi_Ang will be used in DSPROC to rotate an event

      pb(4)=Dsqrt(pb(1)**2+pb(2)**2+pb(3)**2+ampi**2)

      Ax = pb(1)/pb(3)
      Ay = pb(2)/pb(3)

C--------------------
c
c	Vertex

* --- the vertex
*   
*   still, "along the beam" is the THIRD component
*
c     call RNORML(RSIGMA,3)
      CALL RNORMX(RSIGMA,3,PP2RANLUX)

      vx =  0.  +  0.5 * RSIGMA(1)
      vy =  0.  +  0.5 * RSIGMA(2)


c     Target_position  = 0.0              ! "along the beam", in the "hall"
      Target_position  = -25.0            ! "along the beam", in the "hall"

c     Target_thickness = 0.18             ! assuming 0.01 int.length
      Target_thickness = 40.              ! hydrogen target 

c     call RANMAR(RSIGMA,3)
      call PP2RANLUX(RSIGMA,3)
  
      vz = Target_position + Target_thickness*(RSIGMA(1)-0.5)
c ---
      return
      end  
C-------------------    
C
	subroutine pxlab(pb,amtg,amrc,amx,t,px)
c
        double precision pb(4),px(4),px1(4),an(3,3),
     *    amtg,amrc,pcx,ex,ct,st,phi,pi,pbc,anor
        double precision AMX
c
	PARAMETER (ampi = 0.13956755)	
        pi = 3.141592653
c-----------------------------------
c	data ampi /0.13956755/	
c	data pi   /3.141592653/	
c-----------------------------------

	pbc =dsqrt(pb(1)**2+pb(2)**2+pb(3)**2)
	ex = pb(4) + (amtg**2-amrc**2+t)/(2.*amtg)
	pcx = dsqrt( (ex-amx)*(ex+amx) )	
	ct = (t-ampi**2-amx**2+2.*pb(4)*ex)/(2.*pbc*pcx)
	IF (ABS(CT).GT.1.) PRINT*,'Tinv  bad',t
	st = dsqrt(Dabs( (1.-ct)*(1.+ct)) )

	phi = 2.*pi*pp2rnd(-1.)

	px1(4) = ex
	px1(1) = pcx*st*cos(phi)
	px1(2) = pcx*st*sin(phi)
	px1(3) = pcx*ct
c----------rotation  px  to  lab-system-------
	an(1,3) = pb(1)/pbc
	an(2,3) = pb(2)/pbc
	an(3,3) = pb(3)/pbc
c-----------------
	anor = dsqrt( an(2,3)**2+an(1,3)**2 )
	an(1,2) = an(2,3)/anor
	an(2,2) = -an(1,3)/anor
	an(3,2) = 0.
c--------------
	an(1,1) = an(2,2)*an(3,3)
	an(2,1) = -an(1,2)*an(3,3)
	an(3,1) = an(1,2)*an(2,3)-an(2,2)*an(1,3)
c-------------
	if( dsqrt( pb(1)**2+pb(2)**2 )/abs(pb(3)) .le.1.e-10) then
c	call ucopy(px1(1),px(1),8)
c	  do i = 1,8 ! ERROR
	  do i = 1,4
	  px(i) = px1(i)
	  enddo
	return
	endif
c
	px(1) = an(1,1)*px1(1) + an(1,2)*px1(2) + an(1,3)*px1(3) 
	px(2) = an(2,1)*px1(1) + an(2,2)*px1(2) + an(2,3)*px1(3) 
	px(3) = an(3,1)*px1(1) + an(3,2)*px1(2) + an(3,3)*px1(3) 
	px(4) = px1(4)
	return
	end
c-------------------------------------	
	subroutine rexp(t,amx,Wmtarg,Wmrec,q2beam,pb)
c
	real *8  Wmtarg,Wmrec,pb(4)
        real *8  AMX
	common /t_inv/ tt
	common /parexp/ b1,b2,p
	COMMON /TLIM/ TMIN,TMAX
        COMMON /EXS/EMIN1,ER1,EMIN2,ER2,PN
c-------------------------------------
	parameter (AMPI = 0.13956755)
c-------------------------------------

        real *8 s, shs, ex, px, eb, pbc
        real *8 q2beam
        real xx(2)

c       write(0,*) 'rexp: amx = ', amx, ' TMIN,TMAX = ', TMIN,TMAX
	S = q2beam+Wmtarg**2+2.*PB(4)*Wmtarg
	shs = dsqrt(s)
	ex = (s+amx**2-Wmrec**2)/(2.*shs)
	px = dsqrt(ex**2-amx**2)
	eb = (s+q2beam-Wmtarg**2)/(2.*shs)
	pbc = dsqrt(eb**2-q2beam)
C-----------------------
	tminPH = q2beam+amx**2-2.*eb*ex-2.*pbc*px
	tmaxPH = q2beam+amx**2-2.*eb*ex+2.*pbc*px
C----------------------------
c	TMINF = TMAXPH - TMAX
c	TMAXF = TMAXPH - TMIN
	TMINF =        - TMAX
	TMAXF =        - TMIN
c
        A=P*(EXP(B1*TMAXF)-EXP(B1*TMINF))
        B=(1.-P)*(EXP(B2*TMAXF)-EXP(B2*TMINF))
        PN=A/(A+B)

c        write(0,*) ' tmin, tmax = ', tmin, tmax
c        WRITE(0,*)  'b1,b2,p=',b1,b2,p

        EMIN1=EXP(B1*TMINF)
        EMAX1=EXP(B1*TMAXF)
        EMIN2=EXP(B2*TMINF)
        EMAX2=EXP(B2*TMAXF)
        ER1=EMAX1-EMIN1
        ER2=EMAX2-EMIN2
c
c       SP= pp2rnd(-1)
        CALL PP2RANLUX(xx(1),1)
        SP = xx(1)

	IF( SP.LE.PN) THEN
            E=EMIN1+ER1*SP/PN
            T=(1./B1)*ALOG(E)
        ENDIF
	      
        IF( SP.GT.PN .AND. PN.LT.1.) THEN
             E=EMIN2+ER2*( SP - PN )/(1. - PN )
             T=(1./B2)*ALOG(E)
        ENDIF

	t = t + tmaxph

c	WRITE(0,*) 'rexp:T',t

	tt = t
        RETURN
        END
c----------------------------------------------

        subroutine tdist1(t,amx,Wmtarg,Wmrec, q2beam, pb)
c
	real *8  Wmtarg,Wmrec,pb(4)
        real *8  AMX
	common /t_inv/ tt
	common /parexp/ b1,b2,p
	COMMON /TLIM/ TMIN,TMAX
        COMMON /EXS/EMIN1,ER1,EMIN2,ER2,PN
c-------------------------------------
	parameter (AMPI = 0.13956755)
c-------------------------------------
        common / tprime_generated/ t1
	common / init_tdist1/ init_t, vmax1

        real *8 s, shs, ex, px, eb, pbc
        real *8 q2beam
        real xx(2)
c         write(0,*) 'rexp: amx = ', amx, ' TMIN,TMAX = ', TMIN,TMAX
c         write(0,*) 'tgen: wmtarg, wmrec = ', wmtarg, wmrec
c         write(0,*) 'tgen: amx = ', amx
c         write(0,*) 'tgen: pb(4) = ', pb(4)

	S = q2beam+Wmtarg**2+2.*PB(4)*Wmtarg
	shs = dsqrt(s)
	ex = (s+amx**2-Wmrec**2)/(2.*shs)
	px = dsqrt(ex**2-amx**2)
	eb = (s+q2beam-Wmtarg**2)/(2.*shs)
	pbc = dsqrt(eb**2-q2beam)
C-----------------------
	tminPH = q2beam+amx**2-2.*eb*ex-2.*pbc*px
	tmaxPH = q2beam+amx**2-2.*eb*ex+2.*pbc*px
C----------------------------

c--------- generating itself --------------------
	if( init_t.eq. 0) then
	  init_t = 1
	  vmax1 = 0.
	  write(0,*) ' before initialising user-defined t-distrib'
	  do k = 1, 100000

               CALL PP2RANLUX(xx(1),1)

	    t1 = tmin + (tmax-tmin)*xx(1)
	    fun = tfunc1(t1)
	    if( fun .gt. vmax1 ) then
	      vmax1 = fun
	      t1_vmax = t1
	      endif
	  enddo
	  write(0,*) ' after init, vmax = ', vmax1,' t_vmax= ',t1_vmax 
	    vmax1 = 1.5*vmax1
          endif
   
 33       CALL PP2RANLUX(xx(1),1)       
          t1 = tmin + (tmax-tmin)*xx(1)       
          fun = tfunc1(t1)

          CALL PP2RANLUX(xx(1),1)
	  y  = vmax1 * xx(1)
	  if(y .gt.  fun) goto 33
c---------------------------------
	t = -t1 + tmaxph
c	write(0,*) 'tgen: t1= ', t1,  '  tmaxph = ', tmaxph
	tt = t

        RETURN
        END
c
C========= User-defined   beam energy profile ================
c 
       real function efunc1(eb1)
       common  /beam_en_type/ itype_ebeam

	if( itype_ebeam .eq. 1) then
	 efunc1 = 0.3*exp( -0.5*( (eb1-188.)/1.2 )**2 ) + 0.5 * exp( -0.5*( (eb1-190.)/2.0 )**2 ) 
	endif

	return
	end

C========= User-defined   tprime distribution  ================

	real function tfunc1(t1)
	common  /tdist_type/ itype_tdist

	if( itype_tdist .eq. 1) then
	    tfunc1 = exp(- 3. * t1)
	endif

	return
	end
c
C-------------------
      real*8 function traj_Psdv(s,s1,t)
      real*8 s,s1,t, z, alf

      z = dlog(s/s1)
      alf = 1.08 + 0.25 *t !   M.Albrow et al. MAN/HEP/2010/1 June 8, 2010

      traj_Psdv = dexp( 5.*t ) * dexp( 2.*z *( alf - 1.) )


      return
      end
c
C-------------------
      real*8  function traj_p(s,s1,t)
      real*8 s,s1,t, z, alf

      z = dlog(s/s1)
c     alf = 1.0 + 0.40 *t
C      alf = 1.2 + 0.25 *t
      alf = 1.08+ 0.25 *t       ! M.Albrow et al. MAN/HEP/2010/1 June 08,2010

      traj_p = dexp( 0.15*t ) * dexp( 2.*z *( alf - 1.) )
c     traj_p =  exp( 8.*t )

      return
      end
C-------------------
      real*8 function traj_r(s,s1,t)
      real*8 s,s1,t, alf, z

      z = dlog(s/s1)

      alf = 0.5 + 1.1*t ! Cox-Forshaw reggeon ?
      traj_r = dexp( 5.*t ) * dexp( 2.*z *( alf - 1.) )
c     traj_r =  exp( 8.*t )

      return
      end
	
C--------------------------------------------------------------------------------
C    Library routines from compassPWAssh/pwa/utils/dlib.f
C
C-------------------
        SUBROUTINE BOOST1(AM,P3,pi,pf)
        DOUBLE PRECISION P3(4),PI(4),PF(4),A,B,AM
        A=0.
        DO I=1,3
          A=A+P3(I)*PI(I)
        END DO
        PF(4)=(P3(4)*PI(4)+A)/AM
        B=(PF(4)+PI(4))/(P3(4)+AM)
        DO I=1,3
          PF(I)=PI(I)+  P3(I)*B
        ENDDO
        RETURN
        END
C---------------------
        SUBROUTINE BOOST2(AM,P3,pi,pf)
        DOUBLE PRECISION P3(4),PI(4),PF(4),A,B,AM
        A=0.
        DO I=1,3
          A=A-P3(I)*PI(I)
        END DO
        PF(4)=(P3(4)*PI(4)+A)/AM
        B=(PF(4)+PI(4))/(P3(4)+AM)
        DO I=1,3
          PF(I)=PI(I) -  P3(I)*B
        ENDDO
        RETURN
        END
C---------------------
	subroutine vec(p1,p2,p3)
	double precision p1(4),p2(4),p3(4)
	p3(1)=p1(2)*p2(3)-p1(3)*p2(2)
	p3(2)=p1(3)*p2(1)-p1(1)*p2(3)
	p3(3)=p1(1)*p2(2)-p1(2)*p2(1)
	return
	end
C------------------------
	subroutine scl(p1,p2,res)
	double precision p1(4),p2(4)
	res=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
	return
	end
C------------------------	
       real function atan2dmy(x,y)
        real x,y
c        atan2d  = 0.
c        if( y.ne.0. ) atan2dmy = x/y
 
            atan2dmy = atan2(x,y)*180./3.141592654
      
         return
        end
C------------------- 
      real *8 function dams2min(p1, p2)
      real *8 p1(4),p2(4)

      dams2min = (p1(4)-p2(4))**2 - (p1(1)-p2(1))**2 - (p1(2)-p2(2))**2 - (p1(3)-p2(3))**2 

      return
      end        
C-------------------       
      subroutine scatter_pomeron( pa, pa1, b, x,  am )
      real *8 pa(4), pa1(4)
      real *8  x, am, PI2
      real X1(2)
      DATA PI2 /6.2831853071796/

      CALL PP2RANLUX(X1(1),2)

      xx = X1(1)             !
      q2 =-alog(xx)/b        !  Sdv:  
      q  = sqrt(q2)          !  

      phi   = X1(2)*PI2
      pa1(1)= q*cos(phi)
      pa1(2)= q*sin(phi)
      pa1(3)= x*pa(3)
      pa1(4)= dsqrt( pa1(1)**2 + pa1(2)**2 + pa1(3)**2 + am**2 )

      return
      end
C-------------------       
      real *8 function sdams2pl(p1, p2)
      real *8 p1(4),p2(4)

      sdams2pl = (p1(4)+p2(4))**2 - (p1(1)+p2(1))**2 - (p1(2)+p2(2))**2 - (p1(3)+p2(3))**2

      return
      end
C-------------------  
      real*8 function ctd(a, b)
      real*8 a(4),b(4)
      ctd = (a(1)*b(1)+a(2)*b(2)+a(3)*b(3)) / dsqrt((a(1)**2+a(2)**2+a(3)**2)*(b(1)**2+b(2)**2+b(3)**2))
      return
      end
C------------------- 
      real*8 function traj_p1(x,t)
      real*8 x, t, alf, z
      
      z = dlog(1.-x)
c
c     alf = 1.0  + 0.40 *t
c     alf = 1.2  + 0.25 *t
c     alf = 1.08 + 0.25 *t  ! Bialas-Landshoff pomeron
      alf = 1.20 + 0.25 *t  ! Cox-Forshaw pomeron 
 
      traj_p1 = dexp( 5.*t ) * (1-x)**(-2.*( alf - 1.) )
c     traj_p1 = dexp( 5.*t ) * exp( -2.*z* ( alf - 1.) )  ! regularization here ...haven't checked...
c     traj_p1 =  exp( 8.*t )

      return
      end
C------------------- 
      subroutine generate_t(t,amx,Wmtarg,Wmrec,q2beam, pb)

C Sdv ----  !!!!!!  include 'z_usercards.inc'  !!!!!!  ----

      common  /tdist_type/ itype_tdist

      real *8  Wmtarg,Wmrec,pb(4)
      real *8  AMX
      real *8 q2beam

      COMMON /TLIM/ TMIN,TMAX

      common / idebug_mcgen1/ idebug_mcgen1
c     common /tprime_limits/ tprime_min, tprime_max

c     write(0,*) '  tmin ' , tmin, tmax, ' tprime_min,  tprime_max ', tprime_min , tprime_max 
c     write(0,*) '  itype_tdist=   ', itype_tdist

      TMIN = tprime_min_mcgen 
      TMAX = tprime_max_mcgen

      if( idebug_mcgen1 .ne. 0) then
          write(0,*) ' tprime min max mcgen =  ', tprime_min_mcgen, tprime_max_mcgen
      endif

c     write(0,*) ' tmin, tmax = ', tmin, tmax

      if( itype_tdist.eq. 0) then
	        call rexp(t,amx,Wmtarg,Wmrec,q2beam, pb(1))
          endif
	  
      if( itype_tdist.ne. 0) then
	        call tdist1(t,amx,Wmtarg,Wmrec, q2beam, pb(1))
	  endif
	  
      if( idebug_mcgen1 .ne. 0) then
               write(0,*) ' tinv in generate_t =  ', t
          endif

	return
	end

C-----------------
C
c$$$      SUBROUTINE BOOKIN      
c$$$      parameter         (npawc=250000)
c$$$      COMMON/PAWC/  HMEM(npawc)
c$$$C 
c$$$      call hbook1( 10,'Mass of the X system', 500, 0.0, 25., 0.)
c$$$      call hbook1( 11,'x1                  ', 500,-0.0, 1.0, 0.)
c$$$      call hbook1( 12,'x2                  ', 500,-0.0, 1.0, 0.)
c$$$      call hbook1( 19,'Weighd              ', 300,-0.0, 3.E3,0.)
c$$$      call hbook1( 20,'Mass weighted       ', 500, 0.0, 25., 0.)
c$$$      call hbook1( 21,'x1 weighted         ', 500,-0.0, 1.0, 0.)
c$$$      call hbook1( 22,'x2 weighted	   ', 500,-0.0, 1.0, 0.)
c$$$      call hbook1( 27,'Weighd              ', 300,-0.0, 3.E3,0.)
c$$$      call hbook1( 28,'Mass Weighd         ', 500, 0.0, 10., 0.)     
c$$$      call hbook1( 29,'Weighd * Wpipi      ', 300,-0.0, 3.E3,0.)
c$$$      call hbook1( 30,'Mass Weighd * Wpipi ', 500, 0.0, 10., 0.)
c$$$C      
c$$$      call hbook1( 100,'Mass of the X system',500, 0., 5.0,  0.)
c$$$      call hbook1( 101,'Px of particles   ',  600,-3., 3.0,  0.)
c$$$      call hbook1( 102,'Py of particles   ',  600,-3., 3.0,  0.)
c$$$      call hbook1( 103,'Pz of particles   ',  600,-3., 3.0,  0.)
c$$$      call hbook1( 104,'E  of particles   ',  500, 0., 5.0,  0.)
c$$$C      
c$$$      call hbook1( 111,'cos of polar angle in Lab',200,-1.,1.   ,0.)
c$$$      call hbook1( 112,'phi of athim angle in Lab',200, 0.,6.283,0.)
c$$$      call hbook1( 113,'eta of particles      ', 200.,-2., 2.,   0.)
c$$$      call hbook1( 114,'Y   of particles      ', 200.,-2., 2.,   0.)
c$$$      call hbook2( 115,'phi vs eta of particle', 100.,-2., 2., 
c$$$     +                                           100,  0.,6.283, 0.)
c$$$      call hbook1( 116,'Delta phi 12 in Lab   ',200,-6.283,6.283,0.)
c$$$C-    call hbook2( 117,'phi2 vs phi1          ', 100., 0.,6.283, 
c$$$C-   +                                           100,  0.,6.283, 0.)
c$$$C
c$$$      call hbook1( 200,'Mass of the X system ',  500, 0., 5.0,  0.)
c$$$      call hbook1( 201,'Px of the X system   ',  600,-3., 3.0,  0.)
c$$$      call hbook1( 202,'Py of the X system   ',  600,-3., 3.0,  0.)
c$$$      call hbook1( 203,'Pz of the X system   ',  600,-3., 3.0,  0.)
c$$$      call hbook1( 204,'E  of the X system   ',  500, 0., 5.0,  0.)
c$$$      call hbook1( 205,'Pt of the X system   ',  300, 0., 3.0,  0.)  
c$$$      call hbook1( 211,'cos of X system in Lab ',200,-1., 1.   ,0.)
c$$$      call hbook1( 212,'phi of X system in Lab ',200, 0., 6.283,0.)
c$$$      call hbook1( 213,'eta of X system in Lab ',200,-2., 2.0  ,0.)
c$$$      call hbook1( 214,'Y   of X system in Lab ',200,-2., 2.0  ,0.)    
c$$$      call hbook1( 221,'CosGJ of 1st part.in X ',200,-1., 1.0  ,0.)	
c$$$      call hbook1( 222,'PhiTY of 1st part.in X ',200, 0.,6.2832,0.)
c$$$      call hbook2( 223,'pi+/pi-: PhiTY vs CosGJ',100,-1.,1.0,
c$$$     +                                           100, 0.,6.283, 0.) 
c$$$C
c$$$C-    call hbook1(1200,'Mass of like charge pairs ', 350, 0., 3.5, 0.)
c$$$C-    call hbook1(1202,'Py of Like charge pairs   ', 300, 0., 3.0, 0.)
c$$$C-    call hbook1(1201,'Px of Like charge pairs   ', 300, 0., 3.0, 0.)
c$$$C-    call hbook1(1202,'Py of Like charge pairs   ', 300, 0., 3.0, 0.)
c$$$C-    call hbook1(1203,'Pz of Like charge pairs   ', 300, 0., 3.0, 0.)
c$$$C-    call hbook1(1204,'E  of Like charge pairs   ', 300, 0., 3.0, 0.)
c$$$C-    call hbook1(1205,'Pt of Like charge pairs   ', 300, 0., 3.0, 0.)
c$$$C-    call hbook1(1221,'CosGJ, 1st part. Like charg',200,-1., 1.0, 0.)      
c$$$C-    call hbook1(1222,'PhiTY, 1st part. Like charg',200, 0.,6.283,0.)         
c$$$C
c$$$C-    call hbook1(2200,'Mass of Unlike charge pairs',350, 0., 3.5, 0.)
c$$$C-    call hbook1(2201,'Px of UnLike charge pairs  ',600,-3., 3.0, 0.)
c$$$C-    call hbook1(2202,'Py of UnLike charge pairs  ',600,-3., 3.0, 0.)
c$$$C-    call hbook1(2203,'Pz of UnLike charge pairs  ',600,-3., 3.0, 0.)
c$$$C-    call hbook1(2204,'E  of UnLike charge pairs  ',600, 0., 3.0, 0.)
c$$$C-    call hbook1(2205,'Pt of UnLike charge pairs  ',300, 0., 3.0, 0.) 
c$$$C-    call hbook2(2207,'Unlike:  Pt vs pair Mass   ',250, 0.,2.5,
c$$$C-   +                                               100, 0.,2.0,  0.)
c$$$C
c$$$C-    call hbook1(2211,'cos of UnLike charge pairs ',200,-1., 1.   ,0.)
c$$$C-    call hbook1(2212,'phi of UnLike charge pairs ',200, 0., 6.283,0.)
c$$$C-    call hbook1(2213,'eta of UnLike charge pairs ',200,-2., 2.0  ,0.)
c$$$C-    call hbook1(2221,'CosGJ, 1st part. unLike charg',200,-1., 1.0,0.)      
c$$$C-    call hbook1(2214,'Y   of UnLike charge pairs ',200,-2., 2.0  ,0.)  
c$$$C-    call hbook1(2221,'CosGJ, 1st part. unLike charg',200,-1., 1.0,0.)      
c$$$C-    call hbook1(2222,'PhiTY, 1st part. unLike charg',200,0.,6.283,0.) 
c$$$C-    call hbook2(2223,'pi+/pi-: PhiTY vs CosGJ    ',100,-1.,1.0,
c$$$C-   +                                               100, 0.,6.283,0.) 
c$$$C-    call hbook2(2224,'pi+/pi-: CosGJ vs pair Mass',250, 0.,2.5,
c$$$C-   +                                               100,-1.,1.0,  0.) 
c$$$C-    call hbook2(2225,'pi+/pi-: PhiTY vs pair Mass',250, 0.,2.5,
c$$$C-   +                                               100, 0.,6.283,0.)      
c$$$C
c$$$C-    call hbook1(4200,'Mass of pi+pi- system',350, 0., 3.5,  0.)
c$$$C-    call hbook1(4201,'Px  of  pi+pi- system',600,-3., 3.0,  0.)
c$$$C-    call hbook1(4202,'Py  of  pi+pi- system',600,-3., 3.0,  0.)
c$$$C-    call hbook1(4203,'Pz  of  pi+pi- system',600,-3., 3.0,  0.)
c$$$C-    call hbook1(4204,'E   of  pi+pi- system',300, 0., 3.0,  0.)
c$$$C-    call hbook1(4205,'Pt  of  pi+pi- system',300, 0., 3.0,  0.) 
c$$$C-    call hbook1(4206,'Delta phi pi+/pi- in Lab',200,-6.283,6.283,0.)
c$$$C-    call hbook2(4207,'pi+pi-:  Pt vs pair Mass   ',250, 0.,2.5,
c$$$C-   + 					      100, 0.,2.0,  0.)      
c$$$C-    call hbook1(4211,'cos of  pi+pi- system',200,-1., 1.0,  0.)
c$$$C-    call hbook1(4212,'phi of  pi+pi- system',200, 0., 6.283,0.)
c$$$C-    call hbook1(4213,'eta of  pi+pi- system',200,-2., 2.0  ,0.)
c$$$C-    call hbook1(4214,'Y   of  pi+pi- system',200,-2., 2.0  ,0.)     
c$$$C     
c$$$C-    call hbook1(4221,'pi+/pi-: CosGJ of pi+',200,-1., 1.0,  0.)     
c$$$C-    call hbook1(4222,'pi+/pi-: PhiTY of pi+',200, 0.,6.283, 0.)	
c$$$C
c$$$C-    call hbook2(4223,'pi+/pi-: PhiTY vs CosGJ    ',100,-1.,1.0,
c$$$C-   + 					      100, 0.,6.283,0.) 
c$$$C-    call hbook2(4224,'pi+/pi-: CosGJ vs pair Mass',250, 0.,2.5,
c$$$C-   + 					      100,-1.,1.0,  0.) 
c$$$C-    call hbook2(4225,'pi+/pi-: PhiTY vs pair Mass',250, 0.,2.5,
c$$$C-   + 					      100, 0.,6.283,0.)
c$$$C                  
c$$$      RETURN
c$$$      END
C
      SUBROUTINE KinPar(P,cosT,phi,eta,Y)      
      implicit double precision (A-H,O-Z)
      real cosT, phi, eta, Y
      dimension P(6)    
      DATA PI2  / 6.283185307179586476 D0 /
C 
      cosTT =      P(3)/dsqrt(P(1)**2+P(2)**2+P(3)**2) 
      phi   =dacos(P(1)/dsqrt(P(1)**2+P(2)**2))
      if (P(2).lt.0.) phi = pi2-phi  
C 
      cosT= cosTT
      eta = 0.5*dlog((1.+cosTT)/(1.-cosTT))
      Y   = 0.5*dlog((P(4)+P(3))/(P(4)-P(3))) 
      RETURN
      END          
C
C
      SUBROUTINE PP2PXPi(P1,PX,Pi,COSGJ,PHITY)
C                                                      S.Sadovsky 15.05.2012
C                                                      corrected by S.Evdikimov
C     Calculation of the COS_GJ and PHI_TY of particle Pi in the PX CM system:
C     P1 - 5-momentum (3-mom + E + mass) of the first proton in LAB
C     PX - 5-momentum of the central X-system in LAB
C     Pi - 5-momentum of the secondary particle in LAB
C     ------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)          ! Sdv
      real*8 COSGJ, PHITY, cosfi
      DIMENSION P1(6),PX(6),Pi(6)
      DIMENSION XX(6),YY(6),ZZ(6),Pc(6),Pd(6),Yc(6)
      DATA PI2 / 6.283185307179586476 D0 /
C      
      CALL VAB  (P1,PX,YY)
      CALL ARTUR(P1,PX,ZZ) 
C      
      CALL VAB  (YY,ZZ,XX)
      CALL ARTUR(Pi,PX,Pc) 
C           
      Pcm=DSQRT(SAB(Pc,Pc))
      COSGJ=SAB(Pc,ZZ)/(Pcm*DSQRT(SAB(ZZ,ZZ)))
C
      PcX  = SAB(Pc,XX)/DSQRT(SAB(XX,XX))
      PcY  = SAB(Pc,YY)/DSQRT(SAB(YY,YY))
      COSFI= PcX/DSQRT(PcX*PcX+PcY*PcY)
      PHITY= DACOS(COSFI)
      if (PcY.lt.0.) PHITY=PI2-PHITY
C      
      RETURN
      END
C
C
      SUBROUTINE ARTUR(PL,P0,PC)
C
C     Lorenc transformation PL from LAB to the CM system:     
C     PL(5) - 5-momentum (3-mom + E + mass) of particle in LAB
C     P0(5) - 5-momentum (3-mom + E + mass) of decay particle in LAB
C     PC(5) - 5-momentum (3-mom + E + mass) of particle in CM of decay particle
C     --------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)    ! ART
      DIMENSION PC(5),P0(5),PL(5)
C
      AM0 =  P0(5)
      EL4 = (PL(4)*P0(4) - PL(3)*P0(3) - PL(2)*P0(2) - PL(1)*P0(1))/AM0
      DO 10 I = 1,3
   10 PC(I) = PL(I) - P0(I)*(EL4+PL(4))/(P0(4)+AM0)
      PC(4) = EL4
c      PC(5) = PL(5)
      RETURN
      END
C
C
      SUBROUTINE ARTURS(PC,P0,PL)
C
C     Lorenc transformation PC from CM to the LAB system:   
C     PC(5) - 5-momentum (3-mom + E + mass) of particle in CM of decay particle
C     P0(5) - 5-momentum (3-mom + E + mass) of decay particle in LAB
C     PL(5) - 5-momentum (3-mom + E + mass) of the particle in LAB
C     --------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PC(5),P0(5),PL(5)
C
      AM0 =  P0(5)
      EL4 = (PC(4)*P0(4) + PC(3)*P0(3) + PC(2)*P0(2) + PC(1)*P0(1))/AM0
      DO 10 I = 1,3
   10 PL(I) = PC(I) + P0(I)*(EL4+PC(4))/(P0(4)+AM0)
      PL(4) = EL4
      PL(5) = PC(5)
      RETURN
      END
C
C
      FUNCTION SAB(A,B)
C     -------------------
      IMPLICIT REAL*8 (A-H,O-Z)    ! ART
      DIMENSION A(3),B(3)
C
      SAB = A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      RETURN
      END
C
C
      SUBROUTINE VAB(A,B,P)
C     ---------------------
      IMPLICIT REAL*8 (A-H,O-Z)    ! ART
      DIMENSION A(3),B(3),P(3)
C
      P(1) = A(2)*B(3)-A(3)*B(2)
      P(2) = A(3)*B(1)-A(1)*B(3)
      P(3) = A(1)*B(2)-A(2)*B(1)
      RETURN
      END
C
C
      FUNCTION DOT4(P)
C     ----------------
      IMPLICIT REAL*8 (A-H,O-Z)    ! ART
      DIMENSION P(4)
C
      DOT4 = P(4)**2 - P(1)**2 - P(2)**2 - P(3)**2
      RETURN
      END
C
C**************************************************************************
C
      SUBROUTINE STAR3T(PL0,PL1,PL2,PL3)
C.........................................................................
C   TO GENERATE IN CMS ISOTROPIC DECAY OF ONE PARTICLE INTO THREE ONE
C   INPUT  VALUES:   PL0(5),PL1(5),PL2(5),PL3(5) - MASSES OF PARTICLES
C   OUTPUT VALUES:   PL1(1-4),PL2(1-4),PL3(1-4) - 4-VECTORS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL PP2RND
      DIMENSION PL0(5),PL1(5),PL2(5),PL3(5),P1(5),P2(5),P3(5)
      DATA PI2 / 6.283185307179586476 D0 /
C
      T0 = PL0(5) - PL1(5) - PL2(5) - PL3(5)
C
   10 R1 = PP2RND(-1.)
      R2 = PP2RND(-1.)
      IF  (R1.LT.R2) GO TO 20
      R1 = 1. - R1
      R2 = 1. - R2
C
   20 T1 = T0* R1
      T2 = T0*(R2-R1)
      T3 = T0- T1-T2
c
      EE1 = PL1(5) + T1
      EE2 = PL2(5) + T2
      PP1 = DSQRT( T1*(T1 + 2.*PL1(5)) )
      PP2 = DSQRT( T2*(T2 + 2.*PL2(5)) )
      PP3 = DSQRT( T3*(T3 + 2.*PL3(5)) )
C
      IF(PP1.GT.PP2+PP3) GO TO 10
      IF(PP2.GT.PP3+PP1) GO TO 10
      IF(PP3.GT.PP1+PP2) GO TO 10
C
      C1  = 2.*PP2RND(-1.) - 1.
      S1  = DSQRT( 1.-C1*C1 )
      F1  = PI2*PP2RND(-1.)
      CF1 = DCOS(F1)
      SF1 = DSIN(F1)
      P1(1) = PP1*S1*CF1
      P1(2) = PP1*S1*SF1
      P1(3) = PP1*C1
      P1(4) = EE1
C
      C12 =(PP3*PP3 - PP2*PP2 - PP1*PP1)/(2.*PP1*PP2)
      S12 = DSQRT(1.- C12*C12 )
      F12 = PI2*PP2RND(-1.)
      CF12= DCOS(F12)
      SF12= DSIN(F12)
      X2  = PP2*S12*CF12
      Y2  = PP2*S12*SF12
      Z2  = PP2*C12
C
      P2(1) = C1*CF1*X2 - SF1*Y2 + S1*CF1*Z2
      P2(2) = C1*SF1*X2 + CF1*Y2 + S1*SF1*Z2
      P2(3) =-S1*X2              + C1*Z2
      P2(4) = EE2
C
      CALL ARTURS(P1(1),PL0(1),PL1(1))
      CALL ARTURS(P2(1),PL0(1),PL2(1))
C
      DO 30 I=1,4
   30 PL3(I) =PL0(I) - PL1(I) - PL2(I)
C
      RETURN
      END
C
      SUBROUTINE DECAYS(P0,P1,P2)
C...........................................................................
C   TO GENERATE ISOTROPIC (IN C.M.S.) DECAY OF ONE PARTICLE INTO TWO ONES
C   INPUT  VALUES:   P0(1-4) AND P0(5),P1(5),P2(5) - MASSES OF PARTICLES
C   OUTPUT VALUES:   P1(1-4),P2(1-4)
C   THE SUBROUTINE CALLS SUBR. ARTURS TO DO LORENZ TRANSFORMATION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P0(5),P1(5),P2(5),PP1(5)
      REAL PP2RND,X(2)
      DATA PI2 / 6.283185307179586476 D0 /
C
      CALL PP2RANLUX(X(1),2)  
      E1=(P0(5)**2 + P1(5)**2 - P2(5)**2 )/(2*P0(5))
      E2= P0(5) - E1
      P = DSQRT(  E1**2 - P1(5)**2 )
C-    CT= 2.*PP2RND(-1.) -1.
      CT= 2.*X(1) -1.
      ST= DSQRT(1.-CT**2)
C-    FI= PP2RND(-1.)*PI2
      FI= X(2)*PI2
      PP1(1)=P*ST*DCOS(FI)
      PP1(2)=P*ST*DSIN(FI)
      PP1(3)=P*CT
      PP1(4)=E1
      PP1(5)=P1(5)
      CALL ARTURS(PP1(1),P0(1),P1(1))
      DO 1 I=1,4
    1 P2(I) = P0(I) - P1(I)
      RETURN
      END
C
C--------------- S.Evdokimov routines -------------------------
      subroutine star4m(P0,P1,P2,P3,P4)
C      
C     first P0->P1+P234 then P234->P2+P3+P4
C     input: P*(5) - masses of particles, P0(1-4) - 4-momentum
C     output P*(1-5) - 5-momentum
C     IMPORTANT: P2(5)=P3(5)=P4(5)  DO NOT FORGET!!!
C     3-particle phase space approx.calculated only for this case
C     G.I.Kopylov, "Osnovy kinematiki resonansov", str.215, formula 28
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL PP2RND, ENERG
      DIMENSION P0(5),P1(5),P2(5),P3(5),P4(5),P234(5)
      DATA PI / 3.1415926 /
      NCYCL=0
 10   R1=PP2RND(-1.)
      R2=PP2RND(-1.)
C     SA lectures 
      AMU       =(P2(5)+P3(5)+P4(5))**2
      SQRM234   =AMU+R1*((P0(5)-P1(5))**2-AMU)
      ALPHA     =P2(5)/DSQRT(SQRM234)
      ALPHAMIN  =1./3.
C
      RHO       =DSQRT(((P0(5)**2+P1(5)**2-SQRM234)**2)/
     /     (2*P0(5))**2-P1(5)**2)*FUNC(ALPHA)*SQRM234
C
      AMAX      =DSQRT((((P0(5)**2+P1(5)**2-AMU)**2)/
     /     (2*P0(5)))**2-P1(5)**2)*AMU
C
      NCYCL=NCYCL+1
      IF(NCYCL.GT.10000) GOTO 20
      if((AMAX*R2).gt.RHO) goto 10
C
 20   P234(5)=DSQRT(SQRM234)
      call DECAYS(P0,P1,P234)
      call STAR3T(P234,P2,P3,P4)
      return 
      end
C
      function FUNC(ALPHA)
C     G.I.Kopylov, "Osnovy kinematiki resonansov", str.215, formula 28
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FUNC=(1+0.997*ALPHA**2)*((1-ALPHA)**1.5)*
     *     (1-9*ALPHA**2)**2
      return
      end
C
      subroutine star53(P0,P1,P2,P3,P4,P5)
C      
C     first P0->P1+P2345 then P2345->P2+P3+P4+P5
C     input: P*(5) - masses of particles, P0(1-4) - 4-momentum
C     output P*(1-5) - 5-momentum
C     IMPORTANT: P3(5)=P4(5)=P5(5)  DO NOT FORGET!!!
C     cause i'm using star4t only for this case
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL PP2RND, ENERG
      DIMENSION P0(5),P1(5),P2(5),P3(5),P4(5),P5(5),P2345(5),P345(5)
      DATA PI / 3.1415926 /
 10   R1=PP2RND(-1.)
      R2=PP2RND(-1.)
      R3=PP2RND(-1.)
C      
C     SA lectures 
      AMU      = (P2(5)+P3(5)+P4(5)+P5(5))**2
      SQRM2345 = AMU+R1*((P0(5)-P1(5))**2-AMU)
      AMU345   = (P3(5)+P4(5)+P5(5))**2
      SQRM345  = AMU345+R2*((DSQRT(SQRM2345)-P2(5))**2-AMU345)
      ALPHA    = P3(5)/DSQRT(SQRM345)
C
      ALERT=(((P0(5)**2+P1(5)**2-SQRM2345))/
     /     (2*P0(5)))**2-P1(5)**2
c     write(*,*) ALERT, SQRM2345
      RHO       =DSQRT((((P0(5)**2+P1(5)**2-SQRM2345))/
     /     (2*P0(5)))**2-P1(5)**2)*
     *     FUNC(ALPHA)*SQRM345*    ! 3-particle phase space
     *     DSQRT((((SQRM2345+P2(5)**2-SQRM345))/
     /     (2*DSQRT(SQRM2345)))**2-P2(5)**2)/DSQRT(SQRM2345)
C
      AMAX      =DSQRT((((P0(5)**2+P1(5)**2-AMU))/
     /     (2*P0(5)))**2-P1(5)**2)*
     *     *(P0(5)-P1(5))**2*   ! 3-particle phase space
     *     DSQRT(((((P0(5)-P1(5))**2+P2(5)**2-AMU345))/
     /     (2*DSQRT(AMU)))**2-P2(5)**2)/DSQRT(AMU)
C     AMAX=1.
C
      if((AMAX*R3).gt.RHO) goto 10
      write(*,*) RHO ,SQRM2345, AMU345, SQRM345
C
      P2345(5)=DSQRT(SQRM2345)
      P345(5)=DSQRT(SQRM345)
      call DECAYS(P0,P1,P2345)
      call DECAYS(P2345,P2,P345)
      call STAR3T(P345,P3,P4,P5)
      return 
      end
C
      subroutine star5evd(P0,P1,P2,P3,P4,P5)
C
C     first P0->P1+P2345 then P2345->P2+P3+P4+P5
C     input: P*(5) - masses of particles, P0(1-4) - 4-momentum
C     output P*(1-5) - 5-momentum
C     IMPORTANT: P3(5)=P4(5)=P5(5)  DO NOT FORGET!!!
C     cause i'm using star4t only for this case
C     G.I.Kopylov, "Osnovy kinematiki resonansov", str.245, formula 108
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL PP2RND, ENERG
      save NEVENT
      DIMENSION P0(5),P1(5),P2(5),P3(5),P4(5),P5(5),P2345(5)
      DATA PI / 3.1415926 /
      data nevent / 0 /
      NCYCL=0
 10   R1=PP2RND(-1.)
      R2=PP2RND(-1.)
C
C     SA lectures 
      AMU       =(P2(5)+P3(5)+P4(5)+P5(5))**2
      SQRM2345  =AMU+R1*((P0(5)-P1(5))**2-AMU)
C
      ALERT=(((P0(5)**2+P1(5)**2-SQRM2345))/
     /      (2*P0(5)))**2-P1(5)**2
      RHO       =DSQRT((((P0(5)**2+P1(5)**2-SQRM2345))/
     /     (2*P0(5)))**2-P1(5)**2)*
     *     PHASESPACE4(DSQRT(SQRM2345),P2(5),P3(5),P4(5),P5(5))
c      write(*,*) ALERT, SQRM2345
C
      AMAX      =DSQRT((((P0(5)**2+P1(5)**2-AMU)**2)/
     /     (2*P0(5)))**2-P1(5)**2)
C
      NCYCL=NCYCL+1
c     write(*,*) NCYCL
      if(NCYCL.gt.10000) goto 20
      if((AMAX*R2).gt.RHO) goto 10
C
 20   P2345(5)=DSQRT(SQRM2345)
      call DECAYS(P0,P1,P2345)
      call STAR4m(P2345,P2,P3,P4,P5)
      nevent=nevent+1
c     write(*,*) 'nevent=',nevent
      return 
      end
C
      function PHASESPACE4(AM,AM1,AM2,AM3,AM4)
C     G.I.Kopylov, "Osnovy kinematiki resonansov", str.245, formula 108
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      data  AKSI1,    AKSI2,   AKSI3,   AKSI4 /
     /     0.0000,   0.3420, 0.51916, 0.61706 /
      data               C2,      C3,      C4 /
     /               3.1416,  16.340,  43.459 /
      AMU1 = M1
      AMU2 = AMU1+M2
      AMU3 = AMU2+M3
      AMU4 = AMU3+M4
      T4   = AM-AMU4
      T3   = T4*AKSI3
      T2   = T3*AKSI2
      T1   = T2*AKSI1
      AMK4 = T4+AMU4
      AMK3 = T3+AMU3
      AMK2 = T2+AMU2
      AMK1 = T1+AMU1
C
      P2   =DSQRT((2*AMU2+T2+T1)*(2*AMU1+T2+T1)*(2*AM2+T2-T1)*(T2-T1))/
     /     (2*AMU2+2*T2)
      P3   =DSQRT((2*AMU3+T3+T2)*(2*AMU2+T3+T2)*(2*AM3+T3-T2)*(T3-T2))/
     /     (2*AMU3+2*T3)
      P4   =DSQRT((2*AMU4+T4+T3)*(2*AMU3+T4+T3)*(2*AM4+T4-T3)*(T4-T3))/
     /     (2*AMU4+2*T4)
C
      PHASESPACE4=C4*(T4**2)*P4*P3*P2/AM 
      return
      end
C
      subroutine pp2ranlux(X,N)
C
      real X(*),pp2rnd
C
      do I=1,N
         X(I)=pp2rnd(-1.)
      enddo
C      write(*,*) 'pp2ranlux called with N=',n,'x(1)=',x(1)
      return
      end
C
      function Dzero(COSGJ, PHITY)
C     D0 wave distribution
      implicit double precision (A-H,O-Z)
      real*8 Dzero
C
      DzeroVal=(1.5*CosGJ*CosGJ - 0.5)
      DzeroMax=(1.5*1. - 0.5)
      Dzero = DzeroVal*DzeroVal/DzeroMax/DzeroMax
c$$$      write(*,*)'COSGJ, PHITY , Dzero', CosGJ,PHITY,Dzero
      return 
      end
C
      function Dplusminus(COSGJ, PHITY)
C     D0 wave distribution
      implicit double precision (A-H,O-Z)
      real*8 Dplusminus
C
      DplusVal= dsqrt(1.-CosGJ**2)*CosGJ*dcos(PHITY)
      DminusVal=dsqrt(1.-CosGJ**2)*CosGJ*dsin(PHITY)
      D2Max=(1.)*(1.)
      Dplusminus = 0.5*DplusVal*DplusVal/D2Max + 0.5*DminusVal*DminusVal/D2Max
      return 
      end
C
      function UserPolarizationD(COSGJ, PHITY, F2Polarization)
C     D0 wave distribution
      implicit double precision (A-H,O-Z)
      real*8 UserPolarizationD
      real*8 F2Polarization(8)
      data pi / 3.14159265358979/
C
      D0          = F2Polarization(1)
      Dminus      = F2Polarization(2)
      Dplus       = F2Polarization(3)
      D2minus     = F2Polarization(4)
      D2plus      = F2Polarization(5)
      phaseDminus = F2Polarization(6)
      phaseD2minus= F2Polarization(7)
      phaseDplus  = F2Polarization(8)
      if(D0.lt.0) then 
         write(*,*) 'D0 < 0, to be ignored'
         D0=0.
      endif
      if(Dminus.lt.0) then 
         write(*,*) 'Dminus < 0, to be ignored'
         Dminus=0.
      endif
      if(Dplus.lt.0) then 
         write(*,*) 'Dplus < 0, to be ignored'
         Dplus=0.
      endif
      if(D2minus.lt.0) then 
         write(*,*) 'D2minus < 0, to be ignored'
         D2minus=0.
      endif
      if(D2plus.lt.0) then 
         write(*,*) 'D2plus < 0, to be ignored'
         D2plus=0.
      endif
C
      sumDsqr=D0+Dminus+Dplus+D2minus+D2plus
C
      if(sumDsqr.le.0.) then
         UserPolarizationD = 1.
         goto 12
      endif
      SinGJ=dsqrt(1.-COSGJ*COSGJ)
C
      Y0 = dsqrt(5.0/(4.0*pi))*(1.5*CosGJ*CosGJ - 0.5)
      Y0max = dsqrt(5.0/(4.0*pi))*(1.5 - 0.5)
C
      Yminus = -sqrt(2.)*dsqrt(15.0/(8.0*pi))*SinGJ*CosGJ*dcos(PhiTY)
      Yminusmax = sqrt(2.)*dsqrt(15.0/(8.0*pi))
C
      Yplus = -sqrt(2.)*dsqrt(15.0/(8.0*pi))*SinGJ*CosGJ*dsin(PhiTY)
      Yplusmax = sqrt(2.)*dsqrt(15.0/(8.0*pi))
C
      Y2minus= sqrt(2.)*0.25*dsqrt(15.0/(2.0*pi))*SinGJ*SinGJ*dcos(2.*PhiTY)
      Y2minusmax= sqrt(2.)*0.25*dsqrt(15.0/(2.0*pi))
C
      Y2plus= sqrt(2.)*0.25*dsqrt(15.0/(2.0*pi))*SinGJ*SinGJ*dsin(2.*PhiTY)
      Y2plusmax= sqrt(2.)*0.25*dsqrt(15.0/(2.0*pi))
C
      PolarizationD = D0*Y0*Y0 + Dminus*Yminus*Yminus + Dplus*Yplus*Yplus 
     +                  + D2minus*Y2minus*Y2minus + D2plus*Y2plus*Y2plus 
     +                  + 2.*dsqrt(D2minus*D0)*dcos(dabs(phaseD2minus - phaseD0))*Y2minus*Y0
     +                  + 2.*dsqrt(Dminus*D0)*dcos(dabs(phaseDminus - phaseD0))*Y0*Yminus
     +                  + 2.*dsqrt(D2plus*Dplus)*dcos(dabs(phaseD2plus))*Yplus*Y2plus
      PolarizationDmax = D0*Y0max*Y0max + Dminus*Yminusmax*Yminusmax + Dplus*Yplusmax*Yplusmax 
     +                  + D2minus*Y2minusmax*Y2minusmax + D2plus*Y2plusmax*Y2plusmax 
     +                  + dabs(2.*dsqrt(D2minus*D0)*dcos(dabs(phaseD2minus - phaseD0))*Y2minusmax*Y0max)
     +                  + dabs(2.*dsqrt(Dminus*D0)*dcos(dabs(phaseDminus - phaseD0))*Y0max*Yminusmax)
     +                  + dabs(2.*dsqrt(D2plus*Dplus)*dcos(dabs(phaseD2plus))*Yplusmax*Y2plusmax)
c$$$      write(*,*)'COSGJ, PHITY , Dzero', CosGJ,PHITY,Dzero
      UserPolarizationD=PolarizationD/PolarizationDmax
 12   return 
      end
