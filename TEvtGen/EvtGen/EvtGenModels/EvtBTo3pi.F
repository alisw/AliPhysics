C--------------------------------------------------------------------------
C
C Environment:
C      This software is part of the EvtGen package developed jointly
C      for the BaBar and CLEO collaborations.  If you use all or part
C      of it, please give an appropriate acknowledgement.
C
C Copyright Information: See EvtGen/COPYRIGHT
C      Copyright (C) 1998      Caltech, UCSB
C
C Module: EvtBTo3pi.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
C===================================================================
C This package is providing a B -->-- 3pions decay generator
C Its is composed of the following subroutines:
C 
C [*] EVTHowToUse  
C                 This is an How To Use routine where one may find the
C                 implementation of the time dependance: That is to	  
C                 say that it shows how the output of the routine is	
C                 supposed to be used in the mind of the authors.
C
C===================================================================
C [0] EVT3pions
C                 The routine to be called. Note that on the first call
C                 some initialization will be made, in particular the
C                 computation of a normalization factor designed to help
C                 the subsequent time-dependent generation of events.
C                 The normalisation done inside EVT3pions is such that
C                 at the level of the time implementation, the maximum
C                 time-dependant weight is unity : nothing is to be 
C                 computed to generate unity-weight events. The exact
C                 meaning of the above is made explicit in the HowToUse
C                 routine.
C [1] first_step
C                 Generation of the kinematics of the 3 pions 
C                 It uses the function evtranf which is a random number 
C                 generator providing an uniform distribution 
C                 of Real*4 random number between 0 and 1	
C [2] compute
C                 return the amplitudes of the B0     -->-- 3pions
C                                   and of the B0_bar -->-- 3pions
C                 corrected for the generation mechanism
C                 The notations used are the ones of the paper of
C                 A. Snyder and H. Quinn [Phys.Rev.D48 (1993) 2139] 
C [3] BreitWigner
C                 compute the Breit-Wigner of the contributing rho s
C                 taking into account the cosine term linked to the
C                 zero-helicity of the rho s. There is three forms of
C                 Breit-Wigners available. The first one is the simple
C                 non-relativistic form, while the second and third
C                 ones are more involved relativistic expressions which
C                 in addition incorporate the contributions of three
C                 rho resonances [rho(770:1450:1700)]. The parameters
C                 used are the ones resulting from the ALEPH analysis
C                 which might be found in CERN-PPE:97013 (submitted to
C                 Zeitschrift für Physik C). The two parametrizations
C                 of the relativistic Breit-Wigners are the ones of 	
C                 Kuhn-SantaMaria [default] and of Gounaris-Sakurai.            
C                 The default setting is the non-relativistic Breit-
C                 Wigner form. 
C [4] Set_constants
C                 Set the constants values, from the pion mass to the
C                 penguin contributions. It is called by EVT3pions
C
C And other routines which do not deserve comment here.
C===================================================================
c
c      call EvtHowToUse
c      Stop
c      End

      subroutine EvtHowToUse 
      
      Implicit none     
      Real*8  alphaCP
      Integer iset,number,i,j,jset,N_gener,N_asked
      Real*8  p_pi_plus(4),p_pi_minus(4),p_gamma_1(4),p_gamma_2(4)
      Real*8  q_pi_plus(4),q_pi_minus(4),q_gamma_1(4),q_gamma_2(4)
      Real*8  CPweight
      Real*8  Real_B0,Imag_B0,Real_B0bar,Imag_B0bar
      Real*8  t3pions,ttag,t,Weight,Weight_max
      Real*8  AB0,AB0bar,Ainter,xd
      Real*8  m_rhop_2,m_rhom_2
      Real*8  Evtranf,Tag
      Real*8  Kin_moy,Kin_moy2,Kin,Lambda
      
      alphaCP      = 0.4
      Lambda     = Dsin(2.*alphaCP)
      Kin_moy    =0.
      Kin_moy2   =0.
      N_gener    =0
      N_asked    =100000
      
      weight_max = 1.0
      
c xd is not needed by the generator but by the time-dependant 
C sector of the general BaBar Generator.
      xd         = 0.71
      
c run : Simulation of the Anders Ryd Generator

      Do number=1,N_asked ! weighted events as generated here 
      
      If(number.eq.1) then
      iset=10000          ! 10^4 events are used to normalize amplitudes
      Else
      iset=0              ! iset must be reset to zero after the first call 
      End If

c Here is the call to EVT3pions !!!!!!!!!!!!!!!!!!!!!!!!! 
c communication of data is done by argument only <<<<<<<<     
       call EVT3pions(
     + alphaCP,iset,
     + p_pi_plus,p_pi_minus,p_gamma_1,p_gamma_2,
     + Real_B0,Imag_B0,Real_B0bar,Imag_B0bar) 
     
C that is it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
c what follows is just to show how the ouput of EVT3pions
C should be used at the time dependant level.     
c Decay time of the B->-3 pions
 1    Continue
      t3pions=-Dlog(Evtranf())
      If(t3pions.gt.10.) Go to 1
c Decay time of the tagging B
 2    Continue
      ttag   =-Dlog(Evtranf())                
      If(ttag   .gt.10.) Go to 2
c time-difference between the 3pi decay and the tagging decay     
      t=t3pions-ttag      
c select the Tag
      Tag   =evtranf()
C get the relevant quantities            
      AB0   = Real_B0   **2 + Imag_B0   **2
      AB0bar= Real_B0bar**2 + Imag_B0bar**2
      Ainter= Real_B0*Imag_B0bar - Imag_B0*Real_B0bar
c generate acording to the tag 
      If(Tag.gt.0.5) Then
      
c a B0bar tag => at t3pions=ttag the decay is one from a B0      
      Weight= AB0   *Dcos(xd*t/2.D+00)**2
     +      + AB0bar*Dsin(xd*t/2.D+00)**2      
     +      - Ainter*Dsin(xd*t       )          
     
           Else
           
c a B0    tag => at t3pions=ttag the decay is one from a B0bar       
      Weight= AB0bar*Dcos(xd*t/2.D+00)**2
     +      + AB0   *Dsin(xd*t/2.D+00)**2      
     +      + Ainter*Dsin(xd*t       )          
     
      End If
c 	Print*,'weight',weight
      If(Weight.Gt.evtranf()) Then
c	Print*,'unweighted events'
c----------------------------------------------------------------------------         
c unweighted event production
c---------------------------------------------------------------------------- 
      N_gener=N_gener+1  
C here is just a Dalitz plot and a few prints

      m_rhop_2=(p_pi_plus (4)+p_gamma_1(4)+p_gamma_2(4))**2 
      m_rhom_2=(p_pi_minus(4)+p_gamma_1(4)+p_gamma_2(4))**2  
      do j=1,3
      m_rhop_2=m_rhop_2-(p_pi_plus (j)+p_gamma_1(j)+p_gamma_2(j))**2 
      m_rhom_2=m_rhom_2-(p_pi_minus(j)+p_gamma_1(j)+p_gamma_2(j))**2  
      end do
           
c here is a check that weight_max is one

      If(Weight.gt.Weight_max) Then
      Weight_max=Weight
      Print*,' overweighted event found at weight = ',Weight_max           
      End If
      
c here is an example of use of the iset<0 option.
c One does not want to generate a new event but just to
c compute the matrix element of the CP conjugate event:
c the example provided below uses the Kin variable, but
c this kind of facility will be very usefull for other 
c types of applications. For example, to deal with events
c in a likelihood analysis, one needs to compute matrix 
c elements, and it would be nice to be able to do that from
c the BaBar generator itself, to be sure that the Physics
c is the same.

      q_pi_plus (4)=   p_pi_minus(4)
      q_pi_minus(4)=   p_pi_plus (4)
      q_gamma_1 (4)=   p_gamma_1 (4)
      q_gamma_2 (4)=   p_gamma_2 (4)

      do i=1,3
      q_pi_plus (i)= - p_pi_minus(i)
      q_pi_minus(i)= - p_pi_plus (i)
      q_gamma_1 (i)= - p_gamma_1 (i)
      q_gamma_2 (i)= - p_gamma_2 (i)
      end do

c Here is the call to EVT3pions !!!!!!!!!!!!!!!!!!!!!!!!!  
c with the iset<0 option

       jset=-1
           
       call EVT3pions(
     + alphaCP,jset,
     + q_pi_plus,q_pi_minus,q_gamma_1,q_gamma_2,
     + Real_B0,Imag_B0,Real_B0bar,Imag_B0bar) 
c that is it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
C get the relevant quantities            
      AB0   = Real_B0   **2 + Imag_B0   **2
      AB0bar= Real_B0bar**2 + Imag_B0bar**2
      Ainter= Real_B0*Imag_B0bar - Imag_B0*Real_B0bar
      
      If(Tag.lt.0.5) Then
      
c a B0bar tag => at t3pions=ttag the decay is one from a B0      
      CPWeight= AB0   *Dcos(xd*t/2.D+00)**2
     +      + AB0bar*Dsin(xd*t/2.D+00)**2      
     +      - Ainter*Dsin(xd*t       )          
     
           Else
           
c a B0    tag => at t3pions=ttag the decay is one from a B0bar       
      CPWeight= AB0bar*Dcos(xd*t/2.D+00)**2
     +      + AB0   *Dsin(xd*t/2.D+00)**2      
     +      + Ainter*Dsin(xd*t       )          
     
      End If
      
      Kin = (Weight-CPweight)/(Weight+CPweight)/Lambda      
      Kin_moy  = Kin_moy  + Kin
      Kin_moy2 = Kin_moy2 + Kin**2

c----------------------------------------------------------------------------      
      End If    

c end of the loop over events     
      End Do
      
      Print*,'number of unity-weight events generated : ',N_gener
      Print*,'number of trials                        : ',N_asked
      Print*,' measurement of sin(2 alpha)            = ',
     +      Kin_moy/Kin_moy2
      Print*,' true value                             = ',Lambda
      
      End
C===================================================================
      subroutine Evt3pions(
     + alpha_input,iset,
     + p_pi_plus,p_pi_minus,p_gamma_1,p_gamma_2,
     + Real_B0,Imag_B0,Real_B0bar,Imag_B0bar)
c-----------------------------------------------------------------
c	----------------------------------------------------------
c	--- This is the routine to be called by the Main generator
c****************************************************************************
c	--- Inputs are  :
c
c	---	          alpha     : the value of the alpha angle (radians) 
c             Note that the amplitude which are returned incorporate the
c             CP violating phase present in the mixing. 
c
c	--- 		  iset      : to initialize the generator
c	   		              iset should be set to the number
c                                     of events to be used for the pre-run
c                                     devoted to the computation of the
c                                     normalization factor designed to help
c                                     the subsequent time-dependent generation
c                                     of events.
c	   		              After the first call iset should be 
c                                     set to 0.                                
c****************************************************************************    	 		              set to 0.
c       --- Outputs are :
c
c       ---               p_pi_plus : the four momentum of the pi+
c       ---               p_pi_minus: the four momentum of the pi-           
c       ---	          p_gamma_1 : the four momentum of the first  photon   
c       ---		  p_gamma_2 : the four momentum of the second photon
c             Note that : the energy is stored in the fourth component
c                         the values are the ones of the B rest frame
c                         a random rotation has been applied 
c 
c       ---               Real_B0   : The real      part of the amplitude of 
c                                    the B0    ->- 3 pions decay                  
c       ---               Imag_B0   : The imaginary part of the amplitude of 
c                                    the B0    ->- 3 pions decay
c       ---               Real_B0bar: The real      part of the amplitude of 
c                                    the B0bar ->- 3 pions decay
c       ---               Imag_B0bar: The imaginary part of the amplitude of 
c                                    the B0bar ->- 3 pions decay            
c****************************************************************************
c-----------------------------------------------------------------      
      Implicit none
#include "EvtGenModels/EvtBTo3pi.inc"
      Real*8  alpha_input
      Real*8  beta_input
      Integer iset
      Real*8  p_pi_plus(4),p_pi_minus(4),p_gamma_1(4),p_gamma_2(4)
      Real*8  Real_B0,Imag_B0,Real_B0bar,Imag_B0bar
      
c Working quantities
      Integer i,number
      Real*8  p1(5),p2(5),p3(5),Gamma1(5),Gamma2(5) 
      Real*8 factor,factor_max,AB0,AB0bar,Ainter,R1,R2 
      data factor_max/1.D+00/
      Integer ierr
      ierr =0 
c-------------------------------------------------------------------       
      If(iset.eq.0) Then
c------------------------------------------------------------------- 
c     this is the normal mode of operation 
c     First, generate the kinematics  

      p1(5)=M_pip**2
      p2(5)=M_pi0**2
      p3(5)=M_pim**2

 10   continue


      call Evtfirst_step(p1,p2,p3)
      
c     Then, compute the amplitudes 


      Call EvtCompute(p1,p2,p3,
     +    Real_B0,Imag_B0,Real_B0bar,Imag_B0bar,iset,ierr)
      if(ierr.ne.0)  Go To 10
    
c-------------------------------------------------------------------       
      ElseIf(iset.lt.0) Then
c-------------------------------------------------------------------
c     This is an user mode of operation where the kinematics is
c     provided by the user who only wants the corresponding amplitudes
c     to be computed

      Do i=1,4
      p1(i)= p_pi_plus (i)
      p2(i)= p_gamma_1 (i) + p_gamma_2 (i)
      p3(i)= p_pi_minus(i)
      End Do
      p1(5)= M_pip**2
      p2(5)= M_pi0**2
      p3(5)= M_pim**2

	Call EvtCompute(p1,p2,p3,
     +          Real_B0,Imag_B0,Real_B0bar,Imag_B0bar,iset,ierr)
 
      if(ierr.ne.0) Then
	Print*,'the provided kinematics are not physical'
	Print*,'ierr=',ierr
	Print*,'the program will stop'
	Stop
      endif	
c-------------------------------------------------------------------       
      ElseIf(iset.gt.0) Then
c-------------------------------------------------------------------
c     This is the pre-run mode of operation where initializations are
c     performed.
 
      factor_max= 0

C Ghm : beta is not needed for this generation - put a default value
      beta_input = 0.362
C       
      call Evtset_constants(alpha_input, beta_input)
      p1(5)=M_pip**2
      p2(5)=M_pi0**2
      p3(5)=M_pim**2

c     pre-run
      Do number=1,iset
      
 20   continue

      call Evtfirst_step(p1,p2,p3)
 
      Call EvtCompute(p1,p2,p3,Real_B0,Imag_B0,Real_B0bar,
     +  Imag_B0bar,iset,ierr)
         if(ierr.ne.0) Go To 20  
      AB0   = Real_B0   **2 + Imag_B0   **2
      AB0bar= Real_B0bar**2 + Imag_B0bar**2
      Ainter= Real_B0*Imag_B0bar - Imag_B0*Real_B0bar
      R1    =  (AB0-AB0bar) /(AB0+AB0bar)
      R2    =(2.D+00*Ainter)/(AB0+AB0bar)
      factor=(1.D+00+Dsqrt(R1**2+R2**2))*(AB0+AB0bar)/2.D+00
      If(factor.gt.factor_max) factor_max=factor
      
      End Do
c     end of the pre-run 
c     Normalization factor     
      factor_max=1.D+00/Dsqrt(factor_max)
c	Print*,'factor_max',factor_max,p1,p2,p3
c-------------------------------------------------------------------      
      End If
c------------------------------------------------------------------- 
 
      Real_B0   =Real_B0    * factor_max
      Imag_B0   =Imag_B0    * factor_max
      Real_B0bar=Real_B0bar * factor_max
      Imag_B0bar=Imag_B0bar * factor_max
      
      if(iset.lt.0) return

c     P1,p2,p3 ---> random rotation in B rest frame

      Call EvtRotation(p1,1)
      Call EvtRotation(p2,0)
      Call EvtRotation(p3,0)
      
C     Desintegrate the pi_0

      Call EvtGammaGamma(p2,Gamma1,Gamma2)
      
C     Feed the output four vectors

      Do i=1,4
      
      p_pi_plus (i)=p1    (i)      
      p_pi_minus(i)=p3    (i)
      p_gamma_1 (i)=Gamma1(i)
      p_gamma_2 (i)=Gamma2(i)
      
      End Do
      
      Return
      
      End
       
c===================================================================
      subroutine Evtfirst_step(P1,P2,P3)
c-----------------------------------------------------------------
c	----------------------------------------------------------
c	--- This routine generates the 5-vectors P1,P2,P3 
c	--- Associated respectively with the Pi+,Pi0 and Pi-
c	--- 		P1(1) = Px	
c	--- 		P1(2) = Py
c	--- 		P1(3) = Pz
c	--- 		P1(4) = E
c	--- 		P1(5) = M**2
c	----------------------------------------------------------
c 	---	Input Four Vectors                                                            
C 	---	Particle [1] is the pi+                                                       
C 	---	Particle [2] is the pi0                                                       
C 	---	Particle [3] is the pi-                                                       
c	----------------------------------------------------------

c	----------------------------------------------------------      	
c 	--- commons      	
c	----------------------------------------------------------	
#include "EvtGenModels/EvtBTo3pi.inc"
c	----------------------------------------------------------
c 	--- Used Functions 
c	----------------------------------------------------------

	real*8 evtranf

c	----------------------------------------------------------
c 	--- Variables in Argument
c	----------------------------------------------------------

	 real*8 P1(5),P2(5),P3(5)
	 
c	----------------------------------------------------------
c 	--- Local Variables
c	----------------------------------------------------------
	
      	real*8 m12,min_m12, max_m12
      	real*8 m13,min_m13, max_m13
      	real*8 m23,min_m23, max_m23
      	Real*8 cost13,cost12,cost23
      	
      	real*8 p1mom,p2mom,p3mom
      
      	real*8 x, y, z, mass
	integer i
C There is two ways of generating the events:
C The first one used a pole-compensation method to generate the
C events efficiently taking into account the poles due to the
C Breit-Wigners of the rho s. It is activated by setting 
C Phase_Space to .false.     
C The second one generates events according to phase space. It is
C inneficient but allows the exploration of the full Dalitz plot
C in an uniform way. It was found to be usefull fopr some peculiar
C applications. It is activated by setting
C Phase_space to .true.
C Note that in that case, the generation is no longer correct.

	Logical Phase_space	
	data Phase_space/.false./
       	  
c       initialize to avoid warning on Linux.
        mass=0.0

c	----------------------------------------------------------
c	--- Computation
c	----------------------------------------------------------
      	max_m12 = M_B**2 
      	min_m12 = P1(5) + P2(5) 
      
      	max_m13 = M_B**2 
      	min_m13 = P1(5) + P3(5)
      
      	max_m23 = M_B**2 
      	min_m23 = P2(5) + P3(5)
      	

100   	Continue

c	----------------------------------------------------------	
c	--- Generation of the Mass of the Rho(+,-,0) 
c	----------------------------------------------------------

        If(.not.Phase_space) Then
	y = evtranf()*PI - PI/2.
	x = Dtan(y)
	mass = x*Gam_rho/2. +Mass_rho
        End If
	
c	----------------------------------------------------------	
c	--- z is the Flag needed to choose between the generation 
c	--- of a Rho+, a Rho-,or a Rho0
c	----------------------------------------------------------

	z = 3.* evtranf()
	
	if(z.lt.1.) Then

	If(Phase_space) Then
		m12 = evtranf()*(max_m12-min_m12)+min_m12
	Else
		m12 = mass**2
	End If		
		m13 = evtranf()*(max_m13-min_m13)+min_m13
		m23 = MB2 - m12 - m13
	else if (z.lt.2.) Then

	If(Phase_space) Then
		m13 = evtranf()*(max_m13-min_m13)+min_m13
	Else
		m13 =  mass**2
	End If
	 	m12 = evtranf()*(max_m12-min_m12)+min_m12
		m23 = MB2 - m12 - m13
	else

	If(Phase_space) Then
		m23 = evtranf()*(max_m23-min_m23)+min_m23
	Else
		m23 = mass**2
	End If
		m12 = evtranf()*(max_m12-min_m12)+min_m12
      		m13 = MB2 - m12 - m23
	endif

c	----------------------------------------------------------	
c	--- Check that the physics is OK :
c	--- Are the invariant Masses in allowed ranges ?
c	----------------------------------------------------------

      	If(m23.lt.min_m23.or.m23.gt.max_m23) Go to 100
      	If(m13.lt.min_m13.or.m13.gt.max_m13) Go to 100
      	If(m12.lt.min_m12.or.m12.gt.max_m12) Go to 100
      	
c	----------------------------------------------------------
c	--- Are the Cosines of the angles between particles 
c	--- Between -1 and +1 ?
c	----------------------------------------------------------
      
      	P1(4)=(M_B**2+P1(5)-m23)/(2.*M_B)
      	P2(4)=(M_B**2+P2(5)-m13)/(2.*M_B)
      	P3(4)=(M_B**2+P3(5)-m12)/(2.*M_B)
        if (p1(4)**2-P1(5).le.0.0) go to 100
        if (p2(4)**2-P2(5).le.0.0) go to 100
        if (p3(4)**2-P3(5).le.0.0) go to 100
      	p1mom=Dsqrt(p1(4)**2-P1(5))      
      	p2mom=Dsqrt(p2(4)**2-P2(5))
      	p3mom=Dsqrt(p3(4)**2-P3(5))
      	cost13=(2.*p1(4)*p3(4)+P1(5)+p3(5)-m13)/(2.*p1mom*p3mom)
      	cost12=(2.*p1(4)*p2(4)+P1(5)+p2(5)-m12)/(2.*p1mom*p2mom)      
      	cost23=(2.*p2(4)*p3(4)+P2(5)+p3(5)-m23)/(2.*p2mom*p3mom)
      	If(Dabs(cost13).gt.1.) Go to 100
      	If(Dabs(cost12).gt.1.) Go to 100
      	If(Dabs(cost23).gt.1.) Go to 100

c	----------------------------------------------------------
c	--- Filling the 5-vectors P1,P2,P3
c	----------------------------------------------------------      

      	P3(1) = 0
      	P3(2) = 0
      	p3(3) = p3mom

      	P1(3) = p1mom*cost13
      	P1(1) = p1mom*Dsqrt(1.D+00-cost13**2)
      	p1(2) = 0.

      	Do i=1,3
      	P2(i)=-p1(i)-p3(i)
      	End do
12345   continue
	END
	
	
c======================================================================         
      Subroutine EvtCompute(p1,p2,p3,
     +           Real_B0,Imag_B0,Real_B0bar,Imag_B0bar,iset,ierr)                 
c-----------------------------------------------------------------------
      	IMPLICIT None  
#include "EvtGenModels/EvtBTo3pi.inc"

                                                                                                 
      	Real*8 m12, m13, m23, W12, W13, W23, Wtot
      	Real*4 evt_gmas
      	Complex*16 MatBp,MatBm
        Real*8  Real_B0,Imag_B0,Real_B0bar,Imag_B0bar            
        Integer Iset,ierr
      	Real*8     p1(5),p2(5),p3(5)
      	                                              
        COMPLEX*16 Mat_rhop,Mat_rhom,Mat_rho0
        Complex*16 BreitWigner,Mat_1,Mat_2,Mat_3                                         
c	----------------------------------------------------------------      
C 	---	Account for the pole compensation                                                    
c	----------------------------------------------------------------  
        ierr = 0

	if(evt_gmas(p1,p2).lt.0.or.evt_gmas(p1,p3)
     +		.lt.0.or.evt_gmas(p2,p3).lt.0) Then
		ierr=1
		return
	endif	

      	m12 = sqrt(evt_gmas(p1,p2))
      	m13 = sqrt(evt_gmas(p1,p3))
      	m23 = sqrt(evt_gmas(p2,p3))
      
 	W12 = (1./m12)*1./((Mass_rho - m12)**2+(Gam_rho/2.)**2)
      	W13 = (1./m13)*1./((Mass_rho - m13)**2+(Gam_rho/2.)**2)
      	W23 = (1./m23)*1./((Mass_rho - m23)**2+(Gam_rho/2.)**2)
        
	If(iset.ge.0) Then                              
      	Wtot = 1./Dsqrt(W12 + W13 + W23)                                
	else
	Wtot = 1.
	endif 
	                                                                                 
c	----------------------------------------------------------------      
C 	---	Compute Breit-Wigners                                                         
c	----------------------------------------------------------------  

      	Mat_rhop=BreitWigner(p1,p2,p3,ierr)                                            
      	Mat_rhom=BreitWigner(p2,p3,p1,ierr)                                            
      	Mat_rho0=BreitWigner(p1,p3,p2,ierr)   

      	                                  
c	----------------------------------------------------------------      
C 	---	Build up the amplitudes                                                       
c	----------------------------------------------------------------  

        Mat_1 =  Mat_s3*Mat_rhop
      	Mat_2 =  Mat_s4*Mat_rhom
      	Mat_3 =  Mat_s5*Mat_rho0*0.5D+00
                                                                                                                                        
      	MatBp    = (Mat_1+Mat_2+Mat_3) * Wtot 
                                                
      	Mat_1 =  Nat_s3*Mat_rhom
      	Mat_2 =  Nat_s4*Mat_rhop
      	Mat_3 =  Nat_s5*Mat_rho0*0.5D+00
       	                                           
      	MatBm    = (Mat_1+Mat_2+Mat_3) * Wtot
      	
c       Pick up the Real and Imaginary parts

c       changed Real to DBLE and Imag DIMAG to make it complie 
c       with Absoft and g77 (ryd)      	

      	Real_B0    = DBLE(MatBp)
      	Imag_B0    = DIMAG(MatBp) 
      	               	                           
        Real_B0bar = DBLE(MatBm)
        Imag_B0bar = DIMAG(MatBm)
        
      	Return                                                                    
      	End 
      	                                                                                                                    
c======================================================================         
      Function BreitWigner(p1,p2,p3,ierr)
c----------------------------------------------------------------------                                                  
      	IMPLICIT None   
#include "EvtGenModels/EvtBTo3pi.inc"
      	Complex *16 BreitWigner,EvtCRhoF_W 
      	Logical Aleph     
        Integer ierr
      	Data Aleph /.true./                                               

c	---------------------------------------------------------------
c 	---	Input Four Vectors                                                            
c	---------------------------------------------------------------
      	Real*8     p1(5),p2(5),p3(5)                                              
                                           
c	---------------------------------------------------------------
C 	---	intermediate variables                                                        
c	---------------------------------------------------------------
      	Real*8 E12_2,m12_2,beta,gamma,argu,m13_2,costet,coscms,m12                             
      	Real*8 Factor,Num_real,Num_imag                                           
      	Integer i
      	Real *8 p1z,p1zcms12,e1cms12,p1cms12 
                                                                        
c	---------------------------------------------------------------
C 	---	Boost factor                                                                  
c	---------------------------------------------------------------

        BreitWigner=Dcmplx(0,0)  


        ierr  = 0
        	E12_2=(p1(4)+p2(4))**2                                                    
      	m12_2=E12_2                                                               
      	Do i=1,3                                                                  
      		m12_2=m12_2-(p1(i)+p2(i))**2                                              
        End Do                                                                    
      	Argu = 1.D+00 - m12_2 / E12_2                                             

      	If(argu.gt.0.) Then                                                       
      		beta = Dsqrt(Argu)                                                        
        Else                                                       
      		Print *,'Abnormal beta ! Argu  = ',Argu 
      		                               
      		Argu = 0.
      		Beta = 0.                                                                     
        End If                                                     
      	
      	If(m12_2.gt.0.)Then                                                       
      		m12     = Dsqrt(m12_2)                                                    
        Else                                                       
      		Print *,'Abnormal m12  ! m12_2 = ',m12_2 
      		Print*,'p1 = ',p1      		
      		Print*,'p2 = ',p2
      		Print*,'p3 = ',p3                                 
      		Stop                                                                      
        End if 
      	
      	gamma=Dsqrt(E12_2/m12_2)                                                    
c	---------------------------------------------------------------
C 	---	get the cosine in the B CMS                                                   
c	---------------------------------------------------------------

      	m13_2=(p1(4)+p3(4))**2                                                    
      	Do i=1,3                                                                  
      		m13_2=m13_2-(p1(i)+p3(i))**2                                              
      	End Do   
	if(m13_2.lt.0)            Go To 50
        if((p1(4)**2-p1(5)).lt.0) Go To 50
	if((p3(4)**2-p3(5)).lt.0) Go To 50
      	costet= (2.D+00*p1(4)*p3(4)-m13_2+p1(5)+p3(5))                            
     >      /                                                                   
     >        (2.D+00*Dsqrt( (p1(4)**2-p1(5)) * (p3(4)**2-p3(5)) ))              

c	---------------------------------------------------------------
C 	---	get the costet in the 1-2 CMS 
c	---------------------------------------------------------------
      	
      	p1z=dsqrt(P1(4)**2-p1(5))*costet                                               
      	p1zcms12=gamma*(p1z+beta*P1(4))
      	e1cms12 =gamma*(p1(4)+beta*p1z)
      	p1cms12 =Dsqrt(e1cms12**2-p1(5))
      	coscms=p1zcms12/p1cms12   
      	
c	---------------------------------------------------------------      	                 
C 	---	Build the Breit Wigner 
c	---------------------------------------------------------------
      	
      	If(aleph) then
      		BreitWigner= coscms * EvtCRhoF_W( m12_2 )
      	else                                                       
      	 Factor  = 2.D+00* ( (Mass_rho-m12)**2+(0.5D+00*Gam_rho)**2 )              
      	 Factor  = coscms * Gam_rho / Factor                                       
      	 Num_real=  (Mass_rho-m12)                    * Factor                     
      	 Num_imag=                    0.5D+00*Gam_rho * Factor                     
      	 BreitWigner=Dcmplx(Num_real,Num_imag)  
      	End if  

      	Return
 50     continue
	ierr = 2
        Return
      	End
      	                                                                    
c======================================================================         
      Subroutine EvtSet_Constants(balpha, bbeta)   
c----------------------------------------------------------------------                                                     
      IMPLICIT None         
#include "EvtGenModels/EvtBTo3pi.inc"

c	-------------------------------------------------------------	                                                                
c 	---	Basic matrix elements 
c	-------------------------------------------------------------	                                                                
                                                        
      	Complex *16   Mat_Tp0,Mat_Tm0,Mat_T0p,Mat_T0m,Mat_Tpm,Mat_Tmp             
     >             ,Mat_P1,Mat_P0                                               
     >             ,Nat_Tp0,Nat_Tm0,Nat_T0p,Nat_T0m,Nat_Tpm,Nat_Tmp             
     >             ,Nat_P1,Nat_P0, Mat_Ppm, Mat_Pmp                                            
     >		   ,StrongExp                                                                                                                                  
      Real*8 	    balpha, bbeta
      Real*8        calpha,salpha                                               
      Real*8        cbeta,sbeta                                               
      Real*8        StrongPhase 
      Real*8  Rho_gen,Eta_gen,c_factor
c	-------------------------------------------------------------	                                                                                                      
      pi       = 3.141592653 D+00
      alphaCP  = balpha                                               
      calpha   = Dcos(alphaCP)                                             
      salpha   = Dsin(alphaCP)                                         
      betaCP   = bbeta                                               
      cbeta    = Dcos(betaCP)                                             
      sbeta    = Dsin(betaCP)                                         
      M_B      = 5.2794D+00                                           
      M_pip    = 0.13957D0                                                   
      M_pim    = M_pip                                                          
      M_pi0    = 0.134976D0                                                  
      M_Kp     = 0.49368
      Mass_rho = 0.770D0                                                          
      Gam_rho  = 0.150D0
      Mass_Kstarp = 0.8916 
      Gam_Kstarp  = 0.0498
      Mass_Kstar0 = 0.8961 
      Gam_Kstar0  = 0.0505 
      
      MB2         = M_B**2 + M_pip**2 + M_pim**2 + M_pi0**2 
      MC2         = M_B**2 + M_Kp **2 + M_pim**2 + M_pi0**2                        

      StrongPhase = 0.

c	-------------------------------------------------------------	                                                                
      StrongExp = Dcmplx(Dcos(StrongPhase),Dsin(StrongPhase))
      
      rho_gen = 0.D+00
c	if(balpha.ne.0) Then
c      eta_gen = cos(balpha)/sin(balpha)
c	else
c      eta_gen =0.
c	endif	
	eta_gen = 0.5D+00
      c_factor=
     >Sqrt((rho_gen**2+eta_gen**2)/((1.-rho_gen)**2+eta_gen**2))

      Mat_Tp0 = Dcmplx(calpha, -salpha) *1.09D+00        
      Mat_Tm0 = Dcmplx(calpha, -salpha) *1.09D+00                    
      Mat_T0p = Dcmplx(calpha, -salpha) *0.66D+00                    
      Mat_T0m = Dcmplx(calpha, -salpha) *0.66D+00                 
      Mat_Tpm = Dcmplx(calpha, -salpha) *1.00D+00                  
      Mat_Tmp = Dcmplx(calpha, -salpha) *0.47D+00

      Mat_Ppm = -Dcmplx(Dcos(-0.5D0),Dsin(-0.5D0))*0.2D0
      Mat_Pmp =  Dcmplx(Dcos(2.D0),Dsin(2.D0))*0.15D0

      Mat_P1  = 0.5D0*(Mat_Ppm-Mat_Pmp)
      Mat_P0  = 0.5D0*(Mat_Ppm+Mat_Pmp)                                    
      Nat_Tp0 = Dcmplx(calpha, salpha) *1.09D+00                      
      Nat_Tm0 = Dcmplx(calpha, salpha) *1.09D+00                       
      Nat_T0p = Dcmplx(calpha, salpha) *0.66D+00                              
      Nat_T0m = Dcmplx(calpha, salpha) *0.66D+00                     
      Nat_Tpm = Dcmplx(calpha, salpha) *1.00D+00                       
      Nat_Tmp = Dcmplx(calpha, salpha) *0.47D+00                            
      Nat_P1  = Mat_P1                                       
      Nat_P0  = Mat_P0
      
      Mat_Tpm = StrongExp * Mat_Tpm
      Nat_Tpm = StrongExp * Nat_Tpm
                                              
      Mat_S1  = Mat_Tp0 + 2.D+00 * Mat_P1                                       
      Mat_S2  = Mat_T0p - 2.D+00 * Mat_P1                                       
      Mat_S3  = Mat_Tpm + Mat_P1 + Mat_P0                                       
      Mat_S4  = Mat_Tmp - Mat_P1 + Mat_P0                                       
      Mat_S5  = - Mat_Tpm - Mat_Tmp + Mat_Tp0 + Mat_T0p                           
     >        - 2.D+00 * Mat_P0                                                 
                                                                                                                              
      Nat_S1  = Nat_Tp0 + 2.D+00 * Nat_P1                                       
      Nat_S2  = Nat_T0p - 2.D+00 * Nat_P1                                       
      Nat_S3  = Nat_Tpm + Nat_P1 + Nat_P0                                       
      Nat_S4  = Nat_Tmp - Nat_P1 + Nat_P0                                       
      Nat_S5  = - Nat_Tpm - Nat_Tmp + Nat_Tp0 + Nat_T0p                           
     >        - 2.D+00 * Nat_P0      
                                                                        

c=====================================================================           
c=====================================================================


C Gautier 14-Jan-98 : 
C  set matrix elements according to Jerome s calculations
C  (so called Reference set)
c B0    -->-- K*+ pi- Amplitudes (Trees + Penguins)                                    
      MatKstarp = Dcmplx(calpha, -salpha) * Dcmplx( 0.220,0.)
     >          + Dcmplx(cbeta,   sbeta)  * Dcmplx(-1.200,0.) 
c B0    -->-- K*0 pi0 Amplitudes (Trees + Penguins)                                    
      MatKstar0 = Dcmplx(calpha, -salpha) * Dcmplx( 0.015,0.)
     >          + Dcmplx(cbeta,   sbeta)  * Dcmplx( 0.850,0.) 
c B0    -->-- K+ rho- Amplitudes (Trees + Penguins)                                    
      MatKrho   = Dcmplx(calpha, -salpha) * Dcmplx( 0.130,0.)
     >          + Dcmplx(cbeta,   sbeta)  * Dcmplx( 0.160,0.) 
c=====================================================================                  
c B0bar -->-- K*+ pi- Amplitudes (Trees + Penguins)                                    
      NatKstarp = Dcmplx(0.,0.) 
c B0bar -->-- K*0 pi0 Amplitudes (Trees + Penguins)                                    
      NatKstar0 = Dcmplx(0.,0.) 
c B0bar -->-- K+ rho- Amplitudes (Trees + Penguins)                                    
      NatKrho   = Dcmplx(0.,0.)                                     
c=====================================================================                                                                           
c=====================================================================   
c       Print*,'Mat_Tp0',Mat_Tp0
c       Print*,'Mat_T0p',Mat_T0p
c       Print*,'Mat_Tpm',Mat_Tpm
c       Print*,'Mat_Tmp',Mat_Tmp
c       Print*,'Mat_P1',Mat_P1
c       Print*,'Mat_P0',Mat_P0

      Return                                                                    
      End   
      
      
      
C=====================================================================
      	subroutine EvtGammaGamma(P2,PG1Boost,PG2Boost)
c---------------------------------------------------------------------
c	--------------------------------------------------------------
c	---	The aim of this routine is to generate two photons 
c	---	for the decay of the pi0
c	--------------------------------------------------------------
	implicit none
c --- commons

#include "EvtGenModels/EvtBTo3pi.inc"
	
c --- Used Functions 

	Real*8 evtranf	

c --- Variables in Argument

	Real*8 PG1Boost(5), PG2Boost(5)
	Real*8 P2(5)

c --- Local Variables	

	Real*8 EGammaCmsPi0
	Real*8 PGamma1(5), PGamma2(5)
	Real*8 sinThetaRot, cosThetaRot, PhiRot
	Real*8 Beta(4)
	Integer i
	
	Do i=1, 3
          Beta(i) = P2(i)/p2(4)
        Enddo

        EGammaCmsPi0 = Dsqrt(P2(5))/2.
        
        
        cosThetaRot = evtranf()*2.D+00-1.D+00
        sinThetaRot = Dsqrt(1.D+00 - cosThetaRot**2)
        PhiRot = evtranf()*2.D+00*pi
        
        PGamma1(1) = EGammaCmsPi0 * sinThetaRot * Dcos(PhiRot)
        PGamma1(2) = EGammaCmsPi0 * sinThetaRot * Dsin(PhiRot)
        PGamma1(3) = EGammaCmsPi0 * cosThetaRot
        PGamma1(4) = EGammaCmsPi0
        PGamma1(5) = 0.
        
        PGamma2(1) = - PGamma1(1)
        PGamma2(2) = - PGamma1(2)
        PGamma2(3) = - PGamma1(3)
        PGamma2(4) =   PGamma1(4)
        PGamma2(5) = 0.
        
        Call EvtREFREF(BETA,PGamma1,PG1Boost,1)
        Call EvtREFREF(BETA,PGamma2,PG2Boost,1)
        
        end
        
C===================================================================        
       SUBROUTINE EvtREFREF(BETA,PIN,POUT,S)
c------------------------------------------------------------------- 
c	------------------------------------------------------------       
c	---	 GIVEN A FINAL STATE , TO BOOST IT IN ITS C.M.S TAKE
c 	---	 BETA  =-SIGMA(P(K))/SIGMA(E(K))
c 	---	 S     =+1.
c 	---	 TO GO BACK TO THE INITIAL FRAME : SET
c 	---	 S     =-1.
c	------------------------------------------------------------

	implicit None
	Real*8 Beta,Pin,Pout
	Integer s,i
      	dimension BETA(4),PIN(5),POUT(5)                                   
      	Real*8 GAMMA,BETPIN,BETA2
      
        
      	BETA2 =BETA(1)**2+BETA(2)**2+BETA(3)**2
      	IF(BETA2.GT.1.D0) Print*,'Beta**2 .gt.1!!!'
      	IF(BETA2.GT..99999D0) BETA2=0.99999
      	BETA(4)=DSQRT(BETA2)
      	GAMMA=1./DSQRT(1.D0-BETA2)
      	
      	BetPin = 0.
      	Do i=1,3
      		BetPin = Beta(i)*Pin(i) + BetPin
      	Enddo
      	BetPin = s*BetPin	
	
      	DO 1 I=1,3
      	POUT(I)=PIN(I)+S*BETA(I)*GAMMA*
     >      		(PIN(4)+GAMMA/(GAMMA+1.)*BETPIN)
  1   	CONTINUE
      	POUT(4) = GAMMA*(PIN(4)+BETPIN)
      	POUT(5) = PIN(5)
      	            	
      	RETURN
      	END      
                                                                                                                               
C******************************************************************
	function evt_gmas(P1,P2)
C******************************************************************
c	-----------------------------------------------------------	
c	--- 	compute the invariant mass between particle 1 and 2
c	-----------------------------------------------------------
	implicit none
	real*8 p1(5),P2(5)
	real*4 evt_gmas
	integer i
	
	evt_gmas = (P1(4)+P2(4))**2
	
	do i=1,3
		evt_gmas = evt_gmas - (P1(i) + P2(i))**2
	enddo
	
	end      	
      	
C=================================================================== 
      	Subroutine EvtRotation(p,new)      
c------------------------------------------------------------------
      	Implicit None
#include "EvtGenModels/EvtBTo3pi.inc"

      	Real*8 p(5),pstor(3)
      	Real*8 phi2,phi3
      	Real*8 c1,c2,c3,s1,s2,s3
      	Real*8 MatRot(3,3)
        save MatRot
      	Real*8 evtranf
      	Integer new,i,j
      
      	If(new.ne.0) Then
c	---------------------------------------------------------      
C 	--- generate a random rotation 
c	--------------------------------------------------------- 
     
      		phi2=evtranf()*2.*pi
      		phi3=evtranf()*2.*pi
      
                c1=2.D+00*evtranf()-1.D+00
       		c2=Dcos(phi2)
      		c3=Dcos(phi3)
      
                s1=Dsqrt(1.D+00-c1**2)
      		s2=Dsin(phi2)
      		s3=Dsin(phi3)
      		
c	---------------------------------------------------------      
C 	---	Compute the overall rotation matrix
c	---------------------------------------------------------      
      		MatRot(1,1)=c1
      		MatRot(1,2)=s1*c3
      		MatRot(1,3)=s1*s3
      		MatRot(2,1)=-s1*c2
      		MatRot(2,2)=c1*c2*c3-s2*s3
      		MatRot(2,3)=c1*c2*s3+s2*c3
      		MatRot(3,1)=s1*s2
      		MatRot(3,2)=-c1*s2*c3-c2*s3
      		MatRot(3,3)=-c1*s2*s3+c2*c3
      	End If
      	
c	---------------------------------------------------------      	
C 	---	store the input 3-momentum      
c	---------------------------------------------------------
      	Do i=1,3
      		pstor(i)=p(i)
      		p(i) = 0. 
      	End Do
c	---------------------------------------------------------      	
C 	---	Apply the rotation      
c	---------------------------------------------------------      	
      	Do i=1,3
      	  Do j=1,3
      		p(i)=p(i)+MatRot(i,j)*pstor(j)
      	  End do
      	End do
      
      	end
           
	

c======================================================================         
      FUNCTION EvtCRhoF_W( s )
C -----------------------------------------------------------------------
C      
C Author : Andreas Hoecker [Aleph Collaboration, BaBar]
C
C tau --> pi pi0 nu
C
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 EvtCRhoF_W, EvtcBW_KS, EvtcBW_GS

      PARAMETER (AmPi2=0.13956995**2)

      LOGICAL KUHN_SANTA    / .TRUE. / !...type of Breit-Wigner formula
      double precision lambda
	 COMMON /LAM/ lambda
C
C store parameters
C
      IF (KUHN_SANTA) THEN
C...rho(770)            
        AmRho     =  0.7734
        GamRho    =  0.1477
C...rho(1450)
        AmRhoP    =  1.465
        GamRhoP   =  0.696
        beta      =  -0.229
C...rho(1700)
        AmRhoPP   =  1.760
        GamRhoPP  =  0.215
        gamma     =  0.075
C...exponent
        lambda    =  1.0
      ELSE
C...rho(770)            
        AmRho     =  0.7757
        GamRho    =  0.1508
C...rho(1450)
        AmRhoP    =  1.448
        GamRhoP   =  0.503
        beta      =  -0.161
C...rho(1700)
        AmRhoPP   =  1.757
        GamRhoPP  =  0.237
        gamma     =  0.076
C...exponent
        lambda    =  1.0
      ENDIF
C
C init
C
      EvtCRhoF_W = DCMPLX(0.,0.)

      IF (KUHN_SANTA) THEN
C
C Kuehn-Santamaria model
C
        EvtCRhoF_W =  
     >       ( EvtcBW_KS(s,AmRho**2,GamRho)         +     !...BW-rho( 770)
     >       EvtcBW_KS(s,AmRhoP**2,GamRhoP)*(beta)  +     !...BW-rho(1450)
     >       EvtcBW_KS(s,AmRhoPP**2,GamRhoPP)*(gamma) ) / !...BW-rho(1700)
     >       (1.+beta+gamma)
      ELSE
C
C Gounaris-Sakurai model
C
        EvtCRhoF_W =  
     >       ( EvtcBW_GS(s,AmRho**2,GamRho)         + 
     >       EvtcBW_GS(s,AmRhoP**2,GamRhoP)*(beta)  + 
     >       EvtcBW_GS(s,AmRhoPP**2,GamRhoPP)*(gamma) ) / 
     >       (1.+beta+gamma)
      ENDIF

      RETURN
      END
C
C * ================================================================== *
C
C
C * ------------------------------------------------------------- *
C
      FUNCTION EvtcBW_KS( s,Am2,Gam )

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 EvtcBW_KS
      PARAMETER (AmPi2=0.13956995**2)
 
      DOUBLE PRECISION lambda
      COMMON /LAM/ lambda

      EvtcBW_KS = 0.
      IF (s.le.4.*AmPi2) RETURN

      G  =  Gam* (Am2/s)**lambda * ((s-4.*AmPi2)/(Am2-4.*AmPi2))**1.5
      D  =  (Am2-s)**2 + s*G**2
      X  =  Am2*(Am2-s)
      Y  =  Am2*SQRT(s)*G

      EvtcBW_KS = DCMPLX(X/D,Y/D)
      
      RETURN
      END
C
C * ------------------------------------------------------------- *
C
      FUNCTION EvtcBW_GS( s,Am2,Gam )

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 EvtcBW_GS
      PARAMETER (AmPi2=0.13956995**2)

      DOUBLE PRECISION Evtfs,d,N

      DOUBLE PRECISION lambda
      COMMON /LAM/ lambda

      EvtcBW_GS = 0.
      IF (s.le.4.*AmPi2) RETURN

      G  =  Gam* (Am2/s)**lambda * ((s-4.*AmPi2)/(Am2-4.*AmPi2))**1.5
      z1 =  Am2 - s + Evtfs(s,Am2,Gam)
      z2 =  SQRT(s)*G
      z3 =  Am2 + d(Am2)*Gam*SQRT(Am2)

      X  =  z3 * z1
      Y  =  z3 * z2
      N  =  z1**2 + z2**2

      EvtcBW_GS = DCMPLX(X/N,Y/N)
      
      RETURN
      END

C---------------------------------------------------------------------

      FUNCTION Evtfs(s,AmRho2,GamRho)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION k,k_s,k_Am2
      k_s   = k(s)
      k_Am2 = k(AmRho2)

      Evtfs = GamRho * AmRho2 / k_Am2**3 *
     >     ( 
     >          k_s**2 * (h(s)-h(AmRho2)) +
     >          (AmRho2-s) * k_Am2**2 * dh_ds(AmRho2)
     >     )
      RETURN
      END
C --------------
      FUNCTION k(s)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION k
      PARAMETER (AmPi2=0.13956995**2)
      k = 0.5 * SQRT(s - 4.*AmPi2)
      RETURN
      END
C --------------
      FUNCTION h(s)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION k,k_s
      PARAMETER (pi=3.141593)
      PARAMETER (AmPi=0.13957)
      SQRTs = SQRT(s)
      k_s = k(s)
      h = 2./pi * (k_s/SQRTs) * LOG((SQRTs+2.*k_s)/(2.*AmPi))
	RETURN
      END
C --------------
      FUNCTION dh_ds(s)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION k
      PARAMETER (pi=3.141593)
      dh_ds = h(s) * (1./(8.*k(s)**2) - 1./(2.*s)) + 1./(2.*pi*s)
	RETURN
      END
C --------------
      FUNCTION d(AmRho2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION k,k_AmRho2
      PARAMETER (pi=3.141593)
      PARAMETER (AmPi=0.13956995)
      PARAMETER (AmPi2=AmPi**2)
      AmRho    = SQRT(AmRho2)
      k_AmRho2 = k(AmRho2)
      d = 3./pi * AmPi2/k_AmRho2**2 * 
     >         LOG((AmRho+2.*k_AmRho2)/(2.*AmPi)) + 
     >    AmRho/(2.*pi*k_AmRho2) - AmPi2*AmRho/(pi*k_AmRho2**3)
     	 RETURN
      END



