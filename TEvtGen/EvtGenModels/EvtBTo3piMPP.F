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
C Module: EvtBTo3piMPP.F
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
C [*] HowToUse  
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
C [1] first_step_MPP
C                 Generation of the kinematics of the 3 pions 
C                 It uses the function ranf which is a random number 
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
c      Implicit none     
c      Real Hmemor
c      Common/Pawc/Hmemor(1000000)
c      Call Hlimit(1000000)
c      call EvtHowToUse
c      call EvtHowToUse_MPP
c      call EvtHowToUse_P00
c      Stop
c      End

      subroutine EvtHowToUse_MPP 
      
      Implicit none     
      Real*8  alpha
      Integer iset,number,j,N_gener,N_asked


      Real*8  p_pi_minus(4),p_pi_plus_1(4),p_pi_plus_2(4)
      Real*8  Real_Bp,Imag_Bp,Real_Bm,Imag_Bm
      Real*8  Weight,Weight_max
      Real*8  ABp,ABm
      Real*8  m_rho12,m_rho13
      Real*8  Evtranf,Tag
      alpha      = 0.4
      N_gener    = 0
      N_asked    = 100000
      
      weight_max = 1.0
           
c run : Simulation of the Anders Ryd Generator

      Do number=1,N_asked ! weighted events as generated here 
      
      If(number.eq.1) then
      iset=10000          ! 10^4 events are used to normalize amplitudes
      Else
      iset=0              ! iset must be reset to zero after the first call 
      End If

c Here is the call to EVT3pions_MPP  !!!!!!!!!!!!!!!! 
c communication of data is done by argument only <<<<<<<<     
       call EVT3pionsMPP(
     + alpha,iset,
     + p_pi_minus,p_pi_plus_1,p_pi_plus_2,
     + Real_Bp,Imag_Bp,Real_Bm,Imag_Bm) 
     
C that is it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
c select the Tag 
      Tag   =evtranf()
C get the relevant quantities            
      ABp   = Real_Bp   **2 + Imag_Bp   **2
      ABm   = Real_Bp   **2 + Imag_Bm   **2
c generate acording to the tag 
      If(Tag.gt.0.5) Then
      
c a Bm tag =>  the decay is one from a Bp      
      Weight= ABp   
      
      Else
           
c a Bp tag =>  the decay is one from a Bm 
      Weight= ABm
      
      End If
      If(Weight.Gt.evtranf()) Then
c----------------------------------------------------------------------------         
c unweighted event production
c---------------------------------------------------------------------------- 
      N_gener=N_gener+1  
C here is just a Dalitz plot and a few prints

      m_rho12=(p_pi_minus(4)+P_pi_plus_1(4))**2 
      m_rho13=(p_pi_minus(4)+P_pi_plus_2(4))**2 
      do j=1,3
      m_rho12=m_rho12-(p_pi_minus(j)+P_pi_plus_1(j))**2 
      m_rho13=m_rho13-(p_pi_minus(j)+P_pi_plus_2(j))**2  
      end do                 
		     
c here is a check that weight_max is one

      If(Weight.gt.Weight_max) Then
      Weight_max=Weight
      Print*,' overweighted event found at weight = ',Weight_max           
      End If
      
c----------------------------------------------------------------------------      
      End If    

c end of the loop over events     
      End Do
      
      Print*,'number of unity-weight events generated : ',N_gener
      Print*,'number of trials                        : ',N_asked
      
      
      End
C===================================================================
      subroutine Evt3pionsMPP(
     + alpha_input,iset,
     + p_pi_minus,p_pi_plus_1,p_pi_plus_2,
     + Real_Bp,Imag_Bp,Real_Bm,Imag_Bm)
c-----------------------------------------------------------------
c	----------------------------------------------------------
c	--- This is the routine to be called by the Main generator
c           to get the decay of B+ -->-- 3 pions
c       --- AND
c           to get the decay of B- -->-- 3 pions
C       For the sake of clarity, signs refers to B+ decay
c****************************************************************************  	 		              set to 0.
c       --- Outputs are :
c
c       ---               p_pi_minus     : the four momentum of the pi- 
c       ---	          p_pi_plus_1    : the four momentum of the first  pi+
c       ---	          p_pi_plus_2    : the four momentum of the second pi+   
c
c             Note that : the energy is stored in the fourth component
c                         the values are the ones of the B rest frame
c                         a random rotation has been applied 
c 
c       ---               Real_Bp   : The real      part of the amplitude of 
c                                     the B+    ->- 3 pions decay                  
c       ---               Imag_Bp   : The imaginary part of the amplitude of 
c                                     the B+    ->- 3 pions decay
c 
c       ---               Real_Bm   : The real      part of the amplitude of 
c                                     the B-    ->- 3 pions decay                  
c       ---               Imag_Bm   : The imaginary part of the amplitude of 
c      
c****************************************************************************
c-----------------------------------------------------------------      
      Implicit none
#include "EvtGenModels/EvtBTo3pi.inc"
      Real*8  alpha_input
      Integer iset
      Real*8  p_pi_minus(4),p_pi_plus_1(4),p_pi_plus_2(4)
      Real*8  Real_Bp,Imag_Bp,Real_Bm,Imag_Bm
      
c Working quantities
      Integer i,number
      Real*8  p1(5),p2(5),p3(5)
      Real*8  factor_max,ABp,ABm
      Integer ierr
      data factor_max/1.D+00/
      ierr =0
c-------------------------------------------------------------------       
      If(iset.eq.0) Then
c------------------------------------------------------------------- 
c     this is the normal mode of operation 
c     First, generate the kinematics  

      p1(5)= M_pim**2
      p2(5)= M_pip**2
      p3(5)= M_pip**2
          
 10   continue	
      call Evtfirst_step_MPP(p1,p2,p3)
      
c     Then, compute the amplitudes 
      
      Call EvtCompute_MPP(p1,p2,p3,
     +    Real_Bp,Imag_Bp,Real_Bm,Imag_Bm,iset,ierr)
      if(ierr.ne.0) Go To 10
c-------------------------------------------------------------------       
      ElseIf(iset.lt.0) Then
c-------------------------------------------------------------------
c     This is an user mode of operation where the kinematics is
c     provided by the user who only wants the corresponding amplitudes
c     to be computed

      Do i=1,4
      p1(i)= p_pi_minus (i)
      p2(i)= p_pi_plus_1(i) 
      p3(i)= p_pi_plus_2(i) 
      End Do
      p1(5)= M_pim**2
      p2(5)= M_pip**2
      p3(5)= M_pip**2
      
      Call EvtCompute_MPP(p1,p2,p3,
     +    Real_Bp,Imag_Bp,Real_Bm,Imag_Bm,iset,ierr)
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
c    	changed by ryd April 27 1998. this causes an 
c       unaligned access probelm on OSF.
c       (should Beta not be passed as an argument?)    
c      call Evtset_constants(alpha_input,0.362)
      call Evtset_constants(alpha_input,0.362D0)
      p1(5)=M_pim**2
      p2(5)=M_pip**2
      p3(5)=M_pip**2
      
c     pre-run
      Do number=1,iset
      
 20   continue
      call Evtfirst_step_MPP(p1,p2,p3)
 
      Call EvtCompute_MPP(p1,p2,p3,
     +    Real_Bp,Imag_Bp,Real_Bm,Imag_Bm,iset,ierr)     
       if(ierr.ne.0) Go To 20  
      ABp   = Real_Bp   **2 + Imag_Bp   **2
      ABm   = Real_Bm   **2 + Imag_Bm   **2
      
      If(ABp.gt.factor_max) factor_max=ABp
      If(ABm.gt.factor_max) factor_max=ABm
      
      End Do
c     end of the pre-run 

      factor_max=1.D+00/Dsqrt(factor_max)

c-------------------------------------------------------------------      
      End If
c------------------------------------------------------------------- 
    
      Real_Bp   =Real_Bp    * factor_max
      Imag_Bp   =Imag_Bp    * factor_max
      Real_Bm   =Real_Bm    * factor_max
      Imag_Bm   =Imag_Bm    * factor_max
       
      if(iset.lt.0) return	
c     P1,p2,p3 ---> random rotation in B rest frame

      Call EvtRotation(p1,1)
      Call EvtRotation(p2,0)
      Call EvtRotation(p3,0)
           
C     Feed the output four vectors

      Do i=1,4
      
      p_pi_minus     (i)=p1    (i)      
      p_pi_plus_1    (i)=p2    (i) 
      p_pi_plus_2    (i)=p3    (i) 
      
      End Do
      
      Return
      
      End
       
c===================================================================
      subroutine Evtfirst_step_MPP(P1,P2,P3)
c-----------------------------------------------------------------
c	----------------------------------------------------------
c	--- This routine generates the 5-vectors P1,P2,P3 
c	--- Associated respectively with the Pi+ and two Pi0 s 
c	--- 		P1(1) = Px	
c	--- 		P1(2) = Py
c	--- 		P1(3) = Pz
c	--- 		P1(4) = E
c	--- 		P1(5) = M**2
c	----------------------------------------------------------
c 	---	Input Four Vectors                                                            
C 	---	Particle [1] is the pi-                                                       
C 	---	Particle [2] is the pi+                                                      
C 	---	Particle [3] is the pi+                                                     
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

	Logical Phase_space	
	data Phase_space/.false./
       	  
c       initialize to avoid warning on linux
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
c	--- Generation of the Mass of the Rho(0) 
c	----------------------------------------------------------

        If(.not.Phase_space) Then      
	y = evtranf()*PI - PI/2.
	x = Dtan(y)
	mass = x*Gam_rho/2. +Mass_rho
        End If
	
c	----------------------------------------------------------	
c	--- z is the Flag needed to choose between the generation 
c	--- of a Rho0 = pi- pi+[1] or pi- pi+[2]
c	----------------------------------------------------------

	z = evtranf()
	
	if(z.lt..5) Then

		If(Phase_space) Then
		m12 = evtranf()*(max_m12-min_m12)+min_m12
		Else
	 	m12 = mass**2
		End If	
			
	m13 = evtranf()*(max_m13-min_m13)+min_m13
	m23 = MB2 - m12 - m13
	
	else 

		If(Phase_space) Then
		m13 = evtranf()*(max_m13-min_m13)+min_m13
		Else
		m13 =  mass**2
		End If
		
	m12 = evtranf()*(max_m12-min_m12)+min_m12
	m23 = MB2 - m12 - m13

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
      	
      	p1mom=p1(4)**2-P1(5)      
      	p2mom=p2(4)**2-P2(5)
      	p3mom=p3(4)**2-P3(5)
      	If(p1mom.lt.0) Go to 100
      	If(p2mom.lt.0) Go to 100
      	If(p3mom.lt.0) Go to 100
      	p1mom=Dsqrt(p1mom)
      	p2mom=Dsqrt(p2mom)
      	p3mom=Dsqrt(p3mom)
      	
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
      	
	END
	
	
c======================================================================         
      Subroutine EvtCompute_MPP(p1,p2,p3,
     +           Real_Bp,Imag_Bp,Real_Bm,Imag_Bm,iset,ierr)                 
c-----------------------------------------------------------------------
      	IMPLICIT None  
#include "EvtGenModels/EvtBTo3pi.inc"

                                                                                                 
      	Real*8 m12, m13, W12, W13, Wtot
      	Real*8 evt_gmas
      	Complex*16 MatBp,MatBm
        Real*8  Real_Bp,Imag_Bp,Real_Bm,Imag_Bm                                                     
      	Real*8     p1(5),p2(5),p3(5)
        real*8  ASHQ
        Integer iset,ierr
        Data    ASHQ/0.707107 D+00/
      	                                              
        COMPLEX*16 Mat_rhop
        Complex*16 BreitWigner
        
        ierr = 0
c	----------------------------------------------------------------      
C 	---	Account for the pole compensation                                                    
c	----------------------------------------------------------------  
        if(evt_gmas(p1,p2).lt.0.or.evt_gmas(p1,p3)
     +		.lt.0) Then
		ierr=1
		return
	endif	

      	m12 = sqrt(evt_gmas(p1,p2))
      	m13 = sqrt(evt_gmas(p1,p3))
 
 	W12 = (1./m12)*1./((Mass_rho - m12)**2+(Gam_rho/2.)**2)
      	W13 = (1./m13)*1./((Mass_rho - m13)**2+(Gam_rho/2.)**2)
                     
        If(iset.ge.0) Then  
           Wtot = 1.D+00/Dsqrt(W12 + W13)                                        
        else
           Wtot = 1.
        endif
c	----------------------------------------------------------------      
C 	---	Compute Breit-Wigners                                                         
c	----------------------------------------------------------------  
    
      	Mat_rhop = BreitWigner(p1,p2,p3,ierr)
     +   + BreitWigner(p1,p3,p2,ierr)                                                                               
      	                                  
c	----------------------------------------------------------------      
C 	---	Build up the amplitudes                                                       
c	---------------------------------------------------------------- 
c       The factor ASHQ = 1./sqrt(2) is here just to stick to the notations
c       used by Art Snyder and Helen Quinn. It is irrelevant for the genera-
c       tion of events done here. 
                                                                                                                                       
      	MatBp    = Mat_s2 * Mat_rhop * Wtot * ASHQ 
                                                      	                                           
      	MatBm    = Nat_s2 * Mat_rhop * Wtot * ASHQ
      	
c       Pick up the Real and Imaginary parts
c       changed Real to DBLE and Imag DIMAG to make it complie 
c       with Absoft and g77 (ryd)          	

      	Real_Bp    = DBLE(MatBp)
      	Imag_Bp    = DIMAG(MatBp) 
      	               	                           
        Real_Bm    = DBLE(MatBm)
        Imag_Bm    = DIMAG(MatBm)
        
      	Return                                                                    
      	End 
