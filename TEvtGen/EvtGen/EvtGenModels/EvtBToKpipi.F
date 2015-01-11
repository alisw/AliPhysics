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
C Module: EvtBToKpipi.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
C===================================================================
C This package is providing a neutral B -->-- K pi pi decay generator
C Its is composed of the following subroutines:
C 
C [*] HowToUse  
C                 This is an How To Use routine where one may find the
C                 implementation of the time dependance: That is to	  
C                 say that it shows how the output of the routine is	
C                 supposed to be used in the mind of the authors.
C
C===================================================================
C [0] EVTKpipi
C                 The routine to be called. Note that on the first call
C                 some initialization will be made, in particular the
C                 computation of a normalization factor designed to help
C                 the subsequent time-dependent generation of events.
C                 The normalisation done inside EVTKpipi is such that
C                 at the level of the time implementation, the maximum
C                 time-dependant weight is unity : nothing is to be 
C                 computed to generate unity-weight events. The exact
C                 meaning of the above is made explicit in the HowToUse
C                 routine.
C [1] first_step_Kpipi
C                 Generation of the kinematics of the 3 prongs
C                 It uses the function evtranf which is a random number 
C                 generator providing an uniform distribution 
C                 of Real*4 random number between 0 and 1	
C [2] compute
C                 return the amplitudes of the B0 -->-- K+pi-pi0
C                 corrected for the generation mechanism.
c                 Note that this is a Tagging Mode. The CP conjugate 
c                 mode (B0bar -->-- K-pi+pi0) is treated at the same time. 
C [3] BreitWigner
C                 compute the Breit-Wigner of the contributing K* and rho
C                 taking into account the cosine term linked to the
C                 zero-helicity of the spin-1 resonances. There is two 
c                 forms of Breit-Wigners available. The first one is the 
c                 simple non-relativistic form, while the second is the
C                 relativistic expressions.  
C                 The default setting is the relativistic one.
C [4] Set_constants
C                 Set the constants values, from the pion mass to the
C                 penguin contributions. It is called by EVTKpipi
C
C And other routines which do not deserve comment here.
C===================================================================
c      Implicit none     
C      Real Hmemor
C      Common/Pawc/Hmemor(1000000)
C      Call Hlimit(1000000)
C      Call EvtHowToUse_Kpipi
C      Stop
C      End
      

      subroutine EvtHowToUse_Kpipi 
      
      Implicit none     
      Real*8  alphaCP, betaCP
      Integer iset,number,j,N_gener,N_asked

      Real*8  p_K_plus(4),p_pi_minus(4)
      Real*8  p_gamma_1(4),p_gamma_2(4)

      Real*8  Real_B0,Imag_B0,Real_B0bar,Imag_B0Bar
      Real*8  Weight,Weight_max
      Real*8  m_Kstarp,m_Kstar0,m_rhom
      Real*8  Wrong
      Real*8  Evtranf,Tag
      
      alphaCP      = 1.35
      betaCP       = 0.362
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

c Here is the call to EVTKpipi  !!!!!!!!!!!!!!!! 
c communication of data is done by argument only <<<<<<<<     
       call EVTKpipi(
     + alphaCP,betaCP,iset,
     + p_K_plus,p_pi_minus,
     + p_gamma_1,p_gamma_2,
     + Real_B0,Imag_B0,real_B0bar,Imag_B0bar)     
C that is it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
c select the Tag 
      Tag   =evtranf()
c generate acording to the tag 
      If(Tag.gt.0.5) Then
      
c a B0bar tag =>  the decay is one from a B0      
      Weight= Real_B0      **2 + Imag_B0      **2
      
      Else
           
c a B0    tag =>  the decay is one from a B0bar 
      Weight= Real_B0bar   **2 + Imag_B0bar   **2
      
      End If
		      
      If(Weight.Gt.evtranf()) Then
c----------------------------------------------------------------------------         
c unweighted event production
c---------------------------------------------------------------------------- 
      N_gener=N_gener+1  
C here is just a Dalitz plot and a few prints

      m_Kstarp=(p_K_plus   (4)+p_gamma_1(4)+p_gamma_2(4))**2 
      m_rhom  =(p_pi_minus (4)+p_gamma_1(4)+p_gamma_2(4))**2 
      m_Kstar0=(p_K_plus   (4)+p_pi_minus (4))           **2
      do j=1,3
      m_Kstarp=m_Kstarp -(p_K_plus   (j)+p_gamma_1(j)+p_gamma_2(j))**2 
      m_rhom  =m_rhom   -(p_pi_minus (j)+p_gamma_1(j)+p_gamma_2(j))**2
      m_Kstar0=m_Kstar0 -(p_K_plus   (j)+p_pi_minus (j))           **2        
      end do                 

c here is the Dalitz plot when assuming the pion mass for the Kaon
c the energy is redefined
      Wrong =0.  
      do j=1,3
      Wrong = Wrong + p_K_plus(j)**2
      end do     
      Wrong = Dsqrt(Wrong+0.139**2)
 
      m_Kstarp=(Wrong         +p_gamma_1(4)+p_gamma_2(4))**2 
      m_Kstar0=(Wrong         +p_pi_minus (4))           **2
      do j=1,3
      m_Kstarp=m_Kstarp -(p_K_plus   (j)+p_gamma_1(j)+p_gamma_2(j))**2 
      m_Kstar0=m_Kstar0 -(p_K_plus   (j)+p_pi_minus (j))           **2        
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
      subroutine EvtKpipi(
     + alpha_input,beta_input,iset,
     + p_K_plus,p_pi_minus,
     + p_gamma_1,p_gamma_2,
     + Real_B0,Imag_B0,Real_B0bar,Imag_B0bar)
c-----------------------------------------------------------------
c	----------------------------------------------------------
c	--- This is the routine to be called by the Main generator
c           to get the decay of B0    -->-- K+ pi- pi0
c           The decay proceeeds through three channels:
c           a) B0 -->-- K*+ pi-  ; K*+    -->-- K+ pi0
c           b)          K*0 pi0  ; K*0bar -->-- K+ pi-
c           c)          K-  rho+ ; rho+   -->-- pi+ pi0
c           .) The K0 rho0 channel is not implemented since it does 
c              not lead to Kpipi final state, but it is interesting 
c              in itself.
c       It provides at the same time the CP conjugate decay
c                               B0bar -->-- K- pi+ pi0
c****************************************************************************  	 		              
c       --- Outputs are :
c
c       ---               p_K_plus      : the four momentum of the K+ 
c       ---               p_pi_minus    : ........................ pi-
c       ---               p_gamma_1     : ........................ first  gamma of the pi0
c       ---               p_gamma_2     : ........................ second gamma of the pi0
c
c             Note that : the energy is stored in the fourth component
c                         the values are the ones of the B rest frame
c                         a random rotation has been applied 
c 
c       ---               Real_B0      : The real      part of the amplitude of the decay
c                                                                 B0    -->-- K+ pi- pi0
c       ---               Imag_B0      : ... imaginary ..................................
c            similarly      
c       ---               Real_B0bar   : The real      part of the amplitude of the decay 
c                                                                 B0bar -->-- K- pi+ pi0     
c       ---               Imag_B0bar   : ... imaginary ..................................

c****************************************************************************
c-----------------------------------------------------------------      
      Implicit none
#include "EvtGenModels/EvtBTo3pi.inc"
      Real*8  alpha_input, beta_input
      Integer iset
      Real*8  p_K_plus(4),p_pi_minus(4)
      Real*8  p_gamma_1(4),p_gamma_2(4)
      Real*8  Real_B0,Imag_B0,Real_B0bar,Imag_B0Bar
      
c Working quantities
      Integer i,number
      Real*8  p1(5),p2(5),p3(5)
      Real*8  Gamma1(5),Gamma2(5)
      Real*8  factor_max,ABp,ABm
      Integer ierr
	data factor_max/1.D+00/
	ierr =0          
c-------------------------------------------------------------------       
      If(iset.eq.0) Then
c------------------------------------------------------------------- 
c     this is the normal mode of operation 
c     First, generate the kinematics  

      p1(5)= M_Kp **2
      p2(5)= M_pim**2
      p3(5)= M_pi0**2      
          
 10   continue
      call Evtfirst_step_Kpipi(p1,p2,p3)
      
c     Then, compute the amplitudes 
      
      Call EvtCompute_Kpipi(p1,p2,p3,
     +    Real_B0,Imag_B0,Real_B0bar,Imag_B0Bar,iset,ierr)
	if(ierr.ne.0 ) Go To 10
c----------------nedit EvtBto---------------------------------------------------       
      ElseIf(iset.lt.0) Then
c-------------------------------------------------------------------
c     This is an user mode of operation where the kinematics is
c     provided by the user who only wants the corresponding amplitudes
c     to be computed

      Do i=1,4
      p1(i)= p_K_plus  (i)
      p2(i)= p_pi_minus(i)      
      p3(i)= p_gamma_1 (i) + p_gamma_2 (i)
      End Do
      p1(5)= M_Kp **2
      p2(5)= M_pim**2
      p3(5)= M_pi0**2
      
      Call EvtCompute_Kpipi(p1,p2,p3,
     +    Real_B0,Imag_B0,Real_B0bar,Imag_B0Bar,iset,ierr)
      
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
      call Evtset_constants(alpha_input, beta_input)
      p1(5)= M_Kp **2
      p2(5)= M_pim**2
      p3(5)= M_pi0**2      

c     pre-run
      Do number=1,iset
      
  20  continue
      call Evtfirst_step_Kpipi(p1,p2,p3)
 
      Call EvtCompute_Kpipi(p1,p2,p3,
     +    Real_B0,Imag_B0,Real_B0bar,Imag_B0Bar,iset,ierr)     
        if(ierr.ne.0) Go To 20  
      ABp   = Real_B0      **2 + Imag_B0      **2
      ABm   = Real_B0bar   **2 + Imag_B0Bar   **2
      
      If(ABp.gt.factor_max) factor_max=ABp
      If(ABm.gt.factor_max) factor_max=ABm
      
      End Do
c     end of the pre-run 

      factor_max=1.D+00/Dsqrt(factor_max)

c-------------------------------------------------------------------      
      End If
c------------------------------------------------------------------- 
    
      Real_B0   =Real_B0    * factor_max
      Imag_B0   =Imag_B0    * factor_max
      Real_B0bar=Real_B0bar * factor_max
      Imag_B0Bar=Imag_B0Bar * factor_max
      
      if(iset.lt.0) return
c     P1,p2,p3 ---> random rotation in B rest frame

      Call EvtRotation(p1,1)
      Call EvtRotation(p2,0)
      Call EvtRotation(p3,0)
      
C     Desintegrate the pi_0 s
	
      Call EvtGammaGamma(p3,Gamma1,Gamma2)
      
C     Feed the output four vectors

      Do i=1,4
      
      p_K_plus  (i)=p1    (i)
      p_pi_minus(i)=p2    (i)      
      p_gamma_1 (i)=Gamma1(i)
      p_gamma_2 (i)=Gamma2(i)
      End Do

      Return
      
      End
       
c===================================================================
      subroutine Evtfirst_step_Kpipi(P1,P2,P3)
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
C 	---	Particle [1] is the K+                                                       
C 	---	Particle [2] is the pi-                                                       
C 	---	Particle [3] is the pi0                                                      
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
       	  
c	----------------------------------------------------------
c	--- Computation
c	----------------------------------------------------------
      	max_m12 = M_B**2 
      	min_m12 = P1(5) + P2(5) + 2.*Dsqrt(p1(5)*p2(5)) 
      
      	max_m13 = M_B**2 
      	min_m13 = P1(5) + P3(5) + 2.*Dsqrt(p1(5)*p3(5)) 
      
      	max_m23 = M_B**2 
      	min_m23 = P2(5) + P3(5) + 2.*Dsqrt(p2(5)*p3(5)) 
      	
100   	Continue

c	----------------------------------------------------------	
c	--- Generation of the Mass of the Rho or Kstar
c	----------------------------------------------------------


c	----------------------------------------------------------	
c	--- z is the Flag needed to choose between the generation 
c	--- of a K*+, K*0 or rho- resonance
c	----------------------------------------------------------

	z = 3.*evtranf()

        MC2 = M_B**2
	
	If(z.lt.1.) Then
c K*+ production
		If(Phase_space) Then
		m13 = evtranf()*(max_m13-min_m13)+min_m13
		Else
    	        y = evtranf()*PI - PI/2.
	        x = Dtan(y)
	        mass = x*Gam_Kstarp/2. +Mass_Kstarp		
	 	m13 = mass**2
		End If	
			
	m12 = evtranf()*(max_m12-min_m12)+min_m12
	m23 = MC2 - m12 - m13
	
	ElseIf(z.lt.2.) Then
c K*0 production
		If(Phase_space) Then
		m12 = evtranf()*(max_m12-min_m12)+min_m12
		Else
    	        y = evtranf()*PI - PI/2.
	        x = Dtan(y)
	        mass = x*Gam_Kstar0/2. +Mass_Kstar0				
		m12 =  mass**2
		End If
		
	m13 = evtranf()*(max_m13-min_m13)+min_m13
	m23 = MC2 - m12 - m13
	
	Else		
c rho- production
		If(Phase_space) Then
		m23 = evtranf()*(max_m23-min_m23)+min_m23
		Else
    	        y = evtranf()*PI - PI/2.
	        x = Dtan(y)
	        mass = x*Gam_rho/2. +Mass_rho				
		m23 =  mass**2
		End If
		
	m13 = evtranf()*(max_m13-min_m13)+min_m13
	m12 = MC2 - m23 - m13
	
	Endif

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
      Subroutine EvtCompute_Kpipi(p1,p2,p3,
     +           Real_B0,Imag_B0,Real_B0bar,Imag_B0Bar,iset,ierr)                 
c-----------------------------------------------------------------------
      	IMPLICIT None  
#include "EvtGenModels/EvtBTo3pi.inc"

                                                                                                 
      	Real*8 m12, m13, m23, W12, W13, W23, Wtot
      	Real*8 evt_gmas
      	Complex*16 MatB0,MatB0bar,BW12,BW13,BW23
        Real*8  Real_B0,Imag_B0,Real_B0bar,Imag_B0Bar                                                     
      	Real*8     p1(5),p2(5),p3(5)

      	                                              
	Integer ierr,iset
        Complex*16 BrightWagner,BreitWigner
        
	ierr = 0
c	----------------------------------------------------------------      
C 	---	Account for the pole compensation                                                    
c	---------------------------------------------------------------- 
 
      	m12 = evt_gmas(p1,p2)
      	m13 = evt_gmas(p1,p3)
      	m23 = evt_gmas(p2,p3)

	if(m12.lt.0. .or. m13.lt.0. .or. m23.lt.0.) Then
		ierr=1
		Print*,'ierr = ',ierr
		return
	endif	
	
      	m12 = sqrt(m12)
      	m13 = sqrt(m13)
      	m23 = sqrt(m23)
      
 	W12 = 1. / (((Mass_Kstar0 - m12)**2+(Gam_Kstar0/2.)**2)*m12)
      	W13 = 1. / (((Mass_Kstarp - m13)**2+(Gam_Kstarp/2.)**2)*m13)      	
      	W23 = 1. / (((Mass_rho    - m23)**2+(Gam_rho   /2.)**2)*m23)      	
      	
        if(iset.ge.0) Then                   
      	Wtot = 1.D+00/Dsqrt(W12 + W13 + W23)                                        
            else
	Wtot =1.
	Endif                                                                     
c	----------------------------------------------------------------      
C 	---	Compute Breit-Wigners 
c	----------------------------------------------------------------      

        BW13=BrightWagner(p1,p3,p2,Mass_Kstarp,Gam_Kstarp,ierr)
        If(ierr.ne.0) Return
        BW12=BrightWagner(p1,p2,p3,Mass_Kstar0,Gam_Kstar0,ierr)
        If(ierr.ne.0) Return        
c If the rho is to be treated on the same footing as K* ==> use the line below
c        BW23=BrightWagner(p2,p3,p1,Mass_Rho   ,Gam_Rho   ,ierr)
        BW23=BreitWigner(p2,p3,p1,ierr)
        If(ierr.ne.0) Return
        
c	----------------------------------------------------------------                                                                      
c        -              and
C 	---	Build up the amplitudes                                                       
c	---------------------------------------------------------------- 
 
c Here come the relative amplitudes of the three decay channels 
c First, one computes the B0     decay B0    ->- K+ pi- pi0                                                                                                                                        
      	MatB0    = MatKstarp * BW13
     >           + MatKstar0 * BW12                                                     	                                           
     >           + MatKrho   * BW23
c Second, one computes the B0bar decay B0bar ->- K- pi+ pi0                                                                                                                                       
      	MatB0bar = NatKstarp * BW13
     >           + NatKstar0 * BW12                                                      	                                           
     >           + NatKrho   * BW23
          
              	
c       Pick up the Real and Imaginary parts
      	
      	Real_B0    = dreal(MatB0   )*Wtot 
      	Imag_B0    = Imag(MatB0   )*Wtot 
      	               	                           
        Real_B0bar = dreal(MatB0bar)*Wtot
        Imag_B0Bar = Imag(MatB0bar)*Wtot
        
      	Return                                                                    
      	End 
c======================================================================         
      Function BrightWagner(p1,p2,p3,Mass,Width,ierr)
c----------------------------------------------------------------------                                                  
      	IMPLICIT None   
#include "EvtGenModels/EvtBTo3pi.inc"
      	Complex *16 BrightWagner,EvtRBW
        Integer ierr
        Logical RelatBW
        Data    RelatBW/.true./
        Real*8  Mass,Width,Am2Min

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

        BrightWagner=Dcmplx(0,0)  

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
      	
      	If(RelatBW) then
      	        Am2Min       = p1(5) + p2(5) + 2.*Dsqrt( p1(5)*p2(5) )
      		BrightWagner = coscms * EvtRBW(m12_2,Mass**2,Width,Am2Min)
      	Else                                                       
      	 Factor  = 2.D+00* ( (Mass-m12)**2+(0.5D+00*Width)**2 )              
      	 Factor  = coscms * Width / Factor                                       
      	 Num_real=  (Mass-m12)                    * Factor                     
      	 Num_imag=                  0.5D+00*Width * Factor                     
      	 BrightWagner=Dcmplx(Num_real,Num_imag)  
      	End if  

      	Return
 50     continue
	ierr = 2
        Return
      	End
      	                                                                    
      	        
      FUNCTION EvtRBW(s,Am2,Gam,Am2Min)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 EvtRBW

      EvtRBW = 0.
      IF (s.le.Am2Min) RETURN

      G  =  Gam* (Am2/s) * ((s-Am2Min)/(Am2-Am2Min))**1.5
      D  =  (Am2-s)**2 + s*G**2
      X  =  Am2*(Am2-s)
      Y  =  Am2*SQRT(s)*G

      EvtRBW = DCMPLX(X/D,Y/D)
      
      RETURN
      END      	                                                                                                                    
