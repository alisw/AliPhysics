      Complex*16 FUNCTION FPIBEL(W,flpar) 
C**************************************************************************** 
C                       PION FORM FACTOR FROM BELLE 
C
C            Rho(770)+rho(1450)+rho(1700) model is used as in:
C             M.Fujikawa et al., Phys.Rev.D 78 (2008) 072006
C                          -------------------
C       * flpar: = 0 - (all free) fit result for par(1...11) 
C                = 1 - (par(1)=F_pi(0)=1-fixed) fit result for par(2...11) 
C
C       * parameter : par(1)=overall norm. factor (=F_pi(0)) 
C                   : par(2)=rho(770) mass
C                   : par(3)=rho(770) width 
C                   : par(4)=rho(1450) mass
C                   : par(5)=rho(1450) width 
C                   : par(6)=rho(1450) admixture abs. value (|beta|)
C                   : par(7)=rho(1450) admixture phase (phi=arg(beta))
C                   : par(8)=rho(1700) mass
C                   : par(9)=rho(1700) width 
C                   : par(10)=rho(1700) admixture abs. value (|gamma|)
C                   : par(11)=rho(1700) admixture phase (phi3=arg(gamma)) 
C
C       * x:=s=(M_pipi0)^2=W^2 
C
C       * Notice: the following code was extracted (and checked) directly 
C                 from the Belle fit function with minor cosmetic changes 
C                 (this is internal comment of Belle)
C
C       * Notice: the following code was extracted (and checked) directly 
C                 from  tauola_with_belle_fpi.F sent by D. Epifanov and
C                 it is an identical copy of 
C                 Complex*16 FUNCTION FPIBEL(W,flpar) 
C                 used in tauola_with_belle_fpi.F
C
C       * Called: CURR_PIPI0  if (FF2PIRHO = 2  or FF2PIRHO =3) 
C                 in ../RChL-currents/value_parameter.f
C                 FF2PIRHO = 2 is for flpar = 0 
C                 FF2PIRHO = 3 is for flpar = 1  
C
C**************************************************************************** 
      implicit none 
      integer flpar 
      real W 
      real*8 x,par(11),pi,pimas,pi0mas,taumas
      real*8 beta,berho1,berho2,berho3 
      real*8 ps,prho1,prho2,prho3
      real*8 gamma1,gamma2,gamma3 
      real*8 hs,h1,h2,h3,dhds1,dhds2,dhds3 
      real*8 d1,d2,d3,fs1,fs2,fs3 
      real*8 a1,b1,c1,bw_re1,bw_im1
      real*8 a2,b2,c2,bw_re2,bw_im2 
      real*8 a3,b3,c3,bw_re3,bw_im3 
      real*8 sinphi,cosphi,sinphi3,cosphi3 
      complex*16 mbw1,mbw2,mbw3,mff,mephi,mephi3
      parameter(pi=3.141592653589)

      x=dble(W**2)
C------------- PARAMETERS --------------------
C     Particle mass (unit: GeV)

      pimas  = 0.1395702
      pi0mas = 0.13498
      taumas = 1.77699

C---------------------------------------------
C Opt. par. from Table VII of PRD78 (2008) 072006

      if(flpar.eq.0) then 
         par(1) = 1.02 
         par(2) = 0.7749 
         par(3) = 0.1486 
         par(4) = 1.428 
         par(5) = 0.413 
         par(6) = 0.13 
         par(7) = 197 
         par(8) = 1.694 
         par(9) = 0.135 
         par(10)= 0.028 
         par(11)= -3 
      else 
         par(1) = 1.0 
         par(2) = 0.7746 
         par(3) = 0.1481 
         par(4) = 1.446 
         par(5) = 0.434 
         par(6) = 0.15 
         par(7) = 202 
         par(8) = 1.728 
         par(9) = 0.164 
         par(10)= 0.037 
         par(11)= 24 
      endif   

C========= Beta (s, rho(770), rho(1450), rho(1700)) ================== 

       beta  =sqrt( (1.0 - (pimas - pi0mas)**2/x)
     +             *(1.0 - (pimas + pi0mas)**2/x) )

c      ----- rho(770) -----
       berho1=sqrt( ( 1.0 - (pimas - pi0mas)**2/(par(2)*par(2)) )
     +             *( 1.0 - (pimas + pi0mas)**2/(par(2)*par(2)) ) )

c      ----- rho(1450) -----
       berho2=sqrt( ( 1.0 - (pimas - pi0mas)**2/(par(4)*par(4)) )
     +             *( 1.0 - (pimas + pi0mas)**2/(par(4)*par(4)) ) )

c      ----- rho(1700) -----
       berho3=sqrt( ( 1.0 - (pimas - pi0mas)**2/(par(8)*par(8)) )
     +             *( 1.0 - (pimas + pi0mas)**2/(par(8)*par(8)) ) )

C========= Momentum (s, rho(770), rho(1450), rho(1700)) ==============
C
C       ps    : momentum of pi at s
C       prho1 : momentum of pi at s=rho(770)**2
C       prho2 : momentum of pi at s=rho(1450)**2
C       prho3 : momentum of pi at s=rho(1700)**2

       ps  = 0.5 * sqrt(x) * beta
c      ----- rho(770) -----
       prho1 = 0.5* par(2) * berho1
c      ----- rho(1450) -----
       prho2 = 0.5* par(4) * berho2
c      ----- rho(1700) -----
       prho3 = 0.5* par(8) * berho3

C========= Width (rho(770), rho(1450), rho(1700)) ====================

C      ----- rho(770) -----
       gamma1 = par(3) * ( par(2)*par(2)/x ) * (ps/prho1)**3
C      ----- rho(1450) -----
       gamma2 = par(5) * ( par(4)*par(4)/x ) * (ps/prho2)**3
C      ----- rho(1700) -----
       gamma3 = par(9) * ( par(8)*par(8)/x ) * (ps/prho3)**3

C============== h(rho(770), rho(1450), rho(1700)) ====================
C
C       hs : h(s)...s(s=x)
C       h1 : h(s) for s=rho(770)**2
C       h2 : h(s) for s=rho(1450)**2
C       h3 : h(s) for s=rho(1700)**2

        hs = (2.0/pi)*(ps/sqrt(x))*
     +       log( (sqrt(x) + 2.0*ps)/(2.0*pimas) )
C      ----- rho(770) ----- 
        h1 = (2.0/pi)*(prho1/par(2))*
     +       log( (par(2) + 2.0*prho1)/(2.0*pimas) )
C      ----- rho(1450) ----- 
        h2 = (2.0/pi)*(prho2/par(4))*
     +       log( (par(4) + 2.0*prho2)/(2.0*pimas) )
C      ----- rho(1700) ----- 
        h3 = (2.0/pi)*(prho3/par(8))*
     +       log( (par(8) + 2.0*prho3)/(2.0*pimas) )

C============ dhds(rho(770), rho(1450), rho(1700)) ===================
C
C       dhds = dh/ds = h(M_rho)[1/(8*p(M_rho)^2) - 1/(2M_rho^2)
C                                              + 1/(2pi*M_rho^2)
C       dhds1 : dh/ds for M_rho=rho(770)
C       dhds2 : dh/ds for M_rho=rho(1450) 
C       dhds3 : dh/ds for M_rho=rho(1700)

C      ----- rho(770) ----- 
        dhds1 = h1*( 1.0/(8.0*prho1**2) - 1.0/(2.0*par(2)**2) )
     +                       +  1.0/(2.0*pi*par(2)**2)
C      ----- rho(1450) ----- 
        dhds2 = h2*( 1.0/(8.0*prho2**2) - 1.0/(2.0*par(4)**2) )
     +                       +  1.0/(2.0*pi*par(4)**2)
C      ----- rho(1700) ----- 
        dhds3 = h3*( 1.0/(8.0*prho3**2) - 1.0/(2.0*par(8)**2) )
     +                       +  1.0/(2.0*pi*par(8)**2)

C============ d(rho(770), rho(1450), rho(1700)) ====================== 
C
C       d = 3/pi * (m_pi^2/p(M_rho)^2) *
C            ln {(M_rho+2*p(M_rho))/2M_pi} + M_rho/(2pi*p(M_rho))   
C                               - (m_pi^2 * M_rho)/(pi*p(M_rho)^3)       
C       d1 : d for M_rho=rho(770)
C       d2 : d for M_rho=rho(1450) 
C       d2 : d for M_rho=rho(1700) 

C      ----- rho(770) ----- 
        d1 = (3.0/pi)*(pimas**2/prho1**2) 
     +         * log((par(2) + 2.0*prho1)/(2.0*pimas))
     +         + par(2)/(2.0*pi*prho1) 
     +         - ((pimas**2)*par(2))/(pi*prho1**3)
C      ----- rho(1450) ----- 
        d2 = (3.0/pi)*(pimas**2/prho2**2) 
     +         * log((par(4) + 2.0*prho2)/(2.0*pimas))
     +         + par(4)/(2.0*pi*prho2) 
     +         - ((pimas**2)*par(4))/(pi*prho2**3)
C      ----- rho(1700) ----- 
        d3 = (3.0/pi)*(pimas**2/prho3**2) 
     +         * log((par(8) + 2.0*prho3)/(2.0*pimas)) 
     +         + par(8)/(2.0*pi*prho3) 
     +         - ((pimas**2)*par(8))/(pi*prho3**3)

C================ f(s) (rho(770), rho(1450), rho(1700) ) ==========
C 
C       f(s) = gamma-rho * M_rhp**2/p(M_rho)[ p(s)^2(h(s)-h(M_rho)
C                                 + (M_rho^2 -s)p(M_rho)^2 * dh/ds]
C       fs1 : f(s) for M_rho=rho(770)
C       fs2 : f(s) for M_rho=rho(1450) 
C       fs3 : f(s) for M_rho=rho(1700) 
C  
C      ----- rho(770) ----- 
        fs1 = par(3) * (par(2)**2/prho1**3) *
     +        ( ps**2 *(hs - h1) + (par(2)**2 -x)*prho1**2 * dhds1 )
C      ----- rho(1450) ----- 
        fs2 = par(5) * (par(4)**2/prho2**3) *
     +        ( ps**2 *(hs - h2) + (par(4)**2 -x)*prho2**2 * dhds2 )
C      ----- rho(1700) ----- 
        fs3 = par(9) * (par(8)**2/prho3**3) *
     +        ( ps**2 *(hs - h3) + (par(8)**2 -x)*prho3**2 * dhds3 )

C====== BW form G&S model (rho(770) + rho(1450) + rho(1700)) ====== 
C    bw_re - real part ; bw_im - imaginary part 

C      ----- rho(770) -----
       a1=par(2)*par(2) - x + fs1
       b1=sqrt(x)*gamma1
       c1=par(2)*par(2) + d1*par(2)*par(3)

       bw_re1= (a1*c1)/(a1**2 + b1**2) 
       bw_im1= (b1*c1)/(a1**2 + b1**2) 
       mbw1=dcmplx(bw_re1,bw_im1)

c      ----- rho(1450) -----
       a2=par(4)*par(4) - x + fs2
       b2=sqrt(x)*gamma2
       c2=par(4)*par(4)+ d2*par(4)*par(5)

       bw_re2= (a2*c2)/(a2**2 + b2**2) 
       bw_im2= (b2*c2)/(a2**2 + b2**2) 
       mbw2=dcmplx(bw_re2,bw_im2)

c      ----- rho(1700) -----
       a3=par(8)*par(8) - x + fs3
       b3=sqrt(x)*gamma3
       c3=par(8)*par(8)+ d3*par(8)*par(9)

       bw_re3= (a3*c3)/(a3**2 + b3**2) 
       bw_im3= (b3*c3)/(a3**2 + b3**2) 
       mbw3=dcmplx(bw_re3,bw_im3)

C====== admixtures of rho(1450) and rho(1700) ====== 

C --------- exp(i*arg(beta)) ---------- 
       sinphi=sin(par(7)*(pi/180.0))  ! use degree
       cosphi=cos(par(7)*(pi/180.0))  ! use degree
       mephi=dcmplx(cosphi,sinphi)

C --------- exp(i*arg(gamma)) ---------
       sinphi3=sin(par(11)*(pi/180.0))  ! use degree
       cosphi3=cos(par(11)*(pi/180.0))  ! use degree
       mephi3=dcmplx(cosphi3,sinphi3)

C================= Form factor ===================== 

       mff = 1.0/(1.0+par(6)*mephi+par(10)*mephi3) 
       FPIBEL=par(1)*mff*(mbw1+par(6)*mephi*mbw2+par(10)*mephi3*mbw3) 
      RETURN
      END
