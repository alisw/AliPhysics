C **********************************************************************
C                                                  Version 2.0
C     Generator of e+e- pairs produced in            1.10.2002
C                  PbPb collisions at LHC             
C                                            
C     Copyright (c) 2002
C     Authors:          Kai Hencken	     <k.hencken@unibas.ch>
C                       Yuri Kharlov	     <Yuri.Kharlov@cern.ch>
C                       Serguei Sadovsky     <sadovsky@mx.ihep.su>,
C                                            <Serguei.Sadovski@cern.ch>
C
C Permission to use, copy and distribute this software and its
C documentation strictly for non-commercial purposes is hereby granted
C without fee, provided that the above copyright notice appears in all
C copies and that both the copyright notice and this permission notice
C appear in the supporting documentation. The authors make no claims
C about the suitability of this software for any purpose. It is
C provided "as is" without express or implied warranty.
C Any change of the code should be submitted to the authors
C
C     Long write up:   K.Hencken,Yu.Kharlov,S.Sadovsky
C                      Internal ALICE Note 2002-27 
C                      
C **********************************************************************
      subroutine ee_init(Ymin,Ymax,PTmin,PTmax)
C-----------------------------------------------------------------------
C    Generator initialisation
C
C    Input variables:  Ymin      - minimal value of rapidity  \
C                      Ymax      - maximal value of rapidity   | of kinematics 
C                      PTmin     - Pt minimum in MeV/c	       | range
C                      PTmax     - Pt maximum in MeV/c        /
C-----------------------------------------------------------------------
      implicit real*8 (A-H,O-Z)
      external DsdYpY,DsdYmY,DsdXpX,DsdXmX,DsdPhi,DsdXX,DsdYY

      common/parPh / parPhi(7)
      common/parYp / parYpY(6)
      common/parYm / parYmY(6)
      common/parXp / parXpX(5)
      common/parXm / parXmX(4)

      common/eevent/ Xsect2,Dsect2, Xsecttot,Dsecttot, Nevnt
      common/eepars/ Xmin,Xmax,YpYmin,YpYmax,WYpYmax,YmYmin,YmYmax,
     &               Ymed1,Ymed2,sgmY1,sgmY2,XYsect,
     &               Gaus1,Gaus2,Gauss,Exp1,Exp2,Exp3,
     &               XmXmin,XmXmax,XpXmin,XpXmax,Exmx1,Exmx2,Xmed,
     &               Sgm1,AnorX,Icase

      real*8 mass  

      data pi     / 3.141 592 653 589 793 238 462 643d00 /

C
C  Exact differential Cross section:
C - - - - - - - - - - - - - - - - - - -
C
      gm = 2750.0d0             ! Pb gamma factor at the LHC
      mass = 0.5109991d0        ! electron mass 
      call initdiffcross (gm,mass)
C     
C  Cross sections:
C - - - - - - - - - -         
      Nevnt  = 0
      Xsect2 = 0.
      Dsect2 = 0.
      
      Nsd = 1
      Npt = 25                  !  7,25,64
      Eps = 0.00005

      Xmin=    Dlog10(Ptmin)
      Xmax=    Dlog10(Ptmax)
      write(*,*) ' Kinematical limits:'
      write(*,*) '          Xmin,Xmax=',Xmin,Xmax 
      write(*,*) '          Ymin,Ymax=',Ymin,Ymax
      write(*,*) 

      Xsec1  = Dtrint(DsdXX,Nsd,Npt,Eps,Xmin,Xmin,Xmin,Xmax,Xmax,Xmin)
      Xsec2  = Dtrint(DsdXX,Nsd,Npt,Eps,Xmax,Xmax,Xmin,Xmax,Xmax,Xmin)
      XsecX  = Xsec1+Xsec2
      
      Ysec1  = Dtrint(DsdYY,Nsd,Npt,Eps,Ymin,Ymin,Ymin,Ymax,Ymax,Ymin)
      Ysec2  = Dtrint(DsdYY,Nsd,Npt,Eps,Ymax,Ymax,Ymin,Ymax,Ymax,Ymin)
      XsecY  = Ysec1+Ysec2

      IF (Xsec1*Xsec2.le.0.or.Ysec1*Ysec2.le.0.) then
         write(*,*) ' Error: insufficient accuracy of XY-cross sections'  
         write(*,*) ' Xsec1,Xsec2,Xsec=', Xsec1,Xsec2,XsecX
         write(*,*) ' Ysec1,Ysec2,Ysec=', Ysec1,Ysec2,XsecY
      endif
      
      XYsect = XsecX*XsecY
      write(*,*) ' Xsections: Xsec1,Xsec2,XsecX=', Xsec1,Xsec2,XsecX
      write(*,*) '            Ysec1,Ysec2,XsecY=', Ysec1,Ysec2,XsecY
C-    write(*,*) '                  XsecX*YsecY=', XYsect
C
      XYsect = 2371.5239*(82./137.035)**4*XYsect ! Normalization factor 
      write(*,*) ' Normalized:      XsecX*YsecY=', XYsect
      write(*,*) 
C      
C-		   XsecPhi  = Dgauss(DsdPhi,0.0d0,  pi, eps)
C
C  Rapidity initialisation:
C - - - - - - - - - - - - - -  
      if (Ymin.ge.Ymax) then
         write(*,*) 'Wrong values of Ymin,Ymax:',Ymin,Ymax
         stop
      endif	  
CYpY-
C
      YpYmin = 2.*Ymin
      YpYmax = 2.*Ymax
      Yp0    = 0.
      if (YpYmin*YpYmax.gt.0.) then 
         if(YpYmin.gt.0.) Yp0 = YpYmin	           
         if(YpYmin.le.0.) Yp0 = YpYmax
      endif
      WYpYmax = DsdYpY(Yp0)
C     
C-    write(*,*)' YpY: YpYmin,YpYmax,WYpYmax =',YpYmin,YpYmax,WYpYmax

CYmY-
C
      YmYmin = Ymin-Ymax
      YmYmax = Ymax-Ymin
C-    write(*,*) ' YmY:  YmYmin,YmYmax=',YmYmin,YmYmax
C     
      Ymed1  = 0.18
      Ymed2  = 4.00
      
      sgmY1  = 1./dsqrt(2.*parYmY(2))
      sgmY2  = 1./dsqrt(2.*parYmY(4))
C     
      eps = 0.000001
      Gaus1 = Dgauss(DsdYmY ,0.0d0, Ymed1 , eps)
      Gaus2 = Dgauss(DsdYmY ,Ymed1, Ymed2 , eps)
      Expon = Dgauss(DsdYmY ,Ymed2, YmYmax, eps)
      Summ  = Gaus1+Gaus2+Expon
      Gaus2 =(Gaus1+Gaus2)/Summ
      Gaus1 = Gaus1/Summ
      Expon = 1.
C     
C     Phi initialisation:
C   - - - - - - - - - - - 
      Exp0 = parPhi(1)*pi
      Exp1 = parPhi(2)*(1.-dexp(-parPhi(3)*pi))/parPhi(3)
      Exp2 = parPhi(4)*(1.-dexp(-parPhi(5)*pi))/parPhi(5)
      Exp3 = parPhi(6)*(1.-dexp(-parPhi(7)*pi))/parPhi(7)
      Summ = Exp0+Exp1+Exp2+Exp3
C     
      Exp0 = (Exp0+Exp1+Exp2+Exp3)/Summ
      Exp1 = (Exp1+Exp2+Exp3)/Summ   
      Exp2 = (Exp2+Exp3)/Summ 
      Exp3 =  Exp3/Summ
C-    write(*,*) 'Exp0,Exp1,Exp2,Exp3=',Exp0,Exp1,Exp2,Exp3
C
C  Pt initialisation:
C - - - - - - - - - - -        
C
      XmXmin = Xmin-Xmax
      XmXmax = Xmax-Xmin
C     
      XpXmin = 2.*Xmin
      XpXmax = 2.*Xmax
C     
CXmX-
      Exmx1 = parXmX(1)*(1.-dexp(-parXmX(2)*dabs(XmXmax)))/parXmX(2)
      Exmx2 = parXmX(3)*(1.-dexp(-parXmX(4)*dabs(XmXmax)))/parXmX(4)
      Exmx1 = Exmx1/(Exmx1+Exmx2)
      Exmx2 = 1. 
C     
C-    write(*,*) ' XpX:  XpXmin,XpXmax=',XpXmin,XpXmax
C-    write(*,*) ' XmX:  XmXmin,XmXmax=',XmXmin,XmXmax
C-    write(*,*) ' XmX:  Exmx1 ,Exmx2 =',Exmx1 ,Exmx2
C     
CXpX-
      Icase = 2                 !   Gausses and Exponent 
      Xmed  = 0.6
      if (XpXmax.lt.Xmed) then
         Icase = 1              !   Gauss only
         Xmed  = XpXmax
      endif		
C
      if (XpXmin.gt.Xmed) then
         Icase = 3              !   Exponent only
         Xmed  = XpXmin
      endif
C      
      if (Icase.eq.2)     then  !   Gauss and Exponent
         eps = 0.000001
         Gauss = Dgauss(DsdXpX ,XpXmin, Xmed  , eps)
         Expon = Dgauss(DsdXpX ,Xmed  , XpXmax, eps)
         Summ  = Gauss+Expon
         Gauss = Gauss/Summ 
         Expon = Expon/Summ
C-       write(*,*) 'XpX, Case 2: Gauss,Expon=',Gauss,Expon     
      endif                     ! Icase = 2
C
C
      if (Icase.eq.1.or.Icase.eq.2) then !   Gausses only and Icase=2
         
         Sgm1 =       1./dsqrt(2.*parXpX(2))
         
C-       write(*,*) 'XpX, Case 1:  Sgm1=',Sgm1        
      endif
C      
      if (Icase.eq.3.or.Icase.eq.2) then !    Exponent only and Icase=2   
C
         AnorX = 1.-dexp(-parXpX(5)*(XpXmax-Xmed))
C     
C-       write(*,*) 'XpX, Case 3:  AnorX=',AnorX
      endif
         
      return
      end
*
*=======================================================================
      subroutine ee_event(Ymin,Ymax,PTmin,PTmax,
     +                    Ye,Yp,Xe,Xp,Phi,Wtm2)
C------------------------------------------------------------------------------
C    Produce one event 
C
C    Input variables:  Ymin     - minimal value of rapidity  \
C                      Ymax     - maximal value of rapidity   | of kinematics 
C                      PTmin    - Pt minimum in MeV/c	      | range
C                      PTmax    - Pt maximum in MeV/c	     / 
C    Output variables: Ye,Yp    - rapidity of produced e- or e+
C                      Xe,Xp    - log10(pt) of produced e- or e+, pt in MeV/c
C                      Phi      - azymuth angle between e- and e+
C                      Wtm2     - event weight. The sum of these event 
C                                 weights, divided by the total number 
C                                 of generated events, gives the integral 
C                                 cross section of the process of e+e- pair 
C  		                  production in the above mentioned kinematics
C                                 range.
C		                  Sum of the selected event weights, divided 
C                                 by the total number of generated events, 
C                                 gives the integral cross section corresponded 
C                                 to the set of selected events
C------------------------------------------------------------------------------

      implicit real*8 (A-H,O-Z)
      common/parPh / parPhi(7)
      common/parYp / parYpY(6)
      common/parYm / parYmY(6)
      common/parXp / parXpX(5)
      common/parXm / parXmX(4)
      common/eevent/ Xsect2,Dsect2, Xsecttot,Dsecttot, Nevnt
      common/eepars/ Xmin,Xmax,YpYmin,YpYmax,WYpYmax,YmYmin,YmYmax,
     &               Ymed1,Ymed2,sgmY1,sgmY2,XYsect,
     &               Gaus1,Gaus2,Gauss,Exp1,Exp2,Exp3,
     &               XmXmin,XmXmax,XpXmin,XpXmax,Exmx1,Exmx2,Xmed,
     &               Sgm1,AnorX,Icase
      external DsdYpY,DsdYmY,DsdXpX,DsdXmX,DsdPhi 
      data pi     / 3.141 592 653 589 793 238 462 643d00 /

C  Rapidity distributions: 
C
 10   YpY = YpYmin + (YpYmax-YpYmin)*eernd(0)
      Wt  = DsdYpY(YpY)
      if (Wt.lt.eernd(0)*WYpYmax) go to 10
C
CYmY-
      r1 = eernd(0)
      if  (r1.lt.Gaus1) then
 13      call rnorml(rn1)
         YmY = sgmY1*rn1
         if (dabs(YmY).gt.Ymed1) go to 13
      else if  (r1.lt.Gaus2) then
 15      call rnorml(rn2)
         YmY = sgmY2*rn2
         if (dabs(YmY).lt.Ymed1) go to 15
         if (dabs(YmY).gt.Ymed2) go to 15
      else
 20      r2 = eernd(0)
         YmY=-Dlog(r2)/parYmY(6) + Ymed2
         if (YmY.gt.YmYmax) go to 20    
         r2 = eernd(0)
         if (r2.lt.0.5) YmY =-YmY
      endif
C
      Ye = 0.5*(YpY+YmY)
      Yp = 0.5*(YpY-YmY)
      
      if (Ye.lt.Ymin.or.Ye.gt.Ymax) go to 10
      if (Yp.lt.Ymin.or.Yp.gt.Ymax) go to 10      

C      
C Azimuthal angle:       
C 
 30   r1 = eernd(0)
      if (r1.lt.Exp3) then
 31      r2 = eernd(0)
         Phi=-Dlog(r2)/parPhi(7)
         if (Phi.gt.pi) go to 31
         go to 50            
      else if (r1.lt.Exp2) then
 32      r2 = eernd(0)
         Phi=-Dlog(r2)/parPhi(5)
         if (Phi.gt.pi) go to 32
         go to 50  
      else if (r1.lt.Exp1) then 
 33      r2 = eernd(0)
         Phi=-Dlog(r2)/parPhi(3)
         if (Phi.gt.pi) go to 33
         go to 50
      else
         Phi = pi*eernd(0)
      endif
C     
 50   if (eernd(0).gt.0.5) Phi =-Phi
      Phi = Phi+pi
      if (Phi.lt.0.03.or.Phi.gt.2.*pi-0.03) go to 30 
C
C Transverse momentums: 
C     
CXpX-

  60  Jcase = Icase             !       Jcase = Icase, if Icase.ne.2  
C
      if (Icase.eq.2) then      !       Gausses and Exponent       
         Jcase = 3
         if (eernd(0).lt.Gauss) Jcase = 1
      endif                     ! of Icase 2
C           
      if (Jcase.eq.1) then      !       Gauss only
 65      call rnorml(rn1)
         XpX=Sgm1*rn1-parXpX(3)
         if (XpX.lt.XpXmin.or.XpX.gt.Xmed) go to 65
      endif                     ! of Jcase 1	      
C
      if (Jcase.eq.3) then      !       Exponent only
 70      r1  = eernd(0)
         XpX =-Dlog( 1.- r1*AnorX)/parXpX(5)+Xmed 	      
      endif                     ! of Jcase 3 

CXmX-
      r1 = eernd(0)
      if (r1.lt.Exmx1) then
 81      r2 = eernd(0)
         XmX=-Dlog(r2)/parXmX(2)
         if (XmX.gt.dabs(XmXmax)) go to 81           
      else 
 82      r2 = eernd(0)
         XmX=-Dlog(r2)/parXmX(4)
         if (XmX.gt.dabs(XmXmax)) go to 82  
      endif
C
 85   if (eernd(0).gt.0.5) XmX =-XmX
C
      Xe = 0.5*(XpX+XmX)
      Xp = 0.5*(XpX-XmX)
      
      if (Xe.lt.Xmin.or.Xe.gt.Xmax) go to 60
      if (Xp.lt.Xmin.or.Xp.gt.Xmax) go to 60 
C
C --- Good event: 
C
      Pte = 10.d0**Xe
      Ptp = 10.d0**Xp
C
C Event weight if the events would generate uniformly
C      
      Wt = DsdYpY(YpY)*DsdYmY(YmY)*DsdXpX(Xpx)*DsdXmX(XmX)*DsdPhi(Phi)
      Wt = 2.2706950/Wt
C
      Wtm2=(Pte*Ptp)*Wt          ! Transformation factor
C
C --- Exact diff. cross section (Adrian Alscher 1997, Kai Hencken, May '98):  
C
      call Diffcross(Ptp,Yp,Pte,Ye,Phi,Dsigma)
      Wtm2= XYsect*Dsigma*(Pte*Ptp)*Wtm2
C   
      Nevnt  = Nevnt  + 1
      Xsect2 = Xsect2 + Wtm2
      Dsect2 = Dsect2 + Wtm2*Wtm2
C
C     Calculate cross section and its error accumulated so far
C
      Xsecttot = Xsect2/Nevnt
      Dsecttot = dsqrt(Dsect2)/Nevnt
C
      return
      end
c      
C=================================================================
      Double precision function DsdYpY(Y)    
c
c     Y = Yp+Ye
c - - - - - - - - - - - - - - - - - - -
      implicit real*8 (A-H,O-Z)
      common/parYp/   parYpY(6)
C-    data parYpY /  -48.584, 0.11403E-01, 40.649, 0.38920E-04,
C-   +                                     40.283, 0.19238E-01 / 
      data parYpY /  -4.8584, 0.11403E-01, 4.0649, 0.38920E-04,
     +                                     4.0283, 0.19238E-01 / 
c       
      DsdYpY = parYpY(1)*dexp(-parYpY(2)*Y**2)
     +       + parYpY(3)*dexp(-parYpY(4)*Y**4)
     +       + parYpY(5)*dexp(-parYpY(6)*Y**2)
      return
      end
C 
      Double precision function DsdYmY(Y)    
c
c     Y = Yp-Ye
c - - - - - - - - - - - - - - - - - - -
      implicit real*8 (A-H,O-Z)
      common/parYm/   parYmY(6)
C-    data parYmY /   180.,  8., 174.38, 0.17867, 2418.8, 1.3097 /
      data parYmY /   1.80,  8., 1.7438, 0.17867, 24.188, 1.3097 /
C      
      if (abs(Y).lt.0.18) then
C1-      if (abs(Y).lt.0.00) then
         DsdYmY  = parYmY(1)*dexp(-parYmY(2)*abs(Y)**2)      
      else if (abs(Y).lt.4.00) then
         DsdYmY  = parYmY(3)*dexp(-parYmY(4)*abs(Y)**2)
      else
         DsdYmY  = parYmY(5)*dexp(-parYmY(6)*abs(Y))
      endif
         
      return
      end
C      
      Double precision function DsdYY(Yp,Ye)    
c - - - - - - - - - - - - - - - - - - - - - - 
      implicit real*8 (A-H,O-Z)
c
      YpY   = Yp+Ye
      YmY   = Yp-Ye
      DsdYY = DsdYpY(YpY)*DsdYmY(YmY)
      return
      end
C 
C====================================================================
C 
      Double precision function DsdXpX(X)    
c
c     X = Xp+Xe
c - - - - - - - - - - - - - - - - - - -
      implicit real*8 (A-H,O-Z)
      common/parXp/   parXpX(5)

C-    data parXpX /   83.668, 1.2004, 0.47225, 96.951, 2.3814 /
      data parXpX /   8.3668, 1.2004, 0.47225, 9.6951, 2.3814 /
c       
      if (X.lt.0.6)then 
         DsdXpX = parXpX(1)*dexp(-parXpX(2)*(X+parXpX(3))**2)         
      else
         DsdXpX = parXpX(4)*dexp(-parXpX(5)*X)
      endif
      return
      end
C 
      Double precision function DsdXmX(X)    
c
c     X = Xp-Xe
c - - - - - - - - - - - - - - - - - - -
      implicit real*8 (A-H,O-Z)
      common/parXm/   parXmX(4)
c-    data parXmX /   950.38, 39.040, 194.92, 3.5660 /
      data parXmX /   9.5038, 39.040, 1.9492, 3.5660 /
c      
 
      DsdXmX = parXmX(1)*dexp(-parXmX(2)*abs(X))   
     +       + parXmX(3)*dexp(-parXmX(4)*abs(X))
         
      return
      end        
C      
      Double precision function DsdXX(Xp,Xe)    
c - - - - - - - - - - - - - - - - - - - - - - 
      implicit real*8 (A-H,O-Z)
c
      XpX   = Xp+Xe
      XmX   = Xp-Xe
      DsdXX = DsdXpX(XpX)*DsdXmX(XmX)
      return
      end
C====================================================================
C
      Double precision function DsdPhi(Phi)  ! Phi distribution 
c
c     Phi - athimutal angle between e+ and e- (in radians)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      implicit real*8 (A-H,O-Z)
      common/parPh/ parPhi(7)
C                                                         !  For 10**7 events
      data parPhi / 5.4825, 226.00, 11.602,      68.173,  !  accuracy= 0.87 % 
     +			    1.9509, 0.11864E+04, 109.61 / !  N(Wt>20)= 573  
 
C
      data pi     / 3.141 592 653 589 793 238 462 643d00/
C
      DsdPhi  = parPhi(1)
     +        + parPhi(2)*dexp(-parPhi(3)*dabs(Phi-pi)) 
     +        + parPhi(4)*dexp(-parPhi(5)*dabs(Phi-pi))
     +        + parPhi(6)*dexp(-parPhi(7)*dabs(Phi-pi))            
      return
      end
C====================================================================
      SUBROUTINE rnorml(rnd)
*     Random generator of normal distribution
      implicit real*8 (A-H,O-Z)
      PARAMETER (pi=3.141 592 653 589 793 238 462 643d00)
      PARAMETER (pi2=2.*pi)
C      
      u1 = eernd(0)
      u2 = eernd(0)
      p1 = pi2*u1
      p2 = dsqrt(-2.*dlog(u2))
      rnd= dcos(p1)*p2
      RETURN
      END
