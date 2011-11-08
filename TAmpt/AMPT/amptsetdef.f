      SUBROUTINE AMPTSETDEF
c
cgsfs added following line to match C++ call
      double precision xmp, xmu, alpha, rscut2, cutof2
      double precision smearp,smearh,dpcoal,drcoal,ecritl
      character*25 amptvn
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
      COMMON/LUDAT1A/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /AROUT/ IOUT
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
c     parton coalescence radii in case of string melting:
      common /coal/dpcoal,drcoal,ecritl
      common/snn/efrm,npart1,npart2
c     initialization value for parton cascade:
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common /para8/ idpert,npertd,idxsec
      common /rndm3/ iseedp
c     initialization value for hadron cascade:
      COMMON /RUN/ NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy,ikstardcy
clin-6/2009:
c      common/phidcy/iphidcy
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss
      common/embed/iembed,pxqembd,pyqembd,xembd,yembd
      common/popcorn/ipop

      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      SAVE   
c****************

c     flag to select default AMPT or string melting:
      isoft    = 4		! ISOFT (D=1): select Default AMPT or String Melting(see below)
c     read initialization value for hadron cascade:
      NTMAX    = 150		! NTMAX: number of timesteps (D=150), see below
      DT       = 0.2		! DT: timestep in fm (hadron cascade time= DT*NTMAX) (D=0.2)
c     parj(41) and (42) are a and b parameters in Lund string fragmentation:
      PARJ(41) = 0.5		! PARJ(41): parameter a in Lund symmetric splitting function
      PARJ(42) = 0.9      	! PARJ(42): parameter b in Lund symmetric splitting function
c     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and 
c     increase the net-baryon stopping in rapidity (value HIJING is 1):
      ipop      = 1	      	! (D=1,yes;0,no) flag for popcorn mechanism(netbaryon stopping)
c     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn:
      PARJ(5)   = 1.0	      	! PARJ(5) to control BMBbar vs BBbar in popcorn (D=1.0)
c     shadowing flag in HIJING:
      IHPR2(6)  = 1		! shadowing flag (Default=1,yes; 0,no)
c     quenching flag in HIJING:
      IHPR2(4)  = 0		! quenching flag (D=0,no; 1,yes)
c     quenching rate when quenching flag is on (=1.0 GeV/fm):
      HIPR1(14) = 2.0		! quenching parameter -dE/dx (GeV/fm) in case quenching flag=1
c     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV. 
      HIPR1(8) = 2.0		! p0 cutoff in HIJING for minijet productions (D=2.0)
c     read initialization value for parton cascade:
      xmu      = 3.2264d0	! parton screening mass in fm^(-1) (D=3.2264d0), see below
      izpc     = 0		! IZPC: (D=0 forward-angle parton scatterings; 100,isotropic)
      alpha    = 0.3333d0	! alpha in parton cascade
c     quark coalescence radii in momentum and space for string melting:
      dpcoal   = 1d6		! dpcoal in GeV
      drcoal   = 1d6		! drcoal in fm
c     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0):
      ihjsed   = 0		! ihjsed: take HIJING seed from below (D=0)or at runtime(11)
c     2 seeds for random number generators in HIJING/hadron cascade and ZPC:
      nseed    = 53153511	! random seed for HIJING
      iseedp   = 8		! random seed for parton cascade
      iksdcy   = 0		! flag for Ks0 weak decays (D=0,no; 1,yes)
      ikstardcy   = 0		! flag for K* weak decays (D=0,no; 1,yes)
      iphidcy  = 0		! flag for phi decays at end of hadron cascade (D=1,yes; 0,no)
c     flag for OSCAR output for final partons and hadrons:
      ioscar   = 0		! optional OSCAR output (D=0,no; 1,yes; 2&3,more parton info)
clin-5/2008     flag for perturbative treatment of deuterons:
      idpert   = 0		! flag for perturbative deuteron calculation (D=0,no; 1or2,yes)
      npertd   = 1		! integer factor for perturbative deuterons(>=1 & <=10000)
      idxsec   = 1		! choice of cross section assumptions for deuteron reactions
clin-6/2009 To select events that have at least 1 high-Pt minijet parton:
      pttrig   = -7.		! Pt in GeV: generate events with >=1 minijet above this value
      maxmiss  = 1000		! maxmiss (D=1000): maximum # of tries to repeat a HIJING event
      IHPR2(2) = 3		! flag to turn off initial and final state radiation (D=3)
      IHPR2(5) = 1		! flag to turn off Kt kick (D=1)
clin-6/2009 To embed a back-to-back q/qbar pair into each event:
      iembed   = 0		! flag to turn on quark pair embedding (D=0,no; 1to4:yes)
      pxqembd  = 7.             ! Initial Px and Py values (GeV) of the embedded quark (u or d)
      pyqembd  = 0.
      xembd    = 0.             ! Initial x & y values (fm) of the embedded back-to-back q/qbar
      yembd    = 0.
      nsembd   = 1              ! nsembd(D=0), psembd (in GeV),tmaxembd (in radian).
      psembd   = 5.
      tmaxembd = 0.
clin-7/2009 Allow modification of nuclear shadowing:
      ishadow  = 0 		! Flag to enable users to modify shadowing (D=0,no; 1,yes)
      dshadow  = 1.d0		! Factor used to modify nuclear shadowing
c
      RETURN
      END

c$$$%%%%%%%%%% Further explanations:
c$$$ISOFT:  1 Default, 
c$$$        4 String Melting.
c$$$PARJ(41) & (42): 2.2 & 0.5/GeV^2 used for heavy ion (Au+Au, Pb+Pb) collisions,
c$$$        while the HIJING values (0.5 & 0.9/GeV^2) describe well 
c$$$        Nch in pp collisions and are used for d-Au collisions.
c$$$NTMAX:	number of time-steps for hadron cascade. 
c$$$	Use a large value (e.g. 1000) for HBT studies in heavy ion collisions.
c$$$	Using NTMAX=3 effectively turns off hadronic cascade.
c$$$parton screening mass (in 1/fm): its square is inversely proportional to 
c$$$	the parton cross section. Use D=3.2264d0 for 3mb cross section, 
c$$$	and 2.2814d0 for 6mb. Using 1d4 effectively turns off parton cascade.
c$$$ihjsed: if =11, take HIJING random seed at runtime so that 
c$$$	every run may be automatically different (see file 'exec').
c$$$iksdcy: flag for Ks0 weak decays for comparison with data.
c$$$ikstardcy: flag for K* weak decays for comparison with data.
c$$$iphidcy: flag for phi meson decays at the end of hadron cascade for comparison 
c$$$	with data; default is yes; use 0 to turn off these decays. 
c$$$	Note: phi meson decay during hadron cascade is always enabled.
c$$$ioscar:	0 Dafault,
c$$$	1 Write output in the OSCAR format,
c$$$	2 Write out the complete parton information 
c$$$		(ana/parton-initial-afterPropagation.dat)
c$$$        	right after string melting (before parton cascade),
c$$$	3 Write out several more files on parton information (see readme).
c$$$idpert:	flag for perturbative deuteron and antideuteron calculations 
c$$$	with results in ana/ampt_pert.dat:
c$$$	0 No perturbative calculations,
c$$$	1 Trigger a production of NPERTD perturbative deuterons 
c$$$		in each NN collision,	
c$$$	2 Trigger a production of NPERTD perturbative deuterons only in 
c$$$		an NN collision where a conventional deuteron is produced.
c$$$	Note: conventional deuteron calculations are always performed
c$$$		with results in ana/ampt.dat.
c$$$NPERTD:	number of perturbative deuterons produced in each triggered collision;
c$$$	setting it to 0 turns off perturbative deuteron productions.
c$$$idxsec: choose a cross section model for deuteron inelastic/elastic collisions:
c$$$	1: same |matrix element|**2/s (after averaging over initial spins 
c$$$		and isospins) for B+B -> deuteron+meson at the same sqrt(s);
c$$$	2: same |matrix element|**2/s for B+B -> deuteron+meson 
c$$$		at the same sqrt(s)-threshold;
c$$$	3: same |matrix element|**2/s for deuteron+meson -> B+B 
c$$$		at the same sqrt(s);
c$$$ 	4: same |matrix element|**2/s for deuteron+meson -> B+B 
c$$$		at the same sqrt(s)-threshold;
c$$$	1 or 3 also chooses the same cross section for deuteron+meson or baryon
c$$$		elastic collision at the same sqrt(s);
c$$$	2 or 4 also chooses the same cross section for deuteron+meson or baryon
c$$$		elastic collision at the same sqrt(s)-threshold.
c$$$%%%%%%%%%% For jet studies:
c$$$pttrig:	generate events with at least 1 initial minijet parton above this Pt 
c$$$	value, otherwise repeat HIJING event until reaching maxmiss tries;
c$$$	use a negative value to disable this requirement and get normal events.
c$$$maxmiss: maximum number of tries for the repetition of a HIJING event to obtain
c$$$	a minijet above the Pt value of pttrig;	increase maxmiss if some events
c$$$	fail to generate at least 1 initial minijet parton above pttrig. 
c$$$	it is safer to set a large value for high pttrig and/or large b value
c$$$	and/or smaller colliding nuclei.
c$$$IHPR2(2): flag to turn off initial and final state radiation: 
c$$$	0 both radiation off, 1 only final off, 2 only initial off, 3 both on.
c$$$IHPR2(5): flag to turn off Pt kick due to soft interactions: 0 off, 1 on.
c$$$	Setting both IHPR2(2) and IHPR2(5) to zero makes it more likely to 
c$$$	have two high-Pt minijet partons that are close to back-to-back.
c$$$%%%%%%%%%% To embed a back-to-back light q/qbar jet pair 
c$$$%%%%%%%%%%  and a given number of soft pions along each jet into each event:
c$$$iembed: flag to turn on quark pair embedding: 
c$$$        1: on with fixed position(xembd,pembd) and Pt(pxqembd,pyqembd);
c$$$        2: on with fixed position(xembd,pembd) and random azimuthal angle
c$$$         with Pt-magnitude given by sqrt(pxqembd^2+pyqembd^2); 
c$$$        3: on with random position and fixed Pt(pxqembd,pyqembd);
c$$$        4: on with random position and random random azimuthal angle
c$$$         with Pt-magnitude given by sqrt(pxqembd^2+pyqembd^2); 
c$$$	 for iembed=3 or 4: need a position file "embed-jet-xy.txt";
c$$$	Other integers: off.
c$$$pxqembd, pyqembd: sqrt(pxqembd^2+pyqembd^2) > 70MeV/c is required;
c$$$	the embedded quark and antiquark have pz=0.
c$$$xembd, yembd: the embedded quark and antiquark jets have z=0 initially. Note: 
c$$$	the x-axis is defined as the direction along the impact parameter.
c$$$nsembd:	number of soft pions to be embedded with each high-Pt parton
c$$$	in the embedded jet pair.
c$$$psembd: Momentum of each embedded soft pion in GeV.
c$$$tmaxembd: maximum angle(rad) of embedded soft pions relative to high-Pt parton.
c$$$%%%%%%%%%% User modification of nuclear shadowing:
c$$$ishadow: set to 1 to enable users to adjust nuclear shadowing
c$$$	provided the shadowing flag IHPR2(6) is turned on; default value is 0. 
c$$$dshadow: valid when ishadow=1; this parameter modifies the HIJING shadowing
c$$$	parameterization Ra(x,r)==1+fa(x,r) via Ra(x,r)==1+fa(x,r)*dshadow,  
c$$$	so the value of 0.d0 turns off shadowing 
c$$$	and the value of 1.d0 uses the default HIJING shadowing;  
c$$$	currently limited to 0.d0<=dshadow<=1.d0 to make sure Ra(x,r)>0.
