AMPT Users' Guide

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
8/2009 test version v1.25t3/v2.25t3:
     * enable embedding of a given number of soft pions along each parton 
       of the embedded back-to-back jet pair by adding three new 
       input values in file `input.ampt': nsembd, psembd, tmaxembd. 
       nsembd: number of soft pions to be embedded with each high-Pt parton
        in the embedded jet pair;
       psembd: the momentum of each embedded soft pion in GeV/c;
       tmaxembd: maximum angle (in radian) of embedded soft pions relative to 
        the associated high-Pt parton, then the soft pions are generated 
        uniformly in relative polar and azimutal angle in this cone around the 
        high-Pt parton. Note that soft pions are not generated uniformly 
        in solid angle in this cone because that gives a valley at theta=0, 
        unlike the primitive jet-like correlation with a peak at theta=0.
     * enable the embedded jet pair to have random azimuthal angle (still 
        back-to-back) and to take positions according to a user file 
        "embed-jet-xy.txt"; enabled by additional values of iembed:
        iembed= 
         1: on with fixed position(xembd,pembd) and Pt(pxqembd,pyqembd);
         2: on with fixed position(xembd,pembd) and random azimuthal angle
          with Pt-magnitude given by sqrt(pxqembd^2+pyqembd^2); 
         3: on with random position and fixed Pt(pxqembd,pyqembd);
         4: on with random position and random random azimuthal angle
          with Pt-magnitude given by sqrt(pxqembd^2+pyqembd^2); 
          for iembed=3 or 4: need a position file "embed-jet-xy.txt";
         Other integers: turn off embedding.
   The above modification can be found by searching "clin-8/2009". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
7/2009 test version v1.25t2/v2.25t2:
     * added two new input values in file `input.ampt': ishadow and dshadow.  
       Setting ishadow=1 (when shadowing flag is on) enables users to set 
       dshadow (between 0.d0 & 1.d0) to get an intermediate nuclear shadowing 
       between no-shadowing (when dshadow=0.d0) 
       and the default HIJING shadowing (when dshadow=1.d0).
   The above modification can be found by searching "clin-7/2009". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
7/2009 test version v1.25/v2.25:
	corrected the explanations for deuteron cross section assumptions in
	input.ampt and art1f.f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
6/2009 test version v1.25/v2.25 (==v1.25t1/v2.25t1 meaning "test version 1"):
     * added two new input values in file `input.ampt': pttrig and maxmiss.  
       When enabled, event selection takes place so that each event will have 
       at least 1 minijet parton with Pt above pttrig in the initial condition 
       (before parton & hadron cascades).
       For this purpose, the HIJING initial condition is repeated until 
       such an event is generated (but no more than maxmiss times); 
     * added two new input values in file `input.ampt' so that one can turn off
       initial and/or final state radiation as well as the Pt kick 
       in the HIJING event initial condition; with these effects turned off, 
       the selected events using pttrig are more likely to have 
       back-to-back initial minijets;
     * added five new input values in file `input.ampt': iembed, pxqembd, 
       pyqembd, xembd, and yembd. Setting iembed=1 embeds one back-to-back 
       q/qbar jet pair per event at the specified transverse momentum and 
       transverse position before parton cascade starts;
     * changed the format of zpc.dat, including the following: 
       event number is given as the first value, and an event-iteration 
       flag is added as the 2nd value in 1st line of each event in zpc.dat;
       an event-iteration flag other than 0 means the event is being repeated 
       (detailed message is given in the nohup.out file);
     * ana/initial_parton_sm.dat is renamed as 
       ana/parton-initial-afterPropagation.dat, 
       which outputs the complete information of partons that enter the 
       parton cascade; i.e., it gives the minijet gluon information 
       for default runs, or the quark and anti-quark information 
       after string melting (before parton cascade) for string melting runs. 
       This option is activated by setting ioscar to 2.
       It has the same format as ana/ampt.dat except for the 1st line 
       for each event, which gives the following:
       For isoft=1 (default AMPT):
         event number, event-iteration flag, number of partons;
       For isoft=4 (String Melting):
         event number, event-iteration flag, number of partons, 
         number of formed baryons, number of formed mesons, 
         total number of initial particles, 
         total number of initial particles that cannot enter ZPC.
     * added ana/npart-xy.dat, which gives the transverse positions
       of all initial nucleons and their status; format is:
       For each event, the first line gives: 
         event #, event-iteration flag, atomic masses of projectile and target.
       Each of the following four lines gives:
         x,y, sequence number in nucleus (positive for projectile, negative 
         for target), status (0: spectator, 1 or 2: wounded due to 
         elastic collisions, 3: wounded due to inelastic collisions).
     * ioscar=3 now enables the following output files in addition to 
       ana/parton-initial-afterPropagation.dat:
       all parton collision history in ana/parton-collisionsHistory.dat, 
       minijet initial condition in ana/minijet-initial-beforePropagation.dat,
       ana/parton-after-coalescence.dat for String Melting.
     * format of ana/parton-collisionsHistory.dat:
       For each collision, the first line gives: 
         event #, event-iteration flag, sequence numbers of the two partons;
       Each of the following four lines gives:
         PYTHIA particle ID number, three-momentum(Px,Py,Pz), mass, and 
         space-time coordinates(x,y,z,t)
         for parton1&2 before the collision parton1&2 after the collision.
     * format of ana/minijet-initial-beforePropagation.dat:
       For each event, the first line gives: 
         event #, event-iteration flag, atomic masses of projectile and target.
       Each of the following lines gives:
         PYTHIA particle ID number, three-momentum(Px,Py,Pz), mass, and 
         space-time coordinates(x,y,z,t) of one minijet parton at production,
         ID of origin (from 1: projectile, 2: target; 3: independent string).
     * format of ana/parton-after-coalescence.dat:
       For each event, the first line gives: 
         event number, number of partons, number of formed baryons, 
         number of formed mesons, impact-parameter, 
         number of participant nucleons in projectile due to elastic 
         collisions, number of participant nucleons in projectile due to 
         inelastic collisions, and corresponding numbers in target. 
       Each of the following lines gives:
         PYTHIA parton ID number, three-momentum(Px,Py,Pz), mass, 
         sequence number of formed hadron, PYTHIA ID number of formed hadron;
     * `common /para7/ ioscar' is changed to include string melting info;
     * info under "in HJANA1" in file nohup.out is corrected.

   The above modifications can be found by searching "clin-6/2009". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
5/2009 test version v1.24/v2.24:
*  Changes made:
     * corrected a bug on freezeout time (values in `ampt.dat') of particles 
       whose formation time is larger than the hadron cascade termination time;
     * freezeout time of spectator projectile and target nucleons should be ~0
       but it was not correctly updated for string melting runs; now corrected.

   The above modifications can be found by searching "clin-5/2009". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
4/2009 update
     corrected the explanation on the structure of zpc.dat.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3/2009 test version v1.23/v2.23:
*  Changes made:

     * included a subroutine addhad() in amptsub.f to insert 
       user-defined hadrons before the start of the hadron cascade, 
       user must set IDPERT to either 1 or 2, user can set NPERTD to 0 
       if perturbative deuteron productions are not needed;
     * add a comment in input.dat saying that setting NPERTD to 0 
       turns off perturbative deuteron productions;
     * moved the location of the subroutine hbtout() in art1f.f, 
       and as a result there are small changes in hadron freezeout 
       space-time values, e.g., the freezeout space values may change 
       by v*DT with v being the particle momentum, and the freezeout time 
       of spectator projectile or target nucleons change by DT;
     * increased the output accuracy of momentum Pz by 1 more digit
       in ampt.dat and ampt_pert.dat.

   The above modifications can be found by searching "clin-3/2009". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
10/2008 test version v1.22/v2.22:
*  Changes made in the default AMPT model (version 1.22) 
   and the AMPT model with string melting (version 2.22):

     We have included deuteron(d) production and destruction
     in the hadronic cascade via d+M <-> B+B, where M represents a meson 
     including pi, rho, omega and eta; and B represents a baryon 
     including proton, neutron, Delta, N*(1440) and N*(1535), 
     Anti-deuteron processes are included similarly. 
     Elastic collisions of d+M and d+B are also included. 
     The cross sections for d+pi <-> N+N are based on experimental data; 
     for other reactions, 4 different assumptions on their cross sections 
     are available by choosing the value of idxsec in the file `input.ampt'.
     Perturbative production/destruction of deuterons can be
     turned on by using idpert=1 (or 2) and 1<=NPERTD<=10000 in input.ampt,
     and perturbative deuteron results are stored in ana/ampt_pert.dat, 
     which has the same data format as ana/ampt.dat
     (except that ampt_pert.dat does not have the column for deuteron mass
     but has the perturbative probability in the last column).
     Regular (i.e. non-perturbative) results can be found in `ana/ampt.dat'.
     Note that you can use either the regular or perturbative results for
     deuterons; you can not use both (that would be double-counting);
     regular or perturbative results are close but not exactly the same
     due to the biase from triggering perturbative productions.
     A deuteron (anti-deuteron) has the particle ID of 42 (-42) in
     ana/ampt.dat and ana/ampt_pert.dat.

   The above modifications can be found by searching "clin-5/2008", 
     "clin-6/2008", "clin-8/2008" and "clin-9/2008".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
10/01/2008 v1.21/v2.21:
*  Changes made in the default AMPT model (version 1.21) 
   and the AMPT model with string melting (version 2.21):

     We have added an option to write out the complete parton information 
       right after string melting (before parton cascade). 
       Results are stored in ana/initial_parton_sm.dat in the same format 
       as ana/ampt.dat except for the first few lines that provide 
       information for each event. This option is only available for 
       string melting and can be activated by setting ioscar to 2.

     We have added an option to turn off phi meson decays at the end 
       of hadron cascade, i.e., at NT=NTMAX. 
       This option can be activated by setting iphidcy to 0.
       Note that phi decays during hadron cascade are always enabled.

   The above two modifications can be found by searching "clin-5b/2008". 
	
   Note: the following physics extensions are preliminary and under test, 
     therefore they have been disabled in v1.21/v2.21:
     We have included deuteron(d) production and destruction
     in the hadronic cascade via d+M <-> B+B, where M represents a meson 
     including pi, rho, omega and eta; and B represents a baryon 
     including proton, neutron, Delta, N*(1440) and N*(1535), 
     Anti-deuteron processes are included similarly. 
     Elastic collisions of d+M and d+B are also included. 

   The above modifications can be found by searching "clin-5/2008", 
     "clin-6/2008", "clin-8/2008" and "clin-9/2008".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
04/29/2008 v1.13/v2.13:
*  Changes made in the default AMPT model (version 1.13) 
   and the AMPT model with string melting (version 2.13):

     multiple array arguments are checked against out-of-bound in amptsub.f;
     commented out the unused "CALL HJAN1A" in zpc.f;
     a pause statement in ran1() in zpc.f is modified;
     array sizes of ekaon() and sekaon() are increased from 200 to 2000;
     array sizes in /HJJET2/ and /xydr/ are increased from 900 to 150001;
     multiple compound IF statements are broken up in hipyset1.35.f;
     DARWIN added to the list of operating system in `exec';
     added a check on the range of IX, IY, IZ in art1f.f and modified 
       other such checks;
     added a check on the range of npion in art1f.f;
     RAN() renamed to RANART() to avoid conflict with system functions;
     bugs on initializations of xlast() and plast() are fixed;
     the variable ISS is modified to avoid out-of-bound error in EKAON();
     "IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR)" is modified 
       to avoid out-of-bound error in K();
     "DATA NSEED/74769375/" in hijing1.383_ampt.f is commented out;
     PYWIDT() subroutine is modified according to pythia-6115.f to avoid 
       undefined values for variables GGF,GZF,GZPF,ZZF,ZZPF,ZPZPF,API;
     "MDCY(KFPR(ISUB,1),1)" is changed to "MDCY(LUCOMP(KFPR(ISUB,1)),1)"
       to avoid invalid values for the 1st argument of MDCY();
     "if (jscat .ne. 0 .and. next(jscat) .ne. iscat)" is modified 
       to avoid out-of-bound error in next().

   The above modifications can be found by searching "clin-4/2008". 
   They are not found to change ampt.dat from a few tests on a Linux OS.
   We thanks A. Vander Molen and G. Westfall for pointing out these issues. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
2005 v1.12/v2.12:
*  Changes made in the default AMPT model (version 1.12) 
   and the AMPT model with string melting (version 2.12):

     Freezeout time of spectator projectile and target nucleons should be ~0
       but it was not correctly updated in `ampt.dat'; now corrected.

   The above modifications can be found by searching "clin-12/14/03". 

   We have corrected a typo in `input.ampt':
     "IZT (target A number)" is changed to "IZT (target Z number)".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
11/30/2004 v1.11/v2.11:
*  The default AMPT model (version 1.11) and the AMPT model with
   string melting (version 2.11) both use an initialization file 
   `input.ampt'. The analysis directory `ana/' contains the resulting
   data files. The final particle record file is `ana/ampt.dat'.
   The version number of AMPT is written to both the file `ana/version'
   and `nohup.out'.  The AMPT source code has been tested for both the
   f77 and the pgf77 compilers on the Unix, Linux, and OSF1 operating systems.

*  To run the AMPT program, one needs to:
   1. set the initial parameters in `input.ampt'. If one prefers to
      use run-time random number seed, set `ihjsed=11', In this way, every
      run is different even with the same `input.ampt' file.
   2. type `sh exec &' to compile and run the executable `ampt'
      with some general information written in `nohup.out'.

*  Key initial parameters in `input.ampt' are:
   EFRM: sqrt(s_NN) in GeV, e.g. 200 for the maximum RHIC energy.
   NEVNT: the total number of events.
   BMIN, BMAX: the minimum and maximum impact parameter (in fm) 
      for all events with BMAX having an upper limit of HIPR1(34)+HIPR1(35)
      (=19.87 fm for d+Au collisions and 25.60 fm for Au+Au collisions). 
   ISOFT: choice of parton-hadron conversion scenario.
      =1: default AMPT model (version 1.x);
      =4: the AMPT model with string melting (version 2.y).
         Note that values of 2, 3, and 5 have never been used for 
         publications. They are tests of other string melting scenarios:
         =2: a string is decomposed into q+dq+minijet partons instead of 
            using the Lund fragmentation;
         =3: a baryon is decomposed into q+qq instead of 3 quarks;
         =5: same as 4 but partons freeze out according to
            local energy density.
   NTMAX: the number of time-steps for hadron cascade, default(D)=150.
      Note that NTMAX=3 effectively turns off hadron cascade, 
      and a larger value than default is usually necessary 
      for observables at large rapidity or large pseudorapidity.
      We use NTMAX=1000 for HBT studies in central Au+Au
      collisions due to the need for the last interaction points 
      and for LHC calculations due to the longer lifetime of the formed matter.
   DT: value of the time-step (in fm/c) for hadron cascade, D=0.2.
      Note that t_cut=NTMAX*DT is the termination time of hadron cascade. 
   PARJ(41): parameter a in the Lund symmetric splitting function. 
   PARJ(42): parameter b in the Lund symmetric splitting function 
      (in GeV**(-2)). Note that we use default value in HIJING 
      (a=0.5 and b=0.9) for d+Au collisions, 
      and a=2.2 and b=0.5 for collisions of heavy nuclei.
   flag for popcorn mechanism: D=1(Yes) turns on the popcorn mechanism. 
      In general, it increases baryon stopping.
   PARJ(5): controls BMBbar vs. BBbar in the popcorn mechanism, D=1.0. 
   shadowing flag: D=1(Yes) turns on nuclear shadowing. 
   quenching flag: D=0(No) turns off jet quenching 
      since the parton cascade ZPC simulates final-state effects. 
   p0 cutoff: D=2.0 (in GeV/c) for p0 in HIJING for minijet production. 
   parton screening mass: controls the parton cross section, 
      D=3.2264 (in fm**(-1)). Its square is inversely proportional to 
      the parton cross section. Use D=3.2264d0 for 3mb, and 2.2814d0 for 6mb.
   ihjsed: choice of the random number seed, D=0.
      =0: take the `Ran Seed for HIJING' in `input.ampt'
         and disregard the random value generated in the file `exec'.
      =11: take the HIJING random seed at runtime from the file `exec', 
         with the seed written in `nohup.out' and `ana/version'.
   Ran Seed for HIJING: random number seed for HIJING when ihjsed=0.
   Kshort decay flag: depends on the experimental correction procedure, 
      D=0 turns off Kshort decays after the hadron cascade.
      Note that decays of the following resonances and their
      antiparticles are always included: 
      rho, omega, eta, K*, phi, Delta, N*(1440), N*(1535),
      Sigma0 (in order to include its feed down to Lambda). 
   optional OSCAR output: if set to 1, outputs in OSCAR1997A format
      are written in `ana/parton.oscar' and `ana/hadron.oscar'. 
   dpcoal: parton coalescence distance in momentum space (in GeV/c).
   drcoal: parton coalescence distance in coordinate space (in fm).
      dpcoal, drcoal both have D=10**6 for nearest-neighbor coalescence 
      in the AMPT model with string melting. 

*  Key output file are:
   ana/ampt.dat: It contains particle records at hadron kinetic freeze-out, 
      i.e., at the last interaction point. 
      For each event, the first line gives: 
         event number, test number(=1), number of particles in the event, 
         impact parameter, total number of participant nucleons in projectile,
         total number of participant nucleons in target, number of participant 
         nucleons in projectile due to elastic collisions, number of 
         participant nucleons in projectile due to inelastic collisions, 
         and corresponding numbers in target. 
         Note that participant nucleon numbers include nucleons participating 
         in both elastic and inelastic collisions.
      Each of the following lines gives: 
         PYTHIA particle ID number, three-momentum(Px,Py,Pz), mass, and 
         space-time coordinates(x,y,z,t) of one final particle at freeze-out.

   ana/zpc.dat:    similar to `ana/ampt.dat' but for partons at freeze-out.
      The first line of each event gives:
         event number, (event-iteration flag, added in 6/2009),
         number of partons in the event, impact-parameter,
         number of participant nucleons in projectile due to elastic 
         collisions, number of participant nucleons in projectile due to 
         inelastic collisions, and corresponding numbers in target. 
      Each of the following lines gives:
        For isoft=1 (default AMPT):
          PYTHIA particle ID number, three-momentum(Px,Py,Pz), mass, and 
          space-time coordinates(x,y,z,t) of one final parton at freeze-out.
        For isoft=4 (String Melting, format changed in 6/2009):
          PYTHIA particle ID number, three-momentum(Px,Py,Pz), mass, 
          the parent hadron sequence number that the parton comes from,
          the parton sequence number of this parton in the parent hadron
          (1-2 for a meson, 1-3 for a baryon), and freeze-out time (t).

   Note that momenta are in units of GeV/c, mass in GeV/c**2, 
      space in fm, and time in fm/c. 
      If a particle comes from the decay of a resonance which still exists 
      at the termination time of hadron cascade, then its space-time 
      corresponds to the decay point of the parent resonance.
      Also note that the x-axis in AMPT is defined as the direction along 
      the impact parameter, and the z-axis is defined as the beam direction. 


Please do not hesitate to contact us if needed. Have fun!
 
Zi-Wei Lin (linz@ecu.edu)
Che-Ming Ko (ko@comp.tamu.edu)
Bao-An Li (Bao-An_Li@tamu-commerce.edu)
Subrata Pal (spal@tifr.res.in)
Bin Zhang (bzhang@astate.edu)
 
6/25/2009
