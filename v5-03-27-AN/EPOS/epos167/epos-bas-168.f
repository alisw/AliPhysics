
c-----------------------------------------------------------------------
      subroutine aaset(iop)
c-----------------------------------------------------------------------
c     sets parameters and options, initializes counters ...
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incpar'
      include 'epos.incsem'
      common/record/maxrec(2),irecty(30,2)
      common/cfacmss/facmss /cr3pomi/r3pomi,r4pomi/cifset/ifset
      common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon
      character*80 fndat,fnncs,fnIIdat,fnIIncs                 !qgs-II????????
      common/qgsfname/  fndat, fnncs, ifdat, ifncs
      common/qgsIIfname/fnIIdat, fnIIncs, ifIIdat, ifIIncs     !qgs-II????????
      common/ghecsquel/anquasiel,iquasiel
      common/cbincond/nozero,ibmin,ibmax
      !gauss weights
      data (tgss(2,i),i=1,7)/ .3399810436,.8611363116    ,5*0.     /
      data (wgss(2,i),i=1,7)/ .6521451549,.3478548451    ,5*0.     /
      data (tgss(3,i),i=1,7)/ .2386192,.6612094,.9324700  ,4*0.    /
      data (wgss(3,i),i=1,7)/ .4679139,.3607616,.1713245  ,4*0.    /
      data (tgss(5,i),i=1,7)/ .1488743389,.4333953941,.6794095682
     *                       ,.8650633666,.9739065285    ,2*0.     /
      data (wgss(5,i),i=1,7)/ .2955242247,.2692667193,.2190863625
     *                       ,.1494513491,.0666713443    ,2*0.     /
      data (tgss(7,i),i=1,7)/ .9862838,.9284349,.8272013,.6872929
     *                       ,.5152486,.3191124,.1080549           /
      data (wgss(7,i),i=1,7)/ .03511946,.08015809,.1215186,.1572032
     *                       ,.1855384,.2051985,.2152639           /

      if(iop.eq.1)write(ifmt,'(a)')'default settings ...'
      if(iop.eq.1)goto1001
      if(iop.eq.2)goto1002

c  version

      iversn=int(1.67*100) !version number
      iverso=int(1.67*100) !last official version

c  application

      iappl=1          !choice for application (0,1,2,3,4)

c  model

      model=1
      iquasiel=1       !allow (1) or not (0) quasi-elastic event in model 3

c  file names and units

      fnnx='path/epos '      !path epos name
      fnch='zzz.check '         !check-file name
      fnhi='zzz.histo '         !histo-file name
      fndt='zzz.data '          !data-file name
      fncp='zzz.copy '          !copy-file name
      fnii='zzz.initl '       !initl-file name
      fnid='zzz.inidi '       !inidi-file name
      fndr='zzz.inidr '       !inidr-file name
      fnie='zzz.iniev '       !iniev-file name
      fnrj='zzz.inirj '       !inirj-file name
      fncs='zzz.inics '       !inics-file name
      nfnnx=index(fnnx,' ')-1   !length of path epos name
      nfnch=index(fnch,' ')-1   !length of check-file name
      nfnhi=index(fnhi,' ')-1   !length of histo-file name
      nfndt=index(fndt,' ')-1   !length of data-file name
      nfncp=index(fncp,' ')-1   !length of copy-file name
      nfnii=index(fnii,' ')-1   !length of initl-file name
      nfnid=index(fnid,' ')-1   !length of inidi-file name
      nfndr=index(fndr,' ')-1   !length of inidr-file name
      nfnie=index(fnie,' ')-1   !length of iniev-file name
      nfnrj=index(fnrj,' ')-1   !length of inirj-file name
      nfncs=index(fncs,' ')-1   !length of inics-file name

      if(iop.eq.0)then
      ifop=5     !optns-file unit
      ifmt=6     !std output unit
      ifcx=31    !check-file unit (for open)
      ifch=31    !check-file unit (for write)
      ifhi=35    !histo-file unit
      ifdt=51    !data-file unit
      ifcp=52    !copy-file unit
      endif
c  initial seed

      seedi=0d0   !.ne.0.: seed for random number generator: at start program
      seedj=0d0   !.ne.0.: seed for random number generator: for first event


 1001 continue ! ----> entry iop=1 <-----

      call factoriel


c  some basic things

      nevent=1    !number of events
      nfull=-1          !number of full events (ico+hydro+f.o.)(if different from -1)
      nfreeze=1    !number of freeze out events per epos event
      ninicon=1    ! +-number of events per initial condition
                    ! if positive: keep same b, otherwise generate b each time
      engy=-1      !energy
      elab=-1      !energy lab
      ecms=-1      !energy centre of mass
      ekin=-1      !energy kinetic
      pnll=-1      !momentum
      egymin=6.    !minimum energy
      egymax=2.E+06 !maximum energy
      noebin=1     !number of energy bins
      engmin=0     !minimum energy
      engmax=0     !maximum energy
      iologe=-1     !0=linear bins   1=log bins (cms engy)
      iologl=-1     !0=linear bins   1=log bins (Kinetic engy)

c  printout options

      if(iop.eq.0)iprmpt=-2   !prompt (>0) or not (<0)
      ish=0      !1,2,3,4 ...: more and more detailed messages
      irandm=0   !prints all random numbers (1)
      irewch=0   !rewinds check file before each event (1)
      if(iop.eq.0)iecho=1    !verify option for input reading
      modsho=1   !message all modsho events
      idensi=0   !must be 1 when subr xjden1 is used
      ishevt=0   !minimum event number for certain print options
      iwseed=1   !print out seed (1) or not

c  fragmentation and decay parameters

      pud=0.434  !prob u d
      pude=0.434  !prob u d for ee
      pudk=0.434  !prob u d for kinky string
      pudr=0.434 !prob u d for remnant string
      pudd=0.44   !prob u d for break diffractive
      puds=0.14   !s suppression in diquark break
      strcut=0.00012   !cut factor for diffractive string fragmentation
      diqcut=0.12       !cut factor for diffractive string fragmentation
      pdiqua= 0.21   !qq-qqbar probability
      pdiquae=0.09  !qq-qqbar probability for ee
      pdiquak=0.09   !qq-qqbar probability for kinky string
      pdiquar=0.09   !qq-qqbar probability for remnant string
      pdiquad=0.09   !qq-qqbar probability for diffractive remnant string
      ptfra= 0.30  !string break pt
      ptfrae=0.30  !string break pt for ee
      ptfrak=0.30  !string break pt for kinky string in hadronic ia
      ptfrakd=0.30  !string break pt for kinky string diquark break
      ptfrair=0.30 !string break pt for inel remnant string
      ptfradr=0.40 !string break pt for diffr remnant string
      ptfrasr=0.25 !string break del pt for strange breakup of remn string
      ptfraqq=0.40 !string break pt for diquark break
      pbreak=-1.  !break-parameter
      pbreake=-1.  !break-parameter ee
      pbreakk=-1.  !break-parameter for kinky string
      pbreakr=-1. !break-parameter for remnant string
      pbreakd=0.27 !break-parameter for diff remnant string

c  fragmentation and decay options

      ndecay=0   !ndecay suppresses certain resonance decays
                 !0000001: all resonance decays suppressed
                 !0000010: k_short/long (+-20) decays suppressed
                 !0000100: lambda (+-2130) decays suppressed
                 !0001000: sigma (+-1130,+-2230) decays suppressed
                 !0010000: cascade (+-2330,+-1330) decays suppressed
                 !0100000: omega (+-3331) decays suppressed
                 !1000000: pi0 (110) decays suppressed
                 !also several "1"s may be used with obvious meaning
      maxres=99999 !maximum resonance spin
      aouni=0.   !technical parameter for development

c  lepton-nucleon and e+e- options

      iolept=1     !q2-x ditribution calculated (1) or taken from table (<0)
      ydmin=0      ! range of y
      ydmax=0
      qdmin=0      ! range of q**2
      qdmax=0
      themin=0     !minimum scattering angle
      themax=0     !maximum scattering angle
      elomin=0     !minimum energy of outgoing lepton
      elepti=0     !incoming lepton energy
      elepto=0     !outgoing lepton energy
      angmue=3.9645/360.*2*3.14159 !mue angle
      icinpu=0
      itflav=0     ! initial flavor for e+e-
      idisco=0     !deep inelastic contributions
                   !0=all, 1=direct-light, 2=direct-charm, 3=resolved

c  hadron-hadron options +++

      isetcs=2    !option to obtain pomeron parameters
                  ! 0.....determine parameters but do not use Kfit
                  ! 1.....determine parameters and use Kfit
                  ! else..get from table
                  !         should be sufficiently detailed
                  !          say iclegy1=1,iclegy2=99
                  !         table is always done, more or less detailed!!!
                  !and option to use cross ection tables
                  ! 2....tabulation
                  ! else...not
      iclegy1=1   !energy class lower limit ( 1 minimal    1 realistic    )
      iclegy2=99   !energy class upper limit ( 1  option   99 use of table )
      isigma=1    !option for xsection printing (always calculated now)
                  !  0=no, 1=yes

c  hadron-hadron options

      idprojin=1120 !projectile particle id
      idtargin=1120 !target particle id
      idproj=1120 !projectile particle id
      idtarg=1120 !target particle id
      iregge=0    !consider reggeons (1=yes 0=no)
      isopom=1    !consider soft pomerons (1=yes 0=no)
      ishpom=1    !consider semihard pomerons (1=yes 0=no)
      iscreen=1   !consider screening corrections (1=yes 0=no)
      isplit=1    !consider splitting corrections (1=yes 0=no)
      irzptn=0    !recalculate Zptn (1=yes 0=no)  ????????maybe obsolete??????????
      irmdrop=1    !consider droplet for remnants (1=yes 0=no)
      nprmax=10000 !maximum number of pomerons/reggeons
      iemspl=0
      intpol=3     !number of points for interpolation in psfli, psfaz
      ioems=2
      iomega=1    !option for omega calculation (if 2, no diffraction in G)
        !hadron excitation classes (used in psvin)
      icdp=2     !projectile hadron excitation
      icdt=2     !target hadron excitation
              !hadron classes: 1=pions,2=nucleons,3=kaons
      iclpro1=1   !projectile hadron class lower limit
      iclpro2=3   !projectile hadron class upper limit
      icltar1=2   !target hadron class lower limit (should not be change (see epos-sem))
      icltar2=2   !target hadron class upper limit (should not be change (see epos-sem))
      egylow=1.5  !lower limit of lowest energy bin
      delh=0.25   !effective overcriticity for the hard pom (techn param)
      factgam=1         !enhancement factor for gammas

c  hadron-hadron parameters +++

      alpfomi=   8.    !normalization of function fom for Phi^n (z0=alpfomi)
      betfom=    9.d0   !slope of function fom for Phi^n
      gamfom=    5.d0   !Z slope of function fom for Phi^n
      betpom=    0.25    !gluon structure function parameter for the so
      glusea=    0.1   !sea quarks / (sea quarks + gluons) at x -> 0
      r2had(2)=  1.   !r2 proton
      r2hads(2)= 1.06     !diff corr factor proton
      slopom=    0.06    !alpha prime pomeron
      slopoms=   0.175    !alpha prime pomeron dif
      gamhad(2)= 1.1    !gamma proton           increase->  sig up, n up, softPom up
      gamhadsi(2)=-1.   !gamma soft proton (<0 = same as gamhad)
      gamtil=    0.08   !                       increase -> sig up, n up, hard Pom up
      alppar=    0.5  !alpha particip (not 1 !) increase -> sig up, n down, width y down
      alppom=    1.06 !alpha pomeron
      ptsend=    0.8    !string end pt
      ptsendi=   0.05   !string end pt for non excited state
      ptsemx=    2.     !cut factor for gaussian distribution
      ptdiff=    0.27   !pt for diffraction
      q2min=     4.0    !q**2 cutoff
      q2ini=     0.25   !resolution scale
      q2fin=     0.02   !q**2 cutoff timelike      decrease ->  high pt down
      amdrmax=   2.5    !maximum mass leading droplet
      facdif=    0.35   !factor for diffractive profile
      reminv=    0.33     !remnant inversion probability for strangeness
      zrminc=    0.033  !remnant inversion probability string mass evolution
      zrmmax=    10.    !remnant inversion probability system mass evolution
      edmaxi=    1.e12   !defines edmax in epos-sem
      epmaxi=    1.e12   !defines epmax in epos-sem

c  hadron-hadron parameters

      qcdlam=.04           !lambda_qcd squared
      naflav=4          !number of active flavors
      factk=2.         !k-factor value
      alfe=1./137.
      pt2cut=0.         !p_t-cutoff for the born process
      rstrau(1)=1.   !pion !effective ratio of u sea over u sea basis
      rstrad(1)=1.         !effective ratio of d sea over u sea basis
      rstras(1)=0.8         !effective ratio of s sea over u sea basis
      rstrau(2)=1.   !nucl !effective ratio of u sea over u sea basis
      rstrad(2)=1.         !effective ratio of d sea over u sea basis
      rstras(2)=0.8         !effective ratio of s sea over u sea basis
      rstrau(3)=1.   !kaon !effective ratio of u sea over u sea basis
      rstrad(3)=1.         !effective ratio of d sea over u sea basis
      rstras(3)=0.8         !effective ratio of s sea over u sea basis
      rstrau(4)=1.   !j/psi!effective ratio of u sea over u sea basis
      rstrad(4)=1.         !effective ratio of d sea over u sea basis
      rstras(4)=0.8         !effective ratio of s sea over u sea basis
      rstrasi=0.0      !effective ratio of strange sea over u sea increase
      wgtval=0.          !weight for valq - antiq as soft string ends
      wgtsea=1.          !weight for seaq - antiq as soft string ends
      wgtdiq=0.          !weight for diq - antidiq as soft string ends
      wgtqqq(1)=0.17          !weight for seaq - valqq as soft string end for pion
      wgtqqq(2)=0.17          !weight for seaq - valqq as soft string end for nucleon
      wgtqqq(3)=0.17          !weight for seaq - valqq as soft string end for kaon
      wgtqqq(4)=0.17          !weight for seaq - valqq as soft string end for J/Psi
      exmass=0.02         !excitation mass for remnant
      r3pom=0.01      !triple pomeron coupling
        r3pomi=r3pom  !store
      r4pom=0.001     !4-pomeron coupling
        r4pomi=r4pom  !store
      gfactor=0.4     !exp(-gfactor*Z) in front of eps in om5jk (to choose soft or hard)
      gwidth=0.85   !diff relative b-width
!     gamhad(1) defined in psaini
!     gamhad(3) defined in psaini
      gamhadsi(1)=-1.    !gamma soft pion
      gamhadsi(3)=-1.    !gamma soft kaon
      gamhadsi(4)=-1.    !gamma soft charm
      r2had(1)=1.    !r2 pion
      r2had(3)=1.30   !r2 kaon
      r2had(4)=0.     !r2 charm
      r2hads(1)=0.85  !diff corr factor pion
      r2hads(3)=1.15  !diff corr factor kaon
      r2hads(4)=1.25  !diff corr factor kaon
      chad(1)=1.      !c pion
      chad(2)=1.      !c proton
      chad(3)=1.      !c kaon
      chad(4)=1.      !c charm
      wdiff(1)=0.36   !diffractive probability
      wdiff(2)=0.645  !diffractive probability
      wdiff(3)=0.145  !diffractive probability
      wdiff(4)=0.1    !diffractive probability
      alplea(1)=0.75  !alpha leading pion
      alplea(2)=1.25  !alpha leading proton
      alplea(3)=1.    !alpha leading kaon
      alplea(4)=1.    !alpha leading jpsi
      rexndf=1.       !relative value of rexndi compare to rexdif if >0
c following parameters recalculated in xsigma ...
      rexdifi(1)=-0.88  !remnant excitation probability diffractive pion
      rexdifi(2)=-0.81  !remnant excitation probability diffractive proton
      rexdifi(3)=-0.81  !remnant excitation probability diffractive kaon
      rexdifi(4)=1.     !remnant excitation probability diffractive jpsi
      rexndii(1)=0.5   !remnant excitation probability nondiffractive pion
      rexndii(2)=0.81   !remnant excitation probability nondiffractive proton
      rexndii(3)=0.65   !remnant excitation probability nondiffractive kaon
      rexndii(4)=1.   !remnant excitation probability nondiffractive jpsi
c ... up to here.
      xmaxremn=0.6      !maximum momentum fraction allowed to be given for remnant mass determination (very important for multiplicity and xf distri at low energy)
      xmaxdiff=0.18      !maximum momentum fraction allowed to be given for diffractif remnant mass determination (very important for multiplicity and xf distri at low energy)
      rexres=0.5      !remnant excitation probability in nucleus
      zbarfl=1.0      !remnant baryon flow parameter (if negative)
      ptfrab=0.30    !string pt for remnant baryon (in conn with zbarfl)
      alpdif=0.6666     !alpha mass diffractive for cross section and metropolis
      alpdi=0.6666      !alpha mass diffractive
      alpndi=1.05       !alpha mass nondiffractive
      alpsea=0.3      !alpha string end x for sea parton
      alpval=0.3      !alpha string end x for valence parton
      alpdiq=0.3      !alpha string end x for sea diquark
      ammsqq=0.28     !minimum mass string quark quark
      ammsqd=1.08    !minimum mass string quark diquark
      ammsdd=1.88    !minimum mass string diquark diquark
      delrex=.01     !mass to be added to the above masses
      cumpom=0       !cutoff mass for virtual pomerons
      zbcut=1.e-10    !minimum probability for bkmxndif calculation (for Z_p/t)
      alpdro(1)=1.25   !factor for minimum mass of leading droplet
      alpdro(2)=0.6   !pt of leading droplet
      alpdro(3)=alpndi  !alpha mass of leading droplet

      iodiba=0.      ! if iodiba=1, study H-Dibaryon (not used (see ProRef in epos-ems ????)
      bidiba=0.030   ! epsilon of H-Dibaryon

c screening splitting +++

      epscrw=0.18      !overall factor for Z (Zsame,Zother)-> pp xsect .... w_Z
      epscrv=0.0       !low energy value for Zother   ..................... w_M
      epscrp=0.75       !b width param for Z     -> pp xsection ............ w_B
      egyscr=10.       !screening minimum energy -> pp xsection ........... s_M
      epscrs=0.02       !screening power increase soft  -> pp xsctn ........ alp_S
      epscrh=1.75      !screening power increase hard  -> pp xsctn ........ alp_H
      epscrx=0.35       !screening power maxi           -> pp xsctn
      zidinc= 0.2       !inner diquark modif factor increase
      zidmax= 0.6       !inner diquark modif factor max
      zodinc= 0.      !outer diquark modif factor increase
      zodmax= 0.       !outer diquark modif factor max
      zosinc= 0.      !outer strange modif factor increase
      zosmax= 0.       !outer strange modif factor max
      zipinc= 0.4    !inner pt modif factor increase
      zipmax= 0.9      !inner modif factor max
      zopinc= 0.07    !outer pt and end pt modif factor increase
      zopmax= 1.      !outer pt and end pt modif factor max
      zoeinc=0.44         !cutoff timelike evol modif factor increase
      zoemax=4.      !cutoff timelike evol modif factor max

c Reggeon parameters  (not used)

      alpreg=0.734   !alpha_reggeon
      sloreg=0.499   !slope_reggeon
      gamreg=16.46   !gamma_reggeon
      r2reg=0.613    !r^2_reggeon

c  masses

      amhadr(1)=.14            !pion mass
      amhadr(2)=.939           !nucleon mass
      amhadr(3)=.496           !kaon mass
      amhadr(4)=1.868          !d-meson mass
      amhadr(6)=1.116          !lambda mass
      amhadr(5)=2.27           !lambda_c mass
      amhadr(7)=.548           !eta mass
      amhadr(8)=3.5            !J/psi mass
      qcmass=1.6               !c-quark mass
c      amhdibar=2.200           !h-dibaryon mass       !not used any more
      qmass(0)=0.980           !diquark effective bounding energy (for pt distribtions)
      isospin(0)=0
      qmass(1)=0.001           !u quark effective mass (for pt distribtions)
      isospin(1)=1
      qmass(2)=0.003           !d quark effective mass (for pt distribtions)
      isospin(2)=-1
      qmass(3)=0.140           !s quark effective mass (for pt distribtions)
      isospin(3)=0
      qmass(4)=qcmass          !c quark effective mass (for pt distribtions)
      isospin(4)=0
      call idmass(5,qbmass)
      qmass(5)=qbmass          !b quark effective mass (for pt distribtions)
      isospin(5)=0
      call idmass(6,qtmass)
      qmass(6)=qtmass          !t quark effective mass (for pt distribtions)
      isospin(6)=0

c  nucleus-nucleus

      iokoll=0      !fix # of collisions (1)
      laproj=0      !projectile charge number
      maproj=0      !projectile mass number
      latarg=0      !target charge number
      matarg=0      !target mass number
      core=0.34     !hard core distance(2*0.17)
c      ncolmx=100000 !maximum number of collisions
      fctrmx=10     !parameter determining range for density distribution
      bmaxim=10000  !maximum impact parameter
      bminim=0.     !minimum impact parameter
      phimax=2*3.1415927 !maximum phi
      phimin=0      !minimum phi
      ionudi=0      !nuclear diffraction included (1) or not (0)

c rescattering parameters +++

      iorsce=0      !color exchange turned on(1) or off(0)
      iorsdf=3      !droplet formation turned on(>0) or off(0)
      iorshh=0      !other hadron-hadron int. turned on(1) or off(0)
      iocluin=1     !include inwards moving corona segments into clusters (1) or not (0)
      ioquen=0      !jet quenching option (0=no)
      taumin=1.0    !min tau for rescattering
      taurea=0.001    !reaction time (=formation time)
      nsegsu=30     !number of segments per subcluster
      nsegce=4      !number of segments per cell
      ptclu=3  !obove this: string segments escape cluster
      amimfs=1.0    !below this: elastic
      amimel=0.050  !below this: nothing
      delamf=1      !above this: color exch  !cutoffs for kinetic energy
      deuamf=1      !above this: nothing     !mass - minimum mass of pair
      epscri(1)=0.22   !critical energy density
      epscri(3)=0.08   !critical energy density

c rescattering parameters

      amsiac=.8     !interaction mass
      amprif=0      !print option
      delvol=1      !print option
      deleps=1      !print option
      deltau=0.2    !initial delta_tau for space-time evolution
      factau=1.05   !factor for delta_tau for space-time evolution
      numtau=80     !number of tau steps for space-time evolution
      dlzeta=0.5    !delta_zeta for longitudinal droplet binning
      etafac=1.0    !factor determining inner range
      facnuc=0.0    !factor for nuclear size to determine inner range
      hacore=0.8    !hadron compressibility factor
      cepara=0.03   !parameter for excitation for color exchange
      dscale=1.    !scale parameter for hadron-hadron
      iceopt=1      !real color exchange (1) or just excitation (0)

c coupling with urqmd

      iurqmd=0    ! call eposurqmd (1) or not

c initial conditions for hydro

      ispherio=0    !call spherio
      cutico=3      !cutoff parameter for smoothing kernel in epos-ico
                    !  (as small as possible with stable results)
      dssico=0.2    !s step size for string integration in epos-ico
                    !  (as big as possible with stable results)
      icocore=0     !consider core initial condition (1 or 2)
      icotabm=0     !make table for initial condition stuff
      icotabr=0     !read table for initial condition stuff

c  cluster decay

      ioclude=1     !cluster decay option
      amuseg=5.0    !min mass for radial boost
      yradmx= 0.77   !max radial collective boost
      yradmi= 0.15   !max radial collective boost increase
      yradpp=1.2    !max radial collective boost pp basis
      yradpi=0.4    !max radial collective boost pp increase
      yradpx=1.2    !max radial collective boost pp maxi
      facecc=0.45   !eccentricity parameter
      ylongmx=-1.   !max long collective boost ( < 0 -> take from jintpo )
      rcoll=0.0     !radial collective flow param
      bag4rt=0.200  !bag constant ^1/4
      taunll=1.0    !decay time (comoving frame)
      vrad=0.3
      facts=0.3       !gamma-s factor
      factb=1       !gamma-s factor baryons
      factq=1       !gamma-qqbar
      inbxxx=0      !1: fast nboby via hnbxxx   ????? maybe obsolete???????

c  droplet decay initializations

         asuhax(1)=1.134  !two lowest multiplets
         asuhax(2)=1.301  !two lowest multiplets
         asuhax(3)=1.461  !two lowest multiplets
         asuhax(4)=1.673  !two lowest multiplets
         asuhax(5)=0.7700 !two lowest multiplets   rho
         asuhax(6)=0.8920 !two lowest multiplets   K*
         asuhax(7)=1.2320 !two lowest multiplets
         asuhay(1)=0.940  !lowest multiplet
         asuhay(2)=1.200  !lowest multiplet
         asuhay(3)=1.322  !lowest multiplet
         asuhay(4)=1.673  !lowest multiplet
         asuhay(5)=0.1400 !lowest multiplet
         asuhay(6)=0.4977 !lowest multiplet
         asuhay(7)=1.2320 !lowest multiplet

c  droplet specification

      keu=0     !u flavour of droplet
      ked=0     !d flavour of droplet
      kes=0     !s flavour of droplet
      kec=0     !c flavour of droplet
      keb=0     !b flavour of droplet
      ket=0     !t flavour of droplet
      tecm=10   !energy of droplet
      volu=70   !volume of droplet

c  metropolis and grand canonical

      iospec=10 !option for particle species
      iocova=1  !covariant (1) or noncovariant (2) phase space integral
      iopair=2  !single pair (1) or double pair (2) method
      iozero=65 !relative weight of zeros (compared to hadrons)
                ! (-1) nspecs
                ! (-2) nspecs/sqrt(tecm/volu)
      ioflac=1  !test multipl distr without (1) or with (2) flavour conserv
                !  (2 only good for nspecs=3,7)
      iostat=1  !use boltzmann (1) or quantum (0) statistics in hgc-routines
      ioinco=1  !call hnbmin for initial configuration (0)
                !call hgcnbi for initial configuration to generate better
                !initial configuration (1)
                !call hgcnbi for initial configuration to generate optimal
                !initial configuration (2)
      iograc=1  !call hgcaaa in case of ioinco=0 (1)
      epsgc=2.  !required accuracy in hgcaaa 10**(-epsgc)
      iocite=0  !perform counting at metropolis iterations (1) or not (else)
      ioceau=0  !perform counting for exp. autocorrel. time (1) or not (else)
      iociau=0  !perform counting for int. autocorrel. time (1) or not (else)
      ioinct=0  !test grand canonical metropolis sampling (1)
                !to plot results call xhgccc, xhgcfl and xhgcam
      ioinfl=1  !conserve flavor in initial configuration in hgcnbi (1)
                !do not conserve flavor (0)
                !do not conserve flavor and energy (-1)
      iowidn=2  !width of total multiplicity distribution in hgcnbi
                ! sigma_tot -> sigma_tot/iowidn
                ! >0 unnormalized
                ! <0 normalized
      ionlat=2  !determine nlattc ,old method (0)
                !or determine nlattc in hgcnbi as:
                ! (1) max(1.3*<N>,<N>+2*sigma,6)
                ! (2) max(1.5*<N>,<N>+3*sigma,6)
                ! (3) max(2.0*<N>,<N>+4*sigma,6)
      iomom=1   !number of momenta to be changed in hnbodz
      ioobsv=0  !observable for autocorrelation time calculation
                !0: total multiplicity
                !else: particle id for particle species
      iosngl=0  !event # for which counting at metropolis iterations is done
      iorejz=0  !reject pair exchange with only zeros (1) or not (else)
      iompar=4  !parameter for windowing algorithm
      iozinc=0  !if iozevt>0: modifies iozero for testing (in sr hgcnbi)
      iozevt=0  !if >0: modifies iozero for testing (in sr hgcnbi)
      nadd=0    !number of pi0s added to minimum initial configuration
      iterma=-6 !>0: maximum number of iterations
                !<0: - number of iterations per corr time
      iterpr=10000 !iter-increment for printout
      iterpl=1  !iter-increment for plot
      iternc=50 !iterations not counted for analysis
      epsr=1e-4 !required accuracy in sr hnbraw
      keepr=1   !keep most random numbers rnoz in hnbodz (1)
                !  or update all (0)

c  strangelets

      iopenu=1      !option for minimum energy
                    !1: sum of hadron masses
                    !2: bag model curve with minimum at nonzero strangen
      themas=.51225 !parameter theta in berger/jaffe mass formula

c tests

      iotst1=0     !test
      iotst2=0     !test
      iotst3=0     !test
      iotst4=0     !test

c  jpsi

      jpsi=0     !jpsi to be produced (1) or not (0)
      jpsifi=0   !jpsi final state interaction (1) or not (0)
      sigj=0.2   !jpsi nucleon cross section [fm**2]
      taumx=20   !max time for jpsi evolution
      nsttau=100 !time steps for jpsi evolution
      ijphis=0   !fill jpsi histograms (1) or not (0)

c  analysis: intermittency, space-time, droplets, formation time

      ymximi=2   !upper limit for rapidity interval for intermittency analysis
      imihis=0   !fill intermittency histograms (1) or not (0)
      isphis=0   !fill space-time histograms (1) or not (0)
      iologb=0   !0=linear bins   1=log bins
      ispall=1   !xspace: all ptls (1) or only interacting ptls (else)
      wtmini=-3  !tmin in xspace
      wtstep=1   !t-step in xspace
      iwcent=0   !only central point (1) or longitudinal distr (else) in xspace
      iclhis=0   !fill droplet histograms (1) or not (0)
      iwtime=0   !fill formation time histogram (1) or not (else)
      wtimet=100 !max time in wtime
      wtimei=0   !max mass in wtime
      wtimea=1000 !max mass in wtime


c  storing
      maxrec(1)=7
      irecty(1,1)=1
      irecty(2,1)=2
      irecty(3,1)=3
      irecty(4,1)=4
      irecty(5,1)=5
      irecty(6,1)=6
      irecty(7,1)=7
      maxrec(2)=14
      irecty(1,2)=1
      irecty(2,2)=2
      irecty(3,2)=3
      irecty(4,2)=4
      irecty(5,2)=5
      irecty(6,2)=6
      irecty(7,2)=7
      irecty(8,2)=8
      irecty(9,2)=9
      irecty(10,2)=10
      irecty(11,2)=11
      irecty(12,2)=12
      irecty(13,2)=13
      irecty(14,2)=14


c  other

      gaumx=8    !range for gauss distribution
      nclean=1   !clean /cptl/ if nclean > 0
                 !(not for part with istptl<istmax if nclean=1 (do not change analysis)
                 ! for all part with ist.ne.0 if nclean > 1 (particle nb reduce at max))
      istore=0   !0: no storage to data-file
                 !1: epos     standard
                 !2: OSC1997A standart
                 !3: OSC1999A standart
                 !5: calls routine ustore (modifiable by the user)
      ioidch=1   !id choice for storage to data-file
      iframe=0   !frame specification production run
      jframe=0   !frame specification analysis
      kframe=0   !frame specification analysis (2nd frame)
                 ! 1:total
                 !11:nucleon-nucleon
                 !12:target
                 !21:gamma-nucleon
                 !22:lab
                 !32:sphericity
                 !33:thrust
      irescl=1   !momentum rescaling (1) or not (0)
      ifrade=1   !suppression of fragmentation and decay (0)
      idecay=1   !suppression of decay (0)
      jdecay=1   !suppression of cluster decay (0), concerns only ity=60 cluster
      ntrymx=10  !try-again parameter
      istmax=1   !analyse only istptl <= istmax
      irdmpr=0   !random sign for projectile if 1

c  constants

      pi=3.1415927
      pii=1./(2*pi**2)
      hquer=0.197327
      prom=.94
      piom=.14
      ainfin=1e31

c air

      airanxs(1)=14.007
      airznxs(1)=7.
      airwnxs(1)=0.781
      airanxs(2)=15.999
      airznxs(2)=8.
      airwnxs(2)=.21
      airanxs(3)=39.948
      airznxs(3)=18.
      airwnxs(3)=0.009
      airavanxs=airanxs(1)*airwnxs(1)+airanxs(2)*airwnxs(2)
     &         +airanxs(3)*airwnxs(3)
      airavznxs=airznxs(1)*airwnxs(1)+airznxs(2)*airwnxs(2)
     &         +airznxs(3)*airwnxs(3)

c  zero initializations

      ixgeometry=0
      ixtau=0
      iEmsB=0
      iEmsBg=0
      iEmsPm=0
      iEmsPx=0
      iEmsSe=0
      iEmsDr=0
      iEmsRx=0
      iEmsI2=0
      iEmsI1=0
      iSpaceTime=0
      nemsi=0
      facmss=1.
      nstmax=0
      do 6 i=1,99
      prob(i)=0
      do 6 j=1,2
      icbac(i,j)=0
6     icfor(i,j)=0
      imsg=0
      do j=1,mxjerr
       jerr(j)=0
      enddo
      ntevt=0
      nrevt=0
      naevt=0
      nrstr=0
      nrptl=0
      nptlu=0
      do itau=1,mxtau
      volsum(itau)=0
      vo2sum(itau)=0
      nclsum(itau)=0
      do ivol=1,mxvol
      do ieps=1,mxeps
      clust(itau,ivol,ieps)=0
      enddo
      enddo
      enddo
      iutotc=0
      iutote=0
      nopen=0
      nopenr=0
      knxopen=0
      kchopen=0
      khiopen=0
      kdtopen=0
      klgopen=0
      ifdat=0
      ifncs=0
      xpar1=0.
      xpar2=0.
      xpar3=0.
      xpar4=0.
      xpar5=0.
      xpar6=0.
      xpar7=0.
      xpar8=0.
      if(iop.eq.0)khisto=0
      nrclu=0
      nrnody=0
      do n=1,mxnody
      nody(n)=0
      enddo
      nrpri=0
      do n=1,mxpri
      subpri(n)='      '
      ishpri(n)=0
      enddo
      nctcor=0
      ncttim=0
      do n=1,matau
      tauv(n)=0
      enddo
      ncnt=0
      nrnucl(1)=0
      nrnucl(2)=0
      do i=1,mxnucl
      rnucl(i,1)=0
      rnucl(i,2)=0
      rnuclo(i,1)=0
      rnuclo(i,2)=0
      enddo
      inicnt=0
      accept=0.
      reject=0.
      do n=1,matau
      tauv(n)=0.
      enddo
      anquasiel=0.
      nglacc=0
      ifset=1

 1002 continue ! ----> entry iop=2 <----


c  analysis

      xvaria='numptl'
      yvaria='ycmptl'
      normal=11
      xminim=-100
      xmaxim=100
      nrbins=100
      hisfac=1
      do nr=1,mxbins
      do l=1,5
      ar(nr,l)=0
      enddo
      enddo
      nozero=0
      ibmin=1
      ibmax=1e8

      return
      end

c-----------------------------------------------------------------------
      subroutine ustore
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c     modifiable by the user
c-----------------------------------------------------------------------
      include 'epos.inc'

c  count the number of particles to be stored (--> nptevt)

      nptevt=0
      do i=1,nptl
      if(istptl(i).le.istmax)nptevt=nptevt+1
      enddo

c  store event variables:

      write(ifdt,*)nrevt,nptevt,bimevt,phievt,kolevt
     *,pmxevt,egyevt,npjevt,ntgevt

c     nrevt.......... event number
c     nptevt ........ number of (stored!) particles per event
c     bimevt ........ absolute value of impact parameter
c     phievt ........ angle of impact parameter
c     kolevt ........ number of collisions
c     pmxevt ........ reference momentum
c     egyevt ........ pp cm energy (hadron) or string energy (lepton)
c     npjevt ........ number of primary projectile participants
c     ntgevt ........ number of primary target participants
c     npnevt ........ number of primary projectile neutron spectators
c     nppevt ........ number of primary projectile proton spectators
c     ntnevt ........ number of primary target neutron spectators
c     ntpevt ........ number of primary target proton spectators
c     jpnevt ........ number of absolute projectile neutron spectators
c     jppevt ........ number of absolute projectile proton spectators
c     jtnevt ........ number of absolute target neutron spectators
c     jtpevt ........ number of absolute target proton spectators

      do i=1,nptl

      if(istptl(i).le.istmax)then !store events with istptl < istmax

c  store particle variables:

      write(ifdt,*)
     *i,idptl(i),pptl(1,i),pptl(2,i),pptl(3,i),pptl(4,i),pptl(5,i)
     *,iorptl(i),jorptl(i),istptl(i)
     *,xorptl(1,i),xorptl(2,i),xorptl(3,i),xorptl(4,i)
     *,tivptl(1,i),tivptl(2,i)

c     i ............. particle number
c     idptl(i) ...... particle id
c     pptl(1,i) ..... x-component of particle momentum
c     pptl(2,i) ..... y-component of particle momentum
c     pptl(3,i) ..... z-component of particle momentum
c     pptl(4,i) ..... particle energy
c     pptl(5,i) ..... particle mass
c     iorptl(i) ..... particle number of father (if .le. 0 : no father)
c     jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
c     istptl(i) ..... generation flag: last gen. (0) or not (1)
c     ityptl(i) ..... particle type (string, remnant ...)
c     xorptl(1,i) ... x-component of formation point
c     xorptl(2,i) ... y-component of formation point
c     xorptl(3,i) ... z-component of formation point
c     xorptl(4,i) ... formation time
c     tivptl(1,i) ... formation time (always in the pp-cms!)
c     tivptl(2,i) ... destruction time (always in the pp-cms!)

      endif
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine bstora
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c-----------------------------------------------------------------------
      include 'epos.inc'
      common/record/maxrec(2),irecty(30,2)
      character code*8,version*8,frame*4
      code='epos   '
      write(version,'(f5.2,3x)')iversn/100.
      if(iframe.eq. 1)frame='ttcm'
      if(iframe.eq.11)frame='nncm'
      if(iframe.eq.12)frame='targ'
      if(iframe.eq.21)frame='gncm'
      if(iframe.eq.22)frame='lncm'
      ntest=1
      if (istore.eq.2) then     ! OSC1997A
        write (ifdt,'(a)') 'OSC1997A'
        write (ifdt,'(a)') 'final_id_p_x'
        write(ifdt,100) code,version
     *       ,maproj,laproj,matarg,latarg,frame,engy,ntest
 100    format(2(a8,'  '),'(',i3,',',i3,')+(',i3,',',i3,')',
     *       '  ',a4,'  ',e10.4,'  ',i8)
        maxrec(1)=4
        irecty(1,1)=1           !nevt
        irecty(2,1)=2
        irecty(3,1)=3
        irecty(4,1)=4
        maxrec(2)=11
        irecty(1,2)=1           !nr
        irecty(2,2)=2           !id
        irecty(3,2)=3           !px
        irecty(4,2)=4           !py
        irecty(5,2)=5           !pz
        irecty(6,2)=6           !E
        irecty(7,2)=7           !M
        irecty(8,2)=11          !x
        irecty(9,2)=12          !y
        irecty(10,2)=13         !z
        irecty(11,2)=14         !t
      elseif(istore.eq.3) then
 201    format('# ',a)
        write (ifdt,201) 'OSC1999A'
        if(istmax.eq.0) then
          write (ifdt,201) 'final_id_p_x'
        elseif(istmax.ge.2) then
          write (ifdt,201) 'full_event_history'
        endif
 202    format('# ',a8,' ',a8)
        write(ifdt,202) code,version !3rd line
 203    format('# (',i3,',',i3,')+(',i3,',',i3,')',
     *       '  ',a4,'  ',e10.4,'  ',i8)
        write(ifdt,203) maproj,laproj,matarg,latarg,frame,engy,ntest
        maxrec(1)=5
        irecty(1,1)=2           !nevt
        irecty(2,1)=0           !zero
        irecty(3,1)=1           !additional information
        irecty(4,1)=3           !additional information
        irecty(5,1)=4           !additional information
        maxrec(2)=12
        irecty(1,2)=1           !nr
        irecty(2,2)=2           !id
        irecty(3,2)=10          !ist
        irecty(4,2)=3           !px
        irecty(5,2)=4           !py
        irecty(6,2)=5           !pz
        irecty(7,2)=6           !E
        irecty(8,2)=7           !M
        irecty(9,2)=11          !x
        irecty(10,2)=12         !y
        irecty(11,2)=13         !z
        irecty(12,2)=14         !t
                                ! nin nout [optional information]
                                ! ipart id ist px py pz p0 mass x y z t [optional information]

      endif

      end

c-----------------------------------------------------------------------
      subroutine bstore
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c-----------------------------------------------------------------------

      include 'epos.inc'
      common/record/maxrec(2),irecty(30,2)
      common/dimensi/k2(100)

      nptevt=0
      do n=1,nptl
        iok=1               !idcode simple
        if(istptl(n).gt.istmax) then
          iok=0
        endif
      if (iok.eq.1) nptevt=nptevt+1
      enddo
 11   format (i6,' ',$)
 12   format (e12.6,' ',$)
      do i=1,maxrec(1)
        l=irecty(i,1)
        if(l.eq.0)write(ifdt,21) 0
        if(l.eq.1)write(ifdt,11)nrevt
        if(l.eq.2)write(ifdt,11)nptevt
        if(l.eq.3)write(ifdt,12)bimevt
        if(l.eq.4)write(ifdt,12)phievt
        if(l.eq.5)write(ifdt,11)kolevt
        if(l.eq.6)write(ifdt,12)pmxevt
        if(l.eq.7)write(ifdt,12)egyevt
        if(l.eq.8)write(ifdt,11)npjevt
        if(l.eq.9)write(ifdt,11)ntgevt
        if(l.eq.10)write(ifdt,11)npnevt
        if(l.eq.11)write(ifdt,11)nppevt
        if(l.eq.12)write(ifdt,11)ntnevt
        if(l.eq.13)write(ifdt,11)ntpevt
        if(l.eq.14)write(ifdt,11)jpnevt
        if(l.eq.15)write(ifdt,11)jppevt
        if(l.eq.16)write(ifdt,11)jtnevt
        if(l.eq.17)write(ifdt,11)jtpevt
        if(l.eq.20)write(ifdt,12)amproj
        if(l.eq.21)write(ifdt,12)amtarg
      enddo
      write (ifdt,*)            !RETURN
 21   format (i6,' ',$)
 22   format (e12.6,' ',$)
 23   format (i10,' ',$)
      do n=1,nptl
        iok=1                   !idcode simple
        if(istptl(n).gt.istmax) then
          iok=0
        endif
        if (iok.eq.1) then
          id=idptl(n)
          if(istore.eq.2.or.istore.eq.3.or.ioidch.eq.2)then
            id=idtrafo('nxs','pdg',idptl(n))
          endif
          do i=1,maxrec(2)
            l=irecty(i,2)
            if(l.eq.0)write(ifdt,21) 0
            if(l.eq.1)write(ifdt,21) n
            if(l.eq.2)write(ifdt,23) id
            if(l.eq.3.or.l.eq.17)write(ifdt,22) pptl(1,n)
            if(l.eq.4.or.l.eq.17)write(ifdt,22) pptl(2,n)
            if(l.eq.5.or.l.eq.17)write(ifdt,22) pptl(3,n)
            if(l.eq.6.or.l.eq.17)write(ifdt,22) pptl(4,n)
            if(l.eq.7.or.l.eq.17)write(ifdt,22) pptl(5,n)
            if(l.eq.8)write(ifdt,21) iorptl(n)
            if(l.eq.9)write(ifdt,21) jorptl(n)
            if(l.eq.10)write(ifdt,21) istptl(n)
            if(l.eq.11.or.l.eq.18)write(ifdt,22) xorptl(1,n)
            if(l.eq.12.or.l.eq.18)write(ifdt,22) xorptl(2,n)
            if(l.eq.13.or.l.eq.18)write(ifdt,22) xorptl(3,n)
            if(l.eq.14.or.l.eq.18)write(ifdt,22) xorptl(4,n)
            if(l.eq.19)write(ifdt,22) dezptl(n)
            if(l.eq.21)write(ifdt,21) ifrptl(1,n)
            if(l.eq.22)write(ifdt,21) ifrptl(2,n)
            if(l.eq.23)write(ifdt,21) ityptl(n)
            if(l.eq.15) then
              if(iorptl(n).gt.0)then
                write(ifdt,23) idptl(iorptl(n))
              else
                write(ifdt,23) 0
              endif
            endif
            if(l.eq.16) then
              if(jorptl(n).gt.0)then
                write(ifdt,23) idptl(jorptl(n))
              else
                write(ifdt,23) 0
              endif
            endif
          enddo
          write (ifdt,*)        !RETURN
        endif
      enddo
 9999 return
      end


c-----------------------------------------------------------------------
      subroutine bread
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c-----------------------------------------------------------------------

      include 'epos.inc'
      common/record/maxrec(2),irecty(30,2)
      character*255 line
      dimension inptl(mxptl)
      logical info

      info=.false.
      do n=1,mxptl
        inptl(n)=0
      enddo
      read(ifdt,'(a255)',end=1)line
 1    k=1
      nptevt=0
      do i=1,maxrec(1)
        l=irecty(i,1)
        if(l.eq.0)read(line(k:),'(i6)')ldummy !0
        if(l.eq.0) k=k+7
        if(l.eq.1)read(line(k:),'(i6)')ldummy !nrevt
        if(l.eq.1) k=k+7
        if(l.eq.2)read(line(k:),'(i6)')nptevt
        if(l.eq.2) k=k+7
        if(l.eq.3)read(line(k:),'(e12.6)')bimevt
        if(l.eq.3) k=k+13
        if(l.eq.4)read(line(k:),'(e12.6)')phievt
        if(l.eq.4) k=k+13
        if(l.eq.5)read(line(k:),'(i6)')kolevt
        if(l.eq.5) k=k+7
        if(l.eq.6)read(line(k:),'(e12.6)')pmxevt
        if(l.eq.6) k=k+13
        if(l.eq.7)read(line(k:),'(e12.6)')egyevt
        if(l.eq.7) k=k+13
        if(l.eq.8)read(line(k:),'(i6)')npjevt
        if(l.eq.8) k=k+7
        if(l.eq.9)read(line(k:),'(i6)')ntgevt
        if(l.eq.9) k=k+7
        if(l.eq.10)read(line(k:),'(i6)')npnevt
        if(l.eq.10) k=k+7
        if(l.eq.11)read(line(k:),'(i6)')nppevt
        if(l.eq.11) k=k+7
        if(l.eq.12)read(line(k:),'(i6)')ntnevt
        if(l.eq.12) k=k+7
        if(l.eq.13)read(line(k:),'(i6)')ntpevt
        if(l.eq.13) k=k+7
        if(l.eq.14)read(line(k:),'(i6)')jpnevt
        if(l.eq.14) k=k+7
        if(l.eq.15)read(line(k:),'(i6)')jppevt
        if(l.eq.15) k=k+7
        if(l.eq.16)read(line(k:),'(i6)')jtnevt
        if(l.eq.16) k=k+7
        if(l.eq.17)read(line(k:),'(i6)')jtpevt
        if(l.eq.17) k=k+7
        if(l.eq.20)read(line(k:),'(e12.6)')amproj
        if(l.eq.20) k=k+13
        if(l.eq.21)read(line(k:),'(e12.6)')amtarg
        if(l.eq.21) k=k+13
      enddo
      if(nptevt.eq.0)then
        print *,'sorry, '
        print *,'there is no particle number in the event record  '
        stop
      endif
      do n=1,nptevt
        read(ifdt,'(a255)',end=2)line
 2      k=1
        do i=1,maxrec(2)
          l=irecty(i,2)
          if(l.eq.0)read(line(k:),'(i6)') ldummy
          if(l.eq.0) k=k+7
          if(l.eq.1)then
            read(line(k:),'(i6)') nidp
            if(nidp.gt.0.and.nidp.le.mxptl)then
              info=.true.
              inptl(nidp)=n
            endif
            k=k+7
          endif
          if(l.eq.2)read(line(k:),'(i10)') idptl(n)
          if(l.eq.2) k=k+11
          if(l.eq.3.or.l.eq.17)read(line(k:),'(e12.6)') pptl(1,n)
          if(l.eq.3) k=k+13
          if(l.eq.4.or.l.eq.17)read(line(k:),'(e12.6)') pptl(2,n)
          if(l.eq.4) k=k+13
          if(l.eq.5.or.l.eq.17)read(line(k:),'(e12.6)') pptl(3,n)
          if(l.eq.5) k=k+13
          if(l.eq.6.or.l.eq.17)read(line(k:),'(e12.6)') pptl(4,n)
          if(l.eq.6) k=k+13
          if(l.eq.7.or.l.eq.17)read(line(k:),'(e12.6)') pptl(5,n)
          if(l.eq.7) k=k+13
          if(l.eq.8)then
            read(line(k:),'(i6)') iorptl(n)
            k=k+7
            if(info.and.iorptl(n).gt.0)iorptl(n)=inptl(iorptl(n))
          endif
          if(l.eq.9)then
            read(line(k:),'(i6)') jorptl(n)
            k=k+7
            if(info.and.jorptl(n).gt.0)jorptl(n)=inptl(jorptl(n))
          endif
          if(l.eq.10)read(line(k:),'(i6)') istptl(n)
          if(l.eq.10) k=k+7
          if(l.eq.11.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(1,n)
          if(l.eq.11) k=k+13
          if(l.eq.12.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(2,n)
          if(l.eq.12) k=k+13
          if(l.eq.13.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(3,n)
          if(l.eq.13) k=k+13
          if(l.eq.14.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(4,n)
          if(l.eq.14) k=k+13
c     if(i.eq.15)read(line(k:),'(i6)') idiptl(n)
          if(l.eq.15) k=k+7
c     if(i.eq.16)read(line(k:),'(i6)') idjptl(n)
          if(l.eq.16) k=k+7
          if(l.eq.19)read(line(k:),'(e12.6)') dezptl(n)
          if(l.eq.19) k=k+13
c          if(l.eq.21)read(line(k:),'(I6)') ifrptl(1,n)
          if(l.eq.21) k=k+7
c          if(l.eq.22)read(line(k:),'(I6)') ifrptl(2,n)
          if(l.eq.22) k=k+7
          if(l.eq.23)read(line(k:),'(I6)') ityptl(n)
          if(l.eq.23) k=k+7
        enddo
      enddo
      nptl=nptevt
      nevt=1
      end

c-----------------------------------------------------------------------
      subroutine aafinal
c-----------------------------------------------------------------------
c  * calculates xorptl(j,i), representing formation points.
c    (xorptl(j,i),j=1,4) is the 4-vector representing the space-time of
c    creation of particle i.
c-----------------------------------------------------------------------
      include 'epos.inc'
      do i=1,nptl
        if(idptl(i).ne.0.and.istptl(i).le.1)then
          if(    abs(tivptl(1,i)).le.ainfin
     .    .and.abs(xorptl(1,i)).le.ainfin
     .    .and.abs(xorptl(2,i)).le.ainfin
     .    .and.abs(xorptl(3,i)).le.ainfin
     .    .and.abs(xorptl(4,i)).le.ainfin
     .    .and.pptl(5,i).le.ainfin
     .    .and.pptl(4,i).gt.0.)then
            if(ish.ge.4)call alistc('afinal&',i,i)
            t=tivptl(1,i)
            xorptl(1,i)=xorptl(1,i)+pptl(1,i)/pptl(4,i)*(t-xorptl(4,i))
            xorptl(2,i)=xorptl(2,i)+pptl(2,i)/pptl(4,i)*(t-xorptl(4,i))
            xorptl(3,i)=xorptl(3,i)+pptl(3,i)/pptl(4,i)*(t-xorptl(4,i))
            xorptl(4,i)=t
          else
            if(ish.ge.1)then
              if(iorptl(i).gt.0)idior=idptl(iorptl(i))
              write(ifmt,'(a)')
     .        '*** warning (afinal see check file): '
              write(ifch,'(a,i6,i10,i10,i3,1x,7(e7.1,1x))')
     .        '*** warning (afinal): ',
     .        i,idptl(i),idior,ityptl(i),tivptl(1,i), pptl(4,i)
     .        ,pptl(5,i),xorptl(1,i),xorptl(2,i),xorptl(3,i),xorptl(4,i)
            endif
            tivptl(1,i)=2*ainfin
            xorptl(1,i)=2*ainfin
            xorptl(2,i)=2*ainfin
            xorptl(3,i)=2*ainfin
            xorptl(4,i)=2*ainfin
          endif
        endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine afinal
c-----------------------------------------------------------------------
c  does some final calculations, to be called before call aasto.
c  * calculates nptlu, the maximum nptl for all events.
c  * in case of mod(iframe,10) .ne. 1, these vectors are transformed
c    (being originally in the "natural frame",
c  * calculates numbers of spectators:
c    npnevt (number of primary proj neutron spectators)
c    nppevt (number of primary proj proton spectators)
c    ntnevt (number of primary targ neutron spectators)
c    ntpevt (number of primary targ proton spectators)
c    jpnevt (number of absolute proj neutron spectators)
c    jppevt (number of absolute proj proton spectators)
c    jtnevt (number of absolute targ neutron spectators)
c    jtpevt (number of absolute targ proton spectators)
c-----------------------------------------------------------------------

      include 'epos.inc'
      common/geom/rmproj,rmtarg,bmax,bkmx

      double precision pp1,pp2,pp3,pp4,pp5
      call utpri('afinal',ish,ishini,4)

      nptlu=max0(nptl,nptlu)

      if(mod(iframe,10).ne.1)then
      if(iframe.eq.12)then   !targ
      pp1=0d0
      pp2=0d0
      pp3=dsinh(dble(yhaha))
      pp4=dcosh(dble(yhaha))
      pp5=1d0
      else
      stop'transformation not yet defined'
      endif
      endif

      do i=1,nptl
        if(idptl(i).ne.0.and.istptl(i).le.1)then
          if(    abs(tivptl(1,i)).le.ainfin
     .    .and.abs(xorptl(1,i)).le.ainfin
     .    .and.abs(xorptl(2,i)).le.ainfin
     .    .and.abs(xorptl(3,i)).le.ainfin
     .    .and.abs(xorptl(4,i)).le.ainfin
     .    .and.pptl(5,i).le.ainfin
     .    .and.pptl(4,i).gt.0.)then
            if(mod(iframe,10).ne.1)then
              if(iframe.eq.12)then
                call utlob4(-1,pp1,pp2,pp3,pp4,pp5
     .          ,xorptl(1,i),xorptl(2,i),xorptl(3,i),xorptl(4,i))
                call utlob5(-yhaha
     .          ,pptl(1,i), pptl(2,i), pptl(3,i), pptl(4,i), pptl(5,i))
              else
                stop'transformation not yet defined'
              endif
            endif
          else
            tivptl(1,i)=ainfin
            xorptl(1,i)=ainfin
            xorptl(2,i)=ainfin
            xorptl(3,i)=ainfin
            xorptl(4,i)=ainfin
          endif
        endif
      enddo

      if(ish.ge.2)then
        if(model.eq.1)call alistf('EPOS&')
        if(model.eq.2)call alistf('QGSJET01&')
        if(model.eq.3)call alistf('GHEISHA&')
        if(model.eq.4)call alistf('PYTHIA&')
        if(model.eq.5)call alistf('HIJING&')
        if(model.eq.6)call alistf('SIBYLL 2.1&')
        if(model.eq.7)call alistf('QGSJET II&')
      endif

c      if(isto.eq.1)stop
c$$$      call testconex(2)

      npnevt=0
      nppevt=0
      ntnevt=0
      ntpevt=0
      jpnevt=0
      jppevt=0
      jtnevt=0
      jtpevt=0
      if(ish.ge.6)write(ifch,'(/31a1/a/31a1)')('-',l=1,31)
     *,'primary and absolute spectators',('-',l=1,31)
      if(ish.ge.6)write(ifch,'(/a//a/)')'projectile nucleons:'
     *,'     i    id   ior   ist'
      do i=1,maproj
      if(ish.ge.6)write(ifch,'(4i6)')i,idptl(i),iorptl(i),istptl(i)
      io=iorptl(i)
      id=idptl(i)
      is=istptl(i)
      if(io.eq.0.and.id.eq.1220)npnevt=npnevt+1
      if(io.eq.0.and.id.eq.1120)nppevt=nppevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1220)jpnevt=jpnevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1120)jppevt=jppevt+1
      enddo
      if(ish.ge.6)write(ifch,'(/a//a/)')'target nucleons:'
     *,'     i    id   ior   ist'
      do i=maproj+1,maproj+matarg
      if(ish.ge.6)write(ifch,'(4i6)')i,idptl(i),iorptl(i),istptl(i)
      io=iorptl(i)
      id=idptl(i)
      is=istptl(i)
      if(io.eq.0.and.id.eq.1220)ntnevt=ntnevt+1
      if(io.eq.0.and.id.eq.1120)ntpevt=ntpevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1220)jtnevt=jtnevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1120)jtpevt=jtpevt+1
      enddo
      if(ish.ge.6)then
      write(ifch,'(/a/)')'numbers of participants and spectators:'
      write(ifch,'(a,i4,a,i4)')'primary participants:   projectile:'
     *,npjevt,'   target:',ntgevt
      write(ifch,'(a,i4,a,i4)')'primary spectators:     projectile:'
     *,npnevt+nppevt,'   target:',ntnevt+ntpevt
      write(ifch,'(a,i4,a,i4)')
     *'primary spectator neutrons:   projectile:',npnevt
     *,'   target:',ntnevt
      write(ifch,'(a,i4,a,i4)')
     *'primary spectator protons:    projectile:',nppevt
     *,'   target:',ntpevt
      write(ifch,'(a,i4,a,i4)')'absolute spectators:    projectile:'
     *,jpnevt+jppevt,'   target:',jtnevt+jtpevt
      endif

      if(isigma.eq.1.and.ntevt.gt.0)then
        b1=bminim
        b2=min(bmax,bmaxim)
        a=pi*(b2**2-b1**2)
        sigineex=anintine/float(ntevt)*a*10
        sigdifex=anintdiff/float(ntevt)*a*10
        sigsdex=anintsdif/float(ntevt)*a*10
      endif


      if(ish.eq.27)stop'change this????????????' !call hnbcor(2)
      if(imihis.eq.1)call wimi
      if(imihis.eq.1.and.nrevt.eq.nevent)call wimino
      if(isphis.eq.1)call xspace(1)
      if(iclhis.eq.1)call wclu
      if(iclhis.eq.1.and.nrevt.eq.nevent)call wclufi
      if(iwtime.eq.1)call wtime(1)
      if(iwtime.eq.1.and.nrevt.eq.nevent)call wtime(2)

      if(ish.ge.8)call alistc('afinal&',1,nptl)

      call utprix('afinal',ish,ishini,4)
      return
      end

c-----------------------------------------------------------------------
      subroutine bfinal
c-----------------------------------------------------------------------
      include 'epos.inc'
      if(jerr(1).gt.4.and.jerr(1).gt.0.01*nevent)then
        write(ifch,'(3x,70a1)')('#',i=1,70)
        write(ifch,*)'  #   number of events:',nevent
        write(ifch,*)'  #   number of (flav > 9) warnings:',jerr(1)
        write(ifch,*)'  #        (OK when happens rarely)'
        write(ifch,'(3x,70a1)')('#',i=1,70)
      endif
      if(jerr(3).gt.4.and.jerr(3).gt.0.01*jerr(2))then
        write(ifch,'(3x,70a1)')('#',i=1,70)
        write(ifch,*)'  #   number of clusters:',jerr(2)
        write(ifch,*)'  #   number of neg m^2 clusters:',jerr(3)
        write(ifch,*)'  #          (OK when happens rarely)'
        write(ifch,'(3x,70a1)')('#',i=1,70)
      endif
      if(jerr(5).gt.4.and.jerr(5).gt.0.01*jerr(4))then
        write(ifch,'(3x,70a1)')('#',i=1,70)
        write(ifch,*)'  #   number of successful remnant cluster'
     &                       ,' decays:',jerr(4)
        write(ifch,*)'  #   number of unsuccessful remnant cluster'
     &                       ,' decays:',jerr(5)
        write(ifch,*)'  #        (OK when happens rarely)'
        write(ifch,'(3x,70a1)')('#',i=1,70)
      endif
      if(jerr(6).gt.4.and.jerr(6).gt.0.01*nevent)then
        write(ifch,'(3x,70a1)')('#',i=1,70)
        write(ifch,*)'  #   number of events:',nevent
        write(ifch,*)'  #   number of low mass remnant clusters:'
     &                      ,jerr(6)
        write(ifch,*)'  #        (OK when happens rarely)'
        write(ifch,'(3x,70a)')('#',i=1,70)
      endif
      end

c-----------------------------------------------------------------------
      subroutine ainit
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incsem'
      include 'epos.incpar'
      common/cquama/quama
      parameter (nptj=129)
      common /cptj/xptj(nptj),qptj(nptj),wptj(nptj)
      common/geom/rmproj,rmtarg,bmax,bkmx
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat,seedp
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctain/mtain
      double precision rcproj,rctarg
      common/geom1/rcproj,rctarg


      external sptj

      call utpri('ainit ',ish,ishini,4)

      inicnt=inicnt+1
      ntevt=0

      if(inicnt.eq.1)write(ifmt,'(a)')'initializations ...'

      if(seedi.ne.0d0)call ranfst(seedi)
      seedc=seedi
      call aseedi

      if(isphis.eq.1)iframe=11  !nncm
      if(icinpu.ge.1)elepti=engy
ctp060829      if(iopenu.eq.2)call smassi(themas)
      if(iopenu.eq.2.and.ish.eq.19)stop'change this?????????' !call smassp

      if(iappl.eq.5)then
      yhaha=0
      ypjtl=0
      endif

      if(ispherio.ne.0)ndecay=1
      if(ispherio.ne.0)idecay=0
      if(ispherio.ne.0)jdecay=0
      if(ifrade.eq.0)irescl=0
      idtarg=idtargin
      idproj=idprojin
      do 111 iii=1,4
 111  rexdif(iii)=abs(rexdifi(iii))
      if(noebin.gt.1)then
        engy=-1
        ekin=-1
        if(iologe.eq.1)engy=
     *       engmin*(engmax/engmin)**((real(nrebin)-0.5)/noebin)
        if(iologe.eq.0.or.(iologe.lt.0.and.iologl.lt.0))engy=
     *       engmin+(engmax-engmin)/noebin*(nrebin-0.5)
        if(iologl.eq.1)ekin=
     *       engmin*(engmax/engmin)**((real(nrebin-0.5))/noebin)
        if(iologl.eq.0)ekin=
     *       engmin+(engmax-engmin)/noebin*(real(nrebin)-0.5)
        elab=-1
        ecms=-1
        pnll=-1
        if(jpsi.lt.0)then
  11      z=0.19*sqrt(-2*alog(rangen()))*cos(2*pi*rangen())
          engy=abs(z)*engmax
          if(engy.lt.egymin)goto11
        endif
      endif

           if(iappl.le.3)then

        if(idtarg.eq.0)then   !in case of Air target, initialize with Argon nucleus
          if(model.eq.6)then     !no Argon in Sibyll
            latarg=7
            matarg=14
          else
            latarg=20
            matarg=40
          endif
        endif

        if((idproj.ne.1120.and.(laproj.ne.-1.or.maproj.ne.1))
     &    .or.maproj.le.0)
     &  call utstop('Invalid projectile setup !&')
c        if((idtarg.ne.1120.and.(latarg.ne.-1.or.matarg.ne.1))
c     &    .or.matarg.le.0)
c     &  call utstop('Invalid target setup !&')

      if(iabs(idtarg).ne.1120.and.iabs(idtarg).ne.1220.and.idtarg.ne.0)
     &  call utstop('Invalid target !&')
        if((((idtarg.eq.-1120.or.iabs(idtarg).eq.1220)
     &    .and.(latarg.ne.-1.or.matarg.ne.1))
     &    .and.(idtarg.ne.1120.or.latarg.lt.0))
     &    .or.matarg.le.0)
     &  call utstop('Invalid target setup !&')


      call idmass(idproj,amproj)
      call idmass(idtarg,amtarg)
      nre=0
      if(engy.ge.0.)nre=nre+1
      if(pnll.ge.0.)nre=nre+1
      if(elab.ge.0.)nre=nre+1
      if(ekin.ge.0.)nre=nre+1
      if(ecms.ge.0.)nre=nre+1
      if(nre.ne.1)stop'invalid energy definition'
      ifirstghe=0
 101  continue
         if(engy.gt.0.)then
      pnll=sqrt(max(0., ((engy**2-amproj**2-amtarg**2)/2/amtarg)**2
     &                   -amproj**2) )
      elab=sqrt(pnll**2+amproj**2)
      ekin=elab-amproj
      ecms=engy
         elseif(ecms.gt.0.)then
      engy=ecms
      pnll=sqrt(max(0., ((engy**2-amproj**2-amtarg**2)/2/amtarg)**2
     &                   -amproj**2) )
      elab=sqrt(pnll**2+amproj**2)
      ekin=elab-amproj
         elseif(elab.gt.0)then
      pnll=sqrt(max(0.,elab**2-amproj**2))
      engy=sqrt( 2*elab*amtarg+amtarg**2+amproj**2 )
      ecms=engy
      ekin=elab-amproj
         elseif(pnll.gt.0)then
      elab=sqrt(pnll**2+amproj**2)
      engy=sqrt( 2*sqrt(pnll**2+amproj**2)*amtarg+amtarg**2+amproj**2 )
      ecms=engy
      ekin=elab-amproj
         elseif(ekin.gt.0.)then
      elab=ekin+amproj
      pnll=sqrt(max(0.,elab**2-amproj**2))
      engy=sqrt( 2*elab*amtarg+amtarg**2+amproj**2 )
      ecms=engy
         endif

         if(model.eq.3.and.ifirstghe.eq.0)then    !det, trit and alp
           if(maproj.eq.2.and.laproj.eq.1)idproj=17
           if(maproj.eq.3.and.laproj.eq.1)idproj=18
           if(maproj.eq.4.and.laproj.eq.2)idproj=19
           if(idproj.ge.17.and.idproj.le.19)then
             elab=elab*maproj
             call idmass(idproj,amproj)
             maproj=1
             laproj=-1
             ifirstghe=1
             engy=-1
             ecms=-1
             pnll=-1
             ekin=-1
             goto 101
           endif
         endif

      if(pnll.le.0.001)call utstop('ainit: energy too low&')
      if(engy.gt.egymax)call utstop('ainit: energy too high&')
      s=engy**2
      pnullx=utpcm(engy,amproj,amtarg)
      yhaha=alog((sqrt(pnll**2+s)+pnll)/sqrt(s))
      ypjtl=alog((sqrt(pnll**2+amproj**2)+pnll)/amproj)

         elseif(iappl.eq.7)then

      call idmass(idproj,amproj)
         if(elab.gt.0)then
      pnll=sqrt(max(0.,elab**2-amproj**2))
      engy=amproj
      ecms=engy
      ekin=elab-amproj
         elseif(pnll.gt.0)then
      elab=sqrt(pnll**2+amproj**2)
      engy=amproj
      ecms=engy
      ekin=elab-amproj
         elseif(ekin.gt.0.)then
      elab=ekin+amproj
      pnll=sqrt(max(0.,elab**2-amproj**2))
      engy=amproj
      ecms=engy
         else
      engy=amproj
      ecms=amproj
      elab=0.
      pnll=0.
      ekin=0.
        endif

      pnullx=0.
      ypjtl=alog((sqrt(pnll**2+amproj**2)+pnll)/amproj)
      yhaha=ypjtl

         endif

      detap=(ypjtl-yhaha)*etafac
      detat=-yhaha*etafac
      tpro=dcosh(detap)
      zpro=dsinh(detap)
      ttar=dcosh(detat)
      ztar=dsinh(detat)

      egyevt=engy
      pmxevt=pnll

      if(iappl.gt.8)stop'update following statement'
      if(iappl.ge.5.and.iappl.le.8)then
      s=12.**2
      endif

      if(iappl.le.3)then
       if(maproj.gt.1)then
        rpj=1.19*maproj**(1./3.)-1.61*maproj**(-1./3.)
        rmproj=rpj+fctrmx*.54
        rcproj=dble(rpj/cosh(yhaha)*facnuc)
       else
        rmproj=0
        rcproj=dble(0.8/cosh(yhaha)*facnuc)
       endif
       if(matarg.gt.1)then
        rtg=1.19*matarg**(1./3.)-1.61*matarg**(-1./3.)
        rmtarg=rtg+fctrmx*.54
        rctarg=dble(rtg/cosh(yhaha)*facnuc)
       else
        rmtarg=0
        rctarg=dble(0.8/cosh(yhaha)*facnuc)
       endif

      endif

      call iclass(idproj,iclpro)
      call iclass(idtarg,icltar)

         if(inicnt.eq.1)then

          call ranfgt(seedp)                     !not to change the seed ...

      call hdecin(.false.)

      if(iappl.eq.1.or.iappl.ge.5)then
      c=6
      call utquaf(sptj,nptj,xptj,qptj,wptj,0.,.33*c,.66*c,c)
      endif

      if(iappl.ne.2)then
        call hnbspd(iospec)
        ktnbod=0
c        if(model.eq.1)call hnbxxxini
        call hnbpajini
      endif

      if(model.eq.1)then
        if(iclegy2.gt.1)then
          egyfac=(egymax*1.0001/egylow)**(1./float(iclegy2-1))
        else
          egyfac=1.
        endif
      endif

      if(model.eq.1)then
        call conini
      else
        iorsce=0
        iorsdf=0
        iorshh=0
        iorsdf=0
        irescl=0       !rescaling
      endif

      if(iappl.ne.6.and.model.eq.1)call psaini

          call ranfst(seedp)                     ! ... after this initialization

         endif

c$$$      if(idproj.eq.1120)icp=2        !????????????? for what ?
c$$$      if(idproj.eq.-1120)icp=-2
c$$$      if(idproj.eq.120)icp=1
c$$$      if(idproj.eq.-120)icp=-1
c$$$      if(idproj.eq.130)icp=4
c$$$      if(idproj.eq.-130)icp=-4



      if(model.eq.1)then                   !only for epos


      if(iappl.le.3)then
        call emsini(engy,idproj,idtarg)
        call paramini(1)
        if(ish.ge.4)then
          do i=idxD0,idxD1
            write(ifch,'(9(a,f8.4))')
     *      'AlpD:',alpD(i,iclpro,icltar)
     * ,'    AlpDp:',alpDp(i,iclpro,icltar)
     * ,'    AlpDpp:',alpDpp(i,iclpro,icltar)
     * ,'    BetD:',betD(i,iclpro,icltar)
     * ,'    BetDp:',betDp(i,iclpro,icltar)
     * ,'    BetDpp:',betDpp(i,iclpro,icltar)
     * ,'    GamD:',gamD(i,iclpro,icltar)
     * ,'    DelD:',delD(i,iclpro,icltar)
     * ,'    AlpPar:',alppar
     * ,'    bkmxdif:',bmxdif(iclpro,icltar)
          enddo
        endif
      endif

      if(iappl.le.3)then
        bkmxndif=conbmxndif()
        bkmx=conbmx()
        if(ish.ge.3)write(ifch,*)'bkmx,bkmxndif',bkmx,bkmxndif

        if(maproj.gt.1.or.matarg.gt.1)then
        bmax=rmproj+rmtarg
        else
        bmax=bkmx
        endif
      endif

      if(ixtau.eq.1)call xtauev(0)

      if(model.eq.1)then
      if(iEmsB.eq.1)call xEmsB(0,0,0)
      if(iEmsBg.eq.1)call xEmsBg(0,0,0)
      if(iEmsPm.eq.1)call xEmsPm(0,0,0)
      if(iEmsPx.eq.1)call xEmsPx(0,0.,0.,0)
      if(iEmsPBx.eq.1)call xEmsP2(0,0,0,0.,0.,0.,0.,0.,0.)
c      if(iEmsPx.eq.1)call xEmsPxNo(0,0.,0.,0,0)
      if(iEmsSe.eq.1)call xEmsSe(0,0.,0.,0,1)
      if(iEmsSe.eq.1)call xEmsSe(0,0.,0.,0,2)
      if(iEmsDr.eq.1)call xEmsDr(0,0.,0.,0)
      if(iEmsRx.eq.1)call xEmsRx(0,0,0.,0.)
      endif

c G function parameters

      if(iappl.eq.1)then
        call GfunPar(1,1,0.,s,alp,bet,betp,epsp,epst,epss,gamvv)  !b=0
        epszero=epss
        do i=1,nclha
         alpff(i)=   engy**epszero*gamhad(i)
        enddo
        betff(1)=   -alppar+epsp
        betff(2)=   -alppar+epst
      endif

      endif

c additional initialization procedures

      if(model.ne.1)then
        call IniEvtModel
      elseif(iappl.le.3)then
c Cross section calculation
        call xsigma
      endif


      if(idtarg.eq.0)idtarg=1120 !air = nucleus



      if(seedj.ne.0.)call ranfst(seedj)
      seedc=seedj


ccc      call MakeFpartonTable

c$$$      call testconex(1)

      call utprix('ainit ',ish,ishini,4)
      return
      end

c---------------------------------------------------------------------
      subroutine aread
c---------------------------------------------------------------------
c  reads and interprets input commands
c---------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incpar'
      include 'epos.incsem'

      double precision histoweight
      common/chiswei/histoweight
      common/cyield/yield/cifset/ifset/caverg/averg
      common/csigma/sigma
      double precision val,val1,val2
      character*160 line,linex,cline
      data nappl /0/
      common/record/maxrec(2),irecty(30,2)
      common/cfacmss/facmss /cr3pomi/r3pomi,r4pomi
      common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon
      character*80 fndat,fnncs,fnIIdat,fnIIncs                 !qgs-II????????
      common/qgsfname/  fndat, fnncs, ifdat, ifncs
      common/qgsIIfname/fnIIdat, fnIIncs, ifIIdat, ifIIncs     !qgs-II????????
      common/qgsnfname/ nfndat, nfnncs
      common/qgsIInfname/ nfnIIdat, nfnIIncs     !qgs-II????????
      common/ghecsquel/anquasiel,iquasiel
      common/cjjj/jjj,cline
      character cmodel*21
      parameter(mxdefine=40)
      character w1define*100,w2define*100
      common/cdefine/ndefine,l1define(mxdefine),l2define(mxdefine)
     &               ,w1define(mxdefine),w2define(mxdefine)
      common/cbincond/nozero,ibmin,ibmax

	  j=-1
      nhsto=0
      ndefine=0

    1 call utword(line,i,j,1)

          if(line(i:j).eq.'#define')then

      call utword(line,i,j,ne)
      if(line(i:j).eq.'bim3')stop'****** bim3->bim03 ****** '
      if(line(i:j).eq.'bim5')stop'****** bim5->bim05 ****** '
      if(line(i:j).eq.'bim6')stop'****** bim6->bim06 ****** '
      if(ndefine+1.gt.mxdefine)stop'too many `define` statements.      '
      l1=j-i+1
      if(l1.gt.99)stop'`define` argument 1 too long.            '
      w1define(ndefine+1)(1:l1)=line(i:j)
      w1define(ndefine+1)(l1+1:l1+1)=' '
      call utword(line,i,j,ne)
      l2=j-i+1
      if(l2.gt.99)stop'`define` argument 2 too long.            '
      w2define(ndefine+1)(1:l2)=line(i:j)
      w2define(ndefine+1)(l2+1:l2+1)=' '
      ndefine=ndefine+1
      l1define(ndefine)=l1
      l2define(ndefine)=l2

          elseif(line(i:j).eq.'goto')then

      call utword(line,i,j,ne)
      ix=i
      jx=j
      linex=line
      call utword(line,i,j,ne)
      do while(line(i:j).ne.linex(ix:jx))
      call utword(line,i,j,ne)
      enddo
      goto1

           elseif(line(i:j).eq.'application')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'application?'
      call utword(line,i,j,0)
      if(nopen.ne.-1)then       !only first read
      if(line(i:j).eq.'analysis')  iappl=0
      if(line(i:j).eq.'hadron')    iappl=1
      if(line(i:j).eq.'geometry')  iappl=2
      if(line(i:j).eq.'read')      iappl=3
      if(line(i:j).eq.'micro')     iappl=4
      if(line(i:j).eq.'kinky')     iappl=5
      if(line(i:j).eq.'ee')        iappl=6
      if(line(i:j).eq.'decay')     iappl=7
      if(line(i:j).eq.'lepton')    iappl=8
      if(line(i:j).eq.'ee')    then
        naflav=5                ! number of flavors
      endif
      nappl=nappl+1
      if(iappl.gt.0.and.nappl.gt.1)call aaset(1)
      if(iappl.eq.0.and.nappl.gt.1)call aaset(2)
      if(iappl.eq.0)jframe=iframe
      if(iappl.eq.0)kframe=iframe
      if(iappl.eq.1)iframe=0
      if(iappl.eq.2)iframe=0
      if(iappl.eq.4)iframe=1
      if(iappl.eq.5)iframe=1
      if(iappl.eq.6)iframe=1
      if(iappl.eq.7)iframe=0
      if(iappl.eq.8)iframe=21        !gncm
      endif

           elseif(line(i:j).eq.'call')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'subroutine?'
      call utword(line,i,j,0)

      if(nopen.eq.-1)then       !-----only second run

      if(line(i:j).eq.'xConThickProj')call xConThick(1)
      if(line(i:j).eq.'xConThickTarg')call xConThick(2)
      if(line(i:j).eq.'xConNuclDensProj')call xConNuclDens(1)
      if(line(i:j).eq.'xConNuclDensTarg')call xConNuclDens(2)
      if(line(i:j).eq.'xConNuclDensProjTarg')call xConNuclDens(1)
      if(line(i:j).eq.'xConNuclDensProjTarg')call xConNuclDens(2)
      if(line(i:j).eq.'xFom')call xfom
      if(line(i:j).eq.'xGeometry')call xGeometry(2)
      if(line(i:j).eq.'xEpsilon')call xEpsilon(2)
      if(line(i:j).eq.'xRanPt')call xRanPt
      if(line(i:j).eq.'xParGam')call xParGam
      if(line(i:j).eq.'xParGampp')call xParGampp
      if(line(i:j).eq.'xParOmega1xy')call xParOmega1xy
c$$$      if(line(i:j).eq.'xParOmega3xyz')call xParOmega3xyz
      if(line(i:j).eq.'xParPro')call xParPro
      if(line(i:j).eq.'xParPro1')call xParPro1
      if(line(i:j).eq.'xParPomInc')call xParPomInc
      if(line(i:j).eq.'xParPomIncX')call xParPomIncX
      if(line(i:j).eq.'xParPomIncP')call xParPomIncP
      if(line(i:j).eq.'xParPomIncM')call xParPomIncM
      if(line(i:j).eq.'xParPomIncXI')call xParPomIncXI
      if(line(i:j).eq.'xParPomIncPI')call xParPomIncPI
      if(line(i:j).eq.'xParPomIncMI')call xParPomIncMI
      if(line(i:j).eq.'xParPomIncJ')call xParPomIncJ
      if(line(i:j).eq.'xParOmega1')call xParOmega1
c$$$      if(line(i:j).eq.'xParOmega3')call xParOmega3
c$$$      if(line(i:j).eq.'xParOmega5')call xParOmega5
      if(line(i:j).eq.'xParOmegaN')call xParOmegaN
      if(line(i:j).eq.'xParGauss')call xParGauss
      if(line(i:j).eq.'xParSigma')call xParSigma
c$$$      if(line(i:j).eq.'xParSigma2')call xParSigma2
c$$$      if(line(i:j).eq.'xScrD')call xScrD
      if(line(i:j).eq.'xFitD1')call xFitD1
c$$$      if(line(i:j).eq.'xExaD2')call xExaD2
      if(line(i:j).eq.'xbExaD')call xbExaD
c$$$      if(line(i:j).eq.'xbExaD2')call xbExaD2
      if(line(i:j).eq.'xbnExaD')call xbnExaD
c$$$      if(line(i:j).eq.'xbnExaD2')call xbnExaD2
      if(line(i:j).eq.'xFitD2')call xFitD2
      if(line(i:j).eq.'xbParD')call xbParD
c$$$      if(line(i:j).eq.'xParD2')call xParD2
      if(line(i:j).eq.'xGexaJ')call xGexaJ
      if(line(i:j).eq.'xbnParD')call xbnParD
      if(line(i:j).eq.'xsParD')call xsParD
c$$$      if(line(i:j).eq.'xmParD2')call xmParD2
      if(line(i:j).eq.'xyParD')call xyParD
c$$$      if(line(i:j).eq.'xyParD2')call xyParD2
      if(line(i:j).eq.'xParPhi1')call xParPhi1
      if(line(i:j).eq.'xParPhi')call xParPhi
      if(line(i:j).eq.'xParH')call xParH
      if(line(i:j).eq.'xParHPhiInt')call xParHPhiInt
      if(line(i:j).eq.'xParZ')call xParZ
      if(line(i:j).eq.'xtauev')call xtauev(2)
      if(line(i:j).eq.'xspace')call xspace(2)
      if(line(i:j).eq.'gakjto'   )call gakjto
      if(line(i:j).eq.'psaevp')call psaevp
c     if(line(i:j).eq.'pyarea')call pyarea
      if(line(i:j).eq.'xhgcam')call xhgcam(0.,-1)
      if(line(i:j).eq.'xhgcfl')call xhgcfl(0.,0.,0.,-1)
      if(line(i:j).eq.'xhgccc')call xhgccc(-1.)
      if(line(i:j).eq.'xjden1')call xjden1(2,0,0.,0.,0.,0.,0.)
      if(line(i:j).eq.'xjden2')call xjden2(2,0,0.,0.,0.,0.)
c     if(line(i:j).eq.'xjdis' )call xjdis(2,0,0)
      if(line(i:j).eq.'xhgcen')call xhgcen
      if(line(i:j).eq.'xhgcmt')call xhgcmt
      if(line(i:j).eq.'xhgcmu')call xhgcmu
      if(line(i:j).eq.'xhgcmx')call xhgcmx
      if(line(i:j).eq.'xhgcpt')call xhgcpt
      if(line(i:j).eq.'xhgcra')call xhgcra
      if(line(i:j).eq.'xhnben')call xhnben
      if(line(i:j).eq.'xhnbit')call xhnbit
      if(line(i:j).eq.'xhnbmu')call xhnbmu
      if(line(i:j).eq.'xhnbmz')call xhnbmz
      if(line(i:j).eq.'xhnbte')call xhnbte(-1)
      if(line(i:j).eq.'xhnbti')call xhnbti(-1)
      if(model.eq.1)then
      if(line(i:j).eq.'xEmsB' )call xEmsB(2,0,0)
      if(line(i:j).eq.'xEmsBg')call xEmsBg(2,0,0)
      if(line(i:j).eq.'xEmsPm')call xEmsPm(2,0,0)
      if(line(i:j).eq.'xEmsPx')call xEmsPx(2,0.,0.,0)
c      if(line(i:j).eq.'xEmsPx')call xEmsPxNo(2,0.,0.,0,0)
      if(line(i:j).eq.'xEmsSe')call xEmsSe(2,0.,0.,0,1)
      if(line(i:j).eq.'xEmsSe')call xEmsSe(2,0.,0.,0,2)
      if(line(i:j).eq.'xEmsDr')call xEmsDr(2,0.,0.,0)
      if(line(i:j).eq.'xEmsRx')call xEmsRx(2,0,0.,0.)
      if(line(i:j-4).eq.'xEmsP2')then
        read(line(j-1:j-1),*)val
        idh=nint(val)
        read(line(j:j),*)val
        jex=nint(val)
        if(line(j-3:j-2).eq.'PE')
     &       call xEmsP2(2,idh,jex,0.,0.,0.,0.,0.,0.)
        if(line(j-3:j-2).eq.'IB')
     &       call xEmsP2(3,idh,jex,0.,0.,0.,0.,0.,0.)
        if(line(j-3:j-2).eq.'OB')
     &       call xEmsP2(4,idh,jex,0.,0.,0.,0.,0.,0.)
      endif
      if(line(i:j).eq.'xConxyzProj')
     &stop'xConxyzProj->xConNuclDensProj'
      if(line(i:j).eq.'xConxyzTarg')
     &stop'xConxyzTarg->xConNuclDensTarg'
      if(line(i:j).eq.'xConxyzProjTarg')
     &stop'xConxyzProjTarg->xConNuclDensProjTarg'

      endif

      elseif(model.eq.1)then  !first run and epos

      if(line(i:j).eq.'xGeometry')then
       call xGeometry(0)
       ixgeometry=1
      elseif(line(i:j).eq.'xEpsilon')then
       call xEpsilon(0)
      elseif(line(i:j).eq.'xtauev')then
       ixtau=1
      elseif(line(i:j).eq.'xEmsB')then
       iEmsB=1
      elseif(line(i:j).eq.'xEmsBg')then
       iEmsBg=1
      elseif(line(i:j).eq.'xEmsPm')then
       iEmsPm=1
      elseif(line(i:j).eq.'xEmsPx')then
       iEmsPx=1
      elseif(line(i:j-4).eq.'xEmsP2')then
       iEmsPBx=1
      elseif(line(i:j).eq.'xEmsSe')then
       iEmsSe=1
      elseif(line(i:j).eq.'xEmsDr')then
       iEmsDr=1
      elseif(line(i:j).eq.'xEmsRx')then
       iEmsRx=1
      elseif(line(i:j).eq.'xEmsI1')then
       iEmsI1=1
       if(iEmsI1+iEmsI2.eq.1)write(ifhi,'(a)')'newpage zone 3 4 1'
      elseif(line(i:j).eq.'xEmsI2')then
       iEmsI2=1
       if(iEmsI1+iEmsI2.eq.1)write(ifhi,'(a)')'newpage zone 3 4 1'
      elseif(line(i:j).eq.'xSpaceTime')then
       iSpaceTime=1
      elseif(line(i:j).eq.'xxSpaceTime')then
       stop'xxSpaceTime->xSpaceTime.              '
      endif
      endif

           elseif(line(i:j).eq.'decayall')then

      nrnody=0

           elseif(line(i:j).eq.'echo')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'on or off?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'on')iecho=1
      if(line(i:j).eq.'off')iecho=0
      if(line(i:j).ne.'on'.and.line(i:j).ne.'off')stop'invalid option'
      if(nopen.eq.-1)iecho=0

           elseif(line(i:j).eq.'fqgsjet')then              !QGSJet

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'dat')fndat(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dat')nfndat=j-i+1             !length of qgsdat01-file name
      if(linex(ix:jx).eq.'ncs')fnncs(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'ncs')nfnncs=j-i+1             !length of sectnu-file name
      if(nfndat.gt.1)ifdat=1
      if(nfnncs.gt.1)ifncs=2

           elseif(line(i:j).eq.'fqgsjetII')then              !QGSJET-II     !qgs-II????????

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'dat')fnIIdat(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dat')nfnIIdat=j-i+1             !length of qgsjet-II.dat name
      if(linex(ix:jx).eq.'ncs')fnIIncs(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'ncs')nfnIIncs=j-i+1             !length of qgsjet-II.ncs name
      if(nfnIIdat.gt.1)ifIIdat=1
      if(nfnIIncs.gt.1)ifIIncs=2

           elseif(line(i:j).eq.'fname')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'pathnx')fnnx(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'check')fnch(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'histo')fnhi(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'data') fndt(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'copy') fncp(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'initl') fnii(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'initl+') fnii(nfnii+1:nfnii+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidi') fnid(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidi+') fnid(nfnid+1:nfnid+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidr') fndr(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidr+') fndr(nfndr+1:nfndr+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'iniev') fnie(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'iniev+') fnie(nfnie+1:nfnie+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inirj') fnrj(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inirj+') fnrj(nfnrj+1:nfnrj+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inics') fncs(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inics+') fncs(nfncs+1:nfncs+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inihy') fnhy(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'pathnx')nfnnx=j-i+1
      if(linex(ix:jx).eq.'check')nfnch=j-i+1
      if(linex(ix:jx).eq.'histo')nfnhi=j-i+1
      if(linex(ix:jx).eq.'data') nfndt=j-i+1
      if(linex(ix:jx).eq.'copy') nfncp=j-i+1
      if(linex(ix:jx).eq.'initl') nfnii=j-i+1
      if(linex(ix:jx).eq.'initl+')nfnii=nfnii+j-i+1
      if(linex(ix:jx).eq.'inidi') nfnid=j-i+1
      if(linex(ix:jx).eq.'inidi+')nfnid=nfnid+j-i+1
      if(linex(ix:jx).eq.'inidr') nfndr=j-i+1
      if(linex(ix:jx).eq.'inidr+')nfndr=nfndr+j-i+1
      if(linex(ix:jx).eq.'iniev') nfnie=j-i+1
      if(linex(ix:jx).eq.'iniev+')nfnie=nfnie+j-i+1
      if(linex(ix:jx).eq.'inirj') nfnrj=j-i+1
      if(linex(ix:jx).eq.'inirj+')nfnrj=nfnrj+j-i+1
      if(linex(ix:jx).eq.'inics') nfncs=j-i+1
      if(linex(ix:jx).eq.'inihy') nfnhy=j-i+1
      if(linex(ix:jx).eq.'inics+')nfncs=nfncs+j-i+1
      if(linex(ix:jx).eq.'check'.and.fnch(1:nfnch).ne.'none') then
      open(unit=ifcx,file=fnch(1:nfnch),status='unknown')
      kchopen=1
      elseif(linex(ix:jx).eq.'pathnx'.and.fnnx(1:nfnnx).ne.'none')then
      knxopen=1
      elseif(linex(ix:jx).eq.'histo'.and.fnhi(1:nfnhi).ne.'none')then
      open(unit=ifhi,file=fnhi(1:nfnhi),status='unknown')
      khiopen=1
      elseif(linex(ix:jx).eq.'data'.and.fndt(1:nfndt).ne.'none')then
      open(unit=ifdt,file=fndt(1:nfndt),status='unknown')
      kdtopen=1
      write(6,'(a,a)')'Opened data file ',fndt(1:nfndt)
      elseif(linex(ix:jx).eq.'copy'.and.fncp(1:nfncp).ne.'none')then
      open(unit=ifcp,file=fncp(1:nfncp),status='unknown')
      kcpopen=1
      endif

           elseif(line(i:j).eq.'frame')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'frame?'
      call utword(line,i,j,0)
        if(nopen.ne.-1)then       ! event definition only, not analysis
        if(line(i:j).eq.'nucleon-nucleon')then
      if(iappl.eq.0)jframe=11
      if(iappl.eq.1)iframe=11
      if(iappl.eq.3)iframe=11
      if(iappl.gt.3.or.iappl.eq.2)stop'invalid option nucleon-nucleon'
        elseif(line(i:j).eq.'target')then
      if(iappl.eq.0)jframe=12
      if(iappl.eq.1)iframe=12
      if(iappl.eq.3)iframe=12
      if(iappl.gt.3.or.iappl.eq.2)stop'invalid option target'
        elseif(line(i:j).eq.'gamma-nucleon')then
      if(iappl.eq.0)jframe=21
      if((iappl+1)/2.eq.4)iframe=21
      if(iappl.eq.3)iframe=21
      if(iappl.ne.0.and.(iappl+1)/2.ne.4.and.iappl.ne.3)
     *stop'invalid option gamma-nucleon'
        elseif(line(i:j).eq.'lab')then
      if(iappl.eq.0)jframe=22
      if((iappl+1)/2.eq.4)iframe=22
      if(iappl.eq.3)iframe=22
      if(iappl.ne.0.and.(iappl+1)/2.ne.4.and.iappl.ne.3)
     *stop'invalid option lab'
        elseif(line(i:j).eq.'sphericity')then
      if(iappl.eq.0)jframe=32
      if(iappl.ne.0)stop'invalid option sphericity'
        elseif(line(i:j).eq.'thrust')then
      if(iappl.eq.0)jframe=33
      if(iappl.ne.0)stop'invalid option thrust'
        elseif(line(i:j).eq.'breit')then
          if(iappl.ne.0)stop'invalid option breit'
        else
      stop'frame not recognized'
        endif
        endif

           elseif(line(i:j).eq.'frame+')then

      call utword(line,i,j,0)

           elseif(line(i:j).eq.'binning')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'log/lin?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'lin')iologb=0
      if(line(i:j).eq.'log')iologb=1

           elseif(line(i:j).eq.'beginhisto')then

      if(nopen.ne.-1)then
        jjj=j
        cline=line
        call xini
        j=jjj
        line=cline
      endif

           elseif(line(i:j).eq.'endhisto')then

      if(nopen.eq.-1)then
        nhsto=nhsto+1
        call xhis(nhsto)
      endif

           elseif(line(i:j).eq.'noweak')then

      !only used in xini

           elseif(line(i:j).eq.'histogram')then

      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)

           elseif(line(i:j).eq.'plot')then

      stop' "plot" not used any more'

           elseif(line(i:j).eq.'root')then

      stop' "root" not used any more'

           elseif(line(i:j).eq.'idcode')then

      call utword(line,i,j,0)

           elseif(line(i:j).eq.'xpara')then

      call utword(line,i,j,0)
      call utword(line,i,j,0)

           elseif(line(i:j).eq.'histoweight')then

      if(nopen.eq.-1)then
      write(ifhi,'(a,d22.14)')'histoweight ',histoweight
      endif

           elseif(line(i:j).eq.'input')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'input file?'
      call utword(line,i,j,0)
      if(nopen.ge.0)then
       nopen=nopen+1
       if(nopen.gt.9)stop'too many nested input commands'
       open(unit=20+nopen,file=line(i:j),status='old')
       if(iprmpt.eq.1)iprmpt=-1
      endif

           elseif(line(i:j).eq.'nodecays')then

      call utword(line,i,j,0)
      do while (line(i:j).ne.'end')
       if(nrnody.ge.mxnody)then
        write(ifmt,'(a)')'too many nodecays; command ignored'
       else
        nrnody=nrnody+1
        read(line(i:j),*)val
        nody(nrnody)=nint(val)
       endif
      call utword(line,i,j,0)
      enddo

           elseif(line(i:j).eq.'nodecay')then

      call utword(line,i,j,0)
      if(nopen.ne.-1)then
      if(nrnody.ge.mxnody)then
      write(ifmt,'(a)')'too many nodecay commands; command ignored'
      j=160
      i=j+1
      goto1
      endif
      nrnody=nrnody+1
      read(line(i:j),*)val
      nody(nrnody)=nint(val)
      endif

           elseif(line(i:j).eq.'dodecay')then

      call utword(line,i,j,0)
      if(nopen.ne.-1)then
      read(line(i:j),*)val
      idx=nint(val)
      nrn=0
      imv=0
      do while(nrn.lt.nrnody)
       nrn=nrn+1
       if(idx.eq.nody(nrn))then
         nrnody=nrnody-1
         imv=1
       endif
       if(imv.eq.1.and.nrn.le.nrnody)nody(nrn)=nody(nrn+1)
      enddo
      endif

           elseif(line(i:j).eq.'print')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'subroutine?'
      call utword(line,i,j,0)
       if(line(i:j).eq.'*')then
      nrpri=0
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'level?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      ish=nint(val)
       else
      nrpri=nrpri+1
      subpri(nrpri)='                    '
      subpri(nrpri)(1:j-i+1)=line(i:j)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'level?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      ishpri(nrpri)=nint(val)
       endif

           elseif(line(i:j).eq.'printcheck')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'screen or file?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'screen')ifch=ifmt
      if(line(i:j).eq.'file')ifch=ifcx
      if(line(i:j).ne.'screen'.and.line(i:j).ne.'file')
     *write(ifmt,'(a)')'invalid option; command ignored'

           elseif(line(i:j).eq.'prompt')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'on or off or auto?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'on')iprmpt=2
      if(line(i:j).eq.'off')iprmpt=-2
      if(line(i:j).eq.'auto'.and.nopen.eq.0)iprmpt=1
      if(line(i:j).eq.'auto'.and.nopen.eq.1)iprmpt=-1
      if(line(i:j).ne.'on'.and.line(i:j).ne.'off'
     *.and.line(i:j).ne.'auto')stop'invalid option'

           elseif(line(i:j).eq.'return')then

      if(nopen.ne.-1)then
        close(ifop)
        nopen=nopen-1
        if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
ccc      close(20+nopen)
ccc      nopen=nopen-1
ccc      if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
      endif

           elseif(line(i:j).eq.'run')then

           elseif(line(i:j).eq.'runprogram')then

      if(kchopen.eq.0.and.fnch(1:nfnch).ne.'none')then
        open(unit=ifcx,file=fnch(1:nfnch),status='unknown')
        kchopen=1
      endif
      if(khiopen.eq.0.and.fnhi(1:nfnhi).ne.'none')then
        open(unit=ifhi,file=fnhi(1:nfnhi),status='unknown')
        khiopen=1
      endif
      if(kdtopen.eq.0.and.fndt(1:nfndt).ne.'none')then
        open(unit=ifdt,file=fndt(1:nfndt),status='unknown')
        kdtopen=1
        write(6,'(a,a)')'Opened data file ',fndt(1:nfndt)
      endif
      if(kcpopen.eq.0.and.fncp(1:nfncp).ne.'none')then
        open(unit=ifcp,file=fncp(1:nfncp),status='unknown')
        kcpopen=1
      endif
      return

           elseif(line(i:j).eq.'if')then

      call utword(line,i,j,0)
      ix=i
      jx=j
      linex=line
      call utword(line,i,j,0)
      read(line(i:j),*)val1
      call utword(line,i,j,0)
      read(line(i:j),*)val2
      ifset=1
      if(linex(ix:jx).eq.'engy')then
        call idmass(idproj,amproj)
        call idmass(idtarg,amtarg)
        if(engy.gt.0.)then
          xxengy=engy
        elseif(ecms.gt.0.)then
          xxengy=ecms
        elseif(elab.gt.0)then
          xxengy=sqrt( 2*elab*amtarg+amtarg**2+amproj**2 )
        elseif(pnll.gt.0)then
          xxengy=sqrt( 2*sqrt(pnll**2+amproj**2)*amtarg
     *                   +amtarg**2+amproj**2 )
        elseif(ekin.gt.0.)then
          xxelab=ekin+amproj
          xxengy=sqrt( 2*xxelab*amtarg+amtarg**2+amproj**2 )
        endif
        if(xxengy.lt.val1.or.xxengy.gt.val2)ifset=0
      endif

           elseif(line(i:j).eq.'set'.and.ifset.eq.1)then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'p/o-name p/o-value?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'p/o-value?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
c       general
      if(linex(ix:jx).eq.'model') model=nint(val)
      if(linex(ix:jx).eq.'iquasiel') iquasiel=nint(val)
      if(linex(ix:jx).eq.'iversn')iversn=nint(val)
      if(linex(ix:jx).eq.'iappl' )iappl=nint(val)
      if(linex(ix:jx).eq.'nevent')nevent=nint(val)
      if(linex(ix:jx).eq.'nfull') nfull=nint(val)
      if(linex(ix:jx).eq.'nfreeze')nfreeze=nint(val)
      if(linex(ix:jx).eq.'ninicon')ninicon=nint(val)
      if(linex(ix:jx).eq.'egymin' )egymin=sngl(val)
      if(linex(ix:jx).eq.'egymax' )egymax=sngl(val)
c       constants
      if(linex(ix:jx).eq.'ainfin')ainfin=sngl(val)
c       printout options
      if(linex(ix:jx).eq.'iprmpt')iprmpt=nint(val)
      if(linex(ix:jx).eq.  'ish' )ish=nint(val)
      if(linex(ix:jx).eq.'ishsub')ishsub=nint(val)
      if(linex(ix:jx).eq.'irandm')irandm=nint(val)
      if(linex(ix:jx).eq.'irewch')irewch=nint(val)
      if(linex(ix:jx).eq.'iecho ')iecho =nint(val)
      if(linex(ix:jx).eq.'modsho')modsho=nint(val)
      if(linex(ix:jx).eq.'idensi')idensi=nint(val)
      if(linex(ix:jx).eq.'ishevt')ishevt=nint(val)
      if(linex(ix:jx).eq.'iwseed')iwseed=nint(val)
c       fragmentation and decay
      if(linex(ix:jx).eq.'pud'   )pud   =sngl(val)
      if(linex(ix:jx).eq.'pude'  )pude  =sngl(val)
      if(linex(ix:jx).eq.'pudk'  )pudk  =sngl(val)
      if(linex(ix:jx).eq.'puds ' )puds   =sngl(val)
      if(linex(ix:jx).eq.'pudr ' )pudr   =sngl(val)
      if(linex(ix:jx).eq.'pudd ' )pudd   =sngl(val)
      if(linex(ix:jx).eq.'strcut' )strcut   =sngl(val)
      if(linex(ix:jx).eq.'diqcut' )diqcut   =sngl(val)
      if(linex(ix:jx).eq.'pdiqua')pdiqua=sngl(val)
      if(linex(ix:jx).eq.'pdiquae')pdiquae=sngl(val)
      if(linex(ix:jx).eq.'pdiquak')pdiquak=sngl(val)
      if(linex(ix:jx).eq.'pdiquar')pdiquar=sngl(val)
      if(linex(ix:jx).eq.'pdiquad')pdiquad=sngl(val)
      if(linex(ix:jx).eq.'ioptf ')ioptf =nint(val)
      if(linex(ix:jx).eq.'delrex')delrex=sngl(val)
      if(linex(ix:jx).eq.'ndecay')ndecay=nint(val)
      if(linex(ix:jx).eq.'maxres')maxres=nint(val)
      if(linex(ix:jx).eq.'pbreak')pbreak=sngl(val)
      if(linex(ix:jx).eq.'pbreake')pbreake=sngl(val)
      if(linex(ix:jx).eq.'pbreakk')pbreakk=sngl(val)
      if(linex(ix:jx).eq.'pbreakr')pbreakr=sngl(val)
      if(linex(ix:jx).eq.'pbreakd')pbreakd=sngl(val)
      if(linex(ix:jx).eq.'ptfra')ptfra=sngl(val)
      if(linex(ix:jx).eq.'ptfrae')ptfrae=sngl(val)
      if(linex(ix:jx).eq.'ptfrak')ptfrak=sngl(val)
      if(linex(ix:jx).eq.'ptfrakd')ptfrakd=sngl(val)
      if(linex(ix:jx).eq.'ptfrair')ptfrair=sngl(val)
      if(linex(ix:jx).eq.'ptfradr')ptfradr=sngl(val)
      if(linex(ix:jx).eq.'ptfrasr')ptfrasr=sngl(val)
      if(linex(ix:jx).eq.'ptfraqq')ptfraqq=sngl(val)
      if(linex(ix:jx).eq.'ptfrab')ptfrab=sngl(val)
      if(linex(ix:jx).eq.'aouni ')aouni=sngl(val)
c       lepton-nucleon and e+e-
      if(linex(ix:jx).eq.'iolept')iolept=nint(val)
      if(linex(ix:jx).eq.'ydmin')ydmin=sngl(val)
      if(linex(ix:jx).eq.'ydmax')ydmax=sngl(val)
      if(linex(ix:jx).eq.'qdmin')qdmin=sngl(val)
      if(linex(ix:jx).eq.'qdmax')qdmax=sngl(val)
      if(linex(ix:jx).eq.'themin')themin=sngl(val)
      if(linex(ix:jx).eq.'themax')themax=sngl(val)
      if(linex(ix:jx).eq.'elomin')elomin=sngl(val)
      if(linex(ix:jx).eq.'engy' )engy=sngl(val)
      if(linex(ix:jx).eq.'elab' )elab=sngl(val)
      if(linex(ix:jx).eq.'ekin' )ekin=sngl(val)
      if(linex(ix:jx).eq.'ecms' )ecms=sngl(val)
      if(linex(ix:jx).eq.'elepti')elepti=sngl(val)
      if(linex(ix:jx).eq.'elepto')elepto=sngl(val)
      if(linex(ix:jx).eq.'angmue')angmue=sngl(val)
      if(linex(ix:jx).eq.'noebin')noebin=nint(val)
      if(linex(ix:jx).eq.'engmin')engmin=sngl(val)
      if(linex(ix:jx).eq.'engmax')engmax=sngl(val)
      if(linex(ix:jx).eq.'iologe')iologe=nint(val)
      if(linex(ix:jx).eq.'iologl')iologl=nint(val)
      if(linex(ix:jx).eq.'itflav')itflav=nint(val)
      if(linex(ix:jx).eq.'idisco')idisco=nint(val)
c       hadron-hadron
      if(linex(ix:jx).eq. 'pnll' )pnll=sngl(val)
      if(linex(ix:jx).eq.'idproj')idprojin=nint(val)
      idproj=idprojin
      if(linex(ix:jx).eq.'idtarg')idtargin=nint(val)
      idtarg=idtargin
      if(idtarg.eq.0)idtarg=1120
      if(linex(ix:jx).eq.'ptq   ')ptq   =sngl(val)
      if(linex(ix:jx).eq.'rstrau(1)')rstrau(1)=sngl(val)
      if(linex(ix:jx).eq.'rstrad(1)')rstrad(1)=sngl(val)
      if(linex(ix:jx).eq.'rstras(1)')rstras(1)=sngl(val)
      if(linex(ix:jx).eq.'rstrau(2)')rstrau(2)=sngl(val)
      if(linex(ix:jx).eq.'rstrad(2)')rstrad(2)=sngl(val)
      if(linex(ix:jx).eq.'rstras(2)')rstras(2)=sngl(val)
      if(linex(ix:jx).eq.'rstrau(3)')rstrau(3)=sngl(val)
      if(linex(ix:jx).eq.'rstrad(3)')rstrad(3)=sngl(val)
      if(linex(ix:jx).eq.'rstras(3)')rstras(3)=sngl(val)
      if(linex(ix:jx).eq.'rstrau(4)')rstrau(4)=sngl(val)
      if(linex(ix:jx).eq.'rstrad(4)')rstrad(4)=sngl(val)
      if(linex(ix:jx).eq.'rstras(4)')rstras(4)=sngl(val)
      if(linex(ix:jx).eq.'rstrasi')rstrasi=sngl(val)
      if(linex(ix:jx).eq.'wgtval')wgtval=sngl(val)
      if(linex(ix:jx).eq.'wgtsea')wgtsea=sngl(val)
      if(linex(ix:jx).eq.'wgtdiq')wgtdiq=sngl(val)
      if(linex(ix:jx).eq.'wgtqqq(1)')wgtqqq(1)=sngl(val)
      if(linex(ix:jx).eq.'wgtqqq(2)')wgtqqq(2)=sngl(val)
      if(linex(ix:jx).eq.'wgtqqq(3)')wgtqqq(3)=sngl(val)
      if(linex(ix:jx).eq.'wgtqqq(4)')wgtqqq(4)=sngl(val)
      if(linex(ix:jx).eq.'wproj ')wproj =sngl(val)
      if(linex(ix:jx).eq.'wtarg ')wtarg =sngl(val)
      if(linex(ix:jx).eq.'wexcit')wexcit=sngl(val)
c      if(linex(ix:jx).eq.'cutmsq')cutmsq=sngl(val)
      if(linex(ix:jx).eq.'cutmss')cutmss=sngl(val)
      if(linex(ix:jx).eq.'exmass')exmass=sngl(val)
      if(linex(ix:jx).eq.'iregge')iregge=nint(val)
      if(linex(ix:jx).eq.'isopom')isopom=nint(val)
      if(linex(ix:jx).eq.'ishpom')ishpom=nint(val)
      if(linex(ix:jx).eq.'iscreen')iscreen=nint(val)
      if(linex(ix:jx).eq.'isplit')isplit=nint(val)
      if(linex(ix:jx).eq.'irzptn')irzptn=nint(val)
      if(linex(ix:jx).eq.'irmdrop')irmdrop=nint(val)
      if(linex(ix:jx).eq.'nprmax')nprmax=nint(val)
      if(linex(ix:jx).eq.'iemspl')iemspl=nint(val)
      if(linex(ix:jx).eq.'intpol')intpol=nint(val)
      if(linex(ix:jx).eq.'isigma')isigma=nint(val)
      if(linex(ix:jx).eq.'iomega')iomega=nint(val)
      if(linex(ix:jx).eq.'isetcs')isetcs=nint(val)
      if(linex(ix:jx).eq.'iemsb' )iemsb= nint(val)
      if(linex(ix:jx).eq.'iemspm')iemspm=nint(val)
      if(linex(ix:jx).eq.'iemspx')iemspx=nint(val)
      if(linex(ix:jx).eq.'iemsse')iemsse=nint(val)
      if(linex(ix:jx).eq.'iemsdr')iemsdr=nint(val)
      if(linex(ix:jx).eq.'iemsrx')iemsrx=nint(val)
      if(linex(ix:jx).eq.'iemsi2')iemsi2=nint(val)
      if(linex(ix:jx).eq.'iemsi1')iemsi1=nint(val)
      if(linex(ix:jx).eq.'ioems' )ioems= nint(val)
      if(linex(ix:jx).eq.'ispacetime')ispacetime= nint(val)
c       unified parameters
      if(linex(ix:jx).eq.'iclpro1')iclpro1=nint(val)
      if(linex(ix:jx).eq.'iclpro2')iclpro2=nint(val)
      if(linex(ix:jx).eq.'icltar1')icltar1=nint(val)
      if(linex(ix:jx).eq.'icltar2')icltar2=nint(val)
      if(linex(ix:jx).eq.'iclegy1')iclegy1=nint(val)
      if(linex(ix:jx).eq.'iclegy2')iclegy2=nint(val)
      if(linex(ix:jx).eq.'egylow')egylow=sngl(val)
      if(linex(ix:jx).eq.'alppom')alppom=sngl(val)
      if(linex(ix:jx).eq.'gamhad(1)')stop'gamhad(1) not allowed'
      if(linex(ix:jx).eq.'gamhad(2)')gamhad(2)=sngl(val)
      if(linex(ix:jx).eq.'gamhad(3)')stop'gamhad(3) not allowed'
      if(linex(ix:jx).eq.'gamhads(1)')gamhadsi(1)=sngl(val)
      if(linex(ix:jx).eq.'gamhads(2)')gamhadsi(2)=sngl(val)
      if(linex(ix:jx).eq.'gamhads(3)')gamhadsi(3)=sngl(val)
      if(linex(ix:jx).eq.'gamhads(4)')gamhadsi(4)=sngl(val)
      if(linex(ix:jx).eq.'gamtil')gamtil=sngl(val)
      if(linex(ix:jx).eq.'slopom')slopom=sngl(val)
      if(linex(ix:jx).eq.'slopoms')slopoms=sngl(val)
      if(linex(ix:jx).eq.'r2had(1)' )r2had(1)= sngl(val)
      if(linex(ix:jx).eq.'r2had(2)' )r2had(2)= sngl(val)
      if(linex(ix:jx).eq.'r2had(3)' )r2had(3)= sngl(val)
      if(linex(ix:jx).eq.'r2had(4)' )r2had(4)= sngl(val)
      if(linex(ix:jx).eq.'r2hads(1)' )r2hads(1)= sngl(val)
      if(linex(ix:jx).eq.'r2hads(2)' )r2hads(2)= sngl(val)
      if(linex(ix:jx).eq.'r2hads(3)' )r2hads(3)= sngl(val)
      if(linex(ix:jx).eq.'r2hads(4)' )r2hads(4)= sngl(val)
      if(linex(ix:jx).eq.'r3pom'   )r3pom= sngl(val)
      if(linex(ix:jx).eq.'r4pom'   )r4pom= sngl(val)
      if(linex(ix:jx).eq.'egyscr'  )egyscr= sngl(val)
      if(linex(ix:jx).eq.'epscrw'  )epscrw= sngl(val)
      if(linex(ix:jx).eq.'epscrx'  )epscrx= sngl(val)
      if(linex(ix:jx).eq.'epscrs'  )epscrs= sngl(val)
      if(linex(ix:jx).eq.'epscrh'  )epscrh= sngl(val)
      if(linex(ix:jx).eq.'epscrp'  )epscrp= sngl(val)
      if(linex(ix:jx).eq.'epscrv'  )epscrv= sngl(val)
      if(linex(ix:jx).eq.'gfactor' )gfactor= sngl(val)
      if(linex(ix:jx).eq.'gwidth'  )gwidth= sngl(val)
      if(linex(ix:jx).eq.'chad(1)' )chad(1)=  sngl(val)
      if(linex(ix:jx).eq.'chad(2)' )chad(2)=  sngl(val)
      if(linex(ix:jx).eq.'chad(3)' )chad(3)=  sngl(val)
      if(linex(ix:jx).eq.'chad(4)' )chad(4)=  sngl(val)
      if(linex(ix:jx).eq.'wdiff(1)')wdiff(1)=sngl(val)
      if(linex(ix:jx).eq.'wdiff(2)')wdiff(2)=sngl(val)
      if(linex(ix:jx).eq.'wdiff(3)')wdiff(3)=sngl(val)
      if(linex(ix:jx).eq.'wdiff(4)')wdiff(4)=sngl(val)
      if(linex(ix:jx).eq.'facdif')facdif=sngl(val)
      if(linex(ix:jx).eq.'rexndf')   rexndf=sngl(val)
      if(linex(ix:jx).eq.'rexndi(1)')rexndii(1)=sngl(val)
      if(linex(ix:jx).eq.'rexndi(2)')rexndii(2)=sngl(val)
      if(linex(ix:jx).eq.'rexndi(3)')rexndii(3)=sngl(val)
      if(linex(ix:jx).eq.'rexndi(4)')rexndii(4)=sngl(val)
      if(linex(ix:jx).eq.'rexdif(1)')rexdifi(1)=sngl(val)
      if(linex(ix:jx).eq.'rexdif(2)')rexdifi(2)=sngl(val)
      if(linex(ix:jx).eq.'rexdif(3)')rexdifi(3)=sngl(val)
      if(linex(ix:jx).eq.'rexdif(4)')rexdifi(4)=sngl(val)
      if(linex(ix:jx).eq.'rexres')rexres=sngl(val)
      if(linex(ix:jx).eq.'zbarfl')zbarfl=sngl(val)
      if(linex(ix:jx).eq.'alpreg')alpreg=sngl(val)
      if(linex(ix:jx).eq.'sloreg')sloreg=sngl(val)
      if(linex(ix:jx).eq.'gamreg')gamreg=sngl(val)
      if(linex(ix:jx).eq.'r2reg' )r2reg= sngl(val)
c      if(linex(ix:jx).eq.'amhdibar')amhdibar= sngl(val)
      if(linex(ix:jx).eq.'ptdiff')ptdiff=sngl(val)
      if(linex(ix:jx).eq.'alppar')alppar=sngl(val)
      if(linex(ix:jx).eq.'alpsea')alpsea=sngl(val)
      if(linex(ix:jx).eq.'alpval')alpval=sngl(val)
      if(linex(ix:jx).eq.'alpdiq')alpdiq=sngl(val)
      if(linex(ix:jx).eq.'alplea(1)')alplea(1)=sngl(val)
      if(linex(ix:jx).eq.'alplea(2)')alplea(2)=sngl(val)
      if(linex(ix:jx).eq.'alplea(3)')alplea(3)=sngl(val)
      if(linex(ix:jx).eq.'alplea(4)')alplea(4)=sngl(val)
      if(linex(ix:jx).eq.'alpdif')alpdif=sngl(val)
      if(linex(ix:jx).eq.'alpdi')alpdi=sngl(val)
      if(linex(ix:jx).eq.'alpndi')alpndi=sngl(val)
      if(linex(ix:jx).eq.'ammsqq')ammsqq=sngl(val)
      if(linex(ix:jx).eq.'ammsqd')ammsqd=sngl(val)
      if(linex(ix:jx).eq.'ammsdd')ammsdd=sngl(val)
      if(linex(ix:jx).eq.'cumpom')cumpom=sngl(val)
      if(linex(ix:jx).eq.'reminv')reminv=sngl(val)
      if(linex(ix:jx).eq.'ptsend')ptsend=sngl(val)
      if(linex(ix:jx).eq.'ptsendi')ptsendi=sngl(val)
      if(linex(ix:jx).eq.'ptsemx')ptsemx=sngl(val)
      if(linex(ix:jx).eq.'zrminc')zrminc=sngl(val)
      if(linex(ix:jx).eq.'zrmmax')zrmmax=sngl(val)
      if(linex(ix:jx).eq.'edmaxi')edmaxi=sngl(val)
      if(linex(ix:jx).eq.'epmaxi')epmaxi=sngl(val)
      if(linex(ix:jx).eq.'zopinc')zopinc=sngl(val)
      if(linex(ix:jx).eq.'zopmax')zopmax=sngl(val)
      if(linex(ix:jx).eq.'zipinc')zipinc=sngl(val)
      if(linex(ix:jx).eq.'zipmax')zipmax=sngl(val)
      if(linex(ix:jx).eq.'zidinc')zidinc=sngl(val)
      if(linex(ix:jx).eq.'zidmax')zidmax=sngl(val)
      if(linex(ix:jx).eq.'zodinc')zodinc=sngl(val)
      if(linex(ix:jx).eq.'zodmax')zodmax=sngl(val)
      if(linex(ix:jx).eq.'zosinc')zosinc=sngl(val)
      if(linex(ix:jx).eq.'zosmax')zosmax=sngl(val)
      if(linex(ix:jx).eq.'zoeinc')zoeinc=sngl(val)
      if(linex(ix:jx).eq.'zoemax')zoemax=sngl(val)
      if(linex(ix:jx).eq.'zbcut') zbcut=sngl(val)
      if(linex(ix:jx).eq.'xmaxremn')xmaxremn=sngl(val)
      if(linex(ix:jx).eq.'xmaxdiff')xmaxdiff=sngl(val)
      if(linex(ix:jx).eq.'alpdro(1)')alpdro(1)=sngl(val)
      if(linex(ix:jx).eq.'alpdro(2)')alpdro(2)=sngl(val)
      if(linex(ix:jx).eq.'alpdro(3)')alpdro(3)=sngl(val)
      if(linex(ix:jx).eq.'amdrmax')amdrmax=sngl(val)
      if(linex(ix:jx).eq.'iodiba')iodiba=nint(val)
      if(linex(ix:jx).eq.'bidiba')bidiba=sngl(val)

c       hard pomeron parameters
      if(linex(ix:jx).eq.'q2min' )q2min=sngl(val)
      if(linex(ix:jx).eq.'q2ini' )q2ini=sngl(val)
      if(linex(ix:jx).eq.'q2fin' )q2fin=sngl(val)
      if(linex(ix:jx).eq.'betpom')betpom=sngl(val)
      if(linex(ix:jx).eq.'alpfom')alpfomi=sngl(val)
      if(linex(ix:jx).eq.'betfom')betfom=dble(val)
      if(linex(ix:jx).eq.'gamfom')gamfom=dble(val)
      if(linex(ix:jx).eq.'glusea')glusea=sngl(val)
      if(linex(ix:jx).eq.'factk' )factk=sngl(val)
      if(linex(ix:jx).eq.'naflav')naflav=nint(val)
      if(linex(ix:jx).eq.'pt2cut')pt2cut=sngl(val)
      if(linex(ix:jx).eq.'factgam')factgam=sngl(val)
      if(linex(ix:jx).eq.'delh')delh=sngl(val)
c       nucleus-nucleus
      if(linex(ix:jx).eq.'iokoll')iokoll=nint(val)
      if(linex(ix:jx).eq.'laproj')laproj=nint(val)
      if(linex(ix:jx).eq.'maproj')maproj=nint(val)
      if(linex(ix:jx).eq.'latarg')latarg=nint(val)
      if(linex(ix:jx).eq.'matarg')matarg=nint(val)
      if(linex(ix:jx).eq.'core'  )core  =sngl(val)
c      if(linex(ix:jx).eq.'ncolmx')ncolmx=nint(val)
      if(linex(ix:jx).eq.'fctrmx')fctrmx=sngl(val)
      if(linex(ix:jx).eq.'bmaxim')bmaxim=sngl(val)
      if(linex(ix:jx).eq.'bminim')bminim=sngl(val)
      if(linex(ix:jx).eq.'phimax')phimax=sngl(val)
      if(linex(ix:jx).eq.'phimin')phimin=sngl(val)
c       rescattering parameters
      if(linex(ix:jx).eq.'iorsce')iorsce=nint(val)
      if(linex(ix:jx).eq.'iorsdf')iorsdf=nint(val)
      if(linex(ix:jx).eq.'iorshh')iorshh=nint(val)
      if(linex(ix:jx).eq.'iocluin')iocluin=nint(val)
      if(linex(ix:jx).eq.'ioquen')ioquen=nint(val)
      if(linex(ix:jx).eq.'hacore')hacore=sngl(val)
      if(linex(ix:jx).eq.'amimfs')amimfs=sngl(val)
      if(linex(ix:jx).eq.'amimel')amimel=sngl(val)
      if(linex(ix:jx).eq.'cepara')cepara=sngl(val)
      if(linex(ix:jx).eq.'dscale')dscale=sngl(val)
      if(linex(ix:jx).eq.'iceopt')iceopt=nint(val)
      if(linex(ix:jx).eq.'delamf')delamf=sngl(val)
      if(linex(ix:jx).eq.'deuamf')deuamf=sngl(val)
      if(linex(ix:jx).eq.'taurea')taurea=sngl(val)
      if(linex(ix:jx).eq.'nsegsu')nsegsu=nint(val)
      if(linex(ix:jx).eq.'nsegce')nsegce=nint(val)
      if(linex(ix:jx).eq.'ptclu') ptclu=sngl(val)
      if(linex(ix:jx).eq.'epscri(1)')epscri(1)=sngl(val)
      if(linex(ix:jx).eq.'epscri(3)')epscri(3)=sngl(val)
      if(linex(ix:jx).eq.'amsiac')amsiac=sngl(val)
      if(linex(ix:jx).eq.'amprif')amprif=sngl(val)
      if(linex(ix:jx).eq.'delvol')delvol=sngl(val)
      if(linex(ix:jx).eq.'deleps')deleps=sngl(val)
      if(linex(ix:jx).eq.'taumin')taumin=sngl(val)
      if(linex(ix:jx).eq.'deltau')deltau=sngl(val)
      if(linex(ix:jx).eq.'factau')factau=sngl(val)
      if(linex(ix:jx).eq.'numtau')numtau=nint(val)
      if(linex(ix:jx).eq.'dlzeta')dlzeta=sngl(val)
      if(linex(ix:jx).eq.'etafac')etafac=sngl(val)
      if(linex(ix:jx).eq.'facnuc')facnuc=sngl(val)
c       urqmd
      if(linex(ix:jx).eq.'iurqmd')iurqmd=int(val)
c       spherio
      if(linex(ix:jx).eq.'ispherio')ispherio=int(val)
c       ico
      if(linex(ix:jx).eq.'cutico') cutico=sngl(val)
      if(linex(ix:jx).eq.'dssico') dssico=sngl(val)
      if(linex(ix:jx).eq.'icocore')icocore=int(val)
      if(linex(ix:jx).eq.'icotabm')icotabm=int(val)
      if(linex(ix:jx).eq.'icotabr')icotabr=int(val)
c       droplet decay
      if(linex(ix:jx).eq.'dezzer')dezzer=sngl(val)
      if(linex(ix:jx).eq.'ioclude')ioclude=nint(val)
      if(linex(ix:jx).eq.'amuseg')amuseg=sngl(val)
      if(linex(ix:jx).eq.'yradmx')yradmx=sngl(val)
      if(linex(ix:jx).eq.'yradmi')yradmi=sngl(val)
      if(linex(ix:jx).eq.'yradpp')yradpp=sngl(val)
      if(linex(ix:jx).eq.'yradpi')yradpi=sngl(val)
      if(linex(ix:jx).eq.'yradpx')yradpx=sngl(val)
      if(linex(ix:jx).eq.'facecc')facecc=sngl(val)
      if(linex(ix:jx).eq.'rcoll' )rcoll= sngl(val)
      if(linex(ix:jx).eq.'ylongmx' )ylongmx= sngl(val)
      if(linex(ix:jx).eq.'bag4rt')bag4rt=sngl(val)
      if(linex(ix:jx).eq.'taunll')taunll=sngl(val)
c       droplet specification
      if(linex(ix:jx).eq. 'keu'  )keu=nint(val)
      if(linex(ix:jx).eq. 'ked'  )ked=nint(val)
      if(linex(ix:jx).eq. 'kes'  )kes=nint(val)
      if(linex(ix:jx).eq. 'kec'  )kec=nint(val)
      if(linex(ix:jx).eq. 'keb'  )keb=nint(val)
      if(linex(ix:jx).eq. 'ket'  )ket=nint(val)
      if(linex(ix:jx).eq. 'tecm' )tecm=sngl(val)
      if(linex(ix:jx).eq. 'volu' )volu=sngl(val)
      if(linex(ix:jx).eq. 'vrad' )vrad=sngl(val)
      if(linex(ix:jx).eq. 'facts')facts=sngl(val)
      if(linex(ix:jx).eq. 'factb')factb=sngl(val)
      if(linex(ix:jx).eq. 'factq')factq=sngl(val)
      if(linex(ix:jx).eq.'inbxxx')inbxxx=sngl(val)
c       metropolis
      if(linex(ix:jx).eq.'iospec')iospec=nint(val)
      if(linex(ix:jx).eq.'iocova')iocova=nint(val)
      if(linex(ix:jx).eq.'iopair')iopair=nint(val)
      if(linex(ix:jx).eq.'iozero')iozero=nint(val)
      if(linex(ix:jx).eq.'ioflac')ioflac=nint(val)
      if(linex(ix:jx).eq.'iostat')iostat=nint(val)
      if(linex(ix:jx).eq.'ioinco')ioinco=nint(val)
      if(linex(ix:jx).eq.'iograc')iograc=nint(val)
      if(linex(ix:jx).eq.'epsgc' )epsgc=sngl(val)
      if(linex(ix:jx).eq.'iocite')iocite=nint(val)
      if(linex(ix:jx).eq.'ioceau')ioceau=nint(val)
      if(linex(ix:jx).eq.'iociau')iociau=nint(val)
      if(linex(ix:jx).eq.'ioinct')ioinct=nint(val)
      if(linex(ix:jx).eq.'ioinfl')ioinfl=nint(val)
      if(linex(ix:jx).eq.'iowidn')iowidn=nint(val)
      if(linex(ix:jx).eq.'ionlat')ionlat=nint(val)
      if(linex(ix:jx).eq.'iomom')iomom=nint(val)
      if(linex(ix:jx).eq.'ioobsv')ioobsv=nint(val)
      if(linex(ix:jx).eq.'iosngl')iosngl=nint(val)
      if(linex(ix:jx).eq.'iorejz')iorejz=nint(val)
      if(linex(ix:jx).eq.'iompar')iompar=nint(val)
      if(linex(ix:jx).eq.'iozinc')iozinc=nint(val)
      if(linex(ix:jx).eq.'iozevt')iozevt=nint(val)
      if(linex(ix:jx).eq. 'nadd' )nadd=nint(val)
      if(linex(ix:jx).eq.'iterma')iterma=nint(val)
      if(linex(ix:jx).eq.'itermx')stop'STOP: set iterma, not itermx'
      if(linex(ix:jx).eq.'iterpr')iterpr=nint(val)
      if(linex(ix:jx).eq.'iterpl')iterpl=nint(val)
      if(linex(ix:jx).eq.'iternc')iternc=nint(val)
      if(linex(ix:jx).eq.'epsr'  )epsr=sngl(val)
      if(linex(ix:jx).eq.'keepr' )keepr=nint(val)
c       strangelets
      if(linex(ix:jx).eq.'iopenu')iopenu=nint(val)
      if(linex(ix:jx).eq.'themas')themas=sngl(val)
c       tests
      if(linex(ix:jx).eq.'iotst1')iotst1=nint(val)
      if(linex(ix:jx).eq.'iotst2')iotst2=nint(val)
      if(linex(ix:jx).eq.'iotst3')iotst3=nint(val)
      if(linex(ix:jx).eq.'iotst4')iotst4=nint(val)
c       jpsi
      if(linex(ix:jx).eq.'jpsi  ')jpsi  =nint(val)
      if(linex(ix:jx).eq.'jpsifi')jpsifi=nint(val)
      if(linex(ix:jx).eq.'sigj  ')sigj  =sngl(val)
      if(linex(ix:jx).eq.'taumx ')taumx =sngl(val)
      if(linex(ix:jx).eq.'nsttau')nsttau=nint(val)
      if(linex(ix:jx).eq.'ijphis')ijphis=nint(val)
c       analysis: intermittency, space-time, droplets, formation time
      if(linex(ix:jx).eq.'ymximi')ymximi=sngl(val)
      if(linex(ix:jx).eq.'imihis')imihis=nint(val)
      if(linex(ix:jx).eq.'isphis')isphis=nint(val)
      if(linex(ix:jx).eq.'iologb')iologb=nint(val)
      if(linex(ix:jx).eq.'ispall')ispall=nint(val)
      if(linex(ix:jx).eq.'wtmini')wtmini=sngl(val)
      if(linex(ix:jx).eq.'wtstep')wtstep=sngl(val)
      if(linex(ix:jx).eq.'iwcent')iwcent=nint(val)
      if(linex(ix:jx).eq.'iclhis')iclhis=nint(val)
      if(linex(ix:jx).eq.'iwtime')iwtime=nint(val)
      if(linex(ix:jx).eq.'wtimet')wtimet=sngl(val)
      if(linex(ix:jx).eq.'wtimei')wtimei=sngl(val)
      if(linex(ix:jx).eq.'wtimea')wtimea=sngl(val)
c       other
      if(linex(ix:jx).eq.'gaumx ')gaumx =sngl(val)
      if(linex(ix:jx).eq.'nclean')nclean=nint(val)
      if(linex(ix:jx).eq.'istore')istore=nint(val)
      if(linex(ix:jx).eq.'ioidch')ioidch=nint(val)
      if(linex(ix:jx).eq.'iframe')iframe=nint(val)
      if(linex(ix:jx).eq.'jframe')jframe=nint(val)
      if(linex(ix:jx).eq.'labsys')stop'labsys no longer supported'
      if(linex(ix:jx).eq.'irescl')irescl=nint(val)
      if(linex(ix:jx).eq.'ifrade')ifrade=nint(val)
      if(linex(ix:jx).eq.'idecay')idecay=nint(val)
      if(linex(ix:jx).eq.'jdecay')jdecay=nint(val)
      if(linex(ix:jx).eq.'ntrymx')ntrymx=nint(val)
      if(linex(ix:jx).eq.'istmax')istmax=nint(val)
      if(linex(ix:jx).eq.'ionudi')ionudi=nint(val)
      if(linex(ix:jx).eq.'seedi')seedi =val
      if(linex(ix:jx).eq.'seedj')seedj =val
      if(linex(ix:jx).eq.'xvaria')xvaria='      '
      if(linex(ix:jx).eq.'xvaria')xvaria(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'yvaria')yvaria='      '
      if(linex(ix:jx).eq.'yvaria')yvaria(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'normal')normal=nint(val)
      if(linex(ix:jx).eq.'xminim')xminim=sngl(val)
      if(linex(ix:jx).eq.'xmaxim')xmaxim=sngl(val)
      if(linex(ix:jx).eq.'nrbins')nrbins=nint(val)
      if(linex(ix:jx).eq.'hisfac')hisfac=sngl(val)
      if(linex(ix:jx).eq.'xshift')xshift=sngl(val)
      if(linex(ix:jx).eq.'etacut')etacut=sngl(val)
      if(linex(ix:jx).eq.'xpar1' )xpar1=sngl(val)
      if(linex(ix:jx).eq.'xpar2' )xpar2=sngl(val)
      if(linex(ix:jx).eq.'xpar3' )xpar3=sngl(val)
      if(linex(ix:jx).eq.'xpar4' )xpar4=sngl(val)
      if(linex(ix:jx).eq.'xpar5' )xpar5=sngl(val)
      if(linex(ix:jx).eq.'xpar6' )xpar6=sngl(val)
      if(linex(ix:jx).eq.'xpar7' )xpar7=sngl(val)
      if(linex(ix:jx).eq.'xpar8' )xpar8=sngl(val)
      if(linex(ix:jx).eq.'irdmpr')irdmpr=nint(val)
c       frame definitions
      if(linex(ix:jx).eq.'engy'.and.iframe.eq.0)iframe=11
      if(linex(ix:jx).eq.'ecms'.and.iframe.eq.0)iframe=11
      if(linex(ix:jx).eq.'elab'.and.iframe.eq.0)iframe=12
      if(linex(ix:jx).eq.'ekin'.and.iframe.eq.0)iframe=12
      if(linex(ix:jx).eq.'pnll'.and.iframe.eq.0)iframe=12
      if(linex(ix:jx).eq.'noebin'.and.dabs(val-1.d0).gt.1.d-10)iframe=11

           elseif(line(i:j).eq.'set')then

      call utword(line,i,j,0)
      call utword(line,i,j,0)
      ifset=1

           elseif(line(i:j).eq.'stop')then  !same as return

      if(nopen.ne.-1)then
      close(20+nopen)
      nopen=nopen-1
      if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
      endif

           elseif(line(i:j).eq.'stopprogram')then

      close(unit=ifcx)
      close(unit=ifhi)
      close(unit=ifdt)
      stop

           elseif(line(i:j).eq.'EndEposInput')then

      return

           elseif(line(i:j).eq.'string')then

      nstmax=nstmax+1
      ns=nstmax
      icinpu=0
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: prob icbac1 icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      prob(ns)=sngl(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icbac1 icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icbac(ns,1)=nint(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icbac(ns,2)=nint(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icfor(ns,1)=nint(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icfor(ns,2)=nint(val)

           elseif(line(i:j).eq.'kinks')then

      nptl=0
ctp290806      nrow=0
      nel=0
 10   continue
      call utword(line,i,j,0)
      if(line(i:j).eq.'endkinks')goto 12
      nel=nel+1
ctp290806      nrow=1+(nel-1)/4
      nc=mod(nel-1,4)+1
      read(line(i:j),*)a
      if(nc.eq.1)nptl=nptl+1
      if(nc.eq.1)idptl(nptl)=nint(a)
      if(nc.eq.2)pptl(1,nptl)=a
      if(nc.eq.3)pptl(2,nptl)=a
      if(nc.eq.4)then
        pptl(3,nptl)=a
        istptl(nptl)=20
        pptl(4,nptl)=sqrt(pptl(3,nptl)**2+pptl(2,nptl)**2
     $       +pptl(1,nptl)**2)
      endif
      goto10
 12   continue

           elseif(line(i:j).eq.'record')then

      call utworn(line,j,ne)
c      if(ne.eq.0.and.iprmpt.gt.0)
c     *     write(ifmt,'(a)')'kinks: icbac1 icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'event')ir=1
      if(line(i:j).eq.'particle')ir=2
      maxrec(ir)=0
 20   call utworn(line,j,ne)
c      if(ne.eq.0.and.iprmpt.gt.0)
c     *     write(6,'(a)')'<kinks-data (px-py-pz)>? (End=endkinks)'
      call utword(line,i,j,0)
      if(line(i:j).eq.'endrecord')then
         goto 22
      endif
      maxrec(ir)=maxrec(ir)+1
      irecty(maxrec(ir),ir)=0
      if(ir.eq.1)then
        if(line(i:j).eq.'0') irecty(maxrec(ir),ir)=0
        if(line(i:j).eq.'nevt') irecty(maxrec(ir),ir)=1
        if(line(i:j).eq.'nptl') irecty(maxrec(ir),ir)=2
        if(line(i:j).eq.'b') irecty(maxrec(ir),ir)=3
        if(line(i:j).eq.'phi') irecty(maxrec(ir),ir)=4
        if(line(i:j).eq.'ncol') irecty(maxrec(ir),ir)=5
        if(line(i:j).eq.'pmx') irecty(maxrec(ir),ir)=6
        if(line(i:j).eq.'egy') irecty(maxrec(ir),ir)=7
        if(line(i:j).eq.'npj') irecty(maxrec(ir),ir)=8
        if(line(i:j).eq.'ntg') irecty(maxrec(ir),ir)=9
        if(line(i:j).eq.'npn') irecty(maxrec(ir),ir)=10
        if(line(i:j).eq.'npp') irecty(maxrec(ir),ir)=11
        if(line(i:j).eq.'ntn') irecty(maxrec(ir),ir)=12
        if(line(i:j).eq.'ntp') irecty(maxrec(ir),ir)=13
        if(line(i:j).eq.'jpn') irecty(maxrec(ir),ir)=14
        if(line(i:j).eq.'jpp') irecty(maxrec(ir),ir)=15
        if(line(i:j).eq.'jtn') irecty(maxrec(ir),ir)=16
        if(line(i:j).eq.'jtp') irecty(maxrec(ir),ir)=17
        if(line(i:j).eq.'amproj') irecty(maxrec(ir),ir)=20
        if(line(i:j).eq.'amtarg') irecty(maxrec(ir),ir)=21
      else
        if(line(i:j).eq.'0') irecty(maxrec(ir),ir)=0
        if(line(i:j).eq.'i') irecty(maxrec(ir),ir)=1
        if(line(i:j).eq.'id') irecty(maxrec(ir),ir)=2
        if(line(i:j).eq.'p1') irecty(maxrec(ir),ir)=3
        if(line(i:j).eq.'p2') irecty(maxrec(ir),ir)=4
        if(line(i:j).eq.'p3') irecty(maxrec(ir),ir)=5
        if(line(i:j).eq.'p4') irecty(maxrec(ir),ir)=6
        if(line(i:j).eq.'p5') irecty(maxrec(ir),ir)=7
        if(line(i:j).eq.'fa') irecty(maxrec(ir),ir)=8
        if(line(i:j).eq.'mo') irecty(maxrec(ir),ir)=9
        if(line(i:j).eq.'st') irecty(maxrec(ir),ir)=10
        if(line(i:j).eq.'x1') irecty(maxrec(ir),ir)=11
        if(line(i:j).eq.'x2') irecty(maxrec(ir),ir)=12
        if(line(i:j).eq.'x3') irecty(maxrec(ir),ir)=13
        if(line(i:j).eq.'x4') irecty(maxrec(ir),ir)=14
        if(line(i:j).eq.'idfa') irecty(maxrec(ir),ir)=15
        if(line(i:j).eq.'idmo') irecty(maxrec(ir),ir)=16
        if(line(i:j).eq.'p') irecty(maxrec(ir),ir)=17
        if(line(i:j).eq.'x') irecty(maxrec(ir),ir)=18
        if(line(i:j).eq.'dez') irecty(maxrec(ir),ir)=19
        if(line(i:j).eq.'c1') irecty(maxrec(ir),ir)=21
        if(line(i:j).eq.'c2') irecty(maxrec(ir),ir)=22
        if(line(i:j).eq.'ty') irecty(maxrec(ir),ir)=23
      endif
      if(irecty(maxrec(ir),ir).eq.0)then
        write(*,*) 'unknown variable ',line(i:j)
        stop
      endif
      goto 20
 22   continue

           elseif(line(i:j).eq.'switch')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'option on/off?'
      call utword(line,i,j,0)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'on/off?'
        if(line(i:j).eq.'droplet')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' )iorsdf=1
      if(line(i:j).eq.'off')iorsdf=0
        elseif(line(i:j).eq.'cascade')then
      call utword(line,i,j,0)
c      if(line(i:j).eq.'on' )iorsce=1
      if(line(i:j).eq.'on' )iorshh=1
c      if(line(i:j).eq.'off')iorsce=0
      if(line(i:j).eq.'off')iorshh=0
        elseif(line(i:j).eq.'soft')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' )isopom=1
      if(line(i:j).eq.'off')isopom=0
        elseif(line(i:j).eq.'hard')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' )ishpom=1
      if(line(i:j).eq.'off')ishpom=0
        elseif(line(i:j).eq.'splitting')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' )isplit=1
      if(line(i:j).eq.'on' )iscreen=1
      if(line(i:j).eq.'off')isplit=0
      if(line(i:j).eq.'off')iscreen=0
        elseif(line(i:j).eq.'fusion')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' )iorsce=0
      if(line(i:j).eq.'on' )iorsdf=3
      if(line(i:j).eq.'on' )iorshh=0
      if(line(i:j).eq.'off')iorsce=0
      if(line(i:j).eq.'off')iorsdf=0
      if(line(i:j).eq.'off')iorshh=0
        elseif(line(i:j).eq.'urqmd')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' ) iurqmd=1
      if(line(i:j).eq.'off') iurqmd=0
        elseif(line(i:j).eq.'spherio')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' ) ispherio=1
      if(line(i:j).eq.'off') ispherio=0
        elseif(line(i:j).eq.'decay')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' ) ndecay=0
      if(line(i:j).eq.'off') ndecay=1
      if(line(i:j).eq.'on' ) idecay=1
      if(line(i:j).eq.'off') idecay=0
        elseif(line(i:j).eq.'clusterdecay')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' ) jdecay=1
      if(line(i:j).eq.'off') jdecay=0
        elseif(line(i:j).eq.'fragdecay')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'on' ) ifrade=1
      if(line(i:j).eq.'off') ifrade=0
        elseif(line(i:j).eq.'icocore')then
          stop'switch icocore not supported any more.    '
        endif

           elseif(line(i:j).eq.'idchoice')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'nxs')then
        ioidch=1
      elseif(line(i:j).eq.'pdf')then
        ioidch=2
      else
        stop'invalid idchoice.     '
      endif

           elseif(line(i:j).eq.'make')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'icotable')icotabm=1

           elseif(line(i:j).eq.'read')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'icotable')icotabr=1

           elseif(line(i:j).eq.'output')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'epos' )     istore=1
      if(line(i:j).eq.'osc1997a' )  istore=2
      if(line(i:j).eq.'osc1999a' )  istore=3
      if(line(i:j).eq.'ustore' )    istore=5

           elseif(line(i:j).eq.'model')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'epos')then
        model=1
      else
        nij=j-i+1
        if(nij.gt.20)stop'cmodel too small'
        cmodel(1:nij)=line(i:j)
        cmodel(nij+1:nij+1)=' '
        call NumberModel(cmodel,model)
      endif
      if(iappl.ne.1.and.iappl.ne.3.and.model.ne.1)
     &call utstop('Application not possible with this model&')

           elseif(line(i:j).eq.'trigger')then

      call utword(line,i,j,0)
      ntc=1
      if(line(i:j).eq.'or'.or.line(i:j).eq.'contr')then
        call utword(line,i,j,0)
        read(line(i:j),*)ztc
        ntc=nint(ztc)
        call utword(line,i,j,1)
      endif
      do n=1,ntc
        if(n.ne.1)call utword(line,i,j,0)
        call utword(line,i,j,0)
        call utword(line,i,j,0)
      enddo

           elseif(line(i:j).eq.'noerrorbut')then

      call utword(line,i,j,0)

           elseif(line(i:j).eq.'b')then

      nbarray=nbarray+1
      call utword(line,i,j,0)
      read(line(i:j),*)val
      barray(nbarray)=val

           elseif(line(i:j).eq.'message')then

      call utword(line,i,j,0)
      if(nopen.eq.-1)then      !only write in second read
      write(ifmt,'(a,$)')line(i:j)
      endif

           elseif(line(i:j).eq.'endmessage')then

      if(nopen.eq.-1)then      !only write in second read
      write(ifmt,'(a)')' '
      endif

           elseif(line(i:j).eq.'write'.or.line(i:j).eq.'writex')then

      if(line(i:j).eq.'write')ii=1
      if(line(i:j).eq.'writex')ii=2
      call utword(line,i,j,0)
      idol=0
      if(line(i:j).eq.'$')then
       idol=1
       call utword(line,i,j,0)
      endif
      !write: only write in second read; writex: only write in first read
      if(ii.eq.1.and.nopen.eq.-1.or.ii.eq.2.and.nopen.ne.-1)then
       if(j-6.ge.i)then
        do l=i,j-6
         if(line(l:l+6).eq.'$iversn')then
          if(mod(iversn,100).le.9)then
           write(line(l:l+6),'(i1,a,i1,3x)')
     *     int(iversn/100),'.0',mod(iversn,100)
          else
           write(line(l:l+6),'(i1,a,i2,3x)')
     *     int(iversn/100),'.',mod(iversn,100)
          endif
         endif
        enddo
       endif
       if(j-8.ge.i)then
        do l=i,j-8
         if(line(l:l+8).eq.'$xxxyield')then
          write(line(l:l+8),'(f8.6,1x)')yield
         endif
         if(line(l:l+7).eq.'$xxyield')then
          write(line(l:l+8),'(f7.5,1x)')yield
         endif
         if(line(l:l+6).eq.'$xyield')then
          write(line(l:l+8),'(f6.4,1x)')yield
         endif
        enddo
       endif
       if(j-5.ge.i)then
        do l=i,j-5
         if(line(l:l+5).eq.'$yield')then
          if(yield.lt.1.0)then
           write(line(l:l+5),'(f5.3,1x)')yield
          elseif(yield.lt.100.)then
           write(line(l:l+5),'(f5.2,1x)')yield
          elseif(yield.lt.1000.)then
           write(line(l:l+5),'(f5.1,1x)')yield
          elseif(yield.lt.10000.)then
           write(line(l:l+5),'(f6.1)')yield
          else
           write(line(l:l+5),'(i6)')nint(yield)
          endif
         endif
        enddo
        do l=i,j-5
         if(line(l:l+5).eq.'$averg')then
          if(averg.lt.1.0)then
           write(line(l:l+5),'(f5.3,1x)')averg
          elseif(averg.lt.100.)then
           write(line(l:l+5),'(f5.2,1x)')averg
          elseif(averg.lt.1000.)then
           write(line(l:l+5),'(f5.1,1x)')averg
          elseif(averg.lt.10000.)then
           write(line(l:l+5),'(f6.1)')averg
          else
           write(line(l:l+5),'(i6)')nint(averg)
          endif
         endif
        enddo
        do l=i,j-5
         if(line(l:l+5).eq.'$sigma')then
          if(sigma.lt.1.0)then
           write(line(l:l+5),'(f5.3,1x)')sigma
          elseif(sigma.lt.100.)then
           write(line(l:l+5),'(f5.2,1x)')sigma
          elseif(sigma.lt.1000.)then
           write(line(l:l+5),'(f5.1,1x)')sigma
          elseif(sigma.lt.10000.)then
           write(line(l:l+5),'(f6.1)')sigma
          else
           write(line(l:l+5),'(i6)')nint(sigma)
          endif
         endif
        enddo
       endif
       if(idol.eq.0)then
        write(ifhi,'(a)')line(i:j)
       else
        write(ifhi,'(a,a,$)')line(i:j),' '
       endif
      endif

           elseif(line(i:j).eq.'nozero')then

      nozero=1

           elseif(line(i:j).eq.'ibmin')then

      call utword(line,i,j,0)
      read(line(i:j),*)val
      ibmin=nint(val)

           elseif(line(i:j).eq.'ibmax')then

      call utword(line,i,j,0)
      read(line(i:j),*)val
      ibmax=nint(val)

            elseif(line(i:j).eq.'writearray'
     $     .or.line(i:j).eq.'writehisto')then

      if(nopen.eq.-1)then !second run
       ih=0
       if(line(i:j).eq.'writearray') ih=1
       call utword(line,i,j,0)
       if(line(i:j).eq.'s')then
        call utword(line,i,j,0)
        linex=line
        ix=i
        jx=j
        call utword(line,i,j,0)
        if(linex(ix:jx).eq.'inicon')stop'error 060307'
       else
        ioint=0
        iocontr=0
        if(line(i:j).eq.'int')then
         ioint=1
         call utword(line,i,j,0)
        endif
        if(line(i:j).eq.'contr')then
         iocontr=1
         call utword(line,i,j,0)
         read(line(i:j),*)val
         ncontr=nint(val)
         call utword(line,i,j,0)
        endif
        read(line(i:j),*)val
        nco=nint(val)
        if(ih.eq.1)write(ifhi,'(a,i3)')'array',nco
        if(ioint.eq.0)then
         sum=0
         averg=0
         do k=1,nrbins
          if(iocontr.eq.0.and.ionoerr.eq.0)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.1)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.2)then
            ar3=ar(k,3)
            ar4=ar(k,4)
            ar5=ar(k,5)
          else
            ar3=ary(k,ncontr)
            ar4=ardy(k,ncontr)
          endif
          iok=1
          if(k.lt.ibmin.or.k.gt.ibmax)iok=0
          sum=sum+ar3
          averg=averg+ar(k,1)*ar3
          if(nco.eq.2)then
            if(nozero.eq.1.and.ar3.eq.0.)iok=0
            if(iok.eq.1)write(ifhi,'(3e12.4)')ar(k,1),ar3
          elseif(nco.eq.3)then
            if(nozero.eq.1.and.ar3.eq.0..and.ar4.eq.0.)iok=0
            if(iok.eq.1)write(ifhi,'(3e12.4)')ar(k,1),ar3,ar4
          elseif(nco.eq.4)then
          if(nozero.eq.1.and.ar3.eq.0..and.ar4.eq.0..and.ar5.eq.0.)iok=0
            if(iok.eq.1)write(ifhi,'(4e12.4)')ar(k,1),ar3,ar4,ar5
          endif
         enddo
         if(sum.gt.0.)averg=averg/sum
        else
         sum=0.
         sum2=0.
         sum3=0.
         err2=0.
         do k=1,nrbins
          if(iocontr.eq.0.and.ionoerr.eq.0)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.1)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.2)then
            ar3=ar(k,3)
            ar4=ar(k,4)
            ar5=ar(k,5)
          else
            ar3=ary(k,ncontr)
            ar4=ardy(k,ncontr)
          endif
          sum=sum+ar3*(ar(2,1)-ar(1,1))
          if(nco.eq.2)write(ifhi,'(3e12.4)')ar(k,1),sum
          if(ionoerr.eq.0)then
            err2=err2+(ar4*(ar(2,1)-ar(1,1)))**2
            if(nco.eq.3)write(ifhi,'(3e12.4)')ar(k,1),sum,sqrt(err2)
          elseif(ionoerr.eq.1)then
            sum2=sum2+(ar4*(ar(2,1)-ar(1,1)))
            if(nco.eq.3)write(ifhi,'(3e12.4)')ar(k,1),sum,sum2
          elseif(ionoerr.eq.2)then
            sum2=sum2+(ar4*(ar(2,1)-ar(1,1)))
            sum3=sum3+(ar5*(ar(2,1)-ar(1,1)))
            if(nco.eq.3)write(ifhi,'(3e12.4)')ar(k,1),sum,sum2
            if(nco.eq.4)write(ifhi,'(3e12.4)')ar(k,1),sum,sum2,sum3
          endif
         enddo
        endif
        if(ih.eq.1)write(ifhi,'(a)')'  endarray'
       endif
      else !nopen .ge. 0 -- first run
        call utword(line,i,j,0)
        if(line(i:j).eq.'s')then
          call utword(line,i,j,0)
          call utword(line,i,j,0)
        elseif(line(i:j).eq.'int')then
          call utword(line,i,j,0)
        elseif(line(i:j).eq.'contr')then
          call utword(line,i,j,0)
          call utword(line,i,j,0)
        endif
      endif
      nozero=0
      ibmin=1
      ibmax=1e8

           else

      write(ifmt,'(a,a,a)')'command "',line(i:j),'" not found'
      j=160
      stop

           endif

      i=j+1
      goto 1

      end

c-----------------------------------------------------------------------
      subroutine aseed(modus)
c-----------------------------------------------------------------------

      include 'epos.inc'
      double precision seedf
      call utpri('aseed ',ish,ishini,3)

      call ranfgt(seedf)
      if(iwseed.eq.1)then
        if(nrevt.eq.0)then
          write(ifmt,'(a,d27.16)')'seedj:',seedf
        elseif(mod(nrevt,modsho).eq.0)then
          if(modus.eq.1)
     *   write(ifmt,'(a,i10,10x,a,d27.16)')'nrevt:',nrevt,'seedf:',seedf
          if(modus.eq.2)
     *   write(ifmt,'(a,d27.16)')'seed:',seedf
        endif
      endif
      seedc=seedf

      call utprix('aseed ',ish,ishini,3)
      return
      end

c-----------------------------------------------------------------------
      subroutine aseedi
c-----------------------------------------------------------------------

      include 'epos.inc'
      double precision seedini
      call utpri('aseedi',ish,ishini,3)

      call ranfgt(seedini)

      if(ish.ge.1)write(ifmt,'(a,d26.15)')'seedi:',seedini

      call utprix('aseedi',ish,ishini,3)
      return
      end

c$$$c-----------------------------------------------------------------------
c$$$        subroutine aseed(modus)        !Flush ????
c$$$c-----------------------------------------------------------------------
c$$$
c$$$      include 'epos.inc'
c$$$      double precision seedf
c$$$      call utpri('aseed',ish,ishini,4)
c$$$
c$$$      call ranfgt(seedf)
c$$$      if(modus.eq.2)then
c$$$        write(ifmt,'(a,d26.15)')'seed:',seedf
c$$$      elseif(modus.eq.1)then
c$$$        if(mod(nrevt,modsho).eq.0)then
c$$$          write(ifmt,100)'nrevt:',nrevt,'seedf:',seedf
c$$$          call flush(ifmt)
c$$$        endif
c$$$      endif
c$$$      seedc=seedf
c$$$
c$$$  100 format(a,i10,10x,a,d26.15)
c$$$      call utprix('aseed',ish,ishini,4)
c$$$      return
c$$$      end
c$$$
c-----------------------------------------------------------------------
      subroutine astati
c-----------------------------------------------------------------------

      include 'epos.inc'
      common/geom/rmproj,rmtarg,bmax,bkmx
      common/ghecsquel/anquasiel,iquasiel

      call utpri('astati',ish,ishini,1)
      if(ish.ge.1.and.iappl.eq.1.)then
        if(abs(accept+reject).gt.1.e-5)write(ifch,'(a,f9.5)')
     *'EMS acc.rate:',accept/(accept+reject)
        if(antot.ne.0.)write(ifch,*)'initial soft,hard(%)'
     *                           ,ansf/antot*100.
     *                           ,ansh/antot*100.,' of' ,antot
        if(antotf.ne.0.)write(ifch,*)'final soft,hard(%)'
     *                           ,ansff/antotf*100.
     *                           ,anshf/antotf*100.,' of' ,antotf
        if(antotre.ne.0.)write(ifch,*)
     *                  'droplet,string(+d),reson(+d), (had)(%) '
        if(antotre.ne.0.)write(ifch,*)'     '
     *                           ,andropl/antotre*100.
     *                           ,anstrg0/antotre*100.
     *                           ,'(',anstrg1/antotre*100.,') '
     *                           ,anreso0/antotre*100.
     *                           ,'(',anreso1/antotre*100.,') '
        if(antotre.ne.0.)write(ifch,*)'     '
     *             ,' (',anghadr/antotre*100.,')',' of' ,antotre
       if(pp4ini.gt.0.)write(ifch,*)'Energy loss',(pp4ini-pp4max)/pp4ini
      write(ifch,*)'ine cross section:',sigineex
      write(ifch,*)'diffr cross section:',sigdifex
      write(ifch,*)'SD cross section:',sigsdex
c      if(model.eq.3)write(ifch,*)'quasi-elastic cross section:'
c     &,anquasiel/float(ntevt)*a*10
         endif

c$$$      call testconex(3)

      if(iprmpt.le.0)goto1000

      write(ifch,'(77a1)')('-',i=1,77)
      write(ifch,'(a)')'statistics'
      write(ifch,'(77a1)')('-',i=1,77)
      write(ifch,'(a,i6)')'nr of messages:',imsg
      write(ifch,'(a,i8)')'maximum nptl:',nptlu
      write(ifch,'(77a1)')('-',i=1,77)

1000  continue
      call utprix('astati',ish,ishini,1)

      return
      end

c-----------------------------------------------------------------------
      subroutine atitle
c-----------------------------------------------------------------------

      include 'epos.inc'

      write(ifmt,'(67a1/a1,8x,a,f5.2,5x,a,7x,a1/a,22x,a,10x,a1/67a1)')
     *('#',l=1,68),'EPOS',iversn/100.
     *,'K. WERNER, T. PIEROG, S. PORTEBOEUF.','#'
     *,'#','Contact: werner@subatech.in2p3.fr'
     *,('#',l=1,68)
      if(iversn.eq.iverso)return
      write(ifmt,'(a,8x,a,10x,a/a,5x,a,6x,a/67a1)')
     * '#','WARNING: This is a non-official beta version!!!','#'
     *,'#','Do not publish results without contacting the authors.','#'
     *,('#',l=1,67)

      return
      end

c-----------------------------------------------------------------------
      subroutine avehep
c-----------------------------------------------------------------------

      include 'epos.inc'

      call utpri('avehep',ish,ishini,4)


      call utprix('avehep',ish,ishini,4)
      end


c-----------------------------------------------------------------------
      subroutine aepos(nin)
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      double precision eppass,etpass
      common/emnpass/eppass(mamx,4),etpass(mamx,4)
      integer iutime(5)
      call utpri('aepos',ish,ishini,4)

      call timer(iutime)
      timeini=iutime(3)+float(iutime(4))/1000.

      if(ish.ge.2)then
      call alist('start event&',0,0)
      write(ifch,*)'event number:',nrevt+1
      endif

c if random sign for projectile, set it here
      if(irdmpr.ne.0.and.maproj.eq.1)
     &   idproj=idprojin*(1-2*int(rangen()+0.5d0))

c for Air target, set the target nucleus
      if(idtargin.eq.0.and.model.ne.6)then
        call getairmol(latarg,matarg)
        if(ish.ge.2)write(ifch,*)'Air Target, select (Z,A) :'
     &                           ,latarg,matarg
      endif

      if(iappl.eq.4)then
      call amicro
      if(ish.ge.2)call alist('list before int/decays&',1,nptl)
      nevt=1
      nbdky=nptl
      call bjinta(ier)
      if(ier.eq.1)stop'error in bjinta'
      if(ish.ge.2)call alist('list after int/decays&',1,nptl)
      goto1000
      endif

      ntry=0
      if(nin.le.1)bimevt=-1
      ntevt0=ntevt
 3    ntevt=ntevt0
 2    ntevt=ntevt+1
      iret=0
      ntry=ntry+1
      if(iappl.eq.1.or.iappl.eq.2)naevt=naevt+1
 5    nevt=0
      if(nrevt.eq.0)nptly=nptl
      call cleanup
      nptl=0

      minfra=mxptl
      maxfra=0
      if(iappl.eq.1.or.iappl.eq.2)then !---hadron---geometry---
      if(ntry.lt.10000.and.engy.ge.egymin)then    !if no inel scattering -> nothing !
        if(model.eq.1)then
          call emsaaa(iret)
        else
          call emsaaaModel(model,iret)
        endif
        if(iret.gt.0)goto2
      else
        ntevt=ntevt0+1
        iret=0
        nevt=1
        call conre     !define projectile and target (elastic scattering)
        call conwr     !when the MC is suppose to produce something but failed
      endif
      if(iappl.eq.2)then
        nevt=1
        goto 1000
      endif

      elseif(iappl.eq.3) then
        call bread
        if(ish.ge.2)call alist('list after reading&',1,nptl)
        goto 1000

      elseif(iappl.eq.5) then !---kinky---
         nptl=nptly
         nevt=1
         do i=1,nptl
           istptl(i)=20
         enddo

      elseif(iappl.eq.6)then !---ee---

        call timann
        nevt=1

      elseif(iappl.eq.7)then !---decay---

        call conwr
        nevt=1

      elseif(iappl.eq.8)then !---lepton---

        call psadis(iret)
        if(iret.gt.0)goto5
        nevt=1

      endif

      if(nevt.eq.0)stop'************ should not be ***************'

      if(ish.ge.2)call alist('list before fragmentation&',1,nptl)
      nptlx=nptl+1
      if(iappl.ne.2.and.iappl.ne.7.and.nevt.eq.1.and.ifrade.ne.0)then
        call gakfra(iret)
        if(iret.gt.0)goto 3
        maxfra=nptl
        if(ish.ge.2.and.model.eq.1)
     &              call alist('list after fragmentation&',nptlx,nptl)
        if(irescl.eq.1)then
          call utghost(iret)
          if(iret.gt.0)goto 3
        endif
c       nptlx=nptl+1
      endif

      if(ispherio.eq.1.and.irescl.eq.1)then
        call utrsph(iret)
        if(iret.gt.0)goto 3
      endif



      if(iappl.ne.2.and.nevt.eq.1)then
        nbdky=nptl
        call bjinta(ier)
        if(ier.eq.1)goto3
        if(iappl.eq.1.and.irescl.eq.1)then
          call utresc(iret)
          if(iret.gt.0)goto 3
        endif
        if(ish.ge.1)then
          !if(ish.ge.2.and.ifrade.ne.0)
          !&    call alist('list after int/decays&',1,nptl)
          if(iappl.eq.1)then
            numbar=0
            pp4=0.
            do j=1,nptl
              if(istptl(j).eq.0)then
             if(idptl(j).gt. 1000.and.idptl(j).lt. 10000)numbar=numbar+1
             if(idptl(j).lt.-1000.and.idptl(j).gt.-10000)numbar=numbar-1
             if((((idptl(j).eq.1120.or.idptl(j).eq.1220)
     *           .and.idproj.gt.1000).or.(iabs(idptl(j)).gt.100
     *           .and.idproj.lt.1000)).and.pptl(4,j)
     *           .gt.pp4.and.pptl(3,j).gt.0.)pp4=pptl(4,j)
              endif
            enddo
            pp4max=pp4max+pp4
            pp4ini=pp4ini+pptl(4,1)
            nvio=isign(matarg,idtarg)-numbar
            if(iabs(idproj).gt.1000)nvio=nvio+isign(maproj,idproj)
            if(ish.ge.2)write (ifch,*)'- Baryon number conservation : '
     &                  ,nvio,' -'

          endif
        endif
      endif


      if((iappl.eq.1.or.iappl.eq.2).and.nevt.eq.0)then
        if(nin.le.1)bimevt=-1
        goto2
      endif

      if(ifrade.ne.0.and.iappl.eq.2.and.
     $     idproj.eq.1120.and.idtarg.eq.1120)then
       numbar=0
       do j=1,nptl
        if(istptl(j).eq.0)then
         if(idptl(j).gt. 1000.and.idptl(j).lt. 10000)numbar=numbar+1
         if(idptl(j).lt.-1000.and.idptl(j).gt.-10000)numbar=numbar-1
        endif
       enddo
       nvio=maproj+matarg-numbar
       if(nvio.ne.0)then
        call alist('complete list&',1,nptl)
        write(6,'(//10x,a,i3//)')'ERROR: baryon number violation:',nvio
        write(6,'(10x,a//)')
     *        'a complete list has been printed into the check-file'
        stop
       endif
      endif


      ifirst=0
      if(nrevt+1.eq.1)ifirst=1
      if(jpsi.gt.0)then
        npjpsi=0
        do i=1,jpsi
          call jpsifo(npjpsi)
          call jpsian(ifirst)
        enddo
        if(ish.ge.1)call jtauan(0,0)
        if(nrevt+1.eq.nevent)call jpsihi
      endif

      if(ixtau.eq.1)call xtauev(1)

1000  continue
      if(iabs(nin).eq.iabs(ninicon))nrevt=nrevt+1
      nglacc=nglacc+nglevt
      call timer(iutime)
      timefin=iutime(3)+float(iutime(4))/1000.
      call utprix('aepos',ish,ishini,4)
      return
      end

c-----------------------------------------------------------------------
      subroutine cleanup
c-----------------------------------------------------------------------
      include 'epos.inc'
      do i=1,nptl
        do  k=1,5
          pptl(k,i)=0
        enddo
        iorptl(i)  =0
        jorptl(i)  =0
        idptl(i)   =0
        istptl(i)  =0
        tivptl(1,i)=0
        tivptl(2,i)=0
        ifrptl(1,i)=0
        ifrptl(2,i)=0
        ityptl(i)  =0
        iaaptl(i)  =0
        radptl(i)  =0
        dezptl(i)  =0
        itsptl(i)  =0
        rinptl(i)  =-9999
        do  k=1,4
          xorptl(k,i)=0
          ibptl(k,i) =0
        enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine emsaaa(iret)
c-----------------------------------------------------------------------
c  basic EMS routine to determine Pomeron configuration
c-----------------------------------------------------------------------

      include 'epos.inc'
      common/col3/ncol,kolpt

      call utpri('emsaaa',ish,ishini,4)
      if(ish.ge.3)call alist('Determine Pomeron Configuration&',0,0)

      iret=0

1     nptl=0
      call conaa(iret)
      if(iret.gt.0)goto 1001   !no interaction
      if(iappl.eq.2.and.ixgeometry.eq.1)call xGeometry(1)
      if(iappl.eq.2)goto1000
      call conre
      call conwr
      call GfunParK(iret)
      if(iret.gt.0)goto 1    !redo
      if(ionudi.eq.0
     &  .and.(maproj.ne.1.or.matarg.ne.1).and.nglevt.eq.0)goto 1001
      call integom1(iret)
      if(iret.gt.0)goto 1    !redo
      call emsaa(iret)
      if(iret.gt.0)goto 1    !redo
      if(ncol.eq.0)goto 1001 !no interaction

1000  call utprix('emsaaa',ish,ishini,4)
      return

1001  iret=1
      goto 1000

      end


c----------------------------------------------------------------------
      subroutine alist(text,n1,n2)
c----------------------------------------------------------------------
c    ior  jor  i  ifr1  ifr2     id  ist  ity      pt  m  y
c----------------------------------------------------------------------
c       ist                                     ity
c                                  light cluster ........ 19
c   ptl ...  0 1                   soft pom ............. 20-23   25(reggeon)
c   clu ... 10 11                  hard pom low mass .... 30
c   ptn ... 20 21                  proj remnant ......... 40-49
c   str ... 29                     targ remnant ......... 50-59
c   pom ... 30 31 32(virtual)      cluster .............. 60
c   rem ... 40 41                  direct photon ........ 71,72
c----------------------------------------------------------------------
      include 'epos.inc'
      common/cxyzt/xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
     *,optl(mxptl),uptl(mxptl),sptl(mxptl),rptl(mxptl,3)
      parameter(itext=40)
      character  text*40
      dimension pp(5)
      if(n1.gt.n2)return
      imax=itext+1
      do i=itext,1,-1
      if(text(i:i).eq.'&')imax=i
      enddo
      if(imax.gt.1)then
      write(ifch,'(/1x,89a1/1x,a,a,a,90a1)')
     *('#',k=1,89),'############  ',text(1:imax-1),'  '
     *,('#',k=1,74-imax)
      write(ifch,'(1x,89a1/)')('#',k=1,89)
      endif
      if(n1.eq.0.and.n2.eq.0)return
      if(imax.gt.1)then
      write(ifch,'(1x,a,a/1x,89a1)')
     *'   ior   jor     i  ifr1  ifr2       id ist ity',
     *'        pt         m         E         y'
     *,('-',k=1,89)
      endif

      do j=1,5
        pp(j)=0.
      enddo
      nqu=0
      nqd=0
      nqs=0
      do i=n1,n2
        ptptl=pptl(1,i)**2.+pptl(2,i)**2.
        if(ptptl.le.0.d0)then
          ptptl=0
        else
          ptptl=sqrt(ptptl)
        endif
        amtptl=pptl(1,i)**2.+pptl(2,i)**2.+pptl(5,i)**2.
        if(amtptl.le.0.d0)then
          if(idptl(i).lt.10000)then
            call idmass(idptl(i),amtptl)
          endif
          amtptl=sqrt(amtptl*amtptl+pptl(1,i)**2.+pptl(2,i)**2.)
        else
          amtptl=sqrt(amtptl)
        endif
        rap=0.
        if(amtptl.gt.0.d0.and.pptl(4,i).gt.0.)
     &  rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amtptl)
         write(ifch,'(1x,i6,i6,i6,i6,i6,i10,2i3,2x,4(e9.3,1x),$)')
     &          iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &    ,idptl(i),istptl(i),ityptl(i),ptptl,pptl(5,i),pptl(4,i),rap
        write(ifch,*)' '
c        if(istptl(i).ne.12)write(ifch,*)' '
c        if(istptl(i).eq.12)write(ifch,'(1x,3(e9.3,1x))')
c     &           sptl(i),sqrt(uptl(i)-xorptl(1,i)**2)
c     &            ,sqrt(optl(i)-xorptl(2,i)**2)
        if(mod(istptl(i),10).eq.0.and.n1.eq.1.and.n2.eq.nptl)then
          do j=1,4
            pp(j)=pp(j)+pptl(j,i)
          enddo
          if(istptl(i).ne.40.and.istptl(i).ne.30)then
            call idqufl(i,idptl(i),ifl1,ifl2,ifl3)
            nqu=nqu+ifl1
            nqd=nqd+ifl2
            nqs=nqs+ifl3
          endif
        endif
      enddo
      end

c----------------------------------------------------------------------
      subroutine blist(text,n1,n2)
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(itext=40)
      character  text*40
      dimension pp(5)
      if(n1.gt.n2)return
      imax=itext+1
      do i=itext,1,-1
      if(text(i:i).eq.'&')imax=i
      enddo
      if(imax.gt.1)then
      write(ifch,'(/1x,89a1/1x,a,a,a,90a1)')
     *('#',k=1,89),'#############  ',text(1:imax-1),'  '
     *,('#',k=1,74-imax)
      write(ifch,'(1x,89a1/)')('#',k=1,89)
      endif
      if(n1.eq.0.and.n2.eq.0)return
      if(imax.gt.1)then
      write(ifch,'(1x,a,a,a/1x,90a1)')
     *'   ior   jor     i  ifr1   ifr2      id ist ity',
     *'        pt      mass    energy','       rap'
     *,('-',k=1,90)
      endif

      do j=1,5
        pp(j)=0.
      enddo
      nqu=0
      nqd=0
      nqs=0
      do i=n1,n2
        amtptl=pptl(1,i)**2.+pptl(2,i)**2.+pptl(5,i)**2.
        if(amtptl.le.0.d0)then
          if(idptl(i).lt.10000)then
            call idmass(idptl(i),amtptl)
          endif
          amtptl=sqrt(amtptl*amtptl+pptl(1,i)**2.+pptl(2,i)**2.)
        else
          amtptl=sqrt(amtptl)
        endif
        pt=pptl(1,i)**2.+pptl(2,i)**2.
        if(pt.gt.0.)pt=sqrt(pt)
        rap=0.
        if(amtptl.gt.0.d0.and.pptl(4,i).gt.0.)
     &  rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amtptl)
        write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &       ,idptl(i),istptl(i),ityptl(i)
     &       ,pt,pptl(5,i),pptl(4,i),rap
 125  format (1x,i6,i6,i6,i6,i6,i10,2i3,2x,5(e9.3,1x)
     *     ,f9.2,4x,5(e8.2,1x))
        if(mod(istptl(i),10).eq.0.and.n1.eq.1.and.n2.eq.nptl)then
          do j=1,4
            pp(j)=pp(j)+pptl(j,i)
          enddo
          if(istptl(i).ne.40.and.istptl(i).ne.30)then
            call idqufl(i,idptl(i),ifl1,ifl2,ifl3)
            nqu=nqu+ifl1
            nqd=nqd+ifl2
            nqs=nqs+ifl3
          endif
        endif
      enddo
      end

c----------------------------------------------------------------------
      subroutine alistf(text)
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(itext=40)
      character  text*40
      dimension pp(5),erest(5),errp(4)
      n1=1
      n2=nptl
      imax=itext+1
      do i=itext,1,-1
      if(text(i:i).eq.'&')imax=i
      enddo
      if(imax.gt.1)then
      write(ifch,'(/1x,89a1/1x,a,a,a,90a1)')
     *('#',k=1,89),'#############  ',text(1:imax-1),'  '
     *,('#',k=1,74-imax)
      write(ifch,'(1x,75a1/)')('#',k=1,89)
      endif
      if(imax.gt.1)then
      write(ifch,'(1x,a,a,a/1x,124a1)')
     *'   ior   jor        i     ifr1   ifr2         id ist ity',
     *'            px         py         pz         p0       mass',
     *'       rap'
     *,('-',k=1,124)
      endif

      do j=1,4
        pp(j)=0.
        errp(j)=0.
      enddo
      pp(5)=0.
      do i=n1,n2
        if(istptl(i).eq.0)then
        amtptl=pptl(1,i)**2.+pptl(2,i)**2.+pptl(5,i)**2.
        if(amtptl.le.0.d0)then
          if(idptl(i).lt.10000)then
            call idmass(idptl(i),amtptl)
          endif
          amtptl=sqrt(amtptl*amtptl+pptl(1,i)**2.+pptl(2,i)**2.)
        else
          amtptl=sqrt(amtptl)
        endif
        rap=0.
        if(amtptl.gt.0.d0.and.pptl(4,i).gt.0.)
     &  rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amtptl)
        write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &       ,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5),rap
c     &,(xorptl(j,i),j=1,4)
        do j=1,4
          pp(j)=pp(j)+pptl(j,i)
        enddo
        endif
      enddo
 125  format (1x,i6,i6,3x,i6,3x,i6,i6,i12,2i4,4x,5(e10.4,1x)
     *     ,f9.2,4x,4(e8.2,1x))
 126  format (51x,5(e10.4,1x))
 128  format (51x,65('-'))
      pp(5)=(pp(4)-pp(3))*(pp(4)+pp(3))-pp(2)**2-pp(1)**2
      if(pp(5).gt.0.)then
        pp(5)=sqrt(pp(5))
      else
        pp(5)=0.
      endif
      write (ifch,128)
      write (ifch,126) (pp(i),i=1,5)
      erest(1)=0.
      erest(2)=0.
      erest(3)=maproj*pptl(3,1)+matarg*pptl(3,maproj+1)
      erest(4)=maproj*pptl(4,1)
     &        +matarg*pptl(4,maproj+1)
      erest(5)=amproj
      write (ifch,129)  (erest(j),j=1,5)
 129  format (50x,'(',5(e10.4,1x),')')
      do j=1,4
      if(abs(pp(j)).gt.0.d0)errp(j)=100.*(pp(j)-erest(j))/pp(j)
      enddo
      write (ifch,130)  (errp(j),j=1,4)
 130  format (50x,'(',3x,4(f7.2,4x),2x,'err(%))')
      end

c----------------------------------------------------------------------
      subroutine alist2(text,n1,n2,n3,n4)
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(itext=40)
      character  text*40
      if(n1.gt.n2)return
      imax=itext+1
      do i=itext,1,-1
      if(text(i:i).eq.'&')imax=i
      enddo
      write(ifch,'(1x,a,a,a)')
     *'--------------- ',text(1:imax-1),' ---------------  '
      do i=n1,n2
      write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5)
c     &,(xorptl(j,i),j=1,4)
      enddo
      write(ifch,'(1x,a)')'----->'
      do i=n3,n4
      write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5)
c     &,(xorptl(j,i),j=1,4)
      enddo
 125  format (1x,i6,i6,3x,i6,3x,i6,i6,i12,2i4,4x,5(e8.2,1x))
c     *,4x,4(e8.2,1x))
      end

c----------------------------------------------------------------------
      subroutine alistc(text,n1,n2)
c----------------------------------------------------------------------
      include 'epos.inc'
      parameter(itext=40)
      character  text*40
      if(n1.gt.n2)return
      imax=itext+1
      do i=itext,1,-1
      if(text(i:i).eq.'&')imax=i
      enddo
      if(n1.ne.n2)write(ifch,'(1x,a,a,a)')
     *'--------------- ',text(1:imax-1),' ---------------  '
      do i=n1,n2
      write(ifch,130)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5)
     &,(xorptl(j,i),j=1,4),tivptl(1,i),tivptl(2,i)
      enddo
 130  format (1x,i6,i6,3x,i6,3x,i6,i6,i12,2i4,4x,5(e8.2,1x)
     *,4x,6(e8.2,1x))
      end

c-----------------------------------------------------------------------
      subroutine sigmaint(g0,gz,sigdo)
c-----------------------------------------------------------------------
c hadron-hadron cross sections integration
c-----------------------------------------------------------------------
      common /ar3/  x1(7),a1(7)
      include 'epos.inc'
      include 'epos.incsem'
      double precision PhiExact,vvv11,vvv12,vvv21,om1intbi!,PomNbri
     *,vvv22,ww01,ww02,ww11,ww12,ww21,ww22,gz(0:3)
     *,vvv11e,vvv12e,vvv21e,vvv22e,PhiExpo


c       rs=r2had(iclpro)+r2had(icltar)+slopom*log(engy**2)
        rs=r2had(iclpro)+r2had(icltar)+max(slopom,slopoms)*log(engy**2)
     &     +gwidth*(r2had(iclpro)+r2had(icltar))
     &     +bmxdif(iclpro,icltar)/4./0.0389
        rpom=4.*.0389*rs
        e1=exp(-1.)
c        cpt=chad(iclpro)*chad(icltar)

        gz(0)=0.d0
        gz(1)=0.d0
        gz(2)=0.d0
        gz(3)=0.d0

        sigdo=0.
        do i1=1,7
        do m=1,2

          z=.5+x1(i1)*(m-1.5)
          zv1=exp(-z)
          zv2=(e1*z)
          b1=sqrt(-rpom*log(zv1))
          b2=sqrt(-rpom*log(zv2))

          vvv11=max(0.d0,PhiExpo(.5,1.d0,1.d0,engy**2,b1))
          vvv12=max(0.d0,PhiExpo(.5,1.d0,1.d0,engy**2,b2))
          vvv21=max(0.d0,PhiExpo(1.,1.d0,1.d0,engy**2,b1))
          vvv22=max(0.d0,PhiExpo(1.,1.d0,1.d0,engy**2,b2))
c          vvv11=sqrt(vvv21)
c          vvv12=sqrt(vvv21)
          if(isetcs.eq.0)then
            vvv11e=max(0.d0,PhiExact(.5,1.d0,1.d0,engy**2,b1))
            vvv12e=max(0.d0,PhiExact(.5,1.d0,1.d0,engy**2,b2))
            vvv21e=max(0.d0,PhiExact(1.,1.d0,1.d0,engy**2,b1))
            vvv22e=max(0.d0,PhiExact(1.,1.d0,1.d0,engy**2,b2))
            vvv11=0.5d0*(vvv11+min(vvv11,vvv11e))  !to be as close as possible
            vvv12=0.5d0*(vvv12+min(vvv12,vvv12e))  !than the value with
            vvv21=0.5d0*(vvv21+min(vvv21,vvv21e))  !isetcs > 0
            vvv22=0.5d0*(vvv22+min(vvv22,vvv22e))
          endif
          ww11=1.d0-vvv11
          ww12=1.d0-vvv12
          ww21=1.d0-vvv21
          ww22=1.d0-vvv22
          ww01=vvv21-2d0*vvv11+1d0
          ww02=vvv22-2d0*vvv12+1d0

          gz(0)=gz(0)+a1(i1)*(ww01+ww02/z)
          gz(1)=gz(1)+a1(i1)*(ww11+ww12/z)
          gz(2)=gz(2)+a1(i1)*(ww21+ww22/z)
          gz(3)=gz(3)+a1(i1)*(ww11*z+ww12/z*(1.-log(z)))

          phi1=vvv21
          phi2=vvv22
          phi1x=real(min(50d0,exp(om1intbi(b1,2)/dble(r2hads(iclpro)
     &                                               +r2hads(icltar)))))
          phi2x=real(min(50d0,exp(om1intbi(b2,2)/dble(r2hads(iclpro)
     &                                               +r2hads(icltar)))))
          sigdo=sigdo+a1(i1)*(phi1*(phi1x-1.)+phi2*(phi2x-1.)/z)

        enddo
        enddo
        g0=pi*rpom*10./2.                 !common factor (pi*rpom because of b->z, 10 to have mbarn and /2. because z is from -1 to 1 but we need 0 to 1.

      return
      end
c-----------------------------------------------------------------------
      subroutine xsigma
c-----------------------------------------------------------------------
c hadron-hadron and hadron-nucleus cross sections calculation
c b - impact parameter squared (in case of hadron-nucleus interaction);
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incsem'
      double precision gz(0:3),gzp(0:3),GZ0(2)
c Model 2 Common
      COMMON /Q_AREA1/  IA(2),ICZ,ICP
      COMMON /Q_AREA6/  PIQGS,BM,AM
      COMMON /Q_AREA15/ FP(5),RQ(5),CD(5)
      COMMON /Q_AREA7/  RP1
      COMMON /Q_AREA16/ CC(5)
      double precision RP1,FP,RQ,CD,PIQGS,BM,AM,CC,GDP,GDT,GDD

c Model 3 Common
      COMMON/GSECTI/ AIEL(20),AIIN(20),AIFI(20),AICA(20),ALAM,K0FLAG
      INTEGER K0FLAG
      DOUBLE PRECISION AIEL,AIIN,AIFI,AICA,ALAM
      PARAMETER (KKK=3)
      COMMON /CGCOMP/ ACOMP,ZCOMP,WCOMP,KK
      REAL           ACOMP(KKK),ZCOMP(KKK),WCOMP(KKK)

C...Total cross sections in Pythia
c      double precision SIGT
c      COMMON/PYINT7/SIGT(0:6,0:6,0:5)

c Model 5 Common
      COMMON/HIPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)

c theoretical cross sections
      sigcut=0.
      sigtot=0.
      sigela=0.
      sigine=0.
      sigtotold=0.
      sigtotf=0.
      sigelaold=0.
      sigelaf=0.
      sigineold=0.
      sigineaa=0.
      sigdif=0.
      sloela=0.
      sigsd=0.
c simulated cross sections
      sigineex=0.          !calclated in ems if isigma=1
      sigdifex=0.
      sigsdex=0.


      call utpri('xsigma',ish,ishini,4)

      if(model.eq.1)then                        !epos

        if(icltar.ne.2)stop'xsigma: only icltar=2 allowed.'

        call sigmaint(g0,gz,sigdifold)

        sigelaold=g0*gz(0)               !elastic cross section
        rs=g0*0.4091    !=g0/pi/10.*2./4./.0389
        if(gz(1).gt.0d0)sloela=2.*rs*gz(3)/gz(1)

        sigineold=g0*gz(2)               !inelastic pomerons cross-section
        sigtotold=2.*g0*gz(1)                  !tot cross-section
        sigdifold=sigdifold * g0 !xs in mb
        sigcut=sigineold-sigdifold             !cut cross section
        x=engy
        sigtotf=21*x**0.17+44.*x**(-0.8) !fit to data
        sigelaf=30.*(x-1)**(-3)+17*x**(-0.47)+0.3*alog(x)**2 !fit to data
        sigtotfp=sigtotf
        sigelafp=sigelaf
        if(iclpro.eq.1)then        !pi+p
          sigtotf=15.*x**0.15+55.*x**(-1.5)                     !fit to data
          sigelaf=20.*(x-1)**(-3)+6.*x**(-0.4)+0.2*alog(x)**2. !fit to data
        elseif(iclpro.eq.3)then    !K+p
          sigtotf=12.5*x**0.15+35.*x**(-1.5)                     !fit to data
          sigelaf=15.*(x-1)**(-3)+3.*x**(-0.4)+0.2*alog(x)**2 !fit to data
        elseif(iclpro.eq.4)then    !D+p
          sigtotf=0.!12.5*x**0.15+35.*x**(-1.5)                     !fit to data
          sigelaf=0.!15.*(x-1)**(-3)+3.*x**(-0.4)+0.2*alog(x)**2 !fit to data
        endif
        if(engy.lt.20.)then
          sigcoul=max(0.,sigtotf-sigtotold)
        else
          sigcoul=0.
        endif


        sigdif=sigdifold
        sigdelaf=max(0.,(sigelaf+sigineold-sigtotf))
c        sigdelaf=max(0.,(sigelaf-sigelaold))

        if(rexdifi(iclpro).lt.0..or.rexdifi(icltar).lt.0.)then
          sigdela=min(sigdelaf,sigdif)

c calculate rexdif for proton first (always needed)
          if(rexdifi(icltar).lt.0.)then

            if(engy.lt.min(30.,-log(-rexdifi(icltar))/0.015))then !use fit of sigela to get rexdif
              if(iclpro+icltar.ne.4)then !pi or K - p, calculate sigdif for pp
                iclprosave=iclpro
                iclpro=2
                call sigmaint(g0p,gzp,sigdifp)
                siginep=g0p*gzp(2)
                sigdifp=sigdifp * g0p
                sigdelafp=max(0.,(sigelafp+siginep-sigtotfp))
                sigdelap=min(sigdelafp,sigdifp)
                iclpro=iclprosave
              else
                sigdifp=sigdif
                sigdelafp=sigdelaf
                sigdelap=sigdela
              endif
              if(sigdifp.gt.0.)then
                ratioSig=sigdelap/sigdifp
                rexdif(icltar)=1.-sqrt(ratioSig)
                if(rexdif(icltar).ge.exp(-0.015*engy))
     &               rexdif(icltar)=exp(-0.015*engy)
                rexdif(icltar)=max(rexdif(icltar),abs(rexdifi(icltar)))
              else
                rexdif(icltar)=1.
              endif
            else
c         rexdif(icltar)=max(exp(-1.7/engy**0.3),abs(rexdifi(icltar)))  !strong reduction
c        rexdif(icltar)=max(exp(-0.33/engy**0.066),abs(rexdifi(icltar))) !moderate one (constant sig NSD)
              rexdif(icltar)=abs(rexdifi(icltar))
            endif
          else
            rexdif(icltar)=rexdifi(icltar)
          endif

         if(iclpro.ne.2)then  !pi or K rexdif knowing p rexdif
          if(rexdifi(iclpro).lt.0.)then
            if(engy.lt.min(8.,-log(-rexdifi(iclpro))/0.015)
     &         .and.sigdif.gt.0.)then !use fit of sigela to get rexdif
              ratioSig=sigdela/sigdif
              if(abs(1.-rexdif(icltar)).gt.1.e-6)then
                rexdif(iclpro)=1.-ratioSig/(1.-rexdif(icltar))
              else
                rexdif(iclpro)=abs(rexdifi(iclpro))
              endif
              if(rexdif(iclpro).ge.exp(-0.015*engy))
     &             rexdif(iclpro)=exp(-0.015*engy)
            elseif(sigdif.le.0.)then
              rexdif(iclpro)=1.
            else
c          rexdif(iclpro)=max(exp(-1.7/engy**0.3),abs(rexdifi(iclpro))) !strong reduction
c         rexdif(iclpro)=max(exp(-0.33/engy**0.066),abs(rexdifi(iclpro)))  !moderate one (constant sig NSD)
              rexdif(iclpro)=abs(rexdifi(iclpro))
            endif
           else
             rexdif(iclpro)=abs(rexdifi(iclpro))
           endif
         endif
         sigdela=(1.-rexdif(iclpro))
     &           *(1.-rexdif(icltar))  *sigdif

        else
          rexdif(iclpro)=rexdifi(iclpro)
          rexdif(icltar)=rexdifi(icltar)
          sigdela=(1.-rexdif(iclpro))
     &           *(1.-rexdif(icltar))  *sigdif
        endif


        if(rexndf.gt.0..and.iclpro.eq.2)then
          rexndi(iclpro)=rexndf*rexdif(iclpro)
        else
          rexndi(iclpro)=rexndii(iclpro)
        endif
        if(rexndf.gt.0..and.icltar.eq.2)then
          rexndi(icltar)=rexndf*rexdif(icltar)
        else
          rexndi(icltar)=rexndii(icltar)
        endif
        if(ish.ge.2)write(ifch,*)'Xsigma : rexdif/ndi=',rexdif(iclpro)
     &                                                 ,rexdif(icltar)
     &                                                 ,rexndi(iclpro)
     &                                                 ,rexndi(icltar)

        sigsd=( (1.-rexdif(icltar))*rexdif(iclpro)
     &         +(1.-rexdif(iclpro))*rexdif(icltar) )*sigdif
c        if(engy.lt.10.)sigela=max(sigelaf,sigela)
        sigela=sigelaold+sigdela+sigcoul
        sigine=sigineold-sigdela
        sigtot=sigine+sigela
        sigineaa=eposcrse(ekin,maproj,matarg,idtarg)

c        write(ifmt,*)'Rexdif',rexdif(iclpro),rexdif(icltar)

      elseif(model.eq.2)then

        g0=real(PIQGS*RP1/CD(ICZ)*AM**2*10.D0)
        CALL m2XXFZ(0.D0,GZ0)
        gz(1)=GZ0(1)
        gz(2)=GZ0(2)
        gz(3)=0d0
        sigcut=g0*gz(2)/2.               !cut pomerons cross-section
        sigtot=g0*gz(1)                  !tot cross-section
        gz(0)=sigtot-sigcut
        sigela=gz(0)*CC(ICZ)*CC(2)       !elastic cross section
c GDP - projectile diffraction cross section
        GDP=(1.D0-CC(ICZ))*CC(2)*gz(0)
c GDT - target diffraction cross section
        GDT=(1.D0-CC(2))*CC(ICZ)*gz(0)
c  GDD - double diffractive cross section
        GDD=(1.D0-CC(ICZ))*(1.D0-CC(2))*gz(0)
        sigsd=GDT+GDD
        sigdif=sigsd+GDP
        sigine=sigcut+sigdif
        if(gz(1).gt.0.)sloela=2.*rs*gz(3)/gz(1)
        sigdifold=sigtot-sigcut-sigela       !diffractive cross section
        sigineaa=crseModel(2,ekin,maproj,matarg,idtarg)

      elseif(model.eq.3)then

        call m3SIGMA(ekin,idproj,1120,1,1,sigi,sige)
        sigine=sigi
        sigela=sige
        sigcut=sigine
        sigtot=sigine+sigela
        sigdif=sigtot-sigcut-sigela       !diffractive cross section
        call m3SIGMA(ekin,idproj,idtarg,latarg,matarg,sigi,sige)
        sigineaa=sigi

c$$$      elseif(model.eq.4)then
c$$$
c$$$        sigsd=real(SIGT(0,0,2)+SIGT(0,0,3))
c$$$        sigela=real(SIGT(0,0,1))
c$$$        sigcut=real(SIGT(0,0,5))
c$$$        sigtot=real(SIGT(0,0,0))
c$$$        sigine=sigtot-sigela
c$$$        sigdif=sigtot-sigcut-sigela       !diffractive cross section

      elseif(model.eq.5)then

        sigsd=0.
        sigdif=0.
        sigcut=0.
        sigtot=HINT1(13)
        sigine=HINT1(12)
        sigela=sigtot-sigine
        sigineaa=crseModel(5,ekin,maproj,matarg,idtarg)

      elseif(model.eq.6)then                  !for Sibyll

        call m6SIGMA(iclpro,engy,stot,sela,sine,sdifr,slela,Rho)
        sigtot=stot
        sigela=sela
        sigine=sine
        sigdif=sdifr
        sloela=slela
        sigcut=sigtot-sigdif-sigela     ! cut cross section
        sigsd=sigdif/2.
        sigineaa=crseModel(6,ekin,maproj,matarg,idtarg)

      elseif(model.eq.7)then                  !for QGSJET-II

        call m7SIGMA(stot,scut,sine,slela)
        sigtot=stot
        sigcut=scut
        sigine=sine
        sloela=slela
        sigela=sigtot-sigine     ! elastic cross section
        sigdifr=sigine-sigcut
        sigsd=sigdifr
        sigineaa=crseModel(7,ekin,maproj,matarg,idtarg)

      endif

             if(isigma.eq.1)then  !===============!

      if(ish.ge.1)
     *write (ifmt,225)engy,ekin,sigtot,sigtotf,sigtotold
     *,sigine,sigtotf-sigelaf,sigineold
     *,sigela,sigelaf,sigelaold,sigcut,sloela,sigdif,sigsd
     *,sigineaa
      if(ish.ge.1)
     *write (ifch,225)engy,ekin,sigtot,sigtotf,sigtotold
     *,sigine,sigtotf-sigelaf,sigineold
     *,sigela,sigelaf,sigelaold,sigcut,sloela,sigdif,sigsd
     *,sigineaa

      if(.not.(maproj.eq.1.and.matarg.eq.1))then  !(from tabulation) for pA/AA

        sigtotaa=0.
        sigineaa=0.
        sigelaaa=0.
        if(model.eq.1)sigineaa=eposcrse(ekin,maproj,matarg,idtarg)
        if(model.eq.2)sigineaa=crseModel(2,ekin,maproj,matarg,idtarg)
        if(model.eq.3)then
          sigineaa=gheincs
          sigelaaa=0.
          do k=1,kk
            sigelaaa=sigelaaa+AIEL(k)
          enddo
          sigtotaa=sigineaa+sigelaaa
        endif
        if(model.eq.4)sigineaa=crseModel(4,ekin,maproj,matarg,idtarg)
        if(model.eq.5)sigineaa=crseModel(5,ekin,maproj,matarg,idtarg)
        if(model.eq.6)sigineaa=crseModel(6,ekin,maproj,matarg,idtarg)
        if(model.eq.7)sigineaa=crseModel(7,ekin,maproj,matarg,idtarg)
        if(ish.ge.1)
     &  write (ifmt,226)sigtotaa,sigineaa,sigelaaa
        if(ish.ge.4)write (ifch,226)sigtotaa,sigineaa,sigelaaa

             endif  !================!


      endif

 100  continue


225   format(' hadron-proton cross sections for ',f10.2,' GeV',
     *'  (ekin:',f18.3,' GeV)'/
     *4x,'total cross section:           ',f8.2,3x,f8.2,3x,f8.2/
     *4x,'inelastic cross section:       ',f8.2,3x,f8.2,3x,f8.2/
     *4x,'elastic cross section:         ',f8.2,3x,f8.2,3x,f8.2/
     *4x,'cut cross section:             ',f8.2/
     *4x,'elastic slope parameter:       ',f8.2/
     *4x,'diffr. cross section:          ',f8.2,14x,f8.2/
     *4x,'inelastic (tab) cross section: ',f8.2)
 226  format(' hadron/nucleus-hadron/nucleus cross sections'/
     *4x,'total pA/AA cross section:     ',f8.2/
     *4x,'inelastic pA/AA cross section: ',f8.2/
     *4x,'elastic pA/AA cross section:   ',f8.2)

      call utprix('xsigma',ish,ishini,4)

      return
      end

