C      Program Mult_Gen
       SUBROUTINE multgen
      implicit none

CCC   Documentation and Description:
CCC
C     This code is intended to provide a quick means of producing
C     uncorrelated simulated events for event-by-event studies,
C     detector acceptance and efficiency studies, etc.  The
C     user selects the number of events, the one-particle distribution
C     model, the Geant particles to include, the ranges in transverse
C     momentum, pseudorapidity and azimuthal angle, the mean
C     multiplicity for each particle type for the event run, the
C     mean temperature, Rapidity width, etc., and the standard deviations
C     for the event-to-event variation in the model parameters.
C     Note that these events are produced in the c.m. frame only.
C    
C     Anisotropic flow may also be simulated by introducing explicit
C     phi-dependence (azimuthal angle) in the particle distributions.  
C     The assumed model is taken from Poskanzer and Voloshin, Phys. Rev.
C     C58, 1671 (1998), Eq.(1), where we use,
C
C          E d^3N/dp^3 = (1/2*pi*pt)*[d^2N/dpt*dy]
C             * [1 + SUM(n=1,nflowterms){2*Vn*cos[n(phi-PSIr)]}]
C
C     with up to 'nflowterms' (currently set to 6, see file
C     Parameter_values.inc) Fourier components allowed.  Vn are
C     coefficients and PSIr is the reaction plane angle.
C     The algebraic signs of the Vn terms for n=odd are reversed
C     from their input values for particles with rapidity (y) < 0
C     as suggested in Poskanzer and Voloshin.
C     The flow parameters can depend on pt and rapidity (y) according
C     to the model suggested by Art Poskanzer (Feb. 2000) and as
C     defined in the Function Vn_pt_y.
C
C     The user may select either to have the same multiplicity per
C     particle type for each event or to let the multiplicity vary
C     randomly according to a Poisson distribution. In addition, an
C     overall multiplicative scale factor can be applied to each
C     particle ID's multiplicity (same factor applied to each PID).
C     This scaling can vary randomly according to a Gaussian from
C     event-to-event.  This is to simulate trigger acceptance
C     fluctuations.  Similarly the
C     parameters of the one-particle distribution models may either
C     be fixed to the same value for each event or allowed to randomly
C     vary about a specified mean with a specified standard deviation
C     and assuming a Gaussian distribution.
C
C     With respect to the reaction plane and anisotropic flow simulation,
C     the user may select among four options:
C        (1) ignore reaction plane and anisotropic flow effects; all
C            distributions will be azimuthally invariant, on average.
C        (2) assume a fixed reaction plane angle, PSIr, for all events
C            in the run.
C        (3) assume a Gaussian distribution with specified mean and
C            standard deviation for the reaction plane angles for the
C            events in the run.  PSIr is randomly determined for each
C            event.
C        (4) assume uniformly distributed, random values for the reaction  
C            plane angles from 0 to 360 deg., for each event in the run.
C
C     The user may also select the anisotropic flow parameters, Vn,
C     to either be fixed for each event, or to randomly vary from event
C     to event according to a Gaussian distribution where the user must
C     specify the mean and std. dev.  For both cases the input file must
C     list the 'nflowterms' (e.g. 6) values of the mean and Std.dev. for
C     the Vn parameters for all particle ID types included in the run.
C
C     The available list of particles has been increased to permit a
C     number of meson and baryon resonances.  For those with broad widths
C     the code samples the mass distribution for the resonance and outputs
C     the resonance mass for each instance in a special kinematic file
C     (see file unit=9, filename = 'mult_gen.kin').  The resonance shapes
C     are approximately Breit-Wigner and are specific for each resonance 
C     case.  The additional particle/resonances include: rho(+,-,0),
C     omega(0), eta', phi, J/Psi, Delta(-,0,+,++) and K*(+,-,0).  Masses
C     are sampled for the rho, omega, phi, Deltas and D*s. 
C     Refer to SUBR: Particle_prop and Particle_mass for the explicit
C     parameters, resonance shape models, and sampling ranges.
C
C     The input is from a file, named 'mult_gen.in'.  The output is
C     loaded into a file named 'mult_gen.out' which includes run
C     header information, event header information and the EVENT: and
C     TRACK: formats as in the new STAR TEXT Format for event generator
C     input to GSTAR.  A log file, 'mult_gen.log' is also written which
C     may contain error messages.  Normally this file should be empty
C     after a successful run.  These filenames can easily be changed
C     to more suitable names by the script that runs the program or
C     by hand.
C
C     The method for generating random multiplicities and model parameter
C     values involves the following steps:
C        (1) The Poisson or Gaussian distributions are computed and
C            loaded into function f().
C        (2) The distribution f(x') is integrated from xmin to x
C            and saved from x = xmin to x = xmax.  The range and mesh
C            spaces are specified by the user.
C        (3) The integral of f is normalized to unity where 
C            integral[f(x')](at x = xmin) = 0.0
C            integral[f(x')](at x = xmax) = 1.0
C        (4) A random number generator is called which delivers values
C            between 0.0 and 1.0.  
C        (5) We consider the coordinate x (from xmin to xmax) to be
C            dependent on the integral[f].  Using the random number
C            for the selected value of integral[f] the value of x
C            is obtained by interpolation.
C
C     An interpolation subroutine from Rubin Landau, Oregon State Univ.,
C     is used to do this interpolation; it involves uneven mesh point 
C     spacing.
C
C     The method for generating the particle momenta uses the
C     standard random elimination method and involves the following
C     steps:
C
C     For model_type = 1,2,3,4 which are functions of pt,y (see following):
C        (1) The y range is computed using the pseudorapidity (eta)
C            range and includes ample cushioning around the sides
C            along the eta acceptance edges.
C        (2) The transverse momentum (pt) and rapidity (y) are
C            randomly chosen within the specified ranges.
C        (3) The pseudorapidity is computed for this (pt,y) value
C            (and the mass for each pid) and checked against the
C            pseudorapidity acceptance range.
C        (4) If the pseudorapidity is within range then the one-particle
C            model distribution is calculated at this point and its ratio
C            to the maximum value throughout (pt,eta) acceptance region
C            is calculated.
C        (5) Another random number is called and if less than the ratio
C            from step#4 the particle momentum is used; if not, then 
C            another trial value of (pt,y) is obtained.
C        (6) This continues until the required multiplicity for the
C            specific event and particle type has been satisfied.
C        (7) This process is repeated for the requested number of particle
C            types and events.
C
C     For model_type = 5,6 (see following) which are input bin-by-bin
C     in pt,eta:
C        (1) The transverse momentum (pt) and pseudorapidity (eta) are 
C            randomly chosen within the specified ranges.
C        (2) The one-particle model distribution is calculated at this
C            point and its ratio to the maximum value throughout the
C            (pt,eta) region is calculated.
C        (3) Another random number is called and if less than the ratio
C            from step(2) the particle momentum is used; if not then
C            another trial value of (pt,eta) is obtained.
C        (4) This continues until the required multiplicity for the 
C            specific event and particle type has been satisfied.
C        (5) This process is repeated for the requested number of particle
C            types and events. 
C
C     Problematic parameter values are tested, bad input values are checked
C     and in some cases may be changed so that the program will not crash.
C     In some cases the code execution is stopped.
C     Some distributions and/or unusual model parameter values may cause the
C     code to hang up due to the poor performance of the "elimination"
C     method for very strongly peaked distributions.  These are tested for
C     certain problematic values and if necessary these events are aborted.
C     A message, "*** Event No.    2903 ABORTED:" for example is printed
C     in the 'mult_gen.out' file.  Temperatures .le. 0.01 GeV and rapidity
C     width parameters .le. 0.01 will cause the event to abort.
C
C     The input is described below in the 'read' statements and also in
C     the sample input file.  Some additional comments are as follows:
C
C     (1) n_events - Selected number of events in run. Can be anything
C                    .ge. 1.
C     (2) n_pid_type - Number of particle ID types to include in the
C                      particle list. e.g. pi(+) and pi(-) are counted
C                      separately.  The limit is set by parameter npid
C                      in the accompanying include file 'Parameter_values.inc'
C                      and is presently set at 20.
C     (3) model_type - equals 1,2,3,4,5 or 6 so far.  See comments in
C                      Function dNdpty to see what is calculated.
C                      The models included are:
C                    = 1, Factorized mt exponential, Gaussian rapidity model
C                    = 2, Pratt non-expanding, spherical thermal source model
C                    = 3, Bertsch non-expanding spherical thermal source model
C                    = 4, Pratt spherically expanding, thermally equilibrated
C                         source model.
C                    = 5, Factorized pt and eta distributions input bin-by-bin.
C                    = 6, Fully 2D pt,eta distributions input bin-by-bin.
C                         NOTE: model_type = 1-4 are functions of (pt,y)
C                               model_type = 5,6 are functions of (pt,eta)
C     (4) reac_plane_cntrl - Can be either 1,2,3 or 4 where:
C                          = 1 to ignore reaction plane and anisotropic flow,
C                              all distributions will be azimuthally symm.
C                          = 2 to use a fixed reaction plane angle for all
C                              events in the run.
C                          = 3 to assume a randomly varying reaction plane
C                              angle for each event as determined by a
C                              Gaussian distribution.
C                          = 4 to assume a randomly varying reaction plane
C                              for each event in the run as determined by
C                              a uniform distribution from 0 to 360 deg.
C     (5) PSIr_mean, PSIr_stdev - Reaction plane angle mean and Gaussian
C                                 std.dev. (both are in degrees) for cases
C                                 with reac_plane_cntrl = 2 (use mean value)
C                                 and 3.  Note: these are read in regardless
C                                 of the value of reac_plane_cntrl.
C     (6) MultFac_mean, MultFac_stdev - Overall multiplicity scaling factor
C                                       for all PID types; mean and std.dev.;
C                                       for trigger fluctuations event-to-evt.
C     (7) pt_cut_min,pt_cut_max - Range of transverse momentum in GeV/c.
C     (8) eta_cut_min,eta_cut_max - Pseudorapidity range
C     (9) phi_cut_min,phi_cut_max - Azimuthal angular range in degrees.
C    (10) n_stdev_mult - Number of standard deviations about the mean value
C                        of multiplicity to include in the random event-to-
C                        event selection process.  The maximum number of
C                        steps that can be covered is determined by
C                        parameter n_mult_max_steps in the accompanying
C                        include file 'Parameter_values.inc' which is
C                        presently set at 1000, but the true upper limit for
C                        this is n_mult_max_steps - 1 = 999.
C    (11) n_stdev_temp - Same, except for the "Temperature" parameter.
C    (12) n_stdev_sigma- Same, except for the rapidity width parameter.
C    (13) n_stdev_expvel - Same, except for the expansion velocity parameter.
C    (14) n_stdev_PSIr   - Same, except for the reaction plane angle
C    (15) n_stdev_Vn     - Same, except for the anisotropy coefficients, Vn.
C    (16) n_stdev_MultFac - Same, except for the multiplicity scaling factor.
C    (17) n_integ_pts - Number of mesh points to use in the random model
C                       parameter selection process.  The upper limit is
C                       set by parameter nmax_integ in the accompanying
C                       include file 'Parameter_values.inc' which is presently
C                       set at 100, but the true upper limit for n_integ_pts
C                       is nmax_integ - 1 = 99. 
C    (18) n_scan_pts  - Number of mesh points to use to scan the (pt,y)
C                       dependence of the model distributions looking for
C                       the maximum value.  The 2-D grid has
C                       n_scan_pts * n_scan_pts points; no limit to size of
C                       n_scan_pts.
C    (19) irand       - Starting random number seed.
C
C**************************************************************************
C    FOR MODEL_TYPE = 1,2,3 or 4:
C    Input the following 7 lines for each particle type; repeat these
C    set of lines n_pid_type times:
C
C         (a) gpid - Geant Particle ID code number
C         (b) mult_mean,mult_variance_control - Mean multiplicity and
C                                               variance control where:
C             mult_variance_control = 0 for no variance in multiplicity 
C             mult_variance_control = 1 to allow Poisson distribution for
C                                       particle multiplicities for all events.
C             Note that a hard limit exists for the maximum possible
C             multiplicity for a given particle type per event.  This is
C             determined by parameter factorial_max in accompanying include
C             file 'common_facfac.inc' and is presently set at 10000.
C         (c) Temp_mean, Temp_stdev - Temperature parameter mean (in GeV)
C             and standard deviation (Gaussian distribution assumed).
C         (d) sigma_mean, sigma_stdev - Rapidity distribution width (sigma)
C             parameter mean and standard deviation (Gaussian distribution
C             assumed).
C         (e) expvel_mean, expvel_stdev - S. Pratt expansion velocity
C             (in units of c) mean and standard deviation (Gaussian 
C             distribution assumed).
C         (f) Vn_mean(k);  k=1,4  - Anisotropic flow parameters, mean values
C                                   for Fourier component n=1.
C         (g) Vn_stdev(k); k=1,4  - Anisotropic flow parameters, std.dev.
C                                   values for Fourier component n=1.
C
C             Repeat the last two lines of input for remaining Fourier
C             components n=2,3...6.  Include all 6 sets of parameters
C             even if these are not used by the model for Vn(pt,y) (set
C             unused parameter means and std.dev. to 0.0).  List 4 values
C             on every line, even though for n=even the 4th quantity is
C             not used.
C
C**************************************************************************
C    FOR MODEL_TYPE = 5 input the following set of lines for each particle
C                       type; repeat these n_pid_type times.
C
C         (a) gpid - Geant Particle ID code number
C         (b) mult_mean,mult_variance_control - Mean multiplicity and
C                                               variance control where:
C             mult_variance_control = 0 for no variance in multiplicity
C             mult_variance_control = 1 to allow Poisson distribution for
C                                       particle multiplicities for all events.
C         (c) pt_start, eta_start - minimum starting values for pt, eta 
C                                   input for the bin-by-bin distributions.
C         (d) n_pt_bins, n_eta_bins - # input pt and eta bins.
C         (e) delta_pt, pt_bin - pt bin size and function value, repeat for
C                                each pt bin.
C         (f) delta_eta, eta_bin - eta bin size and function value, repeat
C                                  for each eta bin.
C         (g) Vn_mean(k);  k=1,4  - Anisotropic flow parameters, mean values
C                                   for Fourier component n=1.
C         (h) Vn_stdev(k); k=1,4  - Anisotropic flow parameters, std.dev.
C                                   values for Fourier component n=1.
C
C             Repeat the last two lines of input for remaining Fourier
C             components n=2,3...6.  Include all 6 sets of parameters
C             even if these are not used by the model for Vn(pt,y) (set
C             unused parameter means and std.dev. to 0.0).  List 4 values
C             on every line, even though for n=even the 4th quantity is
C             not used.
C
C         NOTE: The pt, eta ranges must fully include the requested ranges
C               in input #4 and 5 above; else the code execution will stop.
C
C         Also, variable bin sizes are permitted for the input distributions.
C
C         Also, this input distribution is used for all events in the run;
C         no fluctuations in this "parent" distribution are allowed from 
C         event-to-event.
C
C**************************************************************************
C    FOR MODEL_TYPE = 6 input the following set of lines for each particle
C                       type; repeat these n_pid_type times.
C
C         (a) gpid - Geant Particle ID code number
C         (b) mult_mean,mult_variance_control - Mean multiplicity and
C                                               variance control where:
C             mult_variance_control = 0 for no variance in multiplicity
C             mult_variance_control = 1 to allow Poisson distribution for
C                                       particle multiplicities for all events.
C         (c) pt_start, eta_start - minimum starting values for pt, eta
C                                   input for the bin-by-bin distributions.
C         (d) n_pt_bins, n_eta_bins - # input pt and eta bins.
C         (e) delta_pt - pt bin size, repeat for each pt bin. 
C         (f) delta_eta - eta bin size, repeat for each eta bin.
C         (g) i,j,pt_eta_bin(i,j) - read pt (index = i) and eta (index = j)
C                                   bin numbers and bin value for full 2D space.
C         (h) Vn_mean(k);  k=1,4  - Anisotropic flow parameters, mean values
C                                   for Fourier component n=1.
C         (i) Vn_stdev(k); k=1,4  - Anisotropic flow parameters, std.dev.
C                                   values for Fourier component n=1.
C
C             Repeat the last two lines of input for remaining Fourier
C             components n=2,3...6.  Include all 6 sets of parameters
C             even if these are not used by the model for Vn(pt,y) (set
C             unused parameter means and std.dev. to 0.0).  List 4 values
C             on every line, even though for n=even the 4th quantity is
C             not used.
C
C         NOTE: The pt, eta ranges must fully include the requested ranges
C               in input #4 and 5 above; else the code execution will stop.
C
C         Also, variable bin sizes are permitted for the input distributions.
C
C         Also, this input distribution is used for all events in the run;
C         no fluctuations in this "parent" distribution are allowed from
C         event-to-event.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Include 'Parameter_values.inc'
      Include 'common_facfac.inc'
      include 'common_dist_bin.inc'
      Common/track/ pout(npid,4,factorial_max)
      real*4 pout
      Common/Geant/geant_mass(200),geant_charge(200),
     1             geant_lifetime(200),geant_width(200)
      real*4 geant_mass,geant_charge,geant_lifetime,geant_width
      Common/Mass/ Mass_integ_save(npid,nmax_integ),
     1             Mass_xfunc_save(npid,nmax_integ)
      real*4 Mass_integ_save,Mass_xfunc_save

      integer n_events, n_pid_type, model_type, n_integ_pts, irand
      integer gpid(npid),mult_mean(npid),mult_variance_control(npid)
      integer i,j,k,pid,mult_min,mult_max,n_mult_steps(npid),ievent
      integer mult_event(npid),n_scan_pts,total_mult,n_vertices
      integer event_abort,status_abort, status
      integer iptbin, ietabin
      integer track_counter,start_v,stop_v
      integer reac_plane_cntrl
      integer jodd,total_mult_inc(npid)

      real*4 pt_cut_min,pt_cut_max,eta_cut_min,eta_cut_max
      real*4 y_cut_min,y_cut_max
      real*4 phi_cut_min,phi_cut_max,n_stdev_mult,n_stdev_temp
      real*4 n_stdev_sigma,n_stdev_expvel,Temp_mean(npid)
      real*4 Temp_stdev(npid),pi,rad,mult_mean_real,mult_stdev
      real*4 mult_min_real,mult_max_real,ran
      real*4 Temp_abort, sigma_abort, bin_value

      real*4 mult_integ(n_mult_max_steps),mult_xfunc(n_mult_max_steps)
      real*4 mult_integ_save(npid,n_mult_max_steps)
      real*4 mult_xfunc_save(npid,n_mult_max_steps)
      real*4 mult_event_real

      real*4 Temp_min,Temp_max,integ(nmax_integ),xfunc(nmax_integ)
      real*4 Temp_integ_save(npid,nmax_integ)
      real*4 Temp_xfunc_save(npid,nmax_integ)
      real*4 Temp_event(npid)

      real*4 sigma_stdev(npid),sigma_mean(npid),sigma_min,sigma_max
      real*4 sigma_integ_save(npid,nmax_integ)
      real*4 sigma_xfunc_save(npid,nmax_integ)
      real*4 sigma_event(npid)

      real*4 expvel_stdev(npid), expvel_mean(npid)
      real*4 expvel_min, expvel_max
      real*4 expvel_integ_save(npid,nmax_integ)
      real*4 expvel_xfunc_save(npid,nmax_integ)
      real*4 expvel_event(npid)

      real*4 masstemp,pxtemp,pytemp,pztemp,pttemp,etatemp,ytemp,phitemp

CCC   Variables associated with anisotropic flow:

      real*4 PSIr_mean, PSIr_stdev, n_stdev_PSIr, n_stdev_Vn
      real*4 PSIr_min, PSIr_max, PSIr_event
      real*4 PSIr_integ_save(nmax_integ)
      real*4 PSIr_xfunc_save(nmax_integ)
      real*4 Vn_mean(nflowterms,4,npid), Vn_stdev(nflowterms,4,npid)
      real*4 Vn_min, Vn_max
      real*4 Vn1_integ_save(nflowterms,npid,nmax_integ)
      real*4 Vn1_xfunc_save(nflowterms,npid,nmax_integ)
      real*4 Vn2_integ_save(nflowterms,npid,nmax_integ)
      real*4 Vn2_xfunc_save(nflowterms,npid,nmax_integ)
      real*4 Vn3_integ_save(nflowterms,npid,nmax_integ)
      real*4 Vn3_xfunc_save(nflowterms,npid,nmax_integ)
      real*4 Vn4_integ_save(nflowterms,npid,nmax_integ)
      real*4 Vn4_xfunc_save(nflowterms,npid,nmax_integ)
      real*4 Vn_event(nflowterms,4,npid), Vn_event_tmp(nflowterms,4)
      real*4 Vn_sum(nflowterms,npid),cosinefac(nflowterms)

CCC   Variables associated with trigger fluctuations:
      real*4 MultFac_mean, MultFac_stdev, n_stdev_MultFac
      real*4 MultFac_min, MultFac_max, MultFac_event
      real*4 MultFac_integ_save(nmax_integ)
      real*4 MultFac_xfunc_save(nmax_integ)
CCC  Open I/O Files:

      open(unit=4,type='old',access='sequential',name='mult_gen.in')
C-->  Commented for AliRoot direct COMMON block access
C      open(unit=7,type='new',access='sequential',name='mult_gen.out')
C      open(unit=8,type='new',access='sequential',name='mult_gen.log')

CCC   NOTE:
CCC   Lines throughout the code which are commented out with "C-OUT" 
CCC   can be activated to output special files (unit=9,10) with just the
CCC   mass,pt,eta,y,phi values listed for all the tracks for all events,
CCC   or the cos[n*(phi - PSI-reaction-plane)] for n=1,2...6.
CCC   These files can be directly loaded into PAW ntuples for analysis.

C-OUT      open(unit=9,type='new',access='sequential',name='mult_gen.kin')
C-OUT      open(unit=10,type='new',access='sequential',name='mult_gen.cos')
C-->  Commented for AliRoot direct COMMON block access
C      open(unit=9,type='new',access='sequential',name='mult_gen.kin')
C      open(unit=10,type='new',access='sequential',name='mult_gen.cos')

CCC   File 'mult_gen.in' is the input file for the run.
CCC   File 'mult_gen.out' is the output file in STAR New TEXT Format.
CCC   File 'mult_gen.log' is a log file for the run.
CCC   File 'mult_gen.kin', if activated, lists track kinematics for PAW ntuples.
CCC   File 'mult_gen.cos', if activated, lists the cos(n*phi) for PAW ntuples.

CCC   Initialize Arrays to Zero:

      do i = 1,npid
         gpid(i) = 0
         mult_mean(i) = 0
         mult_variance_control(i) = 0
         n_mult_steps(i) = 0
         Temp_mean(i) = 0.0
         Temp_stdev(i) = 0.0
         sigma_mean(i) = 0.0
         sigma_stdev(i) = 0.0
         expvel_mean(i) = 0.0
         expvel_stdev(i) = 0.0
         mult_event(i) = 0
         Temp_event(i) = 0.0
         sigma_event(i) = 0.0
         expvel_event(i) = 0.0
         total_mult_inc(i) = 0

         do j = 1,n_mult_max_steps
            mult_integ_save(i,j) = 0.0
            mult_xfunc_save(i,j) = 0.0
         end do

         do j = 1,nmax_integ
            Temp_integ_save(i,j) = 0.0
            Temp_xfunc_save(i,j) = 0.0
            sigma_integ_save(i,j) = 0.0
            sigma_xfunc_save(i,j) = 0.0
            expvel_integ_save(i,j) = 0.0
            expvel_xfunc_save(i,j) = 0.0
            Mass_integ_save(i,j) = 0.0
            Mass_xfunc_save(i,j) = 0.0
         end do
      end do

      do j = 1,nflowterms
         cosinefac(j) = 0.0
         do k = 1,4
            Vn_event_tmp(j,k) = 0.0
         end do
         do i = 1,npid
            Vn_sum(j,i) = 0.0
            do k = 1,4
               Vn_mean(j,k,i)  = 0.0
               Vn_stdev(j,k,i) = 0.0
               Vn_event(j,k,i) = 0.0
            end do
            do k = 1,nmax_integ
               Vn1_integ_save(j,i,k) = 0.0
               Vn1_xfunc_save(j,i,k) = 0.0
               Vn2_integ_save(j,i,k) = 0.0
               Vn2_xfunc_save(j,i,k) = 0.0
               Vn3_integ_save(j,i,k) = 0.0
               Vn3_xfunc_save(j,i,k) = 0.0
               Vn4_integ_save(j,i,k) = 0.0
               Vn4_xfunc_save(j,i,k) = 0.0
            end do
         end do
      end do

      do i = 1,nmax_integ
         PSIr_integ_save(i) = 0.0
         PSIr_xfunc_save(i) = 0.0
         MultFac_integ_save(i) = 0.0
         MultFac_xfunc_save(i) = 0.0
      end do

      do i = 1,n_mult_max_steps
         mult_integ(i) = 0.0
         mult_xfunc(i) = 0.0
      end do

      do i = 1,nmax_integ
         integ(i) = 0.0
         xfunc(i) = 0.0
      end do

      do i = 1,factorial_max
         FACLOG(i) = 0.0
      end do

      do i = 1,60
         geant_mass(i) = 0.0
      end do

      do i = 1,npid
         pt_start(i) = 0.0
         pt_stop(i)  = 0.0
         eta_start(i) = 0.0
         eta_stop(i)  = 0.0
         n_pt_bins(i) = 0
         n_eta_bins(i) = 0
         do j = 1,n_bins_max
            delta_pt(i,j)   = 0.0
            delta_eta(i,j)  = 0.0
            pt_bin(i,j)     = 0.0
            eta_bin(i,j)    = 0.0
            do k = 1,n_bins_max
               pt_eta_bin(i,j,k) = 0.0
            end do
         end do
      end do

CCC  Read Input:

      read(4,*) n_events               ! No. of events to generate
      read(4,*) n_pid_type             ! No. of Geant PID types to include
      read(4,*) model_type             ! Distribution model type (see
CCC                                    ! Function dNdpty for explanation).
      read(4,*) reac_plane_cntrl       ! Reaction plane control option (1-4)
      read(4,*) PSIr_mean,PSIr_stdev   ! Reaction plane angle mean and std.
CCC                                    ! dev., both are in degrees.
      read(4,*) MultFac_mean,MultFac_stdev ! Mult scaling factor mean,stdev.
      read(4,*) pt_cut_min,pt_cut_max  ! Min/Max pt range in GeV/c
      read(4,*) eta_cut_min,eta_cut_max ! Min/Max pseudorapidity range
      read(4,*) phi_cut_min,phi_cut_max ! Min/Max azimuthal angular range (deg)
      read(4,*) n_stdev_mult            ! No.(+/-) standard deviation range
CCC                                     ! for multiplicity
      read(4,*) n_stdev_temp            ! No.(+/-) st.dev. range for Temp.
      read(4,*) n_stdev_sigma           ! No.(+/-) st.dev. range for rapidity
CCC                                     ! width, sigma.
      read(4,*) n_stdev_expvel          ! No.(+/-) st.dev. range for expansion
CCC                                     ! velocity.
      read(4,*) n_stdev_PSIr            ! No.(+/-) st.dev. range for PSIr
      read(4,*) n_stdev_Vn              ! No.(+/-) st.dev. range for anisotropic
CCC                                     ! flow parameters Vn.
      read(4,*) n_stdev_MultFac         ! No.(+/-) st.dev. range for multipli-
CCC                                     ! city scaling factor for all PIDs.
      read(4,*) n_integ_pts             ! No. of integration mesh points to use
CCC                                     ! for random parameter fluctuations.
      read(4,*) n_scan_pts              ! No. of pt and eta mesh points to use
CCC                                     ! in scan for maximum value of dN/dpt*dy
      read(4,*) irand                   ! Random number seed; default=12345.

CCC   Check Validity and Consistency of Input Parameters so far:

      if(n_events .le. 0) n_events = 1
      if(n_pid_type .le. 0) n_pid_type = 1
      if(pt_cut_min .gt. pt_cut_max) then
C-->         write(8,40) pt_cut_min,pt_cut_max
         RETURN
      end if
      if(eta_cut_min .gt. eta_cut_max) then
C-->         write(8,41) eta_cut_min,eta_cut_max
         RETURN
      end if     
      if(phi_cut_min .gt. phi_cut_max) then
C-->         write(8,42) phi_cut_min,phi_cut_max
         RETURN
      end if
40    Format(//10x,'pt_cut_min = ',F7.4,' gt max = ',F7.4,' -STOP')
41    Format(//10x,'eta_cut_min = ',F7.4,' gt max = ',F7.4,' -STOP')
42    Format(//10x,'phi_cut_min = ',F7.4,' gt max = ',F7.4,' -STOP')
      if(n_stdev_mult   .lt. 0.0) n_stdev_mult   = 1.0
      if(n_stdev_temp   .lt. 0.0) n_stdev_temp   = 1.0
      if(n_stdev_sigma  .lt. 0.0) n_stdev_sigma  = 1.0
      if(n_stdev_expvel .lt. 0.0) n_stdev_expvel = 1.0
      if(n_stdev_PSIr   .lt. 0.0) n_stdev_PSIr   = 1.0
      if(n_stdev_Vn     .lt. 0.0) n_stdev_Vn     = 1.0
      if(n_stdev_MultFac .lt. 0.0) n_stdev_MultFac = 1.0
      if(n_integ_pts .le. 0) n_integ_pts = 10
      if(n_scan_pts  .le. 0) n_scan_pts  = 10

      if(irand .le. 0) irand = 12345
      if(n_pid_type .gt. npid) then
C-->         write(8,10) n_pid_type, npid
10       Format(//10x,'No. requested PID types = ',I7,
     1   ', exceeds maximum of ',I7,'; reset')
         n_pid_type = npid
      end if
      if(model_type .lt. 0 .or. model_type .gt. 6) then
C-->         write(8,11) model_type
11       Format(/10x,'model_type = ',I5,' is not allowed; STOP')
         RETURN
      end if
      if(n_integ_pts .gt. nmax_integ) then
C-->         write(8,12) n_integ_pts, nmax_integ
12       Format(/10x,'No. integ. pts = ',I7,
     1   ', exceeds maximum of ',I7,'; reset')
         n_integ_pts = nmax_integ
      end if
      if(reac_plane_cntrl .lt. 1 .OR. reac_plane_cntrl .gt. 4) then
C-->         write(8,*) 'reac_plane_cntrl = ',reac_plane_cntrl,'-STOP'
         RETURN
      end if

CCC   Force the reaction plane angle (mean value) to be within the 
CCC   range 0 to 360 deg.
      PSIr_mean = PSIr_mean - 360.0*float(int(PSIr_mean/360.0))
      if(PSIr_mean .lt. 0.0) PSIr_mean = PSIr_mean + 360.0
      if(PSIr_stdev .lt. 0.0) then
C-->         write(8,*) 'Reaction plane angle Std. dev. is < 0.0 - STOP'
         RETURN
      end if

CCC   Check the multiplicity scaling factor input parameters:
      if(MultFac_mean.lt.0.0 .or. MultFac_stdev.lt.0.0) then
C-->         write(8,48) MultFac_mean, MultFac_stdev
48       Format(//10x,'Multiplicity scaling factor mean or stdev = ',
     1   2F7.4,' - not valid - STOP')
         RETURN
      end if

CCC   FOR MODEL_TYPE = 1,2,3 or 4; 
CCC   Repeat the following lines of input for each particle ID type:
     
      IF(model_type.ge.1 .and. model_type.le.4) then

      do pid = 1,n_pid_type

         read(4,*) gpid(pid)            ! Geant Particle ID code number

         read(4,*) mult_mean(pid), mult_variance_control(pid)
CCC      Mean multiplicity and variance control where:
CCC           mult_variance_control = 0 for no variance in multiplicity
CCC           mult_variance_control = 1 to allow Poisson distribution for 
CCC                                     particle multiplicities for all events.

         read(4,*) Temp_mean(pid), Temp_stdev(pid)
CCC      Temperature parameter mean (in GeV) and standard deviation (Gaussian
CCC      distribution assumed).

         read(4,*) sigma_mean(pid), sigma_stdev(pid)
CCC      Rapidity distribution width (sigma) parameter mean and standard
CCC      deviation (Gaussian distribution assumed).

         read(4,*) expvel_mean(pid), expvel_stdev(pid)
CCC      S. Pratt expansion velocity (in units of c) mean and standard
CCC      deviation (Gaussian distribution assumed).

         do i = 1,nflowterms
            read(4,*) (Vn_mean(i,k,pid) ,k=1,4)
            read(4,*) (Vn_stdev(i,k,pid),k=1,4)
         end do
CCC      Anisotropic flow Fourier component coefficients specifying the
CCC      azimuthal (phi angle) dependence as defined in Phys. Rev. C58,
CCC      1671 (1998), Eq.(1); 'nflowterms' number of parameters are
CCC      allowed (see file 'Parameter_values.inc' where this is currently
CCC      set to 6); these are the mean and Gaussian std. dev. to be used
CCC      if random, Gaussian variation in the Vn values from event-to-event
CCC      are selected (via reac_plane_cntrl).
CCC      The flow parameters are pt and y dependent; see Function Vn_pt_y
CCC      for the full definitions.

CCC      Check Validity and Consistency of Input Parameters

         if(Temp_mean(pid).le.0.0 .or. Temp_stdev(pid).lt.0.0) then
C-->         write(8,45) pid,Temp_mean(pid),Temp_stdev(pid)
45       Format(//10x,'For pid # ',I7,', Temp_mean or Temp_stdev= ',
     1   2F7.4,' - not valid -STOP')
         RETURN
         end if
         if(sigma_mean(pid).le.0.0 .or. sigma_stdev(pid).lt.0.0) then
C-->         write(8,46) pid,sigma_mean(pid),sigma_stdev(pid)
46       Format(//10x,'For pid # ',I7,', sigma_mean or sigma_stdev= ',
     1   2F7.4,' - not valid -STOP')
         RETURN
         end if
         if(expvel_mean(pid).lt.0.0.or.expvel_stdev(pid).lt.0.0) then
C-->         write(8,47) pid,expvel_mean(pid),expvel_stdev(pid)
47       Format(//10x,'For pid #',I7,',expvel_mean or expvel_stdev=',
     1   2F7.4,' - not valid -STOP')
         RETURN
         end if

         do k = 1,4
            do i = 1,nflowterms
               if(Vn_stdev(i,k,pid) .lt. 0.0) then
C-->                  write(8,*) 'For pid#,term#,k= ',pid,i,k
C-->                  write(8,*) 'Vn_stdev = ',Vn_stdev(i,k,pid)
C-->                  write(8,*) 'which is not valid, must be => 0, STOP'
                  RETURN
               end if
            end do
         end do

      end do  !  End of loop over Particle ID types.

      ELSE IF (model_type .eq. 5) then

CCC      Input for Factorized pt, eta bin-by-bin distribution:

         do pid = 1,n_pid_type
            read(4,*) gpid(pid)
            read(4,*) mult_mean(pid), mult_variance_control(pid)
            read(4,*) pt_start(pid), eta_start(pid)
            read(4,*) n_pt_bins(pid), n_eta_bins(pid)
            do i = 1,n_pt_bins(pid)
               read(4,*) delta_pt(pid,i), pt_bin(pid,i)
            end do
            do i = 1,n_eta_bins(pid)
               read(4,*) delta_eta(pid,i), eta_bin(pid,i)
            end do

         do i = 1,nflowterms
            read(4,*) (Vn_mean(i,k,pid) ,k=1,4)
            read(4,*) (Vn_stdev(i,k,pid),k=1,4)
         end do
CCC      Anisotropic flow Fourier component coefficients specifying the
CCC      azimuthal (phi angle) dependence as defined in Phys. Rev. C58,
CCC      1671 (1998), Eq.(1); 'nflowterms' number of parameters are
CCC      allowed (see file 'Parameter_values.inc' where this is currently
CCC      set to 6); these are the mean and Gaussian std. dev. to be used
CCC      if random, Gaussian variation in the Vn values from event-to-event
CCC      are selected (via reac_plane_cntrl).
CCC      The flow parameters are pt and y dependent; see Function Vn_pt_y
CCC      for the full definitions.

CCC      Evaluate grid and find maximum values of pt and eta for input bins:

            pt_stop(pid) = pt_start(pid)
            do i = 1,n_pt_bins(pid)
            pt_stop(pid) = pt_stop(pid) + delta_pt(pid,i)
            pt_bin_mesh(pid,i) = pt_stop(pid)
            end do
            eta_stop(pid) = eta_start(pid)
            do i = 1,n_eta_bins(pid)
            eta_stop(pid) = eta_stop(pid) + delta_eta(pid,i)
            eta_bin_mesh(pid,i) = eta_stop(pid)
            end do

CCC      Check ranges of pt and eta coverage:

            if(pt_cut_min .lt. pt_start(pid)) then
C-->               write(8,50) pt_cut_min,pt_start(pid)
50             Format(//10x,'pt_cut_min = ',F10.7,' .lt. pt_start =',
     1         F10.7,' - STOP')
               RETURN
            end if
            if(pt_cut_max .gt. pt_stop(pid)) then
C-->               write(8,51) pt_cut_max,pt_stop(pid)
51             Format(//10x,'pt_cut_max = ',F10.7,' .gt. pt_stop =',
     1         F10.7,' - STOP')
               RETURN
            end if
            if(eta_cut_min .lt. eta_start(pid)) then
C-->               write(8,52) eta_cut_min,eta_start(pid)
52             Format(//10x,'eta_cut_min = ',F10.7,'.lt.eta_start =',
     1         F10.7,' - STOP')
               RETURN
            end if
            if(eta_cut_max .gt. eta_stop(pid)) then
C-->               write(8,53) eta_cut_max,eta_stop(pid)
53             Format(//10x,'eta_cut_max = ',F10.7,'.gt.eta_stop =',
     1         F10.7,' - STOP')
               RETURN
            end if

         do k = 1,4
            do i = 1,nflowterms
               if(Vn_stdev(i,k,pid) .lt. 0.0) then
C-->                  write(8,*) 'For pid#,term#,k= ',pid,i,k
C-->                  write(8,*) 'Vn_stdev = ',Vn_stdev(i,k,pid)
C-->                  write(8,*) 'which is not valid, must be => 0, STOP'
                  RETURN
               end if
            end do
         end do

         end do ! End loop over particle ID types.

      ELSE IF (model_type .eq. 6) then

CCC      Input for Full 2D (pt,eta) bin-by-bin distribution:

         do pid = 1,n_pid_type
            read(4,*) gpid(pid)
            read(4,*) mult_mean(pid), mult_variance_control(pid)
            read(4,*) pt_start(pid), eta_start(pid)
            read(4,*) n_pt_bins(pid), n_eta_bins(pid)
            do i = 1,n_pt_bins(pid)
               read(4,*) delta_pt(pid,i)
            end do
            do i = 1,n_eta_bins(pid)
               read(4,*) delta_eta(pid,i)
            end do

CCC      Evaluate grid and find maximum values of pt and eta for input bins:

            pt_stop(pid) = pt_start(pid)
            do i = 1,n_pt_bins(pid)
            pt_stop(pid) = pt_stop(pid) + delta_pt(pid,i)
            pt_bin_mesh(pid,i) = pt_stop(pid)
            end do
            eta_stop(pid) = eta_start(pid)
            do i = 1,n_eta_bins(pid)
            eta_stop(pid) = eta_stop(pid) + delta_eta(pid,i)
            eta_bin_mesh(pid,i) = eta_stop(pid)
            end do

CCC   Load 2D bin-by-bin distribution:

            do i = 1,n_pt_bins(pid)*n_eta_bins(pid)
               read(4,*) iptbin,ietabin,bin_value
               pt_eta_bin(pid,iptbin,ietabin) = bin_value
            end do

         do i = 1,nflowterms
            read(4,*) (Vn_mean(i,k,pid) ,k=1,4)
            read(4,*) (Vn_stdev(i,k,pid),k=1,4)
         end do
CCC      Anisotropic flow Fourier component coefficients specifying the
CCC      azimuthal (phi angle) dependence as defined in Phys. Rev. C58,
CCC      1671 (1998), Eq.(1); 'nflowterms' number of parameters are
CCC      allowed (see file 'Parameter_values.inc' where this is currently
CCC      set to 6); these are the mean and Gaussian std. dev. to be used
CCC      if random, Gaussian variation in the Vn values from event-to-event
CCC      are selected (via reac_plane_cntrl).
CCC      The flow parameters are pt and y dependent; see Function Vn_pt_y
CCC      for the full definitions.

CCC      Check ranges of pt and eta coverage:

            if(pt_cut_min .lt. pt_start(pid)) then
C-->               write(8,50) pt_cut_min,pt_start(pid)
               RETURN
            end if
            if(pt_cut_max .gt. pt_stop(pid)) then
C-->               write(8,51) pt_cut_max,pt_stop(pid)
               RETURN
            end if
            if(eta_cut_min .lt. eta_start(pid)) then
C-->               write(8,52) eta_cut_min,eta_start(pid)
               RETURN
            end if
            if(eta_cut_max .gt. eta_stop(pid)) then
C-->               write(8,53) eta_cut_max,eta_stop(pid)
               RETURN
            end if
         do k = 1,4
            do i = 1,nflowterms
               if(Vn_stdev(i,k,pid) .lt. 0.0) then
C-->                  write(8,*) 'For pid#,term#,k= ',pid,i,k
C-->                  write(8,*) 'Vn_stdev = ',Vn_stdev(i,k,pid)
C-->                  write(8,*) 'which is not valid, must be => 0, STOP'
                  RETURN
               end if
            end do
         end do

         end do ! End loop over particle ID types.

      END IF !  End of MODEL_TYPE Options input:
      
CCC   Check some more input parameters; Set constants, and load data tables:

      do pid = 1,n_pid_type
      if(gpid(pid).le.0 .or. gpid(pid).gt.200) then
C-->      write(8,43) pid,gpid(pid)
43    Format(//10x,'For pid # ',I7,', gpid = ',I7,' - not valid -STOP')
      RETURN
      end if
      if(mult_variance_control(pid) .lt. 0  .or.
     1   mult_variance_control(pid) .gt. 1)
     2   mult_variance_control(pid) = 0
      if(mult_mean(pid) .le. 0) then
C-->      write(8,44) pid,mult_mean(pid)
44    Format(//10x,'For pid # ',I7,', mult_mean= ',I7,
     1   ' - not valid -STOP')
      RETURN
      end if
      end do ! End Particle ID input parameter check and verification loop.

      pi = 3.141592654
      rad = 180.0/pi
      Temp_abort  = 0.01
      sigma_abort = 0.01

CCC   Load particle properties array; mass in GeV, etc.:

      Call Particle_prop

CCC   Load log(n!) values into Common/FACFAC/

      Call FACTORIAL

CCC   Set y (rapidity) range, to be used for model_type = 1-4
      if(eta_cut_min .ge. 0.0) then
         y_cut_min = 0.0
      else
         y_cut_min = eta_cut_min
      end if
      if(eta_cut_max .le. 0.0) then
         y_cut_max = 0.0
      else
         y_cut_max = eta_cut_max
      end if

CCC   Obtain integrals of 1D distribution functions of multiplicity
CCC   (Poisson, integer) and parameters (Gaussian, real*4) for the
CCC   particle model distributions, for each particle ID type.
CCC   These integrated quantities are then used with random number
CCC   selection to generate a distribution of multiplicities and
CCC   parameter values for each event in the run.

CCC   Obtain 1D integrals for Gaussian distributions for reaction plane
CCC   angles:

      if(reac_plane_cntrl .eq. 3) then
         if((n_stdev_PSIr*PSIr_stdev).gt. 0.0) then
            PSIr_min = PSIr_mean - n_stdev_PSIr*PSIr_stdev
            PSIr_max = PSIr_mean + n_stdev_PSIr*PSIr_stdev
            Call Gaussian(PSIr_min,PSIr_max,PSIr_mean,PSIr_stdev,
     1           n_integ_pts,integ,xfunc,nmax_integ)
            do i = 1,n_integ_pts + 1
               PSIr_integ_save(i) = integ(i)
               PSIr_xfunc_save(i) = xfunc(i)
            end do
         else
            PSIr_event = PSIr_mean
         end if
      end if
 
CCC   Obtain 1D integral for Gaussian distribution for the trigger fluctuation
CCC   simulations via the overall multiplicity scaling factor.
      if((n_stdev_MultFac*MultFac_stdev).gt.0.0) then
         Call MinMax(MultFac_mean,MultFac_stdev,n_stdev_MultFac,
     1               MultFac_min,MultFac_max)
         Call Gaussian(MultFac_min,MultFac_max,MultFac_mean,
     1                 MultFac_stdev,n_integ_pts,integ,xfunc,nmax_integ)
         do i = 1,n_integ_pts + 1
            MultFac_integ_save(i) = integ(i)
            MultFac_xfunc_save(i) = xfunc(i)
         end do
      else
         MultFac_event = MultFac_mean
      end if

      do pid = 1,n_pid_type

         Call Particle_mass(gpid(pid),pid,n_integ_pts)

         if(mult_variance_control(pid) .ne. 0) then
            mult_mean_real = float(mult_mean(pid))
            mult_stdev = sqrt(mult_mean_real)
            Call MinMax(mult_mean_real,mult_stdev,n_stdev_mult,
     1                  mult_min_real,mult_max_real)
            mult_min = int(mult_min_real)
            mult_max = int(mult_max_real + 1)
            n_mult_steps(pid) = mult_max - mult_min + 1
            if((n_mult_steps(pid) + 1) .gt. n_mult_max_steps) then
C-->               write(8,20) n_mult_steps(pid),n_mult_max_steps
20             Format(//10x,'No. steps in multiplicity integral = ',
     1         I7,' + 1, exceeds max no. of ',I7,'; reset')
            mult_min = mult_mean(pid) - (n_mult_max_steps - 1)/2
            mult_max = mult_mean(pid) + (n_mult_max_steps - 1)/2
            n_mult_steps(pid) = mult_max - mult_min + 1
            end if
            if((mult_max + 1) .gt. factorial_max) then
C-->               write(8,30) mult_max, factorial_max
30             Format(//10x,'In Poisson multiplicity distribution,',
     1         ' max = ',I7,', exceeds limit of ',I7,' - STOP')
               RETURN
            end if

            Call Poisson(mult_min,mult_max,mult_mean(pid),
     1      n_mult_steps(pid),mult_integ,mult_xfunc,n_mult_max_steps)
            do i = 1,n_mult_steps(pid) + 1
               mult_integ_save(pid,i) = mult_integ(i)
               mult_xfunc_save(pid,i) = mult_xfunc(i)
            end do
         else if (mult_variance_control(pid) .eq. 0) then
            mult_event(pid) = mult_mean(pid)
         end if

         if(model_type .le. 4) then
         if(Temp_stdev(pid) .ne. 0.0) then
            Call MinMax(Temp_mean(pid),Temp_stdev(pid),n_stdev_temp,
     1                  Temp_min, Temp_max)
            Call Gaussian(Temp_min,Temp_max,Temp_mean(pid),
     1         Temp_stdev(pid),n_integ_pts,integ,xfunc,nmax_integ)
            do i = 1,n_integ_pts + 1
               Temp_integ_save(pid,i) = integ(i)
               Temp_xfunc_save(pid,i) = xfunc(i)
            end do
         else if(Temp_stdev(pid) .eq. 0.0) then
            Temp_event(pid) = Temp_mean(pid)
         end if

         if(sigma_stdev(pid) .ne. 0.0) then
            Call MinMax(sigma_mean(pid),sigma_stdev(pid),n_stdev_sigma,
     1                  sigma_min, sigma_max)
            Call Gaussian(sigma_min,sigma_max,sigma_mean(pid),
     1         sigma_stdev(pid),n_integ_pts,integ,xfunc,nmax_integ)
            do i = 1,n_integ_pts + 1
               sigma_integ_save(pid,i) = integ(i)
               sigma_xfunc_save(pid,i) = xfunc(i)
            end do
         else if(sigma_stdev(pid) .eq. 0.0) then
            sigma_event(pid) = sigma_mean(pid)
         end if

         if(expvel_stdev(pid) .ne. 0.0) then
          Call MinMax(expvel_mean(pid),expvel_stdev(pid),n_stdev_expvel,
     1                  expvel_min, expvel_max)
            Call Gaussian(expvel_min,expvel_max,expvel_mean(pid),
     1         expvel_stdev(pid),n_integ_pts,integ,xfunc,nmax_integ)
            do i = 1,n_integ_pts + 1
               expvel_integ_save(pid,i) = integ(i)
               expvel_xfunc_save(pid,i) = xfunc(i)
            end do
         else if(expvel_stdev(pid) .eq. 0.0) then
            expvel_event(pid) = expvel_mean(pid)
         end if 
         end if ! End model_type .le. 4 options.

         if(reac_plane_cntrl .gt. 1) then
            do i = 1,nflowterms
             do k = 1,4
               if((n_stdev_Vn*Vn_stdev(i,k,pid)) .gt. 0.0) then
                  Vn_min = Vn_mean(i,k,pid)-n_stdev_Vn*Vn_stdev(i,k,pid)
                  Vn_max = Vn_mean(i,k,pid)+n_stdev_Vn*Vn_stdev(i,k,pid)
                  Call Gaussian(Vn_min,Vn_max,Vn_mean(i,k,pid),
     1            Vn_stdev(i,k,pid),n_integ_pts,integ,xfunc,nmax_integ)
                  if(k.eq.1) then
                  do j = 1,n_integ_pts + 1
                     Vn1_integ_save(i,pid,j) = integ(j)
                     Vn1_xfunc_save(i,pid,j) = xfunc(j)
                  end do
                  else if(k.eq.2) then
                  do j = 1,n_integ_pts + 1
                     Vn2_integ_save(i,pid,j) = integ(j)
                     Vn2_xfunc_save(i,pid,j) = xfunc(j)
                  end do
                  else if(k.eq.3) then
                  do j = 1,n_integ_pts + 1
                     Vn3_integ_save(i,pid,j) = integ(j)
                     Vn3_xfunc_save(i,pid,j) = xfunc(j)
                  end do
                  else if(k.eq.4) then
                  do j = 1,n_integ_pts + 1
                     Vn4_integ_save(i,pid,j) = integ(j)
                     Vn4_xfunc_save(i,pid,j) = xfunc(j)
                  end do
                  end if
               else
                  Vn_event(i,k,pid) = Vn_mean(i,k,pid)
               end if
             end do
            end do
         end if

      end do  !  End of PID Loop:

CCC   Write Run Header Output:

C-->      write(7,200)
C-->      write(7,201) n_events,n_pid_type
C-->      if(model_type .eq. 1) write(7,202)
C-->      if(model_type .eq. 2) write(7,203)
C-->      if(model_type .eq. 3) write(7,204)
C-->      if(model_type .eq. 4) write(7,205)
C-->      if(model_type .eq. 5) write(7,2051)
C-->      if(model_type .eq. 6) write(7,2052)
C-->      write(7,2053) reac_plane_cntrl
C-->      write(7,2054) PSIr_mean, PSIr_stdev
C-->      write(7,2055) MultFac_mean,MultFac_stdev
C-->      write(7,206) pt_cut_min, pt_cut_max
C-->      write(7,207) eta_cut_min, eta_cut_max
C-->      write(7,2071) y_cut_min,y_cut_max
C-->      write(7,208) phi_cut_min, phi_cut_max
C-->      write(7,209) n_stdev_mult,n_stdev_temp,n_stdev_sigma,
C-->     1             n_stdev_expvel
C-->      write(7,2091) n_stdev_PSIr, n_stdev_Vn
C-->      write(7,2092) n_stdev_MultFac
C-->      write(7,210) n_integ_pts
C-->      write(7,211) n_scan_pts
C-->      write(7,212) irand
C-->      write(7,213) (gpid(pid),pid = 1,n_pid_type)
C-->      write(7,214) (mult_mean(pid),pid = 1,n_pid_type)
C-->      write(7,215) (mult_variance_control(pid),pid = 1,n_pid_type)
C-->      if(model_type .le. 4) then
C-->         write(7,216) (Temp_mean(pid),pid = 1,n_pid_type)
C-->         write(7,217) (Temp_stdev(pid),pid = 1,n_pid_type)
C-->         write(7,218) (sigma_mean(pid),pid = 1,n_pid_type)
C-->         write(7,219) (sigma_stdev(pid),pid = 1,n_pid_type)
C-->         write(7,220) (expvel_mean(pid),pid = 1,n_pid_type)
C-->         write(7,221) (expvel_stdev(pid),pid = 1,n_pid_type)
C-->      else if(model_type .ge. 5) then
C-->         write(7,222) (n_pt_bins(pid),pid = 1,n_pid_type)
C-->         write(7,223) (n_eta_bins(pid),pid = 1,n_pid_type)
C-->         write(7,224) (pt_start(pid),pid = 1,n_pid_type)
C-->         write(7,225) (pt_stop(pid),pid = 1,n_pid_type)
C-->         write(7,226) (eta_start(pid),pid = 1,n_pid_type)
C-->         write(7,227) (eta_stop(pid),pid = 1,n_pid_type)
C-->      end if
CCC   Print out flow parameters:
      do pid = 1,n_pid_type
         do i = 1,nflowterms
C-->            write(7,2271) pid,i,(Vn_mean(i,k,pid),k=1,4)
C-->            write(7,2272) pid,i,(Vn_stdev(i,k,pid),k=1,4)
         end do
      end do

200   Format('***  Multiplicity Generator Run Header ***')
201   Format('*    Number of events = ',I7,'; number of Particle ID',
     1       ' types = ',I5)
202   Format('* Factorized mt exponential, Gaussian rapidity model')
203   Format('* Pratt non-expanding, spherical thermal source model')
204   Format('* Bertsch non-expanding spherical thermal source model')
205   Format('* Pratt spherically expanding, thermally equilibrated ',
     1        'source model')
2051  Format('* Factorized pt,eta bin-by-bin Distribution')
2052  Format('* Full 2D pt,eta bin-by-bin Distribution')
2053  Format('* Reaction plane control = ',I5)
2054  Format('* Reaction plane angles - mean and std.dev. = ',2G12.5,
     1          ' (deg)')
2055  Format('* Multiplicity Scaling Factor - mean and std.dev = ',
     1       2G12.5)
206   Format('* Min, Max range in pt              = ', 2G12.5)
207   Format('* Min, Max range in pseudorapidity  = ', 2G12.5)
2071  Format('* Min, Max range in rapdity + cush  = ', 2G12.5)
208   Format('* Min, Max range in azimuthal angle = ', 2G12.5)
209   Format('* No. std. dev. range used for mult and parameters = ',
     1       4F5.2)
2091  Format('* No. std. dev. range for PSIr, Vn  = ', 2G12.5)
2092  Format('* No. std. dev. range for MultFac   = ',G12.5)
210   Format('* No. integration points for parameter variance = ',
     1       I6)
211   Format('* No. of dN/dp(pt,y) scanning grid points to find ',
     1   'maximum = ', I6)
212   Format('* Starting Random Number Seed = ',I10)
213   Format('* Geant PID:          ',10I7)
214   Format('* Mean Multiplicity:  ',10I7)
215   Format('* Mult. Var. Control: ',10I7)
216   Format('* Mean Temperature:   ',10F7.4)
217   Format('* Std. Dev. Temp:     ',10F7.4)
218   Format('* Mean Rap. sigma:    ',10F7.4)
219   Format('* Std. Dev. y-sigma:  ',10F7.4)
220   Format('* Mean expansion vel: ',10F7.4)
221   Format('* Std. Dev. exp. vel: ',10F7.4)
222   Format('* No. input pt bins:  ',10I7)
223   Format('* No. input eta bins: ',10I7)
224   Format('* Input pt start:     ',10F7.4)
225   Format('* Input pt stop:      ',10F7.4)
226   Format('* Input eta start:    ',10F7.4)
227   Format('* Input eta stop:     ',10F7.4)
2271  Format('* Anisotropy Vn mean (pid=',I2,',n=',I1,'):',4G12.5)
2272  Format('* Anisotropy Vn std. (pid=',I2,',n=',I1,'):',4G12.5)

CCC  END Run Header Output

C**************************************************************
C                                                            **
C                     START EVENT LOOP                       **
C                                                            **
C**************************************************************

      DO ievent = 1,n_events

CCC   Compute the Reaction plane angle for this event:
         if(reac_plane_cntrl .eq. 1) then
            PSIr_event = 0.0
         else if(reac_plane_cntrl .eq. 2) then
            PSIr_event = PSIr_mean
         else if(reac_plane_cntrl .eq. 3) then
            if((n_stdev_PSIr*PSIr_stdev) .gt. 0.0) then
               do i = 1,n_integ_pts + 1
                  integ(i) = PSIr_integ_save(i)
                  xfunc(i) = PSIr_xfunc_save(i)
               end do
               Call LAGRNG(ran(irand),integ,PSIr_event,xfunc,
     1              n_integ_pts+1,1,5,n_integ_pts+1,1)
CCC         Ensure that the randomly selected reaction plane angle
CCC         is within the 0 to 360 deg range.
               PSIr_event=PSIr_event-360.0*float(int(PSIr_event/360.0))
               if(PSIr_event .lt. 0.0) PSIr_event = PSIr_event + 360.0
            end if
CCC NOTE: If PSIr_stdev=0.0, PSIr_event already calculated (see preceding)
         else if(reac_plane_cntrl .eq. 4) then
            PSIr_event = 360.0*ran(irand)
         else
            PSIr_event = 0.0
         end if

CCC   Compute the multiplicity scaling factor for simulating trigger
CCC   fluctuations for this event:
         if((n_stdev_MultFac*MultFac_stdev).gt.0.0) then
            do i = 1,n_integ_pts + 1
               integ(i) = MultFac_integ_save(i)
               xfunc(i) = MultFac_xfunc_save(i)
            end do
            Call LAGRNG(ran(irand),integ,MultFac_event,xfunc,
     1         n_integ_pts+1,1,5,n_integ_pts+1,1)
         end if

         do pid = 1,n_pid_type
            if(mult_variance_control(pid) .ne. 0) then
               do i = 1,n_mult_steps(pid) + 1
               mult_integ(i) = mult_integ_save(pid,i)
               mult_xfunc(i) = mult_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),mult_integ,mult_event_real,
     1         mult_xfunc,n_mult_steps(pid)+1,1,5,
     2         n_mult_steps(pid)+1,1)
               mult_event(pid) = mult_event_real
            else if(mult_variance_control(pid) .eq. 0) then
               mult_event(pid) = mult_mean(pid)
            end if
            mult_event(pid) = int(MultFac_event*float(mult_event(pid))
     1                            + 0.5)
CCC   Check each multiplicity wrt upper array size limit:
            if(mult_event(pid).gt.factorial_max) then
C-->               write(8,*) 'For event#',ievent,'and pid#',pid,
C-->     1            'multiplicity=',mult_event(pid),
C-->     2            'which is > array size (factorial_max=',
C-->     3             factorial_max,')-STOP'
               RETURN
            end if

         if(model_type .le. 4) then

            if(Temp_stdev(pid) .ne. 0.0) then
               do i = 1,n_integ_pts + 1
               integ(i) = Temp_integ_save(pid,i)
               xfunc(i) = Temp_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),integ,Temp_event(pid),xfunc,
     1              n_integ_pts+1,1,5,n_integ_pts+1,1)
            end if

            if(sigma_stdev(pid) .ne. 0.0) then
               do i = 1,n_integ_pts + 1
               integ(i) = sigma_integ_save(pid,i)
               xfunc(i) = sigma_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),integ,sigma_event(pid),xfunc,
     1              n_integ_pts+1,1,5,n_integ_pts+1,1)
            end if

            if(expvel_stdev(pid) .ne. 0.0) then
               do i = 1,n_integ_pts + 1
               integ(i) = expvel_integ_save(pid,i)
               xfunc(i) = expvel_xfunc_save(pid,i)
               end do
               Call LAGRNG(ran(irand),integ,expvel_event(pid),xfunc,
     1              n_integ_pts+1,1,5,n_integ_pts+1,1)
            end if
         end if
         
         if(reac_plane_cntrl .gt. 1) then
            do i = 1,nflowterms
             do k = 1,4
               if((n_stdev_Vn*Vn_stdev(i,k,pid)) .gt. 0.0) then
                if(k.eq.1) then
                  do j = 1,n_integ_pts + 1
                     integ(j) = Vn1_integ_save(i,pid,j)
                     xfunc(j) = Vn1_xfunc_save(i,pid,j)
                  end do
                else if (k.eq.2) then
                  do j = 1,n_integ_pts + 1
                     integ(j) = Vn2_integ_save(i,pid,j)
                     xfunc(j) = Vn2_xfunc_save(i,pid,j)
                  end do
                else if (k.eq.3) then
                  do j = 1,n_integ_pts + 1
                     integ(j) = Vn3_integ_save(i,pid,j)
                     xfunc(j) = Vn3_xfunc_save(i,pid,j)
                  end do
                else if (k.eq.4) then
                  do j = 1,n_integ_pts + 1
                     integ(j) = Vn4_integ_save(i,pid,j)
                     xfunc(j) = Vn4_xfunc_save(i,pid,j)
                  end do
                end if
                Call LAGRNG(ran(irand),integ,Vn_event(i,k,pid),xfunc,
     1                 n_integ_pts+1,1,5,n_integ_pts+1,1)
               end if ! (For n_stdev_Vn*Vn_stdev=0, already have Vn_event)
             end do
            end do
         end if

         end do !  End Particle ID setup Loop

         event_abort = 1

         if(model_type .le. 4) then
CCC      Check validity of Temperature and sigma parameters:
         do pid = 1,n_pid_type
            if(Temp_event(pid) .le. Temp_abort) event_abort = -86
            if(sigma_event(pid).le.sigma_abort) event_abort = -86
         end do
         end if

         if(event_abort .eq. 1) then

            total_mult = 0
            do pid = 1,n_pid_type
            total_mult = total_mult + mult_event(pid)
            end do
            n_vertices = 0
            status_abort = 1
            do pid = 1,n_pid_type
CCC   Load Anisotropic flow parameters for this event# and PID type
CCC   into temporary array;
               do i = 1,nflowterms
                  do k = 1,4
                     Vn_event_tmp(i,k) = Vn_event(i,k,pid)
                  end do
               end do
CCC   For the specified Geant PID, multiplicity, model type, temperature,
CCC   rapidity width (sigma), and expansion velocity parameter, generate
CCC   random track list:

            Call track_gen(gpid(pid),mult_event(pid),model_type,
     1      Temp_event(pid),sigma_event(pid),expvel_event(pid),
     2      reac_plane_cntrl,PSIr_event,Vn_event_tmp,
     3      pt_cut_min,pt_cut_max,eta_cut_min,eta_cut_max,
     4      y_cut_min,y_cut_max,
     5      phi_cut_min,phi_cut_max,n_scan_pts,irand,rad,pid,status,
     6      n_integ_pts)
            if(status .eq. -86) status_abort = -86
            end do
         end if

         if(event_abort.eq.1 .and. status_abort.eq.1) then
CCC Event Header Output:
C-->           write(7,230) ievent, total_mult
C-->           write(7,2301) PSIr_event
C-->           write(7,2302) MultFac_event
C-->           write(7,213) (gpid(pid),pid = 1,n_pid_type)
C-->           write(7,231) (mult_event(pid),pid = 1,n_pid_type)
C-->           if(model_type .le. 4) then
C-->              write(7,232) (Temp_event(pid),pid = 1,n_pid_type)
C-->              write(7,233) (sigma_event(pid),pid = 1,n_pid_type)
C-->              write(7,234) (expvel_event(pid),pid = 1,n_pid_type)
C-->           end if

230   Format(/'***  Next Event:  No. ',I7,'; Total # tracks = ',I10)
2301  Format('* Reaction plane angle = ',G12.5,' (deg)')
2302  Format('* Multiplicity Scaling Factor = ',G12.5)
231   Format('* Multiplicity:         ',10I7)
232   Format('* Temperature:          ',10F7.4)
233   Format('* Rapidity Dist. sigma: ',10F7.4)
234   Format('* Expansion Velocity:   ',10F7.4)

CCC   Print out flow parameters for this event:
C-->           do pid = 1,n_pid_type
C-->              do i = 1,nflowterms
C-->                 write(7,2341) pid,i,(Vn_event(i,k,pid),k=1,4)
C-->              end do
C-->           end do
2341  Format('* Anisotropy Vn (pid=',I2,',n=',I1,'):',4G12.5)

C-->           write(7,235) ievent,total_mult,n_vertices
235        Format('EVENT:',3x,3(1x,I6))
C-->           write(7,236)
C-->           write(7,237)
236        Format('*** Tracks for this event ***')
237        Format('* Geant PID #      px          py          pz')
CCC   End Event Header Output:

CCC   Output track kinematics for ievent and pid:
           track_counter = 0
           start_v = 0
           stop_v  = 0
           do pid = 1,n_pid_type
           do i = 1,mult_event(pid)
           track_counter = track_counter + 1
C-->           write(7,240) gpid(pid),(pout(pid,j,i),j=1,3),
C-->     1          track_counter,start_v,stop_v,gpid(pid)
C-OUT           masstemp = pout(pid,4,i)
C-OUT           pxtemp = pout(pid,1,i)
C-OUT           pytemp = pout(pid,2,i)
C-OUT           pztemp = pout(pid,3,i)
C-OUT           Call kinematics2(masstemp,pxtemp,pytemp,pztemp,pttemp,
C-OUT     1                      etatemp,ytemp,phitemp)
C-OUT           write(9,270) gpid(pid),masstemp,pttemp,etatemp,ytemp,phitemp
C-OUT270   Format(1x,I4,5E15.6)
           masstemp = pout(pid,4,i)
           pxtemp = pout(pid,1,i)
           pytemp = pout(pid,2,i)
           pztemp = pout(pid,3,i)
           Call kinematics2(masstemp,pxtemp,pytemp,pztemp,pttemp,
     1                      etatemp,ytemp,phitemp)
C-->           write(9,270) gpid(pid),masstemp,pttemp,etatemp,ytemp,phitemp
270   Format(1x,I4,5E15.6)
C-OUT  Compute the cos(n*phi) Fourier component projections of the
C-OUT  azimuthal distributions for each PID type:
           total_mult_inc(pid) = total_mult_inc(pid) + 1
           jodd = 1
           do j = 1,nflowterms
              if(jodd.eq.1) then
                 if(ytemp.ge.0.0) then
                    cosinefac(j) =
     1                 cos(float(j)*(phitemp-PSIr_event)/rad)
                 else if(ytemp.lt.0.0) then
                    cosinefac(j) =
     1                -cos(float(j)*(phitemp-PSIr_event)/rad)
                 end if
              else if(jodd.eq.(-1)) then
                 cosinefac(j) =
     1              cos(float(j)*(phitemp-PSIr_event)/rad)
              end if
              jodd = -jodd
              Vn_sum(j,pid) = Vn_sum(j,pid) + cosinefac(j)
           end do
C-->           write(10,260) gpid(pid),(cosinefac(j),j=1,nflowterms)
260        Format(1x,I3,6E12.4)
C-OUT

           end do
           end do
240        Format('TRACK:',1x,I6,3(1x,G12.5),4(1x,I6))
CCC   End track kinematics output.

         else
C-->           write(7,250) ievent, event_abort, status_abort
C-->           write(8,250) ievent, event_abort, status_abort
250        Format(/'*** Event No. ',I7,' ABORTED: event_abort,',
     1     'status_abort = ',2I7)
         end if

      END DO !  End Event Loop

CCC   Output Anisotropic flow check sums to terminal:
      do pid = 1,n_pid_type 
C-->         write(6,*) 'Total incl # part type (',gpid(pid),
C-->     1              ') = ',total_mult_inc(pid)
         do j = 1,nflowterms
            Vn_sum(j,pid) = Vn_sum(j,pid)/total_mult_inc(pid)
C-->            write(6,*) 'Flow term #',j,': Check sum = ',
C-->     1                  Vn_sum(j,pid)
         end do
      end do

CCC   Output File Termination:
      ievent = -999
      total_mult = 0
      n_vertices = 0
C-->      write(7,241)
C-->      write(7,235) ievent,total_mult,n_vertices
241   Format(/'***  End of File  ***')

      Close(Unit=4)
      Close(Unit=7)
      Close(Unit=8)
C-OUT      Close(Unit=9)
C-OUT      Close(Unit=10)
      Close(Unit=9)
      Close(Unit=10)
      RETURN
      END

      Subroutine track_gen(gpid,N,model_type,T,sigma,expvel,
     1     reac_plane_cntrl,PSIr,Vn,
     2     pt_cut_min,pt_cut_max,eta_cut_min,eta_cut_max,
     3     y_cut_min,y_cut_max,
     4     phi_cut_min,phi_cut_max,n_scan_pts,irand,rad,pid,status,
     5     n_integ_pts)
      implicit none
      Include 'common_facfac.inc'
      Include 'Parameter_values.inc'
      Common/track/ pout(npid,4,factorial_max)
      real*4 pout
      Common/Geant/geant_mass(200),geant_charge(200),
     1             geant_lifetime(200),geant_width(200)
      real*4 geant_mass,geant_charge,geant_lifetime,geant_width

      real*4 T,sigma,expvel,pt_cut_min,pt_cut_max,eta_cut_min
      real*4 eta_cut_max,phi_cut_min,phi_cut_max,mass
      real*4 dpt,deta,facmax,pt,eta,y,rapidity,dNdpty,dNdp
      real*4 pt_trial,eta_trial,y_trial,ran,rad,phi
      real*4 y_cut_min,y_cut_max,pseudorapidity
      real*4 PSIr,Vn(nflowterms,4),dphi,facmax_phi,dNdphi
      real*4 Vnptmax,Vnymax,Vn_pt_y,Vntmp(nflowterms)
      real*4 masstmp,Mass_function

      integer gpid,N,model_type,npt,neta,n_scan_pts,ipt,ieta
      integer Track_count,irand,i,j,pid,status
      integer reac_plane_cntrl,iphi
      integer n_integ_pts

      do i = 1,factorial_max
      do j = 1,4
         pout(pid,j,i) = 0.0
      end do
      end do

      mass = geant_mass(gpid)
      npt  = n_scan_pts
      neta = n_scan_pts

CCC  Determine maximum value of model distribution in (pt,eta) range:
    
      dpt = (pt_cut_max - pt_cut_min)/float(npt - 1)
      deta = (eta_cut_max - eta_cut_min)/float(neta - 1)
      facmax = 0.0
      do ipt = 1,npt
         pt = pt_cut_min + dpt*float(ipt - 1)
         do ieta = 1,neta
            eta = eta_cut_min + deta*float(ieta - 1)
            y   = rapidity(mass,pt,eta)
            dNdp = dNdpty(1.0,pt,eta,y,mass,T,sigma,expvel,
     1                    model_type,1,pid)
            if(dNdp .gt. facmax) facmax = dNdp
         end do
      end do

CCC   If dNdp always underflows exp() range, then facmax will stay 
CCC   equal to 0.0, thus causing a divide by 0.0 error below. 
CCC   Check for this and Return if this is the case.  This event will
CCC   be aborted in this instance.

      if(facmax .eq. 0.0) then
         status = -86
         Return
      else
         status = 1
      end if

CCC   Determine the maximum values of the azimuthal model distribution
CCC   in phi, for case with reaction plane dependence and anisotropic flow
CCC   Find the absolute maximum possible value given the pt and y dependences
CCC   and assuming all Fourier components add with maximum coherence.

      if(reac_plane_cntrl .gt. 1) then
         facmax_phi = 1.0
         do i = 1,nflowterms
            if(i.eq.(2*(i/2))) then ! Select even Fourier components:
               Vnptmax = max(
     1                 abs(Vn(i,1)+Vn(i,4)*pt_cut_min
     3                    +Vn(i,2)*pt_cut_min*pt_cut_min),
     2                 abs(Vn(i,1)+Vn(i,4)*pt_cut_max
     4                    +Vn(i,2)*pt_cut_max*pt_cut_max))
               Vnymax  = max(
     1                 exp(-Vn(i,3)*y_cut_min*y_cut_min),
     2                 exp(-Vn(i,3)*y_cut_max*y_cut_max))
               if(y_cut_min.lt.0.0 .and. y_cut_max.gt.0.0) then
                  Vnymax  = max(Vnymax,1.0)
               end if
            else ! Select odd Fourier components
               Vnptmax = max(
     1                 abs(Vn(i,1)+Vn(i,2)*pt_cut_min),
     2                 abs(Vn(i,1)+Vn(i,2)*pt_cut_max))
               Vnymax  = max(
     1                 abs(Vn(i,3)+Vn(i,4)*abs(y_cut_min**3)),
     2                 abs(Vn(i,3)+Vn(i,4)*abs(y_cut_max**3)))
               if(y_cut_min.lt.0.0 .and. y_cut_max.gt.0.0) then
                  Vnymax  = max(Vnymax,abs(Vn(i,3)))
               end if
            end if
            facmax_phi = facmax_phi + 2.0*Vnptmax*Vnymax
         end do
CCC   Check facmax_phi; if 0, then event will be aborted.  
         if(facmax_phi.eq.0.0) then
            status = -86
            Return
         else
            status = 1
         end if
      end if

CCC Start Random Track Selection:

      Track_count = 0

100   pt_trial = ran(irand)*(pt_cut_max - pt_cut_min)+pt_cut_min
      if(model_type.ge.1 .and. model_type.le.4) then
         y_trial = ran(irand)*(y_cut_max - y_cut_min)+y_cut_min
         eta_trial = pseudorapidity(mass,pt_trial,y_trial)
         if(eta_trial.lt.eta_cut_min .or. eta_trial.gt.eta_cut_max)
     1      go to 100
      else if (model_type.eq.5 .or. model_type.eq.6) then
         eta_trial=ran(irand)*(eta_cut_max-eta_cut_min)+eta_cut_min
         y_trial = rapidity(mass,pt_trial,eta_trial)
      end if
      dNdp = dNdpty(1.0,pt_trial,eta_trial,y_trial,mass,T,sigma,
     1       expvel,model_type,1,pid)/facmax
      if(ran(irand) .le. dNdp) then
         Track_count = Track_count + 1

CCC   Determine phi angle:
         if(reac_plane_cntrl .eq. 1) then
            phi = (ran(irand)*(phi_cut_max - phi_cut_min) +
     1            phi_cut_min)/rad
         else if(reac_plane_cntrl .gt. 1) then
            do i = 1,nflowterms
               Vntmp(i) = Vn_pt_y(i,Vn(i,1),Vn(i,2),Vn(i,3),Vn(i,4),
     1                    pt_trial,y_trial)
            end do
101         phi = ran(irand)*(phi_cut_max - phi_cut_min)+phi_cut_min
            dNdphi = 1.0
            do i = 1,nflowterms
               dNdphi=dNdphi+2.0*Vntmp(i)*cos(float(i)*(phi-PSIr)/rad)
            end do
            if(ran(irand).gt.(dNdphi/facmax_phi)) then
               go to 101
            else
               phi = phi/rad
            end if
         end if

         masstmp = Mass_function(gpid,pid,n_integ_pts,irand)
         Call kinematics(masstmp,pt_trial,y_trial,phi,Track_count,pid)
         if(Track_count .lt. N) then
            go to 100
         else if(Track_count .ge. N) then
            Return
         end if

      else
         go to 100
      end if

      END

      Real*4 Function rapidity(m,pt,eta)
      implicit none
      real*4 m,pt,eta,theta,pz,E,pi,expR

      pi = 3.141592654
      theta = 2.0*ATAN(expR(-eta))
      if(abs(theta - pi/2.0) .lt. 0.000001) then
         pz = 0.0
      else
         pz = pt/tan(theta)
      end if
      E = sqrt(pt*pt + pz*pz + m*m)
      rapidity = 0.5*log((E+pz)/(E-pz))
      Return
      END

      Real*4 Function pseudorapidity(m,pt,y)
      implicit none

CCC   This Function computes the pseudorapidty for a given mass, pt, y:

      real*4 m,pt,y,mt,coshy,E,pzmag,pz,pmag,expR

      if(y.eq.0.0) then
         pseudorapidity = 0.0
         Return
      end if
      mt = sqrt(m*m + pt*pt)
      coshy = 0.5*(expR(y) + expR(-y))
      E = mt*coshy
      pzmag = sqrt(abs(E*E - mt*mt))
      if(y.eq.0.0) then
         pz = 0.0
      else
         pz = (y/abs(y))*pzmag
      end if
      pmag = sqrt(pt*pt + pz*pz)
      if(pt.ne.0.0) then
         pseudorapidity = 0.5*log((pmag+pz)/(pmag-pz))
      else if (pt.eq.0.0) then
         pseudorapidity = 86.0
C-->         write(8,10)
10       Format(10x,'Function pseudorapidity called with pt=0')
      end if
      Return
      END

      Subroutine kinematics(m,pt,y,phi,index,pid)
      implicit none

CCC  This SUBR takes the input particle mass (m), pt, y and phi and
CCC  computes E, px, py, pz and loads the momenta into the index-th
CCC  row of array pout(,,) in Common/track/ .

      integer index,pid
      Include 'common_facfac.inc'
      Include 'Parameter_values.inc'
      Common/track/ pout(npid,4,factorial_max)
      real*4 pout

      real*4 m,pt,y,phi,mt,coshy,E,pzmag,px,py,pz,expR

      mt = sqrt(m*m + pt*pt)
      coshy = 0.5*(expR(y) + expR(-y))
      E = mt*coshy
      pzmag = sqrt(abs(E*E - mt*mt))
      if(y.eq.0.0) then
         pz = 0.0
      else
         pz = (y/abs(y))*pzmag
      end if
      px = pt*cos(phi)
      py = pt*sin(phi)

      pout(pid,1,index) = px
      pout(pid,2,index) = py
      pout(pid,3,index) = pz
      pout(pid,4,index) = m

      Return
      END

      Subroutine kinematics2(m,px,py,pz,pt,eta,y,phi)
      implicit none

CCC  This SUBR takes the input particle mass (m), px,py,pz and
CCC  computes pt,eta,y,phi:

      real*4 m,px,py,pz,pt,eta,y,phi,pi,E,E0

      pi = 3.141592654
      pt = sqrt(px*px + py*py)
      E  = sqrt(m*m + pz*pz + pt*pt)
      y  = 0.5*log((E + pz)/(E - pz))
      E0 = sqrt(pz*pz + pt*pt)
      if(pt.eq.0.0) then
         eta = -86.0
      else
         eta = 0.5*log((E0 + pz)/(E0 - pz))
      end if
      if(px.eq.0.0 .and. py.eq.0.0) then
         phi = 0.0
      else
         phi = atan2(py,px)
      end if
      phi = (180.0/pi)*phi
      if(phi .lt. 0.0) phi = phi + 360.0
      Return
      END

      Real*4 Function dNdpty(A,pt,eta,y,m,T,sigma,vel,control,ktl,pid)
      implicit none

      real*4 A,pt,eta,y,m,T,sigma,vel
      real*4 pi,mt,coshy,E,ptot,FAC,gamma,yp,sinhyp,coshyp
      real*4 FAC2,FAC3,expR
      integer control, ktl,pid,pt_index,eta_index,index_locate
      Include 'Parameter_values.inc'
      Include 'common_dist_bin.inc'

CCC   Calculates dN/dp^3 using several models:
CCC
CCC      control = 1,  Humanic factorized model
CCC              = 2,  Pratt non-expanding spherical thermal source
CCC              = 3,  Bertsch non-expanding spherical thermal source
CCC              = 4,  Pratt spherical expanding thermally equilibrated
CCC                    source.
CCC              = 5,  Factorized pt, eta bin-by-bin distribution.
CCC              = 6,  Full 2D pt, eta bin-by-bin distribution.
CCC
CCC      ktl     = 0,  to return value of dN/dp^3
CCC      ktl     = 1,  to return value of dN/dpt*dy

      pi = 3.141592654
      mt = sqrt(pt*pt + m*m)
      coshy = 0.5*(expR(y) + expR(-y))
      E = mt*coshy
      ptot = sqrt(E*E - m*m)
      if(ktl .eq. 0) then
         FAC = 2.0*pi*pt*E
      else if(ktl .eq. 1) then
         FAC = 1.0
      end if

      if(control .eq. 1) then
         dNdpty = A*pt*expR(-mt/T)*expR(-y*y/(2.0*sigma*sigma))
         dNdpty = dNdpty/FAC

      else if(control .eq. 2) then
         dNdpty = A*pt*E*expR(-E/T)
         dNdpty = dNdpty/FAC

      else if(control .eq. 3) then
         dNdpty = A*pt*E/(expR(E/T) - 1.0)
         dNdpty = dNdpty/FAC

      else if(control .eq. 4) then
         gamma = 1.0/sqrt(1.0 - vel*vel)
         yp = gamma*vel*ptot/T
         sinhyp = 0.5*(expR(yp) - expR(-yp))
         coshyp = 0.5*(expR(yp) + expR(-yp))
         if(yp .ne. 0.0) then
            FAC2 = sinhyp/yp
         else
            FAC2 = 1.0
         end if
         FAC3 = FAC2+ (T/(gamma*E))*(FAC2 - coshyp)
         dNdpty = A*pt*E*expR(-gamma*E/T)*FAC3
         dNdpty = dNdpty/FAC

      else if(control .eq. 5) then
         pt_index = index_locate(pid,pt,1)
         eta_index = index_locate(pid,eta,2)
         dNdpty = A*pt_bin(pid,pt_index)*eta_bin(pid,eta_index)
         dNdpty = dNdpty/FAC

      else if(control .eq. 6) then
         pt_index = index_locate(pid,pt,1)
         eta_index = index_locate(pid,eta,2)
         dNdpty = A*pt_eta_bin(pid,pt_index,eta_index)
         dNdpty = dNdpty/FAC

      else
         dNdpty = -86.0

      end if

      return
      END

      Integer Function index_locate(pid,arg,kind)
      implicit none

      Include 'Parameter_values.inc'
      Include 'common_dist_bin.inc'

      integer pid,kind,ibin
      real*4 arg

CCC   This Function locates the pt or eta bin number corresponding to the
CCC   input bin mesh, the requested value of pt or eta, for the current
CCC   value of particle type.
CCC
CCC   If kind = 1, then pt bin number is located.
CCC   If kind = 2, then eta bin number is located.

      if(kind .eq. 1) then
         do ibin = 1,n_pt_bins(pid)
            if(arg.le.pt_bin_mesh(pid,ibin)) then
            index_locate = ibin
            Return
            end if
         end do
         index_locate = n_pt_bins(pid)
C-->         write(8,10) pid,arg
10       Format(//10x,'In Function index_locate, for pid = ',I5,
     1   'pt  =',E15.6,' is out of range - use last bin #')
         Return

      else if(kind .eq. 2) then

         do ibin = 1,n_eta_bins(pid)
            if(arg.le.eta_bin_mesh(pid,ibin)) then
            index_locate = ibin
            Return
            end if
         end do
         index_locate = n_eta_bins(pid)
C-->         write(8,11) pid,arg
11       Format(//10x,'In Function index_locate, for pid = ',I5,
     1   'eta =',E15.6,' is out of range - use last bin #')
         Return

      end if
      END

      Real*4 Function expR(x)
      implicit none
      real*4 x
      if(x .gt. 69.0) then
C-->         write(8,10) x
         STOP
      else if(x .lt. -69.0) then
         expR = 0.0
      else
         expR = exp(x)
      end if
10    Format(///10x,'Func. expR(x) called with x = ',E15.6,
     1    ', gt 69.0 - STOP')
      Return
      END

      SUBROUTINE LAGRNG (X,ARG,Y,VAL,NDIM,NFS,NPTS,MAXARG,MAXFS)
	IMPLICIT REAL*4(A-H,O-Z)
C
C     LAGRANGE INTERPOLATION,UNEQUALLY SPACED POINTS
C     ROUTINE OBTAINED FROM R. LANDAU, UNIV. OF OREGON.
C     ARG=VECTOR OF INDEPENDENT VARIABLE CONTAINING MAXARG VALUES.
C     VAL=MATRIX OF FUNCTION VALUES CORRESPONDING TO ARG. (MAXFS
C         FUNCTIONS AT MAXARG VALUES.)
C     X  =VALUE OF INDEP. VARIABLE FOR WHICH INTERPOLATION IS DESIRED.
C     Y  =VECTOR OF MAXFS FUNCTION VALUES RESULTING FROM SIMUL. INTERP.
C     NDIM=NUMBER OF ARG VALUES TO BE USED. (NDIM.LE.MAXARG)
C     NFS=NUMBER OF FUNCTIONS SIMUL. INTERP (NFS.LE.MAXFS)
C     NPTS=NUMBER OF POINTS USED IN INTERPOLATION. (NPTS=2,3,4,5,6)
C
      DIMENSION ARG(MAXARG), VAL(MAXFS,MAXARG), Y(MAXFS)
C
C     -----FIND X0, THE CLOSEST POINT TO X.
C
      NI=1
      NF=NDIM
   10 IF ((X.LE.ARG(NI)).OR.(X.GE.ARG(NF))) GO TO 30
      IF ((NF-NI+1).EQ.2) GO TO 70
      NMID=(NF+NI)/2
      IF (X.GT.ARG(NMID)) GO TO 20
      NF=NMID
      GO TO 10
   20 NI=NMID
      GO TO 10
C
C     ------ X IS ONE OF THE TABLULATED VALUES.
C
   30 IF (X.LE.ARG(NI)) GO TO 60
      NN=NF
   40 NUSED=0
      DO 50 N=1,NFS
   50 Y(N)=VAL(N,NN)
      RETURN
   60 NN=NI
      GO TO 40
C
C     ------- 2 PTS LEFT, CHOOSE SMALLER ONE.
C
   70 N0=NI
      NN=NPTS-2
      GO TO (110,100,90,80), NN
   80 CONTINUE
      IF (((N0+3).GT.NDIM).OR.((N0-2).LT.1)) GO TO 90
      NUSED=6
      GO TO 130
   90 CONTINUE
      IF ((N0+2).GT.NDIM) GO TO 110
      IF ((N0-2).LT.1) GO TO 100
      NUSED=5
      GO TO 130
  100 CONTINUE
      IF (((N0+2).GT.NDIM).OR.((N0-1).LT.1)) GO TO 110
      NUSED=4
      GO TO 130
  110 IF ((N0+1).LT.NDIM) GO TO 120
C
C     ------N0=NDIM, SPECIAL CASE.
C
      NN=NDIM
      GO TO 40
  120 NUSED=3
      IF ((N0-1).LT.1) NUSED=2
  130 CONTINUE
C
C     ------AT LEAST 2 PTS LEFT.
C
      Y0=X-ARG(N0)
      Y1=X-ARG(N0+1)
      Y01=Y1-Y0
      C0=Y1/Y01
      C1=-Y0/Y01
      IF (NUSED.EQ.2) GO TO 140
C
C     ------AT LEAST 3 PTS.
C
      YM1=X-ARG(N0-1)
      Y0M1=YM1-Y0
      YM11=Y1-YM1
      CM1=-Y0*Y1/Y0M1/YM11
      C0=C0*YM1/Y0M1
      C1=-C1*YM1/YM11
      IF (NUSED.EQ.3) GO TO 160
C
C     ------AT LEAST 4 PTS
C
      Y2=X-ARG(N0+2)
      YM12=Y2-YM1
      Y02=Y2-Y0
      Y12=Y2-Y1
      CM1=CM1*Y2/YM12
      C0=C0*Y2/Y02
      C1=C1*Y2/Y12
      C2=-YM1*Y0*Y1/YM12/Y02/Y12
      IF (NUSED.EQ.4) GO TO 180
C
C     ------AT LEAST 5 PTS.
C
      YM2=X-ARG(N0-2)
      YM2M1=YM1-YM2
      YM20=Y0-YM2
      YM21=Y1-YM2
      YM22=Y2-YM2
      CM2=YM1*Y0*Y1*Y2/YM2M1/YM20/YM21/YM22
      CM1=-CM1*YM2/YM2M1
      C0=-C0*YM2/YM20
      C1=-C1*YM2/YM21
      C2=-C2*YM2/YM22
      IF (NUSED.EQ.5) GO TO 200
C
C     ------AT LEAST 6 PTS.
C
      Y3=X-ARG(N0+3)
      YM23=Y3-YM2
      YM13=Y3-YM1
      Y03=Y3-Y0
      Y13=Y3-Y1
      Y23=Y3-Y2
      CM2=CM2*Y3/YM23
      CM1=CM1*Y3/YM13
      C0=C0*Y3/Y03
      C1=C1*Y3/Y13
      C2=C2*Y3/Y23
      C3=YM2*YM1*Y0*Y1*Y2/YM23/YM13/Y03/Y13/Y23
      GO TO 220
  140 CONTINUE
      DO 150 N=1,NFS
  150 Y(N)=C0*VAL(N,N0)+C1*VAL(N,N0+1)
      GO TO 240
  160 CONTINUE
      DO 170 N=1,NFS
  170 Y(N)=CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)
      GO TO 240
  180 CONTINUE
      DO 190 N=1,NFS
  190 Y(N)=CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)+C2*VAL(N,N0+2)
      GO TO 240
  200 CONTINUE
      DO 210 N=1,NFS
  210 Y(N)=CM2*VAL(N,N0-2)+CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)+C
     12*VAL(N,N0+2)
      GO TO 240
  220 CONTINUE
      DO 230 N=1,NFS
  230 Y(N)=CM2*VAL(N,N0-2)+CM1*VAL(N,N0-1)+C0*VAL(N,N0)+C1*VAL(N,N0+1)+C
     12*VAL(N,N0+2)+C3*VAL(N,N0+3)
  240 RETURN
C
      END

      Subroutine FACTORIAL
      implicit none

      integer n
      Include 'common_facfac.inc'
      real*4 FN

CCC   Computes the natural log of n! for n = 0,1,2...(factorial_max -1)
CCC   and puts the result in array FACLOG().
C
CCC   FACLOG(1) = log(0!) = 0.0
CCC   FACLOG(2) = log(1!) = 0.0
CCC   FACLOG(n+1) = log(n!)

      FACLOG(1) = 0.0
      FACLOG(2) = 0.0
      FN = 1.0
      do n = 3,factorial_max
      FN = FN + 1.0
      FACLOG(n) = FACLOG(n-1) + log(FN)
      end do
      Return
      END

      Subroutine MinMax(mean,stdev,n_stdev,min,max)
      implicit none

CCC   Computes range of integration for random number selections

      real*4 mean,stdev,n_stdev,min,max

      min = mean - n_stdev*stdev
      if(min .lt. 0.0) min = 0.0
      max = mean + n_stdev*stdev
      Return
      END

      Subroutine Poisson(min,max,mean,nsteps,integ,xfunc,ndim)
      implicit none

CCC   Computes Poisson distribution from n = min to max;
CCC   Integrates this distribution and records result at each step in
CCC      array integ();
CCC   Records the coordinates in array xfunc().

      integer min,max,mean,nsteps,ndim,i,n
      Include 'Parameter_values.inc'
      real*4 mean_real,mean_real_log,expR
      real*4 integ(ndim)
      real*4 xfunc(ndim)
      real*4 Poisson_dist(n_mult_max_steps)
      Include 'common_facfac.inc'

CCC Initialize arrays to zero:

      do i = 1,ndim
         integ(i) = 0.0
         xfunc(i) = 0.0
         Poisson_dist(i) = 0.0
      end do

      mean_real = float(mean)
      mean_real_log = log(mean_real) 

CCC   Compute Poisson distribution from n = min to max:
      do i = 1,nsteps
      n = (i - 1) + min
      Poisson_dist(i) = expR(-mean_real + float(n)*mean_real_log
     1      - FACLOG(n+1))
      end do

CCC   Integrate the Poisson distribution:
      integ(1) = 0.0
      do i = 2,nsteps
      integ(i) = 0.5*(Poisson_dist(i-1) + Poisson_dist(i)) + integ(i-1)
      end do

CCC   Normalize the integral to unity:
      do i = 1,nsteps
      integ(i) = integ(i)/integ(nsteps)
      end do

CCC   Fill xfunc array:
      do i = 1,nsteps
      xfunc(i) = float(i - 1 + min)
      end do

CCC  Extend integ() and xfunc() by one additional mesh point past the
CCC  end point in order to avoid a bug in the Lagrange interpolation
CCC  subroutine that gives erroneous interpolation results within the
CCC  the last mesh bin.

      integ(nsteps + 1) = integ(nsteps) + 0.01
      xfunc(nsteps + 1) = xfunc(nsteps)

      Return
      END

      Subroutine Gaussian(min,max,mean,stdev,npts,integ,xfunc,ndim)
      implicit none

CCC   Compute Gaussian distribution from x = min to max at npts;
CCC   Integrate this distribution and record result at each mesh in
CCC      array integ();
CCC   Record the coordinates in array xfunc().

      Include 'Parameter_values.inc'
      integer npts,ndim,i
      real*4 min,max,mean,stdev,integ(ndim),xfunc(ndim)
      real*4 dm,x,Gauss_dist(nmax_integ),FAC1,FAC2,pi,expR

CCC   Initialize arrays to zero:
      do i = 1,ndim
         integ(i) = 0.0
         xfunc(i) = 0.0
         Gauss_dist(i) = 0.0
      end do

      pi = 3.141592654
      FAC1 = 1.0/(sqrt(2.0*pi)*stdev)
      FAC2 = 2.0*stdev*stdev
      dm = (max - min)/float(npts - 1)

CCC   Compute normalized Gaussian distribution:
      do i = 1,npts
      x = min + dm*float(i-1)
      xfunc(i) = x
      Gauss_dist(i) = FAC1*expR(-((x - mean)**2)/FAC2)
      end do

CCC   Integrate Gaussian distribution over specified range
      integ(1) = 0.0
      do i = 2,npts
      integ(i) = 0.5*(Gauss_dist(i-1) + Gauss_dist(i))*dm + integ(i-1)
      end do

CCC   Normalize integral to unity:
      do i = 1,npts
      integ(i) = integ(i)/integ(npts)
      end do

CCC  Extend integ() and xfunc() by one mesh point to avoid Lagrange
CCC  interpolation subroutine bug:
      integ(npts + 1) = integ(npts) + 0.01
      xfunc(npts + 1) = xfunc(npts)

      Return
      END

      Subroutine Particle_prop
      implicit none

      Common/Geant/geant_mass(200),geant_charge(200),
     1             geant_lifetime(200),geant_width(200)

CCC   Local Variable Type Declarations:
      integer i
      real*4 geant_mass,geant_charge,geant_lifetime,geant_width
      real*4 hbar
      Parameter(hbar = 6.5821220E-25) ! hbar in GeV*seconds

CCC   Fill arrays geant_mass(),geant_charge(),geant_lifetime(),
CCC   geant_width() with particle mass in GeV, charge in units of
CCC   |e|, mean lifetime in seconds, and width in GeV, where
CCC   width*lifetime = hbar.  For lifetimes listed with values of
CCC   1.0E+30 (i.e. stable particles) the width is set to 0.000.
CCC   Row # = Geant Particle ID # code (see Geant Manual 3.10,
CCC   User's Guide, CONS 300-1 and 300-2).  Note, rows 151 and higher
CCC   are used for resonances.  These assume the masses and widths
CCC   specific to the models used to represent the invariant mass
CCC   distributions and therefore may differ slightly from the Particle
CCC   Data Group's listing.

CCC   Initialize arrays to zero:
      do i = 1,200
         geant_mass(i) = 0.0
         geant_charge(i) = 0.0
         geant_lifetime(i) = 0.0
         geant_width(i) = 0.0
      end do

CCC   Set Particle Masses in GeV:
      geant_mass(1)  = 0.0           ! Gamma
      geant_mass(2)  = 0.000511      ! Positron
      geant_mass(3)  = 0.000511      ! Electron
      geant_mass(4)  = 0.0           ! Neutrino
      geant_mass(5)  = 0.105659      ! Muon +
      geant_mass(6)  = 0.105659      ! Muon -
      geant_mass(7)  = 0.134693      ! Pion 0
      geant_mass(8)  = 0.139567      ! Pion +
      geant_mass(9)  = 0.139567      ! Pion -
      geant_mass(10) = 0.49767       ! Kaon 0 Long
      geant_mass(11) = 0.493667      ! Kaon +
      geant_mass(12) = 0.493667      ! Kaon -
      geant_mass(13) = 0.939573      ! Neutron
      geant_mass(14) = 0.93828       ! Proton
      geant_mass(15) = 0.93828       ! Antiproton
      geant_mass(16) = 0.49767       ! Kaon 0 Short
      geant_mass(17) = 0.5488        ! Eta
      geant_mass(18) = 1.11560       ! Lambda
      geant_mass(19) = 1.18936       ! Sigma +
      geant_mass(20) = 1.19246       ! Sigma 0
      geant_mass(21) = 1.19734       ! Sigma -
      geant_mass(22) = 1.31490       ! Xi 0
      geant_mass(23) = 1.32132       ! Xi -
      geant_mass(24) = 1.67245       ! Omega
      geant_mass(25) = 0.939573      ! Antineutron
      geant_mass(26) = 1.11560       ! Antilambda
      geant_mass(27) = 1.18936       ! Antisigma -
      geant_mass(28) = 1.19246       ! Antisigma 0
      geant_mass(29) = 1.19734       ! Antisigma +
      geant_mass(30) = 1.3149        ! Antixi 0
      geant_mass(31) = 1.32132       ! Antixi +
      geant_mass(32) = 1.67245       ! Antiomega +
      geant_mass(33) = 1.7842        ! Tau +
      geant_mass(34) = 1.7842        ! Tau -
      geant_mass(35) = 1.8694        ! D+
      geant_mass(36) = 1.8694        ! D-
      geant_mass(37) = 1.8647        ! D0
      geant_mass(38) = 1.8647        ! Anti D0
      geant_mass(39) = 1.9710        ! F+, now called Ds+
      geant_mass(40) = 1.9710        ! F-, now called Ds-
      geant_mass(41) = 2.2822        ! Lambda C+
      geant_mass(42) = 80.8000       ! W+
      geant_mass(43) = 80.8000       ! W-
      geant_mass(44) = 92.9000       ! Z0
      geant_mass(45) = 1.877         ! Deuteron
      geant_mass(46) = 2.817         ! Tritium
      geant_mass(47) = 3.755         ! Alpha
      geant_mass(48) = 0.0           ! Geantino
      geant_mass(49) = 2.80923       ! He3
      geant_mass(50) = 0.0           ! Cherenkov
      geant_mass(151) = 0.783        ! rho +
      geant_mass(152) = 0.783        ! rho -
      geant_mass(153) = 0.783        ! rho 0
      geant_mass(154) = 0.782        ! omega 0
      geant_mass(155) = 0.95750      ! eta'
      geant_mass(156) = 1.0194       ! phi
      geant_mass(157) = 3.09693      ! J/Psi
      geant_mass(158) = 1.232        ! Delta -
      geant_mass(159) = 1.232        ! Delta 0
      geant_mass(160) = 1.232        ! Delta +
      geant_mass(161) = 1.232        ! Delta ++
      geant_mass(162) = 0.89183      ! K* +
      geant_mass(163) = 0.89183      ! K* -
      geant_mass(164) = 0.89610      ! K* 0

CCC   Set Particle Charge in |e|:
      geant_charge( 1) =  0      ! Gamma
      geant_charge( 2) =  1      ! Positron
      geant_charge( 3) = -1      ! Electron
      geant_charge( 4) =  0      ! Neutrino
      geant_charge( 5) =  1      ! Muon+
      geant_charge( 6) = -1      ! Muon-
      geant_charge( 7) =  0      ! Pion0
      geant_charge( 8) =  1      ! Pion+
      geant_charge( 9) = -1      ! Pion-
      geant_charge(10) =  0      ! Kaon 0 long
      geant_charge(11) =  1      ! Kaon+
      geant_charge(12) = -1      ! Kaon-
      geant_charge(13) =  0      ! Neutron
      geant_charge(14) =  1      ! Proton
      geant_charge(15) = -1      ! Antiproton
      geant_charge(16) =  0      ! Kaon 0 short
      geant_charge(17) =  0      ! Eta
      geant_charge(18) =  0      ! Lambda
      geant_charge(19) =  1      ! Sigma+
      geant_charge(20) =  0      ! Sigma0
      geant_charge(21) = -1      ! Sigma-
      geant_charge(22) =  0      ! Xi 0
      geant_charge(23) = -1      ! Xi -
      geant_charge(24) = -1      ! Omega
      geant_charge(25) =  0      ! Antineutron
      geant_charge(26) =  0      ! Antilambda
      geant_charge(27) = -1      ! Anti-Sigma -
      geant_charge(28) =  0      ! Anti-Sigma 0
      geant_charge(29) =  1      ! Anti-Sigma +
      geant_charge(30) =  0      ! AntiXi 0
      geant_charge(31) =  1      ! AntiXi +
      geant_charge(32) =  1      ! Anti-Omega +
      geant_charge(33) =  1      ! Tau +
      geant_charge(34) = -1      ! Tau - 
      geant_charge(35) =  1      ! D+ 
      geant_charge(36) = -1      ! D- 
      geant_charge(37) =  0      ! D0 
      geant_charge(38) =  0      ! Anti D0 
      geant_charge(39) =  1      ! F+, now called Ds+ 
      geant_charge(40) = -1      ! F-, now called Ds- 
      geant_charge(41) =  1      ! Lambda C+ 
      geant_charge(42) =  1      ! W+ 
      geant_charge(43) = -1      ! W- 
      geant_charge(44) =  0      ! Z0
      geant_charge(45) =  1      ! Deuteron
      geant_charge(46) =  1      ! Triton
      geant_charge(47) =  2      ! Alpha
      geant_charge(48) =  0      ! Geantino (Fake particle)
      geant_charge(49) =  2      ! He3
      geant_charge(50) =  0      ! Cerenkov (Fake particle)
      geant_charge(151) =  1     ! rho+
      geant_charge(152) = -1     ! rho-
      geant_charge(153) =  0     ! rho 0
      geant_charge(154) =  0     ! omega 0
      geant_charge(155) =  0     ! eta'
      geant_charge(156) =  0     ! phi
      geant_charge(157) =  0     ! J/Psi
      geant_charge(158) = -1     ! Delta -
      geant_charge(159) =  0     ! Delta 0
      geant_charge(160) =  1     ! Delta +
      geant_charge(161) =  2     ! Delta ++
      geant_charge(162) =  1     ! K* +
      geant_charge(163) = -1     ! K* -
      geant_charge(164) =  0     ! K* 0

CCC   Set Particle Lifetimes in seconds:
      geant_lifetime( 1) = 1.0E+30       ! Gamma
      geant_lifetime( 2) = 1.0E+30       ! Positron
      geant_lifetime( 3) = 1.0E+30       ! Electron
      geant_lifetime( 4) = 1.0E+30       ! Neutrino
      geant_lifetime( 5) = 2.19703E-06   ! Muon+
      geant_lifetime( 6) = 2.19703E-06   ! Muon-
      geant_lifetime( 7) = 8.4E-17       ! Pion0
      geant_lifetime( 8) = 2.603E-08     ! Pion+
      geant_lifetime( 9) = 2.603E-08     ! Pion-
      geant_lifetime(10) = 5.16E-08      ! Kaon 0 long
      geant_lifetime(11) = 1.237E-08     ! Kaon+
      geant_lifetime(12) = 1.237E-08     ! Kaon-
      geant_lifetime(13) = 889.1         ! Neutron
      geant_lifetime(14) = 1.0E+30       ! Proton
      geant_lifetime(15) = 1.0E+30       ! Antiproton
      geant_lifetime(16) = 8.922E-11     ! Kaon 0 short
      geant_lifetime(17) = 5.53085E-19   ! Eta
      geant_lifetime(18) = 2.632E-10     ! Lambda
      geant_lifetime(19) = 7.99E-11      ! Sigma+
      geant_lifetime(20) = 7.40E-20      ! Sigma0
      geant_lifetime(21) = 1.479E-10     ! Sigma-
      geant_lifetime(22) = 2.90E-10      ! Xi 0
      geant_lifetime(23) = 1.639E-10     ! Xi -
      geant_lifetime(24) = 8.22E-11      ! Omega
      geant_lifetime(25) = 889.1         ! Antineutron
      geant_lifetime(26) = 2.632E-10     ! Antilambda
      geant_lifetime(27) = 7.99E-11      ! Anti-Sigma -
      geant_lifetime(28) = 7.40E-20      ! Anti-Sigma 0
      geant_lifetime(29) = 1.479E-10     ! Anti-Sigma +
      geant_lifetime(30) = 2.90E-10      ! AntiXi 0
      geant_lifetime(31) = 1.639E-10     ! AntiXi +
      geant_lifetime(32) = 8.22E-11      ! Anti-Omega +
      geant_lifetime(33) = 0.303E-12     ! Tau +
      geant_lifetime(34) = 0.303E-12     ! Tau -
      geant_lifetime(35) = 10.62E-13     ! D+
      geant_lifetime(36) = 10.62E-13     ! D-
      geant_lifetime(37) = 4.21E-13      ! D0
      geant_lifetime(38) = 4.21E-13      ! Anti D0
      geant_lifetime(39) = 4.45E-13      ! F+, now called Ds+
      geant_lifetime(40) = 4.45E-13      ! F-, now called Ds-
      geant_lifetime(41) = 1.91E-13      ! Lambda C+
      geant_lifetime(42) = 2.92E-25      ! W+
      geant_lifetime(43) = 2.92E-25      ! W-
      geant_lifetime(44) = 2.60E-25      ! Z0
      geant_lifetime(45) = 1.0E+30       ! Deuteron
      geant_lifetime(46) = 1.0E+30       ! Triton
      geant_lifetime(47) = 1.0E+30       ! Alpha
      geant_lifetime(48) = 1.0E+30       ! Geantino (Fake particle)
      geant_lifetime(49) = 1.0E+30       ! He3
      geant_lifetime(50) = 1.0E+30       ! Cerenkov (Fake particle)
      geant_lifetime(151) = 3.72E-24     ! rho +
      geant_lifetime(152) = 3.72E-24     ! rho -
      geant_lifetime(153) = 3.72E-24     ! rho 0
      geant_lifetime(154) = 7.84E-23     ! omega 0
      geant_lifetime(155) = 3.16E-21     ! eta'
      geant_lifetime(156) = 1.49E-22     ! phi
      geant_lifetime(157) = 9.68E-21     ! J/Psi
      geant_lifetime(158) = 9.27E-24     ! Delta -, Based on gamma = 71 MeV
      geant_lifetime(159) = 9.27E-24     ! Delta 0, -same-
      geant_lifetime(160) = 9.27E-24     ! Delta +, -same-
      geant_lifetime(161) = 9.27E-24     ! Delta ++,-same- 
      geant_lifetime(162) = 1.322E-23    ! K* +
      geant_lifetime(163) = 1.322E-23    ! K* -
      geant_lifetime(164) = 1.303E-23    ! K* 0

CCC   Set Particle Widths in GeV:
      do i = 1,200
         if(geant_lifetime(i) .le. 0.0) then
            geant_width(i) = 0.0
         else if(geant_lifetime(i) .ge. 1.0E+30) then
            geant_width(i) = 0.0
         else 
            geant_width(i) = hbar/geant_lifetime(i)
         end if
      end do

      Return
      END

      Real*4 Function Vn_pt_y(n,V1,V2,V3,V4,pt,y)
      implicit none

CCC   Description:  This function computes the pt,y dependent flow
CCC                 parameters Vn(pt,y) for the requested Fourier
CCC                 component, n = 1-6, pt (GeV/c) and y (rapidity).
CCC
CCC   Input:    n    = Fourier component, 1,2...6
CCC             V1   = Constant coefficient in pt dependent term
CCC             V2   = Coefficient of pt(pt^2) dependence for n odd (even).
CCC             V3   = Coefficient of y dependence; constant for n=odd,
CCC                    and inverse range squared for Gaussian for n=even.
CCC             V4   = Coefficient of y^3 dependence for n=odd; 
CCC                    pt dependence for n = even.
CCC             pt   = Transverse momentum (GeV/c)
CCC             y    = Rapidity relative to y(C.M.)
CCC
CCC   Output:   Vn_pt_y = Vn(pt,y) based on the model suggested by
CCC                       Art Poskanzer (LBNL, Feb. 2000)
CCC             Vn_pt_y = (V1 + V4*pt + V2*pt*pt)*exp(-V3*y*y) ; n=even
CCC             Vn_pt_y = (V1 + V2*pt)*sign(y)*(V3 + V4*abs(y**3)); n=odd

CCC   Local Variable Type Declarations:

      integer n
      real*4  V1,V2,V3,V4,pt,y,signy

      if(n .eq. (2*(n/2))) then
         Vn_pt_y = (V1 + V4*pt + V2*pt*pt)*exp(-V3*y*y)
      else
         if(y.ge.0.0) then
            signy = 1.0
         else if(y.lt.0.0) then
            signy = -1.0
         end if
         Vn_pt_y = (V1 + V2*pt)*signy*(V3 + V4*abs(y**3))
      end if
      Return
      END

      Subroutine Particle_mass(gpid,pid_index,npts)
      implicit none

CCC   Description:  This subroutine computes the mass distributions for
C     included resonances at 'npts' number of mesh points in mass from the
C     lower mass limit to an upper mass limit (listed below), integrates
C     this mass distribution, normalizes the integral to 1.0, and saves
C     the mass steps and integral function in the arrays in Common/Mass/.
C     The mass distribution integral is then randomly sampled in a later
C     step in order to get the specific resonance mass instances.
C     For non-resonant particles (i.e. either stable or those with negligible
C     widths) this subroutine returns without doing anything, leaving the
C     arrays in Common/Mass/ set to zero.  This subroutine is called for
C     a specific PID index, corresponding to the input list of particle
C     types.
C
C     Input:   gpid       = Geant particle ID code number, see SUBR:
C                           Particle_prop for listing.
C              pid_index  = Particle type array index, determined by input 
C                           particle list.
C              npts       = n_integ_pts in calling program; is the number
C                           of mass mesh points used to load the mass
C                           distribution integral.  Note that one extra
C                           mesh point is added to handle the bug in the
C                           Lagrange interpolator, LAGRNG.
C
C     Output:  Mass_integ_save( , ) - mass distribution integral
C              Mass_xfunc_save( , ) - mass distribution mesh
C              These are in Common/Mass/.

CCC   Include files and common blocks:
      Include 'Parameter_values.inc'
      Common/Geant/geant_mass(200),geant_charge(200),
     1             geant_lifetime(200),geant_width(200)
      real*4 geant_mass,geant_charge,geant_lifetime,geant_width
      Common/Mass/ Mass_integ_save(npid,nmax_integ),
     1             Mass_xfunc_save(npid,nmax_integ)
      real*4 Mass_integ_save,Mass_xfunc_save

CCC   Local Variable Type Declarations:
      integer gpid,pid_index,npts,i
      real*4 dist(nmax_integ),dM,M0,Gamma,GammaS
      real*4 M,Mpi,MK,MN,R_Delta,Jinc,qR
      real*4 Mrho_low,Mrho_high,Momega_low,Momega_high
      real*4 Mphi_low,Mphi_high,MDelta_low,MDelta_high
      real*4 MKstar_low,MKstar_high
      real*4 kcm,k0cm,Ecm,E0cm,beta,beta0,Epicm,ENcm,redtotE

CCC   Set Fixed parameter values:
      Parameter(Mpi = 0.1395675)     ! Charged pion mass (GeV)
      Parameter(MK  = 0.493646)      ! Charged kaon mass (GeV)
      Parameter(MN  = 0.938919)      ! Nucleon average mass (GeV)
      Parameter(R_Delta = 0.81)      ! Delta resonance range parameter 
      Parameter(Mrho_low = 0.28   )  ! Lower rho mass limit
      Parameter(Mrho_high = 1.200)   ! Upper rho mass limit (GeV)
      Parameter(Momega_low = 0.75)   ! Lower omega mass limit (GeV)
      Parameter(Momega_high = 0.81)  ! Upper omega mass limit (GeV)
      Parameter(Mphi_low = 1.009)    ! Lower phi mass limit (GeV)
      Parameter(Mphi_high = 1.029)   ! Upper phi mass limit (GeV)
      Parameter(MDelta_low = 1.1   ) ! Lower Delta mass limit 
      Parameter(MDelta_high = 1.400) ! Upper Delta mass limit (GeV)
      Parameter(MKstar_low = 0.74)   ! Lower Kstar mass limit (GeV)
      Parameter(MKstar_high = 1.04)  ! Upper Kstar mass limit (GeV)

CCC   Check npts:
      if(npts.le.1) Return

CCC   Load mass distribution for rho-meson:
      if(gpid.ge.151 .and. gpid.le.153) then
         do i = 1,nmax_integ
            dist(i) = 0.0
         end do
         M0 = geant_mass(gpid)
         Gamma = geant_width(gpid)
         dM = (Mrho_high - Mrho_low)/float(npts-1)
         do i = 1,npts
            M = Mrho_low + dM*float(i-1)
            Mass_xfunc_save(pid_index,i) = M
            kcm = sqrt(0.25*M*M - Mpi*Mpi)
            dist(i) = kcm/((M-M0)**2 + 0.25*Gamma*Gamma)
         end do

CCC   Load mass distribution for omega-meson:
      else if(gpid.eq.154) then
         do i = 1,nmax_integ
            dist(i) = 0.0
         end do
         M0 = geant_mass(gpid)
         Gamma = geant_width(gpid)
         dM = (Momega_high - Momega_low)/float(npts-1)
         do i = 1,npts
            M = Momega_low + dM*float(i-1)
            Mass_xfunc_save(pid_index,i) = M
            GammaS = Gamma*((M/M0)**3)
            dist(i) = M0*M0*Gamma/((M0*M0 - M*M)**2
     1              + M*M*GammaS*GammaS)
         end do

CCC   Load mass distribution for phi-meson:
      else if(gpid.eq.156) then
         do i = 1,nmax_integ
            dist(i) = 0.0
         end do
         M0 = geant_mass(gpid)
         Gamma = geant_width(gpid)
         dM = (Mphi_high - Mphi_low)/float(npts-1)
         k0cm = sqrt(0.25*M0*M0 - MK*MK)
         E0cm = sqrt(k0cm*k0cm  + MK*MK)
         beta0 = k0cm/E0cm
         do i = 1,npts
            M = Mphi_low + dM*float(i-1)
            Mass_xfunc_save(pid_index,i) = M
            kcm = sqrt(0.25*M*M - MK*MK)
            Ecm = sqrt(kcm*kcm  + MK*MK)
            beta = kcm/Ecm
            dist(i) = 0.25*Gamma*Gamma*((beta/beta0)**3)/
     1                ((M - M0)**2 + 0.25*Gamma*Gamma)
         end do

CCC   Load mass distribution for Delta resonances:
      else if(gpid.ge.158 .and. gpid.le.161) then
         do i = 1,nmax_integ
            dist(i) = 0.0
         end do
         M0 = geant_mass(gpid)
         Gamma = geant_width(gpid)
         dM = (MDelta_high - MDelta_low)/float(npts-1)
         do i = 1,npts
            M = MDelta_low + dM*float(i-1)
            Mass_xfunc_save(pid_index,i) = M
            kcm = ((M*M + Mpi*Mpi - MN*MN)**2)/(4.0*M*M) - Mpi*Mpi
            kcm = sqrt(abs(kcm))
            Epicm = sqrt(kcm*kcm + Mpi*Mpi)
            ENcm  = sqrt(kcm*kcm + MN*MN)
            redtotE = Epicm*ENcm/(Epicm + ENcm)
            Jinc = kcm/redtotE
            qR = kcm*R_Delta/Mpi
            GammaS = 2.0*qR*qR*qR*Gamma/(1.0 + qR*qR)
            dist(i) = (Jinc/(kcm*kcm))*GammaS*GammaS/
     1                ((M - M0)**2 + 0.25*GammaS*GammaS)
         end do

CCC   Load mass distribution for K*(892) resonances:
      else if(gpid.ge.162 .and. gpid.le.164) then
         do i = 1,nmax_integ
            dist(i) = 0.0
         end do
         M0 = geant_mass(gpid)
         Gamma = geant_width(gpid)
         dM = (MKstar_high - MKstar_low)/float(npts-1)
         k0cm = ((M0*M0 + Mpi*Mpi - MK*MK)**2)/(4.0*M0*M0) - Mpi*Mpi
         k0cm = sqrt(k0cm)
         do i = 1,npts
            M = MKstar_low + dM*float(i-1)
            Mass_xfunc_save(pid_index,i) = M
            kcm = ((M*M + Mpi*Mpi - MK*MK)**2)/(4.0*M*M) - Mpi*Mpi
            kcm = sqrt(kcm)
            qR = kcm/k0cm
            GammaS = 2.0*qR*qR*qR*Gamma/(1.0 + qR*qR)
            dist(i) = GammaS*GammaS*M0*M0/
     1                ((M*M - M0*M0)**2 + GammaS*GammaS*M0*M0)
         end do

CCC  Load additional resonance mass distributions here:

      else
         Return        ! Return for Geant PID types without mass distribution
      end if

CCC   Integrate mass distribution from mass_low to mass_high:

      Mass_integ_save(pid_index,1) = 0.0
      do i = 2,npts
         Mass_integ_save(pid_index,i) = 0.5*(dist(i-1) + dist(i))*dM
     1      + Mass_integ_save(pid_index,i-1)
      end do

CCC   Normalize this integral such that the last point is 1.00:
      if(Mass_integ_save(pid_index,npts) .ne. 0.0) then
         do i = 1,npts
         Mass_integ_save(pid_index,i) = 
     1   Mass_integ_save(pid_index,i)/Mass_integ_save(pid_index,npts)
         end do
      end if

CCC   Extend integ() and xfunc() by one mesh point to avoid Lagrange
CCC  interpolation subroutine bug:
      Mass_integ_save(pid_index,npts+1) =
     1   Mass_integ_save(pid_index,npts) + 0.01
      Mass_xfunc_save(pid_index,npts+1) = 
     1   Mass_xfunc_save(pid_index,npts)

      Return  
      END

      Real*4 Function Mass_function(gpid,pid_index,npts,irand)
      implicit none

CCC   Description:  For resonance particles which have mass distributions
C     this function randomly samples the distributions in Common/Mass/
C     and returns a particle mass in GeV in 'Mass_function'.
C     For non-resonant particles this function returns the Geant mass
C     listed in SUBR: Particle_prop.
C
C     Input:   gpid       = Geant particle ID code number, see SUBR:
C                           Particle_prop for listing.
C              pid_index  = Particle type array index, determined by input
C                           particle list.
C              npts       = n_integ_pts in calling program.  Is the number
C                           of mass mesh points for the arrays
C                           in Common/Mass/ minus 1.
C              irand      = random number generator argument/seed
C
C     Output:  Mass_function = particle or resonance mass (GeV)

CCC   Include files and common blocks:
      Include 'Parameter_values.inc'
      Common/Geant/geant_mass(200),geant_charge(200),
     1             geant_lifetime(200),geant_width(200)
      real*4 geant_mass,geant_charge,geant_lifetime,geant_width
      Common/Mass/ Mass_integ_save(npid,nmax_integ),
     1             Mass_xfunc_save(npid,nmax_integ)
      real*4 Mass_integ_save,Mass_xfunc_save

CCC   Local Variable Type Declarations:
      integer gpid,pid_index,npts,irand,i
      real*4 integ(nmax_integ),xfunc(nmax_integ)
      real*4 ran,masstmp

      if(Mass_integ_save(pid_index,npts) .ne. 0.0) then
         do i = 1,npts+1
            integ(i) = Mass_integ_save(pid_index,i)
            xfunc(i) = Mass_xfunc_save(pid_index,i)
         end do
         Call LAGRNG(ran(irand),integ,masstmp,xfunc,
     1               npts+1, 1, 5, npts+1, 1)
         Mass_function = masstmp
      else
         Mass_function = geant_mass(gpid)
      end if

      Return
      END
*
*      real*4 function ran(i)
*      implicit none
*      integer i
*      real*4 r
*      Call ranlux2(r,1,i)
*      ran = r
*      Return
*      END
*
*      Include 'ranlux2.F'






