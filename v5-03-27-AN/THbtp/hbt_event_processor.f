CCC ********************************************************************
CCC  Modifications for ALICE-ROOT application at CERN - 12/15/2000
CCC  1.  Removed all Fortran Data Structures in favor of labelled
CCC      common blocks.  The syntax of the structure_variable to
CCC      common variable name change is:
CCC      In STAR code    In ALICE code
CCC        A(i).B     ==>  A_B(i)
CCC        A(i).B(j)  ==>  A_B(j,i)
CCC  2.  All remaining references in the comments and write statements
CCC      to the data structures are interpreted as applying to the
CCC      new common variables.
CCC  3.  The UNIX random number generator, ran(), was replaced with
CCC      a function which calls a modified version of the CERNLIB
CCC      random number generator, ranlux, herein called ranlux2.  The
CCC      latter allows a user supplied seed value.
CCC  4.  Increased the following array sizes for the LHC Pb+Pb design
CCC      multiplicity criteria of dN_{ch}/dy = 8000 assuming 80% of
CCC      this is pi(+) and pi(-), which are the largest populations
CCC      that this code is required to process.
CCC
CCC      The following are for the max number of tracks that can be
CCC      stored in each sector without overflow:
CCC         common_mesh.inc - increased max_trk_save from 30 to 150
CCC         common_sec_track.inc  - inc. max_trk_sec from 30 to 150
CCC         common_sec_track2.inc - inc. max_trk_sec2 from 30 to 150
CCC
CCC      The following determine the maximum number of tracks that can
CCC      be processed in an event:
CCC         common_track.inc - increased trk_maxlen from 6500 to 25000
CCC         common_track2.inc - increased trk2_maxlen from 6500 to 25000
CCC
C
C     DESCRIPTION OF METHOD:
C     =====================
C
C          This program produces relativistic heavy-ion collision
C     events in which the particle momenta for selected particle ID
C     types and for selected kinematic acceptance ranges are randomly
C     adjusted such that specified one-body distributions and two-body
C     correlation functions are reproduced.  The input to the code may
C     be a set of events from any STAR event generator, so long as the
C     format is in the STAR Geant (GSTAR) text format standard (see
C     STAR Note #235).  The basic method is similar to that of Metropolis
C     et al. and is fully described in Ray and Hoffmann, Phys. Rev. C 54,
C     2582 (1996).  Briefly the steps in the algorithm include:
C
C        (1) For an initial particle distribution of specified particle
C            ID types (maximum of two types allowed) and momentum
C            acceptance ranges [given in terms of transverse momentum 
C            (p_T), azimuthal angle (phi) and pseudorapidity (eta)] the
C            momentum vector of one particle is randomly shifted within
C            a specified range from its initial value.  The shifts are
C            done for px, py and pz independently.
C
C        (2) New one-body and two-body histograms, as well as the two-body
C            correlation function are calculated.
C
C        (3) If the random momentum shift results in an improved overall
C            chi-square, obtained by comparison with a specified reference 
C            for the one-body distribution and the two-body correlation
C            model, then the new momentum vector is retained.  If not,
C            then the vector is restored to its starting value.
C
C        (4) Steps 1-3 are repeated for each accepted particle in the
C            event.
C
C        (5) The entire process, steps 1-4, is repeated until either a
C            satisfactory fit to the model distributions are obtained or 
C            the maximum number of iterations is reached.
C
C        (6) Once the iterative process is complete, the input event file
C            is copied directly to an output event file where the adjusted
C            momentum values for the accepted tracks replace that in the
C            input file.  The event output file is in the GSTAR standard
C            text format.  This event output file may be processed again
C            by this code in order to generate correlations for other 
C            particle types or for different kinematic ranges.  The file
C            is suitable for input into the STAR version of Geant-3, called
C            GSTAR (STAR Note 235).
C
C     In order to reduce cpu demand the particle momenta are sorted into
C     momentum space cells or sectors.  In forming particle pairs only those
C     particles in the same, or adjacent cells are used.  For large events
C     this vastly reduces the required cpu time.  The penalty is that the
C     coding becomes more complicated.  Much of the present code is devoted
C     to the necessary bookeeping chores that must be done in order to
C     determine which cell the tracks are in and if they move to new cells
C     during the track adjustment procedure.  Information about the 
C     momentum space cells are contained in the data structure /sec_trk_map/.
C
C     The sector size must therefore be scaled with the specified correlation
C     range.  All particles will be paired with all possible partners out
C     to Q's equal to the smallest dimension of the momentum space sectors.
C     Particle pairs with Q greater than this sector dimension will suffer
C     reduced acceptance, finally being completely cut-off for Q ~> 2 times
C     the diagonal length thru a sector. 
C
C          In order to generate momentum correlations for particle types
C     having low multiplicity it is necessary for the user to supply this
C     code with an artificially enhanced multiplicity along with a track-
C     write-output fractional acceptance factor (see input variable
C     'trk_accep').  For example, if the user wants to generate HBT
C     correlations for K0-shorts but the assumed multiplicity is too
C     low for the present algorithm to work, the user may increase the
C     input K0-short multiplicity, for instance by a factor of 5, then
C     run the code and set trk_accep = 1/5 = 0.2 in order to randomly
C     reject about 80% of the K0-shorts.  The track rejection is done
C     after the track adjustment procedure is completed.  This procedure
C     should preserve the built-in correlations in the final output
C     event file, although with reduced statistics.
C          Another approach for handling low multiplicity events would
C     be to combine particles from several events, carry out the track
C     adjustment procedure, then separate the tracks into their original
C     events.  This method must insure that no bias is included due to
C     the order of processing the tracks from the first, second, etc.
C     events.  This latter method, once proven, could be used for
C     the low multiplicity particles from any event generator.  For
C     the present version of the code the low multiplicity HBT studies
C     must utilize a Monte Carlo multiplicity generator.
C
C          The code may also be used to read in a previously correlated 
C     set of events and either compute the reference histograms or read in
C     the references, and then compute the correlations for each event and
C     the inclusive correlations without doing the track momentum adjustment
C     procedure.  This feature may be used, for example, to study the 
C     correlations that result in one set of coordinates for events fitted
C     to correlations with respect to a different set of coordinates.  For
C     example, fit the correlations to the Y-K-P form and then evaluate
C     the side-out-long correlations, or vice-versa.
C 
C     TWO-BODY REFERENCE HISTOGRAMS:
C     =============================
C
C          In order to calculate the correlations, an uncorrelated two-body
C     reference spectrum is needed.  The program will calculate this
C     quantity by forming pairs of particles from different events in the
C     input file.  For the particle ID type(s) and momentum acceptance
C     the code forms all possible pairs (given the cell substructure) by 
C     mixing particles from event#1 with those in event#2, then particles
C     from event#2 are mixed with particles from event#3, then events 3
C     and 4 are mixed, and so on.  In this way ample statistics may be
C     achieved for the reference distributions.  These reference distributions
C     can be written out to file and re-used in subsequent runs.  Since
C     all events in the input event file are used in generating the 
C     reference distribution, it is imperative that these events be physically
C     similar.
C
C     ONE-BODY REFERENCE HISTOGRAMS:
C     =============================
C
C          Inclusive sums of the accepted particles for all events in the
C     input event file determine the one-body reference distributions.
C     These are used to constrain the momentum vector shifts.  Although
C     the one-body distributions from realistic event generators are fully
C     three-dimensional, the present code is restricted to only work with
C     the one-dimensional projections onto p_T, phi and eta.  In other words,
C     the p_T distribution used in this code is formed by integrating 
C     the particle distributions over (phi,eta) over the momentum acceptance
C     range.  No particle distribution models are built into the code.
C     The one-body reference distributions are either read-in or determined
C     from the events in the input event text file.
C
C     TWO-BODY CORRELATION MODELS:
C     ===========================
C
C          The code permits both 1-dimensional and 3-dimensional two-body
C     correlation models.  These may be fitted separately or simultaneously.
C     The source may include a mixture of incoherent and coherent components;
C     Coulomb corrections are also included.  The general form assumed 
C     [see Cramer and Kadija, Phys. Rev. C 53, 908 (1996)] is:
C
C       C2 = 1 + lambda*(b**2) + 2.0*sqrt(lambda)*[1 - sqrt(lambda)]*b 
C
C     where lambda is the usual chaoticity parameter.  The third term in
C     this equation may be turned on or off.  Values of lambda < 1.0 may
C     be used without the third term being included.  For 1-dimensional
C     functions b is given by:
C
C       b = exp(- Q**2 * R**2 / 2)
C
C     where Q is either the invariant 4-momentum difference, the total
C     4-momentum difference (i.e. time-like + space-like) or the
C     3-vector momentum difference.  The 3-dimensional functions may be
C     of the Bertsch-Pratt ``side-out-longitudinal'' form given by:
C
C       b = exp[(- Qside**2 * Rside**2 - Qout**2 * Rout**2
C                - Qlong**2 * Rlong**2)/2]
C
C     where the ``out-long'' cross term is omitted.  The 3D function may 
C     also be in the Yano-Koonin-Podgoretski (YKP) form given by (for
C     pairs in the A+A c.m. frame):
C
C       b = exp[(- Qperp**2 * Rperp**2 - Qparallel**2 * Rparallel**2
C                - Qtime**2 * Rtime**2)/2]
C
C     where
C                Qperp = transverse momentum difference
C                Qparallel = Qlong = p_{1z} - p_{2z}
C                Qtime = E_1 - E_2
C
C     The Coulomb correction may be omitted, or included in one of 3 ways:
C
C          (1) Point source Gamow factor
C          (2) Finite source NA35 model (see Eq.(5) in Z. Phys. C73, 443
C              (1997)) where
C              
C              Coulomb correction = 1 + [G(q) - 1]*exp(-q/Q0) 
C
C              and G(q) is the Gamow factor and q is the 3-momentum
C              vector difference.
C          (3) Finite source, Pratt integrated Coulomb wave function
C              integrals, interpolated for a specific source radius
C              from John Cramer's tables for discrete source radii.
C
C     These Coulomb correction factors multiply the above correlation
C     function to give the total correlation function to be fitted for
C     charged, like pairs. For charged, unlike pairs only the Coulomb
C     (attractive) correlation function is used.
C
C     BINNING:
C     =======
C
C          Several types of binning are done in the code. The one-body
C     distributions are binned in p_t, phi and eta.  The full momentum
C     space is subdivided into cells or sectors.  The 1D and 3D two-body
C     distributions are binned into fine and coarse bins, where the fine
C     bins occur at the smaller momentum difference range and the coarse
C     bins at the larger.  For the 3D distributions the (1,1,1) coarse
C     bin must coincide with the 3D fine bins.
C
C     SUMMARY OF EXTERNAL FILES:
C     =========================
C
C File Unit#    File Name                Description
C ----------------------------------------------------------------------------
C     1      hbt_parameters.in          Input switches and controls
C     2      event_text.in              Event generator output in GSTAR text
C     3      event_line.flags           Generated tmp file, input line flags
C     4      event_tracks.select        Generated tmp file, accep. tracks flg.
C     7      hbt_log.out                Generated log file - error reports
C     8      hbt_simulation.out         Generated main output file
C     9      hbt_pair_reference.hist    Generated pair ref. histograms
C     10     event_hbt_text.out         Gen. correlated event text file
C     11     hbt_singles_reference.hist Gen. one-body ref. histograms
C     12     event_text_aux.in          Tmp copy of event_text.in per event
C     14     event_tracks_aux.select    Tmp copy of event_tracks.select/event
C     21-27  cpp_*.dat (*=06,08...18)   Like pair Pratt Coulomb corrections.
C     31-37  cpm_*.dat (*=06,08...18)   Unlike pair Pratt Coulomb corrects.
C ----------------------------------------------------------------------------
C
C     Source of Data for ALICE Application:
C     ====================================
C
C File Unit#    File Name                For ALICE File or Struc?
C ----------------------------------------------------------------------------
C     1      hbt_parameters.in          Call AliHbtp_ function
C     2      event_text.in              Call AliHbtp_ function
C     3      event_line.flags           File not used
C     4      event_tracks.select        File not used
C     7      hbt_log.out                File used as is
C     8      hbt_simulation.out         File used as is
C     9      hbt_pair_reference.hist    File used as is
C     10     event_hbt_text.out         Call AliHbtp_ function
C     11     hbt_singles_reference.hist File used as is
C     12     event_text_aux.in          File not used
C     14     event_tracks_aux.select    File not used
C     21-27  cpp_*.dat (*=06,08...18)   Files are used as is
C     31-37  cpm_*.dat (*=06,08...18)   Files are used as is
C ----------------------------------------------------------------------------
C
C     DESCRIPTION OF INPUT PARAMETERS AND SWITCHES (FILE: hbt_parameters.in):
C     ======================================================================
C
C     Control Switches:
C     ================
C
C        ref_control       = 1 to read reference histograms from input files
C                          = 2 to compute reference histograms by track 
C                              mixing from event pairs in the event input file.
C
C        switch_1d         = 0 to not compute the 1D two-body correlations.
C                          = 1 to compute this using Q-invariant
C                          = 2 to compute this using Q-total
C                          = 3 to compute this using Q-3-vector
C
C        switch_3d         = 0 to not compute the 3D two-body correlations.
C                          = 1 to compute this using the side-out-long form
C                          = 2 to compute this using the YKP form.
C
C        switch_type       = 1 to fit only the like pair correlations
C                          = 2 to fit only the unlike pair correlations
C                          = 3 to fit both the like and unlike pair correl.
C
C        switch_coherence  = 0 for purely incoherent source (but can have
C                              lambda < 1.0)
C                          = 1 for mixed incoherent and coherent source
C
C        switch_coulomb    = 0 no Coulomb correction
C                          = 1 Point source Gamow correction only
C                          = 2 NA35 finite source size correction
C                          = 3 Pratt/Cramer finite source size correction;
C                              interpolated from input tables.
C
C        switch_fermi_bose =  1 Boson pairs
C                          = -1 Fermion pairs
C
C        trk_accep         = 1.0 all adjusted tracks are written out
C                          < 1.0 only this fraction, on average, of the
C                            adjusted tracks are written out.  Used for
C                            low multiplicity events.
C
C        print_full        = 0 for standard, minimal output
C                          = 1 for full, comprehensive (large) output for
C                              each event.
C
C        print_sector_data = 0 std. sector occupancy data printed out
C                          = 1 to print sector occupancy and overflow info.
C
C     Particle ID and Search Parameters:
C     =================================
C
C        n_pid_types       = 1 or 2 only, # particle types to correlate
C 
C        pid(1), pid(2)    = Geant-3 particle ID code numbers
C
C        deltap            = maximum range for random momentum shifts in
C                            GeV/c; px,py,pz independent; Def = 0.1 GeV/c.
C
C        maxit             = maximum # allowed iterations thru all the
C                            tracks for each event. Default = 50.
C                            If maxit=0, then calculate the correlations 
C                            for the input set of events without doing the
C                            track adjustment procedure.
C
C        delchi            = min % change in total chi-square for which
C                            the track adjustment iterations may stop,
C                            Default = 0.1%.
C
C        irand             = initial random # seed, default = 12345
C
C     Source Function Parameters:
C     ==========================
C
C        lambda            = Chaoticity
C
C        R_1d              = Spherical source model radius (fm)
C
C        Rside,Rout,Rlong  = Non-spherical Bertsch-Pratt source model (fm)
C
C        Rperp,Rparallel,R0= Non-spherical Yano-Koonin-Podgoretski source
C                            model (fm).
C
C        Q0                = NA35 Coulomb parameter for finite source size
C                            in (GeV/c) - iff switch_coulomb = 2
C                          = Spherical Coulomb source radius in (fm) iff
C                            switch_coulomb = 3, used to interpolate the
C                            input Pratt/Cramer discrete Coulomb source
C                            radii tables.
C
C     One-body pT, phi, eta Acceptance Bins:
C     =====================================
C
C        n_pt_bins, pt_min, pt_max  = # pt bins,  min/max pt accept. (GeV/c)
C
C        n_phi_bins,phi_min,phi_max = # phi bins, min/max phi accept. (deg.)
C
C        n_eta_bins,eta_min,eta_max = # eta bins, min/max eta accept. 
C
C                                     [NOTE: For each the maximum # of bins
C                                            must be .le. 100]
C
C     Two-body 1D and 3D Correlation Bins:
C     ===================================
C
C        n_1d_fine,  binsize_1d_fine   = # and size (GeV/c), 1D - fine mesh
C
C        n_1d_coarse,binsize_1d_coarse = # and size (GeV/c), 1D - coarse mesh
C
C        n_3d_fine,  binsize_3d_fine   = # and size (GeV/c), 3D - fine mesh
C
C        n_3d_coarse,binsize_3d_coarse = # and size (GeV/c), 3D - coarse mesh
C
C                                      [NOTE: The maximum # of 1D bins (fine
C                                             + coarse) must be .le. 100;
C                                             The maximum # of 3D bins (either
C                                             fine or coarse) must be .le.10).
C                                             For both 1D and 3D there must be
C                                             at least 1 fine bin and 1 coarse
C                                             bin.]
C        n_3d_fine_project             = # of 3D-fine bins to integrate over
C                                        to form 1D projections.  This value
C                                        must be .le. n_3d_fine.
C
C     Momentum Space Track-Sector Cells:
C     =================================
C
C        n_px_bins,px_min,px_max = #, min,max px bins (GeV/c)
C
C        n_py_bins,py_min,py_max = #, min,max py bins (GeV/c)
C
C        n_pz_bins,pz_min,pz_max = #, min,max pz bins (GeV/c)
C
C                                  [NOTE: The maximum number of total sectors,
C                                         equal to the product of the x-y-z
C                                         number of cells must be .le.
C                                         sec_maxlen which is defined in the
C                                         /sec_trk_map/ data structure.]
C
C     Relative Chi-Square Weights:
C     ===========================
C
C        chisq_wt_like_1d          = 1D, like pairs
C        chisq_wt_unlike_1d        = 1D, unlike pairs
C        chisq_wt_like_3d_fine     = 3D, like pairs, fine mesh
C        chisq_wt_unlike_3d_fine   = 3D, unlike pairs, fine mesh
C        chisq_wt_like_3d_coarse   = 3D, like pairs, coarse mesh
C        chisq_wt_unlike_3d_coarse = 3D, unlike pairs, coarse mesh
C        chisq_wt_hist1_1 = summed pt, phi, eta 1-body dist., PID#1
C        chisq_wt_hist1_2 = summed pt, phi, eta 1-body dist., PID#2
C
C
C     FORMAT for hbt_singles_reference.hist:
C     =====================================
C
C     The output content for the one-body reference histograms is:
C
C     Line 1:  n_pid_types,pid(1),pid(2)
C          2:  n_pt_bins,pt_min,pt_max
C          3:  n_phi_bins,phi_min,phi_max
C          4:  n_eta_bins,eta_min,eta_max
C          5:  n_part_used_1_ref,n_part_used_2_ref
C
C     Then for PID #1:      (href1_pt_1(i),i=1,n_pt_bins)
C     (One entry per line)  (href1_phi_1(i),i=1,n_phi_bins)
C                           (href1_eta_1(i),i=1,n_eta_bins)
C
C     and for PID #2:       (href1_pt_2(i),i=1,n_pt_bins)
C     (One entry per line)  (href1_phi_2(i),i=1,n_phi_bins)
C                           (href1_eta_2(i),i=1,n_eta_bins)
C
C
C     FORMAT for hbt_pair_reference.hist:
C     ==================================
C
C     The output content for the two-body reference histograms is:
C
C     Line 1:  n_pid_types,pid(1),pid(2)
C          2:  n_pt_bins,pt_min,pt_max
C          3:  n_phi_bins,phi_min,phi_max
C          4:  n_eta_bins,eta_min,eta_max
C          5:  switch_1d,switch_3d,switch_type
C          6:  n_1d_fine,n_1d_coarse,n_3d_fine,n_3d_coarse
C          7:  binsize_1d_fine,binsize_1d_coarse,
C              binsize_3d_fine,binsize_3d_coarse
C          8:  num_pairs_like_ref,num_pairs_unlike_ref
C
C   The pair distributions (with one entry per line) are:
C
C   1D, like pairs:             (href_like_1d(i),i=1,n_1d_total)
C
C   1D, unlike pairs:           (href_unlike_1d(i),i=1,n_1d_total)
C
C   3D, like pairs, fine mesh:   href_like_3d_fine(i,j,k) ; (i,(j,(k,...)))
C
C   3D, like pairs, coarse mesh: href_like_3d_coarse(i,j,k) ; (i,(j,(k,...)))
C
C   3D, unlike, fine mesh:       href_unlike_3d_fine(i,j,k) ; (i,(j,(k,...)))
C
C   3D, unlike, coarse mesh:     href_unlike_3d_coarse(i,j,k) ; (i,(j,(k,...)))
C
C*************************************************************************
C*************************************************************************

      SUBROUTINE CTEST
      implicit none

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'
      Include 'common_coulomb.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'
      Include 'common_particle.inc'
      
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' '
      
      write(*,*) 'Input data in Fort'
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '      PARAMETERS'
      write(*,*) ' '
      write(*,*) 'ref_control',ref_control
      write(*,*) 'switch_1d',switch_1d 
      write(*,*) 'switch_3d',switch_3d
      write(*,*) 'switch_type',switch_type
      write(*,*) 'switch_coherence',switch_coherence
      write(*,*) 'switch_coulomb',switch_coulomb
      write(*,*) 'switch_fermi_bose',switch_fermi_bose
      write(*,*) 'trk_accep',trk_accep
      write(*,*) 'print_full',print_full
      write(*,*) 'print_sector_data',print_sector_data
      write(*,*) 'n_pid_types',n_pid_types
      write(*,*) 'pid(1)', pid(1)
      write(*,*) 'pid(2)', pid(2)
      write(*,*) 'maxit',maxit
      write(*,*) 'irand',irand
      write(*,*) 'n_part_1_trk', n_part_1_trk      
      write(*,*) 'n_part_2_trk ', n_part_2_trk      
      write(*,*) 'n_part_tot_trk ', n_part_tot_trk      
      write(*,*) 'n_part_used_1_trk ', n_part_used_1_trk      
      write(*,*) 'n_part_used_2_trk',  n_part_used_2_trk     
      write(*,*) 'n_part_1_trk2',    n_part_1_trk2   
      write(*,*) 'n_part_2_trk2', n_part_2_trk2      
      write(*,*) 'n_part_tot_trk2',   n_part_tot_trk2    
      write(*,*) 'n_part_used_1_trk2',  n_part_used_1_trk2     
      write(*,*) 'n_part_used_2_trk2',  n_part_used_2_trk2     
      write(*,*) 'n_part_used_1_ref', n_part_used_1_ref      
      write(*,*) 'n_part_used_2_ref ',  n_part_used_2_ref     
      write(*,*) 'n_part_used_1_inc',  n_part_used_1_inc     
      write(*,*) 'n_part_used_2_inc',     n_part_used_2_inc  
      write(*,*) 'num_pairs_like',   num_pairs_like    
      write(*,*) 'num_pairs_unlike',  num_pairs_unlike     
      write(*,*) 'num_pairs_like_ref ', num_pairs_like_ref      
      write(*,*) 'num_pairs_like_inc ',  num_pairs_like_inc     
      write(*,*) 'num_pairs_unlike_inc', num_pairs_unlike_inc      
      write(*,*) 'event_line_counter', event_line_counter      
      write(*,*) 'file10_line_counter ',file10_line_counter       
      write(*,*) 'lambda',lambda
      write(*,*) 'R_1d ',R_1d 
      write(*,*) 'Rside',Rside
      write(*,*) 'Rout  ',  Rout
      write(*,*) 'Rlong  ', Rlong 
      write(*,*) 'Rperp  ',  Rperp
      write(*,*) 'Rparallel  ', Rparallel 
      write(*,*) 'R0 ',  R0
      write(*,*) 'Q0 ',  Q0
      write(*,*) 'deltap',deltap
      write(*,*) 'delchi',delchi
      write(*,*) 'pi ', pi 
      write(*,*) 'rad ', rad 
      write(*,*) 'hbc ', hbc 
      write(*,*) 'chisq_wt_like_1d ', chisq_wt_like_1d 
      write(*,*) 'chisq_wt_unlike_1d ',  chisq_wt_unlike_1d
      write(*,*) 'chisq_wt_like_3d_fine  ',chisq_wt_like_3d_fine  
      write(*,*) 'chisq_wt_unlike_3d_fine ',  chisq_wt_unlike_3d_fine
      write(*,*) 'chisq_wt_like_3d_coarse ', chisq_wt_like_3d_coarse
      write(*,*) 'chisq_wt_unlike_3d_coarse',chisq_wt_unlike_3d_coarse 
      write(*,*) 'chisq_wt_hist1_1 ', chisq_wt_hist1_1 
      write(*,*) 'chisq_wt_hist1_2 ', chisq_wt_hist1_2 
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '         MESH       '  
      write(*,*) ' '
      write(*,*) ' n_pt_bins ',  n_pt_bins
      write(*,*) ' pt_min ', pt_min 
      write(*,*) ' pt_max  ', pt_max  
      write(*,*) ' n_phi_bins ', n_phi_bins 
      write(*,*) ' phi_min ', phi_min 
      write(*,*) ' phi_max ', phi_max 
      write(*,*) ' n_eta_bins ', n_eta_bins 
      write(*,*) ' eta_min', eta_min 
      write(*,*) ' eta_max', eta_max 
      write(*,*) ' n_1d_fine ', n_1d_fine 
      write(*,*) ' binsize_1d_fine ',binsize_1d_fine  
      write(*,*) ' n_1d_coarse ',n_1d_coarse  
      write(*,*) ' binsize_1d_coarse ', binsize_1d_coarse 
      write(*,*) ' n_3d_fine  ',n_3d_fine   
      write(*,*) ' binsize_3d_fine ',binsize_3d_fine  
      write(*,*) ' n_3d_coarse ', n_3d_coarse 
      write(*,*) ' binsize_3d_coarse ', binsize_3d_coarse 
      write(*,*) ' n_3d_fine_project ',n_3d_fine_project  
      write(*,*) ' n_px_bins ',n_px_bins  
      write(*,*) ' px_min ',px_min  
      write(*,*) ' px_max ', px_max 
      write(*,*) ' n_py_bins ', n_py_bins 
      write(*,*) ' py_min ',  py_min
      write(*,*) ' py_max',  py_max
      write(*,*) ' n_pz_bins ',n_pz_bins  
      write(*,*) ' pz_min ', pz_min 
      write(*,*) ' pz_max ', pz_max 
      write(*,*) '  '
   
      End



      SUBROUTINE HBTPROCESSOR
      implicit none


      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'
      Include 'common_coulomb.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'
      Include 'common_particle.inc'

CCC   Set Data I/O control for ALICE or Standalone application
C     ALICE = 1   ! This is for the ALICE AliRoot application
CCC   ALICE = 0   ! This is for the standalone application

CCC   Initialize error code for ALICE application:
      errorcode = 0

CCC   Open Output Files:

      open(unit=7,status='unknown',access='sequential',
     1     file='hbt_log.out')
      open(unit=8,status='unknown',access='sequential',
     1     file='hbt_simulation.out')

CCC   Initialize Arrays and Data Structures:
      If(ALICE .eq. 1) then
C     In fact we not need to call initialization, 
C     because we can easily assume that is already done
         Call alihbtp_initialize
      Else If (ALICE .eq. 0) Then
         Call initialize
      End If

      Write(6,100)
CCC   Read Input Controls and Parameters:
      Call read_data(1)
      If(errorcode .eq. 1) Return
      
CCC   Setup values and check input parameter ranges and consistency:
      Call set_values
      If(errorcode .eq. 1) Return

CCC   Produce Basic Output File Header:
      Call write_data(1,0)
      If(errorcode .eq. 1) Return
      
      
      Write(6,101)
CCC   Read Event Input file and fill flag files:
      Call read_data(2)
      If(errorcode .eq. 1) Return
      
      Write(6,102)
CCC   Get the Reference Histograms and write out if new calculation:
      Call getref_hist
      If(errorcode .eq. 1) Return
      Call write_data(3,0)
      If(errorcode .eq. 1) Return
      Write(6,103)

      Write(6,104)
CCC   Compute the correlation model and print out:
      Call correl_model
      Call write_data(4,0)
      If(errorcode .eq. 1) Return

      Write(6,105)
CCC   Carry out the Track Adjustment/Correlation Fitting Procedure:
      Call correlation_fit
      Write(6,106)

CCC   Final Output of Inclusive Quantities:
      Call write_data(6,0)
      If(errorcode .eq. 1) Return
      
CCC   Close Output Files:
      close(unit=7)
      close(unit=8)

CCC   Formats:
100   Format(5x,'Read Input Controls, Setup values, check input:')
101   Format(5x,'Read Event Input file and fill flag files:')
102   Format(5x,'Get the Reference Histograms:')
103   Format(5x,'Finished with Reference Histograms:')
104   Format(5x,'Compute the correlation model:')
105   Format(5x,'Start Track Adjustment/Correlation Fitting Procedure:')
106   Format(5x,'Finished with Track Fitting Procedure:')

      Return
      END

C-------------------------------------------------------------------


      subroutine initialize
      implicit none

CCC   This subroutine sets all arrays and structures to zero:

      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'
      Include 'common_coulomb.inc'
      Include 'common_event_summary.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'
      Include 'common_particle.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,j,k

      do i = 1,trk_maxlen
         trk_id(i)            = 0
         trk_px_sec(i)        = 0
         trk_py_sec(i)        = 0
         trk_pz_sec(i)        = 0
         trk_sector(i)        = 0
         trk_flag(i)          = 0
         trk_out_flag(i)      = 0
         trk_merge_flag(i)    = 0
         trk_ge_pid(i)        = 0
         trk_start_vertex(i)  = 0
         trk_stop_vertex(i)   = 0
         trk_event_line(i)    = 0
         trk_px(i)            = 0.0
         trk_py(i)            = 0.0
         trk_pz(i)            = 0.0
         trk_E(i)             = 0.0
         trk_pt(i)            = 0.0
         trk_phi(i)           = 0.0
         trk_eta(i)           = 0.0
      end do
      
      do i = 1,trk2_maxlen
         trk2_id(i)            = 0
         trk2_px_sec(i)        = 0
         trk2_py_sec(i)        = 0
         trk2_pz_sec(i)        = 0
         trk2_sector(i)        = 0
         trk2_flag(i)          = 0
         trk2_out_flag(i)      = 0
         trk2_merge_flag(i)    = 0
         trk2_ge_pid(i)        = 0
         trk2_start_vertex(i)  = 0
         trk2_stop_vertex(i)   = 0
         trk2_event_line(i)    = 0
         trk2_px(i)            = 0.0
         trk2_py(i)            = 0.0
         trk2_pz(i)            = 0.0
         trk2_E(i)             = 0.0
         trk2_pt(i)            = 0.0
         trk2_phi(i)           = 0.0
         trk2_eta(i)           = 0.0
      end do

      do i = 1,sec_maxlen
         stm_sec_id(i)         = 0
         stm_n_trk_sec(i)      = 0
         stm_flag(i)           = 0
         do j = 1,max_trk_sec
            stm_track_id(j,i)  = 0
         end do
      end do

      do i = 1,sec_maxlen2
         stm2_sec_id(i)         = 0
         stm2_n_trk_sec(i)      = 0
         stm2_flag(i)           = 0
         do j = 1,max_trk_sec2
            stm2_track_id(j,i)  = 0
         end do
      end do

      do i = 1,part_maxlen
         part_id(i)       = 0
         part_charge(i)   = 0
         part_mass(i)     = 0.0
         part_lifetime(i) = 0.0
      end do

      do i = 1,max_trk_save
         old_sec_trkid(i) = 0
         new_sec_trkid(i) = 0
      end do

      do i = 1,max_h_1d
         hist_like_1d(i)   = 0
         hist_unlike_1d(i) = 0
         htmp_like_1d(i)   = 0
         htmp_unlike_1d(i) = 0
         href_like_1d(i)   = 0
         href_unlike_1d(i) = 0
         hinc_like_1d(i)   = 0
         hinc_unlike_1d(i) = 0
         hist1_pt_1(i)     = 0
         hist1_pt_2(i)     = 0
         hist1_phi_1(i)    = 0
         hist1_phi_2(i)    = 0
         hist1_eta_1(i)    = 0
         hist1_eta_2(i)    = 0
         htmp1_pt_1(i)     = 0
         htmp1_pt_2(i)     = 0
         htmp1_phi_1(i)    = 0
         htmp1_phi_2(i)    = 0
         htmp1_eta_1(i)    = 0
         htmp1_eta_2(i)    = 0
         href1_pt_1(i)     = 0
         href1_pt_2(i)     = 0
         href1_phi_1(i)    = 0
         href1_phi_2(i)    = 0
         href1_eta_1(i)    = 0
         href1_eta_2(i)    = 0
         hinc1_pt_1(i)     = 0
         hinc1_pt_2(i)     = 0
         hinc1_phi_1(i)    = 0
         hinc1_phi_2(i)    = 0
         hinc1_eta_1(i)    = 0
         hinc1_eta_2(i)    = 0
      end do

      do i = 1,max_h_3d
      do j = 1,max_h_3d
      do k = 1,max_h_3d
         hist_like_3d_fine(i,j,k)     = 0
         hist_unlike_3d_fine(i,j,k)   = 0
         hist_like_3d_coarse(i,j,k)   = 0
         hist_unlike_3d_coarse(i,j,k) = 0
         htmp_like_3d_fine(i,j,k)     = 0
         htmp_unlike_3d_fine(i,j,k)   = 0
         htmp_like_3d_coarse(i,j,k)   = 0
         htmp_unlike_3d_coarse(i,j,k) = 0
         href_like_3d_fine(i,j,k)     = 0
         href_unlike_3d_fine(i,j,k)   = 0
         href_like_3d_coarse(i,j,k)   = 0
         href_unlike_3d_coarse(i,j,k) = 0
         hinc_like_3d_fine(i,j,k)     = 0
         hinc_unlike_3d_fine(i,j,k)   = 0
         hinc_like_3d_coarse(i,j,k)   = 0
         hinc_unlike_3d_coarse(i,j,k) = 0
      end do
      end do
      end do 

      do i = 1,max_c2_1d
         c2mod_like_1d(i)   = 0.0
         c2mod_unlike_1d(i) = 0.0
         c2fit_like_1d(i)   = 0.0
         c2fit_unlike_1d(i) = 0.0
         c2err_like_1d(i)   = 0.0
         c2err_unlike_1d(i) = 0.0
      end do

      do i = 1,max_c2_3d
      do j = 1,max_c2_3d
      do k = 1,max_c2_3d
         c2mod_like_3d_fine(i,j,k)     = 0.0
         c2mod_unlike_3d_fine(i,j,k)   = 0.0
         c2mod_like_3d_coarse(i,j,k)   = 0.0
         c2mod_unlike_3d_coarse(i,j,k) = 0.0
         c2fit_like_3d_fine(i,j,k)     = 0.0
         c2fit_unlike_3d_fine(i,j,k)   = 0.0
         c2fit_like_3d_coarse(i,j,k)   = 0.0
         c2fit_unlike_3d_coarse(i,j,k) = 0.0
         c2err_like_3d_fine(i,j,k)     = 0.0
         c2err_unlike_3d_fine(i,j,k)   = 0.0
         c2err_like_3d_coarse(i,j,k)   = 0.0
         c2err_unlike_3d_coarse(i,j,k) = 0.0
      end do
      end do
      end do

      do i = 1,max_c2_coul
         c2_coul_like(i)   = 0.0
         c2_coul_unlike(i) = 0.0
         q_coul(i)         = 0.0
      end do

      do i = 1,max_events
         num_iter(i)                      = 0.0
         n_part_used_1_store(i)           = 0.0
         n_part_used_2_store(i)           = 0.0
         num_sec_flagged_store(i)         = 0.0
         frac_trks_out(i)                 = 0.0
         frac_trks_flag(i)                = 0.0
         chisq_like_1d_store(i)           = 0.0
         chisq_unlike_1d_store(i)         = 0.0
         chisq_like_3d_fine_store(i)      = 0.0
         chisq_unlike_3d_fine_store(i)    = 0.0
         chisq_like_3d_coarse_store(i)    = 0.0
         chisq_unlike_3d_coarse_store(i)  = 0.0
         chisq_hist1_1_store(i)           = 0.0
         chisq_hist1_2_store(i)           = 0.0
         chisq_total_store(i)             = 0.0
      end do

      Return
      END

C---------------------------------------------------------------------


      subroutine set_values
      implicit none

CCC   This subroutine sets parameters based on the main input.
CCC   The consistency of the input parameters and controls is
CCC   checked.  Any problems are reported in the Log File,
CCC   'hbt_log.out'.  Most inconsistencies or array size limit
CCC   overflows will cause the code execution to STOP.

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'
      Include 'common_coulomb.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'
      Include 'common_particle.inc'

CCC   Local Variable Type Declarations:

      integer*4 iphistep, ptmaxsteps, iptstep

      real*4 px1,py1,pz1,E1, pt1,phi1
      real*4 px2,py2,pz2,E2
      real*4 pt_step,phi_step
      real*4 pxstepmin, pxstepmax, pystepmin, pystepmax

CCC  Check Input Controls:

      if(ref_control .lt. 1 .or. ref_control .gt. 2) then
         write(7,101) ref_control
         errorcode = 1
         Return
      end if

      if(switch_1d .lt. 0 .or. switch_1d .gt. 3) then
         write(7,102) switch_1d
         errorcode = 1
         Return
      end if

      if(switch_3d .lt. 0 .or. switch_3d .gt. 2) then
         write(7,103) switch_3d
         errorcode = 1
         Return
      end if

      if(switch_type .lt. 1 .or. switch_type .gt. 3) then
         write(7,104) switch_type
         errorcode = 1
         Return
      end if

      if(switch_coherence .lt. 0 .or. switch_coherence .gt. 1) then
         write(7,105) switch_coherence
         errorcode = 1
         Return
      end if

      if(switch_coulomb .lt. 0 .or. switch_coulomb .gt. 3) then
         write(7,106) switch_coulomb
         errorcode = 1
         Return
      end if

      if(switch_fermi_bose.ne.-1 .and. switch_fermi_bose.ne.1) then
         write(7,107) switch_fermi_bose
         errorcode = 1
         Return
      end if 

      if(n_pid_types .lt. 1 .or. n_pid_types .gt. 2) then
         write(7,108) n_pid_types
         errorcode = 1
         Return
      end if

      if(switch_type .ge. 2 .and. n_pid_types .eq. 1) then
         write(7,109) switch_type, n_pid_types
         errorcode = 1
         Return
      end if

      if(n_pid_types .eq. 1) then
         if(pid(1).gt.0 .and. pid(2).gt.0) then
            write(7,1091) pid(1),pid(2)
            errorcode = 1
            Return
         end if
      end if

      if(pid(1).eq.0 .and. pid(2).eq.0) then
         write(7,1092)
         errorcode = 1
         Return
      end if

      if(n_pid_types .eq. 2) then
         if(pid(1).gt.0.and.pid(2).gt.0.and.pid(1).ne.pid(2))then
            continue
         else
            write(7,1093) pid(1), pid(2)
            errorcode = 1
            Return
         end if
      end if

      if(pid(1).gt.0.and.pid(2).gt.0.and.pid(1).eq.pid(2))then
         write(7,1094) pid(1), pid(2)
         errorcode = 1
         Return
      end if

      if(trk_accep .le. 0.0) then
         write(7,10941) trk_accep
         errorcode = 1
         Return

      end if

CCC   Check Input Parameters:

      if(deltap .le. 0.0) deltap = 0.1
      if(maxit  .lt. 0  ) maxit = 50
      if(delchi .lt. 0.0) delchi = 0.1
      if(irand  .le. 0  ) irand = 12345

CCC   Check Coulomb source radius in range for Pratt type Coulomb correction.

      if(switch_coulomb .eq. 3 .and. (Q0 .lt. coulradmin .or.
     1   Q0 .gt. coulradmax)) then
         write(7,132) Q0
         errorcode = 1
         Return
      end if

CCC   Load the Pratt type Coulomb correction if this form is selected.

      if(switch_coulomb .eq. 3 .and. (Q0 .ge. coulradmin .and.
     1   Q0 .le. coulradmax)) then
         Call read_data(6)
      end if

CCC   Check and determine the one-body distribution's binning:

      if(n_pt_bins .lt. 1 .or. n_pt_bins .gt. max_h_1d) then
         write(7,110) n_pt_bins
         errorcode = 1
         Return
      end if

      if(n_phi_bins .lt. 1 .or. n_phi_bins .gt. max_h_1d) then
         write(7,111) n_phi_bins
         errorcode = 1
         Return
      end if

      if(n_eta_bins .lt. 1 .or. n_eta_bins .gt. max_h_1d) then
         write(7,112) n_eta_bins
         errorcode = 1
         Return
      end if

      if(pt_min .gt. pt_max .or. pt_min .lt. 0.0) then
         write(7,113) pt_min, pt_max
         errorcode = 1
         Return
      end if

      if(phi_min.gt.phi_max .or. phi_min.lt.0.0 .or.
     1   phi_max.gt.360.0) then
         write(7,114) phi_min, phi_max
         errorcode = 1
         Return
      end if

      if(eta_min .gt. eta_max) then
         write(7,115) eta_min, eta_max
         errorcode = 1
         Return
      end if

      pt_bin_size  = (pt_max  - pt_min )/float(n_pt_bins)
      phi_bin_size = (phi_max - phi_min)/float(n_phi_bins)
      eta_bin_size = (eta_max - eta_min)/float(n_eta_bins)

CCC   Check and determine the two-body distribution's binning:

      n_1d_total = n_1d_fine + n_1d_coarse
      n_3d_total = n_3d_fine + n_3d_coarse - 1

      if(switch_1d .gt. 0) then
      if(n_1d_fine .lt. 1) then
         write(7,116) n_1d_fine
         errorcode = 1
         Return
      end if

      if(n_1d_coarse .lt. 1) then
         write(7,117) n_1d_coarse
         errorcode = 1
         Return
      end if

      if(n_1d_total .gt. max_h_1d) then
         write(7,118) n_1d_total
         errorcode = 1
         Return
      end if

      qmid_1d = binsize_1d_fine  *float(n_1d_fine)
      qmax_1d = binsize_1d_coarse*float(n_1d_coarse) + qmid_1d
      end if

      if(switch_3d .gt. 0) then
      if(n_3d_fine .lt. 1 .or. n_3d_fine .gt. max_h_3d) then
         write(7,119) n_3d_fine
         errorcode = 1
         Return
      end if

      if(n_3d_coarse .lt. 1 .or. n_3d_coarse .gt. max_h_3d) then
         write(7,120) n_3d_coarse
         errorcode = 1
         Return
      end if

      qmid_3d = binsize_3d_fine  *float(n_3d_fine)
      qmax_3d = binsize_3d_coarse*float(n_3d_coarse)

      if(abs(qmid_3d - binsize_3d_coarse) .gt. 0.00001) then
         write(7,121) qmid_3d, binsize_3d_coarse
         errorcode = 1
         Return
      end if

      if(n_3d_fine_project .gt. n_3d_fine) then
         write(7,1211) n_3d_fine_project, n_3d_fine
         n_3d_fine_project = n_3d_fine
      end if

      if(n_3d_fine_project .lt. 1) then
         write(7,1212) n_3d_fine_project
         n_3d_fine_project = 1
      end if
      end if

CCC   Check and determine Track-Sector Parameters:

      if(n_px_bins .lt. 1) then
         write(7,122) n_px_bins
         errorcode = 1
         Return
      end if

      if(n_py_bins .lt. 1) then
         write(7,123) n_py_bins
         errorcode = 1
         Return
      end if

      if(n_pz_bins .lt. 1) then
         write(7,124) n_pz_bins
         errorcode = 1
         Return
      end if

      n_sectors = n_px_bins * n_py_bins * n_pz_bins
      if(n_sectors .gt. sec_maxlen) then
         write(7,125) n_sectors
         errorcode = 1
         Return
      end if
 
      if(n_sectors .gt. sec_maxlen2 .and. ref_control .eq. 2) then
         write(7,1251) n_sectors
         errorcode = 1
         Return
      end if

      if(trk_maxlen .ne. trk2_maxlen .and. ref_control .eq. 2) then
         write(7,1252)
         errorcode = 1
         Return
      end if

      if(max_trk_save .ne. max_trk_sec .or.
     1   max_trk_save .ne. max_trk_sec2 .or.
     2   max_trk_sec  .ne. max_trk_sec2) then
         write(7,12521) max_trk_save,max_trk_sec,max_trk_sec2
         errorcode = 1
         Return
      end if

      delpx = (px_max - px_min)/float(n_px_bins)
      delpy = (py_max - py_min)/float(n_py_bins)
      delpz = (pz_max - pz_min)/float(n_pz_bins)

CCC   Check that the Track-Sector Grid includes the acceptance range:
CCC   The Track-Sector Grid is a 3D {px,py,pz} box, while the acceptance
CCC   is defined in cylindrical {pt,phi,eta} coordinates.
CCC
CCC   Check the z-momentum components:

      if(eta_min .ge. 0.0) then
         Call Hbtp_kin(px1,py1,pz1,E1,pt_min,0.0,eta_min,0.14,2)
         Call Hbtp_kin(px2,py2,pz2,E2,pt_max,0.0,eta_max,0.14,2)
      else if(eta_max .le. 0.0) then
         Call Hbtp_kin(px1,py1,pz1,E1,pt_max,0.0,eta_min,0.14,2)
         Call Hbtp_kin(px2,py2,pz2,E2,pt_min,0.0,eta_max,0.14,2)
      else if(eta_min .le. 0.0 .and. eta_max .ge. 0.0) then
         Call Hbtp_kin(px1,py1,pz1,E1,pt_max,0.0,eta_min,0.14,2)
         Call Hbtp_kin(px2,py2,pz2,E2,pt_max,0.0,eta_max,0.14,2)
      end if

      if(pz1 .lt. pz_min .or. pz2 .gt. pz_max) then
         write(7,126) pz1,pz_min,pz2,pz_max
         errorcode = 1
         Return
      end if

CCC   Check the x,y-momentum components by scanning over the perimeter
CCC   of the (pt,phi) acceptance domain space with 100 trial grid points.
CCC   The overall required px_min, px_max, py_min and py_max to cover the
CCC   acceptance by the track-sectors is determined.  These values are
CCC   then compared with the min/max px and py ranges for the track-sectors.
CCC
CCC   Divide the pt and phi acceptance ranges into 24 equal steps:

      pt_step  = (pt_max  - pt_min)/24.0
      phi_step = (phi_max - phi_min)/24.0
      pxstepmax = -1000.
      pxstepmin =  1000.
      pystepmax = -1000.
      pystepmin =  1000.
      phi1 = phi_min - phi_step
      do iphistep = 1,25
         phi1 = phi1 + phi_step
         ptmaxsteps = 2
         if(iphistep.eq.1 .or. iphistep.eq.25) ptmaxsteps = 25
         pt1 = pt_min - pt_step
         do iptstep = 1,ptmaxsteps
            if(iphistep.eq.1 .or. iphistep.eq.25) then
               pt1 = pt1 + pt_step
            else if(iphistep.gt.1 .and. iphistep.lt.25) then
               if(iptstep.eq.1) pt1 = pt_min
               if(iptstep.eq.2) pt1 = pt_max
            end if
            Call Hbtp_kin(px1,py1,pz1,E1,pt1,phi1,0.0,0.14,2)
               if(px1.gt.pxstepmax) pxstepmax = px1
               if(px1.lt.pxstepmin) pxstepmin = px1
               if(py1.gt.pystepmax) pystepmax = py1
               if(py1.lt.pystepmin) pystepmin = py1
         end do
      end do

      if(pxstepmin .lt. px_min .or. pxstepmax .gt. px_max) then
         write(7,127) pxstepmin,px_min,pxstepmax,px_max
         errorcode = 1
         Return
      end if

      if(pystepmin .lt. py_min .or. pystepmax .gt. py_max) then
         write(7,128) pystepmin,py_min,pystepmax,py_max
         errorcode = 1
         Return
      end if

CCC   Load Geant Particle Data:
      Call Hbtp_particle_prop

CCC   Check Requested Particle ID Numbers:

      if(n_pid_types.eq.1 .and. pid(1).le.0 .and. pid(2).le.0) then
         write(7,131) pid(1),pid(2)
         errorcode = 1
         Return
      end if

CCC   Initialize Masses to 0.0
    
      mass1 = 0.0
      mass2 = 0.0

      if(n_pid_types .eq. 1 .and. pid(1) .ne. 0) then
         if(pid(1) .lt. 1 .or. pid(1) .gt. part_maxlen) then
            write(7,129) pid(1)
            errorcode = 1
            Return
         else
            mass1 = part_mass(pid(1))
         end if
      else if(n_pid_types .eq. 1 .and. pid(2) .ne. 0) then
         if(pid(2) .lt. 1 .or. pid(2) .gt. part_maxlen) then
            write(7,130) pid(2)
            errorcode = 1
            Return
         else
            mass2 = part_mass(pid(2))
         end if
      else if(n_pid_types .eq. 2) then
         if(pid(1) .lt. 1 .or. pid(1) .gt. part_maxlen) then
            write(7,129) pid(1)
            errorcode = 1
            Return
         else
            mass1 = part_mass(pid(1))
         end if
         if(pid(2) .lt. 1 .or. pid(2) .gt. part_maxlen) then
            write(7,130) pid(2)
            errorcode = 1
            Return
         else
            mass2 = part_mass(pid(2))
         end if
      end if

CCC   Set Math Constants:

      pi = 3.141592654
      hbc = 0.19732891
      rad = 180.0/pi

CCC   FORMATS:

101   Format(5x,'ref_control = ',I5,'Out of Range - STOP')
102   Format(5x,'switch_1d   = ',I5,'Out of Range - STOP')
103   Format(5x,'switch_3d   = ',I5,'Out of Range - STOP')
104   Format(5x,'switch_type = ',I5,'Out of Range - STOP')
105   Format(5x,'switch_coherence = ',I5,'Out of Range - STOP')
106   Format(5x,'switch_coulomb = ',I5,'Out of Range - STOP')
107   Format(5x,'switch_fermi_bose = ',I5,'Out of Range - STOP')
108   Format(5x,'n_pid_types = ',I5,'Out of Range - STOP')
109   Format(5x,'switch_type & n_pid_types = ',2I5,
     1       'Incompatible - STOP')
1091  Format(5x,'For n_pid_types=1, cannot have both PID#1,2 = ',
     1       2I5,' .ne.0 - STOP')
1092  Format(5x,'Both PID# 1 and 2 = 0, - STOP')
1093  Format(5x,'For n_pid_types=2, PID#1,2 = ',2I5,
     1      ' are incorrect - STOP')
1094  Format(5x,'Both PID# 1,2 = ',2I5,' are equal - STOP')
10941 Format(5x,'Track acceptance output frac .le. 0.0 - STOP')
132   Format(5x,'Coulomb source radius = ',E12.4,' - For Pratt ',
     1       'Correction, Out of Range - STOP')
110   Format(5x,'# pt bins  = ',I5,'Out of Range - STOP')
111   Format(5x,'# phi bins = ',I5,'Out of Range - STOP')
112   Format(5x,'# eta bins = ',I5,'Out of Range - STOP')
113   Format(5x,'pt min/max accept. range = ',2E12.4,' is wrong-STOP')
114   Format(5x,'phi min/max accept. range = ',2E12.4,' is wrong-STOP')
115   Format(5x,'eta min/max accept. range = ',2E12.4,' is wrong-STOP')
116   Format(5x,'# 1d fine mesh for C2 = ',I5,' .lt.1 - STOP')
117   Format(5x,'# 1d coarse mesh for C2 = ',I5,' .lt.1 - STOP')
118   Format(5x,'Total # 1d mesh for C2 = ',I5,' .gt.max_h_1d - STOP')
119   Format(5x,'# 3d fine mesh for C2 = ',I5,'Out of Range - STOP')
120   Format(5x,'# 3d coarse mesh for C2 = ',I5,'Out of Range - STOP')
121   Format(5x,'3D C2 fine range & coarse grid = ',2E12.4,
     1       'Not Equal - STOP')
1211  Format(5x,'# 3D fine bins projected = ',I5,
     1       ' TOO BIG - reduce to n_3d_fine = ',I5)
1212  Format(5x,'# 3D fine bins projected = ',I5,
     1       ' Set to 1')
122   Format(5x,'#track-sector px bins = ',I5,' .lt.1 - STOP')
123   Format(5x,'#track-sector py bins = ',I5,' .lt.1 - STOP')
124   Format(5x,'#track-sector pz bins = ',I5,' .lt.1 - STOP')
125   Format(5x,'Total # trk-sec = ',I5,' .gt.sec_maxlen - STOP')
1251  Format(5x,'Total # trk-sec = ',I5,' .gt.sec_maxlen2 for ',
     1       'Reference calc. - STOP')
1252  Format(5x,'trk_maxlen .ne. trk2_maxlen for Ref. Calc. - STOP')
12521 Format(5x,'max_trk_save,max_trk_sec,max_trk_sec2 = ',
     1       3I5,' are not all equal - STOP')
126   Format(5x,'pz accept. not covered by sectors-STOP:',4E12.4)
127   Format(5x,'px accept. not covered by sectors-STOP:',4E12.4)
128   Format(5x,'py accept. not covered by sectors-STOP:',4E12.4)
131   Format(5x,'Particle ID values = ',2I5,' not valid - STOP')
129   Format(5x,'Particle ID value #1 = ',I5,' not valid - STOP')
130   Format(5x,'Particle ID value #2 = ',I5,' not valid - STOP')

      Return
      END

C----------------------------------------------------------------------


      subroutine read_data(mode)
      implicit none

CCC   This subroutine does all the data input associated with all input
C     files.  Some diagnostic output is printed here if errors occur
C     during the file reading.  Two auxiliary output files, which tag
C     the events input tracks are written out.
C
C     The type of input is controlled by the value of 'mode'
C     where:
C     (The following mostly applies to the standalone application
C     that reads from files and writes temporary scratch files.
C     This is the ALICE=0 mode.)
C
C        mode = 1, read the parameter and switches input file
C
C        mode = 2, scan the event text file and write out two
C                  auxiliary output/tag files; select and mark
C                  accepted tracks to use.
C
C        mode = 3, read the reference pair and one-body histograms
C
C        mode = 4, read the next event from the event text file, 
C                  'event_text.in,' and load the accepted tracks
C                  into the 'trk' data structure.
C
C        mode = 5, same as mode=4, except the accepted tracks are
C                  loaded into the 'trk2' data structure.
C
C        mode = 6, read the input Coulomb correction tables and
C                  interpolate for the requested source radius, arrays
C                  in common/coulomb/ are filled for like and unlike
C                  charged pairs.
C
C        mode = 7, read the next event from the event text file, 
C                  'event_text.in,' and load the accepted tracks
C                  into the 'trk' data structure. Then copy the event
C                  data in 'event_text.in' to 'event_text_aux.in' and
C                  from 'event_tracks.select' to 'event_tracks_aux.select'
C
C        mode = 8, read contents of 'event_text_aux.in' using flag values
C                  in 'event_tracks_aux.select' and copy into
C                  'event_hbt_text.out' (i.e. the main event output file)
C                  where the momentum values for accepted tracks are
C                  overwritten with the adjusted (correlated) parameters
C                  in the 'trk' data structure.
C
C                  If trk_accep .lt. 1.0, then only write this fraction
C                  of the final tracks, as determined by a random number
C                  throw.
C
C     Summary of Files:
C     ----------------
C
C File Unit #        Filename               Description
C ---------------------------------------------------------------------------
C     1       hbt_parameters.in          Input switches, parameters
C     2       event_text.in              Event Gen output, GSTAR text format
C     3       event_line.flags           Event file line flags
C     4       event_tracks.select        Event file selected tracks
C     7       hbt_log.out                Log and error messages
C     8       hbt_simulation.out         Full Output
C     9       hbt_pair_reference.hist    Reference pair histograms
C    10       event_hbt_text.out         Updated/correlated event text file
C    11       hbt_singles_reference.hist Reference one-body histograms
C    12       event_text_aux.in          Tmp. copy of 'event_text.in'/event
C    14       event_tracks_aux.select    Tmp. copy 'event_tracks.select'/event
C   21-27     cpp_*.dat (*=06,08,...18)  Like pair Pratt Coul. Correct
C   31-37     cpm_*.dat (*=06,08,...18)  Unlike pair Pratt Coul. Correct
C ---------------------------------------------------------------------------
C

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'
      Include 'common_coulomb.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_particle.inc'

      integer LNBLNK  

CCC   Local Variable Type Declarations:

      real*4 px,py,pz,E,pt,phi,eta,mass
      real*4 acheck(10), function(20)
      real*4 hbtpran

      integer*4 i,j,k,mode,flag,flag4,flag0,ntracks
      integer*4 ge_pid,tid,start_v,stop_v,eg_pid
      integer*4 ref_check,pidok,accepok,check(13)
      integer*4 event_counter,track_counter
      integer*4 track_counter_1,track_counter_2

      character*5 evg_label,event_line,vertex_line,track_line,dummy
      character*5 gener_line
      character*80 comment_event_label
      character*87 vertex_label
      character*93 gener_label

      parameter (event_line  = 'EVENT')
      parameter (vertex_line = 'VERTE')
      parameter (track_line  = 'TRACK')
      parameter (gener_line  = 'GENER')
      parameter (flag4 = 4)
      parameter (flag0 = 0)
C     ALICE USE ONLY 
      CHARACTER*100 CHROOT  
      CHARACTER*100 FILNAM
      INTEGER*4 LNROOT
      LOGICAL EXISTS
      CHROOT=' '
C     
      
CCC   Begin Data Input Options:

C------------------------
      IF (mode.eq.1) THEN              ! Read Input parameters from File#1
C------------------------
      
CCC   For standalone version (ALICE = 0), read parameters from
CCC   File Unit=1, 'hbt_parameters.in'
CCC   For ALICE-ROOT version (ALICE=1) load parameters from Call to C++ funct
      If(ALICE .eq. 1) then
         Call AliHbtp_SetParameters
      Else If(ALICE .eq. 0) Then

         open(unit=1,status='old',access='sequential',
     1        file='hbt_parameters.in')

CCC   Read Control Switches:  (See Main program listing for complete
CCC                            description of input parameters)

      read(1,*) ref_control
      read(1,*) switch_1d
      read(1,*) switch_3d
      read(1,*) switch_type
      read(1,*) switch_coherence
      read(1,*) switch_coulomb
      read(1,*) switch_fermi_bose
      read(1,*) trk_accep
      read(1,*) print_full
      read(1,*) print_sector_data

CCC   Read Parameters:

      read(1,*) n_pid_types
      read(1,*) pid(1),pid(2)
      read(1,*) deltap
      read(1,*) maxit
      read(1,*) delchi
      read(1,*) irand

CCC   Read Source Parameters:

      read(1,*) lambda
      read(1,*) R_1d
      read(1,*) Rside, Rout, Rlong
      read(1,*) Rperp, Rparallel, R0
      read(1,*) Q0

CCC   Read one-body {pt,phi,eta} bins:

      read(1,*) n_pt_bins ,pt_min ,pt_max
      read(1,*) n_phi_bins,phi_min,phi_max
      read(1,*) n_eta_bins,eta_min,eta_max

CCC   Read two-body 1D and 3D bins:

      read(1,*) n_1d_fine,   binsize_1d_fine
      read(1,*) n_1d_coarse, binsize_1d_coarse
      read(1,*) n_3d_fine,   binsize_3d_fine
      read(1,*) n_3d_coarse, binsize_3d_coarse
      read(1,*) n_3d_fine_project

CCC   Read momentum space track sector bins in {px,py,pz}:

      read(1,*) n_px_bins,px_min,px_max
      read(1,*) n_py_bins,py_min,py_max
      read(1,*) n_pz_bins,pz_min,pz_max

CCC   Relative Chi-Square weights for track adjustment fitting:

      read(1,*) chisq_wt_like_1d
      read(1,*) chisq_wt_unlike_1d
      read(1,*) chisq_wt_like_3d_fine
      read(1,*) chisq_wt_unlike_3d_fine
      read(1,*) chisq_wt_like_3d_coarse
      read(1,*) chisq_wt_unlike_3d_coarse
      read(1,*) chisq_wt_hist1_1
      read(1,*) chisq_wt_hist1_2

      Close(unit=1)
      End If  !  ALICE Data I/O Option

C-----------------------------
      ELSE IF (mode.eq.2) THEN    
C-----------------------------

C     Open event generator text file, 'event_text.in,' and read it,
C     noting each type of line input.  Write out a file called 
C     'event_line.flags' which identifies the type of information on
C     each line where:
C
C        'EVENT:'  lines are assigned flag = 1
C        'VERTEX:' lines are assigned flag = 2
C        'TRACK:'  lines are assigned flag = 3
C        'GENER:'  lines are assigned flag = 5
C        All other lines are assigned flag = 0

      If(ALICE .eq. 0) Then
      open(unit=2,status='old',access='sequential',
     1     file='event_text.in')
      open(unit=3,status='unknown',access='sequential',
     1     file='event_line.flags')

CCC   Set Event Counter:

      event_counter = 0
30    read(2,10,err=20,end=25) evg_label
10    Format(A)
      if(evg_label .eq. event_line) then
         event_counter = event_counter + 1
         flag = 1
      else if(evg_label .eq. vertex_line) then
         flag = 2
      else if(evg_label .eq. track_line) then
         flag = 3
      else if(evg_label .eq. gener_line) then
         flag = 5
      else
         flag = 0
      end if

      write(3,11) flag
11    Format(1x,I1)
      go to 30         ! Return to S.N. 30 and read next line in file
20    write(7,12) event_counter
12    Format(5x,'Read error in event_text.in at event# ',I5,' - STOP')
      Stop
25    Continue
      Close(unit=2)
      Close(unit=3)
      End If  !  ALICE Data I/O Option

C     Next, re-open the 'event_text.in' and 'event_line.flags' files
C     again and read thru the entire files.  For each track, check its'
C     particle ID and kinematics (pt,phi,eta) with respect to the
C     selected particle ID type(s) for the run and the acceptances.
C     Fill another file called, 'event_tracks.select,' which is the same
C     as 'event_line.flags' except that the accepted tracks are marked
C     with flag = 4.
C
C     NOTE:  Assume all vertices in 'event_text.in' are at microscopic
C            distances (fermis) such that all particles in the event
C            file are considered as primaries.  Also for each event
C            the code will only accept tracks up to the limit imposed
C            by trk_maxlen in the 'trk' table.

      If(ALICE .eq. 1) Then
CCC   For ALICE application do the following:
CCC      Store number of events in 'n_events'
CCC      Count number accepted tracks in each event, check wrt trk_maxlen
CCC      Mark accepted tracks in all events

         Call AliHbtp_GetNumberEvents(n_events)
         do i = 1,n_events
            Call AliHbtp_SetActiveEventNumber(i)
            track_counter = 0
            Call AliHbtp_GetNumberTracks(ntracks)
            do j = 1,ntracks
               Call AliHbtp_GetTrack(j,flag,px,py,pz,ge_pid)
               eg_pid = ge_pid
CCC   Check if this track's particle ID is one to be used

         pidok = 0
         accepok = 0
         if(pid(1).gt.0 .and. eg_pid.eq.pid(1)) pidok = 1
         if(pid(2).gt.0 .and. eg_pid.eq.pid(2)) pidok = 1
         if(pidok.eq.1 .and. eg_pid.le.part_maxlen) then
            mass = part_mass(eg_pid)
            Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
            if(pt.ge.pt_min .and. pt.le.pt_max .and.
     1         phi.ge.phi_min .and. phi.le.phi_max .and.
     2         eta.ge.eta_min .and. eta.le.eta_max) then
               if(track_counter .lt. trk_maxlen) then
                  accepok = 1
               else
                  write(7,62) trk_maxlen, event_counter
62                Format(5x,'#tracks exceeds trk_maxlen = ',
     1                  I6,' for event#',I4)
               end if
            end if
         end if

         if(pidok.eq.1 .and. accepok.eq.1) then
            track_counter = track_counter + 1
C	    write(*,*) '  FFF: 1 calling PutTrack j = ',j
            Call AliHbtp_PutTrack(j,flag4,px,py,pz,ge_pid)
         else
C	    write(*,*) '  FFF: 2 calling PutTrack j = ',j
            Call AliHbtp_PutTrack(j,flag0,px,py,pz,ge_pid)
         end if
         end do
         end do
      
      Else If(ALICE .eq. 0) Then

      open(unit=2,status='old',access='sequential',
     1     file='event_text.in')
      open(unit=3,status='old',access='sequential',
     1     file='event_line.flags')
      open(unit=4,status='unknown',access='sequential',
     1     file='event_tracks.select')

CCC   Set Event Counter:

      event_counter = 0
40    read(3,11,err=45,end=50) flag
      if(flag.eq.1) then
         event_counter = event_counter + 1
         track_counter = 0
      end if

      if(flag.ne.3) then
         read(2,10) dummy
         write(4,11) flag
      else if(flag.eq.3) then
         read(2,41,err=46,end=50) ge_pid,px,py,pz,tid,start_v,
     1                            stop_v,eg_pid
41    Format(7x,I6,3(1x,G12.5),4(1x,I6))

CCC   Check if the 'event_text.in' file includes non-zero PID
CCC   values for the variable 'eg_pid'.  If this is zero, then
CCC   use the ge_pid value.
         if(eg_pid.eq.0 .and. ge_pid.ne.0) eg_pid = ge_pid

CCC   Check if this track's particle ID is one to be used

         pidok = 0
         accepok = 0
         if(pid(1).gt.0 .and. eg_pid.eq.pid(1)) pidok = 1
         if(pid(2).gt.0 .and. eg_pid.eq.pid(2)) pidok = 1
         if(pidok.eq.1 .and. eg_pid.le.part_maxlen) then
            mass = part_mass(eg_pid)
            Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
            if(pt.ge.pt_min .and. pt.le.pt_max .and.
     1         phi.ge.phi_min .and. phi.le.phi_max .and.
     2         eta.ge.eta_min .and. eta.le.eta_max) then
               if(track_counter .lt. trk_maxlen) then
                  accepok = 1
               else
                  write(7,621) trk_maxlen, event_counter
621                Format(5x,'#tracks exceeds trk_maxlen = ',
     1                  I6,' for event#',I4)
               end if
            end if
         end if

         if(pidok.eq.1 .and. accepok.eq.1) then
            track_counter = track_counter + 1
            write(4,11) flag4
         else
            write(4,11) flag
         end if
    
      end if    ! End Flag=3 options

      go to 40  ! Return to S.N. 40 and read next record
45    write(7,60) event_counter
60    Format(5x,'Read error in event_line.flags at event#',I5,
     1       ' - STOP')
      Stop
46    write(7,61) event_counter
61    Format(5x,'Read error in event_text.in (2nd pass) at event#',I5,
     1       ' - STOP')
      Stop
50    Continue
      
      n_events = event_counter - 1  !  Set # events in event_text.in file
C                                   !  This is one less than the counter
C                                   !  value since the last 'EVENT:' line is
C                                   !  used to mark the End-of-File.

      Close(unit=2)
      Close(unit=3)
      Close(unit=4)

      End If    !  ALICE Data I/O Option

C-----------------------------
      ELSE IF(mode.eq.3) THEN
C-----------------------------

C     Read the reference histograms for pairs, then for singles for one
C     or two particle ID types.  Check switches, bins and mesh information
C     to be sure the input reference histograms are compatible with the
C     present run conditions.

      open(unit=9,status='old',access='sequential',
     1     file='hbt_pair_reference.hist')

      read(9,*) (check(i),i=1,3)
      read(9,*)  check(4),acheck(1),acheck(2)
      read(9,*)  check(5),acheck(3),acheck(4)
      read(9,*)  check(6),acheck(5),acheck(6)
      read(9,*) (check(i),i=7,9)
      read(9,*) (check(i),i=10,13)
      read(9,*) (acheck(i),i=7,10)
      read(9,*) num_pairs_like_ref, num_pairs_unlike_ref

CCC   Determine if the Input Reference pair histograms are compatible
CCC   with the present run parameters:

      ref_check = 1
      if(check(1)  .ne. n_pid_types ) ref_check = 0
      if(check(2)  .ne. pid(1)      ) ref_check = 0
      if(check(3)  .ne. pid(2)      ) ref_check = 0
      if(check(4)  .ne. n_pt_bins   ) ref_check = 0
      if(check(5)  .ne. n_phi_bins  ) ref_check = 0
      if(check(6)  .ne. n_eta_bins  ) ref_check = 0
      if(check(7)  .ne. switch_1d   ) ref_check = 0
      if(check(8)  .ne. switch_3d   ) ref_check = 0
      if(check(9)  .ne. switch_type ) ref_check = 0
      if(check(10) .ne. n_1d_fine   ) ref_check = 0
      if(check(11) .ne. n_1d_coarse ) ref_check = 0
      if(check(12) .ne. n_3d_fine   ) ref_check = 0
      if(check(13) .ne. n_3d_coarse ) ref_check = 0

      if(abs(acheck( 1) - pt_min )         .gt. 0.000001) ref_check = 0
      if(abs(acheck( 2) - pt_max )         .gt. 0.000001) ref_check = 0
      if(abs(acheck( 3) - phi_min )        .gt. 0.000001) ref_check = 0
      if(abs(acheck( 4) - phi_max )        .gt. 0.000001) ref_check = 0
      if(abs(acheck( 5) - eta_min )        .gt. 0.000001) ref_check = 0
      if(abs(acheck( 6) - eta_max )        .gt. 0.000001) ref_check = 0
      if(abs(acheck( 7) - binsize_1d_fine) .gt. 0.000001) ref_check = 0
      if(abs(acheck( 8) - binsize_1d_coarse).gt.0.000001) ref_check = 0
      if(abs(acheck( 9) - binsize_3d_fine) .gt. 0.000001) ref_check = 0
      if(abs(acheck(10) - binsize_3d_coarse).gt.0.000001) ref_check = 0

      if(ref_check .eq. 0) then
         write(7,70)
70    Format(5x,'Reference Pair Histogram Parameters not consistent',
     1       ' with present run conditions - STOP')
         errorcode = 1
         Return
      else if(ref_check .eq. 1) then

         if(switch_1d.gt.0 .and. n_1d_total.gt.0) then
            if(switch_type.eq.1 .or. switch_type.eq.3) then
               read(9,*) (href_like_1d(i),i=1,n_1d_total)
            end if
            if(switch_type.eq.2 .or. switch_type.eq.3) then
               read(9,*) (href_unlike_1d(i),i=1,n_1d_total)
            end if
         end if   ! End 1D input option

         if(switch_3d.gt.0) then
            if(switch_type.eq.1 .or. switch_type.eq.3) then

               if(n_3d_fine.gt.0) then
                  do i = 1,n_3d_fine
                  do j = 1,n_3d_fine
                  do k = 1,n_3d_fine
                     read(9,*) href_like_3d_fine(i,j,k)
                  end do
                  end do
                  end do
               end if

               if(n_3d_coarse.gt.0) then
                  do i = 1,n_3d_coarse
                  do j = 1,n_3d_coarse
                  do k = 1,n_3d_coarse
                     read(9,*) href_like_3d_coarse(i,j,k)
                  end do
                  end do
                  end do
               end if

            end if

            if(switch_type.eq.2 .or. switch_type.eq.3) then

               if(n_3d_fine.gt.0) then
                  do i = 1,n_3d_fine
                  do j = 1,n_3d_fine
                  do k = 1,n_3d_fine
                     read(9,*) href_unlike_3d_fine(i,j,k)
                  end do
                  end do
                  end do
               end if

               if(n_3d_coarse.gt.0) then
                  do i = 1,n_3d_coarse
                  do j = 1,n_3d_coarse
                  do k = 1,n_3d_coarse
                     read(9,*) href_unlike_3d_coarse(i,j,k)
                  end do
                  end do
                  end do
               end if

            end if

         end if     ! End 3D input option
      end if        ! End Reference Check OK/Not OK Option

      Close(unit=9)

CCC   Next read the one-body histograms for 1 or 2 particle ID types:

      open(unit=11,status='old',access='sequential',
     1     file='hbt_singles_reference.hist')

      read(11,*) (check(i),i=1,3)
      read(11,*)  check(4),acheck(1),acheck(2)
      read(11,*)  check(5),acheck(3),acheck(4)
      read(11,*)  check(6),acheck(5),acheck(6)
      read(11,*) n_part_used_1_ref, n_part_used_2_ref

CCC   Determine if Reference one-body histograms are compatible with
CCC   the present run conditions.

      ref_check = 1
      if(check(1) .ne. n_pid_types)   ref_check = 0
      if(check(2) .ne. pid(1)     )   ref_check = 0
      if(check(3) .ne. pid(2)     )   ref_check = 0
      if(check(4) .ne. n_pt_bins  )   ref_check = 0
      if(check(5) .ne. n_phi_bins )   ref_check = 0
      if(check(6) .ne. n_eta_bins )   ref_check = 0

      if(abs(acheck(1) - pt_min  ).gt.0.000001) ref_check = 0
      if(abs(acheck(2) - pt_max  ).gt.0.000001) ref_check = 0
      if(abs(acheck(3) - phi_min ).gt.0.000001) ref_check = 0
      if(abs(acheck(4) - phi_max ).gt.0.000001) ref_check = 0
      if(abs(acheck(5) - eta_min ).gt.0.000001) ref_check = 0
      if(abs(acheck(6) - eta_max ).gt.0.000001) ref_check = 0

      if(ref_check .eq. 0) then
         write(7,71)
71       Format(5x,'Reference One-Body Histogram parameters not ',
     1          'consistent with current run - STOP')
         errorcode = 1
         Return
      else if(ref_check .eq. 1) then

         if(pid(1).gt.0) then
            read(11,*) (href1_pt_1(i) ,i=1,n_pt_bins)
            read(11,*) (href1_phi_1(i),i=1,n_phi_bins)
            read(11,*) (href1_eta_1(i),i=1,n_eta_bins)
         end if

         if(pid(2).gt.0) then
            read(11,*) (href1_pt_2(i) ,i=1,n_pt_bins)
            read(11,*) (href1_phi_2(i),i=1,n_phi_bins)
            read(11,*) (href1_eta_2(i),i=1,n_eta_bins)
         end if

      end if    ! End one-body reference histogram input

      Close(unit=11)

C-----------------------------
      ELSE IF(mode.eq.4) THEN
C-----------------------------

CCC   Read the next event from 'event_text.in' and load accepted tracks
C     into the 'trk' data structure using the flag information about each
C     line type in the file 'event_tracks.select'.
C
C     For this mode to run successfully the calling program must:
C       (1) initially set the event_line_counter = 0
C       (2) open the 'event_text.in' and 'event_tracks.select' files
C           as units 2 and 4, respectively.
C       (3) Close units 2 and 4 when finished.

CCC   Initialize accepted track counters for this new event:

      track_counter   = 0     ! Counts all accepted tracks
      track_counter_1 = 0     ! Counts all accepted tracks of type pid(1)
      track_counter_2 = 0     ! Counts all accepted tracks of type pid(2)

      If(ALICE .eq. 1) Then
         Call AliHbtp_GetNumberTracks(ntracks)
         do i = 1,ntracks
            Call AliHbtp_GetTrack(i,flag,px,py,pz,ge_pid)
            eg_pid = ge_pid
            if(flag.eq.flag4) then
         track_counter = track_counter + 1

         if(eg_pid.eq.pid(1) .and. pid(1).gt.0) then
            track_counter_1 = track_counter_1 + 1
         end if

         if(eg_pid.eq.pid(2) .and. pid(2).gt.0) then
            track_counter_2 = track_counter_2 + 1
         end if

         mass = part_mass(eg_pid)
         Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
         trk_ge_pid(track_counter)        = eg_pid
         trk_px(track_counter)            = px
         trk_py(track_counter)            = py
         trk_pz(track_counter)            = pz
         trk_id(track_counter)            = track_counter
         trk_start_vertex(track_counter)  = 0
         trk_stop_vertex(track_counter)   = 0
         trk_event_line(track_counter)    = 0
         trk_flag(track_counter)          = 0
         trk_px_sec(track_counter)        = 0
         trk_py_sec(track_counter)        = 0
         trk_pz_sec(track_counter)        = 0
         trk_sector(track_counter)        = 0
         trk_out_flag(track_counter)      = 0
         trk_merge_flag(track_counter)    = 0
         trk_E(track_counter)             = E
         trk_pt(track_counter)            = pt
         trk_phi(track_counter)           = phi
         trk_eta(track_counter)           = eta
        end if
       end do
         n_part_1_trk   = track_counter_1
         n_part_2_trk   = track_counter_2
         n_part_tot_trk = track_counter

      Else If(ALICE .eq. 0) Then

80    read(4,11,err=81,end=82) flag
      event_line_counter = event_line_counter + 1
      
      if(flag .ne. 4) then
         read(2,10,err=83,end=82) dummy
      else if(flag .eq. 4) then
         read(2,41) ge_pid,px,py,pz,tid,start_v,stop_v,eg_pid

CCC   Check if the 'event_text.in' file includes non-zero PID
CCC   values for the variable 'eg_pid'.  If this is zero, then
CCC   use the ge_pid value.
         if(eg_pid.eq.0 .and. ge_pid.ne.0) eg_pid = ge_pid

         track_counter = track_counter + 1
      
         if(eg_pid.eq.pid(1) .and. pid(1).gt.0) then
            track_counter_1 = track_counter_1 + 1
         end if

         if(eg_pid.eq.pid(2) .and. pid(2).gt.0) then
            track_counter_2 = track_counter_2 + 1
         end if

         mass = part_mass(eg_pid)
         Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
         trk_ge_pid(track_counter)        = eg_pid
         trk_px(track_counter)            = px
         trk_py(track_counter)            = py
         trk_pz(track_counter)            = pz
         trk_id(track_counter)            = track_counter
         trk_start_vertex(track_counter)  = start_v
         trk_stop_vertex(track_counter)   = stop_v
         trk_event_line(track_counter)    = event_line_counter
         trk_flag(track_counter)          = 0
         trk_px_sec(track_counter)        = 0
         trk_py_sec(track_counter)        = 0
         trk_pz_sec(track_counter)        = 0
         trk_sector(track_counter)        = 0
         trk_out_flag(track_counter)      = 0
         trk_merge_flag(track_counter)    = 0
         trk_E(track_counter)             = E
         trk_pt(track_counter)            = pt
         trk_phi(track_counter)           = phi
         trk_eta(track_counter)           = eta
      end if

      if(flag.ne.1) then
         go to 80   ! Return to S.N. 80 and read next record in file
      else if(flag.eq.1) then
         n_part_1_trk   = track_counter_1
         n_part_2_trk   = track_counter_2
         n_part_tot_trk = track_counter
      end if

82    Return
81    write(7,84)
84    Format(5x,'Read error from file event_tracks.select for mode=4',
     1       ' - STOP')
      Stop
83    write(7,85)
85    Format(5x,'Read error from file event_text.in for mode=4',
     1       ' - STOP')
      Stop
      End If  !  ALICE Data I/O Option

C-----------------------------
      ELSE IF(mode.eq.5) THEN
C-----------------------------

CCC   Read the next event from 'event_text.in' and load accepted tracks
C     into the 'trk2' data structure using the flag information about each
C     line type in the file 'event_tracks.select'.
C
C     For this mode to run successfully the calling program must:
C       (1) initially set the event_line_counter = 0
C       (2) open the 'event_text.in' and 'event_tracks.select' files
C           as units 2 and 4, respectively.
C       (3) Close units 2 and 4 when finished.

CCC   Initialize accepted track counters for this new event:

      track_counter   = 0     ! Counts all accepted tracks
      track_counter_1 = 0     ! Counts all accepted tracks of type pid(1)
      track_counter_2 = 0     ! Counts all accepted tracks of type pid(2)

      If(ALICE .eq. 1) Then
         Call AliHbtp_GetNumberTracks(ntracks)
         do i = 1,ntracks
            Call AliHbtp_GetTrack(i,flag,px,py,pz,ge_pid)
            eg_pid = ge_pid
            if(flag.eq.flag4) then
         track_counter = track_counter + 1

         if(eg_pid.eq.pid(1) .and. pid(1).gt.0) then
            track_counter_1 = track_counter_1 + 1
         end if

         if(eg_pid.eq.pid(2) .and. pid(2).gt.0) then
            track_counter_2 = track_counter_2 + 1
         end if

         mass = part_mass(eg_pid)
         Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
         trk2_ge_pid(track_counter)        = eg_pid
         trk2_px(track_counter)            = px
         trk2_py(track_counter)            = py
         trk2_pz(track_counter)            = pz
         trk2_id(track_counter)            = track_counter
         trk2_start_vertex(track_counter)  = 0
         trk2_stop_vertex(track_counter)   = 0
         trk2_event_line(track_counter)    = 0
         trk2_flag(track_counter)          = 0
         trk2_px_sec(track_counter)        = 0
         trk2_py_sec(track_counter)        = 0
         trk2_pz_sec(track_counter)        = 0
         trk2_sector(track_counter)        = 0
         trk2_out_flag(track_counter)      = 0
         trk2_merge_flag(track_counter)    = 0
         trk2_E(track_counter)             = E
         trk2_pt(track_counter)            = pt
         trk2_phi(track_counter)           = phi
         trk2_eta(track_counter)           = eta
        end if
       end do
         n_part_1_trk2   = track_counter_1
         n_part_2_trk2   = track_counter_2
         n_part_tot_trk2 = track_counter

      Else If(ALICE.eq.0) Then

90    read(4,11,err=91,end=92) flag
      event_line_counter = event_line_counter + 1
      
      if(flag .ne. 4) then
         read(2,10,err=93,end=92) dummy
      else if(flag .eq. 4) then
         read(2,41) ge_pid,px,py,pz,tid,start_v,stop_v,eg_pid

CCC   Check if the 'event_text.in' file includes non-zero PID
CCC   values for the variable 'eg_pid'.  If this is zero, then
CCC   use the ge_pid value.
         if(eg_pid.eq.0 .and. ge_pid.ne.0) eg_pid = ge_pid

         track_counter = track_counter + 1
      
         if(eg_pid.eq.pid(1) .and. pid(1).gt.0) then
            track_counter_1 = track_counter_1 + 1
         end if

         if(eg_pid.eq.pid(2) .and. pid(2).gt.0) then
            track_counter_2 = track_counter_2 + 1
         end if

         mass = part_mass(eg_pid)
         Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
         trk2_ge_pid(track_counter)        = eg_pid
         trk2_px(track_counter)            = px
         trk2_py(track_counter)            = py
         trk2_pz(track_counter)            = pz
         trk2_id(track_counter)            = track_counter
         trk2_start_vertex(track_counter)  = start_v
         trk2_stop_vertex(track_counter)   = stop_v
         trk2_event_line(track_counter)    = event_line_counter
         trk2_flag(track_counter)          = 0
         trk2_px_sec(track_counter)        = 0
         trk2_py_sec(track_counter)        = 0
         trk2_pz_sec(track_counter)        = 0
         trk2_sector(track_counter)        = 0
         trk2_out_flag(track_counter)      = 0
         trk2_merge_flag(track_counter)    = 0
         trk2_E(track_counter)             = E
         trk2_pt(track_counter)            = pt
         trk2_phi(track_counter)           = phi
         trk2_eta(track_counter)           = eta
      end if

      if(flag.ne.1) then
         go to 90   ! Return to S.N. 90 and read next record in file
      else if(flag.eq.1) then
         n_part_1_trk2   = track_counter_1
         n_part_2_trk2   = track_counter_2
         n_part_tot_trk2 = track_counter
      end if

92    Return
91    write(7,94)
94    Format(5x,'Read error from file event_tracks.select for mode=5',
     1       ' - STOP')
      Stop
93    write(7,95)
95    Format(5x,'Read error from file event_text.in for mode=5',
     1       ' - STOP')
      Stop

      End If   !  ALICE Data I/O Option

C-----------------------------
      ELSE IF(mode.eq.6) THEN
C-----------------------------

CCC   Read finite source size Coulomb pair correlation corrections and
CCC   interpolate to requested source radius and save the results for q,
CCC   like and unlike pairs in common/coulomb/.

      if(switch_coulomb.eq.3 .and. Q0.ge.coulradmin .and.
     1   Q0.le.coulradmax) then

CCC   Initially, read and interpolate like pair Coulomb corrections:
C     ALICE

      If(ALICE .eq. 1) then
      
        CALL GETENVF('ALICE_ROOT',CHROOT)
        LNROOT = LNBLNK(CHROOT)
        
        IF(LNROOT.LE.0) THEN 
         PRINT*,'**********************************'
         PRINT*,'*        HBT PROCESSOR           *'
         PRINT*,'*        -----------             *'
         PRINT*,'*   DATA File  not found         *'
         PRINT*,'*         Program STOP           *'
         PRINT*,'*   Check ALICE_ROOT environment *'
         PRINT*,'*           variable             *'
         PRINT*,'**********************************'        
         errorcode = 1
         return
        ENDIF

        FILNAM=CHROOT(1:LNROOT)//'/data/cpp_06.dat' 
	open(unit=21,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpp_08.dat' 
	open(unit=22,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpp_10.dat' 
	open(unit=23,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpp_12.dat' 
	open(unit=24,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpp_14.dat' 
	open(unit=25,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpp_16.dat' 
	open(unit=26,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpp_18.dat' 
	open(unit=27,status='old',access='sequential',
     1       file=FILNAM)
	
      ELSE
        open(unit=21,status='old',access='sequential',
     1       file='cpp_06.dat')
        open(unit=22,status='old',access='sequential',
     1       file='cpp_08.dat')
        open(unit=23,status='old',access='sequential',
     1       file='cpp_10.dat')
        open(unit=24,status='old',access='sequential',
     1       file='cpp_12.dat')
        open(unit=25,status='old',access='sequential',
     1       file='cpp_14.dat')
        open(unit=26,status='old',access='sequential',
     1       file='cpp_16.dat')
        open(unit=27,status='old',access='sequential',
     1       file='cpp_18.dat')
      ENDIF
      

      do i = 1,max_c2_coul
         do j = 1,ncoulradsteps
            read(20+j,*) q_coul(i), function(j)
         end do
         Call AliHbtp_interp(coulradmin,coulradmax,coulradstep,
     1        ncoulradsteps,function,20,Q0,c2_coul_like(i))
      end do

      close(unit=21)
      close(unit=22)
      close(unit=23)
      close(unit=24)
      close(unit=25)
      close(unit=26)
      close(unit=27)

CCC   Next read and interpolate the unlike pair Coulomb corrections:

      If(ALICE .eq. 1) then
        FILNAM=CHROOT(1:LNROOT)//'/data/cpm_06.dat' 
	open(unit=31,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpm_08.dat' 
	open(unit=32,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpm_10.dat' 
	open(unit=33,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpm_12.dat' 
	open(unit=34,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpm_14.dat' 
	open(unit=35,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpm_16.dat' 
	open(unit=36,status='old',access='sequential',
     1       file=FILNAM)

        FILNAM=CHROOT(1:LNROOT)//'/data/cpm_18.dat' 
	open(unit=37,status='old',access='sequential',
     1       file=FILNAM)
    
      else
        open(unit=31,status='old',access='sequential',
     1       file='cpm_06.dat')
        open(unit=32,status='old',access='sequential',
     1       file='cpm_08.dat')
        open(unit=33,status='old',access='sequential',
     1       file='cpm_10.dat')
        open(unit=34,status='old',access='sequential',
     1       file='cpm_12.dat')
        open(unit=35,status='old',access='sequential',
     1       file='cpm_14.dat')
        open(unit=36,status='old',access='sequential',
     1       file='cpm_16.dat')
        open(unit=37,status='old',access='sequential',
     1       file='cpm_18.dat')
      EndIf
      
      do i = 1,max_c2_coul
         do j = 1,ncoulradsteps
            read(30+j,*) q_coul(i), function(j)
         end do
         Call AliHbtp_interp(coulradmin,coulradmax,coulradstep, 
     1        ncoulradsteps,function,20,Q0,c2_coul_unlike(i))
      end do

      close(unit=31)
      close(unit=32)
      close(unit=33)
      close(unit=34)
      close(unit=35)
      close(unit=36)
      close(unit=37)

CCC   Convert the input q values which are in MeV/c, to GeV/c:

      do i = 1,max_c2_coul
         q_coul(i) = 0.001*q_coul(i)
      end do

      end if   


C----------------------------
      ELSE IF(mode.eq.7) THEN
C----------------------------

CCC   Read next event from 'event_text.in', load accepted tracks into 'trk'
CCC   data structure using the flag information in the file
CCC   'event_tracks.select', copy contents of 'event_text.in' and
CCC   'event_tracks.select', for this one event only, into temporary files
CCC   'event_text_aux.in' and 'event_tracks_aux.select', respectively.
C
C     For this mode to run successfully the calling program must:
C       (1) initially set the event_line_counter = 0
C       (2) open the 'event_text.in' and 'event_tracks.select' files
C           as units 2 and 4, respectively.
C       (3) Close units 2 and 4 when finished.

CCC   Initialize accepted track counters for this new event:

      track_counter   = 0     ! Counts all accepted tracks
      track_counter_1 = 0     ! Counts all accepted tracks of type pid(1)
      track_counter_2 = 0     ! Counts all accepted tracks of type pid(2)

      If(ALICE .eq. 1) Then
         Call AliHbtp_GetNumberTracks(ntracks)
         do i = 1,ntracks
            Call AliHbtp_GetTrack(i,flag,px,py,pz,ge_pid)
            eg_pid = ge_pid
            if(flag.eq.flag4) then
         track_counter = track_counter + 1

         if(eg_pid.eq.pid(1) .and. pid(1).gt.0) then
            track_counter_1 = track_counter_1 + 1
         end if

         if(eg_pid.eq.pid(2) .and. pid(2).gt.0) then
            track_counter_2 = track_counter_2 + 1
         end if

         mass = part_mass(eg_pid)
         Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
         trk_ge_pid(track_counter)        = eg_pid
         trk_px(track_counter)            = px
         trk_py(track_counter)            = py
         trk_pz(track_counter)            = pz
         trk_id(track_counter)            = track_counter
         trk_start_vertex(track_counter)  = 0
         trk_stop_vertex(track_counter)   = 0
         trk_event_line(track_counter)    = 0
         trk_flag(track_counter)          = 0
         trk_px_sec(track_counter)        = 0
         trk_py_sec(track_counter)        = 0
         trk_pz_sec(track_counter)        = 0
         trk_sector(track_counter)        = 0
         trk_out_flag(track_counter)      = 0
         trk_merge_flag(track_counter)    = 0
         trk_E(track_counter)             = E
         trk_pt(track_counter)            = pt
         trk_phi(track_counter)           = phi
         trk_eta(track_counter)           = eta
        end if
       end do
         n_part_1_trk   = track_counter_1
         n_part_2_trk   = track_counter_2
         n_part_tot_trk = track_counter

      Else If(ALICE .eq. 0) Then

CCC   Open temporary files:

      open(unit=12,status='unknown',access='sequential',
     1     file='event_text_aux.in')
      open(unit=14,status='unknown',access='sequential',
     1     file='event_tracks_aux.select')

100   read(4,11,err=101,end=102) flag
      event_line_counter = event_line_counter + 1
      write(14,11) flag
      if(flag.eq.1) then
         read(2,10,err=103,end=102) comment_event_label
         write(12,10) comment_event_label
      else if(flag .eq. 2) then
         read(2,10,err=103,end=102) vertex_label
         write(12,10) vertex_label
      else if(flag .eq. 3) then
         read(2,10,err=103,end=102) comment_event_label
         write(12,10) comment_event_label
      else if(flag .eq. 5) then
         read(2,10,err=103,end=102) gener_label
         write(12,10) gener_label
      else if(flag .eq. 4) then
         read(2,41) ge_pid,px,py,pz,tid,start_v,stop_v,eg_pid
         write(12,41) ge_pid,px,py,pz,tid,start_v,stop_v,eg_pid

CCC   Check if the 'event_text.in' file includes non-zero PID
CCC   values for the variable 'eg_pid'.  If this is zero, then
CCC   use the ge_pid value.
         if(eg_pid.eq.0 .and. ge_pid.ne.0) eg_pid = ge_pid

         track_counter = track_counter + 1

         if(eg_pid.eq.pid(1) .and. pid(1).gt.0) then
            track_counter_1 = track_counter_1 + 1
         end if

         if(eg_pid.eq.pid(2) .and. pid(2).gt.0) then
            track_counter_2 = track_counter_2 + 1
         end if

         mass = part_mass(eg_pid)
         Call Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,1)
         trk_ge_pid(track_counter)        = eg_pid
         trk_px(track_counter)            = px
         trk_py(track_counter)            = py
         trk_pz(track_counter)            = pz
         trk_id(track_counter)            = track_counter
         trk_start_vertex(track_counter)  = start_v
         trk_stop_vertex(track_counter)   = stop_v
         trk_event_line(track_counter)    = event_line_counter
         trk_flag(track_counter)          = 0
         trk_px_sec(track_counter)        = 0
         trk_py_sec(track_counter)        = 0
         trk_pz_sec(track_counter)        = 0
         trk_sector(track_counter)        = 0
         trk_out_flag(track_counter)      = 0
         trk_merge_flag(track_counter)    = 0
         trk_E(track_counter)             = E
         trk_pt(track_counter)            = pt
         trk_phi(track_counter)           = phi
         trk_eta(track_counter)           = eta
      else
         read(2,10,err=103,end=102) comment_event_label
         write(12,10) comment_event_label
      end if

      if(flag.ne.1) then
         go to 100   ! Return to S.N. 100 and read next record in file
      else if(flag.eq.1) then
         n_part_1_trk   = track_counter_1
         n_part_2_trk   = track_counter_2
         n_part_tot_trk = track_counter
      end if

102   Close(unit=12)
      Close(unit=14)
      Return
101   write(7,104)
104   Format(5x,'Read error from file event_tracks.select for mode=7',
     1       ' - STOP')
      Stop
103   write(7,105)
105   Format(5x,'Read error from file event_text.in for mode=7',
     1       ' - STOP')
      Stop

      End If  !  ALICE Data I/O Option

C----------------------------
      ELSE IF(mode.eq.8) THEN
C----------------------------

CCC   Read contents of 'event_text_aux.in' using the flag values in
CCC   tmp. file 'event_tracks_aux.select' and copy this into the final
CCC   output event file, 'event_hbt_text.out', where the momentum values
CCC   of the accepted tracks in the initial input event file are replaced
CCC   with the adjusted/correlated values obtained from the 'trk' table.
C
C     For this to work successfully the calling program must:
C      (1) initially set the event_line_counter = 0
C      (2) open the 'event_hbt_text.out' file as unit = 10
C      (3) Close unit 10 when finished

CCC   Initialize accepted track counters:

      track_counter = 0

      If(ALICE .eq. 1) Then
         Call AliHbtp_GetNumberTracks(ntracks)
         do i = 1,ntracks
            Call AliHbtp_GetTrack(i,flag,px,py,pz,ge_pid)
            if(flag.eq.flag4) then
              track_counter = track_counter + 1
              if(trk_accep .ge. 1.000  .or. (trk_accep .lt. 1.00
     1           .and. hbtpran(irand) .le. trk_accep)) then
C                 write(*,*) '  FFF: 3 calling PutTrack i = ',i
                 Call AliHbtp_PutTrack(i,flag,
     1                trk_px(track_counter),
     2                trk_py(track_counter),
     3                trk_pz(track_counter),
     4                ge_pid)
              end if
            end if
         end do

      Else If(ALICE .eq. 0) Then

CCC   Open temporary, auxiliary files:

      open(unit=12,status='old',access='sequential',
     1     file='event_text_aux.in')
      open(unit=14,status='old',access='sequential',
     1     file='event_tracks_aux.select')

120   read(14,11,err=121,end=122) flag
      file10_line_counter = file10_line_counter + 1
      if(flag.eq.1) then
         read(12,10,err=123,end=122) comment_event_label
         write(10,10) comment_event_label
      else if(flag .eq. 2) then
         read(12,10,err=123,end=122) vertex_label
         write(10,10) vertex_label
      else if(flag .eq. 3) then
         read(12,10,err=123,end=122) comment_event_label
         write(10,10) comment_event_label
      else if(flag .eq. 5) then
         read(12,10,err=123,end=122) gener_label
         write(10,10) gener_label
      else if(flag .eq. 4) then
         read(12,41,err=123,end=122)
     1      ge_pid,px,py,pz,tid,start_v,stop_v,eg_pid
         track_counter = track_counter + 1
         if(tid.eq.0) tid = trk_id(track_counter)
         if(trk_event_line(track_counter).eq.file10_line_counter)then
            if(trk_accep .ge. 1.000  .or. (trk_accep .lt. 1.00
     1         .and. hbtpran(irand) .le. trk_accep)) then
            write(10,841)ge_pid                   ,
     1                   trk_px(track_counter)    ,
     2                   trk_py(track_counter)    ,
     3                   trk_pz(track_counter)    ,
     4                   tid                      ,
     5                   start_v                  ,
     6                   stop_v                   ,
     7                   trk_ge_pid(track_counter)
841         Format('TRACK:',1x,I6,3(1x,G12.5),4(1x,I6))
            end if
         else
            write(7,127)
            write(7,126) track_counter, trk_event_line(track_counter),
     1                file10_line_counter
127         Format(5x,'Track table rows and Event file line count ',
     1          'out-of-synch. - STOP')
126         Format(5x,'track_counter, trk().event_line,',
     1          'file10_line_counter = ',3I10)
            Stop
         end if
      else
         read(12,10,err=123,end=122) comment_event_label
         write(10,10) comment_event_label
      end if

      if(flag .ne. 1) go to 120  !  Return to S.N. 120 and read next record

122   Close(unit=12,status='delete')
      Close(unit=14,status='delete')
      Return
121   write(7,124)
124   Format(5x,'Read error from file event_tracks_aux.select',
     1       ' for mode = 8 - STOP')
      Stop
123   write(7,125)
125   Format(5x,'Read error from file event_text_aux.in',
     1       ' for mode = 8 - STOP')
      Stop
     
      End If  !  ALICE Data I/O Option

C-----------------
      END IF     ! End of read_data mode selection options
C-----------------

      Return
      END

C-----------------------------------------------------------------------


      subroutine getref_hist
      implicit none

CCC   This subroutine controls the reading or calculation and output
CCC   of the several reference histograms.  These include:
CCC      (a) the one-body {pt,phi,eta} 1D distributions for 1 or 2
CCC          particle ID types.
CCC      (b) the two-body pair-wise histograms for like and unlike
CCC          pairs; for 1D and/or 3D fine mesh and 3D coarse mesh
CCC          distributions.

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'
      Include 'common_particle.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,ipt,iphi,ieta,sign_toggle

      if(ref_control .eq. 1) then

CCC   read pair and one-body reference histograms:
         Call read_data(3)
      else if(ref_control .eq. 2) then

CCC   calculate the pair and one-body histograms:
CCC   Open event and flag files:

      If(ALICE .eq. 0) Then
         open(unit=2,status='old',access='sequential',
     1        file='event_text.in')
         open(unit=4,status='old',access='sequential',
     1        file='event_tracks.select')
      End If

CCC   Initialize counters:

         n_part_used_1_ref    = 0
         n_part_used_2_ref    = 0
         num_pairs_like_ref   = 0
         num_pairs_unlike_ref = 0
         event_line_counter   = 0

CCC   Read event header lines (no tracks are in this part):
      If(ALICE .eq. 0) Then
         Call read_data(4)
      End If

CCC   Set toggle switch to alternate between loading event tracks into
CCC   table 'trk' and table 'trk2':
         sign_toggle = 1

CCC   Start Event Loop:
C         write(*,*) 'REF HISTO N Ev = ', n_events
         do i = 1,n_events
           If(ALICE .eq. 1) Then
              Call AliHbtp_SetActiveEventNumber(i)
           End If
           if(sign_toggle .eq. 1) then   ! Put tracks into 'trk'
             Call read_data(4)
             Call tindex(1,0)
             Call stm_build(1,0,0)
             if(pid(1) .gt. 0) then
               Call histog1(1,0,1,pid(1),0.0,0.0,0.0)
               n_part_used_1_ref = n_part_used_1_ref + n_part_used_1_trk

               do ipt = 1,n_pt_bins
               href1_pt_1(ipt) = href1_pt_1(ipt) + hist1_pt_1(ipt)
               if (href1_pt_1(ipt) .lt. 0 ) then
                 write(*,*) 'href1_pt_1 bin ',ipt,'is less then 0'
               endif 
               end do

               do iphi = 1,n_phi_bins
               href1_phi_1(iphi) = href1_phi_1(iphi) + hist1_phi_1(iphi)
               if (href1_phi_1(iphi) .lt. 0 ) then
                 write(*,*) 'href1_phi_1 bin ',iphi,'is less then 0'
               endif 
               end do

               do ieta = 1,n_eta_bins
               href1_eta_1(ieta) = href1_eta_1(ieta) + hist1_eta_1(ieta)
               if (href1_eta_1(ieta) .lt. 0 ) then
                 write(*,*) 'href1_eta_1 bin ',ieta,'is less then 0'
               endif 
               end do
             end if

             if(pid(2) .gt. 0) then
               Call histog1(1,0,2,pid(2),0.0,0.0,0.0)
               n_part_used_2_ref = n_part_used_2_ref + n_part_used_2_trk

               do ipt = 1,n_pt_bins
               href1_pt_2(ipt) = href1_pt_2(ipt) + hist1_pt_2(ipt)
               end do

               do iphi = 1,n_phi_bins
               href1_phi_2(iphi) = href1_phi_2(iphi) + hist1_phi_2(iphi)
               end do

               do ieta = 1,n_eta_bins
               href1_eta_2(ieta) = href1_eta_2(ieta) + hist1_eta_2(ieta)
               end do
             end if

           else if(sign_toggle .eq. (-1)) then  ! Put tracks into 'trk2'
             Call read_data(5)
             Call tindex(2,0)
             Call stm_build(2,0,0)
             if(pid(1) .gt. 0) then
               Call histog1(4,0,1,pid(1),0.0,0.0,0.0)
               n_part_used_1_ref = n_part_used_1_ref +n_part_used_1_trk2

               do ipt = 1,n_pt_bins
               href1_pt_1(ipt) = href1_pt_1(ipt) + hist1_pt_1(ipt)
               end do

               do iphi = 1,n_phi_bins
               href1_phi_1(iphi) = href1_phi_1(iphi) + hist1_phi_1(iphi)
               end do

               do ieta = 1,n_eta_bins
               href1_eta_1(ieta) = href1_eta_1(ieta) + hist1_eta_1(ieta)
               end do
             end if

             if(pid(2) .gt. 0) then
               Call histog1(4,0,2,pid(2),0.0,0.0,0.0)
               n_part_used_2_ref = n_part_used_2_ref +n_part_used_2_trk2

               do ipt = 1,n_pt_bins
               href1_pt_2(ipt) = href1_pt_2(ipt) + hist1_pt_2(ipt)
               end do

               do iphi = 1,n_phi_bins
               href1_phi_2(iphi) = href1_phi_2(iphi) + hist1_phi_2(iphi)
               end do

               do ieta = 1,n_eta_bins
               href1_eta_2(ieta) = href1_eta_2(ieta) + hist1_eta_2(ieta)
               end do
             end if

           end if   !  End read and load to trk or trk2 option

           sign_toggle = -sign_toggle

           if(i .gt. 1) then    !  Compute 2-body reference histograms
              Call histog2(4,0,0,0,0,0.0,0.0,0.0,0.0)
              num_pairs_like_ref = num_pairs_like_ref
     1           + n_part_used_1_trk * n_part_used_1_trk2
     2           + n_part_used_2_trk * n_part_used_2_trk2
              num_pairs_unlike_ref = num_pairs_unlike_ref
     1           + n_part_used_1_trk * n_part_used_2_trk2
     2           + n_part_used_2_trk * n_part_used_1_trk2
C              write(*,*) 'num_pairs_like_ref',num_pairs_like_ref
C              write(*,*) 'num_pairs_unlike_ref',num_pairs_unlike_ref
           end if

        end do   !  End of Event Loop

CCC   Write out the pair and one-body reference Histograms:
      Call write_data(2,0)

      If(ALICE .eq. 0) Then
      Close(unit=2)
      Close(unit=4)
      End If

      end if  !  End Reference Histogram read/calculate option

      Return
      END

C----------------------------------------------------------------------


      subroutine AliHbtp_interp(rmin,rmax,rstep,nrsteps,function,
     1                          ndim,r,answer)
      implicit none

CCC   This routine interpolates the function values and puts the result
CCC   into 'answer'.  It uses 2,3 or 4 mesh points which must be equally
CCC   spaced.  The method uses the Lagrange interpolation formulas given
CCC   in Abramowitz and Stegun, ``Handbook of Mathematical Functions,''
CCC   (Dover Publications, New York, 1970); pages 878-879.

CCC   Definition of Variables in the Argument List:

CCC   rmin     = lower limit of independent variable for input function
CCC   rmax     = upper limit of independent variable for input function
CCC   rstep    = step size of independent variable
CCC   nrsteps  = (redundant) # of input steps
CCC   function(ndim) = Array of function values to be interpolated
CCC   ndim     = array dimension size in calling program
CCC   r        = coordinate value of independent variable where interpolation
CCC              is needed.
CCC   answer   = interpolated value

CCC   The algorithm will use the maximum number of points in the
CCC   interpolation, up to a maximum of 4

CCC   If the requested coordinate value, r, is out-of-range, then
CCC   'answer' is returned with a 0.0 value.

CCC   Local Variable Type Declarations:

      integer*4 ndim, nrsteps, ik

      real*4 rmin,rmax,rstep,r,answer,rshift,p
      real*4 function(ndim),w1,w2,w3,w4

CCC   Check Mesh:

      if(abs(((rmax-rmin)/float(nrsteps-1))-rstep).gt.0.00001) then
         write(7,10) rmin,rmax,rstep,nrsteps
10       Format(5x,'Interp mesh inconsistent:',3E12.5,I5,
     1          ' - STOP')
      Return
      end if

CCC   Check range:
     
      if(r .lt. rmin .or. r .gt. rmax) then
         write(7,11) rmin,rmax,r
11       Format(5x,'Interp called with r out-of-range =',3E12.5)
         answer = 0.0
         Return
      end if

CCC   Begin interpolation:

      if(nrsteps .eq. 2) then
         p = (r - rmin)/rstep
         answer = (1.0 - p)*function(1) + p*function(2)
      else if(nrsteps .eq. 3) then
         p = (r - (rmin + rstep))/rstep
         answer = 0.5*p*(p-1.0)*function(1) + (1.0 - p*p)
     1            *function(2) + 0.5*p*(p+1.0)*function(3)
      else if(nrsteps .ge. 4) then
         rshift = r - rmin

         if(rshift .le. rstep) then
            ik = 2
            p = (rshift - rstep)/rstep
         else if(rshift .ge. (rmax - rstep - rmin)) then
            ik = nrsteps - 2
            p = (rshift - (rmax - rmin - 2.0*rstep))/rstep
         else
            ik = int(rshift/rstep + 1.000001)
            if(ik .le. 1) ik = 2
            if(ik .ge. (nrsteps-1)) ik = nrsteps - 2
            p = (rshift - float(ik-1)*rstep)/rstep
         end if

         w1 = -p*(p-1.0)*(p-2.0)/6.0
         w2 = (p*p-1.0)*(p-2.0)/2.0
         w3 = -p*(p+1.0)*(p-2.0)/2.0
         w4 = p*(p*p-1.0)/6.0

         answer = w1*function(ik-1) + w2*function(ik)
     1          + w3*function(ik+1) + w4*function(ik+2)
      end if  !  End # interplation points option

      Return
      END

C--------------------------------------------------------------------


      subroutine Hbtp_particle_prop
      implicit none

CCC   Fill particle properties table /particle/ with Geant 3 particle ID
CCC   numbers, charge (in units of |e|), mass in GeV/c**2 and lifetime
CCC   in seconds.  See the Geant 3.15 Manual User's Guide, pages: CONS
CCC   300-1 and -2.

      Include 'common_particle.inc'

CCC   Local Variable Type Declarations:

      integer*4 i

      do i = 1,part_maxlen
      part_id(i) = i
      end do

CCC   Set Particle Masses:

      part_mass( 1) = 0.0                   ! Gamma
      part_mass( 2) = 0.00051099906         ! Positron
      part_mass( 3) = 0.00051099906         ! Electron
      part_mass( 4) = 0.0                   ! Neutrino
      part_mass( 5) = 0.105658389           ! Muon+
      part_mass( 6) = 0.105658389           ! Muon-
      part_mass( 7) = 0.1349743             ! Pion0
      part_mass( 8) = 0.1395679             ! Pion+
      part_mass( 9) = 0.1395679             ! Pion-
      part_mass(10) = 0.497671              ! Kaon 0 long
      part_mass(11) = 0.493646              ! Kaon+
      part_mass(12) = 0.493646              ! Kaon-
      part_mass(13) = 0.93956563            ! Neutron
      part_mass(14) = 0.93827231            ! Proton
      part_mass(15) = 0.93827231            ! Antiproton
      part_mass(16) = 0.497671              ! Kaon 0 short
      part_mass(17) = 0.54745               ! Eta
      part_mass(18) = 1.11563               ! Lambda
      part_mass(19) = 1.18937               ! Sigma+
      part_mass(20) = 1.19255               ! Sigma0
      part_mass(21) = 1.197465              ! Sigma-
      part_mass(22) = 1.31485               ! Xi 0
      part_mass(23) = 1.32133               ! Xi -
      part_mass(24) = 1.67243               ! Omega
      part_mass(25) = 0.93956563            ! Antineutron
      part_mass(26) = 1.11563               ! Antilambda
      part_mass(27) = 1.18937               ! Anti-Sigma -
      part_mass(28) = 1.19255               ! Anti-Sigma 0
      part_mass(29) = 1.197465              ! Anti-Sigma +
      part_mass(30) = 1.31485               ! AntiXi 0
      part_mass(31) = 1.32133               ! AntiXi +
      part_mass(32) = 1.67243               ! Anti-Omega +
      part_mass(33) = 0.0
      part_mass(34) = 0.0
      part_mass(35) = 0.0
      part_mass(36) = 0.0
      part_mass(37) = 0.0
      part_mass(38) = 0.0
      part_mass(39) = 0.0
      part_mass(40) = 0.0
      part_mass(41) = 0.0
      part_mass(42) = 0.0
      part_mass(43) = 0.0
      part_mass(44) = 0.0
      part_mass(45) = 1.875613              ! Deuteron
      part_mass(46) = 2.80925               ! Triton
      part_mass(47) = 3.727417              ! Alpha
      part_mass(48) = 0.0                   ! Geantino (Fake particle)
      part_mass(49) = 2.80923               ! He3
      part_mass(50) = 0.0                   ! Cerenkov (Fake particle)

CCC   Set Particle Charge:

      part_charge( 1) =  0      ! Gamma
      part_charge( 2) =  1      ! Positron
      part_charge( 3) = -1      ! Electron
      part_charge( 4) =  0      ! Neutrino
      part_charge( 5) =  1      ! Muon+
      part_charge( 6) = -1      ! Muon-
      part_charge( 7) =  0      ! Pion0
      part_charge( 8) =  1      ! Pion+
      part_charge( 9) = -1      ! Pion-
      part_charge(10) =  0      ! Kaon 0 long
      part_charge(11) =  1      ! Kaon+
      part_charge(12) = -1      ! Kaon-
      part_charge(13) =  0      ! Neutron
      part_charge(14) =  1      ! Proton
      part_charge(15) = -1      ! Antiproton
      part_charge(16) =  0      ! Kaon 0 short
      part_charge(17) =  0      ! Eta
      part_charge(18) =  0      ! Lambda
      part_charge(19) =  1      ! Sigma+
      part_charge(20) =  0      ! Sigma0
      part_charge(21) = -1      ! Sigma-
      part_charge(22) =  0      ! Xi 0
      part_charge(23) = -1      ! Xi -
      part_charge(24) = -1      ! Omega
      part_charge(25) =  0      ! Antineutron
      part_charge(26) =  0      ! Antilambda
      part_charge(27) = -1      ! Anti-Sigma -
      part_charge(28) =  0      ! Anti-Sigma 0
      part_charge(29) =  1      ! Anti-Sigma +
      part_charge(30) =  0      ! AntiXi 0
      part_charge(31) =  1      ! AntiXi +
      part_charge(32) =  1      ! Anti-Omega +
      part_charge(33) =  0
      part_charge(34) =  0
      part_charge(35) =  0
      part_charge(36) =  0
      part_charge(37) =  0
      part_charge(38) =  0
      part_charge(39) =  0
      part_charge(40) =  0
      part_charge(41) =  0
      part_charge(42) =  0
      part_charge(43) =  0
      part_charge(44) =  0
      part_charge(45) =  1      ! Deuteron
      part_charge(46) =  1      ! Triton
      part_charge(47) =  2      ! Alpha
      part_charge(48) =  0      ! Geantino (Fake particle)
      part_charge(49) =  2      ! He3
      part_charge(50) =  0      ! Cerenkov (Fake particle)

CCC   Set Particle Lifetimes:

      part_lifetime( 1) = 1.0E+30       ! Gamma
      part_lifetime( 2) = 1.0E+30       ! Positron
      part_lifetime( 3) = 1.0E+30       ! Electron
      part_lifetime( 4) = 1.0E+30       ! Neutrino
      part_lifetime( 5) = 2.19703E-06   ! Muon+
      part_lifetime( 6) = 2.19703E-06   ! Muon-
      part_lifetime( 7) = 8.4E-17       ! Pion0
      part_lifetime( 8) = 2.603E-08     ! Pion+
      part_lifetime( 9) = 2.603E-08     ! Pion-
      part_lifetime(10) = 5.16E-08      ! Kaon 0 long
      part_lifetime(11) = 1.237E-08     ! Kaon+
      part_lifetime(12) = 1.237E-08     ! Kaon-
      part_lifetime(13) = 889.1         ! Neutron
      part_lifetime(14) = 1.0E+30       ! Proton
      part_lifetime(15) = 1.0E+30       ! Antiproton
      part_lifetime(16) = 8.922E-11     ! Kaon 0 short
      part_lifetime(17) = 5.53085E-19   ! Eta
      part_lifetime(18) = 2.632E-10     ! Lambda
      part_lifetime(19) = 7.99E-11      ! Sigma+
      part_lifetime(20) = 7.40E-20      ! Sigma0
      part_lifetime(21) = 1.479E-10     ! Sigma-
      part_lifetime(22) = 2.90E-10      ! Xi 0
      part_lifetime(23) = 1.639E-10     ! Xi -
      part_lifetime(24) = 8.22E-11      ! Omega
      part_lifetime(25) = 889.1         ! Antineutron
      part_lifetime(26) = 2.632E-10     ! Antilambda
      part_lifetime(27) = 7.99E-11      ! Anti-Sigma -
      part_lifetime(28) = 7.40E-20      ! Anti-Sigma 0
      part_lifetime(29) = 1.479E-10     ! Anti-Sigma +
      part_lifetime(30) = 2.90E-10      ! AntiXi 0
      part_lifetime(31) = 1.639E-10     ! AntiXi +
      part_lifetime(32) = 8.22E-11      ! Anti-Omega +
      part_lifetime(33) = 0.0
      part_lifetime(34) = 0.0
      part_lifetime(35) = 0.0
      part_lifetime(36) = 0.0
      part_lifetime(37) = 0.0
      part_lifetime(38) = 0.0
      part_lifetime(39) = 0.0
      part_lifetime(40) = 0.0
      part_lifetime(41) = 0.0
      part_lifetime(42) = 0.0
      part_lifetime(43) = 0.0
      part_lifetime(44) = 0.0
      part_lifetime(45) = 1.0E+30       ! Deuteron
      part_lifetime(46) = 1.0E+30       ! Triton
      part_lifetime(47) = 1.0E+30       ! Alpha
      part_lifetime(48) = 1.0E+30       ! Geantino (Fake particle)
      part_lifetime(49) = 1.0E+30       ! He3
      part_lifetime(50) = 1.0E+30       ! Cerenkov (Fake particle)

      Return
      END

C----------------------------------------------------------------


      subroutine correl_model
      implicit none

CCC   This subroutine computes the requested 2-body model correlation
CCC   function which is to be fitted by the track adjustment procedure.
CCC   The model values are calculated on the requested fine and coarse
CCC   mesh grid in momentum space.  The model values are computed at the
CCC   mid point of each 1D bin or at the center of each 3D cell.  This
CCC   could be refined at a later date to correspond to the integral of
CCC   the model function over the bin width (cell volume) divided by the
CCC   the bin width (cell volume).

C     The model includes the following options which are selected by the
C     'switch*' parameters in common/parameters/:
C
C        switch_1d:          1D model as function of either Qinvar, Qtotal
C                            or Q-vector
C        switch_3d:          3D model as function of either the Bertsch-Pratt
C                            side-out-long kinematics (but no cross term) or
C                            the Yano-Koonin-Podgoretski perp-parallel-time
C                            kinematics.
C        switch_type:        Like and/or Unlike particles
C        switch_coherence:   Purely incoherent source or a mixed incoherent-
C                            coherent source.
C        switch_coulomb:     Either (a) no Coulomb correction, (b) Gamow
C                            factor, (c) NA35 parametrization, or (d) Pratt
C                            Coulomb wave function integration for finite
C                            size, spherical source.
C        switch_fermi_bose:  Fermion or boson identical pairs.

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_correlations.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,j,k

      real*4 R_1dsq, Rsidesq, Routsq, Rlongsq
      real*4 Rperpsq, Rparallelsq, R0sq
      real*4 sqrtlambda,fermi_bose_sign,coulomb_factor,coherence_fac
      real*4 q,q1,q2,q3
      real*4 b,b1,b2,b3
      real*4 massavg

CCC   Set Constants:

      sqrtlambda  = sqrt(abs(lambda))
      R_1dsq      = ((R_1d     /hbc)**2)/2.0
      Rsidesq     = ((Rside    /hbc)**2)/2.0
      Routsq      = ((Rout     /hbc)**2)/2.0
      Rlongsq     = ((Rlong    /hbc)**2)/2.0
      Rperpsq     = ((Rperp    /hbc)**2)/2.0
      Rparallelsq = ((Rparallel/hbc)**2)/2.0
      R0sq        = ((R0       /hbc)**2)/2.0

      fermi_bose_sign = float(switch_fermi_bose)
      coherence_fac   = switch_coherence*2.0*sqrtlambda*
     1                  (1.0 - sqrtlambda)

CCC   Determine average particle pair mass for Coulomb correction:

      massavg = 0.14
      if(mass1.eq.0.0 .and. mass2.gt.0.0) massavg = mass2
      if(mass1.gt.0.0 .and. mass2.eq.0.0) massavg = mass1
      if(mass1.gt.0.0 .and. mass2.gt.0.0) massavg = 0.5*(mass1+mass2)
    
CCC   Compute 1D correlation model arrays:

      If(switch_1d .ge. 1) then
        If(n_1d_fine .gt. 0) then        ! Fill the 1D fine mesh bins
          q = -0.5*binsize_1d_fine
          do i = 1,n_1d_fine
            q = q + binsize_1d_fine
            b = exp(-q*q*R_1dsq)

            if(switch_type.eq.1 .or. switch_type.eq.3) then
              c2mod_like_1d(i) = 1.0 + fermi_bose_sign*(lambda
     1                         *b*b + coherence_fac*b)
              if(switch_coulomb.eq.0) then
                coulomb_factor = 1.0
              else if(switch_coulomb.gt.0) then
                Call coulomb(switch_coulomb,q,1,massavg,Q0,
     1                       coulomb_factor)
              end if
              c2mod_like_1d(i) = coulomb_factor*c2mod_like_1d(i)
            end if

            if(switch_type.eq.2 .or. switch_type.eq.3) then
             c2mod_unlike_1d(i) = 1.0 
             if(switch_coulomb.eq.0) then
               coulomb_factor = 1.0
             else if(switch_coulomb.gt.0) then
               Call coulomb(switch_coulomb,q,-1,massavg,Q0,
     1                      coulomb_factor)
             end if
             c2mod_unlike_1d(i) = coulomb_factor*c2mod_unlike_1d(i)
            end if

          end do    ! End of 1D fine mesh filling do-loop
        end if      ! End of 1D fine mesh option

        If(n_1d_coarse .gt. 0) then      ! Fill the 1D coarse mesh bins
          q = qmid_1d -0.5*binsize_1d_coarse
          do i = n_1d_fine + 1, n_1d_total
            q = q + binsize_1d_coarse
            b = exp(-q*q*R_1dsq)

            if(switch_type.eq.1 .or. switch_type.eq.3) then
              c2mod_like_1d(i) = 1.0 + fermi_bose_sign*(lambda
     1                         *b*b + coherence_fac*b)
              if(switch_coulomb.eq.0) then
                coulomb_factor = 1.0
              else if(switch_coulomb.gt.0) then
                Call coulomb(switch_coulomb,q,1,massavg,Q0,
     1                       coulomb_factor)
              end if
              c2mod_like_1d(i) = coulomb_factor*c2mod_like_1d(i)
            end if

            if(switch_type.eq.2 .or. switch_type.eq.3) then
             c2mod_unlike_1d(i) = 1.0 
             if(switch_coulomb.eq.0) then
               coulomb_factor = 1.0
             else if(switch_coulomb.gt.0) then
               Call coulomb(switch_coulomb,q,-1,massavg,Q0,
     1                      coulomb_factor)
             end if
             c2mod_unlike_1d(i) = coulomb_factor*c2mod_unlike_1d(i)
            end if

          end do    ! End of 1D coarse mesh filling do-loop
        end if      ! End of 1D coarse mesh option
      end if        ! End of 1D option

CCC   Compute 3D correlation model arrays:

      If(switch_3d .ge. 1) Then
       If(n_3d_fine .gt. 0) then     !  Fill the 3D fine mesh bins
        q1 = -0.5*binsize_3d_fine
        do i = 1,n_3d_fine
         q1 = q1 + binsize_3d_fine
         if(switch_3d.eq.1) b1=exp(-q1*q1*Rsidesq)
         if(switch_3d.eq.2) b1=exp(-q1*q1*Rperpsq)

          q2 = -0.5*binsize_3d_fine
          do j = 1,n_3d_fine
           q2 = q2 + binsize_3d_fine
           if(switch_3d.eq.1) b2=exp(-q2*q2*Routsq)
           if(switch_3d.eq.2) b2=exp(-q2*q2*Rparallelsq)

            q3 = -0.5*binsize_3d_fine
            do k = 1,n_3d_fine
             q3 = q3 + binsize_3d_fine
             if(switch_3d.eq.1) b3=exp(-q3*q3*Rlongsq)
             if(switch_3d.eq.2) b3=exp(-q3*q3*R0sq)

             b = b1*b2*b3
             if(switch_3d.eq.1) q = sqrt(q1*q1+q2*q2+q3*q3)
             if(switch_3d.eq.2) q = sqrt(q1*q1+q2*q2)

             if(switch_type.eq.1 .or. switch_type.eq.3) then
              c2mod_like_3d_fine(i,j,k) = 1.0 + fermi_bose_sign*(lambda
     1                         *b*b + coherence_fac*b)
              if(switch_coulomb.eq.0) then
                coulomb_factor = 1.0
              else if(switch_coulomb.gt.0) then
                Call coulomb(switch_coulomb,q,1,massavg,Q0,
     1                       coulomb_factor)
              end if
              c2mod_like_3d_fine(i,j,k) =
     1           coulomb_factor*c2mod_like_3d_fine(i,j,k)
             end if

             if(switch_type.eq.2 .or. switch_type.eq.3) then
              c2mod_unlike_3d_fine(i,j,k) = 1.0 
              if(switch_coulomb.eq.0) then
                coulomb_factor = 1.0
              else if(switch_coulomb.gt.0) then
                Call coulomb(switch_coulomb,q,-1,massavg,Q0,
     1                       coulomb_factor)
              end if
              c2mod_unlike_3d_fine(i,j,k) =
     1           coulomb_factor*c2mod_unlike_3d_fine(i,j,k)
             end if

            end do
           end do
          end do     ! End of 3D Fine Mesh Filling do-loops
         end if      ! End 3D fine mesh option

       If(n_3d_coarse .gt. 0) then     !  Fill the 3D coarse mesh bins
        q1 = -0.5*binsize_3d_coarse
        do i = 1,n_3d_coarse
         q1 = q1 + binsize_3d_coarse
         if(switch_3d.eq.1) b1=exp(-q1*q1*Rsidesq)
         if(switch_3d.eq.2) b1=exp(-q1*q1*Rperpsq)

          q2 = -0.5*binsize_3d_coarse
          do j = 1,n_3d_coarse
           q2 = q2 + binsize_3d_coarse
           if(switch_3d.eq.1) b2=exp(-q2*q2*Routsq)
           if(switch_3d.eq.2) b2=exp(-q2*q2*Rparallelsq)

            q3 = -0.5*binsize_3d_coarse
            do k = 1,n_3d_coarse
             q3 = q3 + binsize_3d_coarse
             if(switch_3d.eq.1) b3=exp(-q3*q3*Rlongsq)
             if(switch_3d.eq.2) b3=exp(-q3*q3*R0sq)

             b = b1*b2*b3
             if(switch_3d.eq.1) q = sqrt(q1*q1+q2*q2+q3*q3)
             if(switch_3d.eq.2) q = sqrt(q1*q1+q2*q2)

             if(switch_type.eq.1 .or. switch_type.eq.3) then
              c2mod_like_3d_coarse(i,j,k) = 1.0+fermi_bose_sign*(lambda
     1                         *b*b + coherence_fac*b)
              if(switch_coulomb.eq.0) then
                coulomb_factor = 1.0
              else if(switch_coulomb.gt.0) then
                Call coulomb(switch_coulomb,q,1,massavg,Q0,
     1                       coulomb_factor)
              end if
              c2mod_like_3d_coarse(i,j,k) =
     1           coulomb_factor*c2mod_like_3d_coarse(i,j,k)
             end if

             if(switch_type.eq.2 .or. switch_type.eq.3) then
              c2mod_unlike_3d_coarse(i,j,k) = 1.0 
              if(switch_coulomb.eq.0) then
                coulomb_factor = 1.0
              else if(switch_coulomb.gt.0) then
                Call coulomb(switch_coulomb,q,-1,massavg,Q0,
     1                       coulomb_factor)
              end if
              c2mod_unlike_3d_coarse(i,j,k) =
     1           coulomb_factor*c2mod_unlike_3d_coarse(i,j,k)
             end if

            end do
           end do
          end do     ! End of 3D Coarse Mesh Filling do-loops
          c2mod_like_3d_coarse(1,1,1) = 0.0
          c2mod_unlike_3d_coarse(1,1,1) = 0.0
        end if       ! End of 3D Coarse Mesh Option

      End If         ! End of 3D Option

      Return
      END

C----------------------------------------------------------------------


      subroutine coulomb(control,q,sign,mass,Q0,factor)
      implicit none

CCC   Compute Coulomb correction to the two-body correlation functions
C     for like and unlike charges for particles of the same mass.
C     Three methods are allowed:
C
C     If control = 1, Then use the gamow factor for point sources
C     If control = 2, Then use the NA35 finite source size empirical
C                     correction factor from eq.(5) in Z. Phys. C73,
C                     443 (1997).
C     If control = 3, Then use the Pratt finite source size, numerically
C                     integrated Coulomb correction factor with inter-
C                     polated tables.
C
C     Other parameters in the argument list are:
C
C        q     = 3-vector momentum difference for track pair, in GeV/c.
C        sign  = algebraic sign of the charge product for the track pair.
C        mass  = particle mass in GeV, it is assumed that both particles
C                have the same mass, e.g. pi+ and pi-, but not K+ and pi-.
C        Q0    = NA35 parameter in GeV/c if control = 2
C              = Source radius in fm if control = 3
C       factor = Multiplicative Coulomb correction result which is
C                calculated here and returned to the calling program.
C

      Include 'common_coulomb.inc'

CCC   Local Variable Type Declarations:

      integer*4 control, sign

      real*4 pi,q,mass,Q0,factor,alpha,eta,eta2pi
      real*4 gamow
      parameter (pi = 3.141592654)
      parameter (alpha = 0.00729735)

CCC   Compute Gamow factor for control options 1 and 2:

      if(control.eq.1 .or. control.eq.2) then
         if(q .le. 0.001) then
            if(sign .gt. 0) gamow = 0.0
            if(sign .lt. 0) gamow = 86.0
         else
            eta = sign*mass*alpha/q
            eta2pi = 2.0*pi*eta
            gamow = eta2pi/(exp(eta2pi) - 1.0)
         end if
      end if

CCC   Compute Coulomb Correction factor for options 1, 2 and 3:

      if(control .eq. 1) then
         factor = gamow
      else if(control .eq. 2) then
         factor = 1.0 + (gamow - 1.0)*exp(-q/Q0)
      else if(control .eq. 3) then

         if(q .le. q_coul(1)) then
            if(sign .gt. 0) factor = c2_coul_like(1)
            if(sign .lt. 0) factor = c2_coul_unlike(1)
         else if(q .ge. q_coul(max_c2_coul - 1)) then
            if(sign .gt. 0) factor = c2_coul_like(max_c2_coul - 1)
            if(sign .lt. 0) factor = c2_coul_unlike(max_c2_coul - 1)
         else
            if(sign .gt. 0) then
               Call LAGRNG1(q,q_coul,factor,c2_coul_like,
     1                     max_c2_coul,1,5,max_c2_coul,1)
            else if(sign .lt. 0) then
               Call LAGRNG1(q,q_coul,factor,c2_coul_unlike,
     1                     max_c2_coul,1,5,max_c2_coul,1)
            end if
         end if
      end if  ! END Coulomb correction evaluation, control selection opt.

      Return
      END

C---------------------------------------------------------------------


      SUBROUTINE LAGRNG1 (X,ARG,Y,VAL,NDIM,NFS,NPTS,MAXARG,MAXFS)
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

C-------------------------------------------------------------------


      subroutine Hbtp_kin(px,py,pz,E,pt,phi,eta,mass,control)
      implicit none

CCC   Four-momentum kinematics conversion:
C
C        If control = 1, use input {px,py,pz,mass} to calculate
C                        {E,pt,phi,eta}
C        If control = 2, use input {pt,phi,eta,mass} to calculate
C                        {px,py,pz,E}
C
C        Units:  Momentum are in GeV/c
C                Energy and mass are in GeV
C                Angles are in degrees

CCC   Local Variable Type Declarations:

      integer*4 control

      real*4 px,py,pz,E,pt,phi,eta,mass
      real*4 theta,pi,rad,pcut,x,y
      parameter (pi = 3.141592654)
      parameter (pcut = 0.000001)

      rad = 180.0/pi

      If(control .eq. 1) Then   !   Use {px,py,pz,mass} --> {E,pt,phi,eta}
         pt = sqrt(px*px + py*py)
         E  = sqrt(pt*pt + pz*pz + mass*mass)

CCC   Compute azimuthal angle phi; treat pt = 0.0 and py = 0.0 cases
CCC   separate.

         if(pt .le. pcut) then
            phi = 0.0
         else if(pt.gt.pcut.and.abs(py).le.pcut.and.px.lt.0.0) then
            phi = pi
         else
            phi = atan2(py,px)
         end if
         if(phi .lt. 0.0) phi = phi + 2.0*pi
         phi = phi*rad

CCC   Compute pseudorapidity:

         if(pt.le.pcut .and. abs(pz).le.pcut) then
            eta = 0.0
         else if(pt.le.pcut .and. abs(pz).gt.pcut) then
            eta = 0.5*log((E+pz)/(E-pz))   !  Use beam rapidity
         else
            theta = atan2(pt,pz)
            eta   = -log(tan(theta/2.0))
         end if

      Else If(control .eq. 2) Then  ! Use {pt,phi,eta,mass} --> {E,px,py,pz}

         px = pt*cos(phi/rad)
         py = pt*sin(phi/rad)
         if(abs(eta) .le. pcut) then
            pz = 0.0
         else
            x = exp(-eta)
            y = atan(x)
            theta = 2.0*y
C            theta = 2.0*atan(exp(-eta))
            pz = pt/tan(theta)
         end if

         E = sqrt(pt*pt + pz*pz + mass*mass)

      End If   !  End control options

      Return
      END

C----------------------------------------------------------------------


      subroutine qdiff(px1,py1,pz1,E1,px2,py2,pz2,E2,qinvar2,qtotal2,
     1           qvector2,qside2,qout2,qlong2,qperp2,qtime2)
      implicit none

CCC   This subroutine computes the various relative momenta for given
CCC   input 4-momentum for particles 1 and 2.  The subr: returns the
CCC   square of the momentum. All energy and momenta are in GeV.
C
C     Input 4-momentum for particle 1: {px1,py1,pz1,E1}
C     Input 4-momentum for particle 2: {px2,py2,pz2,E2}
C
CCC   Computed Momentum Differences are the following:
C
C        qinvar2   = Q-invariant**2
C        qtotal2   = Q-Total**2      (space**2 + time**2)
C        qvector2  = 3-momentum vector difference squared
C        qside2    = q-side**2 of Bertsch-Pratt 3D source models
C        qout2     = q-out**2  of Bertsch-Pratt 3D source models
C        qlong2    = q-long**2 of Bertsch-Pratt 3D source models
C        qperp2    = q-perpendicular**2 of YKP 3D source models
C        qparallel2= q-parallel**2 of YKP 3D source models
C                  = qlong2 (not assigned a separate variable name)
C        qtime2    = q-time-like**2 of YKP 3D source models
    
CCC   Local Variable Type Declarations:

      real*4 px1,py1,pz1,E1,px2,py2,pz2,E2
      real*4 qinvar2,qtotal2,qvector2
      real*4 qside2,qout2,qlong2
      real*4 qperp2,qtime2
      real*4 px12sq,py12sq,pz12sq,E12sq

      px12sq = (px1 - px2)**2
      py12sq = (py1 - py2)**2
      pz12sq = (pz1 - pz2)**2
      E12sq  = (E1  -  E2)**2
        
      qvector2   = px12sq + py12sq + pz12sq
      qinvar2    = qvector2 - E12sq
      qtotal2    = qvector2 + E12sq

      qlong2     = pz12sq
      qout2      = (px1*px1 - px2*px2 + py1*py1 - py2*py2)
      qout2      = qout2*qout2/((px1+px2)**2 + (py1+py2)**2)
      qperp2     = px12sq + py12sq
      qside2     = qperp2 - qout2
      qtime2     = E12sq
      if (qside2 .lt. 0) then
C        write(*,*) 'qside2 is less then 0', qside2
C	write(*,*) ' qperp2, qout2', qperp2, qout2
C	write(*,*) ' px1,py1,pz1,E1 ',px1,py1,pz1,E1
C	write(*,*) ' px2,py2,pz2,E2 ',px2,py2,pz2,E2
	qside2 = 0.0
      end if
      Return
      END

C----------------------------------------------------------------------- 


      subroutine mean_rms(a,ndim,npts,mean,rms)
      implicit none

CCC   Calculate the mean and standard deviation (rms) for input
C     distribution a() for npts number of values.
C     ndim = dimension of array a() in calling program.

CCC   Local Variable Type Declarations:

      integer*4 ndim, npts, i
      real*4 a(ndim), mean, rms, sum_mean, sum_rms

      if(npts .le. 0) then
         mean = 0.0
         rms  = 0.0
         return
      else if(npts .eq. 1) then
         mean = a(1)
         rms  = 0.0
         return
      else
         sum_mean = 0.0
         sum_rms  = 0.0
         do i = 1,npts
         sum_mean = sum_mean + a(i)
         end do
         mean = sum_mean/float(npts)

         do i = 1,npts
         sum_rms  = sum_rms + (a(i) - mean)**2
         end do
         rms = sqrt((sum_rms)/float(npts - 1))
         return
      end if

      END

C-----------------------------------------------------------------------


      subroutine tindex(mode,track_id)
      implicit none

CCC   This subroutine locates tracks in {px,py,pz} sectors
C     and sets the sector index numbers
C     in track table 'trk' and/or 'trk2', depending on the value of 'mode'.
C     If a track's momentum is out of the sector ranges, then the track
C     will be assigned to, and counted in the nearest sector cell on the
C     edge or corner.
C
C     If mode = 1, apply to tracks in table 'trk'
C     If mode = 2, apply to tracks in table 'trk2'
C
C     If track_id = 0, then do this for all tracks
C     If track_id = i (where i.gt.0) then do this for track row i only

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,track_id,mode

C-----------------------
      If(mode.eq.1) Then
C-----------------------

       If(track_id .eq. 0) Then
         do i = 1,n_part_tot_trk
         trk_px_sec(i) = int(((trk_px(i) - px_min)/delpx)+1.00001)
            if(trk_px_sec(i) .lt.1) trk_px_sec(i) = 1
            if(trk_px_sec(i) .gt. n_px_bins) trk_px_sec(i) = n_px_bins
         trk_py_sec(i) = int(((trk_py(i) - py_min)/delpy)+1.00001)
            if(trk_py_sec(i) .lt.1) trk_py_sec(i) = 1
            if(trk_py_sec(i) .gt. n_py_bins) trk_py_sec(i) = n_py_bins
         trk_pz_sec(i) = int(((trk_pz(i) - pz_min)/delpz)+1.00001)
            if(trk_pz_sec(i) .lt.1) trk_pz_sec(i) = 1
            if(trk_pz_sec(i) .gt. n_pz_bins) trk_pz_sec(i) = n_pz_bins
         trk_sector(i) = trk_px_sec(i) + (trk_py_sec(i) - 1)*
     1      n_px_bins + (trk_pz_sec(i) - 1)*n_px_bins*n_py_bins
         end do

       Else If(track_id .gt. 0) Then
         i = track_id
         trk_px_sec(i) = int(((trk_px(i) - px_min)/delpx)+1.00001)
            if(trk_px_sec(i) .lt.1) trk_px_sec(i) = 1
            if(trk_px_sec(i) .gt. n_px_bins) trk_px_sec(i) = n_px_bins
         trk_py_sec(i) = int(((trk_py(i) - py_min)/delpy)+1.00001)
            if(trk_py_sec(i) .lt.1) trk_py_sec(i) = 1
            if(trk_py_sec(i) .gt. n_py_bins) trk_py_sec(i) = n_py_bins
         trk_pz_sec(i) = int(((trk_pz(i) - pz_min)/delpz)+1.00001)
            if(trk_pz_sec(i) .lt.1) trk_pz_sec(i) = 1
            if(trk_pz_sec(i) .gt. n_pz_bins) trk_pz_sec(i) = n_pz_bins
         trk_sector(i) = trk_px_sec(i) + (trk_py_sec(i) - 1)*
     1      n_px_bins + (trk_pz_sec(i) - 1)*n_px_bins*n_py_bins
       End If

C-----------------------------
      Else If (mode.eq.2) Then
C-----------------------------

       If(track_id .eq. 0) Then
         do i = 1,n_part_tot_trk2
         trk2_px_sec(i) = int(((trk2_px(i) - px_min)/delpx)+1.00001)
            if(trk2_px_sec(i) .lt.1) trk2_px_sec(i) = 1
            if(trk2_px_sec(i) .gt. n_px_bins) trk2_px_sec(i) = n_px_bins
         trk2_py_sec(i) = int(((trk2_py(i) - py_min)/delpy)+1.00001)
            if(trk2_py_sec(i) .lt.1) trk2_py_sec(i) = 1
            if(trk2_py_sec(i) .gt. n_py_bins) trk2_py_sec(i) = n_py_bins
         trk2_pz_sec(i) = int(((trk2_pz(i) - pz_min)/delpz)+1.00001)
            if(trk2_pz_sec(i) .lt.1) trk2_pz_sec(i) = 1
            if(trk2_pz_sec(i) .gt. n_pz_bins) trk2_pz_sec(i) = n_pz_bins
         trk2_sector(i) = trk2_px_sec(i) + (trk2_py_sec(i) - 1)*
     1      n_px_bins + (trk2_pz_sec(i) - 1)*n_px_bins*n_py_bins
         end do

       Else If(track_id .gt. 0) Then
         i = track_id
         trk2_px_sec(i) = int(((trk2_px(i) - px_min)/delpx)+1.00001)
            if(trk2_px_sec(i) .lt.1) trk2_px_sec(i) = 1
            if(trk2_px_sec(i) .gt. n_px_bins) trk2_px_sec(i) = n_px_bins
         trk2_py_sec(i) = int(((trk2_py(i) - py_min)/delpy)+1.00001)
            if(trk2_py_sec(i) .lt.1) trk2_py_sec(i) = 1
            if(trk2_py_sec(i) .gt. n_py_bins) trk2_py_sec(i) = n_py_bins
         trk2_pz_sec(i) = int(((trk2_pz(i) - pz_min)/delpz)+1.00001)
            if(trk2_pz_sec(i) .lt.1) trk2_pz_sec(i) = 1
            if(trk2_pz_sec(i) .gt. n_pz_bins) trk2_pz_sec(i) = n_pz_bins
         trk2_sector(i) = trk2_px_sec(i) + (trk2_py_sec(i) - 1)*
     1      n_px_bins + (trk2_pz_sec(i) - 1)*n_px_bins*n_py_bins
       End If

C------------
      End If
C------------

      Return
      END

C------------------------------------------------------------------------


      subroutine stm_build(mode,track_index,old_sector)
      implicit none

CCC   This subroutine fills or updates the track-sector information
C     table sec_trk_map or, for the reference calculations, it fills
C     sec_trk_map2.  These track-sector tables contain the information
C     about the track occupancy, status, etc. for all of the {px,py,pz}
C     sectors.
C
C     For Mode = 1:
C
C        If track_index = 0, then fill information for all tracks in 'trk',
C                            into table 'stm'
C        If track_index = i (where i.gt.0) then fill only the track-sector
C                           information for track i in 'trk', into table 'stm'
C                           Also for this case, if old_sector .ne. 0, then
C                           remove the track # and ID information for this
C                           old sector # from table stm
C
C     For Mode = 2: Fill information for all tracks in 'trk2' into table 'stm2'

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,j,mode,track_index,old_sector,row
      integer*4 temp(max_trk_sec)

C------------------------
      IF (mode.eq.1) Then
C------------------------

       If (track_index .eq. 0) Then
         do i = 1,sec_maxlen
            stm_sec_id(i) = 0
            stm_n_trk_sec(i) = 0
            do j = 1,max_trk_sec
               stm_track_id(j,i) = 0
            end do
            stm_flag(i) = 0
         end do
         do i = 1,n_sectors
         stm_sec_id(i) = i
         end do
         do i = 1,n_part_tot_trk
         if(trk_flag(i) .eq. 0) then
         row = trk_sector(i)
         stm_n_trk_sec(row) = stm_n_trk_sec(row) + 1
         if(stm_n_trk_sec(row) .le. max_trk_sec) then
            stm_track_id(stm_n_trk_sec(row),row) = trk_id(i)
            stm_flag(row) = 0
            trk_flag(i)   = 0
         else
            stm_n_trk_sec(row) = stm_n_trk_sec(row) - 1
            stm_flag(row) = 1
            trk_flag(i)   = 1
            trk_sector(i) = 0
         end if
         end if
         end do

       Else If (track_index .gt. 0) Then

         if(old_sector .ne. 0) then
CCC   Remove track from old sector:
            j = 0
            do i = 1,stm_n_trk_sec(old_sector)
               if(stm_track_id(i,old_sector) .ne. track_index) then
                  j = j + 1
                  temp(j) = stm_track_id(i,old_sector)
               end if
            end do
            stm_n_trk_sec(old_sector) = j
            do i = 1,max_trk_sec
               stm_track_id(i,old_sector) = 0
            end do
            do i = 1,stm_n_trk_sec(old_sector)
               stm_track_id(i,old_sector) = temp(i)
            end do
         end if
CCC   Update with new sector location of track:
         i = track_index
         if(trk_flag(i) .eq. 0) then
         row = trk_sector(i)
         stm_n_trk_sec(row) = stm_n_trk_sec(row) + 1
         if(stm_n_trk_sec(row) .le. max_trk_sec) then
            stm_track_id(stm_n_trk_sec(row),row) = trk_id(i)
            stm_flag(row) = 0
            trk_flag(i)   = 0
         else
            stm_n_trk_sec(row) = stm_n_trk_sec(row) - 1
            stm_flag(row) = 1
            trk_flag(i)   = 1
            trk_sector(i) = 0
         end if
         end if
       End If

C-----------------------------
      Else If (mode.eq.2) Then
C-----------------------------

       If (track_index .eq. 0) Then
         do i = 1,sec_maxlen2
            stm2_sec_id(i) = 0
            stm2_n_trk_sec(i) = 0
            do j = 1,max_trk_sec2
               stm2_track_id(j,i) = 0
            end do
            stm2_flag(i) = 0
         end do
         do i = 1,n_sectors
         stm2_sec_id(i) = i
         end do
         do i = 1,n_part_tot_trk2
         if(trk2_flag(i) .eq. 0) then
         row = trk2_sector(i)
         stm2_n_trk_sec(row) = stm2_n_trk_sec(row) + 1
         if(stm2_n_trk_sec(row) .le. max_trk_sec2) then
            stm2_track_id(stm2_n_trk_sec(row),row) = trk2_id(i)
            stm2_flag(row) = 0
            trk2_flag(i)   = 0
         else
            stm2_n_trk_sec(row) = stm2_n_trk_sec(row) - 1
            stm2_flag(row) = 1
            trk2_flag(i)   = 1
            trk2_sector(i) = 0
         end if
         end if
         end do
       end if

C------------
      End If   !  End mode = 1,2 selection options
C------------

      Return
      END

C-----------------------------------------------------------------------


      subroutine sec_index(index,nbins,index_min,index_max)
      implicit none

CCC   Calculate track-sector neighboring bins and min->max range:

CCC   Local Variable Type Declarations:

      integer*4 index,nbins,index_min,index_max

      index_min = index - 1
      if(index_min .lt. 1) index_min = 1
      index_max = index + 1
      if(index_max .gt. nbins) index_max = nbins

      Return
      END

C-----------------------------------------------------------------------


      subroutine dist_range(mode,ntracks_out,ntracks_flagged)
      implicit none

CCC   Determine if tracks are out of acceptance range in pt, phi and eta,
C     and, if so, then set the 'out_flag' variable in the track table 'trk'

CCC   For Mode = 1, use track table 'trk'
CCC   For Mode = 2, use track table 'trk2'

CCC   Count the number of flagged tracks, i.e. trk(i).flag = 1, and "out
C     of acceptance range" tracks, i.e. trk(i).out_flag = 1, for both
C     particle ID types.  Determine the number of tracks to use in the
C     correlation fit for each particle ID type.

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,mode,ntracks_out,ntracks_flagged

C------------------------
      If (mode.eq.1) Then
C------------------------

      do i = 1,trk_maxlen
      trk_out_flag(i) = 0
      end do

      ntracks_flagged = 0

      do i = 1,n_part_tot_trk
      if(trk_flag(i) .eq. 0) then
      if(trk_pt(i) .lt. pt_min .or. trk_pt(i) .gt. pt_max)
     1   trk_out_flag(i)=1
      if(trk_phi(i).lt.phi_min .or. trk_phi(i).gt.phi_max)
     1   trk_out_flag(i)=1
      if(trk_eta(i).lt.eta_min .or. trk_eta(i).gt.eta_max)
     1   trk_out_flag(i)=1
      else if(trk_flag(i) .eq. 1) then
      ntracks_flagged = ntracks_flagged + 1
      end if
      end do

      ntracks_out = 0
      do i = 1,n_part_tot_trk
      if(trk_out_flag(i) .eq. 1) ntracks_out = ntracks_out + 1
      end do

      n_part_used_1_trk = 0
      n_part_used_2_trk = 0
      do i = 1,n_part_tot_trk
         if(trk_flag(i) .eq. 0) then
            if(trk_ge_pid(i) .eq. pid(1)) then
               n_part_used_1_trk = n_part_used_1_trk + 1
            else if(trk_ge_pid(i) .eq. pid(2)) then
               n_part_used_2_trk = n_part_used_2_trk + 1
            end if
         end if
      end do 

C-----------------------------
      Else If (mode.eq.2) Then
C-----------------------------

      do i = 1,trk2_maxlen
      trk2_out_flag(i) = 0
      end do

      ntracks_flagged = 0

      do i = 1,n_part_tot_trk2
      if(trk2_flag(i) .eq. 0) then
      if(trk2_pt(i) .lt. pt_min .or. trk2_pt(i) .gt. pt_max)
     1   trk2_out_flag(i)=1
      if(trk2_phi(i).lt.phi_min .or. trk2_phi(i).gt.phi_max)
     1   trk2_out_flag(i)=1
      if(trk2_eta(i).lt.eta_min .or. trk2_eta(i).gt.eta_max)
     1   trk2_out_flag(i)=1
      else if(trk2_flag(i) .eq. 1) then
      ntracks_flagged = ntracks_flagged + 1
      end if
      end do

      ntracks_out = 0
      do i = 1,n_part_tot_trk2
      if(trk2_out_flag(i) .eq. 1) ntracks_out = ntracks_out + 1
      end do

      n_part_used_1_trk2 = 0
      n_part_used_2_trk2 = 0
      do i = 1,n_part_tot_trk2
         if(trk2_flag(i) .eq. 0) then
            if(trk2_ge_pid(i) .eq. pid(1)) then
               n_part_used_1_trk2 = n_part_used_1_trk2 + 1
            else if(trk2_ge_pid(i) .eq. pid(2)) then
               n_part_used_2_trk2 = n_part_used_2_trk2 + 1
            end if
         end if
      end do 

C------------
      End If  !  End mode = 1,2 selection option
C------------

      Return
      END

C--------------------------------------------------------------------


      subroutine histog1(mode,itrack,pid_index,pidnum,pt_save,
     1                   phi_save,eta_save)
      implicit none

CCC   This subroutine computes and/or updates the one-body histograms
C     according to the selected 'mode' and for the requested particle
C     ID type.

C     Note:  If the track momentum is out-of-range in {pt,phi,eta},
C            then it is ignored.  The {pt,phi,eta} dependences for 
C            the 1-dimensional histogramming are treated independently.
C            It is therefore possible for the sum of the number of
C            particles in the pt, phi and eta one-body, 1D histograms
C            to be unequal.

CCC   Mode = 1,  Fill the one-body histograms (hist1*) for selected
C                particle ID type, for the initial input distribution,
C                using the momenta in 'trk'
C
C     Mode = 2,  Remove particle 'itrack' from temporary one-body hist-
C                ogram (htmp1*) for selected particle ID type, using
C                momentum values given by pt_save, phi_save, eta_save.
C
C     Mode = 3,  Add particle 'itrack' to the temporary one-body hist-
C                ogram (htmp1*) for selected particle ID type, using
C                momentum values in track table 'trk'. 
C
C     Mode = 4,  Fill the one-body histograms (hist1*) for selected
C                particle ID type, for the initial input distribution,
C                using the momenta in 'trk2'
C
C     itrack =   track index # for the removed or added track for mode =
C                2 or 3, respectively.
C
C     pid_index = 1 or 2 for the first or second particle ID type, and
C                 for filling/update either hist1*1 or hist1*2, and
C                 similarly for htmp1*1 or htmp1*2
C
C     pidnum    =   Geant particle ID # for the track(s) to be filled or
C                updated.
C
C     pt_save  = Removed track's pt value.
C
C     phi_save = Removed track's phi value.
C
C     eta_save = Removed track's eta value.

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'

CCC   Local Variable Type Declarations:

      integer*4 mode,itrack,i,pid_index,pidnum,index
      integer*4 trk_counter,trk2_counter

      real*4    pt_save,phi_save,eta_save
      real*4    delpt,delphi,deleta

C-------------------------
      If (mode.eq.1) Then
C-------------------------

CCC   Fill one-body histograms for requested particle ID from table 'trk'
CCC   Initialize necessary arrays to zero:

         if(pid_index .eq. 1) then
            do i = 1,max_h_1d
               hist1_pt_1(i)  = 0
               hist1_phi_1(i) = 0
               hist1_eta_1(i) = 0
            end do
         else if(pid_index .eq. 2) then
            do i = 1,max_h_1d
               hist1_pt_2(i)  = 0
               hist1_phi_2(i) = 0
               hist1_eta_2(i) = 0
            end do
         end if

         trk_counter = 0

         do i = 1,n_part_tot_trk
            if(trk_ge_pid(i).eq.pidnum .and. trk_flag(i).eq.0) then
               trk_counter = trk_counter + 1
               delpt  = trk_pt(i)   - pt_min
               delphi = trk_phi(i)  - phi_min
               deleta = trk_eta(i)  - eta_min

               index = int((delpt/pt_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) hist1_pt_1(index) =
     1                               hist1_pt_1(index) + 1
                  if(pid_index.eq.2) hist1_pt_2(index) =
     1                               hist1_pt_2(index) + 1
               end if

               index = int((delphi/phi_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) hist1_phi_1(index) =
     1                               hist1_phi_1(index) + 1
                  if(pid_index.eq.2) hist1_phi_2(index) =
     1                               hist1_phi_2(index) + 1
               end if

               index = int((deleta/eta_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) hist1_eta_1(index) =
     1                               hist1_eta_1(index) + 1
                  if(pid_index.eq.2) hist1_eta_2(index) =
     1                               hist1_eta_2(index) + 1
               end if

            end if
         end do

         if(pid_index .eq. 1) n_part_used_1_trk = trk_counter
         if(pid_index .eq. 2) n_part_used_2_trk = trk_counter

C--------------------------------
      Else If (mode .eq. 2) Then
C--------------------------------

CCC   Remove track # 'itrack' from fitting histograms in htmp1*,
CCC   use pt_save, phi_save, eta_save for the old momentum values
CCC   in order to determine which bins to decrement.

         if(trk_ge_pid(itrack).eq.pidnum.and.trk_flag(itrack).eq.0)then
               delpt  = pt_save  - pt_min
               delphi = phi_save - phi_min
               deleta = eta_save - eta_min

               index = int((delpt/pt_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) htmp1_pt_1(index) =
     1                               htmp1_pt_1(index) - 1
                  if(pid_index.eq.2) htmp1_pt_2(index) =
     1                               htmp1_pt_2(index) - 1
               end if

               index = int((delphi/phi_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) htmp1_phi_1(index) =
     1                               htmp1_phi_1(index) - 1
                  if(pid_index.eq.2) htmp1_phi_2(index) =
     1                               htmp1_phi_2(index) - 1
               end if

               index = int((deleta/eta_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) htmp1_eta_1(index) =
     1                               htmp1_eta_1(index) - 1
                  if(pid_index.eq.2) htmp1_eta_2(index) =
     1                               htmp1_eta_2(index) - 1
               end if

         end if

C--------------------------------
      Else If (mode .eq. 3) Then
C--------------------------------

CCC   Add track # 'itrack' to fitting histograms in htmp1*,
CCC   use pt, phi and eta values in track table 'trk(itrack)'
CCC   for the new/added track position.

         if(trk_ge_pid(itrack).eq.pidnum.and.trk_flag(itrack).eq.0)then
               delpt  = trk_pt(itrack)  - pt_min
               delphi = trk_phi(itrack) - phi_min
               deleta = trk_eta(itrack) - eta_min

               index = int((delpt/pt_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) htmp1_pt_1(index) =
     1                               htmp1_pt_1(index) + 1
                  if(pid_index.eq.2) htmp1_pt_2(index) =
     1                               htmp1_pt_2(index) + 1
               end if

               index = int((delphi/phi_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) htmp1_phi_1(index) =
     1                               htmp1_phi_1(index) + 1
                  if(pid_index.eq.2) htmp1_phi_2(index) =
     1                               htmp1_phi_2(index) + 1
               end if

               index = int((deleta/eta_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) htmp1_eta_1(index) =
     1                               htmp1_eta_1(index) + 1
                  if(pid_index.eq.2) htmp1_eta_2(index) =
     1                               htmp1_eta_2(index) + 1
               end if

         end if

C------------------------------
      Else If (mode.eq.4) Then
C------------------------------

CCC   Fill one-body histograms for requested particle ID from table 'trk2'
CCC   Initialize necessary arrays to zero:

         if(pid_index .eq. 1) then
            do i = 1,max_h_1d
               hist1_pt_1(i)  = 0
               hist1_phi_1(i) = 0
               hist1_eta_1(i) = 0
            end do
         else if(pid_index .eq. 2) then
            do i = 1,max_h_1d
               hist1_pt_2(i)  = 0
               hist1_phi_2(i) = 0
               hist1_eta_2(i) = 0
            end do
         end if

         trk2_counter = 0

         do i = 1,n_part_tot_trk2
            if(trk2_ge_pid(i).eq.pidnum .and. trk2_flag(i).eq.0) then
               trk2_counter = trk2_counter + 1
               delpt  = trk2_pt(i)   - pt_min
               delphi = trk2_phi(i)  - phi_min
               deleta = trk2_eta(i)  - eta_min

               index = int((delpt/pt_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) hist1_pt_1(index) =
     1                               hist1_pt_1(index) + 1
                  if(pid_index.eq.2) hist1_pt_2(index) =
     1                               hist1_pt_2(index) + 1
               end if

               index = int((delphi/phi_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) hist1_phi_1(index) =
     1                               hist1_phi_1(index) + 1
                  if(pid_index.eq.2) hist1_phi_2(index) =
     1                               hist1_phi_2(index) + 1
               end if

               index = int((deleta/eta_bin_size) + 0.99999)
               if(index.ge.1 .and. index.le.max_h_1d) then
                  if(pid_index.eq.1) hist1_eta_1(index) =
     1                               hist1_eta_1(index) + 1
                  if(pid_index.eq.2) hist1_eta_2(index) =
     1                               hist1_eta_2(index) + 1
               end if

            end if
         end do

         if(pid_index .eq. 1) n_part_used_1_trk2 = trk2_counter
         if(pid_index .eq. 2) n_part_used_2_trk2 = trk2_counter

C------------
      End If   !   End Mode = 1,2,3,4 Selection Options
C------------

      Return
      END

C-----------------------------------------------------------------------

      subroutine histog2(mode,itrack,px_sec_save,py_sec_save,
     1           pz_sec_save,px_save,py_save,pz_save,E_save)
      implicit none

CCC   This subroutine computes and/or updates the two-body histograms
C     according to the selected 'mode' and for the necessary particle
C     ID type(s).

C     Mode = 1,  Fill the two-body histograms (hist*) for all particles
C                in table 'trk' for like and unlike pairs, for 1D and/or
C                3D fine and 3D coarse mesh bins.
C
C     Mode = 2,  Remove all old track pairs for 'itrack' particle from 
C                all htmp* histograms, for particles in table 'trk', for
C                like and unlike pairs, for 1D and/or 3D fine and 3D coarse
C                mesh bins; using the saved momentum and track values.
C
C     Mode = 3,  Add all new track pairs for 'itrack' particle to
C                all htmp* histograms, for particles in table 'trk', for
C                like and unlike pairs, for 1D and/or 3D fine and 3D coarse
C                mesh bins; using the values in table 'trk(itrack).*'
C
C     Mode = 4,  Fill and accumulate reference histograms (href*) for all
C                particle pairs from tables 'trk' and 'trk2', for like and
C                unlike pairs, for 1D and/or 3D fine and 3D coarse
C                mesh bins.
C
C     itrack =   single track index in table 'trk' for pairs to be removed
C                (mode = 2) or added (mode = 3).

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'

CCC   Local Variable Type Declarations:

      integer*4 mode,itrack,i,j,k,jx,jy,jz
      integer*4 jsec,trkj_sector,imin,imax,njtrks
      integer*4 index1,index2,index3
      integer*4 findex1,findex2,findex3
      integer*4 ipxmin,ipymin,ipzmin
      integer*4 ipxmax,ipymax,ipzmax
      integer*4 trki_pid,trkj_pid
      integer*4 px_sec_save,py_sec_save,pz_sec_save

      real*4    qinvar2,qtotal2,qvector2
      real*4    qside2, qout2,  qlong2
      real*4    qperp2, qtime2
      real*4    qdiff1, qdiff2, qdiff3
      real*4    px_save,py_save,pz_save,E_save

      If (mode .eq. 1) Then     ! Full hist* filling; initialize arrays to zero.

        do i = 1,max_h_1d
          hist_like_1d(i)   = 0
          hist_unlike_1d(i) = 0
        end do

        do i = 1,max_h_3d
        do j = 1,max_h_3d
        do k = 1,max_h_3d
          hist_like_3d_fine(i,j,k)     = 0
          hist_unlike_3d_fine(i,j,k)   = 0
          hist_like_3d_coarse(i,j,k)   = 0
          hist_unlike_3d_coarse(i,j,k) = 0
        end do
        end do
        end do

      End If

CCC   Select # of particles to loop over for each 'mode':

      If (mode .eq. 1) Then
         imin = 2
         imax = n_part_tot_trk
      Else If (mode .eq. 2 .or. mode .eq. 3) Then
         imin = itrack
         imax = itrack
      Else If (mode .eq. 4) Then
         imin = 1
         imax = n_part_tot_trk
      End If

C------------------------------------------------------
CCC   Begin Primary Loop over particles in Table 'trk':
C------------------------------------------------------

      do i = imin,imax
      if(trk_flag(i) .eq. 0) then
        trki_pid = trk_ge_pid(i)
        if(mode.eq.2) then
          Call sec_index(px_sec_save,n_px_bins,ipxmin,ipxmax)
          Call sec_index(py_sec_save,n_py_bins,ipymin,ipymax)
          Call sec_index(pz_sec_save,n_pz_bins,ipzmin,ipzmax)
        else
          Call sec_index(trk_px_sec(i),n_px_bins,ipxmin,ipxmax)
          Call sec_index(trk_py_sec(i),n_py_bins,ipymin,ipymax)
          Call sec_index(trk_pz_sec(i),n_pz_bins,ipzmin,ipzmax)
        end if

CCC   Begin Loop over neighboring sectors for track # 'i':

        do jx = ipxmin,ipxmax
        do jy = ipymin,ipymax
        do jz = ipzmin,ipzmax
          trkj_sector = jx + (jy-1)*n_px_bins
     1                + (jz-1)*n_px_bins*n_py_bins
          njtrks = 0
          if(mode.le.3) njtrks = stm_n_trk_sec(trkj_sector)
          if(mode.eq.4) njtrks = stm2_n_trk_sec(trkj_sector)
          if(njtrks .gt. 0) then

CCC   Begin Secondary Loop over particles in selected sectors in tables
CCC   'trk' or 'trk2':
            do jsec = 1,njtrks
              if(mode.le.3) j = stm_track_id(jsec,trkj_sector)
              if(mode.eq.4) j = stm2_track_id(jsec,trkj_sector)
              if((mode.eq.1 .and. j.lt.i .and. trk_flag(j).eq.0)
     1       .or.(mode.eq.2 .and. j.ne.i .and. trk_flag(j).eq.0)
     2       .or.(mode.eq.3 .and. j.ne.i .and. trk_flag(j).eq.0)
     3       .or.(mode.eq.4 .and. trk2_flag(j).eq.0)) then

CCC   Obtain 1D and 3D relative momenta:

                if(mode.eq.1 .or. mode.eq.3) then
                  trkj_pid = trk_ge_pid(j)
                  Call qdiff(trk_px(i),trk_py(i),trk_pz(i),trk_E(i),
     1                       trk_px(j),trk_py(j),trk_pz(j),trk_E(j),
     2              qinvar2,qtotal2,qvector2,qside2,qout2,qlong2,
     3              qperp2,qtime2)
                else if(mode.eq.2) then
                  trkj_pid = trk_ge_pid(j)
                  Call qdiff(px_save,py_save,pz_save,E_save,
     1                       trk_px(j),trk_py(j),trk_pz(j),trk_E(j),
     2              qinvar2,qtotal2,qvector2,qside2,qout2,qlong2,
     3              qperp2,qtime2)
                else if(mode.eq.4) then
                  trkj_pid = trk2_ge_pid(j)
                  Call qdiff(trk_px(i),trk_py(i),trk_pz(i),trk_E(i),
     1              trk2_px(j),trk2_py(j),trk2_pz(j),trk2_E(j),
     2              qinvar2,qtotal2,qvector2,qside2,qout2,qlong2,
     3              qperp2,qtime2)
                end if

C-----------------------------------------------
CCC   Fill and/or Update 1D two-body Histograms:
C-----------------------------------------------

                if(switch_1d .gt. 0) then

                  if(switch_1d .eq. 1) then
                    qdiff1 = sqrt(qinvar2)
                  else if(switch_1d .eq. 2) then
                    qdiff1 = sqrt(qtotal2)
                  else if(switch_1d .eq. 3) then
                    qdiff1 = sqrt(qvector2)
                  else
                    qdiff1 = sqrt(qvector2)
                  end if

                  if(qdiff1 .le. qmid_1d) then
                    index1 = int((qdiff1/binsize_1d_fine)+0.99999)
                    if(index1 .eq. 0) index1 = 1
                  else if(qdiff1.gt.qmid_1d.and.qdiff1.le.qmax_1d) then
                    index1 = int(((qdiff1-qmid_1d)/binsize_1d_coarse)
     1                     + 0.99999)
                    if(index1.eq.0) index1 = 1
                    index1 = index1 + n_1d_fine
                  else
                    index1 = -86
                  end if

                  if(index1.ge.1.and.index1.le.n_1d_total) then
                    if((trki_pid.eq.trkj_pid).and.(switch_type.eq.1
     1                  .or. switch_type.eq.3)) then
                      if(mode.eq.1) then
                        hist_like_1d(index1) = hist_like_1d(index1) + 1
                      else if(mode.eq.2) then
                        htmp_like_1d(index1) = htmp_like_1d(index1) - 1
                      else if(mode.eq.3) then
                        htmp_like_1d(index1) = htmp_like_1d(index1) + 1
                      else if(mode.eq.4) then
                        href_like_1d(index1) = href_like_1d(index1) + 1
                      end if

                    else if((trki_pid.ne.trkj_pid).and.(switch_type.eq.2
     1                       .or. switch_type.eq.3)) then
                      if(mode.eq.1) then
                      hist_unlike_1d(index1) = hist_unlike_1d(index1)+1
                      else if(mode.eq.2) then
                      htmp_unlike_1d(index1) = htmp_unlike_1d(index1)-1
                      else if(mode.eq.3) then
                      htmp_unlike_1d(index1) = htmp_unlike_1d(index1)+1
                      else if(mode.eq.4) then
                      href_unlike_1d(index1) = href_unlike_1d(index1)+1
                      end if

                    end if
                  end if
                end if     ! End 1D Histogram Fill and/or Update.

C-----------------------------------------------
CCC   Fill and/or Update 3D Two-Body Histograms:
C-----------------------------------------------

                if(switch_3d .gt. 0) then
                  if(switch_3d .eq. 1) then
                    qdiff1 = sqrt(qside2)
                    qdiff2 = sqrt(qout2)
                    qdiff3 = sqrt(qlong2)
                  else if(switch_3d .eq. 2) then
                    qdiff1 = sqrt(qperp2)
                    qdiff2 = sqrt(qtime2)
                    qdiff3 = sqrt(qlong2)
                  else
                    qdiff1 = sqrt(qperp2)
                    qdiff2 = sqrt(qtime2)
                    qdiff3 = sqrt(qlong2)
                  end if

                  if(qdiff1 .le. qmid_3d) then
                    findex1 = int((qdiff1/binsize_3d_fine)+0.99999)
                    if(findex1 .eq. 0) findex1 = 1
                    index1 = 1
                  else if(qdiff1.gt.qmid_3d.and.qdiff1.le.qmax_3d) then
                    index1 = int((qdiff1/binsize_3d_coarse)+0.99999)
                    if(index1.eq.1) index1 = 2
                    findex1 = 0
                  else
                    index1 = -86
                    findex1 = -86
                  end if

                  if(qdiff2 .le. qmid_3d) then
                    findex2 = int((qdiff2/binsize_3d_fine)+0.99999)
                    if(findex2 .eq. 0) findex2 = 1
                    index2 = 1
                  else if(qdiff2.gt.qmid_3d.and.qdiff2.le.qmax_3d) then
                    index2 = int((qdiff2/binsize_3d_coarse)+0.99999)
                    if(index2.eq.1) index2 = 2
                    findex2 = 0
                  else
                    index2 = -86
                    findex2 = -86
                  end if

                  if(qdiff3 .le. qmid_3d) then
                    findex3 = int((qdiff3/binsize_3d_fine)+0.99999)
                    if(findex3 .eq. 0) findex3 = 1
                    index3 = 1
                  else if(qdiff3.gt.qmid_3d.and.qdiff3.le.qmax_3d) then
                    index3 = int((qdiff3/binsize_3d_coarse)+0.99999)
                    if(index3.eq.1) index3 = 2
                    findex3 = 0
                  else
                    index3 = -86
                    findex3 = -86
                  end if

                  if((index1.ge.1.and.index1.le.n_3d_coarse).and.
     1               (index2.ge.1.and.index2.le.n_3d_coarse).and.
     2               (index3.ge.1.and.index3.le.n_3d_coarse)) then

                   if((index1+index2+index3).eq.3) then

                    if(findex1.ge.1.and.findex1.le.n_3d_fine.and.
     1                 findex2.ge.1.and.findex2.le.n_3d_fine.and.
     2                 findex3.ge.1.and.findex3.le.n_3d_fine) then

                     if((trki_pid.eq.trkj_pid).and.(switch_type.eq.1
     1                   .or. switch_type.eq.3)) then

                      if(mode.eq.1) then
                       hist_like_3d_fine(findex1,findex2,findex3) =
     1                 hist_like_3d_fine(findex1,findex2,findex3) +1
                      else if(mode.eq.2) then
                       htmp_like_3d_fine(findex1,findex2,findex3) =
     1                 htmp_like_3d_fine(findex1,findex2,findex3) -1
                      else if(mode.eq.3) then
                       htmp_like_3d_fine(findex1,findex2,findex3) =
     1                 htmp_like_3d_fine(findex1,findex2,findex3) +1
                      else if(mode.eq.4) then
                       href_like_3d_fine(findex1,findex2,findex3) =
     1                 href_like_3d_fine(findex1,findex2,findex3) +1
                      end if

                     else if((trki_pid.ne.trkj_pid).and.(switch_type
     1                    .eq.2 .or. switch_type.eq.3)) then

                      if(mode.eq.1) then
                       hist_unlike_3d_fine(findex1,findex2,findex3) =
     1                 hist_unlike_3d_fine(findex1,findex2,findex3) +1
                      else if(mode.eq.2) then
                       htmp_unlike_3d_fine(findex1,findex2,findex3) =
     1                 htmp_unlike_3d_fine(findex1,findex2,findex3) -1
                      else if(mode.eq.3) then
                       htmp_unlike_3d_fine(findex1,findex2,findex3) =
     1                 htmp_unlike_3d_fine(findex1,findex2,findex3) +1
                      else if(mode.eq.4) then
                       href_unlike_3d_fine(findex1,findex2,findex3) =
     1                 href_unlike_3d_fine(findex1,findex2,findex3) +1
                      end if

                     end if
                    end if

                   else if((index1+index2+index3).gt.3) then

                     if((trki_pid.eq.trkj_pid).and.(switch_type.eq.1
     1                   .or. switch_type.eq.3)) then

                      if(mode.eq.1) then
                       hist_like_3d_coarse(index1,index2,index3) =
     1                 hist_like_3d_coarse(index1,index2,index3) +1
                      else if(mode.eq.2) then
                       htmp_like_3d_coarse(index1,index2,index3) =
     1                 htmp_like_3d_coarse(index1,index2,index3) -1
                      else if(mode.eq.3) then
                       htmp_like_3d_coarse(index1,index2,index3) =
     1                 htmp_like_3d_coarse(index1,index2,index3) +1
                      else if(mode.eq.4) then
                       href_like_3d_coarse(index1,index2,index3) =
     1                 href_like_3d_coarse(index1,index2,index3) +1
                      end if

                     else if((trki_pid.ne.trkj_pid).and.(switch_type
     1                    .eq.2 .or. switch_type.eq.3)) then

                      if(mode.eq.1) then
                       hist_unlike_3d_coarse(index1,index2,index3) =
     1                 hist_unlike_3d_coarse(index1,index2,index3) +1
                      else if(mode.eq.2) then
                       htmp_unlike_3d_coarse(index1,index2,index3) =
     1                 htmp_unlike_3d_coarse(index1,index2,index3) -1
                      else if(mode.eq.3) then
                       htmp_unlike_3d_coarse(index1,index2,index3) =
     1                 htmp_unlike_3d_coarse(index1,index2,index3) +1
                      else if(mode.eq.4) then
                       href_unlike_3d_coarse(index1,index2,index3) =
     1                 href_unlike_3d_coarse(index1,index2,index3) +1
                      end if

                     end if
                    end if    ! End 3D Fine/Coarse Grid
                   end if
                  end if      ! End 3D Histogram Filling and/or Update

               end if
               end do         ! End Secondary Track Loop

            end if
            end do
            end do
            end do            ! End Neighboring Sector Loop

         end if
      end do                  ! End Primary Track Loop

      Return
      END

C-----------------------------------------------------------------------


      subroutine correlation_fit
      implicit none

CCC   This subroutine carries out the track momentum adjustment
CCC   procedure in order to fit the requested model correlation
CCC   function and the input one-body {pt,phi,eta} distributions.

CCC   The method used is similar to the Metropolis method.  Briefly:
C
C     1. The accepted tracks for each event in the 'event_text.in'
C        input file are loaded into the 'trk' data structure table.
C     2. The sector information, histograms, C2 and initial chi-square
C        are computed.
C     3. Each track momentum is randomly shifted within a specified
C        range, one track at a time, the histograms are updated, and
C        a new C2 and chi-square are computed.
C     4. If the new track momentum is acceptable and if it results in a
C        smaller value of chi-square, then this shifted momentum is
C        retained, if not then the original momentum value is restored.
C     5. This is done for all particles in the track table for the event.
C     6. The entire process is repeated either until the maximum # of
C        iterations is reached, or until the chi-square improvement with
C        each iteration diminishes sufficiently.
C     7. After completing the event loop, summary information is gathered
C        and inclusive correlation functions and one-body distributions
C        are calculated.

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'
      Include 'common_event_summary.inc'

      Include 'common_track.inc'
      Include 'common_track2.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'
      Include 'common_particle.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,j,k,ievent,niter,ntracks_out,nev
      integer*4 ntracks_flagged,track_status,pid_index
      integer*4 px_sec_save,py_sec_save,pz_sec_save

      real*4    px_save,py_save,pz_save,E_save
      real*4    pt_save,phi_save,eta_save,mass
      real*4    chisq_like_1d,chisq_unlike_1d
      real*4    chisq_like_3d_fine,chisq_unlike_3d_fine
      real*4    chisq_like_3d_coarse,chisq_unlike_3d_coarse
      real*4    chisq_hist1_1,chisq_hist1_2
      real*4    chisq_total,chisq_total_oldvalue,chisq_total_newvalue
      real*4    hbtpran

CCC   Initialize counters:

      n_part_used_1_inc    = 0
      n_part_used_2_inc    = 0
      num_pairs_like_inc   = 0
      num_pairs_unlike_inc = 0
      event_line_counter   = 0
      file10_line_counter  = 0

CCC   Open event input, track selection flags and event output files:

      If(ALICE .eq. 0) Then
      open(unit=2,status='old',access='sequential',
     1     file='event_text.in')
      open(unit=4,status='old',access='sequential',
     1     file='event_tracks.select')
      open(unit=10,status='unknown',access='sequential',
     1     file='event_hbt_text.out')
CCC   Read/Write event header from/to I/O event text files
         Call read_data(7)
         Call read_data(8)
      End If

C-------------------------------------
C  Begin Event Loop, 
C
      do ievent = 2, n_events + 1   
C-------------------------------------

        If(ALICE .eq. 1) Then
           Call AliHbtp_SetActiveEventNumber(ievent-1)
C           write(*,*) 'NEXT EVENT:', ievent
        End If
        Call read_data(7)
        if(n_part_tot_trk .gt. 0) then

          write(6,98)
          Call tindex(1,0)  !  Fill initial track-sector info.
          Call stm_build(1,0,0)  ! Fill initial sector info.
          Call dist_range(1,ntracks_out,ntracks_flagged)
          num_pairs_like = (n_part_used_1_trk*(n_part_used_1_trk-1))/2
     1                   + (n_part_used_2_trk*(n_part_used_2_trk-1))/2
          num_pairs_unlike = n_part_used_1_trk*n_part_used_2_trk
          num_pairs_like_inc = num_pairs_like_inc + num_pairs_like
          num_pairs_unlike_inc = num_pairs_unlike_inc + num_pairs_unlike
          n_part_used_1_inc = n_part_used_1_inc + n_part_used_1_trk
          n_part_used_2_inc = n_part_used_2_inc + n_part_used_2_trk
          
C          write (*,*) 'num_pairs_like = ',num_pairs_like
C          write (*,*) 'num_pairs_unlike = ',num_pairs_unlike
C          write (*,*) 'num_pairs_like_inc = ',num_pairs_like_inc
c          write (*,*) 'num_pairs_unlike_inc = ',num_pairs_unlike_inc
c          write (*,*) 'n_part_used_1_inc = ',n_part_used_1_inc
C          write (*,*) 'n_part_used_2_inc = ',n_part_used_2_inc
          
          if(pid(1).gt.0) Call histog1(1,0,1,pid(1),0.,0.,0.)
          if(pid(2).gt.0) Call histog1(1,0,2,pid(2),0.,0.,0.)
          Call histog2(1,0,0,0,0,0.0,0.0,0.0,0.0)
          Call correl_fit(1)
          Call chisquare(1,chisq_like_1d,chisq_unlike_1d,
     1                   chisq_like_3d_fine,chisq_unlike_3d_fine,
     2                   chisq_like_3d_coarse,chisq_unlike_3d_coarse,
     3                   chisq_hist1_1,chisq_hist1_2)
          chisq_total = chisq_wt_like_1d    *chisq_like_1d
     1                + chisq_wt_unlike_1d  *chisq_unlike_1d
     2         + chisq_wt_like_3d_fine      *chisq_like_3d_fine
     3         + chisq_wt_unlike_3d_fine    *chisq_unlike_3d_fine
     4         + chisq_wt_like_3d_coarse    *chisq_like_3d_coarse
     5         + chisq_wt_unlike_3d_coarse  *chisq_unlike_3d_coarse
     6         + chisq_wt_hist1_1           *chisq_hist1_1
     7         + chisq_wt_hist1_2           *chisq_hist1_2
          chisq_total_oldvalue = chisq_total
          Call hist1_copy(1)
          Call hist2_copy(1)

          niter = 0
1000      Continue       ! Starting Point for Track Shift Iteration Loop:
          niter = niter + 1

          if(niter.eq.1) then
             write(8,99)
             write(8,98)
             write(8,99)
          end if
98        Format(5x,'**  START NEXT EVENT  **')
99        Format(5x,'************************')
          write(8,100) ievent,niter,chisq_total
100       Format(10x,'Event#+1,Iteration# and Chi-Sq = ',2I5,E11.4)

          IF(maxit .eq. 0) GO TO 1001   ! Option to compute correlations
C                                       ! for input events.

C-------------------------------------
C  Begin Track Adjustment Loop:

          do i = 1,n_part_tot_trk   
C-------------------------------------

          if(trk_flag(i) .eq. 0) then

CCC   Save initial track parameters (those that might change):

          px_save      = trk_px(i)
          py_save      = trk_py(i)
          pz_save      = trk_pz(i)
          E_save       = trk_E(i)
          pt_save      = trk_pt(i)
          phi_save     = trk_phi(i)
          eta_save     = trk_eta(i)
          px_sec_save  = trk_px_sec(i)
          py_sec_save  = trk_py_sec(i)
          pz_sec_save  = trk_pz_sec(i)
          old_sec_save = trk_sector(i)

CCC   Save the sector values that the track is initially located in:

          old_sec_ntrk = stm_n_trk_sec(trk_sector(i))
          old_sec_flag = stm_flag(trk_sector(i))
          do k = 1,stm_n_trk_sec(trk_sector(i))
            old_sec_trkid(k) = stm_track_id(k,trk_sector(i))
          end do

CCC   Determine the particle ID type:

          if(trk_ge_pid(i).eq.pid(1) .and. pid(1).gt.0) then
            pid_index = 1
          else if(trk_ge_pid(i).eq.pid(2).and.pid(2).gt.0) then
            pid_index = 2
          else
            pid_index = 1
          end if

CCC   Randomly shift track momentum vector and compute new kinematics:

          trk_px(i) = trk_px(i) + deltap*(2.0*hbtpran(irand) - 1.0)
          trk_py(i) = trk_py(i) + deltap*(2.0*hbtpran(irand) - 1.0)
          trk_pz(i) = trk_pz(i) + deltap*(2.0*hbtpran(irand) - 1.0)
          mass = part_mass(trk_ge_pid(i))
          Call Hbtp_kin(trk_px(i),trk_py(i),trk_pz(i),trk_E(i),
     1                    trk_pt(i),trk_phi(i),trk_eta(i),mass,1)
          Call tindex(1,i)
          new_sec_save = trk_sector(i)

CCC   Determine if track has been shifted to a new sector, and if so,
CCC   whether this overfills this new sector.  If all is well, then
CCC   update histograms.  If not, then restore track parameters to their
CCC   initial values prior to shifting.  Keep the status of track(i) in
CCC   'track_status', where a value of 0 means the track is OK to use.
CCC
CCC   The Logical steps are the following:
CCC
C     IF(new track position is in same sector) THEN
C       o Remove old track position from htmp1*, htmp* using old saved values.
C       o Add new track position to htmp1*, htmp* using values in 'trk'
C         (Sector information is unchanged)
C     ELSE IF(new track position is in a different sector) THEN
C       IF(# tracks in new sector is still OK, with the new track) THEN 
C         o Save values of new sector before trk was shifted into it.
C         o Remove old trk position from htmp1*, htmp* using old saved values
C         o Add new trk position to htmp1*, htmp* using values in trk
C         o Update sector information in stm
C       ELSE IF(# tracks in new sector becomes too many with new trk) THEN
C         o Restore track parameters to pre-shifted values
C         o Set track_status = 1, indicating the track could not be moved
C       END IF
C     END IF

          track_status = 0
          if(old_sec_save .eq. new_sec_save) then
            Call histog1(2,i,pid_index,pid(pid_index),pt_save,
     1                   phi_save,eta_save)
            Call histog2(2,i,px_sec_save,py_sec_save,pz_sec_save,
     1           px_save,py_save,pz_save,E_save)

            Call histog1(3,i,pid_index,pid(pid_index),0.,0.,0.)
            Call histog2(3,i,0,0,0,0.0,0.0,0.0,0.0)

          else if(old_sec_save .ne. new_sec_save) then

            if(stm_n_trk_sec(new_sec_save) .lt. max_trk_sec) then
              new_sec_ntrk = stm_n_trk_sec(new_sec_save)
              new_sec_flag = stm_flag(new_sec_save)
              if(new_sec_ntrk .gt. 0) then
                do k = 1,new_sec_ntrk
                  new_sec_trkid(k) = stm_track_id(k,new_sec_save)
                end do
              end if

              Call histog1(2,i,pid_index,pid(pid_index),pt_save,
     1                   phi_save,eta_save)
              Call histog2(2,i,px_sec_save,py_sec_save,pz_sec_save,
     1           px_save,py_save,pz_save,E_save)

              Call histog1(3,i,pid_index,pid(pid_index),
     1           0.,0.,0.)
              Call histog2(3,i,0,0,0,0.0,0.0,0.0,0.0)

              Call stm_build(1,i,old_sec_save)
           
            else if(stm_n_trk_sec(new_sec_save) .ge. max_trk_sec) then
     
              track_status = 1
              trk_px(i)      = px_save
              trk_py(i)      = py_save
              trk_pz(i)      = pz_save
              trk_E(i)       = E_save
              trk_pt(i)      = pt_save
              trk_phi(i)     = phi_save
              trk_eta(i)     = eta_save
              trk_px_sec(i)  = px_sec_save
              trk_py_sec(i)  = py_sec_save
              trk_pz_sec(i)  = pz_sec_save
              trk_sector(i)  = old_sec_save

            end if
          end if        ! End Histogram and Sector Update

CCC   If the track was succesfully shifted then compute C2 and determine
C     if the chi-square decreases (improves) or increases.  If it improves,
C     then store the new chi-square value and keep the 1- and 2-body
C     histograms in hist1* and hist*, repsectively.  If chi-square
C     increases (worsens), then restore the track parameters to the
C     pre-shifted values, restore the histograms and if a new sector was
C     involved, then restore both the old and new sector values.
C
C     The Logical steps are the following:
C
C     IF(new track position is OK, (i.e. track_status = 0)) Then
C       o Compute C2 using htmp*
C       o Compute chi-square and save
C       IF(chi-square improves) Then
C         o Replace previous (best) chi-square with new value
C         o Update histograms, i.e. copy htmp1* -> hist1* and
C                                   copy htmp*  -> hist*
C       ELSE IF(chi-square increases) Then
C         o Restore track parameters
C         o Restore histograms, i.e. copy hist1* -> htmp1* and
C                                    copy hist*  -> htmp*
C         IF(new sector was used) Then
C           o Restore old sector values to pre-shifted values
C           o Restore new sector values to pre-shifted values
C         END IF
C       END IF
C     END IF

          If(track_status .eq.0) Then
            Call correl_fit(2)                   
            Call chisquare(2,chisq_like_1d,chisq_unlike_1d,
     1                     chisq_like_3d_fine,chisq_unlike_3d_fine,
     2                     chisq_like_3d_coarse,chisq_unlike_3d_coarse,
     3                     chisq_hist1_1,chisq_hist1_2)
            chisq_total_newvalue =
     1         chisq_wt_like_1d           *chisq_like_1d
     2       + chisq_wt_unlike_1d         *chisq_unlike_1d
     3       + chisq_wt_like_3d_fine      *chisq_like_3d_fine
     4       + chisq_wt_unlike_3d_fine    *chisq_unlike_3d_fine
     5       + chisq_wt_like_3d_coarse    *chisq_like_3d_coarse
     6       + chisq_wt_unlike_3d_coarse  *chisq_unlike_3d_coarse
     7       + chisq_wt_hist1_1           *chisq_hist1_1
     8       + chisq_wt_hist1_2           *chisq_hist1_2

            if(chisq_total_newvalue .lt. chisq_total_oldvalue) then
              chisq_total_oldvalue = chisq_total_newvalue
              Call hist1_copy(2)
              Call hist2_copy(2)
            else if(chisq_total_newvalue.ge.chisq_total_oldvalue) then
              trk_px(i)      = px_save
              trk_py(i)      = py_save
              trk_pz(i)      = pz_save
              trk_E(i)       = E_save
              trk_pt(i)      = pt_save
              trk_phi(i)     = phi_save
              trk_eta(i)     = eta_save
              trk_px_sec(i)  = px_sec_save
              trk_py_sec(i)  = py_sec_save
              trk_pz_sec(i)  = pz_sec_save
              trk_sector(i)  = old_sec_save
              Call hist1_copy(1)
              Call hist2_copy(1)

              If(old_sec_save .ne. new_sec_save) then

                stm_n_trk_sec(old_sec_save) = old_sec_ntrk
                stm_flag(old_sec_save)      = old_sec_flag
                do k = 1,max_trk_sec
                  stm_track_id(k,old_sec_save) = 0
                end do
                do k = 1,old_sec_ntrk
                  stm_track_id(k,old_sec_save) = old_sec_trkid(k)
                end do

                stm_n_trk_sec(new_sec_save) = new_sec_ntrk
                stm_flag(new_sec_save)      = new_sec_flag
                do k = 1,max_trk_sec
                  stm_track_id(k,new_sec_save) = 0
                end do
                do k = 1,new_sec_ntrk
                  stm_track_id(k,new_sec_save) = new_sec_trkid(k)
                end do

              end if
            end if
          end if         ! End Chi-Square Determination
        end if           ! End valid tracks option
        end do           ! End Track Shift Loop

CCC   Check chi-square for this iteration --
C        Best, current chi-square value is in 'chisq_total_oldvalue'
C        Chi-square value at the beginning of the iteration loop is in
C        'chisq_total'.

        If(abs(200.0*(chisq_total_oldvalue - chisq_total)/
     1    (chisq_total_oldvalue + chisq_total)) .lt. delchi) then
          write(8,101)
101       Format(/5x,'Chi-Sq reduced .lt. delchi % on last iteration',
     1           ' - Stop Search')
          go to 1001
        End If
        If (niter .gt. maxit) Then
          write(8,102)
102       Format(/5x,'Max # Search Iterations Reached - Abort track ',
     1           'Adj. process')
          go to 1001
        End If
        chisq_total = chisq_total_oldvalue
        go to 1000

1001    Continue

CCC   Finished Track Adjustment Iteration Loop for event # 'ievent'

        if((ievent - 1) .le. max_events) then
          Call dist_range(1,ntracks_out,ntracks_flagged)
          num_iter(ievent-1) = float(niter)
          n_part_used_1_store(ievent-1) = float(n_part_used_1_trk)
          n_part_used_2_store(ievent-1) = float(n_part_used_2_trk)
          n_part_tot_store(ievent-1) = float(n_part_tot_trk)
          frac_trks_out(ievent-1)=float(ntracks_out)/
     1                            float(n_part_tot_trk)
          frac_trks_flag(ievent-1) = 
     1         float(ntracks_flagged)/float(n_part_tot_trk)
        end if
       
        Call correl_fit(1)
        Call chisquare(1,chisq_like_1d,chisq_unlike_1d,
     1                 chisq_like_3d_fine,chisq_unlike_3d_fine,
     2                 chisq_like_3d_coarse,chisq_unlike_3d_coarse,
     3                 chisq_hist1_1,chisq_hist1_2)
        chisq_total = chisq_wt_like_1d    *chisq_like_1d
     1              + chisq_wt_unlike_1d  *chisq_unlike_1d
     2       + chisq_wt_like_3d_fine      *chisq_like_3d_fine
     3       + chisq_wt_unlike_3d_fine    *chisq_unlike_3d_fine
     4       + chisq_wt_like_3d_coarse    *chisq_like_3d_coarse
     5       + chisq_wt_unlike_3d_coarse  *chisq_unlike_3d_coarse
     6       + chisq_wt_hist1_1           *chisq_hist1_1
     7       + chisq_wt_hist1_2           *chisq_hist1_2

        if((ievent - 1) .le. max_events) then
          chisq_like_1d_store(ievent-1)          = chisq_like_1d
          chisq_unlike_1d_store(ievent-1)        = chisq_unlike_1d
          chisq_like_3d_fine_store(ievent-1)     = chisq_like_3d_fine
          chisq_unlike_3d_fine_store(ievent-1)   = chisq_unlike_3d_fine
          chisq_like_3d_coarse_store(ievent-1)   = chisq_like_3d_coarse
          chisq_unlike_3d_coarse_store(ievent-1) =chisq_unlike_3d_coarse
          chisq_hist1_1_store(ievent-1)          = chisq_hist1_1
          chisq_hist1_2_store(ievent-1)          = chisq_hist1_2
          chisq_total_store(ievent-1)            = chisq_total

CCC   Count # sectors with stm().flag = 1, indicating that too many 
C     tracks were attempted to be loaded into that sector.

          num_sec_flagged_store(ievent-1) = 0.0
          do k = 1,n_sectors
            if(stm_flag(k) .eq. 1) then
              num_sec_flagged_store(ievent-1) = 
     1             num_sec_flagged_store(ievent-1) + 1.0
            end if
          end do
        end if

        Call hist1_incl_sum
        Call hist2_incl_sum
        if(print_full .eq. 1) Call write_data(5,ievent-1)

      end if    !  End event-with-tracks processing.
        Call read_data(8)

C-------------------------------
      end do    ! End Event Loop
C-------------------------------

CCC   Compute Correlation Functions for the Inclusive Histograms

      Call correl_fit(3)

CCC   Compute Mean and Std. dev of event monitor and summary quantities:

      if(n_events .le. max_events) then
         nev = n_events
      else
         nev = max_events
      end if

      Call mean_rms(num_iter,nev,nev,niter_mean,niter_rms)
      Call mean_rms(n_part_used_1_store,nev,nev,npart1_mean,npart1_rms)
      Call mean_rms(n_part_used_2_store,nev,nev,npart2_mean,npart2_rms)
      Call mean_rms(n_part_tot_store,nev,nev,npart_tot_mean,
     1              npart_tot_rms)
      Call mean_rms(num_sec_flagged_store,nev,nev,
     1              nsec_flag_mean,nsec_flag_rms)
      Call mean_rms(frac_trks_out,nev,nev,
     1              frac_trks_out_mean,frac_trks_out_rms)
      Call mean_rms(frac_trks_flag,nev,nev,
     1              frac_trks_flag_mean,frac_trks_flag_rms)
      Call mean_rms(chisq_like_1d_store,nev,nev,
     1              chi_l1d_mean,chi_l1d_rms)
      Call mean_rms(chisq_unlike_1d_store,nev,nev,
     1              chi_u1d_mean,chi_u1d_rms)
      Call mean_rms(chisq_like_3d_fine_store,nev,nev,
     1              chi_l3f_mean,chi_l3f_rms)
      Call mean_rms(chisq_unlike_3d_fine_store,nev,nev,
     1              chi_u3f_mean,chi_u3f_rms)
      Call mean_rms(chisq_like_3d_coarse_store,nev,nev,
     1              chi_l3c_mean,chi_l3c_rms)
      Call mean_rms(chisq_unlike_3d_coarse_store,nev,nev,
     1              chi_u3c_mean,chi_u3c_rms)
      Call mean_rms(chisq_hist1_1_store,nev,nev,
     1              chi_1_1_mean, chi_1_1_rms)
      Call mean_rms(chisq_hist1_2_store,nev,nev,
     1              chi_1_2_mean, chi_1_2_rms)
      Call mean_rms(chisq_total_store,nev,nev,
     1              chi_tot_mean, chi_tot_rms)

      If(ALICE .eq. 0) Then
      Close(unit=2)
      Close(unit=4)
      Close(unit=10)
      End If

      Return
      END

C------------------------------------------------------------------------


      subroutine hist1_copy(mode)
      implicit none

CCC   Copy 1-body histograms if:
CCC
CCC   mode = 1, then copy hist1* -> htmp1*
CCC   mode = 2, then copy htmp1* -> hist1*

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'

CCC   Local Variable Type Declarations:

      integer*4  mode, i

C---------------------------
      If(mode .eq. 1) Then       !  Copy hist1* -> htmp1*
C---------------------------

        if(pid(1) .gt. 0) then
          do i = 1,n_pt_bins
            htmp1_pt_1(i)  = hist1_pt_1(i)
          end do
          do i = 1,n_phi_bins
            htmp1_phi_1(i) = hist1_phi_1(i)
          end do
          do i = 1,n_eta_bins
            htmp1_eta_1(i) = hist1_eta_1(i)
          end do
        end if

        if(pid(2) .gt. 0) then
          do i = 1,n_pt_bins     
            htmp1_pt_2(i)  = hist1_pt_2(i)
          end do
          do i = 1,n_phi_bins
            htmp1_phi_2(i) = hist1_phi_2(i)
          end do
          do i = 1,n_eta_bins
            htmp1_eta_2(i) = hist1_eta_2(i)
          end do
        end if
     
C--------------------------------
      Else If (mode .eq. 2) Then   ! Copy htmp1* -> hist1*
C--------------------------------

        if(pid(1) .gt. 0) then
          do i = 1,n_pt_bins     
            hist1_pt_1(i)  = htmp1_pt_1(i)
          end do
          do i = 1,n_phi_bins
            hist1_phi_1(i) = htmp1_phi_1(i)
          end do
          do i = 1,n_eta_bins
            hist1_eta_1(i) = htmp1_eta_1(i)
          end do
        end if

        if(pid(2) .gt. 0) then
          do i = 1,n_pt_bins     
            hist1_pt_2(i)  = htmp1_pt_2(i)
          end do
          do i = 1,n_phi_bins
            hist1_phi_2(i) = htmp1_phi_2(i)
          end do
          do i = 1,n_eta_bins
            hist1_eta_2(i) = htmp1_eta_2(i)
          end do
        end if 

C------------
      End If
C------------

      Return
      END

C----------------------------------------------------------------------

      subroutine hist2_copy(mode)
      implicit none

CCC   Copy 2-body histograms if:
CCC
CCC   mode = 1, then copy hist* -> htmp*
CCC   mode = 2, then copy htmp* -> hist*

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'

CCC   Local Variable Type Declarations:

      integer*4  mode, i,j,k

C---------------------------
      If (mode .eq. 1) Then    !   Copy hist* -> htmp*
C---------------------------

        if(switch_1d.gt.0 .and. n_1d_total.gt.0) then
          if(switch_type.eq.1 .or. switch_type.eq.3) then
            do i = 1,n_1d_total
              htmp_like_1d(i) = hist_like_1d(i)
            end do
          end if
          if(switch_type.eq.2 .or. switch_type.eq.3) then
            do i = 1,n_1d_total
              htmp_unlike_1d(i) = hist_unlike_1d(i)
            end do
          end if
        end if     ! End 1D histogram copy

        if(switch_3d.gt.0) then
          if(switch_type.eq.1 .or. switch_type.eq.3) then

            if(n_3d_fine .gt. 0) then
              do i = 1,n_3d_fine
              do j = 1,n_3d_fine
              do k = 1,n_3d_fine
                htmp_like_3d_fine(i,j,k) = hist_like_3d_fine(i,j,k)
              end do
              end do
              end do
            end if

            if(n_3d_coarse .gt. 0) then
              do i = 1,n_3d_coarse
              do j = 1,n_3d_coarse
              do k = 1,n_3d_coarse
                htmp_like_3d_coarse(i,j,k) = hist_like_3d_coarse(i,j,k)
              end do
              end do
              end do
            end if

          end if

          if(switch_type.eq.2 .or. switch_type.eq.3) then

            if(n_3d_fine .gt. 0) then
              do i = 1,n_3d_fine
              do j = 1,n_3d_fine
              do k = 1,n_3d_fine
                htmp_unlike_3d_fine(i,j,k) = hist_unlike_3d_fine(i,j,k)
              end do
              end do
              end do
            end if

            if(n_3d_coarse .gt. 0) then
              do i = 1,n_3d_coarse
              do j = 1,n_3d_coarse
              do k = 1,n_3d_coarse
               htmp_unlike_3d_coarse(i,j,k)=hist_unlike_3d_coarse(i,j,k)
              end do
              end do
              end do
            end if

          end if
        end if     ! End 3D histogram copy

C--------------------------------
      Else If (mode .eq. 2) Then   !   Copy htmp* -> hist*
C--------------------------------

        if(switch_1d.gt.0 .and. n_1d_total.gt.0) then
          if(switch_type.eq.1 .or. switch_type.eq.3) then
            do i = 1,n_1d_total
              hist_like_1d(i) = htmp_like_1d(i)
            end do
          end if
          if(switch_type.eq.2 .or. switch_type.eq.3) then
            do i = 1,n_1d_total
              hist_unlike_1d(i) = htmp_unlike_1d(i)
            end do
          end if
        end if     ! End 1D histogram copy

        if(switch_3d.gt.0) then
          if(switch_type.eq.1 .or. switch_type.eq.3) then

            if(n_3d_fine .gt. 0) then
              do i = 1,n_3d_fine
              do j = 1,n_3d_fine
              do k = 1,n_3d_fine
                hist_like_3d_fine(i,j,k) = htmp_like_3d_fine(i,j,k)
              end do
              end do
              end do
            end if

            if(n_3d_coarse .gt. 0) then
              do i = 1,n_3d_coarse
              do j = 1,n_3d_coarse
              do k = 1,n_3d_coarse
                hist_like_3d_coarse(i,j,k) = htmp_like_3d_coarse(i,j,k)
              end do
              end do
              end do
            end if

          end if

          if(switch_type.eq.2 .or. switch_type.eq.3) then

            if(n_3d_fine .gt. 0) then
              do i = 1,n_3d_fine
              do j = 1,n_3d_fine
              do k = 1,n_3d_fine
                hist_unlike_3d_fine(i,j,k) = htmp_unlike_3d_fine(i,j,k)
              end do
              end do
              end do
            end if

            if(n_3d_coarse .gt. 0) then
              do i = 1,n_3d_coarse
              do j = 1,n_3d_coarse
              do k = 1,n_3d_coarse
               hist_unlike_3d_coarse(i,j,k)=htmp_unlike_3d_coarse(i,j,k)
              end do
              end do
              end do
            end if

          end if
        end if     ! End 3D histogram copy

C-------------
      End If       ! End mode selection options
C-------------

      Return
      END

C-----------------------------------------------------------------------


      subroutine hist1_incl_sum
      implicit none

CCC   Sum 1-body histograms for each event into inclusive totals, where
CCC   hinc1* = SUM[hist1*]

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'

CCC   Local Variable Type Declarations:

      integer*4 i

      if(pid(1) .gt. 0) then
        do i = 1,n_pt_bins
          hinc1_pt_1(i) = hinc1_pt_1(i) + hist1_pt_1(i)
        end do
        do i = 1,n_phi_bins
          hinc1_phi_1(i) = hinc1_phi_1(i) + hist1_phi_1(i)
        end do
        do i = 1,n_eta_bins
          hinc1_eta_1(i) = hinc1_eta_1(i) + hist1_eta_1(i)
        end do
      end if

      if(pid(2) .gt. 0) then
        do i = 1,n_pt_bins
          hinc1_pt_2(i) = hinc1_pt_2(i) + hist1_pt_2(i)
        end do
        do i = 1,n_phi_bins
          hinc1_phi_2(i) = hinc1_phi_2(i) + hist1_phi_2(i)
        end do
        do i = 1,n_eta_bins
          hinc1_eta_2(i) = hinc1_eta_2(i) + hist1_eta_2(i)
        end do
      end if

      Return
      END


C------------------------------------------------------------------------


      subroutine hist2_incl_sum
      implicit none

CCC   Sum 2-body histograms for each event into inclusive totals, where
CCC   hinc* = SUM[hist*]

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,j,k

      if(switch_1d.gt.0 .and. n_1d_total.gt.0) then
        if(switch_type.eq.1 .or. switch_type.eq.3) then
          do i = 1,n_1d_total
            hinc_like_1d(i) = hinc_like_1d(i) + hist_like_1d(i)
          end do
        end if
        if(switch_type.eq.2 .or. switch_type.eq.3) then
          do i = 1,n_1d_total
            hinc_unlike_1d(i) = hinc_unlike_1d(i) + hist_unlike_1d(i)
          end do
        end if
      end if      ! End 1D Inclusive Histogram Sum

      if(switch_3d.gt.0) then
        if(switch_type.eq.1 .or. switch_type.eq.3) then

          if(n_3d_fine .gt. 0) then
            do i = 1,n_3d_fine
            do j = 1,n_3d_fine
            do k = 1,n_3d_fine
              hinc_like_3d_fine(i,j,k) = hinc_like_3d_fine(i,j,k)
     1                                 + hist_like_3d_fine(i,j,k)
            end do
            end do
            end do
          end if

          if(n_3d_coarse .gt. 0) then
            do i = 1,n_3d_coarse
            do j = 1,n_3d_coarse
            do k = 1,n_3d_coarse
              hinc_like_3d_coarse(i,j,k) = hinc_like_3d_coarse(i,j,k)
     1                                   + hist_like_3d_coarse(i,j,k)
            end do
            end do
            end do
          end if

        end if

        if(switch_type.eq.2 .or. switch_type.eq.3) then

          if(n_3d_fine .gt. 0) then
            do i = 1,n_3d_fine
            do j = 1,n_3d_fine
            do k = 1,n_3d_fine
              hinc_unlike_3d_fine(i,j,k) = hinc_unlike_3d_fine(i,j,k)
     1                                   + hist_unlike_3d_fine(i,j,k)
            end do
            end do
            end do
          end if

          if(n_3d_coarse .gt. 0) then
            do i = 1,n_3d_coarse
            do j = 1,n_3d_coarse
            do k = 1,n_3d_coarse
             hinc_unlike_3d_coarse(i,j,k) = hinc_unlike_3d_coarse(i,j,k)
     1                                    + hist_unlike_3d_coarse(i,j,k)
            end do
            end do
            end do
          end if

        end if
      end if      ! End 3D Inclusive Histogram Sum

      Return
      END

C--------------------------------------------------------------------------


      subroutine correl_fit(mode)
      implicit none

CCC   This subroutine calculates the 2-body correlation function with
CCC   errors for the cases:
CCC
CCC      (1) 1D and/or 3D fine and coarse grid distributions
CCC      (2) like pairs and/or unlike pairs
CCC
CCC   Uses the signal and reference histograms.  The input parameter
CCC   'mode' selects which histograms to use.
CCC
CCC      Mode = 1, use hist*
CCC      Mode = 2, use htmp*
CCC      Mode = 3, use hinc*

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'

CCC   Local Variable Type Declarations:

      integer*4 mode,i,j,k

CCC   Initialize correlation functions and error arrays to zero: 

      do i = 1,max_c2_1d
         c2fit_like_1d(i)   = 0.0
         c2fit_unlike_1d(i) = 0.0
         c2err_like_1d(i)   = 0.0
         c2err_unlike_1d(i) = 0.0
      end do

      do i = 1,max_c2_3d
      do j = 1,max_c2_3d
      do k = 1,max_c2_3d
         c2fit_like_3d_fine(i,j,k)     = 0.0
         c2fit_unlike_3d_fine(i,j,k)   = 0.0
         c2fit_like_3d_coarse(i,j,k)   = 0.0
         c2fit_unlike_3d_coarse(i,j,k) = 0.0
         c2err_like_3d_fine(i,j,k)     = 0.0
         c2err_unlike_3d_fine(i,j,k)   = 0.0
         c2err_like_3d_coarse(i,j,k)   = 0.0
         c2err_unlike_3d_coarse(i,j,k) = 0.0
      end do
      end do
      end do

CCC   Compute 1D Correlation Functions and Errors:

      if(switch_1d .gt. 0) then
        if(switch_type.eq.1 .or. switch_type.eq.3) then

          if(mode .eq. 1) then
            Call c2_1d(hist_like_1d,href_like_1d,c2fit_like_1d,
     1           c2err_like_1d,max_h_1d,max_c2_1d,n_1d_total,
     2           num_pairs_like,num_pairs_like_ref)
          else if (mode .eq. 2) then
            Call c2_1d(htmp_like_1d,href_like_1d,c2fit_like_1d,
     1           c2err_like_1d,max_h_1d,max_c2_1d,n_1d_total,
     2           num_pairs_like,num_pairs_like_ref)
          else if (mode .eq. 3) then
            Call c2_1d(hinc_like_1d,href_like_1d,c2fit_like_1d,
     1           c2err_like_1d,max_h_1d,max_c2_1d,n_1d_total,
     2           num_pairs_like_inc,num_pairs_like_ref)
          end if

        end if

        if(switch_type.eq.2 .or. switch_type.eq.3) then

          if(mode .eq. 1) then
            Call c2_1d(hist_unlike_1d,href_unlike_1d,c2fit_unlike_1d,
     1           c2err_unlike_1d,max_h_1d,max_c2_1d,n_1d_total,
     2           num_pairs_unlike,num_pairs_unlike_ref)
          else if (mode .eq. 2) then
            Call c2_1d(htmp_unlike_1d,href_unlike_1d,c2fit_unlike_1d,
     1           c2err_unlike_1d,max_h_1d,max_c2_1d,n_1d_total,
     2           num_pairs_unlike,num_pairs_unlike_ref)
          else if (mode .eq. 3) then
            Call c2_1d(hinc_unlike_1d,href_unlike_1d,c2fit_unlike_1d,
     1           c2err_unlike_1d,max_h_1d,max_c2_1d,n_1d_total,
     2           num_pairs_unlike_inc,num_pairs_unlike_ref)
          end if
        end if
      end if        ! End 1D correlations

CCC   Compute 3D Correlation Functions and Errors:

      if(switch_3d .gt. 0) then
        if(switch_type.eq.1 .or. switch_type.eq.3) then

          if(mode .eq. 1) then
            Call c2_3d(hist_like_3d_fine,href_like_3d_fine,
     1           c2fit_like_3d_fine,c2err_like_3d_fine,
     2           max_h_3d,max_c2_3d,n_3d_fine,
     3           num_pairs_like,num_pairs_like_ref)
            Call c2_3d(hist_like_3d_coarse,href_like_3d_coarse,
     1           c2fit_like_3d_coarse,c2err_like_3d_coarse,
     2           max_h_3d,max_c2_3d,n_3d_coarse,
     3           num_pairs_like,num_pairs_like_ref)
          else if(mode .eq. 2) then
            Call c2_3d(htmp_like_3d_fine,href_like_3d_fine,
     1           c2fit_like_3d_fine,c2err_like_3d_fine,
     2           max_h_3d,max_c2_3d,n_3d_fine,
     3           num_pairs_like,num_pairs_like_ref)
            Call c2_3d(htmp_like_3d_coarse,href_like_3d_coarse,
     1           c2fit_like_3d_coarse,c2err_like_3d_coarse,
     2           max_h_3d,max_c2_3d,n_3d_coarse,
     3           num_pairs_like,num_pairs_like_ref)
          else if(mode .eq. 3) then
            Call c2_3d(hinc_like_3d_fine,href_like_3d_fine,
     1           c2fit_like_3d_fine,c2err_like_3d_fine,
     2           max_h_3d,max_c2_3d,n_3d_fine,
     3           num_pairs_like_inc,num_pairs_like_ref)
            Call c2_3d(hinc_like_3d_coarse,href_like_3d_coarse,
     1           c2fit_like_3d_coarse,c2err_like_3d_coarse,
     2           max_h_3d,max_c2_3d,n_3d_coarse,
     3           num_pairs_like_inc,num_pairs_like_ref)
          end if

        end if

        if(switch_type.eq.2 .or. switch_type.eq.3) then

          if(mode .eq. 1) then
            Call c2_3d(hist_unlike_3d_fine,href_unlike_3d_fine,
     1           c2fit_unlike_3d_fine,c2err_unlike_3d_fine,
     2           max_h_3d,max_c2_3d,n_3d_fine,
     3           num_pairs_unlike,num_pairs_unlike_ref)
            Call c2_3d(hist_unlike_3d_coarse,href_unlike_3d_coarse,
     1           c2fit_unlike_3d_coarse,c2err_unlike_3d_coarse,
     2           max_h_3d,max_c2_3d,n_3d_coarse,
     3           num_pairs_unlike,num_pairs_unlike_ref)
          else if(mode .eq. 2) then
            Call c2_3d(htmp_unlike_3d_fine,href_unlike_3d_fine,
     1           c2fit_unlike_3d_fine,c2err_unlike_3d_fine,
     2           max_h_3d,max_c2_3d,n_3d_fine,
     3           num_pairs_unlike,num_pairs_unlike_ref)
            Call c2_3d(htmp_unlike_3d_coarse,href_unlike_3d_coarse,
     1           c2fit_unlike_3d_coarse,c2err_unlike_3d_coarse,
     2           max_h_3d,max_c2_3d,n_3d_coarse,
     3           num_pairs_unlike,num_pairs_unlike_ref)
          else if(mode .eq. 3) then
            Call c2_3d(hinc_unlike_3d_fine,href_unlike_3d_fine,
     1           c2fit_unlike_3d_fine,c2err_unlike_3d_fine,
     2           max_h_3d,max_c2_3d,n_3d_fine,
     3           num_pairs_unlike_inc,num_pairs_unlike_ref)
            Call c2_3d(hinc_unlike_3d_coarse,href_unlike_3d_coarse,
     1           c2fit_unlike_3d_coarse,c2err_unlike_3d_coarse,
     2           max_h_3d,max_c2_3d,n_3d_coarse,
     3           num_pairs_unlike_inc,num_pairs_unlike_ref)
          end if
        end if
      end if     ! End 3D correlations

      Return
      END


C-----------------------------------------------------------------------


      subroutine c2_1d(h,href,c2,c2err,maxh,maxc2,n,num_pairs_sig,
     1                 num_pairs_bkg)
      implicit none

CCC   Computes the two-body correlation function for 1D distributions.
CCC   Errors are also computed.
CCC
CCC   Description of Input Variables in Argument List:
C
C     h(maxh)       = signal histogram (numerator)
C     href(maxh)    = background histogram (denominator)
C     c2(maxc2)     = correlation function = a/b
C     c2err(maxc2)  = correlation function error
C     maxh          = dimension of histogram arrays
C     maxc2         = dimension of correlation function array.
C     n             = # bins to use
C     num_pairs_sig = # pairs used in signal histogram
C     num_pairs_bkg = # pairs used in background histogram
C

CCC   Local Variable Type Declarations:

      integer*4 maxh,maxc2,n,num_pairs_sig,num_pairs_bkg
      integer*4 h(maxh), href(maxh)
      integer*4 k

      real*4    c2(maxc2), c2err(maxc2)
      real*4    a,a_error,b,b_error

      do k = 1,n
        if(href(k).le.0 .or. h(k).le.0) then
          c2(k) = 0.0
          c2err(k) = 1.0
        else
          a = float(h(k))/float(num_pairs_sig)
          a_error = sqrt(float(h(k)))/float(num_pairs_sig)
          b = float(href(k))/float(num_pairs_bkg)
          b_error = sqrt(float(href(k)))/float(num_pairs_bkg)
          c2(k) = a/b
          c2err(k) = c2(k)*sqrt((a_error/a)**2 + (b_error/b)**2)
        end if
      end do

      Return
      END

C-----------------------------------------------------------------------


      subroutine c2_3d(h,href,c2,c2err,maxh,maxc2,n,num_pairs_sig,
     1                 num_pairs_bkg)
      implicit none

CCC   Computes the two-body correlation function for 3D distributions.
CCC   Errors are also computed.
CCC
CCC   Description of Input Variables in Argument List:
C
C     h(maxh,maxh,maxh)         = 3D signal histogram (numerator)
C     href(maxh,maxh,maxh))     = 3D background histogram (denominator)
C     c2(maxc2,maxc2,maxc2)     = 3D correlation function = a/b
C     c2err(maxc2,maxc2,maxc2)  = 3D correlation function error
C     maxh          = dimension of 3D histogram arrays
C     maxc2         = dimension of 3D correlation function array.
C     n             = # bins to use
C     num_pairs_sig = # pairs used in signal histogram
C     num_pairs_bkg = # pairs used in background histogram
C

CCC   Local Variable Type Declarations:

      integer*4 maxh,maxc2,n,num_pairs_sig,num_pairs_bkg
      integer*4 h(maxh,maxh,maxh), href(maxh,maxh,maxh)
      integer*4 i,j,k

      real*4    c2(maxc2,maxc2,maxc2), c2err(maxc2,maxc2,maxc2)
      real*4    a,a_error,b,b_error

      do i = 1,n
      do j = 1,n
      do k = 1,n
        if(href(i,j,k).le.0 .or. h(i,j,k).le.0) then
          c2(i,j,k) = 0.0
          c2err(i,j,k) = 1.0
        else
          a = float(h(i,j,k))/float(num_pairs_sig)
          a_error = sqrt(float(h(i,j,k)))/float(num_pairs_sig)
          b = float(href(i,j,k))/float(num_pairs_bkg)
          b_error = sqrt(float(href(i,j,k)))/float(num_pairs_bkg)
          c2(i,j,k) = a/b
          c2err(i,j,k) = c2(i,j,k)*sqrt((a_error/a)**2 + (b_error/b)**2)
        end if
      end do
      end do
      end do

      Return
      END


C-------------------------------------------------------------------------


      subroutine chisquare(mode,chisq_like_1d,chisq_unlike_1d,
     1           chisq_like_3d_fine,chisq_unlike_3d_fine,
     2           chisq_like_3d_coarse,chisq_unlike_3d_coarse,
     3           chisq_hist1_1,chisq_hist1_2)
      implicit none

CCC   This subroutine calculates the chi-squares for the following:
C       o Like pair   1D 2-body correlations
C       o Unlike pair 1D 2-body correlations
C       o Like pair   3D, Fine Mesh   2-body correlations
C       o Unlike pair 3D, Fine Mesh   2-body correlations
C       o Like pair   3D, Coarse Mesh 2-body correlations
C       o Unlike pair 3D, Coarse Mesh 2-body correlations
C       o One-body 1D {pt,phi,eta} (summed) distributions for PID#1
C       o One-body 1D {pt,phi,eta} (summed) distributions for PID#2
C        
C         (where the separate chi-squares for the 1D pt, phi and eta
C          one-body distributions are added and only the sum is returned.)
C
C     'Mode' determines which one-body histogram is compared to the
C     reference histogram, where:
C
C     If mode = 1, then hist1* are used
C     If mode = 2, then htmp1* are used
C
C     The one-body reference histograms used in the chi-square calculation
C     are in href1*

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'

CCC   Local Variable Type Declarations:

      integer*4 i,j,k,mode

      real*4    chisq_like_1d, chisq_unlike_1d
      real*4    chisq_like_3d_fine,chisq_unlike_3d_fine
      real*4    chisq_like_3d_coarse,chisq_unlike_3d_coarse
      real*4    chisq_hist1_1,chisq_hist1_2

      real*4    n1fac,n2fac   ! # part 1(2) used/# part 1(2) used in Ref.
      real*4    avgerrsq_pt_1, avgerrsq_phi_1, avgerrsq_eta_1
      real*4    avgerrsq_pt_2, avgerrsq_phi_2, avgerrsq_eta_2
      real*4    avgerrsq_pt_ref_1,avgerrsq_phi_ref_1,avgerrsq_eta_ref_1
      real*4    avgerrsq_pt_ref_2,avgerrsq_phi_ref_2,avgerrsq_eta_ref_2
      real*4    chisq1

CCC   Initialize all chi-square values to zero:

      chisq_like_1d          = 0.0
      chisq_unlike_1d        = 0.0
      chisq_like_3d_fine     = 0.0
      chisq_unlike_3d_fine   = 0.0
      chisq_like_3d_coarse   = 0.0
      chisq_unlike_3d_coarse = 0.0
      chisq_hist1_1          = 0.0
      chisq_hist1_2          = 0.0

      if(switch_1d .gt. 0) then
        if(switch_type.eq.1 .or. switch_type.eq.3) then
          do i = 1,n_1d_total
            if(c2fit_like_1d(i) .ne. 0.0) then
              chisq_like_1d = chisq_like_1d + ((c2fit_like_1d(i)
     1          - c2mod_like_1d(i))/c2err_like_1d(i))**2
            end if
          end do
        end if
        if(switch_type.eq.2 .or. switch_type.eq.3) then
          do i = 1,n_1d_total
            if(c2fit_unlike_1d(i) .ne. 0.0) then
              chisq_unlike_1d = chisq_unlike_1d + ((c2fit_unlike_1d(i)
     1          - c2mod_unlike_1d(i))/c2err_unlike_1d(i))**2
            end if
          end do
        end if
      end if        ! End 1D correlation function, chi-square option

      if(switch_3d .gt. 0) then
        if(switch_type.eq.1 .or. switch_type.eq.3) then

          if(n_3d_fine .gt. 0) then
            do i = 1,n_3d_fine
            do j = 1,n_3d_fine
            do k = 1,n_3d_fine
              if(c2fit_like_3d_fine(i,j,k).ne.0.0) then
              chisq_like_3d_fine = chisq_like_3d_fine 
     1          + ((c2fit_like_3d_fine(i,j,k)
     2          -   c2mod_like_3d_fine(i,j,k))
     3             /c2err_like_3d_fine(i,j,k))**2
              end if
            end do
            end do
            end do
          end if

          if(n_3d_coarse .gt. 0) then
            do i = 1,n_3d_coarse
            do j = 1,n_3d_coarse
            do k = 1,n_3d_coarse
              if((i+j+k).gt.3) then
              if(c2fit_like_3d_coarse(i,j,k).ne.0.0) then
              chisq_like_3d_coarse = chisq_like_3d_coarse 
     1          +((c2fit_like_3d_coarse(i,j,k)
     2          -  c2mod_like_3d_coarse(i,j,k))
     3            /c2err_like_3d_coarse(i,j,k))**2
              end if
              end if
            end do
            end do
            end do
          end if

        end if

        if(switch_type.eq.2 .or. switch_type.eq.3) then

          if(n_3d_fine .gt. 0) then
            do i = 1,n_3d_fine
            do j = 1,n_3d_fine
            do k = 1,n_3d_fine
              if(c2fit_unlike_3d_fine(i,j,k).ne.0.0) then
              chisq_unlike_3d_fine = chisq_unlike_3d_fine 
     1          + ((c2fit_unlike_3d_fine(i,j,k)
     2          -   c2mod_unlike_3d_fine(i,j,k))
     3             /c2err_unlike_3d_fine(i,j,k))**2
              end if
            end do
            end do
            end do
          end if

          if(n_3d_coarse .gt. 0) then
            do i = 1,n_3d_coarse
            do j = 1,n_3d_coarse
            do k = 1,n_3d_coarse
              if((i+j+k).gt.3) then
              if(c2fit_unlike_3d_coarse(i,j,k).ne.0.0) then
              chisq_unlike_3d_coarse = chisq_unlike_3d_coarse 
     1          +((c2fit_unlike_3d_coarse(i,j,k)
     2          -  c2mod_unlike_3d_coarse(i,j,k))
     3            /c2err_unlike_3d_coarse(i,j,k))**2
              end if
              end if
            end do
            end do
            end do
          end if

        end if
      end if         !  End of 3D Correlation Function, Chi-Square Option

CCC   Obtain chi-squares for one-body distributions

      if(pid(1) .gt. 0) then
        n1fac = float(n_part_used_1_trk)/float(n_part_used_1_ref)
        avgerrsq_pt_1      = float(n_part_used_1_trk)/float(n_pt_bins)
        avgerrsq_phi_1     = float(n_part_used_1_trk)/float(n_phi_bins)
        avgerrsq_eta_1     = float(n_part_used_1_trk)/float(n_eta_bins)
        avgerrsq_pt_ref_1  = float(n_part_used_1_ref)/float(n_pt_bins)
        avgerrsq_phi_ref_1 = float(n_part_used_1_ref)/float(n_phi_bins)
        avgerrsq_eta_ref_1 = float(n_part_used_1_ref)/float(n_eta_bins)
      end if

      if(pid(2) .gt. 0) then
        n2fac = float(n_part_used_2_trk)/float(n_part_used_2_ref)
        avgerrsq_pt_2      = float(n_part_used_2_trk)/float(n_pt_bins)
        avgerrsq_phi_2     = float(n_part_used_2_trk)/float(n_phi_bins)
        avgerrsq_eta_2     = float(n_part_used_2_trk)/float(n_eta_bins)
        avgerrsq_pt_ref_2  = float(n_part_used_2_ref)/float(n_pt_bins)
        avgerrsq_phi_ref_2 = float(n_part_used_2_ref)/float(n_phi_bins)
        avgerrsq_eta_ref_2 = float(n_part_used_2_ref)/float(n_eta_bins)
      end if

      if(pid(1) .gt. 0) then
        if(mode .eq. 1) then

          chisq_hist1_1 = 
     1      chisq1(hist1_pt_1,href1_pt_1,max_h_1d,avgerrsq_pt_1,
     2             avgerrsq_pt_ref_1,n1fac,n_pt_bins)
     3     +chisq1(hist1_phi_1,href1_phi_1,max_h_1d,avgerrsq_phi_1,
     4             avgerrsq_phi_ref_1,n1fac,n_phi_bins)
     5     +chisq1(hist1_eta_1,href1_eta_1,max_h_1d,avgerrsq_eta_1,
     6             avgerrsq_eta_ref_1,n1fac,n_eta_bins)

        else if(mode .eq. 2) then

          chisq_hist1_1 =
     1      chisq1(htmp1_pt_1,href1_pt_1,max_h_1d,avgerrsq_pt_1,
     2             avgerrsq_pt_ref_1,n1fac,n_pt_bins)
     3     +chisq1(htmp1_phi_1,href1_phi_1,max_h_1d,avgerrsq_phi_1,
     4             avgerrsq_phi_ref_1,n1fac,n_phi_bins)
     5     +chisq1(htmp1_eta_1,href1_eta_1,max_h_1d,avgerrsq_eta_1,
     6             avgerrsq_eta_ref_1,n1fac,n_eta_bins)

        end if
      end if     ! End pid(1) one-body histogram chi-square calculation

      if(pid(2) .gt. 0) then
        if(mode .eq. 1) then

          chisq_hist1_2 = 
     1      chisq1(hist1_pt_2,href1_pt_2,max_h_1d,avgerrsq_pt_2,
     2             avgerrsq_pt_ref_2,n2fac,n_pt_bins)
     3     +chisq1(hist1_phi_2,href1_phi_2,max_h_1d,avgerrsq_phi_2,
     4             avgerrsq_phi_ref_2,n2fac,n_phi_bins)
     5     +chisq1(hist1_eta_2,href1_eta_2,max_h_1d,avgerrsq_eta_2,
     6             avgerrsq_eta_ref_2,n2fac,n_eta_bins)

        else if(mode .eq. 2) then

          chisq_hist1_2 =
     1      chisq1(htmp1_pt_2,href1_pt_2,max_h_1d,avgerrsq_pt_2,
     2             avgerrsq_pt_ref_2,n2fac,n_pt_bins)
     3     +chisq1(htmp1_phi_2,href1_phi_2,max_h_1d,avgerrsq_phi_2,
     4             avgerrsq_phi_ref_2,n2fac,n_phi_bins)
     5     +chisq1(htmp1_eta_2,href1_eta_2,max_h_1d,avgerrsq_eta_2,
     6             avgerrsq_eta_ref_2,n2fac,n_eta_bins)

        end if
      end if     ! End pid(2) one-body histogram chi-square calculation

      Return
      END

C----------------------------------------------------------------------


      real*4 function chisq1(h,href,maxh,herravgsq,hreferravgsq,
     1                       numfac,nbins)
      implicit none

CCC   Compute chi-square for 1D histogram h(), with respect to the
CCC   reference histogram, href().
C
C     h(maxh)       = 1D histogram array
C     href(maxh)    = 1D reference histogram array
C     maxh          = dimension of histogram arrays
C     herravgsq     = average error squared in histogram h's bins
C     hreferravgsq  = average error squared in ref. hist. href's bins
C     numfac        = ratio of total number of entries in h to that
C                     in href
C     nbins         = # bins to use in chi-square sum, starting at array
C                     element 1,2,... nbins (where nbins .le. maxh)
C
C     The chi-square value is returned in chisq1

CCC   Local Variable Type Declarations:

      integer*4 maxh, nbins, i
      integer*4 h(maxh),href(maxh)

      real*4 herravgsq,hreferravgsq,numfac,numfacsq
      real*4 herrsq,hreferrsq

      chisq1 = 0.0
      numfacsq = numfac*numfac

      do i = 1,nbins
        if(h(i) .gt. 0) then
          herrsq = float(h(i))
        else
          herrsq = herravgsq
        end if

        if(href(i) .gt. 0) then
          hreferrsq = float(href(i))
        else
          hreferrsq = hreferravgsq
        end if

        chisq1 = chisq1 + ((float(h(i)) - numfac*float(href(i)))**2)
     1           /(herrsq + numfacsq*hreferrsq)
      end do

      Return
      END

C-----------------------------------------------------------------------


      Subroutine write_data(mode,ievent)
      implicit none

CCC   This subroutine writes the main output file, 'hbt_simulation.out'
C     on File Unit 8.  File Unit 8 is opened and closed by the main
C     program.
C
C     Also, the computed 1- and 2-body reference histograms are printed
C     out from this subroutine on File Units 11 and 9, respectively.  These
C     files are opened/closed here.
C
C     Output content determined by input parameter 'mode', where:
C
C     Mode           Description of Output
C    -----   -----------------------------------------------------------
C      1     basic output file header
C            input and derived quantities
C
C      2     reference histograms (1 and 2-body)
C            saved to separate I/O File Unit=11,9 respectively
C
C      3     reference histogram output
C 
C      4     correlation model
C
C      5     correlation fit and one-body distributions
C            for each event, optional output 
C
C      6     inclusive one-body distributions and inclusive
C            correlation fit; projection onto 1D axes.
C 

      Include 'common_parameters.inc'
      Include 'common_mesh.inc'
      Include 'common_histograms.inc'
      Include 'common_correlations.inc'
      Include 'common_coulomb.inc'
      Include 'common_event_summary.inc'

      Include 'common_track.inc'
      Include 'common_sec_track.inc'
      Include 'common_sec_track2.inc'

CCC   Local Variable Type Declarations:

      integer*4 mode,i,j,k,ievent,ref_print,nev
     
      real*4 nfac1,nfac2,ref_error
      real*4 c2mod_proj1(max_c2_3d)
      real*4 c2mod_proj2(max_c2_3d)
      real*4 c2mod_proj3(max_c2_3d)
      real*4 c2fit_proj1(max_c2_3d)
      real*4 c2fit_proj2(max_c2_3d)
      real*4 c2fit_proj3(max_c2_3d)
      real*4 c2err_proj1(max_c2_3d)
      real*4 c2err_proj2(max_c2_3d)
      real*4 c2err_proj3(max_c2_3d)

C-------------------------------------------
      If(mode.eq.1) Then   !Basic Output Header
C-------------------------------------------

      write(8,100)          
      write(8,101)    
      write(8,100)          
C     write(8,102) n_events          
      write(8,103) n_pid_types,pid(1),mass1,pid(2),mass2         
      write(8,104) ref_control
      write(8,105) switch_1d         
      write(8,106) switch_3d         
      write(8,107) switch_type         
      write(8,108) switch_coherence         
      write(8,109) switch_coulomb         
      write(8,110) switch_fermi_bose         
      write(8,1101) trk_accep
      write(8,111) print_full,print_sector_data         
C     write(8,112) n_part_used_1_ref,n_part_used_2_ref          
C     write(8,113) n_part_used_1_inc,n_part_used_2_inc         
C     write(8,114) num_pairs_like_ref,num_pairs_unlike_ref         
C     write(8,115) num_pairs_like_inc,num_pairs_unlike_inc         
      write(8,116) lambda         
      write(8,117) R_1d          
      write(8,118) Rside,Rout,Rlong         
      write(8,119) Rperp,Rparallel,R0        
      write(8,120) Q0         
      write(8,121) irand          
      write(8,122) maxit          
      write(8,123) deltap         
      write(8,124) delchi         
      write(8,125) chisq_wt_like_1d         
      write(8,126) chisq_wt_unlike_1d         
      write(8,127) chisq_wt_like_3d_fine         
      write(8,128) chisq_wt_unlike_3d_fine         
      write(8,129) chisq_wt_like_3d_coarse         
      write(8,130) chisq_wt_unlike_3d_coarse         
      write(8,131) chisq_wt_hist1_1         
      write(8,132) chisq_wt_hist1_2         
      write(8,133)          
      write(8,134) n_pt_bins,pt_bin_size,pt_min,pt_max         
      write(8,135) n_phi_bins,phi_bin_size,phi_min,phi_max          
      write(8,136) n_eta_bins,eta_bin_size,eta_min,eta_max          
      write(8,137)          
      write(8,138) n_px_bins,delpx,px_min,px_max         
      write(8,139) n_py_bins,delpy,py_min,py_max          
      write(8,140) n_pz_bins,delpz,pz_min,pz_max          
      write(8,141) n_sectors         
      write(8,142)          
      write(8,143) n_1d_fine,n_1d_coarse,n_1d_total          
      write(8,144) binsize_1d_fine,binsize_1d_coarse         
      write(8,145) qmid_1d,qmax_1d         
      write(8,146)          
      write(8,147) n_3d_fine,n_3d_coarse,n_3d_total         
      write(8,148) binsize_3d_fine,binsize_3d_coarse         
      write(8,149) qmid_3d,qmax_3d
      write(8,150) n_3d_fine_project

CCC   Formats for Mode=1 Output         

 100  format(   15x,50('*'))      
 101  format(   15x,'*****',7x,'HBT CORRELATION SIMULATION',7x,'*****')
 102  format(///15x,'Number of Events in Event Text Input File=',I5)
 103  format(  /15x,'#PID types=',I2,'  PID#,mass=',I2,F8.5,
     1'  PID#,mass=',I2,F8.5)
 104  format(/  15x,'Reference Spectra Selection Option=',I2)
 105  format(// 15x,'Control Switches: Switch_1d         =',I2)
 106  format(   15x,'                  Switch_3d         =',I2)
 107  format(   15x,'                  Switch_type       =',I2)
 108  format(   15x,'                  Switch_coherence  =',I2)
 109  format(   15x,'                  Switch_coulomb    =',I2)
 110  format(   15x,'                  Switch_fermi_bose =',I2)
1101  format(   15x,'                  trk_accep         =',F10.7)
 111  format(/  15x,'Print Options: Full=',I2,'  Sectors=',I2)
 112  format(// 15x,
     1'Number particles used in Reference, for PID types=',2I5)
 113  format( 15x,
     1'Number particles used in Inclusive, for PID types=',2I5)
 114  format(/  15x,
     1'Number pairs used in Reference, like and unlike=',2I5)
 115  format(  15x,
     1'Number pairs used in Inclusive, like and unlike=',2I5)
 116  format(//  15x,'Correlation Model Parameters: Chaoticity =',F8.5)
 117  format(    15x,'1D Spherical Source Radius=',F8.4)
 118  format(    15x,'Bertsch-Pratt R-side,out,long=',3F8.4)
 119  format(    15x,'YKP R-perp,parallel,time=',3F8.4)
 120  format(    15x,'Coulomb parameter=',F8.5)
 121  format(//  15x,'Iteration Controls: Random # seed        =',I10)
 122  format(    15x,'                    Max # iterations     =',I5)
 123  format(    15x,'                    Momentum Shift Range =',F8.5)
 124  format(    15x,'                    Min % Chi-Sq limit   =',F8.5)
 125  format(//  15x,'CHI-Sq Weights: correl,like,  1d       =',F8.5)
 126  format(    15x,'                correl,unlike,1d       =',F8.5)
 127  format(    15x,'                correl,like,  3d_fine  =',F8.5)
 128  format(    15x,'                correl,unlike,3d_fine  =',F8.5)
 129  format(    15x,'                correl,like,  3d_coarse=',F8.5)
 130  format(    15x,'                correl,unlike,3d_coarse=',F8.5)
 131  format(    15x,'                1-body, PID#1          =',F8.5)
 132  format(    15x,'                1-body, PID#2          =',F8.5)
 133  format(//  15x,'Momentum Space Acceptance Range and 1D Bins:')
 134  format(    15x,'#bins,bin size,min,max for pt =',I5,3F8.5)
 135  format(    15x,'#bins,bin size,min,max for phi=',I5,3F8.4)
 136  format(    15x,'#bins,bin size,min,max for eta=',I5,3F8.4)
 137  format(//  15x,'Momentum Space Sectors:')
 138  format(    15x,'#sectors,sectorsize,min,max for px=',I5,3F8.4)
 139  format(    15x,'#sectors,sectorsize,min,max for py=',I5,3F8.4)
 140  format(    15x,'#sectors,sectorsize,min,max for pz=',I5,3F8.4)
 141  format(    15x,'Total Number of Sectors           =',I5)
 142  format(//  15x,'2-Body Correlations, 1-D Grid:')
 143  format(    15x,'#bins_fine,coarse,total     =',3I5)
 144  format(    15x,'bin size - fine, coarse     =',2F8.5)
 145  format(    15x,'Q mid point, Q maximum      =',2F8.5)
 146  format(//  15x,'2-Body Correlations, 3-D Grid:')
 147  format(    15x,'#bins - fine, coarse, total =',3I5)
 148  format(    15x,'bin size - fine, coarse     =',2F8.5)
 149  format(    15x,'Q mid point, Q maximum      =',2F8.5)
 150  format(    15x,'# 3D fine bin projected     =',I5)

CCC   END mode=1 Output and Formats

C-----------------------------
      Else If(mode.eq.2) Then	!Store 2- and 1-body Ref. Histograms
C-----------------------------


      open(unit=9,status='unknown',access='sequential',
     1     file='hbt_pair_reference.hist')
      open(unit=11,status='unknown',access='sequential',
     1     file='hbt_singles_reference.hist')

C     Write Pair Reference Hist:

      write(9,201) n_pid_types,pid(1),pid(2)
      write(9,202) n_pt_bins,pt_min,pt_max
      write(9,202) n_phi_bins,phi_min,phi_max
      write(9,202) n_eta_bins,eta_min,eta_max
      write(9,201) switch_1d,switch_3d,switch_type
      write(9,203) n_1d_fine,n_1d_coarse,n_3d_fine,n_3d_coarse
      write(9,204) binsize_1d_fine,binsize_1d_coarse,
     1             binsize_3d_fine,binsize_3d_coarse
      write(9,201) num_pairs_like_ref,num_pairs_unlike_ref
 201  format(2x,3I10)
 202  format(2x,I10,2E15.6)
 203  format(2x,4I10)
 204  format(2x,4E15.6)
 205  format(2x,I20)

      if(switch_1d.gt.0.and.n_1d_total.gt.0) then
     	if(switch_type.eq.1.or.switch_type.eq.3) then
          write(9,205) (href_like_1d(i),i=1,n_1d_total)
        end if
        if(switch_type.eq.2.or.switch_type.eq.3) then
          write(9,205) (href_unlike_1d(i),i=1,n_1d_total)
        endif
      endif	!End 1D Ref. Hist. Output

      if(switch_3d.gt.0) then
	if(switch_type.eq.1.or.switch_type.eq.3) then

	  if(n_3d_fine.gt.0) then
            do i=1,n_3d_fine
            do j=1,n_3d_fine
            do k=1,n_3d_fine
               write(9,205) href_like_3d_fine(i,j,k)
            enddo
            enddo
            enddo
          endif

          if(n_3d_coarse.gt.0) then
            do i=1,n_3d_coarse
            do j=1,n_3d_coarse
            do k=1,n_3d_coarse
               write(9,205) href_like_3d_coarse(i,j,k)
            enddo
            enddo
            enddo
          endif

        end if

        if(switch_type.eq.2.or.switch_type.eq.3) then

          if(n_3d_fine.gt.0) then     	
            do i=1,n_3d_fine
            do j=1,n_3d_fine
            do k=1,n_3d_fine
               write(9,205) href_unlike_3d_fine(i,j,k)
            enddo
            enddo
            enddo
          endif

          if(n_3d_coarse.gt.0) then     	
            do i=1,n_3d_coarse
            do j=1,n_3d_coarse
            do k=1,n_3d_coarse
               write(9,205) href_unlike_3d_coarse(i,j,k)
            enddo
            enddo
            enddo
          endif

        endif
      endif	!End 3D Reference Histograms Output

CC    Write One-Body  - singles histograms:
      
      write(11,201) n_pid_types,pid(1),pid(2) 
      write(11,202) n_pt_bins,pt_min,pt_max
      write(11,202) n_phi_bins,phi_min,phi_max
      write(11,202) n_eta_bins,eta_min,eta_max
      write(11,201) n_part_used_1_ref,n_part_used_2_ref

      if(pid(1).gt.0) then
        write(11,205)(href1_pt_1(i),i=1,n_pt_bins)
        write(11,205)(href1_phi_1(i),i=1,n_phi_bins)
        write(11,205)(href1_eta_1(i),i=1,n_eta_bins)
      endif
  

      if(pid(2).gt.0) then
        write(11,205)(href1_pt_2(i),i=1,n_pt_bins)
        write(11,205)(href1_phi_2(i),i=1,n_phi_bins)
        write(11,205)(href1_eta_2(i),i=1,n_eta_bins)
      endif
      
      close(unit=9)
      close(unit=11)

CCC   END mode=2 Reference Histogram Output

C-----------------------------
      Else If(mode.eq.3) Then	!Print out the Reference Histograms
C-----------------------------

      write(8,300)
      write(8,301)
      write(8,302) n_pt_bins,pt_min,pt_max
      write(8,303) n_phi_bins,phi_min,phi_max
      write(8,304) n_eta_bins,eta_min,eta_max
      write(8,305) n_part_used_1_ref,n_part_used_2_ref

      write(8,306)
      do i=1,n_pt_bins
         write(8,307) i,href1_pt_1(i),href1_pt_2(i)
      enddo

      write(8,308)
      do i=1,n_phi_bins
         write(8,307) i,href1_phi_1(i),href1_phi_2(i)
      enddo

      write(8,309)
      do i=1,n_eta_bins
         write(8,307) i,href1_eta_1(i),href1_eta_2(i)
      enddo

      write(8,310)
      write(8,311) n_1d_fine,n_1d_coarse
      write(8,312) binsize_1d_fine,binsize_1d_coarse
      write(8,313) n_3d_fine,n_3d_coarse
      write(8,314) binsize_3d_fine,binsize_3d_coarse
      write(8,315) num_pairs_like_ref,num_pairs_unlike_ref

      if(switch_1d.gt.0.and.n_1d_total.gt.0) then
        write(8,316)
            do i=1,n_1d_total
               write(8,307) i,href_like_1d(i),href_unlike_1d(i)
            enddo
      endif	!End Print Out of 2-body, 1D reference histogram

      if(switch_3d.gt.0.and.n_3d_fine.gt.0) then
        write(8,317)
        do i=1,n_3d_fine
        do j=1,n_3d_fine
        do k=1,n_3d_fine
           write(8,318) i,j,k,href_like_3d_fine(i,j,k),
     1                  href_unlike_3d_fine(i,j,k)
        enddo
        enddo
        enddo
      endif	     !End Print Out of 2-Body, 3D-Fine Mesh Ref. Hist.


      if(switch_3d.gt.0.and.n_3d_coarse.gt.0) then
        write(8,319)
        do i=1,n_3d_coarse
        do j=1,n_3d_coarse
        do k=1,n_3d_coarse
           write(8,318) i,j,k,href_like_3d_coarse(i,j,k),
     1                  href_unlike_3d_coarse(i,j,k)
        enddo
        enddo
        enddo
      endif	    !End Print Out of 2-Body, 3D-Coarse Mesh Ref. Hist.

CCC   Formats for mode=3 Output:

 300  format(///5x,15('*'),'REFERENCE HISTOGRAMS',15('*'))
 301  format(//15x,'ONE-BODY REFERENCE DISTRIBUTIONS:')
 302  format(/ 15x,'PT BINS: (#,min,max)=',I5,2F8.4)
 303  format(  15x,'PHI BINS:(#,min,max)=',I5,2F8.4)
 304  format(  15x,'ETA BINS:(#,min,max)=',I5,2F8.4)
 305  format(  15x,'Number particles used in Ref, PID type1,2=',2I8)
 306  format(/  9x,'PT',10x,'BIN#',5x,'PID-1',5x,'PID-2')
 308  format(/  9x,'PHI',9x,'BIN#',5x,'PID-1',5x,'PID-2')
 309  format(/  9x,'ETA',9x,'BIN#',5x,'PID-1',5x,'PID-2')
 307  format(  20x,I5,2I10)
 310  format(///15x,'TWO-BODY REFERENCE DISTRIBUTIONS:')
 311  format(/  15x,'#BINS FOR 1D-Fine and Coarse Grid =',2I5)
 312  format(   15x,'BIN SIZES FOR 1D-Fine and Coarse  =',2F8.5)
 313  format(   15x,'#BINS FOR 3D-Fine and Coarse Grid =',2I5)
 314  format(   15x,'BIN SIZES FOR 3D-Fine and Coarse  =',2F8.5)
 315  format(   15x,'Number of Like and Unlike Pairs For Ref. = ', 
     1             2I10)
 316  format(/5x,'2-BODY, 1D',6x,'BIN#',5x,'LIKE',5x,'UNLIKE')
 317  format(/3x,'2-BODY, 3D-FINE',2x,'BIN:i',4x,'j',4x,
     1     'k',5x,'LIKE',5x,'UNLIKE')
 318  format(   20x,3I5,2I10)
 319  format(/2x,'2-BODY, 3D-COARSE',1x,'BIN:i',4x,'j',4x,
     1     'k',5x,'LIKE',5x,'UNLIKE')

CCC   END mode=3  Output and Formats

C----------------------------
      Else If(mode.eq.4) Then	!Print Correlation Function Model
C----------------------------

      write(8,400)
      write(8,311) n_1d_fine,n_1d_coarse
      write(8,312) binsize_1d_fine,binsize_1d_coarse
      write(8,313) n_3d_fine,n_3d_coarse
      write(8,314) binsize_3d_fine,binsize_3d_coarse


      if(switch_1d.gt.0.and.n_1d_total.gt.0) then
        write(8,316)
            do i=1,n_1d_total
               write(8,407) i,c2mod_like_1d(i),c2mod_unlike_1d(i)
            enddo
      endif	!End Print Out of 2-body, 1D Model Correction Functions

      if(switch_3d.gt.0.and.n_3d_fine.gt.0) then
        write(8,317)
        do i=1,n_3d_fine
        do j=1,n_3d_fine
        do k=1,n_3d_fine
          write(8,418) i,j,k,c2mod_like_3d_fine(i,j,k),
     1                 c2mod_unlike_3d_fine(i,j,k)
        enddo
        enddo
        enddo
      endif  !End Print Out of 2-Body, 3D-Fine mesh Model Correl. Function 

      if(switch_3d.gt.0.and.n_3d_coarse.gt.0) then
        write(8,319)
        do i=1,n_3d_coarse
        do j=1,n_3d_coarse
        do k=1,n_3d_coarse
          write(8,418) i,j,k,c2mod_like_3d_coarse(i,j,k),
     1                 c2mod_unlike_3d_coarse(i,j,k)
        enddo
        enddo
        enddo
      endif  !End Print Out of 2-Body, 3D-Coarse Model Correlation Function

      if(switch_coulomb.eq.3) then ! Print interpolated Pratt Model
CC                                 ! Coulomb Correction for Finite
CC                                 ! Source Radius Q0
        write(8,401)Q0
        write(8,402)
        do i=1,max_c2_coul
           write(8,403) i,q_coul(i),c2_coul_like(i),
     1                  c2_coul_unlike(i)
        enddo
      endif	    

CCC   Additional Formats for Mode=4 Output:

 400  format(///5x,15('*'),'MODEL CORRELATION FUNCTIONS',15('*'))
 401  format(///15x,'COULOMB SOURCE RADIUS FOR PRATT MODEL=',F8.4)
 402  format(// 5x,'q-bin',2x,'q',4x,'C2_coul_like',2x,
     1             'C2_coul_unlike')
 403  format(   5x,I5,3E15.6)
 407  format(  20x,I5,2F10.7)
 418  format(   20x,3I5,2F10.7)
 
CCC   END MODE = 4  OUTPUT AND FORMATS

C------------------------------
      Else If(mode.eq.5) Then   ! Optional Output for 1- and 2-Body Fits
C------------------------------ ! for each event.
C      write(*,*)'Event ', ievent
C      write(*,*)'chisq_total_store(event) ',chisq_total_store(ievent)
      write(8,500) ievent
      write(8,501) n_part_1_trk,n_part_2_trk,n_part_tot_trk
      write(8,502) n_part_used_1_trk, n_part_used_2_trk
      write(8,503) num_pairs_like, num_pairs_unlike

CCC   Output one-body distributions for event:

      write(8,504)
      if(pid(1) .gt. 0) then
         nfac1 = float(n_part_used_1_trk)/float(n_part_used_1_ref)
         write(8,505) nfac1

         write(8,507)
         do i = 1,n_pt_bins
           ref_print = int(nfac1*float(href1_pt_1(i)))
           ref_error = nfac1*sqrt(float(href1_pt_1(i)))
           write(8,510) i,hist1_pt_1(i),ref_print,ref_error
         end do

         write(8,508)
         do i = 1,n_phi_bins
           ref_print = int(nfac1*float(href1_phi_1(i)))
           ref_error = nfac1*sqrt(float(href1_phi_1(i)))
           write(8,510) i,hist1_phi_1(i),ref_print,ref_error
         end do

         write(8,509)
         do i = 1,n_eta_bins
           ref_print = int(nfac1*float(href1_eta_1(i)))
           ref_error = nfac1*sqrt(float(href1_eta_1(i)))
           write(8,510) i,hist1_eta_1(i),ref_print,ref_error
         end do

      end if    !  End PID # 1, One-Body Distribution Output

      if(pid(2) .gt. 0) then
         nfac2 = float(n_part_used_2_trk)/float(n_part_used_2_ref)
         write(8,506) nfac2

         write(8,507)
         do i = 1,n_pt_bins
           ref_print = int(nfac2*float(href1_pt_2(i)))
           ref_error = nfac2*sqrt(float(href1_pt_2(i)))
           write(8,510) i,hist1_pt_2(i),ref_print,ref_error
         end do

         write(8,508)
         do i = 1,n_phi_bins
           ref_print = int(nfac2*float(href1_phi_2(i)))
           ref_error = nfac2*sqrt(float(href1_phi_2(i)))
           write(8,510) i,hist1_phi_2(i),ref_print,ref_error
         end do

         write(8,509)
         do i = 1,n_eta_bins
           ref_print = int(nfac2*float(href1_eta_2(i)))
           ref_error = nfac2*sqrt(float(href1_eta_2(i)))
           write(8,510) i,hist1_eta_2(i),ref_print,ref_error
         end do

      end if    !  End PID # 2, One-Body Distribution Output

CCC   Output Two-Body Correlation Functions for Event:

      write(8,520)
      if(switch_1d.gt.0 .and. n_1d_total.gt.0) then
         write(8,530)
         write(8,521)
         write(8,522)
         do i = 1,n_1d_total
            write(8,523) i,c2mod_like_1d(i),c2fit_like_1d(i),
     1           c2err_like_1d(i),c2mod_unlike_1d(i),
     2           c2fit_unlike_1d(i),c2err_unlike_1d(i)
         end do
      end if    !  End 1D Correlation Model and Fit Output 

      if(switch_3d.gt.0 .and. n_3d_fine.gt.0) then
         write(8,531)
         write(8,524)
         write(8,525)
         do i = 1,n_3d_fine
         do j = 1,n_3d_fine
         do k = 1,n_3d_fine
            write(8,526) i,j,k,c2mod_like_3d_fine(i,j,k),
     1      c2fit_like_3d_fine(i,j,k),c2err_like_3d_fine(i,j,k),
     2      c2mod_unlike_3d_fine(i,j,k),c2fit_unlike_3d_fine(i,j,k),
     3      c2err_unlike_3d_fine(i,j,k)
         end do
         end do
         end do
      end if    !  End 3D Fine Mesh Correlation Model and Fit Output

      if(switch_3d.gt.0 .and. n_3d_coarse.gt.0) then
         write(8,532)
         write(8,524)
         write(8,525)
         do i = 1,n_3d_coarse
         do j = 1,n_3d_coarse
         do k = 1,n_3d_coarse
            write(8,526) i,j,k,c2mod_like_3d_coarse(i,j,k),
     1      c2fit_like_3d_coarse(i,j,k),c2err_like_3d_coarse(i,j,k),
     2      c2mod_unlike_3d_coarse(i,j,k),c2fit_unlike_3d_coarse(i,j,k),
     3      c2err_unlike_3d_coarse(i,j,k)
         end do
         end do
         end do
      end if    !  End 3D Coarse Mesh Correlation Model and Fit Output

CCC   Output Event Summary and Chi-Square Information for Event:

      write(8,539) ievent
      write(8,540) num_iter(ievent)
      write(8,541) n_part_used_1_store(ievent),
     1             n_part_used_2_store(ievent)
      write(8,5411) n_part_tot_store(ievent)
      write(8,542) num_sec_flagged_store(ievent)
      write(8,543) frac_trks_out(ievent),frac_trks_flag(ievent)
      write(8,544) chisq_like_1d_store(ievent),
     1             chisq_unlike_1d_store(ievent)
      write(8,545) chisq_like_3d_fine_store(ievent),
     1             chisq_unlike_3d_fine_store(ievent)
      write(8,546) chisq_like_3d_coarse_store(ievent),
     1             chisq_unlike_3d_coarse_store(ievent)
      write(8,547) chisq_hist1_1_store(ievent),
     1             chisq_hist1_2_store(ievent)
      write(8,548) chisq_total_store(ievent)

CCC   Formats for Mode = 5 Output:

500   Format(///5x,5('*'),'Fitted 1-Body Distributions and ',
     1      'Correlations for Event #',I5,5('*'))
501   Format(//15x,'Number of Particles of PID types 1,2,total = ',
     1      3I5)
502   Format(  15x,'Number of Particles of PID types 1,2 Used  = ',
     1      2I5)
503   Format(  15x,'Number of Pairs Used - Like and Unlike = ',2I10)
504   Format(//5x,'Fitted and Normalized Reference One-Body ',
     1      'Distributions')
505   Format( /10x,'Particle Type 1: Reference Scale Factor = ',E12.5)
506   Format( /10x,'Particle Type 2: Reference Scale Factor = ',E12.5)
507   Format(/2x,' PT: BIN#',5x,'hist1',7x,'href1-scaled',2x,
     1       'ref-err-scaled')
508   Format(/2x,'PHI: BIN#',5x,'hist1',7x,'href1-scaled',2x,
     1       'ref-err-scaled')
509   Format(/2x,'ETA: BIN#',5x,'hist1',7x,'href1-scaled',2x,
     1       'ref-err-scaled')
510   Format(7x,I4,3x,I7,8x,I7,7x,F10.5)
520   Format(//5x,'Model and Fitted Correlations')
530   Format(//21x,'One-Dimensional Fit - Fine & Coarse Mesh')
531   Format(//25x,'Three-Dimensional Fit - Fine Mesh')
532   Format(//24x,'Three-Dimensional Fit - Coarse Mesh')
521   Format(/1x,'BIN',13x,'LIKE PAIRS',27x,'UNLIKE PAIRS')
522   Format(8x,'MOD',9x,'FIT',9x,'ERR',11x,'MOD',9x,'FIT',9x,
     1       'ERR',/)
523   Format(1x,I3,3E12.4,2x,3E12.4)
524   Format(/2x,'BINS',12x,'LIKE PAIRS',24x,'UNLIKE PAIRS')
525   Format(1x,' i j k',4x,'MOD',8x,'FIT',8x,'ERR',10x,'MOD',
     1      8x,'FIT',8x,'ERR',/)
526   Format(1x,3I2,3E11.4,2x,3E11.4)
539   Format(///10x,'Event and Chi-Square Summary for Event #',I5)
540   Format( //15x,'Number of Iterations             =',F10.2)
541   Format(   15x,'# Particles Used for PID Types1,2=',2F10.2)
5411  Format(   15x,'Total # Particles in track table =',F10.2)
542   Format(   15x,'# Sectors Flagged                =',F10.2)
543   Format(   15x,'Frac Trks Out of Accep., Flagged =',2E11.4)
544   Format(   15x,'Chi-Sq: 1D - Like & Unlike       =',2E11.4)
545   Format(   15x,'Chi-Sq: 3D - Fine -Like & Unlike =',2E11.4)
546   Format(   15x,'Chi-Sq: 3D - Coarse-Like &Unlike =',2E11.4)
547   Format(   15x,'Chi-Sq: One-Body Dist. PID# 1&2  =',2E11.4)
548   Format(   15x,'Chi-Sq: Total Weighted           =',E11.4)

CCC   End Mode = 5 Output and Formats
      
C------------------------------
      Else If(mode.eq.6) Then	! Inclusive 1 & 2 Body Output
C------------------------------
      write(8,600) n_events
      write(8,601) n_part_used_1_inc,n_part_used_2_inc
      write(8,602) num_pairs_like_inc,num_pairs_unlike_inc

      write(8,603) 
      if(pid(1).gt.0) then
C    Division by zero check      
        IF (n_part_used_1_ref .LE. 0) THEN
         PRINT*,'************************************'
         PRINT*,'*        HBT PROCESSOR             *'
         PRINT*,'* Number of particles selected for *'
         PRINT*,'*   processing is less or equal    *'
         PRINT*,'* !!!!!!!!  ZER0  !!!!!!!!!!       *'
         PRINT*,'* unable to proceed                *'
         PRINT*,'* EXITING FORTRAN                  *'
         PRINT*,'*                                  *'
         PRINT*,'* HINT: broad the parameter regions*'
         PRINT*,'* OR/AND number of particles OR/AND*'
         PRINT*,'* number of events                 *'
         PRINT*,'************************************'
         WRITE(7,5481)
5481     FORMAT(5x,'Number of particles selected for processing is',
     1             ' less or equal 0',
     2             ' - STOP')
         errorcode = 1
         Return   
        END IF
        nfac1=float(n_part_used_1_inc)/float(n_part_used_1_ref)
        write(8,604) nfac1
        write(8,605)
        do i = 1,n_pt_bins
          ref_print=int(nfac1*float(href1_pt_1(i)))
          ref_error=nfac1*sqrt(float(href1_pt_1(i)))
          write(8,510) i,hinc1_pt_1(i),ref_print,ref_error
        enddo

        write(8,606)
        do i = 1,n_phi_bins
          ref_print=int(nfac1*float(href1_phi_1(i)))
          ref_error=nfac1*sqrt(float(href1_phi_1(i)))
          write(8,510) i,hinc1_phi_1(i),ref_print,ref_error
        enddo

        write(8,607)
        do i = 1,n_eta_bins
          ref_print=int(nfac1*float(href1_eta_1(i)))
          ref_error=nfac1*sqrt(float(href1_eta_1(i)))
          write(8,510) i,hinc1_eta_1(i),ref_print,ref_error
        enddo

      endif	!END PID #1  One-BODY INCL. DISTRIBUTION OUTPUT

      if(pid(2).gt.0) then
        nfac2=float(n_part_used_2_inc)/float(n_part_used_2_ref)
        write(8,608) nfac2
        write(8,605)
        do i = 1,n_pt_bins
          ref_print=int(nfac2*float(href1_pt_2(i)))
          ref_error=nfac2*sqrt(float(href1_pt_2(i)))
          write(8,510) i,hinc1_pt_2(i),ref_print,ref_error
        enddo

        write(8,606)
        do i = 1,n_phi_bins
          ref_print=int(nfac2*float(href1_phi_2(i)))
          ref_error=nfac2*sqrt(float(href1_phi_2(i)))
          write(8,510) i,hinc1_phi_2(i),ref_print,ref_error
        enddo

        write(8,607)
        do i = 1,n_eta_bins
          ref_print=int(nfac2*float(href1_eta_2(i)))
          ref_error=nfac2*sqrt(float(href1_eta_2(i)))
          write(8,510) i,hinc1_eta_2(i),ref_print,ref_error
        enddo

      endif	!END PID #2  One-BODY INCl. DISTRIBUTION OUTPUT

CC    OUTPUT TWO-BODY INCLUSIVE HISTOGRAMS:

      write(8,660)
      if(switch_1d.gt.0.and.n_1d_total.gt.0) then
         write(8,316)
         do i=1,n_1d_total
            write(8,307) i,hinc_like_1d(i),hinc_unlike_1d(i)
         end do
      end if   ! End Print out of 2-Body, 1D Inclusive Histograms

      if(switch_3d.gt.0.and.n_3d_fine.gt.0) then
         write(8,317)
         do i=1,n_3d_fine
         do j=1,n_3d_fine
         do k=1,n_3d_fine
            write(8,318) i,j,k,hinc_like_3d_fine(i,j,k),
     1                         hinc_unlike_3d_fine(i,j,k)
        enddo
        enddo
        enddo
      endif   ! End Print out of 2-Body, 3D-Fine Inclusive Histograms

      if(switch_3d.gt.0.and.n_3d_coarse.gt.0) then
         write(8,319)
         do i=1,n_3d_coarse
         do j=1,n_3d_coarse
         do k=1,n_3d_coarse
            write(8,318) i,j,k,hinc_like_3d_coarse(i,j,k),
     1                         hinc_unlike_3d_coarse(i,j,k)
        enddo
        enddo
        enddo
      endif   ! End Print out of 2-Body, 3D-Coarse Inclusive Histograms

CC    OUTPUT TWO-BODY INCL.CORRELATION FUNCTIONS FOR EVENT
      
      write(8,620)
      if(switch_1d.gt.0.and.n_1d_total.gt.0) then
        write(8,530)
        write(8,521)
        write(8,522)
        do i=1,n_1d_total
           write(8,523) i,c2mod_like_1d(i),c2fit_like_1d(i),
     1                  c2err_like_1d(i),c2mod_unlike_1d(i),
     2                  c2fit_unlike_1d(i),c2err_unlike_1d(i)
        enddo
      endif

      if(switch_3d.gt.0.and.n_3d_fine.gt.0) then
        write(8,531)
        write(8,524)
        write(8,525)
        do i=1,n_3d_fine
        do j=1,n_3d_fine
        do k=1,n_3d_fine
           write(8,526) i,j,k,c2mod_like_3d_fine(i,j,k),
     1                       c2fit_like_3d_fine(i,j,k),
     2                       c2err_like_3d_fine(i,j,k),
     3                       c2mod_unlike_3d_fine(i,j,k),
     4                       c2fit_unlike_3d_fine(i,j,k),
     5                       c2err_unlike_3d_fine(i,j,k)

        enddo
        enddo
        enddo
      endif

      if(switch_3d.gt.0.and.n_3d_coarse.gt.0) then
        write(8,532)
        write(8,524)
        write(8,525)
        do i=1,n_3d_coarse
        do j=1,n_3d_coarse
        do k=1,n_3d_coarse
           write(8,526) i,j,k,c2mod_like_3d_coarse(i,j,k),
     1                       c2fit_like_3d_coarse(i,j,k),
     2                       c2err_like_3d_coarse(i,j,k),
     3                       c2mod_unlike_3d_coarse(i,j,k),
     4                       c2fit_unlike_3d_coarse(i,j,k),
     5                       c2err_unlike_3d_coarse(i,j,k)

        enddo
        enddo
        enddo
      endif

CCC   Compute and Print 1D projections of 3D fine mesh C2 model,
CCC   fit and errors for like and unlike pairs.

      if(switch_3d .gt. 0 .and. n_3d_fine .gt. 0) then
         if(switch_type .eq. 1 .or. switch_type .eq. 3) then
            Call c2_3d_projected(hinc_like_3d_fine,
     1         href_like_3d_fine,c2mod_like_3d_fine,
     2         c2mod_proj1,c2mod_proj2,c2mod_proj3,
     3         c2fit_proj1,c2fit_proj2,c2fit_proj3,
     4         c2err_proj1,c2err_proj2,c2err_proj3,
     5         max_h_3d, max_c2_3d, n_3d_fine,
     6         n_3d_fine_project,num_pairs_like_inc,
     7         num_pairs_like_ref)
            write(8,650)
            write(8,651)
            write(8,657)
            do i = 1,n_3d_fine
             write(8,658) i,c2mod_proj1(i),c2fit_proj1(i),c2err_proj1(i)
            end do
            write(8,652)
            write(8,657)
            do i = 1,n_3d_fine
             write(8,658) i,c2mod_proj2(i),c2fit_proj2(i),c2err_proj2(i)
            end do
            write(8,653)
            write(8,657)
            do i = 1,n_3d_fine
             write(8,658) i,c2mod_proj3(i),c2fit_proj3(i),c2err_proj3(i)
            end do
         end if       !  End Like pair output

         if(switch_type .eq. 2 .or. switch_type .eq. 3) then
            Call c2_3d_projected(hinc_unlike_3d_fine,
     1         href_unlike_3d_fine,c2mod_unlike_3d_fine,
     2         c2mod_proj1,c2mod_proj2,c2mod_proj3,
     3         c2fit_proj1,c2fit_proj2,c2fit_proj3,
     4         c2err_proj1,c2err_proj2,c2err_proj3,
     5         max_h_3d, max_c2_3d, n_3d_fine,
     6         n_3d_fine_project,num_pairs_unlike_inc,
     7         num_pairs_unlike_ref)
            write(8,654)
            write(8,657)
            do i = 1,n_3d_fine
             write(8,658) i,c2mod_proj1(i),c2fit_proj1(i),c2err_proj1(i)
            end do
            write(8,655)
            write(8,657)
            do i = 1,n_3d_fine
             write(8,658) i,c2mod_proj2(i),c2fit_proj2(i),c2err_proj2(i)
            end do
            write(8,656)
            write(8,657)
            do i = 1,n_3d_fine
             write(8,658) i,c2mod_proj3(i),c2fit_proj3(i),c2err_proj3(i)
            end do
         end if       !  End Unlike pair output
      end if          !  End 1D projections


CCC   EVENT AND CHISQ SUMMARY INFORMATION:

      if(n_events.le.max_events) then
        nev=n_events
      else
        nev=max_events
      endif
 
      write(8,621)
      write(8,622)

      do i=1,nev
        write(8,623) i,num_iter(i),n_part_used_1_store(i),
     1               n_part_used_2_store(i),
     2               num_sec_flagged_store(i),
     3               frac_trks_out(i),frac_trks_flag(i),
     4               chisq_total_store(i)
      enddo

      write(8,6231) trk_maxlen
      write(8,6232)
      do i=1,nev
        write(8,6233) i,n_part_tot_store(i)
      end do

      write(8,624)
      do i=1,nev
        write(8,625) i,chisq_like_1d_store(i),
     1               chisq_unlike_1d_store(i)
      enddo


      write(8,626)
      do i=1,nev
        write(8,625) i,chisq_like_3d_fine_store(i),
     1               chisq_unlike_3d_fine_store(i)
      enddo


      write(8,627)
      do i=1,nev
        write(8,625) i,chisq_like_3d_coarse_store(i),
     1               chisq_unlike_3d_coarse_store(i)
      enddo


      write(8,628)
      do i=1,nev
        write(8,625) i,chisq_hist1_1_store(i),
     1                 chisq_hist1_2_store(i)
      enddo

CCC   Output the Mean and RMS values for the Event Loop:

      write(8,629)
      write(8,630)
      write(8,631) niter_mean,niter_rms
      write(8,632) npart1_mean,npart1_rms
      write(8,633) npart2_mean,npart2_rms
      write(8,6331) npart_tot_mean, npart_tot_rms
      write(8,634) nsec_flag_mean,nsec_flag_rms
      write(8,635) frac_trks_out_mean,frac_trks_out_rms
      write(8,636) frac_trks_flag_mean,frac_trks_flag_rms
      write(8,637) chi_l1d_mean,chi_l1d_rms
      write(8,638) chi_u1d_mean,chi_u1d_rms
      write(8,639) chi_l3f_mean,chi_l3f_rms
      write(8,640) chi_u3f_mean,chi_u3f_rms
      write(8,641) chi_l3c_mean,chi_l3c_rms
      write(8,642) chi_u3c_mean,chi_u3c_rms
      write(8,643) chi_1_1_mean,chi_1_1_rms
      write(8,644) chi_1_2_mean,chi_1_2_rms
      write(8,645) chi_tot_mean,chi_tot_rms

CCC   FORMATS FOR MODE = 6 OUTPUT

 600  format(/// 2x,'FITTED 1-BODY DIST. AND CORRELATIONS ',
     1      'FOR INCLUSIVE SUM OF',I5,' EVENTS')
 601  format(// 15x,'Inclusive # Particles USED of PID ', 
     1      'types 1,2=',2I8)
 602  format(   15x,'Inclusive # of pairs used; like/unlike=',
     1       2I10)
 603  format(//  5x,'Inclusive and Normalized Reference ',
     1       'One-Body Distributions')
 604  format(/  10x,'Inclusive: Particle Type 1 - Reference ',
     1       'Scale Factor=',E12.5)
 605  format(/   2x,' PT: BIN#',5x,'hinc1',7x,'href1-scaled',2x,
     1      'ref-err-scaled')
 606  format(/   2x,'PHI: BIN#',5x,'hinc1',7x,'href1-scaled',2x,
     1      'ref-err-scaled')
 607  format(/   2x,'ETA: BIN#',5x,'hinc1',7x,'href1-scaled',2x,
     1      'ref-err-scaled')
 608  format(/  10x,'Inclusive: Particle Type 2 - ',
     1      'Reference Scale Factor=',E12.5)
 620  format(//  5x,'MODEL AND INCLUSIVE FITTED CORRELATIONS')
 621  format(// 15x,'Event and Chi-Square Summary Lists')
 622  format(/   3x,'event',2x,'#iter',3x,'#PID1',4x,'#PID2',3x,
     1      '#sec-flg',3x,'frac-out',4x,'frac-flg',3x,'CHISQ-TOT')
 623  format(3x,I5,2x,F6.0,2(1x,F8.0),1x,F9.0,
     1           2(1x,F11.8),1x,E11.4)
6231  format(/5x,'Max# tracks allowed in track table = ',I8)
6232  format(/5x,'event',4x,'Tot# trks')
6233  format(5x,I5,F12.2)
 624  format(/5x,'event',4x,'CHI_l1d',8x,'CHI_u1d')
 626  format(/5x,'event',4x,'CHI_l3f',8x,'CHI_u3f')
 627  format(/5x,'event',4x,'CHI_l3c',8x,'CHI_u3c')
 628  format(/5x,'event',4x,'CHI_1-1',8x,'CHI_1-2')
 625  format(5x,I5,2E15.6)
 629  format(// 10x,'Event and Chi-Square Summary - ',
     1     'Mean and RMS Values')
 630  format(/  14x,'Quantity',15x,'Mean',11x,'RMS')
 631  format(    5x,'Number of Iterations      ',2E15.6)
 632  format(    5x,'#PID Type 1               ',2E15.6)
 633  format(    5x,'#PID Type 2               ',2E15.6)
6331  format(    5x,'Tot # Tracks in Table     ',2E15.6)
 634  format(    5x,'#Sectors Flagged          ',2E15.6)
 635  format(    5x,'Frac. Trks Out of Accept. ',2E15.6)
 636  format(    5x,'Frac. Trks Flagged        ',2E15.6)
 637  format(    5x,'CHISQ like 1D             ',2E15.6)
 638  format(    5x,'CHISQ unlike 1D           ',2E15.6)
 639  format(    5x,'CHISQ like 3D Fine        ',2E15.6)
 640  format(    5x,'CHISQ unlike 3D Fine      ',2E15.6)
 641  format(    5x,'CHISQ like 3D Coarse      ',2E15.6)
 642  format(    5x,'CHISQ unlike 3D Coarse    ',2E15.6)
 643  format(    5x,'CHISQ 1 Body #1           ',2E15.6)
 644  format(    5x,'CHISQ 1 Body #2           ',2E15.6)
 645  format(    5x,'CHISQ Total               ',2E15.6)
 650  format(//10x ,'Inclusive Three-Dimensional Projected Fits -',
     1       ' Fine Mesh')
 651  format( /25x ,'Like Pairs - Axis #1 ')
 652  format( /25x ,'Like Pairs - Axis #2 ')
 653  format( /25x ,'Like Pairs - Axis #3 ')
 654  format( /25x ,'Unlike Pairs - Axis #1 ')
 655  format( /25x ,'Unlike Pairs - Axis #2 ')
 656  format( /25x ,'Unlike Pairs - Axis #3 ')
 657  format(    2x,'BIN#',3x,'Model',8x,'Fit',8x,'Error')
 658  format(    3x,I3,3E12.4)
 660  format(//  5x,'INCLUSIVE TWO-BODY HISTOGRAMS')

CCC   END MODE = 6 OUTPUT AND FORMATS

C----------------
      END IF
C----------------

      Return
      END 

C-----------------------------------------------------------------------


      subroutine c2_3d_projected(h,href,c2mod,
     1           c2mod_proj1,c2mod_proj2,c2mod_proj3,
     2           c2fit_proj1,c2fit_proj2,c2fit_proj3,
     3           c2err_proj1,c2err_proj2,c2err_proj3,
     4           maxh,maxc2,n,n_proj,num_pairs_sig,
     5           num_pairs_bkg)

      implicit none

CCC   This Subroutine computes the projected two-body correlation
CCC   function for 3D distributions - fine mesh only; for both the
CCC   correlation model (weighted with the reference histogram) and
CCC   the inclusive correlation fit.
CCC
CCC   Description of Input Variables in the Argument List:
CCC
CCC      h(maxh,maxh,maxh)         = 3D fine mesh inclusive signal histog.
CCC      href(maxh,maxh,maxh)      = 3D fine mesh inclusive background hist.
CCC      c2mod(maxc2,maxc2,maxc2)  = 3D fine mesh correlation model
CCC      maxh       = Dimension of 3D fine mesh histogram arrays
CCC      maxc2      = Dimension of 3D fine mesh correlation function arrays
CCC      n          = Number of bins to use
CCC      n_proj     = Number of bins to integrate in (i,j) to project onto (k)
CCC      num_pairs_sig = # pairs used in signal histogram
CCC      num_pairs_bkg = # pairs used in background histogram
CCC
CCC   Description of Output quantities:
CCC
CCC      c2mod_proj1,2,3(maxc2) = Reference histogram weighted 1D projections
CCC                               of C2 model function along {1,2,3} axes.
CCC      c2fit_proj1,2,3(maxc2) = Fitted 3D correlation function projected
CCC                               onto {1,2,3} axes.
CCC      c2err_proj1,2,3(maxc2) = Error in fitted 3D correlation function 
CCC                               projected onto {1,2,3} axes.

CCC   Local Variable Type Declarations:

      integer*4 maxh,maxc2,n,n_proj,num_pairs_sig,num_pairs_bkg
      integer*4 h(maxh,maxh,maxh),href(maxh,maxh,maxh)
      integer*4 i,j,k

      real*4 c2mod(maxc2,maxc2,maxc2)
      real*4 c2mod_proj1(maxc2),c2mod_proj2(maxc2),c2mod_proj3(maxc2)
      real*4 c2fit_proj1(maxc2),c2fit_proj2(maxc2),c2fit_proj3(maxc2)
      real*4 c2err_proj1(maxc2),c2err_proj2(maxc2),c2err_proj3(maxc2)
      real*4 a,a_error,b,b_error
      real*4 sum1n,sum1d,sum2n,sum2d,sum3n,sum3d

CCC   Initialize arrays to zero:

      do i = 1,maxc2
         c2mod_proj1(i) = 0.0
         c2mod_proj2(i) = 0.0
         c2mod_proj3(i) = 0.0
         c2fit_proj1(i) = 0.0
         c2fit_proj2(i) = 0.0
         c2fit_proj3(i) = 0.0
         c2err_proj1(i) = 0.0
         c2err_proj2(i) = 0.0
         c2err_proj3(i) = 0.0
      end do

CCC   Project Reference spectra (histogram) weighted model correlation:

      do i = 1,n
         sum1n = 0.0
         sum1d = 0.0
         sum2n = 0.0
         sum2d = 0.0
         sum3n = 0.0
         sum3d = 0.0
         do j = 1,n_proj
         do k = 1,n_proj
            sum1n = sum1n + c2mod(i,j,k)*float(href(i,j,k))
            sum1d = sum1d +              float(href(i,j,k))
            sum2n = sum2n + c2mod(k,i,j)*float(href(k,i,j))
            sum2d = sum2d +              float(href(k,i,j))
            sum3n = sum3n + c2mod(j,k,i)*float(href(j,k,i))
            sum3d = sum3d +              float(href(j,k,i))
         end do
         end do
         if(sum1d .le. 0.0) then
            c2mod_proj1(i) = 0.0
         else
            c2mod_proj1(i) = sum1n/sum1d
         end if
         if(sum2d .le. 0.0) then
            c2mod_proj2(i) = 0.0
         else
            c2mod_proj2(i) = sum2n/sum2d
         end if
         if(sum3d .le. 0.0) then
            c2mod_proj3(i) = 0.0
         else
            c2mod_proj3(i) = sum3n/sum3d
         end if
      end do

CCC   Calculate and Project the fitted correlation functions:

      do i = 1,n
         sum1n = 0.0
         sum1d = 0.0
         sum2n = 0.0
         sum2d = 0.0
         sum3n = 0.0
         sum3d = 0.0
         do j = 1,n_proj
         do k = 1,n_proj
            sum1n = sum1n + float(h(i,j,k))
            sum1d = sum1d + float(href(i,j,k))
            sum2n = sum2n + float(h(k,i,j))
            sum2d = sum2d + float(href(k,i,j))
            sum3n = sum3n + float(h(j,k,i))
            sum3d = sum3d + float(href(j,k,i))
         end do
         end do
         if(sum1n .le. 0.0 .or. sum1d .le. 0.0) then
            c2fit_proj1(i) = 0.0
            c2err_proj1(i) = 1.0
         else
            a = sum1n/float(num_pairs_sig)
            a_error = sqrt(sum1n)/float(num_pairs_sig)
            b = sum1d/float(num_pairs_bkg)
            b_error = sqrt(sum1d)/float(num_pairs_bkg)
            c2fit_proj1(i) = a/b
            c2err_proj1(i) = c2fit_proj1(i)*sqrt((a_error/a)**2
     1                       + (b_error/b)**2)
         end if
         if(sum2n .le. 0.0 .or. sum2d .le. 0.0) then
            c2fit_proj2(i) = 0.0
            c2err_proj2(i) = 1.0
         else
            a = sum2n/float(num_pairs_sig)
            a_error = sqrt(sum2n)/float(num_pairs_sig)
            b = sum2d/float(num_pairs_bkg)
            b_error = sqrt(sum2d)/float(num_pairs_bkg)
            c2fit_proj2(i) = a/b
            c2err_proj2(i) = c2fit_proj2(i)*sqrt((a_error/a)**2
     1                       + (b_error/b)**2)
         end if
         if(sum3n .le. 0.0 .or. sum3d .le. 0.0) then
            c2fit_proj3(i) = 0.0
            c2err_proj3(i) = 1.0
         else
            a = sum3n/float(num_pairs_sig)
            a_error = sqrt(sum3n)/float(num_pairs_sig)
            b = sum3d/float(num_pairs_bkg)
            b_error = sqrt(sum3d)/float(num_pairs_bkg)
            c2fit_proj3(i) = a/b
            c2err_proj3(i) = c2fit_proj3(i)*sqrt((a_error/a)**2
     1                       + (b_error/b)**2)
         end if
      end do

      Return
      END

C----------------------------------------------------------------------

C>>>>>>>>>>>>>>   Piotr, this needs to be replaced with Ali random
C>>>>>>>>>>>>>>>   number generator

*      real*4 function hbtpran(i)
*      implicit none
*      integer i
*      real*4 r
*      Call ranhbtp(r,1,i)
*      hbtpran = r
*      Return
*      END

*      Include 'ranlux2.f'
C----------------------------------------------------------------------




