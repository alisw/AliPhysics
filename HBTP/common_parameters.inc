CCC   common/parameters/
C
C     This common block includes the following data:
C
C          o Run control switches
C          o Correlation model selection switches
C          o Correlation model parameters
C          o Numbers of particles and pairs of particles
C          o Track momentum adjustment parameters
C          o Constants and Chi-Square values

      common/parameters/
     1   ALICE,errorcode,
     2   ref_control,print_full,print_sector_data,n_pid_types,pid(2),
     3   n_events,switch_1d,switch_3d,switch_type,switch_coherence,
     4   switch_coulomb,switch_fermi_bose,n_part_1_trk,n_part_2_trk,
     5   n_part_tot_trk,n_part_1_trk2,n_part_2_trk2,n_part_tot_trk2,
     6   n_part_used_1_trk,n_part_used_2_trk,n_part_used_1_trk2,
     7   n_part_used_2_trk2,n_part_used_1_ref,n_part_used_2_ref,
     8   n_part_used_1_inc,n_part_used_2_inc,num_pairs_like,
     9   num_pairs_unlike,num_pairs_like_ref,num_pairs_unlike_ref,
     1   num_pairs_like_inc,num_pairs_unlike_inc,event_line_counter,
     2   maxit,irand,file10_line_counter,
     3   lambda,R_1d,Rside,Rout,Rlong,Rperp,Rparallel,R0,Q0,deltap,
     4   delchi,pi,rad,hbc,chisq_wt_like_1d,chisq_wt_unlike_1d,
     5   chisq_wt_like_3d_fine,chisq_wt_unlike_3d_fine,
     6   chisq_wt_like_3d_coarse,chisq_wt_unlike_3d_coarse,
     7   chisq_wt_hist1_1,chisq_wt_hist1_2,mass1,mass2,trk_accep

CCC   Variable Type Declarations:

CCC   ALICE specific variables, DATA I/O option and return error code

      integer*4 ALICE                  ! To I/O from files or AliRoot strucs
      integer*4 errorcode              ! Processor error code to stop
CCC                                      calculation and return to calling 
CCC                                      program.     
CCC   Run Control and Model Switches:

      integer*4 ref_control            ! Select Reference Histograms
      integer*4 print_full             ! Full print out option - each event
      integer*4 print_sector_data      ! Print sector overflow diagnostics
      integer*4 n_pid_types            ! # particle ID types to correlate
      integer*4 pid                    ! Geant particle ID #s, max of 2 types
      integer*4 n_events               ! # events in input event text file
      integer*4 switch_1d              ! Include 1D correlations
      integer*4 switch_3d              ! Include 3D correlations
      integer*4 switch_type            ! For like, unlike or both PID pairs
      integer*4 switch_coherence       ! To include incoh/coher mixed source
      integer*4 switch_coulomb         ! Coulomb correction selection options
      integer*4 switch_fermi_bose      ! For fermions or bosons

CCC   Numbers of particles and pairs:

      integer*4 n_part_1_trk           ! Total # PID #1 in 'trk', all flags
      integer*4 n_part_2_trk           ! Total # PID #2 in 'trk', all flags
      integer*4 n_part_tot_trk         ! Total # all part. in 'trk', all flgs
      integer*4 n_part_used_1_trk      ! # PID#1, used (flag=0) in 'trk'
      integer*4 n_part_used_2_trk      ! # PID#2, used (flag=0) in 'trk'

      integer*4 n_part_1_trk2          ! Total # PID #1 in 'trk2', all flags
      integer*4 n_part_2_trk2          ! Total # PID #2 in 'trk2', all flags
      integer*4 n_part_tot_trk2        ! Total # all part. in 'trk2', all flgs
      integer*4 n_part_used_1_trk2     ! # PID#1, used (flag=0) in 'trk2'
      integer*4 n_part_used_2_trk2     ! # PID#2, used (flag=0) in 'trk2'

      integer*4 n_part_used_1_ref      ! # PID#1, used (flag=0) in Reference
      integer*4 n_part_used_2_ref      ! # PID#2, used (flag=0) in Reference
      integer*4 n_part_used_1_inc      ! # PID#1, used (flag=0) in Inclusive 
      integer*4 n_part_used_2_inc      ! # PID#2, used (flag=0) in Inclusive

      integer*4 num_pairs_like         ! # like pairs used (flag=0) in fit 
      integer*4 num_pairs_unlike       ! # unlike pairs used (flag=0) in fit
      integer*4 num_pairs_like_ref     ! # like pairs used (flag=0) in Ref. 
      integer*4 num_pairs_unlike_ref   ! # unlike pairs used (flag=0) in Ref. 
      integer*4 num_pairs_like_inc     ! # like pairs used (flag=0) in Incl. 
      integer*4 num_pairs_unlike_inc   ! # unlike pairs used (flag=0) in Incl. 

CCC   Counters:

      integer*4 event_line_counter     ! Input event text file line counter
      integer*4 file10_line_counter    ! Output, correlated event text file
CCC                                    !    line counter

CCC   Correlation Model Parameters:

      real*4    lambda                 ! Chaoticity parameter
      real*4    R_1d                   ! Spherical source radius (fm)
      real*4    Rside                  ! 3D Bertsch-Pratt source 'side' R (fm)
      real*4    Rout                   ! 3D Bertsch-Pratt source 'out'  R (fm)
      real*4    Rlong                  ! 3D Bertsch-Pratt source 'long' R (fm)
      real*4    Rperp                  ! 3D YKP source transverse radius  (fm)
      real*4    Rparallel              ! 3D YKP source longitudinal radius(fm)
      real*4    R0                     ! 3D YKP source emission time durat(fm)
      real*4    Q0                     ! NA35 Coulomb parameter (GeV/c) or
CCC                                    ! Coul radius for Pratt finite src (fm)

CCC   Search Control Parameters:

      integer*4 maxit                  ! Max # iterations in track adjustment
      real*4    deltap                 ! Max limit for x,y,z momt shifts(GeV/c)
      real*4    delchi                 ! Min% change in Chi-Sq to stop iterat.

CCC   Chi-Square Values:

      real*4    chisq_wt_like_1d          ! 1D, Like pairs
      real*4    chisq_wt_unlike_1d        ! 1D, Unlike pairs
      real*4    chisq_wt_like_3d_fine     ! 3D, Like pairs, Fine Mesh
      real*4    chisq_wt_unlike_3d_fine   ! 3D, Unlike pairs, Fine Mesh
      real*4    chisq_wt_like_3d_coarse   ! 3D, Like pairs, Coarse Mesh
      real*4    chisq_wt_unlike_3d_coarse ! 3D, Unlike pairs, Coarse Mesh
      real*4    chisq_wt_hist1_1          ! One-body, particle ID type #1
      real*4    chisq_wt_hist1_2          ! One-body, particle ID type #2

CCC   Particle Masses:

      real*4    mass1, mass2           ! Particle ID# 1 and 2 masses (GeV)

CCC   Constants:

      integer*4 irand                  ! Random # starting seed (Def=12345)
      real*4    pi, hbc                ! PI and h-bar-c (GeV*fm)
      real*4    rad                    ! radian = 180.0/pi
C     parameter (pi = 3.141592654)
C     parameter (hbc = 0.19732891)

CCC  Random Track Selection Fraction, for low multiplicity particles

      real*4 trk_accep                 ! ranges from 0.0 -> 1.0

CCC   END common/parameters/

