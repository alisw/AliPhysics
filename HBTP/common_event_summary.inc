CCC   common/event_summary/
C
C     This common block includes the following data:
C
C          o Event processing summary and diagnostic information
C          o Chi-Square values

      integer*4 max_events            ! Limit on stored values, this is NOT
C                                     ! a limit on the number of events that
C                                     ! can be processed in the input event
C                                     ! text file.
      parameter (max_events = 10000)

      common/event_summary/
     1num_iter(max_events), niter_mean, niter_rms,
     2n_part_used_1_store(max_events), npart1_mean, npart1_rms,
     3n_part_used_2_store(max_events), npart2_mean, npart2_rms,
     4n_part_tot_store(max_events), npart_tot_mean, npart_tot_rms,
     5num_sec_flagged_store(max_events), nsec_flag_mean, nsec_flag_rms,
     6frac_trks_out(max_events), frac_trks_out_mean,frac_trks_out_rms,
     7frac_trks_flag(max_events),frac_trks_flag_mean,frac_trks_flag_rms,
     8chisq_like_1d_store(max_events),         chi_l1d_mean,chi_l1d_rms,
     9chisq_unlike_1d_store(max_events),       chi_u1d_mean,chi_u1d_rms,
     1chisq_like_3d_fine_store(max_events),    chi_l3f_mean,chi_l3f_rms,
     2chisq_unlike_3d_fine_store(max_events),  chi_u3f_mean,chi_u3f_rms,
     3chisq_like_3d_coarse_store(max_events),  chi_l3c_mean,chi_l3c_rms,
     4chisq_unlike_3d_coarse_store(max_events),chi_u3c_mean,chi_u3c_rms,
     5chisq_hist1_1_store(max_events), chi_1_1_mean, chi_1_1_rms,
     6chisq_hist1_2_store(max_events), chi_1_2_mean, chi_1_2_rms,
     7chisq_total_store(max_events), chi_tot_mean, chi_tot_rms

CCC   Variable Type Declarations:

CCC   Event processing summary and diagnostic information:

      real*4 num_iter, niter_mean, niter_rms
      real*4 n_part_used_1_store, npart1_mean, npart1_rms
      real*4 n_part_used_2_store, npart2_mean, npart2_rms
      real*4 n_part_tot_store, npart_tot_mean, npart_tot_rms
      real*4 num_sec_flagged_store, nsec_flag_mean, nsec_flag_rms
      real*4 frac_trks_out, frac_trks_out_mean,frac_trks_out_rms
      real*4 frac_trks_flag,frac_trks_flag_mean,frac_trks_flag_rms

CCC   Chi-Square values:

      real*4 chisq_like_1d_store,         chi_l1d_mean,chi_l1d_rms
      real*4 chisq_unlike_1d_store,       chi_u1d_mean,chi_u1d_rms
      real*4 chisq_like_3d_fine_store,    chi_l3f_mean,chi_l3f_rms
      real*4 chisq_unlike_3d_fine_store,  chi_u3f_mean,chi_u3f_rms
      real*4 chisq_like_3d_coarse_store,  chi_l3c_mean,chi_l3c_rms
      real*4 chisq_unlike_3d_coarse_store,chi_u3c_mean,chi_u3c_rms
      real*4 chisq_hist1_1_store, chi_1_1_mean, chi_1_1_rms
      real*4 chisq_hist1_2_store, chi_1_2_mean, chi_1_2_rms
      real*4 chisq_total_store, chi_tot_mean, chi_tot_rms

CCC   END common/event_summary/
