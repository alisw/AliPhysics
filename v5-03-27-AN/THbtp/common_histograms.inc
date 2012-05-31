CCC   common/histograms/
C
C     This common block includes the following data:
C
C          o One- and Two-Body Histograms:
C   
C            One-Body:  hist1* - Originally filled and updated after each
C                                iteration where Chi-Square is improved.
C                       htmp1* - Working copy used within the track adjustment
C                                iteration loop
C                       href1* - Fixed reference distributions obtained from
C                                event text input file for all events.
C                       hinc1* - Inclusive event sum of all final, fitted
C                                one-body distributions
C
C                                (*) denotes pt,phi,eta and PID#1,2
C                                NOTE: The pt,phi,eta distributions are
C                                      projections of the full, 3D dist.
C                                      onto these axes.
C
C            Two-Body:  hist* - Originally filled and updated after each
C                               iteration where Chi-Square is improved.
C                       htmp* - Working copy used within the track adjustment
C                               iteration loop
C                       href* - Fixed reference distributions obtained from
C                               event text input file for all events.
C                       hinc* - Inclusive event sum of all final, fitted
C                               two-body distributions
C
C                               (*) denotes arrays for the following:
C                                   - Like pairs,   1D distributions
C                                   - Unlike pairs, 1D distributions           
C                                   - Like pairs,   3D dist., Fine Mesh
C                                   - Unlike pairs, 3D dist., Fine Mesh
C                                   - Like pairs,   3D dist., Coarse Mesh
C                                   - Unlike pairs, 3D dist., Coarse Mesh

      integer*4 max_h_1d, max_h_3d
      parameter (max_h_1d = 100)
      parameter (max_h_3d = 10 )

      common/histograms/
     1 hist1_pt_1(max_h_1d),hist1_phi_1(max_h_1d),hist1_eta_1(max_h_1d),
     2 hist1_pt_2(max_h_1d),hist1_phi_2(max_h_1d),hist1_eta_2(max_h_1d),
     3 htmp1_pt_1(max_h_1d),htmp1_phi_1(max_h_1d),htmp1_eta_1(max_h_1d),
     4 htmp1_pt_2(max_h_1d),htmp1_phi_2(max_h_1d),htmp1_eta_2(max_h_1d),
     5 href1_pt_1(max_h_1d),href1_phi_1(max_h_1d),href1_eta_1(max_h_1d),
     6 href1_pt_2(max_h_1d),href1_phi_2(max_h_1d),href1_eta_2(max_h_1d),
     7 hinc1_pt_1(max_h_1d),hinc1_phi_1(max_h_1d),hinc1_eta_1(max_h_1d),
     8 hinc1_pt_2(max_h_1d),hinc1_phi_2(max_h_1d),hinc1_eta_2(max_h_1d),
     9 hist_like_1d(max_h_1d),hist_unlike_1d(max_h_1d),
     1 htmp_like_1d(max_h_1d),htmp_unlike_1d(max_h_1d),
     2 href_like_1d(max_h_1d),href_unlike_1d(max_h_1d),
     3 hinc_like_1d(max_h_1d),hinc_unlike_1d(max_h_1d),
     4 hist_like_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     5 hist_unlike_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     6 hist_like_3d_coarse(max_h_3d,max_h_3d,max_h_3d),
     7 hist_unlike_3d_coarse(max_h_3d,max_h_3d,max_h_3d),
     4 htmp_like_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     5 htmp_unlike_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     6 htmp_like_3d_coarse(max_h_3d,max_h_3d,max_h_3d),
     7 htmp_unlike_3d_coarse(max_h_3d,max_h_3d,max_h_3d),
     4 href_like_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     5 href_unlike_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     6 href_like_3d_coarse(max_h_3d,max_h_3d,max_h_3d),
     7 href_unlike_3d_coarse(max_h_3d,max_h_3d,max_h_3d),
     4 hinc_like_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     5 hinc_unlike_3d_fine(max_h_3d,max_h_3d,max_h_3d),
     6 hinc_like_3d_coarse(max_h_3d,max_h_3d,max_h_3d),
     7 hinc_unlike_3d_coarse(max_h_3d,max_h_3d,max_h_3d)

CCC   Variable Type Declarations:

CCC   One-Body Distributions:

      integer*4 hist1_pt_1,hist1_phi_1,hist1_eta_1
      integer*4 hist1_pt_2,hist1_phi_2,hist1_eta_2
      integer*4 htmp1_pt_1,htmp1_phi_1,htmp1_eta_1
      integer*4 htmp1_pt_2,htmp1_phi_2,htmp1_eta_2
      integer*4 href1_pt_1,href1_phi_1,href1_eta_1
      integer*4 href1_pt_2,href1_phi_2,href1_eta_2
      integer*4 hinc1_pt_1,hinc1_phi_1,hinc1_eta_1
      integer*4 hinc1_pt_2,hinc1_phi_2,hinc1_eta_2

CCC   Two-Body Distributions 1D:

      integer*4 hist_like_1d, hist_unlike_1d
      integer*4 htmp_like_1d, htmp_unlike_1d
      integer*4 href_like_1d, href_unlike_1d
      integer*4 hinc_like_1d, hinc_unlike_1d

CCC   Two-Body Distributions 3D:

      integer*4 hist_like_3d_fine
      integer*4 hist_unlike_3d_fine
      integer*4 hist_like_3d_coarse
      integer*4 hist_unlike_3d_coarse
      integer*4 htmp_like_3d_fine
      integer*4 htmp_unlike_3d_fine
      integer*4 htmp_like_3d_coarse
      integer*4 htmp_unlike_3d_coarse
      integer*4 href_like_3d_fine
      integer*4 href_unlike_3d_fine
      integer*4 href_like_3d_coarse
      integer*4 href_unlike_3d_coarse
      integer*4 hinc_like_3d_fine
      integer*4 hinc_unlike_3d_fine
      integer*4 hinc_like_3d_coarse
      integer*4 hinc_unlike_3d_coarse

CCC   END common/histograms/
